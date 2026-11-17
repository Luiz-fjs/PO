# Proposta de Solução e Implementação: CVRP por Método Exato e Simulated Annealing

## 1. Visão geral

O CVRP será resolvido por duas abordagens, aplicadas às mesmas instâncias e no mesmo computador:

1. **Método exato**: modelo de Programação Linear Inteira Mista (PLIM) escrito em AMPL e resolvido pelo Gurobi via *branch-and-cut*, com limite de tempo. Ele retorna a solução ótima ou, quando o tempo acaba, a melhor solução viável (limitante superior, UB) e o melhor limitante inferior (LB), o que permite calcular o **gap de otimalidade**.
2. **Metaheurística**: *Simulated Annealing* (SA) com vizinhanças múltiplas, penalização de capacidade, listas granulares de vizinhos e reaquecimento. O SA é avaliado pelo **desvio em relação ao ótimo ou BKS** publicado no CVRPLIB.

A hipótese do trabalho, vinda da revisão bibliográfica, é a seguinte: o método exato prova a otimalidade nas instâncias pequenas, mas perde rapidamente a capacidade de fechar o gap (e até de encontrar soluções viáveis) à medida que *n* cresce. Já o SA mantém desvios baixos em tempo controlado em todas as faixas.

---

## 2. Instâncias de teste

Todas as instâncias vêm do CVRPLIB (<https://galgos.inf.puc-rio.br/cvrplib/en/cvrp/instances>, consultado em 16/11/2026), que fornece coordenadas, demandas, capacidade e a melhor solução conhecida (BKS). A distância é euclidiana arredondada ao inteiro mais próximo (convenção TSPLIB `EUC_2D`). A exceção é o E-n13-k4, cuja matriz de distâncias é dada explicitamente (`EXPLICIT / LOWER_ROW`).

Atendendo ao pedido do professor, foram acrescentadas duas faixas acima de 200 clientes: **muito grandes** (X com 300 a 1000 clientes, o maior tamanho do conjunto X) e **massiva** (Leuven1, do conjunto AGS, com 3000 clientes).

| Faixa | Instância | Clientes (n) | Q | Frota | Ótimo/BKS | Ótimo provado? |
|---|---|---:|---:|---|---:|:---:|
| Pequena | E-n13-k4 | 12 | 6000 | = 4 | 247 | sim |
| Pequena | E-n22-k4 | 21 | 6000 | = 4 | 375 | sim |
| Pequena | E-n23-k3 | 22 | 4500 | = 3 | 569 | sim |
| Pequena | E-n30-k3 | 29 | 4500 | = 3 | 534 | sim |
| Média | A-n32-k5 | 31 | 100 | = 5 | 784 | sim |
| Média | A-n45-k6 | 44 | 100 | = 6 | 944 | sim |
| Média | A-n60-k9 | 59 | 100 | = 9 | 1.354 | sim |
| Média | A-n80-k10 | 79 | 100 | = 10 | 1.763 | sim |
| Grande | X-n101-k25 | 100 | 206 | ≥ 25 | 27.591 | sim |
| Grande | X-n153-k22 | 152 | 144 | ≥ 22 | 21.220 | sim |
| Grande | X-n200-k36 | 199 | 402 | ≥ 36 | 58.578 | sim |
| **Muito grande** | **X-n303-k21** | 302 | 794 | ≥ 21 | 21.736 | não |
| **Muito grande** | **X-n502-k39** | 501 | 13 | ≥ 39 | 69.226 | não |
| **Muito grande** | **X-n1001-k43** | 1000 | 131 | ≥ 43 | 72.332 | não |
| **Massiva** | **Leuven1** | 3000 | 25 | ≥ 203 | 192.848 | não |

Fontes: n, Q, BKS e status de otimalidade vêm da tabela de instâncias do CVRPLIB. Os arquivos `.vrp` e `.sol` estão em `data/raw/`. O custo de cada BKS foi recalculado a partir das rotas do `.sol` com o leitor implementado (`src/cvrp_io.py`) e bateu exatamente com o valor publicado nas 15 instâncias, o que valida o leitor e o cálculo de distâncias.

**Tamanho da frota.** Os conjuntos clássicos (E, A) usam exatamente K veículos, como indicam os comentários dos arquivos ("No of trucks: 5"). No conjunto X a frota é livre, com mínimo ⌈Σq/Q⌉ (Uchoa et al., 2017). Essa distinção importa na prática:

- No E-n30-k3 com frota livre, o solver encontrou uma solução de custo 503, menor que o ótimo publicado de 534, que usa 3 veículos.
- O BKS do X-n101-k25 usa 26 rotas, e não 25.

O modelo trata os dois casos com os parâmetros `Kmin` e `Kmax`.

**Para os exatos ainda maiores.** Leuven1 entra **apenas no SA**. O modelo compacto teria cerca de 9 milhões de variáveis binárias, inviável nesta máquina (8 GB de RAM). O CVRPLIB ainda lista os conjuntos AGS completo (até 30.000 clientes) e XL, caso se queira estender a faixa massiva.

---

## 3. Método exato: formulação matemática

### 3.1 Notação

- $V = \{0, 1, \dots, n\}$: nós, sendo 0 o depósito; $V_c = V \setminus \{0\}$: clientes.
- $q_i > 0$: demanda do cliente $i$ (com $q_0 = 0$); $Q$: capacidade do veículo.
- $c_{ij}$: distância entre $i$ e $j$.
- $K_{\min}, K_{\max}$: limites da frota ($=K$ nos conjuntos E/A; $\lceil \sum_i q_i / Q \rceil$ e $n$ no conjunto X).
- $A = \{(i,j) \in V \times V : i \neq j,\ (i = 0 \lor j = 0 \lor q_i + q_j \le Q)\}$: arcos viáveis. Os arcos entre dois clientes cuja soma de demandas excede $Q$ são eliminados *a priori*.

**Variáveis**

- $x_{ij} \in \{0,1\}$: vale 1 se algum veículo percorre o arco $(i,j)$.
- $f_{ij} \ge 0$: carga a bordo do veículo ao percorrer $(i,j)$, definida para $j \neq 0$, pois o veículo volta vazio ao depósito.

### 3.2 Modelo: fluxo de veículos com dois índices e fluxo de commodity única (Gavish & Graves, 1978)

$$
\min \sum_{(i,j) \in A} c_{ij}\, x_{ij}
$$

sujeito a

$$
\begin{aligned}
&\sum_{j:(i,j)\in A} x_{ij} = 1 && \forall i \in V_c && (1)\ \text{um sucessor}\\
&\sum_{i:(i,j)\in A} x_{ij} = 1 && \forall j \in V_c && (2)\ \text{um predecessor}\\
&\sum_{j:(0,j)\in A} x_{0j} = \sum_{i:(i,0)\in A} x_{i0} && && (3)\ \text{saídas = retornos}\\
&K_{\min} \le \sum_{j:(0,j)\in A} x_{0j} \le K_{\max} && && (4)\ \text{frota}\\
&\sum_{h:(h,i)\in A} f_{hi} - \sum_{j:(i,j)\in A,\, j\neq 0} f_{ij} = q_i && \forall i \in V_c && (5)\ \text{conservação de fluxo}\\
&q_j\, x_{ij} \le f_{ij} \le (Q - q_i)\, x_{ij} && \forall (i,j) \in A,\ j \neq 0 && (6)\ \text{acoplamento e capacidade}\\
&x_{ij} \in \{0,1\},\quad f_{ij} \ge 0
\end{aligned}
$$

**Interpretação.**

- As restrições (1) e (2) garantem que cada cliente seja visitado exatamente uma vez.
- As restrições (3) e (4) garantem que as rotas comecem e terminem no depósito e respeitem a frota.
- A restrição (5) faz cada cliente "retirar" sua demanda do veículo.
- A restrição (6) só permite fluxo em arcos usados e limita a carga em $Q$.

Como $q_i > 0$ para todo cliente (verificado nas 15 instâncias), um ciclo formado apenas por clientes exigiria que a carga diminuísse a cada nó e voltasse ao valor inicial, o que é impossível. Assim, (5) e (6) eliminam subciclos e impõem a capacidade ao mesmo tempo.

**Tamanho do modelo.** São $|A| \approx n^2$ variáveis binárias, $\approx n^2$ contínuas e $\approx 2n^2 + 3n$ restrições. Na instância X-n1001-k43 isso dá cerca de 2 milhões de variáveis.

### 3.3 Por que fluxo de commodity única e não MTZ

Também foi implementada a formulação clássica com restrições MTZ de carga (Miller et al., 1960) com o *lifting* de Kara et al. (2004), em `ampl/cvrp_mtz.mod`. As duas foram comparadas no conjunto E com 120 s de limite (arquivo `results/formulation_comparison.txt`):

| Instância | Ótimo | LP MTZ | LP fluxo | Tempo MTZ | Tempo fluxo |
|---|---:|---:|---:|---:|---:|
| E-n13-k4 | 247 | 226,0 (gap 8,5%) | 238,3 (gap 3,5%) | 0,6 s | 0,6 s |
| E-n22-k4 | 375 | 308,5 (17,7%) | 330,6 (11,8%) | > 120 s (LB 348) | **5,0 s** |
| E-n23-k3 | 569 | 515,0 (9,5%) | 489,8 (13,9%) | 19,9 s | **1,0 s** |
| E-n30-k3 | 534 | 347,0 (35,0%) | 433,6 (18,8%) | > 120 s (UB 601, LB 363) | > 120 s (UB 555, LB 470) |

Essa comparação foi feita com o HiGHS, antes da licença acadêmica. A formulação de fluxo tem relaxação linear mais forte em 3 das 4 instâncias e resolve bem mais rápido, ao custo de cerca de 2× mais variáveis. Por isso ela é o modelo principal (`ampl/cvrp.mod`).

### 3.4 Implementação e configuração

| Item | Configuração |
|---|---|
| Linguagem de modelagem | AMPL (versão 20261109, via `amplpy` 0.18.1), licença AMPL for Academics |
| Solver | Gurobi 13.0.3 (*branch-and-cut*), licença AMPL for Academics |
| Limite de tempo | 3600 s por instância |
| Critério de otimalidade | UB − LB < 1 (custos inteiros), ou seja, `mip:gapabs=0.999` |
| Threads | 1 por instância; 4 instâncias em paralelo |
| Máquina | Apple M2, 8 núcleos, 8 GB de RAM, macOS 26 |
| Validação | as rotas retornadas pelo solver são reconstruídas, conferidas (cada cliente exatamente uma vez, capacidade) e o custo é recalculado de forma independente |

**Métricas do método exato** (por instância):

- solução encontrada (UB);
- limitante inferior (LB);
- **gap** = (UB − LB)/UB;
- **desvio** = (UB − BKS)/BKS;
- tempo de CPU;
- status: ótimo, viável no limite de tempo ou sem solução.

---

## 4. Metaheurística: Simulated Annealing

O SA (Kirkpatrick et al., 1983) parte de uma solução inicial e, a cada iteração, gera um vizinho. O vizinho é aceito se melhora a função de avaliação, ou com probabilidade $e^{-\Delta/T}$ se piora. A temperatura $T$ é reduzida aos poucos. O desenho abaixo incorpora as conclusões da revisão: o desempenho depende da vizinhança e do resfriamento (Arce-Ortega et al., 2025), e o reaquecimento ajuda a escapar de regiões pouco promissoras.

### 4.1 Representação

A solução é uma lista de rotas, cada uma com a sequência de clientes, a carga acumulada e o custo da rota. Também é mantido o vetor `rota[i]` / `posição[i]`, que permite avaliar movimentos em O(1).

### 4.2 Solução inicial

A solução inicial vem do algoritmo de economias de Clarke & Wright (1964), na versão paralela.

- **Frota livre (conjunto X e Leuven1):** a solução de Clarke & Wright é usada diretamente.
- **Frota fixa (conjuntos E e A):** se o algoritmo gerar mais de K rotas, as menores são dissolvidas e seus clientes reinseridos nas K restantes. A capacidade violada é tratada pela penalidade da seção 4.4.

### 4.3 Vizinhanças

A cada iteração, um dos movimentos abaixo é sorteado (probabilidades iguais, a calibrar):

| Movimento | Descrição |
|---|---|
| *Relocate* | Remove o cliente $i$ e o reinsere após um cliente $j$, na mesma rota ou em outra. |
| *Swap* | Troca as posições de $i$ e $j$, na mesma rota ou entre rotas. |
| *2-opt* (intra) | Inverte o trecho entre $i$ e $j$ de uma mesma rota. |
| *2-opt\** (inter) | Troca as "caudas" de duas rotas a partir de $i$ e $j$. |

**Listas granulares.** Nas instâncias com centenas ou milhares de clientes, avaliar pares $(i,j)$ aleatórios desperdiça iterações com clientes distantes. Por isso, $j$ é sorteado entre os $k=20$ vizinhos mais próximos de $i$, a ideia de vizinhança granular de Toth & Vigo (2003). Isso mantém o custo de cada iteração constante e independente de *n*, o que é o que viabiliza X-n1001 e Leuven1.

**Frota.** Com frota livre, um movimento *relocate* para uma rota vazia permite abrir uma rota nova, e rotas esvaziadas são removidas. Com frota fixa, movimentos que esvaziem uma rota são proibidos.

### 4.4 Função de avaliação

A busca pode passar por soluções inviáveis, avaliadas por

$$
F(s) = \sum_{r \in s} \text{custo}(r) + \lambda \sum_{r \in s} \max\{0,\ \text{carga}(r) - Q\}.
$$

O peso $\lambda$ é adaptativo. Ele começa igual à distância média entre um cliente e seus vizinhos granulares. A cada 100 iterações, $\lambda$ é multiplicado por 1,2 se menos de 50% das soluções visitadas forem viáveis, e dividido por 1,2 caso contrário. A melhor solução registrada é sempre **viável**.

### 4.5 Esquema de temperatura

| Parâmetro | Definição | Valores para calibração |
|---|---|---|
| $T_0$ | $T_0 = -\bar\Delta^+ / \ln p_0$, sendo $\bar\Delta^+$ a piora média de 1000 movimentos aleatórios a partir da solução inicial | $p_0 \in \{0{,}3;\ 0{,}5\}$ |
| Resfriamento | geométrico, $T_{k+1} = \alpha T_k$ | $\alpha \in \{0{,}95;\ 0{,}98;\ 0{,}99\}$ |
| Iterações por temperatura | $L = \beta \cdot n$ | $\beta \in \{1, 5, 10\}$ |
| Reaquecimento | se a melhor solução não melhorar em $R$ patamares consecutivos, $T \leftarrow 0{,}5\,T_0$ e a busca recomeça da melhor solução | $R \in \{20, 50\}$ |
| Parada | tempo limite $t_{\max} = \min(3600,\ \max(60,\ n))$ segundos, ou 10 reaquecimentos sem melhora | — |

O tempo do SA nunca passa do limite dado ao método exato (3600 s), o que torna a comparação justa. Os parâmetros serão calibrados num subconjunto que não entra nos resultados finais: A-n33-k5, A-n54-k7 e X-n125-k30, com 3 execuções por configuração. A configuração vencedora é então fixada para todas as instâncias.

### 4.6 Pseudocódigo

```
s ← ClarkeWright(); s* ← s se viável; T ← T0; sem_melhora ← 0
enquanto tempo < t_max e reaquecimentos < 10:
    repetir L vezes:
        m ← sorteia movimento; (i, j) ← i aleatório, j ∈ vizinhos_granulares(i)
        Δ ← avaliação incremental de m em F
        se Δ < 0 ou rand() < exp(−Δ/T): aplica m em s
        se s viável e custo(s) < custo(s*): s* ← s; sem_melhora ← 0
        ajusta λ a cada 100 iterações
    T ← α·T; sem_melhora ← sem_melhora + 1
    se sem_melhora ≥ R: T ← 0,5·T0; s ← s*; sem_melhora ← 0; reaquecimentos += 1
retorna s*
```

### 4.7 Protocolo experimental do SA

- 10 execuções independentes por instância (sementes 1 a 10), 1 thread, mesma máquina do método exato.
- Métricas:
  - melhor custo, custo médio e desvio-padrão;
  - **desvio** = (custo − BKS)/BKS, para o melhor e para a média;
  - tempo médio até a melhor solução e tempo total.
- Cada solução final é validada pelo mesmo verificador usado no método exato (`check_routes`).
- Comparação final por faixa: desvio do SA × gap/desvio do exato × tempo.

---

## 5. Resultados do método exato

Configuração: AMPL 20261109 + Gurobi 13.0.3, modelo `ampl/cvrp.mod`, limite de 3600 s por instância, 1 thread, 4 instâncias em paralelo num Apple M2 com 8 GB de RAM. Dados brutos em `results/exact/*.json` e `results/exact_results.csv`.

**Tabela 1 – Resultados do método exato por instância**

| Faixa | Instância | n | Ótimo/BKS | Solução (UB) | Limite inferior (LB) | Gap (%) | Desvio BKS (%) | Tempo (s) | Status | Variáveis | Restrições |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Pequena | E-n13-k4 | 12 | 247 | 247 | 247 | 0,00 | 0,00 | 0,3 | Ótimo | 300 | 326 |
| Pequena | E-n22-k4 | 21 | 375 | 375 | 375 | 0,00 | 0,00 | 4,9 | Ótimo | 903 | 947 |
| Pequena | E-n23-k3 | 22 | 569 | 569 | 569 | 0,00 | 0,00 | 0,5 | Ótimo | 970 | 1.016 |
| Pequena | E-n30-k3 | 29 | 534 | 534 | 534 | 0,00 | 0,00 | 1.089,8 | Ótimo | 1.707 | 1.767 |
| Média | A-n32-k5 | 31 | 784 | 784 | 784 | 0,00 | 0,00 | 770,3 | Ótimo | 1.953 | 2.017 |
| Média | A-n45-k6 | 44 | 944 | 985 | 902 | 8,43 | 4,34 | 3.600,6 | Limite de tempo | 3.916 | 4.006 |
| Média | A-n60-k9 | 59 | 1.354 | 1.370 | 1.285 | 6,20 | 1,18 | 3.600,2 | Limite de tempo | 7.021 | 7.141 |
| Média | A-n80-k10 | 79 | 1.763 | 1.801 | 1.689 | 6,22 | 2,16 | 3.600,2 | Limite de tempo | 12.561 | 12.721 |
| Grande | X-n101-k25 | 100 | 27.591 | 27.996 | 26.276 | 6,14 | 1,47 | 3.600,2 | Limite de tempo | 20.100 | 20.302 |
| Grande | X-n153-k22 | 152 | 21.220 | 22.486 | 19.625 | 12,72 | 5,97 | 3.600,5 | Limite de tempo | 44.940 | 45.246 |
| Grande | X-n200-k36 | 199 | 58.578 | 61.434 | 57.079 | 7,09 | 4,88 | 3.600,7 | Limite de tempo | 79.401 | 79.801 |
| Muito grande | X-n303-k21 | 302 | 21.736* | 23.706 | 20.148 | 15,01 | 9,06 | 3.601,0 | Limite de tempo | 182.710 | 183.316 |
| Muito grande | X-n502-k39 | 501 | 69.226* | 73.902 | 67.071 | 9,24 | 6,75 | 3.645,3 | Limite de tempo | 502.503 | 503.507 |
| Muito grande | X-n1001-k43 | 1000 | 72.332* | 1.374.831 | 0 | — | 1.800,72 | 3.611,5 | Limite de tempo | 2.001.000 | 2.003.002 |

\* Melhor solução conhecida (BKS), sem otimalidade provada na literatura (CVRPLIB). Nas demais, o valor é o ótimo.
Gap = (UB − LB)/UB. Desvio = (UB − Ótimo/BKS)/(Ótimo/BKS). No X-n1001-k43 o Gurobi não obteve limitante inferior útil (LB = 0) em 3600 s, por isso o gap não é informado.

**Análise.**

- **Pequenas:** o método exato prova a otimalidade em todas as instâncias. As de 12 a 22 clientes levam menos de 5 s, mas o E-n30-k3 já exige 1.090 s.
- **Médias:** a otimalidade só é provada na A-n32-k5, em 770 s. Da A-n45-k6 em diante, todas param no limite de tempo, com gap entre 6,2% e 8,4% e solução de 1,2% a 4,3% acima do ótimo.
- **Grandes (100 a 199 clientes):** sempre há solução viável, com desvio de 1,5% a 6,0% e gap de 6,1% a 12,7%.
- **Muito grandes:**
  - X-n303-k21 e X-n502-k39 terminam com desvio de 6,8% a 9,1% em relação ao BKS.
  - No X-n1001-k43, o modelo tem cerca de 2 milhões de variáveis. O Gurobi não obtém limitante inferior em 1 hora, e a solução encontrada é cerca de 19 vezes o BKS.
- **Comparação com o HiGHS:** numa rodada preliminar com o mesmo modelo e o mesmo limite de tempo (`results/highs_partial/`), o HiGHS não encontrou solução viável para X-n502-k39 e X-n1001-k43, e para X-n303-k21 parou em 218.766, cerca de 10 vezes o BKS.
- **Conclusão:** os resultados confirmam a hipótese do trabalho. O método exato é adequado até cerca de 30 clientes e perde rapidamente a capacidade de provar o ótimo à medida que *n* cresce, o que abre espaço para o Simulated Annealing nas faixas maiores.


---

## 6. Reprodução

```bash
python3 -m venv .venv && .venv/bin/pip install amplpy highspy numpy
.venv/bin/python -m amplpy.modules install highs gurobi
.venv/bin/python src/make_dat.py              # gera data/dat/*.dat
BACKEND=ampl SOLVER=gurobi TL=3600 ./run_all_exact.sh   # roda todas as instâncias
.venv/bin/python src/make_table.py            # gera results/tabela_exato.md e .csv
```

Para rodar uma instância isolada no próprio AMPL, use `ampl/cvrp.run`.

**Licença AMPL.** A licença demo instalada pelo `amplpy` limita o modelo a 2000 variáveis e 2000 restrições. Os resultados foram obtidos com a licença **AMPL for Academics**, gratuita para estudantes (<https://portal.ampl.com>), que remove o limite e libera o Gurobi. Essa licença é temporária (a atual expira em 16/12/2026). Para ativar, use `.venv/bin/python -m amplpy.modules activate <UUID>`, sem colocar o UUID em arquivos do projeto. Sem a licença, use `SOLVER=highs`.

---

## Referências

- Arce-Ortega, D. V. et al. (2025). A simulated annealing approach for the homogeneous capacitated vehicle routing problem. *Mathematics*, 13(19), 3209. https://doi.org/10.3390/math13193209
- Arnold, F., Gendreau, M., & Sörensen, K. (2019). Efficiently solving very large-scale routing problems. *Computers & Operations Research*, 107, 32–42. https://doi.org/10.1016/j.cor.2019.03.006
- Clarke, G., & Wright, J. W. (1964). Scheduling of vehicles from a central depot to a number of delivery points. *Operations Research*, 12(4), 568–581. https://doi.org/10.1287/opre.12.4.568
- Gavish, B., & Graves, S. C. (1978). The travelling salesman problem and related problems. Working Paper OR-078-78, Operations Research Center, MIT. http://hdl.handle.net/1721.1/5363
- Gurobi Optimization, LLC. *Gurobi Optimizer Reference Manual*, versão 13.0. https://docs.gurobi.com
- Huangfu, Q., & Hall, J. A. J. (2017). Parallelizing the dual revised simplex method. *Mathematical Programming Computation*, 10(1), 119–142. https://doi.org/10.1007/s12532-017-0130-5 (HiGHS)
- Kara, I., Laporte, G., & Bektas, T. (2004). A note on the lifted Miller–Tucker–Zemlin subtour elimination constraints for the capacitated vehicle routing problem. *European Journal of Operational Research*, 158(3), 793–795. https://doi.org/10.1016/S0377-2217(03)00377-1
- Kirkpatrick, S., Gelatt, C. D., & Vecchi, M. P. (1983). Optimization by simulated annealing. *Science*, 220(4598), 671–680. https://doi.org/10.1126/science.220.4598.671
- Miller, C. E., Tucker, A. W., & Zemlin, R. A. (1960). Integer programming formulation of traveling salesman problems. *Journal of the ACM*, 7(4), 326–329. https://doi.org/10.1145/321043.321046
- Toth, P., & Vigo, D. (2003). The granular tabu search and its application to the vehicle-routing problem. *INFORMS Journal on Computing*, 15(4), 333–346. https://doi.org/10.1287/ijoc.15.4.333.24890
- Uchoa, E. et al. (2017). New benchmark instances for the capacitated vehicle routing problem. *European Journal of Operational Research*, 257(3), 845–858. https://doi.org/10.1016/j.ejor.2016.08.012
- CVRPLIB. Capacitated Vehicle Routing Problem Library. https://galgos.inf.puc-rio.br/cvrplib/ (acesso em 16/11/2026).
- Demais referências conforme a revisão bibliográfica (Pessoa et al., 2020; Errami et al., 2024; Ilhan, 2020, 2021; Fitzpatrick et al., 2024).
