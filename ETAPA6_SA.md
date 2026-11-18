# Etapa 6 – Simulated Annealing (implementação em Python)

Código: `src/sa.py` (algoritmo), `src/run_sa.py` (experimentos), `src/calibrate_sa.py`
(calibração), `src/make_table_sa.py` (tabela comparativa, teste de Wilcoxon e figuras).
Python puro (sem dependências para o SA; `scipy` e `matplotlib` só para tabela/figuras).

## Como rodar (no mesmo computador do método exato)

```bash
# 1) calibração (baixar antes A-n33-k5, A-n54-k7, X-n125-k30 .vrp/.sol do CVRPLIB para data/raw/)
.venv/bin/python src/calibrate_sa.py -j 6 --time-limit 60 --seeds 1-3        # ~1 h

# 2) bateria final: 15 instâncias x 10 sementes, t_max = min(3600, max(60, n)) s
.venv/bin/python src/run_sa.py -j 6 --params results/calibration/best_params.json   # ~3-4 h
#    (pode ser interrompido e retomado: execuções prontas são puladas)

# 3) tabela exato x SA (Markdown + LaTeX), Wilcoxon e figuras em results/figs/
.venv/bin/python src/make_table_sa.py
```

Um `-j` igual ao número de núcleos de desempenho do M2 (4) deixa a comparação de tempo mais justa
com o exato, que rodou 4 instâncias em paralelo.

## O que mudou em relação à proposta (Etapa 4) — descrever assim no artigo

| Item | Proposta | Implementado | Motivo (observado nos testes) |
|---|---|---|---|
| Resfriamento | geométrico por patamares, α ∈ {0,95; 0,98; 0,99}, L = βn | geométrico **em função do tempo**: T(τ) = T_s·(T_f/T_s)^τ, τ = fração decorrida do ciclo | com α por patamar, T caía de T0 a ~0 em < 1 s e o SA virava descida gulosa; amarrar ao tempo garante que todo o orçamento seja usado no recozimento |
| Reaquecimento | após R patamares sem melhora, T ← 0,5·T0 | o orçamento é dividido em C ciclos; cada ciclo recomeça da melhor solução com T_s = ρ·T0 | mesmo efeito, controlado pelo tempo |
| T0 | piora média de F | piora média só da **distância** | a penalidade inflava T0 em instâncias apertadas |
| Vizinhanças | relocate, swap, 2-opt, 2-opt* | + **or-opt** (bloco de 2–3 clientes, com inversão opcional) e abertura de rota nova (frota livre) | ganho claro em A e X |
| λ | ×/÷ 1,2 a cada 100 it. | ×/÷ 1,05 a cada max(100, n) it. | o ajuste brusco "esquentava" a busca mesmo com T baixo |
| Fim de ciclo | — | **reparo + polimento**: descida estrita com λ alto (≤ 5% do ciclo) | nas X grandes (carga média 98,6% da capacidade) o fim do ciclo terminava inviável |

Parâmetros padrão (até a calibração): p0 = 0,05; T_f = 0,01·T0; C = 5 ciclos; ρ = 0,1; k = 20 vizinhos.

## Prévia (20 s por execução, 2 sementes, máquina de 2 núcleos – NÃO são os resultados finais)

| Instância | Exato (3600 s): desvio | SA (20 s): desvio médio |
|---|---:|---:|
| E-n13 … A-n32 | 0,00% (ótimo) | 0,00% (ótimo) |
| A-n45-k6 | 4,34% | 1,22% |
| A-n60-k9 | 1,18% | 0,07% |
| A-n80-k10 | 2,16% | 0,34% |
| X-n101-k25 | 1,47% | 0,29% |
| X-n153-k22 | 5,97% | 2,63% |
| X-n200-k36 | 4,88% | 1,12% |
| X-n303-k21 | 9,06% | 2,36% |
| X-n502-k39 | 6,75% | 1,47% |
| X-n1001-k43 | 1800,72% | 5,37% |
| Leuven1 | — | 3,76% |

Toda solução final do SA passa por `check_routes` (cada cliente uma vez, capacidade) e o custo é
recalculado do zero.
