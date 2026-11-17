# CVRP - formulação de fluxo de veículos com dois índices (arcos direcionados)
# e fluxo de commodity única (Gavish & Graves, 1978).
#
# Nó 0 = depósito; 1..n = clientes.

param n integer > 0;                    # número de clientes
set C := 1..n;                          # clientes
set N := 0..n;                          # clientes + depósito

param Q > 0;                            # capacidade do veículo
param q {N} >= 0;                       # demanda (q[0] = 0)
param c {N, N} >= 0;                    # distância (EUC_2D arredondada)

# Frota: conjuntos E e A usam exatamente K veículos (Kmin = Kmax = K);
# conjunto X tem frota livre (Kmin = ceil(demanda/Q), Kmax = n).
param Kmin integer > 0 default ceil(sum {i in C} q[i] / Q);
param Kmax integer >= Kmin default n;

# Arcos viáveis: arco entre dois clientes só existe se q[i] + q[j] <= Q.
set A := {i in N, j in N: i != j and (i = 0 or j = 0 or q[i] + q[j] <= Q)};

var x {A} binary;                       # 1 se algum veículo percorre o arco (i,j)
var f {(i, j) in A: j != 0} >= 0;       # carga a bordo no arco (i,j); volta vazio ao depósito

minimize TotalDistance:
    sum {(i, j) in A} c[i, j] * x[i, j];

# (1)-(2) cada cliente tem exatamente um sucessor e um predecessor
subject to LeaveOnce {i in C}:
    sum {(i, j) in A} x[i, j] = 1;

subject to EnterOnce {j in C}:
    sum {(i, j) in A} x[i, j] = 1;

# (3) todo veículo que sai do depósito retorna a ele
subject to DepotBalance:
    sum {(0, j) in A} x[0, j] = sum {(i, 0) in A} x[i, 0];

# (4) tamanho da frota
subject to FleetSize:
    Kmin <= sum {(0, j) in A} x[0, j] <= Kmax;

# (5) conservação de fluxo: cada cliente retém exatamente sua demanda
subject to FlowBalance {i in C}:
    sum {(h, i) in A} f[h, i] - sum {(i, j) in A: j != 0} f[i, j] = q[i];

# (6)-(7) acoplamento fluxo-arco: elimina subciclos e impõe a capacidade
subject to FlowLower {(i, j) in A: j != 0}:
    f[i, j] >= q[j] * x[i, j];

subject to FlowUpper {(i, j) in A: j != 0}:
    f[i, j] <= (Q - q[i]) * x[i, j];
