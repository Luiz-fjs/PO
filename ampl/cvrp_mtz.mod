# Formulação ALTERNATIVA (usada só na comparação de força da relaxação).
# CVRP - formulação de fluxo de veículos com dois índices (arcos direcionados)
# e restrições MTZ de carga com lifting de Desrochers-Laporte (Kara et al., 2004).
#
# Nó 0 = depósito; 1..n = clientes.

param n integer > 0;                    # número de clientes
set C := 1..n;                          # clientes
set N := 0..n;                          # clientes + depósito

param Q > 0;                            # capacidade do veículo
param q {N} >= 0;                       # demanda (q[0] = 0)
param c {N, N} >= 0;                    # distância (EUC_2D arredondada)

# Número mínimo de veículos: ceil(demanda total / Q). A frota é ilimitada,
# convenção do CVRPLIB para os conjuntos X e AGS (Uchoa et al., 2017).
param Kmin integer > 0 default ceil(sum {i in C} q[i] / Q);
param Kmax integer >= Kmin default n;

# Arcos viáveis: arco entre dois clientes só existe se q[i] + q[j] <= Q.
set A := {i in N, j in N: i != j and (i = 0 or j = 0 or q[i] + q[j] <= Q)};

var x {A} binary;                       # 1 se o veículo vai de i para j
var u {i in C} >= q[i], <= Q;           # carga acumulada ao sair do cliente i

minimize TotalDistance:
    sum {(i, j) in A} c[i, j] * x[i, j];

# Cada cliente tem exatamente um sucessor e um predecessor.
subject to LeaveOnce {i in C}:
    sum {(i, j) in A} x[i, j] = 1;

subject to EnterOnce {j in C}:
    sum {(i, j) in A} x[i, j] = 1;

# Todo veículo que sai do depósito retorna a ele.
subject to DepotBalance:
    sum {(0, j) in A} x[0, j] = sum {(i, 0) in A} x[i, 0];

subject to FleetSize:
    Kmin <= sum {(0, j) in A} x[0, j] <= Kmax;

# MTZ de carga com lifting: elimina subciclos sem o depósito e impõe a capacidade.
subject to LoadMTZ {(i, j) in A: i != 0 and j != 0}:
    u[j] >= u[i] + q[j] - Q * (1 - x[i, j]) + (Q - q[i] - q[j]) * x[j, i];
