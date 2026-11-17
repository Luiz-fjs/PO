#!/bin/zsh
# Roda o método exato em todas as instâncias, 4 em paralelo, 1 thread cada.
#   BACKEND=ampl (padrão): tudo via AMPL; exige licença Community Edition
#                          para instâncias com mais de 2000 variáveis.
#   BACKEND=auto         : AMPL no conjunto E, API do HiGHS nos demais.
cd "$(dirname "$0")/src"
TL=${TL:-3600}
BACKEND=${BACKEND:-ampl}
SOLVER=${SOLVER:-gurobi}
mkdir -p ../results/logs
printf '%s\n' X-n1001-k43 X-n502-k39 X-n303-k21 X-n200-k36 X-n153-k22 X-n101-k25 \
  A-n80-k10 A-n60-k9 A-n45-k6 A-n32-k5 E-n30-k3 E-n23-k3 E-n22-k4 E-n13-k4 |
  xargs -P 4 -I{} sh -c "../.venv/bin/python run_exact.py --backend $BACKEND --solver $SOLVER --time-limit $TL {} > ../results/logs/{}.log 2>&1"
