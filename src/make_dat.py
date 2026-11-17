"""Gera os arquivos .dat do AMPL a partir das instâncias CVRPLIB em data/raw."""
from __future__ import annotations

import sys

from cvrp_io import DAT_DIR, Instance, read_instance
from instances import BY_NAME, EXACT_INSTANCES


def fleet_lines(inst: Instance) -> list[str]:
    k = BY_NAME[inst.name].fixed_fleet
    return [f"param Kmin := {k};", f"param Kmax := {k};"] if k else []


def write_dat(inst: Instance) -> None:
    DAT_DIR.mkdir(parents=True, exist_ok=True)
    n = inst.n_customers
    lines = [f"# {inst.name} (CVRPLIB) - BKS = {inst.bks}",
             f"param n := {n};",
             f"param Q := {inst.capacity};",
             *fleet_lines(inst),
             "param q :=",
             *(f"  {i} {inst.demand[i]}" for i in range(n + 1)),
             ";",
             "param c :="]
    # formato (i, j, valor) por linha: simples e robusto para matrizes grandes
    for i in range(n + 1):
        row = inst.dist[i]
        lines.append(" ".join(f"[{i},{j}] {row[j]}" for j in range(n + 1) if j != i))
    lines.append(";")
    (DAT_DIR / f"{inst.name}.dat").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    names = sys.argv[1:] or [i.name for i in EXACT_INSTANCES]
    for name in names:
        write_dat(read_instance(name))
        print(f"wrote {DAT_DIR / (name + '.dat')}")
