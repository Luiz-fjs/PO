"""Compara o SA com 1 ciclo (results/sa) e com 5 ciclos (results/sa_5cycles).

Pareia as execuções pela semente e aplica o teste de Wilcoxon por instância,
nas médias das 15 instâncias e em todos os pares. Valida as rotas de 5 ciclos.

Uso:
    .venv/bin/python src/compare_cycles.py > results/cycles_comparison.txt
"""
from __future__ import annotations

import json
import statistics as st

from scipy.stats import wilcoxon

from cvrp_io import ROOT, check_routes, read_instance, routes_cost
from instances import INSTANCES


def load_costs(folder: str, name: str, inst, validate: bool) -> dict[int, int]:
    costs = {}
    for p in (ROOT / "results" / folder / name).glob("seed*.json"):
        r = json.loads(p.read_text())
        if validate:
            check_routes(inst, r["routes"])
            assert routes_cost(inst, r["routes"]) == r["cost"]
        costs[r["seed"]] = r["cost"]
    return costs


def main() -> None:
    pairs, inst_means = [], []
    print(f"{'instance':12s} {'1c dev%':>8s} {'5c dev%':>8s} {'5c<1c':>5s} {'5c>1c':>5s} "
          f"{'1c best':>8s} {'5c best':>8s} {'p':>6s}")
    for info in INSTANCES:
        inst = read_instance(info.name)
        c1 = load_costs("sa", info.name, inst, validate=False)
        c5 = load_costs("sa_5cycles", info.name, inst, validate=True)
        seeds = sorted(set(c1) & set(c5))
        a, b = [c1[s] for s in seeds], [c5[s] for s in seeds]
        dev = lambda v: (v - inst.bks) / inst.bks * 100
        pairs += [(dev(x), dev(y)) for x, y in zip(a, b)]
        inst_means.append((dev(st.mean(a)), dev(st.mean(b))))
        p = f"{wilcoxon(a, b).pvalue:.3f}" if a != b else "-"
        print(f"{info.name:12s} {dev(st.mean(a)):8.2f} {dev(st.mean(b)):8.2f} "
              f"{sum(y < x for x, y in zip(a, b)):5d} {sum(y > x for x, y in zip(a, b)):5d} "
              f"{min(a):8d} {min(b):8d} {p:>6s}")
    x, y = zip(*inst_means)
    print(f"\nmean deviation over instances: 1 cycle {st.mean(x):.3f}%  5 cycles {st.mean(y):.3f}%")
    print("Wilcoxon, 15 instance means:", wilcoxon(x, y, zero_method="zsplit"))
    nz = [(u, v) for u, v in pairs if u != v]
    print(f"Wilcoxon, seed pairs ({len(nz)} non-tied of {len(pairs)}):",
          wilcoxon([u for u, _ in nz], [v for _, v in nz]))


if __name__ == "__main__":
    main()
