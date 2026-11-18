"""Calibração dos parâmetros do SA (fatorial completo), como na seção 4.5 da proposta.

As instâncias de calibração NÃO entram nos resultados finais:
A-n33-k5, A-n54-k7 e X-n125-k30 (baixar .vrp e .sol do CVRPLIB para data/raw/).

Uso:
    python src/calibrate_sa.py -j 6 --time-limit 60 --seeds 1-3

Saídas:
    results/calibration/runs.csv         uma linha por (configuração, instância, semente)
    results/calibration/ranking.csv      configurações ordenadas pelo desvio médio
    results/calibration/best_params.json melhor configuração (usar em run_sa.py --params)
"""
from __future__ import annotations

import argparse
import csv
import itertools
import json
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from statistics import mean

from cvrp_io import ROOT, read_instance
from run_sa import fleet_of, parse_seeds
from sa import SAParams, solve

CALIB_INSTANCES = ["A-n33-k5", "A-n54-k7", "X-n125-k30"]
GRID = {
    "p0": [0.02, 0.05, 0.2],
    "tf_ratio": [1e-2, 1e-3],
    "n_cycles": [1, 5, 10],
    "reheat_factor": [0.1, 0.5],
}
OUT = ROOT / "results" / "calibration"


def run(cfg: dict, name: str, seed: int, tl: float) -> dict:
    inst = read_instance(name)
    res = solve(inst, fleet_of(name), SAParams(seed=seed, time_limit=tl, **cfg))
    return dict(**cfg, instance=name, seed=seed, cost=res["cost"], bks=inst.bks,
                dev=100 * (res["cost"] - inst.bks) / inst.bks, t_best=res["time_to_best"])


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("-j", "--jobs", type=int, default=max(1, (os.cpu_count() or 2) - 2))
    ap.add_argument("--time-limit", type=float, default=60)
    ap.add_argument("--seeds", default="1-3")
    a = ap.parse_args()

    from cvrp_io import RAW_DIR
    names = [n for n in CALIB_INSTANCES if (RAW_DIR / f"{n}.vrp").exists() and (RAW_DIR / f"{n}.sol").exists()]
    missing = sorted(set(CALIB_INSTANCES) - set(names))
    if missing:
        print("AVISO: instâncias de calibração ausentes (rode download_calibration.py):", missing)
    if not names:
        raise SystemExit("nenhuma instância de calibração em data/raw/")
    keys = list(GRID)
    configs = [dict(zip(keys, vals)) for vals in itertools.product(*GRID.values())]
    tasks = [(c, n, s) for c in configs for n in names for s in parse_seeds(a.seeds)]
    print(f"{len(configs)} configurações x {len(names)} instâncias -> {len(tasks)} execuções")
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    with ProcessPoolExecutor(max_workers=a.jobs) as ex:
        futs = [ex.submit(run, c, n, s, a.time_limit) for c, n, s in tasks]
        for k, f in enumerate(as_completed(futs), 1):
            rows.append(f.result())
            if k % 20 == 0:
                print(f"{k}/{len(tasks)}", flush=True)

    with open(OUT / "runs.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    ranking = []
    for c in configs:
        devs = [r["dev"] for r in rows if all(r[k] == c[k] for k in keys)]
        ranking.append(dict(**c, mean_dev=mean(devs), max_dev=max(devs)))
    ranking.sort(key=lambda r: (r["mean_dev"], r["max_dev"]))
    with open(OUT / "ranking.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(ranking[0]))
        w.writeheader()
        w.writerows(ranking)
    best = {k: ranking[0][k] for k in keys}
    (OUT / "best_params.json").write_text(json.dumps(best, indent=2))
    print("melhor configuração:", best, f"desvio médio {ranking[0]['mean_dev']:.2f}%")


if __name__ == "__main__":
    main()
