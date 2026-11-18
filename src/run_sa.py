"""Executa o Simulated Annealing em cada instância com várias sementes (Etapa 6).

Uso:
    python src/run_sa.py                         # todas as instâncias, sementes 1..10
    python src/run_sa.py -j 4 --seeds 1-10 A-n32-k5 X-n101-k25
    python src/run_sa.py --time-limit 30 ...     # sobrescreve t_max = min(3600, max(60, n))

Cada execução é salva em results/sa/<instância>/seed<k>.json (execuções já
existentes são puladas, então o script pode ser interrompido e retomado).
Parâmetros do SA: os padrões de SAParams, ou um JSON com --params
(ex.: results/calibration/best_params.json, gerado por calibrate_sa.py).
"""
from __future__ import annotations

import argparse
import json
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import fields
from pathlib import Path

from cvrp_io import ROOT, read_instance
from instances import BY_NAME, INSTANCES
from sa import SAParams, solve

OUT_DIR = ROOT / "results" / "sa"


def fleet_of(name: str) -> int | None:
    if name in BY_NAME:
        return BY_NAME[name].fixed_fleet
    return int(name.rsplit("-k", 1)[1]) if name[0] in "EA" else None


def parse_seeds(text: str) -> list[int]:
    out: list[int] = []
    for part in text.split(","):
        if "-" in part:
            a, b = part.split("-")
            out += range(int(a), int(b) + 1)
        else:
            out.append(int(part))
    return out


def run_one(name: str, seed: int, overrides: dict, out_dir: str) -> str:
    out = Path(out_dir) / name / f"seed{seed}.json"
    if out.exists():
        return f"{name} seed {seed}: já existe"
    inst = read_instance(name)
    params = SAParams(seed=seed, **overrides)
    res = solve(inst, fleet_of(name), params)
    res.update(instance=name, seed=seed, n=inst.n_customers, bks=inst.bks)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(res))
    dev = 100 * (res["cost"] - inst.bks) / inst.bks if inst.bks else float("nan")
    return (f"{name} seed {seed}: {res['cost']} (desvio {dev:.2f}%) "
            f"t_best={res['time_to_best']:.1f}s t={res['time']:.1f}s")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("instances", nargs="*")
    ap.add_argument("--seeds", default="1-10")
    ap.add_argument("-j", "--jobs", type=int, default=max(1, (os.cpu_count() or 2) // 2))
    ap.add_argument("--time-limit", type=float, default=None)
    ap.add_argument("--params", type=Path, default=None, help="JSON com parâmetros do SA")
    ap.add_argument("--out", type=Path, default=OUT_DIR)
    a = ap.parse_args()

    names = a.instances or [i.name for i in INSTANCES]
    overrides: dict = {}
    if a.params:
        valid = {f.name for f in fields(SAParams)}
        overrides = {k: v for k, v in json.loads(a.params.read_text()).items()
                     if k in valid and k not in ("seed", "time_limit")}
    if a.time_limit is not None:
        overrides["time_limit"] = a.time_limit

    # instâncias maiores primeiro, para equilibrar a fila
    jobs = sorted(((n, s) for n in names for s in parse_seeds(a.seeds)),
                  key=lambda x: -read_instance(x[0], with_matrix=False).n_customers)
    print(f"{len(jobs)} execuções, {a.jobs} em paralelo, parâmetros extras: {overrides}")
    with ProcessPoolExecutor(max_workers=a.jobs) as ex:
        futs = [ex.submit(run_one, n, s, overrides, str(a.out)) for n, s in jobs]
        for f in as_completed(futs):
            print(f.result(), flush=True)


if __name__ == "__main__":
    main()
