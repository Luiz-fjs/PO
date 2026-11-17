"""Executa o modelo exato (ampl/cvrp.mod) em cada instância e salva os resultados.

Uso:
    .venv/bin/python src/run_exact.py --time-limit 3600 [--backend ampl|highspy] [instâncias...]

backend "ampl"    : AMPL (amplpy) + solver (--solver, padrão gurobi) lendo ampl/cvrp.mod e data/dat/<inst>.dat.
backend "highspy" : o MESMO modelo montado direto na API do HiGHS. Serve apenas de
                    reserva enquanto a licença AMPL é a demo (limite de 2000
                    variáveis/restrições); os resultados devem ser refeitos no AMPL.
"""
from __future__ import annotations

import argparse
import json
import math
import time
from pathlib import Path

from cvrp_io import DAT_DIR, ROOT, Instance, check_routes, read_instance, routes_cost
from instances import BY_NAME, EXACT_INSTANCES

MODEL = ROOT / "ampl" / "cvrp.mod"
SOLVER = "gurobi"   # sobrescrito por --solver
OUT_DIR = ROOT / "results" / "exact"


def arcs(inst: Instance) -> list[tuple[int, int]]:
    q, Q, n = inst.demand, inst.capacity, inst.n_customers
    return [(i, j) for i in range(n + 1) for j in range(n + 1)
            if i != j and (i == 0 or j == 0 or q[i] + q[j] <= Q)]


def succ_to_routes(succ: dict[int, list[int]]) -> list[list[int]]:
    routes = []
    for start in succ.get(0, []):
        route, cur = [], start
        while cur != 0:
            route.append(cur)
            cur = succ[cur][0]
        routes.append(route)
    return routes


def solve_ampl(inst: Instance, time_limit: float, threads: int) -> dict:
    from amplpy import AMPL

    ampl = AMPL()
    ampl.read(str(MODEL))
    ampl.read_data(str(DAT_DIR / f"{inst.name}.dat"))
    ampl.option["solver"] = SOLVER
    ampl.option[f"{SOLVER}_options"] = (
        f"outlev=0 timelim={time_limit} threads={threads} "
        "mip:gap=0 mip:gapabs=0.999 mip:bestbound=1 mip:return_gap=3"
    )
    n_vars = int(ampl.get_value("_nvars"))
    n_cons = int(ampl.get_value("_ncons"))
    t0 = time.perf_counter()
    ampl.solve()
    wall = time.perf_counter() - t0
    status = ampl.get_value("solve_result")
    solve_msg = ampl.get_value("solve_message")
    obj = ampl.get_objective("TotalDistance")
    has_sol = status in ("solved", "limit") and ampl.get_value("solve_result_num") < 500
    ub = obj.value() if has_sol else None
    try:
        lb = ampl.get_value("TotalDistance.bestbound")
    except Exception:
        lb = None
    routes = None
    if has_sol:
        x = ampl.get_variable("x").get_values().to_dict()
        succ: dict[int, list[int]] = {}
        for (i, j), v in x.items():
            if v > 0.5:
                succ.setdefault(int(i), []).append(int(j))
        routes = succ_to_routes(succ)
    return dict(backend=f"ampl+{SOLVER}", status=status, message=solve_msg, ub=ub, lb=lb,
                time=float(ampl.get_value("_solve_elapsed_time")), wall=wall,
                n_vars=n_vars, n_cons=n_cons, routes=routes)


def solve_highspy(inst: Instance, time_limit: float, threads: int) -> dict:
    """Mesmo modelo de ampl/cvrp.mod (fluxo de commodity única) na API do HiGHS."""
    import highspy
    import numpy as np

    n, Q, q = inst.n_customers, inst.capacity, inst.demand
    k_fixed = BY_NAME[inst.name].fixed_fleet
    kmin, kmax = (k_fixed, k_fixed) if k_fixed else (inst.min_vehicles, n)
    A = arcs(inst)
    idx = {a: k for k, a in enumerate(A)}
    nx = len(A)
    farcs = [a for a in A if a[1] != 0]
    fidx = {a: nx + k for k, a in enumerate(farcs)}
    ncol = nx + len(farcs)
    cost = [float(inst.dist[i][j]) for (i, j) in A] + [0.0] * len(farcs)
    lower = [0.0] * ncol
    upper = [1.0] * nx + [highspy.kHighsInf] * len(farcs)
    integ = [highspy.HighsVarType.kInteger] * nx + [highspy.HighsVarType.kContinuous] * len(farcs)

    rows_lo, rows_hi, starts, index, value = [], [], [0], [], []

    def add_row(coefs: dict[int, float], lo: float, hi: float) -> None:
        for k, v in coefs.items():
            index.append(k); value.append(v)
        starts.append(len(index)); rows_lo.append(lo); rows_hi.append(hi)

    out = {i: [] for i in range(n + 1)}
    inc = {i: [] for i in range(n + 1)}
    for (i, j), k in idx.items():
        out[i].append(k); inc[j].append(k)
    for i in range(1, n + 1):
        add_row({k: 1.0 for k in out[i]}, 1.0, 1.0)            # LeaveOnce
    for j in range(1, n + 1):
        add_row({k: 1.0 for k in inc[j]}, 1.0, 1.0)            # EnterOnce
    bal = {k: 1.0 for k in out[0]}
    for k in inc[0]:
        bal[k] = bal.get(k, 0.0) - 1.0
    add_row(bal, 0.0, 0.0)                                     # DepotBalance
    add_row({k: 1.0 for k in out[0]}, float(kmin), float(kmax))  # FleetSize
    fin = {i: [] for i in range(n + 1)}
    fout = {i: [] for i in range(n + 1)}
    for (i, j), k in fidx.items():
        fout[i].append(k); fin[j].append(k)
    for i in range(1, n + 1):                                  # FlowBalance
        coefs = {k: 1.0 for k in fin[i]}
        coefs.update({k: -1.0 for k in fout[i]})
        add_row(coefs, float(q[i]), float(q[i]))
    for (i, j), k in fidx.items():
        add_row({k: 1.0, idx[(i, j)]: -float(q[j])}, 0.0, highspy.kHighsInf)          # FlowLower
        add_row({k: 1.0, idx[(i, j)]: -float(Q - q[i])}, -highspy.kHighsInf, 0.0)     # FlowUpper

    lp = highspy.HighsLp()
    lp.num_col_, lp.num_row_ = ncol, len(rows_lo)
    lp.col_cost_, lp.col_lower_, lp.col_upper_ = np.array(cost), np.array(lower), np.array(upper)
    lp.row_lower_, lp.row_upper_ = np.array(rows_lo), np.array(rows_hi)
    lp.a_matrix_.format_ = highspy.MatrixFormat.kRowwise
    lp.a_matrix_.start_ = np.array(starts, dtype=np.int32)
    lp.a_matrix_.index_ = np.array(index, dtype=np.int32)
    lp.a_matrix_.value_ = np.array(value)
    lp.integrality_ = integ

    h = highspy.Highs()
    h.setOptionValue("output_flag", False)
    h.setOptionValue("time_limit", float(time_limit))
    h.setOptionValue("threads", threads)
    h.setOptionValue("mip_rel_gap", 0.0)
    h.setOptionValue("mip_abs_gap", 0.999)   # custos inteiros: gap < 1 prova o ótimo
    h.passModel(lp)
    t0 = time.perf_counter()
    h.run()
    wall = time.perf_counter() - t0
    info = h.getInfo()
    model_status = h.modelStatusToString(h.getModelStatus())
    has_sol = info.primal_solution_status == 2   # kSolutionStatusFeasible
    ub = info.objective_function_value if has_sol else None
    lb = info.mip_dual_bound
    routes = None
    if has_sol:
        col = h.getSolution().col_value
        succ: dict[int, list[int]] = {}
        for (i, j), k in idx.items():
            if col[k] > 0.5:
                succ.setdefault(i, []).append(j)
        routes = succ_to_routes(succ)
    return dict(backend="highspy", status=model_status, message=model_status, ub=ub, lb=lb,
                time=h.getRunTime(), wall=wall, n_vars=ncol, n_cons=len(rows_lo), routes=routes)


def run(name: str, backend: str, time_limit: float, threads: int) -> dict:
    inst = read_instance(name)
    solver = solve_ampl if backend == "ampl" else solve_highspy
    res = solver(inst, time_limit, threads)
    res.update(instance=name, tier=BY_NAME[name].tier, n=inst.n_customers,
               bks=inst.bks, bks_is_optimal=BY_NAME[name].bks_is_optimal,
               time_limit=time_limit, threads=threads)
    # validação independente da solução retornada pelo solver
    if res["routes"] is not None:
        check_routes(inst, res["routes"])
        recomputed = routes_cost(inst, res["routes"])
        assert abs(recomputed - res["ub"]) < 0.5, (recomputed, res["ub"])
        res["ub"] = recomputed
        res["n_routes"] = len(res["routes"])
    if res["lb"] is not None and (math.isinf(res["lb"]) or math.isnan(res["lb"])):
        res["lb"] = None
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    (OUT_DIR / f"{name}.json").write_text(json.dumps(res, indent=1))
    return res


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("instances", nargs="*")
    ap.add_argument("--backend", choices=["ampl", "highspy", "auto"], default="ampl",
                    help="auto: AMPL no conjunto E (cabe na licença demo), highspy nos demais")
    ap.add_argument("--time-limit", type=float, default=3600)
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--solver", choices=["gurobi", "highs", "cplex"], default="gurobi")
    args = ap.parse_args()
    SOLVER = args.solver
    names = args.instances or [i.name for i in EXACT_INSTANCES]
    for name in names:
        backend = args.backend
        if backend == "auto":
            backend = "ampl" if name.startswith("E-") else "highspy"
        r = run(name, backend, args.time_limit, args.threads)
        print(f"{name:12s} {r['backend']:10s} status={r['status']:<10} UB={r['ub']} "
              f"LB={r['lb']} t={r['time']:.1f}s vars={r['n_vars']} cons={r['n_cons']}", flush=True)
