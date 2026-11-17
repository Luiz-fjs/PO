"""Agrega results/exact/*.json numa tabela (Markdown e CSV)."""
from __future__ import annotations

import csv
import json

from cvrp_io import ROOT
from instances import EXACT_INSTANCES

OUT_DIR = ROOT / "results"


def fmt(v, nd=0):
    if v is None:
        return "—"
    return f"{v:,.{nd}f}".replace(",", "X").replace(".", ",").replace("X", ".")


def status_pt(r: dict) -> str:
    if r["ub"] is not None and r["lb"] is not None and r["ub"] - r["lb"] < 1 - 1e-6:
        return "Ótimo"
    if r["ub"] is not None:
        return "Viável (limite de tempo)"
    if "memory" in str(r.get("message", "")).lower():
        return "Sem solução (memória)"
    return "Sem solução viável"


def main() -> None:
    rows = []
    for info in EXACT_INSTANCES:
        p = ROOT / "results" / "exact" / f"{info.name}.json"
        if not p.exists():
            continue
        r = json.loads(p.read_text())
        ub, lb, bks = r["ub"], r["lb"], r["bks"]
        lb = None if lb is None else max(0.0, lb)
        gap = (ub - lb) / ub * 100 if ub is not None and lb is not None else None
        dev = (ub - bks) / bks * 100 if ub is not None else None
        rows.append(dict(
            faixa=info.tier, instancia=info.name, n=r["n"], bks=bks,
            bks_tipo="ótimo" if info.bks_is_optimal else "BKS",
            ub=ub, lb=lb, gap=gap, desvio=dev, tempo=r["time"], status=status_pt(r),
            rotas=r.get("n_routes"), variaveis=r["n_vars"], restricoes=r["n_cons"],
            backend=r["backend"]))

    with open(OUT_DIR / "exact_results.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)

    head = ("| Faixa | Instância | n | Ótimo/BKS | Solução (UB) | Limite inferior (LB) "
            "| Gap (%) | Desvio BKS (%) | Tempo (s) | Status | Variáveis | Restrições |")
    lines = [head, "|" + "---|" * 12]
    for r in rows:
        star = "" if r["bks_tipo"] == "ótimo" else "*"
        lines.append(
            f"| {r['faixa']} | {r['instancia']} | {r['n']} | {fmt(r['bks'])}{star} | {fmt(r['ub'])} "
            f"| {fmt(r['lb'])} | {fmt(r['gap'], 2)} | {fmt(r['desvio'], 2)} | {fmt(r['tempo'], 1)} "
            f"| {r['status']} | {fmt(r['variaveis'])} | {fmt(r['restricoes'])} |")
    (OUT_DIR / "tabela_exato.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
