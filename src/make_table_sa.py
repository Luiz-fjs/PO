"""Tabela comparativa método exato x Simulated Annealing (Etapa 6).

Lê results/exact_results.csv (gerado por make_table.py) e results/sa/<inst>/seed*.json.
Gera:
    results/tabela_comparacao.md    tabela em Markdown
    results/tabela_comparacao.tex   tabela em LaTeX (para o template SBPO)
    results/comparacao.csv          dados agregados
    results/figs/desvio_por_n.pdf   desvio (%) x n, exato e SA
    results/figs/convergencia_<inst>.pdf
    results/figs/rotas_<inst>.pdf
e imprime o teste de Wilcoxon pareado (desvio do SA médio x desvio do exato).
"""
from __future__ import annotations

import csv
import json
from statistics import mean, pstdev

from cvrp_io import ROOT, read_instance
from instances import INSTANCES

RES = ROOT / "results"
SA_DIR = RES / "sa"
OUT = RES
FIGS = OUT / "figs"


def br(v, nd=2):
    if v is None:
        return "—"
    return f"{v:,.{nd}f}".replace(",", "X").replace(".", ",").replace("X", ".")


def load_exact() -> dict[str, dict]:
    p = RES / "exact_results.csv"   # sempre lido de results/
    out = {}
    if p.exists():
        for r in csv.DictReader(p.open()):
            out[r["instancia"]] = r
    return out


def load_sa(name: str) -> list[dict]:
    d = SA_DIR / name
    return [json.loads(p.read_text()) for p in sorted(d.glob("seed*.json"))] if d.exists() else []


def main() -> None:
    exact = load_exact()
    rows = []
    for info in INSTANCES:
        runs = load_sa(info.name)
        if not runs:
            continue
        bks = runs[0]["bks"]
        costs = [r["cost"] for r in runs]
        ex = exact.get(info.name)
        ex_ub = float(ex["ub"]) if ex and ex["ub"] not in ("", "None") else None
        ex_gap = float(ex["gap"]) if ex and ex["gap"] not in ("", "None") else None
        rows.append(dict(
            faixa=info.tier, instancia=info.name, n=runs[0]["n"], bks=bks,
            bks_otimo=info.bks_is_optimal,
            ex_ub=ex_ub, ex_gap=ex_gap,
            ex_dev=(ex_ub - bks) / bks * 100 if ex_ub is not None else None,
            ex_t=float(ex["tempo"]) if ex else None,
            ex_status=ex["status"] if ex else "—",
            sa_best=min(costs), sa_mean=mean(costs), sa_std=pstdev(costs),
            sa_dev_best=(min(costs) - bks) / bks * 100,
            sa_dev_mean=(mean(costs) - bks) / bks * 100,
            sa_tbest=mean(r["time_to_best"] for r in runs),
            sa_t=mean(r["time"] for r in runs),
            runs=len(runs),
        ))

    with open(OUT / "comparacao.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    # ---------- Markdown ----------
    hdr = ("| Faixa | Instância | n | Ótimo/BKS | Exato: UB | Exato: gap (%) | Exato: desvio (%) | Exato: tempo (s) "
           "| SA: melhor | SA: média ± dp | SA: desvio melhor (%) | SA: desvio médio (%) | SA: t. até melhor (s) | SA: tempo (s) |")
    md = [hdr, "|" + "---|" * 14]
    for r in rows:
        star = "" if r["bks_otimo"] else "*"
        md.append(
            f"| {r['faixa']} | {r['instancia']} | {r['n']} | {br(r['bks'], 0)}{star} | {br(r['ex_ub'], 0)} | "
            f"{br(r['ex_gap'])} | {br(r['ex_dev'])} | {br(r['ex_t'], 1)} | {br(r['sa_best'], 0)} | "
            f"{br(r['sa_mean'], 1)} ± {br(r['sa_std'], 1)} | {br(r['sa_dev_best'])} | {br(r['sa_dev_mean'])} | "
            f"{br(r['sa_tbest'], 1)} | {br(r['sa_t'], 1)} |")
    md.append("")
    md.append("\\* BKS sem otimalidade provada. Desvio = (custo − Ótimo/BKS)/(Ótimo/BKS). "
              f"SA: {rows[0]['runs']} execuções independentes por instância.")

    # ---------- Wilcoxon ----------
    paired = [(r["ex_dev"], r["sa_dev_mean"]) for r in rows if r["ex_dev"] is not None]
    stat_txt = ""
    try:
        from scipy.stats import wilcoxon
        if len(paired) >= 5:
            st = wilcoxon([a for a, _ in paired], [b for _, b in paired])
            stat_txt = (f"Teste de Wilcoxon pareado (desvio do exato x desvio médio do SA, "
                        f"{len(paired)} instâncias): W = {st.statistic:.1f}, p = {st.pvalue:.4f}")
    except ImportError:
        pass
    if stat_txt:
        md += ["", stat_txt]
        print(stat_txt)
    (OUT / "tabela_comparacao.md").write_text("\n".join(md) + "\n")

    # ---------- LaTeX (compacta, cabe na largura do template SBPO) ----------
    tex = [r"\begin{table}[ht]", r"\centering", r"\scriptsize", r"\setlength{\tabcolsep}{3pt}",
           r"\caption{Comparação entre o método exato (AMPL + Gurobi, 3600 s) e o SA "
           f"({rows[0]['runs']} execuções).}}",
           r"\label{tab:comparacao}",
           r"\begin{tabular}{lrrrrrrrrrr}", r"\hline",
           r" & & & \multicolumn{3}{c}{Exato} & \multicolumn{5}{c}{Simulated Annealing} \\",
           r"Instância & $n$ & Ót./BKS & UB & Gap\% & $t$(s) & Melhor & Média & Dev.m\% & Dev.µ\% & $t^*$(s) \\",
           r"\hline"]
    tier = None
    for r in rows:
        if r["faixa"] != tier:
            if tier is not None:
                tex.append(r"\hline")
            tier = r["faixa"]
        star = "" if r["bks_otimo"] else "$^*$"
        tex.append(
            f"{r['instancia']} & {r['n']} & {br(r['bks'], 0)}{star} & {br(r['ex_ub'], 0)} & "
            f"{br(r['ex_gap'])} & {br(r['ex_t'], 0)} & {br(r['sa_best'], 0)} & {br(r['sa_mean'], 1)} & "
            f"{br(r['sa_dev_best'])} & {br(r['sa_dev_mean'])} & {br(r['sa_tbest'], 1)} \\\\")
    tex += [r"\hline", r"\end{tabular}", r"\end{table}"]
    (OUT / "tabela_comparacao.tex").write_text("\n".join(tex).replace("µ", r"$\mu$") + "\n")

    figures(rows)
    print("\n".join(md))


def figures(rows: list[dict]) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    FIGS.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})

    # desvio x n
    fig, ax = plt.subplots(figsize=(6, 3.2))
    ex = [(r["n"], r["ex_dev"]) for r in rows if r["ex_dev"] is not None and r["ex_dev"] < 100]
    ax.plot([a for a, _ in ex], [b for _, b in ex], "o-", color="#c0504d", label="Exato (UB em 3600 s)")
    ax.plot([r["n"] for r in rows], [r["sa_dev_mean"] for r in rows], "s-", color="#1f77b4",
            label="SA (média)")
    ax.fill_between([r["n"] for r in rows], [r["sa_dev_best"] for r in rows],
                    [(r["sa_mean"] + r["sa_std"] - r["bks"]) / r["bks"] * 100 for r in rows],
                    color="#1f77b4", alpha=0.15, label="SA (melhor a média + dp)")
    ax.set_xscale("log")
    ax.set_xlabel("número de clientes $n$ (escala log)")
    ax.set_ylabel("desvio do ótimo/BKS (%)")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(FIGS / "desvio_por_n.pdf")
    fig.savefig(FIGS / "desvio_por_n.png", dpi=200)
    plt.close(fig)

    # convergência e rotas da melhor execução numa instância média e numa grande
    for name in ("A-n80-k10", "X-n303-k21"):
        runs = load_sa(name)
        if not runs:
            continue
        best = min(runs, key=lambda r: r["cost"])
        fig, ax = plt.subplots(figsize=(4.5, 2.8))
        for r in runs:
            t, c = zip(*r["trace"])
            ax.step(t, c, where="post", lw=0.8, alpha=0.6)
        ax.axhline(best["bks"], ls="--", color="k", lw=0.8, label="ótimo/BKS")
        ax.set_xscale("symlog", linthresh=1)
        ax.set_xlabel("tempo (s)")
        ax.set_ylabel("melhor custo viável")
        ax.set_title(name)
        ax.legend(frameon=False)
        fig.tight_layout()
        fig.savefig(FIGS / f"convergencia_{name}.pdf")
        plt.close(fig)

        inst = read_instance(name, with_matrix=False)
        if inst.coords:
            fig, ax = plt.subplots(figsize=(4, 4))
            cmap = plt.get_cmap("tab20")
            for k, route in enumerate(best["routes"]):
                pts = [inst.coords[0]] + [inst.coords[c] for c in route] + [inst.coords[0]]
                ax.plot([p[0] for p in pts], [p[1] for p in pts], "-o", ms=2, lw=0.8, color=cmap(k % 20))
            ax.plot(*inst.coords[0], "ks", ms=6)
            ax.set_title(f"{name}: custo {best['cost']} ({len(best['routes'])} rotas)")
            ax.set_xticks([]); ax.set_yticks([])
            fig.tight_layout()
            fig.savefig(FIGS / f"rotas_{name}.pdf")
            plt.close(fig)


if __name__ == "__main__":
    import argparse
    from pathlib import Path

    ap = argparse.ArgumentParser()
    ap.add_argument("--sa-dir", type=Path, default=SA_DIR)
    ap.add_argument("--out", type=Path, default=OUT)
    a = ap.parse_args()
    SA_DIR, OUT, FIGS = a.sa_dir, a.out, a.out / "figs"
    OUT.mkdir(parents=True, exist_ok=True)
    main()
