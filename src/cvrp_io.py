"""Leitura de instâncias CVRPLIB (formato TSPLIB) e utilidades comuns.

Suporta EDGE_WEIGHT_TYPE = EUC_2D (distância euclidiana arredondada ao inteiro
mais próximo, convenção TSPLIB/CVRPLIB) e EXPLICIT/LOWER_ROW (caso do E-n13-k4).
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RAW_DIR = ROOT / "data" / "raw"
DAT_DIR = ROOT / "data" / "dat"


@dataclass
class Instance:
    name: str
    capacity: int
    demand: list[int]                      # índice 0 = depósito
    dist: list[list[int]]                  # matriz (n+1) x (n+1)
    coords: list[tuple[float, float]] | None = None
    bks: int | None = None                 # custo da melhor solução conhecida
    bks_routes: list[list[int]] = field(default_factory=list)

    @property
    def n_customers(self) -> int:
        return len(self.demand) - 1

    @property
    def min_vehicles(self) -> int:
        # limitante trivial: ceil(demanda total / capacidade)
        return math.ceil(sum(self.demand) / self.capacity)


def nint(value: float) -> int:
    """Arredondamento ao inteiro mais próximo (TSPLIB): floor(x + 0.5)."""
    return int(math.floor(value + 0.5))


def _parse_header_value(line: str) -> str:
    return line.split(":", 1)[1].strip()


def read_instance(name: str, raw_dir: Path = RAW_DIR, with_matrix: bool = True) -> Instance:
    lines = (raw_dir / f"{name}.vrp").read_text().splitlines()
    header: dict[str, str] = {}
    coords: dict[int, tuple[float, float]] = {}
    demand: dict[int, int] = {}
    weights: list[int] = []
    depot = 1
    section = None
    for raw in lines:
        line = raw.strip()
        if not line or line == "EOF":
            continue
        upper = line.upper()
        if upper.startswith("NODE_COORD_SECTION"):
            section = "coord"; continue
        if upper.startswith("DEMAND_SECTION"):
            section = "demand"; continue
        if upper.startswith("EDGE_WEIGHT_SECTION"):
            section = "weight"; continue
        if upper.startswith("DEPOT_SECTION"):
            section = "depot"; continue
        if upper.startswith("DISPLAY_DATA_SECTION"):
            section = "skip"; continue
        if ":" in line and section is None:
            key = line.split(":", 1)[0].strip().upper()
            header[key] = _parse_header_value(line)
            continue
        parts = line.split()
        if section == "coord":
            coords[int(parts[0])] = (float(parts[1]), float(parts[2]))
        elif section == "demand":
            demand[int(parts[0])] = int(parts[1])
        elif section == "weight":
            weights.extend(int(float(p)) for p in parts)
        elif section == "depot":
            if int(parts[0]) > 0:
                depot = int(parts[0])

    dim = int(header["DIMENSION"])
    capacity = int(header["CAPACITY"])
    if depot != 1:
        raise ValueError(f"{name}: depósito diferente do nó 1 não suportado")
    # nó TSPLIB k -> índice k-1 (depósito = 0)
    dem = [demand[k] for k in range(1, dim + 1)]
    edge_type = header.get("EDGE_WEIGHT_TYPE", "").upper()

    dist: list[list[int]] = []
    pts = None
    if edge_type == "EUC_2D":
        pts = [coords[k] for k in range(1, dim + 1)]
        if with_matrix:
            dist = [[nint(math.hypot(xi - xj, yi - yj)) for (xj, yj) in pts] for (xi, yi) in pts]
    elif edge_type == "EXPLICIT":
        fmt = header.get("EDGE_WEIGHT_FORMAT", "").upper()
        if fmt != "LOWER_ROW":
            raise ValueError(f"{name}: formato {fmt} não suportado")
        dist = [[0] * dim for _ in range(dim)]
        it = iter(weights)
        for i in range(1, dim):          # LOWER_ROW: linha i tem i elementos (sem diagonal)
            for j in range(i):
                dist[i][j] = dist[j][i] = next(it)
    else:
        raise ValueError(f"{name}: EDGE_WEIGHT_TYPE {edge_type} não suportado")

    inst = Instance(name=name, capacity=capacity, demand=dem, dist=dist, coords=pts)
    sol_path = raw_dir / f"{name}.sol"
    if sol_path.exists():
        inst.bks, inst.bks_routes = read_solution(sol_path)
    return inst


def read_solution(path: Path) -> tuple[int, list[list[int]]]:
    """Lê arquivo .sol do CVRPLIB (clientes numerados 1..n, depósito implícito)."""
    routes: list[list[int]] = []
    cost = None
    for line in path.read_text().splitlines():
        line = line.strip()
        if line.lower().startswith("route"):
            routes.append([int(v) for v in line.split(":", 1)[1].split()])
        elif line.lower().startswith("cost"):
            cost = int(float(line.split()[1]))
    return cost, routes


def routes_cost(inst: Instance, routes: list[list[int]]) -> int:
    total = 0
    for r in routes:
        path = [0] + r + [0]
        total += sum(inst.dist[a][b] for a, b in zip(path, path[1:]))
    return total


def check_routes(inst: Instance, routes: list[list[int]]) -> None:
    """Valida cobertura única e capacidade (levanta exceção se inválida)."""
    seen = sorted(c for r in routes for c in r)
    if seen != list(range(1, inst.n_customers + 1)):
        raise ValueError(f"{inst.name}: clientes ausentes ou repetidos")
    for r in routes:
        load = sum(inst.demand[c] for c in r)
        if load > inst.capacity:
            raise ValueError(f"{inst.name}: rota com carga {load} > {inst.capacity}")
