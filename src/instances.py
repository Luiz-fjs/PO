"""Conjunto de instâncias do experimento, agrupadas por faixa de tamanho.

BKS e status de otimalidade conferidos em
https://galgos.inf.puc-rio.br/cvrplib/en/cvrp/instances (consulta em 16/11/2026).
"""
from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class InstanceInfo:
    name: str
    tier: str
    bks_is_optimal: bool

    @property
    def fixed_fleet(self) -> int | None:
        # Conjuntos clássicos (E, A): frota fixa em K = sufixo "-kK" do nome.
        # Conjuntos X e AGS: frota livre (Uchoa et al., 2017).
        if self.name[0] in "EA":
            return int(self.name.rsplit("-k", 1)[1])
        return None


INSTANCES = [
    InstanceInfo("E-n13-k4", "Pequena", True),
    InstanceInfo("E-n22-k4", "Pequena", True),
    InstanceInfo("E-n23-k3", "Pequena", True),
    InstanceInfo("E-n30-k3", "Pequena", True),
    InstanceInfo("A-n32-k5", "Média", True),
    InstanceInfo("A-n45-k6", "Média", True),
    InstanceInfo("A-n60-k9", "Média", True),
    InstanceInfo("A-n80-k10", "Média", True),
    InstanceInfo("X-n101-k25", "Grande", True),
    InstanceInfo("X-n153-k22", "Grande", True),
    InstanceInfo("X-n200-k36", "Grande", True),
    InstanceInfo("X-n303-k21", "Muito grande", False),
    InstanceInfo("X-n502-k39", "Muito grande", False),
    InstanceInfo("X-n1001-k43", "Muito grande", False),
    InstanceInfo("Leuven1", "Massiva", False),
]

# Leuven1 (3000 clientes) gera ~9 milhões de variáveis binárias no modelo
# compacto; fica restrita ao Simulated Annealing.
EXACT_INSTANCES = [i for i in INSTANCES if i.tier != "Massiva"]
BY_NAME = {i.name: i for i in INSTANCES}
