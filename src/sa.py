"""Simulated Annealing para o CVRP (Etapa 6).

Implementação em Python puro, seguindo a seção 4 da proposta (Proposta.md):

* solução inicial: economias de Clarke & Wright (versão paralela);
* vizinhanças: relocate, swap, 2-opt (intra-rota) e 2-opt* (entre rotas),
  com o segundo cliente sorteado entre os k vizinhos mais próximos do primeiro
  (vizinhança granular, Toth & Vigo, 2003);
* função de avaliação F(s) = distância + lambda * excesso de carga, com lambda
  adaptativo (a busca pode atravessar soluções inviáveis; a melhor solução
  guardada é sempre viável);
* resfriamento geométrico, L = beta * n iterações por temperatura e
  reaquecimento a partir da melhor solução após R patamares sem melhora.

Todos os movimentos são avaliados de forma incremental: O(1) para os movimentos
entre rotas e O(tamanho da rota) para os intra-rota, de modo que o custo de uma
iteração não depende de n.

Uso isolado:
    python src/sa.py A-n32-k5 --seed 1 --time-limit 60
"""
from __future__ import annotations

import argparse
import heapq
import math
import random
import time
from dataclasses import asdict, dataclass

from cvrp_io import Instance, check_routes, read_instance, routes_cost


@dataclass
class SAParams:
    p0: float = 0.05                # prob. inicial de aceitar uma piora média -> T0
    tf_ratio: float = 1e-2          # temperatura final de cada ciclo: Tf = tf_ratio * T0
    n_cycles: int = 5               # ciclos de recozimento no orçamento de tempo
    reheat_factor: float = 0.1      # reaquecimento: T <- reheat_factor * T0
    beta: float = 1.0               # T é atualizada a cada max(100, beta*n) iterações
    granular_k: int = 20            # tamanho da lista de vizinhos
    time_limit: float | None = None  # None -> min(3600, max(60, n))
    seed: int = 1
    lam_every: int | None = None    # ajuste de lambda a cada X iterações (None -> max(100, n))
    lam_factor: float = 1.05
    target_feasible: float = 0.5
    new_route_prob: float = 0.02    # só com frota livre: abre rota nova
    polish_stall: int = 20          # polimento termina após polish_stall*n movimentos sem melhora


# operadores sorteados com probabilidades iguais:
# 0 relocate, 1 swap, 2 2-opt/2-opt*, 4 or-opt (3 = rota nova, só frota livre)
OPS = (0, 1, 2, 4)


def default_time_limit(n: int) -> float:
    return float(min(3600, max(60, n)))


# --------------------------------------------------------------------------
# Construção
# --------------------------------------------------------------------------
def granular_neighbors(d: list[list[int]], n: int, k: int) -> list[list[int]]:
    k = min(k, n - 1)
    neigh: list[list[int]] = [[]]
    custs = range(1, n + 1)
    for i in custs:
        row = d[i]
        neigh.append([j for j in heapq.nsmallest(k + 1, custs, key=row.__getitem__) if j != i][:k])
    return neigh


def clarke_wright(inst: Instance, neigh: list[list[int]] | None = None) -> list[list[int]]:
    """Economias de Clarke & Wright (1964), versão paralela.

    Com n grande, só os pares (i, j) da lista granular são considerados.
    """
    d, q, Q, n = inst.dist, inst.demand, inst.capacity, inst.n_customers
    d0 = d[0]
    if neigh is None:
        pairs = ((i, j) for i in range(1, n + 1) for j in range(i + 1, n + 1))
    else:
        pairs = {(min(i, j), max(i, j)) for i in range(1, n + 1) for j in neigh[i]}
    savings = sorted(((d0[i] + d0[j] - d[i][j], i, j) for i, j in pairs), reverse=True)
    route_of = list(range(n + 1))
    routes = {i: [i] for i in range(1, n + 1)}
    load = {i: q[i] for i in range(1, n + 1)}
    for s, i, j in savings:
        if s <= 0:
            break
        ri, rj = route_of[i], route_of[j]
        if ri == rj or load[ri] + load[rj] > Q:
            continue
        Ri, Rj = routes[ri], routes[rj]
        if Ri[-1] == i:
            A = Ri
        elif Ri[0] == i:
            A = Ri[::-1]
        else:
            continue
        if Rj[0] == j:
            B = Rj
        elif Rj[-1] == j:
            B = Rj[::-1]
        else:
            continue
        routes[ri] = A + B
        load[ri] += load[rj]
        for c in B:
            route_of[c] = ri
        del routes[rj], load[rj]
    return list(routes.values())


def fix_fleet(inst: Instance, routes: list[list[int]], K: int) -> list[list[int]]:
    """Ajusta a solução para exatamente K rotas (conjuntos E e A).

    Rotas excedentes (as de menor carga) são dissolvidas e seus clientes
    reinseridos na posição mais barata, preferindo posições que respeitem a
    capacidade; a capacidade eventualmente violada é tratada pela penalidade.
    Se houver menos de K rotas, a maior é dividida ao meio.
    """
    d, q, Q = inst.dist, inst.demand, inst.capacity
    routes = [r[:] for r in routes]
    while len(routes) > K:
        routes.sort(key=lambda r: sum(q[c] for c in r))
        orphan = routes.pop(0)
        for c in sorted(orphan, key=lambda c: -q[c]):
            best = None
            for r in routes:
                over = sum(q[x] for x in r) + q[c] > Q
                path = [0] + r + [0]
                for p in range(len(path) - 1):
                    a, b = path[p], path[p + 1]
                    key = (over, d[a][c] + d[c][b] - d[a][b])
                    if best is None or key < best[0]:
                        best = (key, r, p)
            _, r, p = best
            r.insert(p, c)
    while len(routes) < K:
        routes.sort(key=len)
        big = routes.pop()
        h = len(big) // 2
        routes += [big[:h], big[h:]]
    return routes


# --------------------------------------------------------------------------
# Simulated Annealing
# --------------------------------------------------------------------------
class SimulatedAnnealing:
    def __init__(self, inst: Instance, fixed_fleet: int | None, params: SAParams):
        self.inst = inst
        self.d, self.q, self.Q = inst.dist, inst.demand, inst.capacity
        self.n = inst.n_customers
        self.K = fixed_fleet
        self.P = params
        self.rng = random.Random(params.seed)
        self.neigh = granular_neighbors(self.d, self.n, params.granular_k)

    # ---------------- estado ----------------
    def _load_state(self, routes: list[list[int]]) -> None:
        n = self.n
        self.routes = [r[:] for r in routes if r]
        self.rid = [0] * (n + 1)
        self.pos = [0] * (n + 1)
        self.load = []
        for r, seq in enumerate(self.routes):
            for p, c in enumerate(seq):
                self.rid[c] = r
                self.pos[c] = p
            self.load.append(sum(self.q[c] for c in seq))
        Q = self.Q
        self.cost = routes_cost(self.inst, self.routes)
        self.excess = sum(l - Q for l in self.load if l > Q)

    def _rcost(self, seq: list[int]) -> int:
        if not seq:
            return 0
        d = self.d
        s = d[0][seq[0]] + d[seq[-1]][0]
        for a, b in zip(seq, seq[1:]):
            s += d[a][b]
        return s

    def _set_route(self, r: int, seq: list[int]) -> None:
        self.routes[r] = seq
        rid, pos, q, Q = self.rid, self.pos, self.q, self.Q
        new = 0
        for p, c in enumerate(seq):
            rid[c] = r
            pos[c] = p
            new += q[c]
        old = self.load[r]
        self.excess += (new - Q if new > Q else 0) - (old - Q if old > Q else 0)
        self.load[r] = new

    def _drop_empty(self) -> None:
        empties = [r for r, seq in enumerate(self.routes) if not seq]
        for r in sorted(empties, reverse=True):
            last = len(self.routes) - 1
            if r != last:
                self.routes[r] = self.routes[last]
                self.load[r] = self.load[last]
                for c in self.routes[r]:
                    self.rid[c] = r
            self.routes.pop()
            self.load.pop()

    # ---------------- movimentos ----------------
    # Cada avaliação devolve (delta_custo, delta_excesso, movimento) ou None.
    def _ex(self, l: int) -> int:
        return l - self.Q if l > self.Q else 0

    def eval_move(self, op: int, i: int, j: int):
        d, rid, pos, routes, load, q, ex = self.d, self.rid, self.pos, self.routes, self.load, self.q, self._ex
        r1, p1 = rid[i], pos[i]
        R1 = routes[r1]
        L1 = len(R1)
        a = R1[p1 - 1] if p1 else 0
        b = R1[p1 + 1] if p1 + 1 < L1 else 0

        if op == 3:  # nova rota com i sozinho (frota livre)
            if L1 == 1:
                return None
            dc = d[0][i] + d[i][0] - (d[a][i] + d[i][b] - d[a][b])
            dex = ex(load[r1] - q[i]) - ex(load[r1])
            return dc, dex, ("new", r1, p1)

        r2, p2 = rid[j], pos[j]
        R2 = routes[r2]
        L2 = len(R2)

        if op == 4:  # or-opt: bloco de 2-3 clientes a partir de i vai para perto de j
            sl = 2 if self.rng.random() < 0.5 else 3
            if p1 + sl > L1:
                return None
            seg = R1[p1:p1 + sl]
            if self.rng.random() < 0.5:
                seg = seg[::-1]
            if r1 == r2:
                if j in seg:
                    return None
                seq = R1[:p1] + R1[p1 + sl:]
                k = seq.index(j)
                k = k + 1 if self.rng.random() < 0.5 else k
                new = seq[:k] + seg + seq[k:]
                if new == R1:
                    return None
                return self._rcost(new) - self._rcost(R1), 0, ("set", r1, new)
            if sl == L1 and self.K is not None:
                return None
            b2 = R1[p1 + sl] if p1 + sl < L1 else 0
            if self.rng.random() < 0.5:
                c, e, ins = j, (R2[p2 + 1] if p2 + 1 < L2 else 0), p2 + 1
            else:
                c, e, ins = (R2[p2 - 1] if p2 else 0), j, p2
            f, l = seg[0], seg[-1]
            dc = (d[c][f] + d[l][e] - d[c][e]
                  - (d[a][R1[p1]] + d[R1[p1 + sl - 1]][b2] - d[a][b2]))
            sq = 0
            for x in seg:
                sq += q[x]
            l1, l2 = load[r1], load[r2]
            dex = ex(l1 - sq) + ex(l2 + sq) - ex(l1) - ex(l2)
            return dc, dex, ("oropt", r1, p1, r2, ins, seg)

        if r1 == r2:  # ---------- movimentos intra-rota ----------
            if op == 0:    # relocate: i para logo depois/antes de j
                seq = R1[:p1] + R1[p1 + 1:]
                k = seq.index(j)
                k = k + 1 if self.rng.random() < 0.5 else k
                new = seq[:k] + [i] + seq[k:]
            elif op == 1:  # swap
                new = R1[:]
                new[p1], new[p2] = j, i
            else:          # 2-opt: torna i e j adjacentes invertendo um trecho
                lo, hi = (p1, p2) if p1 < p2 else (p2, p1)
                if hi == lo + 1:
                    return None
                u, v = R1[lo], R1[hi]
                if self.rng.random() < 0.5:   # (u,su),(v,sv) -> (u,v),(su,sv)
                    su = R1[lo + 1]
                    sv = R1[hi + 1] if hi + 1 < L1 else 0
                    dc = d[u][v] + d[su][sv] - d[u][su] - d[v][sv]
                    return dc, 0, ("rev", r1, lo + 1, hi)
                pu = R1[lo - 1] if lo else 0      # (pu,u),(pv,v) -> (pu,pv),(u,v)
                pv = R1[hi - 1]
                dc = d[pu][pv] + d[u][v] - d[pu][u] - d[pv][v]
                return dc, 0, ("rev", r1, lo, hi - 1)
            if new == R1:
                return None
            return self._rcost(new) - self._rcost(R1), 0, ("set", r1, new)

        # ---------- movimentos entre rotas ----------
        l1, l2 = load[r1], load[r2]
        if op == 0:  # relocate i para a rota de j
            if L1 == 1 and self.K is not None:
                return None
            if self.rng.random() < 0.5:
                c, e, ins = j, (R2[p2 + 1] if p2 + 1 < L2 else 0), p2 + 1
            else:
                c, e, ins = (R2[p2 - 1] if p2 else 0), j, p2
            dc = d[c][i] + d[i][e] - d[c][e] - (d[a][i] + d[i][b] - d[a][b])
            qi = q[i]
            dex = ex(l1 - qi) + ex(l2 + qi) - ex(l1) - ex(l2)
            return dc, dex, ("reloc", r1, p1, r2, ins)
        if op == 1:  # swap entre rotas
            c = R2[p2 - 1] if p2 else 0
            e = R2[p2 + 1] if p2 + 1 < L2 else 0
            dc = (d[a][j] + d[j][b] - d[a][i] - d[i][b]
                  + d[c][i] + d[i][e] - d[c][j] - d[j][e])
            t = q[j] - q[i]
            dex = ex(l1 + t) + ex(l2 - t) - ex(l1) - ex(l2)
            return dc, dex, ("swap", r1, p1, r2, p2)
        # 2-opt*: troca de caudas entre as rotas de i e j
        sj = R2[p2 + 1] if p2 + 1 < L2 else 0
        pre1 = 0
        for c in R1[:p1 + 1]:
            pre1 += q[c]
        pre2 = 0
        for c in R2[:p2 + 1]:
            pre2 += q[c]
        if self.rng.random() < 0.5:
            # A: R1[:p1+1] + R2[p2+1:] | R2[:p2+1] + R1[p1+1:]
            dc = d[i][sj] + d[j][b] - d[i][b] - d[j][sj]
            n1, n2 = pre1 + l2 - pre2, pre2 + l1 - pre1
            kind = "2optA"
        else:
            # B: R1[:p1+1] + rev(R2[:p2+1]) | rev(R1[p1+1:]) + R2[p2+1:]
            if b == 0 and sj == 0 and self.K is not None:
                return None
            dc = d[i][j] + d[b][sj] - d[i][b] - d[j][sj]
            n1 = pre1 + pre2
            n2 = l1 + l2 - n1
            kind = "2optB"
        dex = ex(n1) + ex(n2) - ex(l1) - ex(l2)
        return dc, dex, (kind, r1, p1, r2, p2)

    def apply(self, mv, dc: int) -> None:
        kind = mv[0]
        routes = self.routes
        if kind == "set":
            self._set_route(mv[1], mv[2])
        elif kind == "rev":
            _, r, lo, hi = mv
            R = routes[r]
            self._set_route(r, R[:lo] + R[lo:hi + 1][::-1] + R[hi + 1:])
        elif kind == "new":
            _, r1, p1 = mv
            R1 = routes[r1]
            i = R1[p1]
            self._set_route(r1, R1[:p1] + R1[p1 + 1:])
            routes.append([])
            self.load.append(0)
            self._set_route(len(routes) - 1, [i])
        elif kind == "oropt":
            _, r1, p1, r2, ins, seg = mv
            R1, R2 = routes[r1], routes[r2]
            self._set_route(r1, R1[:p1] + R1[p1 + len(seg):])
            self._set_route(r2, R2[:ins] + seg + R2[ins:])
            if not routes[r1]:
                self._drop_empty()
        else:
            _, r1, p1, r2, p2 = mv
            R1, R2 = routes[r1], routes[r2]
            if kind == "reloc":
                i = R1[p1]
                self._set_route(r1, R1[:p1] + R1[p1 + 1:])
                self._set_route(r2, R2[:p2] + [i] + R2[p2:])
            elif kind == "swap":
                n1, n2 = R1[:], R2[:]
                n1[p1], n2[p2] = R2[p2], R1[p1]
                self._set_route(r1, n1)
                self._set_route(r2, n2)
            elif kind == "2optA":
                self._set_route(r1, R1[:p1 + 1] + R2[p2 + 1:])
                self._set_route(r2, R2[:p2 + 1] + R1[p1 + 1:])
            else:  # 2optB
                self._set_route(r1, R1[:p1 + 1] + R2[:p2 + 1][::-1])
                self._set_route(r2, R1[p1 + 1:][::-1] + R2[p2 + 1:])
            if not routes[r1] or not routes[r2]:
                self._drop_empty()
        self.cost += dc

    def _random_move(self):
        rng = self.rng
        i = rng.randint(1, self.n)
        if self.K is None and rng.random() < self.P.new_route_prob:
            return self.eval_move(3, i, 0)
        j = rng.choice(self.neigh[i])
        return self.eval_move(OPS[rng.randrange(4)], i, j)

    def descent(self, lam: float, stall: int, t_end: float) -> None:
        """Reparo + polimento: descida estrita (T = 0) com penalidade alta.

        Aceita só movimentos que reduzem F; lambda é multiplicado por 10 a cada
        rodada em que a solução continua inviável. Termina após `stall`
        movimentos seguidos sem melhora numa solução viável.
        """
        rng, n, neigh = self.rng, self.n, self.neigh
        randint, choice, randrange = rng.randint, rng.choice, rng.randrange
        eval_move, apply_move = self.eval_move, self.apply
        perf = time.perf_counter
        k = 0
        while True:
            no_imp = 0
            while no_imp < stall:
                k += 1
                if k % 5000 == 0 and perf() >= t_end:
                    return
                i = randint(1, n)
                m = eval_move(OPS[randrange(4)], i, choice(neigh[i]))
                if m is not None and m[0] + lam * m[1] < 0:
                    apply_move(m[2], m[0])
                    no_imp = 0
                else:
                    no_imp += 1
            if self.excess == 0 or time.perf_counter() >= t_end:
                return
            lam *= 10

    # ---------------- laço principal ----------------
    def solve(self, initial: list[list[int]] | None = None, log=None) -> dict:
        P, rng, n = self.P, self.rng, self.n
        t_start = time.perf_counter()
        t_max = P.time_limit if P.time_limit is not None else default_time_limit(n)

        if initial is None:
            cw_neigh = None if n <= 1000 else granular_neighbors(self.d, n, max(50, P.granular_k))
            initial = clarke_wright(self.inst, cw_neigh)
            if self.K is not None:
                initial = fix_fleet(self.inst, initial, self.K)
        self._load_state(initial)
        initial_cost = self.cost
        t_init = time.perf_counter() - t_start

        # lambda inicial: distância média até os vizinhos granulares
        d = self.d
        lam = sum(d[i][j] for i in range(1, n + 1) for j in self.neigh[i]) / sum(
            len(self.neigh[i]) for i in range(1, n + 1))
        lam = max(lam, 1.0)
        lam0 = lam

        # T0 = -media(pioras) / ln(p0), estimada com 1000 movimentos aleatórios
        worse = []
        for _ in range(1000):
            m = self._random_move()
            if m is not None and m[0] > 0:   # só a variação de distância
                worse.append(m[0])
        T0 = (-(sum(worse) / len(worse)) / math.log(P.p0)) if worse else 1.0
        Tf = P.tf_ratio * T0
        L = max(100, int(P.beta * n))       # iterações entre atualizações de T
        # tempo restante após a construção, com folga para o polimento do último ciclo
        cycle_dur = 0.98 * max(1e-3, t_max - (time.perf_counter() - t_start)) / P.n_cycles
        cycle = 0
        cycle_t0 = time.perf_counter()
        Ts = T0                               # temperatura no início do ciclo
        T = T0
        log_ratio = math.log(Tf / Ts)

        best_routes = [r[:] for r in self.routes] if self.excess == 0 else None
        best_cost = self.cost if self.excess == 0 else math.inf
        t_best = time.perf_counter() - t_start
        trace = [(round(t_best, 3), best_cost)] if best_routes else []

        it = 0
        feas_count = 0
        exp, rand = math.exp, rng.random
        lam_every = P.lam_every or max(100, n)
        eval_move, apply_move = self.eval_move, self.apply
        randint, choice, randrange = rng.randint, rng.choice, rng.randrange
        neigh, free, pnew = self.neigh, self.K is None, P.new_route_prob
        perf = time.perf_counter
        while True:
            for _ in range(L):
                it += 1
                i = randint(1, n)
                if free and rand() < pnew:
                    m = eval_move(3, i, 0)
                else:
                    m = eval_move(OPS[randrange(4)], i, choice(neigh[i]))
                if m is not None:
                    dc, dex, mv = m
                    dF = dc + lam * dex
                    if dF <= 0 or rand() < exp(-dF / T):
                        apply_move(mv, dc)
                        if self.excess == 0 and self.cost < best_cost:
                            best_cost = self.cost
                            best_routes = [r[:] for r in self.routes]
                            t_best = perf() - t_start
                            trace.append((round(t_best, 3), best_cost))
                if self.excess == 0:
                    feas_count += 1
                if it % lam_every == 0:
                    if feas_count < P.target_feasible * lam_every:
                        lam *= P.lam_factor
                    else:
                        lam /= P.lam_factor
                    feas_count = 0
            now = perf()
            if now - t_start >= t_max:
                break
            # resfriamento geométrico parametrizado pelo tempo:
            # T(tau) = Ts * (Tf/Ts)^tau, tau = fração decorrida do ciclo
            tau = (now - cycle_t0) / cycle_dur
            if tau < 1.0:
                T = Ts * exp(log_ratio * tau)
                continue
            # fim do ciclo: reparo/polimento da solução corrente e
            # reaquecimento a partir da melhor solução viável
            # polimento limitado a 5% do ciclo (e ao tempo total)
            self.descent(max(lam, lam0) * 10, P.polish_stall * n,
                         min(t_start + t_max, perf() + 0.05 * cycle_dur))
            if self.excess == 0 and self.cost < best_cost:
                best_cost = self.cost
                best_routes = [r[:] for r in self.routes]
                t_best = perf() - t_start
                trace.append((round(t_best, 3), best_cost))
            now = perf()
            cycle += 1
            if log:
                log(f"  ciclo {cycle}: it={it} best={best_cost} atual={self.cost} exc={self.excess} "
                    f"T={T:.3g} lam={lam:.3g} t={now - t_start:.1f}s")
            if best_routes is not None:
                self._load_state(best_routes)
            Ts = P.reheat_factor * T0
            T = Ts
            log_ratio = math.log(Tf / Ts)
            cycle_t0 = now

        total = time.perf_counter() - t_start
        if best_routes is None:
            raise RuntimeError(f"{self.inst.name}: SA não encontrou solução viável")
        check_routes(self.inst, best_routes)
        assert routes_cost(self.inst, best_routes) == best_cost
        return dict(cost=int(best_cost), routes=best_routes, n_routes=len(best_routes),
                    initial_cost=int(initial_cost), time_to_best=t_best, time=total,
                    time_init=t_init, iterations=it, cycles=cycle + 1, T0=T0, Tf=Tf,
                    trace=trace, params=asdict(P))


def solve(inst: Instance, fixed_fleet: int | None, params: SAParams, log=None) -> dict:
    return SimulatedAnnealing(inst, fixed_fleet, params).solve(log=log)


if __name__ == "__main__":
    from instances import BY_NAME

    ap = argparse.ArgumentParser()
    ap.add_argument("instance")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--time-limit", type=float, default=None)
    ap.add_argument("--p0", type=float, default=0.05)
    ap.add_argument("--tf-ratio", type=float, default=1e-2)
    ap.add_argument("--cycles", type=int, default=5)
    ap.add_argument("--reheat-factor", type=float, default=0.1)
    a = ap.parse_args()
    inst = read_instance(a.instance)
    info = BY_NAME.get(a.instance)
    K = info.fixed_fleet if info else (int(a.instance.rsplit("-k", 1)[1]) if a.instance[0] in "EA" else None)
    params = SAParams(seed=a.seed, time_limit=a.time_limit, p0=a.p0, tf_ratio=a.tf_ratio,
                      n_cycles=a.cycles, reheat_factor=a.reheat_factor)
    res = solve(inst, K, params, log=print)
    gap = f"{100 * (res['cost'] - inst.bks) / inst.bks:.2f}%" if inst.bks else "-"
    print(f"{a.instance}: custo={res['cost']} (BKS {inst.bks}, desvio {gap}) rotas={res['n_routes']} "
          f"inicial={res['initial_cost']} t_best={res['time_to_best']:.1f}s t={res['time']:.1f}s "
          f"it={res['iterations']} ({res['iterations'] / res['time']:.0f} it/s)")
