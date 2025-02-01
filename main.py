import matplotlib.pyplot as plt
from IPython import display

import networkx as nx
import numpy as np
import pandas as pd
import time

from typing import List
from qiskit import QuantumCircuit, transpile
from qiskit.circuit import ParameterVector
from qiskit.quantum_info import SparsePauliOp
from qiskit_aer import AerSimulator

from itertools import product

from graph import graphs
from ansatz import (
    build_ansatz,
    build_with_gates,
    gates
)
from qitevolver import QITEvolver
from hamiltonian import build_maxcut_hamiltonian

class MaxCutSolver:

    backend = AerSimulator()
    shots = 100_000
    challenges = ['base', 'balanced', 'connected']

    def solve(self):
    # Loop through each type of ansatz, and try it on every graph,
    # summing up the scores. Then output it into a csv (?)

        for length in range(1, 2):
            for gate_perm in product(gates, repeat=length):
                # Perm is a permuatation of gates to build the ansatz
                for graph in graphs:
                    # Now we have a graph to test the ansatz on
                    ansatz = build_with_gates(graph, gate_perm) # Builds ansatz with given gate permutation
                    ham = build_maxcut_hamiltonian(graph)

                    qite_evolver = QITEvolver(ham, ansatz)
                    qite_evolver.evolve(num_steps=40, lr = 0.1, verbose = True)

                    qite_evolver.plot_convergence() # Comment out to remove graphs

                    # Run on the backend after optimizing ansatz
                    optimized_state = ansatz.assign_parameters(qite_evolver.param_vals[-1])
                    optimized_state.measure_all()
                    counts = self.backend.run(optimized_state, shots = self.shots)

                    # Find the sampled bitstring with the largest cut value
                    cut_vals = sorted(((bs, self.compute_cut_size(graph, bs)) for bs in counts), key=lambda t: t[1])
                    best_bs = cut_vals[-1][0]

                    # most_likely_soln = ""   # TODO: CHANGE TO MOST LIKELY

                    (best_brute, best_bal, best_con), (XS_brute, XS_balanced, XS_connected) = self.calculate_best_score(self, graph)

                    sum_counts = 0
                    for bs in counts:
                        if bs in XS_brute:
                            sum_counts += counts[bs]

                    sum_balanced_counts = 0
                    for bs in counts:
                        if bs in XS_balanced:
                            sum_balanced_counts += counts[bs]

                    sum_connected_counts = 0
                    for bs in counts:
                        if bs in XS_connected:
                            sum_connected_counts += counts[bs]

                    final_score = self.final_score(graph, XS_brute, XS_balanced, XS_connected, counts, self.shots, ansatz, 'base')
                    print("Final Score: " + final_score + "for test ", gate_perm)

    def final_score(self, graph, XS_brut, XS_balanced, XS_connected, counts,shots,ansatz,challenge):

        if(challenge=='base'):
            sum_counts = 0
            for bs in counts:
                if bs in XS_brut:
                    sum_counts += counts[bs]
        elif(challenge=='balanced'):
            sum_balanced_counts = 0
            for bs in counts:
                if bs in XS_balanced:
                    sum_balanced_counts += counts[bs]
            sum_counts = sum_balanced_counts
        elif(challenge=='connected'):
            sum_connected_counts = 0
            for bs in counts:
                if bs in XS_connected:
                    sum_connected_counts += counts[bs]
            sum_counts = sum_connected_counts


        transpiled_ansatz = transpile(ansatz, basis_gates = ['cx','rz','sx','x'])
        cx_count = transpiled_ansatz.count_ops()['cx']
        score = (4*2*graph.number_of_edges())/(4*2*graph.number_of_edges() + cx_count) * sum_counts/shots

        return np.round(score,5)



    def calculate_best_score(self, graph) -> tuple[int, int, int]:
        verbose = False

        G = graph
        n = len(G.nodes())
        w = np.zeros([n, n])
        for i in range(n):
            for j in range(n):
                temp = G.get_edge_data(i, j, default=0)
                if temp != 0:
                    w[i, j] = 1.0
        if verbose:
            print(w)

        best_cost_brute = 0
        best_cost_balanced = 0
        best_cost_connected = 0

        for b in range(2**n):
            x = [int(t) for t in reversed(list(bin(b)[2:].zfill(n)))]

            # Create subgraphs based on the partition
            subgraph0 = G.subgraph([i for i, val in enumerate(x) if val == 0])
            subgraph1 = G.subgraph([i for i, val in enumerate(x) if val == 1])

            bs = "".join(str(i) for i in x)

            # Check if subgraphs are not empty
            if len(subgraph0.nodes) > 0 and len(subgraph1.nodes) > 0:
                cost = 0
                for i in range(n):
                    for j in range(n):
                        cost = cost + w[i, j] * x[i] * (1 - x[j])
                if best_cost_brute < cost:
                    best_cost_brute = cost
                    xbest_brute = x
                    XS_brut = []
                if best_cost_brute == cost:
                    XS_brut.append(bs)

                outstr = "case = " + str(x) + " cost = " + str(cost)

                if (len(subgraph1.nodes)-len(subgraph0.nodes))**2 <= 1:
                    outstr += " balanced"
                    if best_cost_balanced < cost:
                        best_cost_balanced = cost
                        xbest_balanced = x
                        XS_balanced = []
                    if best_cost_balanced == cost:
                        XS_balanced.append(bs)

                if nx.is_connected(subgraph0) and nx.is_connected(subgraph1):
                    outstr += " connected"
                    if best_cost_connected < cost:
                        best_cost_connected = cost
                        xbest_connected = x
                        XS_connected = []
                    if best_cost_connected == cost:
                        XS_connected.append(bs)
                if verbose:
                    print(outstr)

        return (best_cost_brute, best_cost_balanced, best_cost_connected), (XS_brut, XS_balanced, XS_connected)
        
    def compute_cut_size(self, graph, bitstring):
        """
        Get the cut size of the partition of ``graph`` described by the given
        ``bitstring``.
        """
        cut_sz = 0
        for (u, v) in graph.edges:
            if bitstring[u] != bitstring[v]:
                cut_sz += 1
        return cut_sz

if __name__ == '__main__':
    solver = MaxCutSolver()
    solver.solve()
