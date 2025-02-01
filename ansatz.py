from qiskit.circuit import ParameterVector
import networkx as nx
from qiskit import QuantumCircuit

gates = {0: 'CX', 1: 'CY', 2: 'CZ', 3: 'RX', 4: 'RY', 5: 'RZ'}

CX = 0
CY = 1
CZ = 2
RX = 3
RY = 4
RZ = 5

def build_ansatz(graph: nx.Graph) -> QuantumCircuit:

    ansatz = QuantumCircuit(graph.number_of_nodes())
    ansatz.h(range(graph.number_of_nodes()))    # Apply hadamard to all qubits (putting in superposition)

    theta = ParameterVector(r"$\theta$", graph.number_of_edges())
    for t, (u, v) in zip(theta, graph.edges):
        ansatz.cz(u, v)
        ansatz.rz(t, v)
        ansatz.cy(u, v)

    return ansatz

def build_with_gates(graph: nx.Graph, gates: list[int]):
    ansatz = QuantumCircuit(graph.number_of_nodes())
    ansatz.h(range(graph.number_of_nodes()))

    theta = ParameterVector(r"$\theta$", graph.number_of_nodes())
    for t, (u, v) in zip(theta, graph.edges):
        for gate in gates:
            if gate == 0:
                ansatz.cx(u, v)
            elif gate == 1:
                ansatz.cy(u, v)
            elif gate == 2:
                ansatz.cz(u, v)
            elif gate == 3:
                ansatz.rx(t, v)
            elif gate == 4:
                ansatz.ry(t, v)
            elif gate == 5:
                ansatz.rz(t, v)
    
    return ansatz