from qiskit.circuit import ParameterVector
import networkx as nx
from qiskit import QuantumCircuit

gates = {
    0: 'CX_UV', 1: 'CY_UV', 2: 'CZ_UV', 
    3: 'CX_VU', 4: 'CY_VU', 5: 'CZ_VU',
    6: 'RX_U', 7: 'RY_U', 8: 'RZ_U',
    9: 'RX_V', 10: 'RY_V', 11: 'RZ_V'
}

# Deprecated TODO: CHECK IF THIS IS STILL NEEDED
# CX = 0
# CY = 1
# CZ = 2
# RX = 3
# RY = 4
# RZ = 5

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
                ansatz.cx(v, u)
            elif gate == 4:
                ansatz.cy(v, u)
            elif gate == 5:
                ansatz.cz(v, u)
            elif gate == 6:
                ansatz.rx(t, u)
            elif gate == 7:
                ansatz.ry(t, u)
            elif gate == 8:
                ansatz.rz(t, u)
            elif gate == 9:
                ansatz.rx(t, v)
            elif gate == 10:
                ansatz.ry(t, v)
            elif gate == 11:
                ansatz.rz(t, v)            
    
    return ansatz