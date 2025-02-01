from qiskit.quantum_info import SparsePauliOp
import networkx as nx


def build_maxcut_hamiltonian(graph: nx.Graph) -> SparsePauliOp:
    """
    Build the MaxCut Hamiltonian for the given graph H = (|E|/2)*I - (1/2)*Σ_{(i,j)∈E}(Z_i Z_j)
    """
    num_qubits = len(graph.nodes)
    edges = list(graph.edges())
    num_edges = len(edges)

    pauli_terms = ["I"*num_qubits] # start with identity
    coeffs = [-num_edges / 2]

    for (u, v) in edges: # for each edge, add -(1/2)*Z_i Z_j
        z_term = ["I"] * num_qubits
        z_term[u] = "Z"
        z_term[v] = "Z"
        pauli_terms.append("".join(z_term))
        coeffs.append(0.5)

    return SparsePauliOp.from_list(list(zip(pauli_terms, coeffs)))

def build_2nd_hamiltonian(graph: nx.Graph, penalty_weight: float) -> SparsePauliOp:
    """
    Build a custom Hamiltonian for the given graph with the following cost function:
    H = (1/2) * Σ_{(i,j) ∈ E} (x_i - x_j)^2 - λ * (0.5 * num_values_in_E - Σ_i x_i)^2
    
    Args:
        graph: NetworkX graph representing the problem.
        num_values_in_E: Total number of values in the set E (used in the penalty term).
        penalty_weight: The penalty weight λ to control the strength of the second term.
    
    Returns:
        SparsePauliOp representing the Hamiltonian.
    """
    num_qubits = len(graph.nodes)
    edges = list(graph.edges())
    num_edges = len(edges)

    pauli_terms = []  # List to store the Pauli strings
    coeffs = []  # List to store the coefficients

    # First term: (1/2) * Σ_{(i,j) ∈ E} (x_i - x_j)^2
    for (u, v) in edges:
        z_term = ["I"] * num_qubits
        z_term[u] = "Z"
        z_term[v] = "Z"
        pauli_terms.append("".join(z_term))
        coeffs.append(0.5)

    # Second term: penalty term - λ * (0.5 * num_values_in_E - Σ_i x_i)^2
    # This can be written as a sum of X operators for each qubit
    penalty_pauli_term = ["I"] * num_qubits
    for i in range(num_qubits):
        penalty_pauli_term[i] = "X"  # Apply X to all qubits for the sum Σ_i x_i
    pauli_terms.append("".join(penalty_pauli_term))
    coeffs.append(-penalty_weight * (0.5 * num_edges) ** 2)

    # Create the Hamiltonian using SparsePauliOp
    return SparsePauliOp.from_list(list(zip(pauli_terms, coeffs)))