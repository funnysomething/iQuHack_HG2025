# Adi's algorithm
theta = ParameterVector(r"$\theta$", graph.number_of_edges())
    for t, (u, v) in zip(theta, graph.edges):
        ansatz.ry(t, v)
        ansatz.cx(u, v)

    beta = ParameterVector('b', graph.number_of_nodes())
    for b, n in zip(beta, graph.nodes):
        ansatz.ry(b, n)

# Brandon's
theta = ParameterVector(r"$\theta$", graph.number_of_edges())
    for t, (u, v) in zip(theta, graph.edges):
        ansatz.ry(t, v)
        ansatz.cx(u, v)

    beta = ParameterVector('b', graph.number_of_nodes())
    for b, n in zip(beta, graph.nodes):
        ansatz.rz(b, n)     
# Change ry to rz ^
