import networkx as nx
import numpy as np
import groupEigenLs
import sys

"""
glee.py
-------
Geometric Laplacian Eigenmap Embedding.
"""


def laplacian_degrees(graph):
    """Return the Laplacian matrix and the diagonal degree matrix."""
    adj = nx.to_scipy_sparse_matrix(graph, format='csr')
    degs = sparse.diags(adj.sum(axis=1).flatten(), [0], adj.shape, format='csr')
    return degs - adj, degs


def eigenmaps(graph, dim, method='glee', return_vals=False):
    """Return the Eigenmap embedding of the graph.
    Params
    ------
    graph (nx.Graph): the graph to embed.
    method (str): 'glee' for GLEE (default), or 'eigen' for original
    Laplacian Eigenmaps.
    dim (int): number of embedding dimensions.
    return_vals (bool): whether to return the eigenvalues. Default False.
    """


    lapl, degs = laplacian_degrees(graph)
    shape = (graph.order(), graph.order())
    invdegs = sparse.diags(1 / degs.diagonal(), 0, shape, format='csr')

    # added by mkrahn to calc the least sig eigenvalue
    if int(sys.argv[-1]) == -1:
        eigenval, _ = np.linalg.eig(lapl.todense())
        eigenval = np.abs(eigenval)
        eigenval = np.sort(eigenval)
        #temp = list(map(np.sum,[eigenval[0:x] for x in range (1,len(eigenval))])) 
        eigenval = eigenval / eigenval[-1]
        dim = np.where(eigenval >= 0.7)[0][0]
        #dim = min(dim,12)
    
    eigenval, _ = np.linalg.eig(lapl.todense())
    eigenval = eigenval[:dim]
    ls = groupEigenLs.groupLs(eigenval)


    if dim > graph.order() - 1 or dim < 2:
        raise ValueError('dim must be grater than 0 and less than graph.order()')

    eig_fun = {
        'glee': {
            True: np.linalg.eigh,
            False: lambda m: sparse.linalg.eigsh(
                m, k=dim, which='LM', return_eigenvectors=True)},
        'eigen': {
            True: np.linalg.eig,
            # In the following we need to compute dim+1 eigenvectors so
            # that after deleting the zero eigenvalue we get exactly dim.
            False: lambda m: sparse.linalg.eigs(
                m, k=dim+1, which='SM', return_eigenvectors=True)}
    }
    matrix_fun = {
        'glee': {True: lambda: lapl.A, False: lambda: lapl},
        'eigen': {True: lambda: invdegs.A.dot(lapl.A),
                  False: lambda: invdegs.dot(lapl)}
    }
    post_process = {
        'glee': lambda vals, vecs: \
        (vals, vecs.dot(np.diag(np.sqrt(vals)))),
        'eigen': lambda vals, vecs: \
        (vals[1:], vecs[:, 1:].dot(
            np.diag([1/np.sqrt(v.dot(degs.A).dot(v))
                     for v in vecs[:, 1:].T])))
    }

    is_full = dim is None or dim > graph.order() - 3
    eig = eig_fun[method][is_full]
    matrix = matrix_fun[method][is_full]()
    vals, vecs = eig(matrix)
    # All eigenvalues are guaranteed to be non-negative, but sometimes the
    # zero eigenvalues can come very close to zero but negative, so we take
    # the absolute value as we need to take sqrt.
    vals, vecs = np.abs(vals.real), vecs.real
    vals, vecs = post_process[method](vals, vecs)

    indices = np.argsort(vals)
    indices = indices[::-1] if method == 'glee' else indices
    vals = vals[indices][:dim]
    vecs = (vecs.T[indices].T)[:, :dim]
    return (vecs, vals) if return_vals else (vecs, ls)


from scipy import sparse
import os

dirname = str(os.path.dirname(os.path.abspath(__file__)))

adja1 = np.loadtxt(os.path.join(dirname, 'adja1.txt'), delimiter = ',')
adja2 = np.loadtxt(os.path.join(dirname, 'adja2.txt'), delimiter = ',')

G1 = nx.convert_matrix.from_numpy_matrix(adja1)
G2 = nx.convert_matrix.from_numpy_matrix(adja2)

dimension = int(sys.argv[-1])

embeddedGraph1, groupLS1 = eigenmaps(G1, dimension, 'glee')
embeddedGraph2, groupLS2 = eigenmaps(G2, dimension, 'glee')


with open(os.path.join(dirname, 'embeddedG1.txt'), 'w') as f:
    np.savetxt(f, embeddedGraph1)

with open(os.path.join(dirname, 'groupG1.txt'), 'w') as f:
    np.savetxt(f, groupLS1)

with open(os.path.join(dirname, 'embeddedG2.txt'), 'w') as f:
    np.savetxt(f, embeddedGraph2)

with open(os.path.join(dirname, 'groupG2.txt'), 'w') as f:
    np.savetxt(f, groupLS2)


