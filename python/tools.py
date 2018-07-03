import networkx as nx
import numpy as np
from scipy.optimize import brentq

def build_graph_from_adjacency(inadjacency):
    """
    Takes an adjacency_list like: "23 41 18" or 18 times  "23 41 1"   (edge from 23 --> 41)
    possibly having multiple edges and build a graph with no multiple edges but weigths representing how many of them there are
    Necessary in case of using algorithms that do not accept MultiGraphs. E.g. eigenvector centrality.    
    """
    
    adjacency_list=open(inadjacency,'r')
    edges={}
    for row in adjacency_list:
        a=row.split()
        e=(a[0],a[1])
        w=int(a[2])
        if(e not in edges):edges[e]=w
        else:edges[e]+=w
    G=nx.DiGraph()
    for e in edges: G.add_edge(e[0],e[1],weight=edges[e])
    adjacency_list.close()

    return G


def shift_rank(ranks):
    '''
    Shifts all scores so that the minimum is in zero and the others are all positive
    '''
    min_r=min(ranks)
    N=len(ranks)
    for i in range(N): ranks[i]=ranks[i]-min_r
    return ranks    

def btl(A,tol):
    N = np.shape(A)[0]
    g = np.random.rand(N)
    wins = np.array(np.sum(A,axis=1)).flatten();
    matches = np.array(A+np.transpose(A));
    totalMatches = np.array(np.sum(matches,axis=0)).flatten()
    g_prev = np.random.rand(N)
    eps = 1e-6
    while np.linalg.norm(g-g_prev) > tol:
        g_prev = g
        for i in range(N):
            if totalMatches[i]>0:
                q = np.divide(matches[i,:],g_prev[i]+g_prev)
                q[i] = 0
                g[i] = (wins[i]+eps)/np.sum(q)
            else:
                g[i] = 0
        g = g/np.sum(g)
    return np.log(g)

def eqs39(beta,s,A):
    N = np.shape(A)[0]
    x = 0
    for i in range(N):
        for j in range(N):
            if A[i,j] == 0:
                continue
            else:
                x += (s[i]-s[j]) * ( A[i,j] - (A[i,j]+A[j,i]) / (1+np.exp(-2*beta*(s[i]-s[j]))) )
    return x

def get_optimal_temperature(ranks,A):    
    return brentq(eqs39,0.01,20,args=(ranks,A))