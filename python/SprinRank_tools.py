import networkx as nx
import numpy as np

def SpringRank_planted_network(N,beta,alpha,K,prng,l0=0.5,l1=1.):
    '''

    Uses the SpringRank generative model to build a directed, possibly weigthed and having self-loops, network.
    Can be used to generate benchmarks for hierarchical networks

    Steps:
        1. Generates the scores (default is factorized Gaussian)
        2. Extracts A_ij entries (network edges) from Poisson distribution with average related to SpringRank energy

    INPUT:

        N=# of nodes
        beta= inverse temperature, controls noise
        alpha=controls prior's variance
        K=E/N  --> average degree, controls sparsity
        l0=prior spring's rest length 
        l1=interaction spring's rest lenght

    OUTPUT:
        G: nx.DiGraph()         Directed (possibly weighted graph, there can be self-loops)
        
    '''
    G=nx.DiGraph()

    scores=prng.normal(l0,1./np.sqrt(alpha*beta),N)  # planted scores ---> uses factorized Gaussian
    for i in range(N):G.add_node(i,score=scores[i])

    #  ---- Fixing sparsity i.e. the average degree  ---- 
    Z=0.
    for i in range(N):
        for j in range(N):  
            Z+=np.exp(-0.5*beta*np.power(scores[i]-scores[j]-l1,2))
    c=float(K*N)/Z        
    #  --------------------------------------------------

    # ----  Building the graph   ------------------------ 
    for i in range(N):
        for j in range(N):

            H_ij=0.5*np.power((scores[i]-scores[j]-l1),2)
            lambda_ij=c*np.exp(-beta*H_ij)

            A_ij=prng.poisson(lambda_ij,1)[0]

            if A_ij>0:G.add_edge(i,j,weight=A_ij)

    return G     
