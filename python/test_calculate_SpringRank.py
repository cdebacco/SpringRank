'''
Example of SPringRank usage

From a given network, it extracts the SpringRank scores

'''

import networkx as nx
import numpy as np
import SpringRank_tools as sr
import tools as tl

network='US_CS'

alpha=0.
l0=1.
l1=1.

'''
Builds graph and Adjacency matrix from input network
'''
inadjacency='../data/'+network+'_adjacency.dat'

G=tl.build_graph_from_adjacency(inadjacency)

nodes=list(G.nodes())			#  determines the order of the entries of matrix A

A=nx.to_numpy_matrix(G,nodelist=nodes)


'''
Extracts SpringRank
'''
rank=sr.SpringRank(A,alpha=alpha,l0=l0,l1=l1)

rank=tl.shift_rank(rank)   # (optional) shifts so that the min is in zero and the others are positive

'''
Order results so that the first node is the highest-ranked one
'''
X=[(nodes[i],rank[i]) for i in range(G.number_of_nodes())]
X= sorted(X, key=lambda tup: tup[1],reverse=True)
'''
Prints results
'''
print('SpringRank scores:')
outfile='../data/'+network+'_SpringRank_'+'a'+str(alpha)+'_l0_'+str(l0)+'_l1_'+str(l1)+'.dat'
outf=open(outfile,'w')

for i in range(G.number_of_nodes()):
	outf.write("{} {}\n".format(X[i][0],X[i][1]))
	# print nodes[i],rank[i]
	print(X[i][0],X[i][1])
print('Results saved in:', outfile)
outf.close()