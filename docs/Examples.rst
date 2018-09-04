.. _examples:

########
Examples
########

Load a graph from a list of edges
---------------------------------

In many cases, network structure is described by a list of edges.
In this example, we consider a standard file format composed of three columns, where
the first and second columns are the pairs of adjacent nodes, and third column is the weight of the edge between them, e.g.,  
 
.. code-block:: bash

   source,target,weight
   ACE,@BD,1
   BDF,@BD,1
   CEG,@BD,1
   DFH,@BD,1
   EGI,@BD,1
   FHJ,@BD,1
   GIK,@BD,1
   HJL,@BD,1
   IKM,@BD,1
   @BD,ACE,1
   BDF,ACE,1
   CEG,ACE,1
   DFH,ACE,1
   EGI,ACE,1
   FHJ,ACE,1
   GIK,ACE,1
   HJL,ACE,1
   IKM,ACE,1
   @BD,BDF,1
   ACE,BDF,1
   CEG,BDF,1
   DFH,BDF,1
   EGI,BDF,1
   FHJ,BDF,1
   GIK,BDF,1
   HJL,BDF,1
   IKM,BDF,1
   @BD,CEG,1
   ACE,CEG,1
   BDF,CEG,1
   DFH,CEG,1
   EGI,CEG,1
   FHJ,CEG,1
   GIK,CEG,1
   HJL,CEG,1
   IKM,CEG,1
   @BD,DFH,1
   ACE,DFH,1
   BDF,DFH,1
   CEG,DFH,1
   EGI,DFH,1
   FHJ,DFH,1
   GIK,DFH,1
   HJL,DFH,1
   IKM,DFH,1
   @BD,EGI,1
   ACE,EGI,1
   BDF,EGI,1
   CEG,EGI,1
   DFH,EGI,1
   @BD,FHJ,1
   ACE,FHJ,1
   BDF,FHJ,1
   CEG,FHJ,1
   DFH,FHJ,1
   @BD,GIK,1
   ACE,GIK,1
   BDF,GIK,1
   CEG,GIK,1
   DFH,GIK,1
   @BD,HJL,1
   ACE,HJL,1
   BDF,HJL,1
   CEG,HJL,1
   DFH,HJL,1
   @BD,IKM,1
   ACE,IKM,1
   BDF,IKM,1
   CEG,IKM,1
   DFH,IKM,1
   KMO,JLN,1
   LNP,JLN,1
   MOQ,JLN,1
   NPR,JLN,1
   OQS,JLN,1
   PRT,JLN,1
   QSU,JLN,1
   RTV,JLN,1
   SUW,JLN,1
   JLN,KMO,1
   LNP,KMO,1
   MOQ,KMO,1
   NPR,KMO,1
   OQS,KMO,1
   PRT,KMO,1
   QSU,KMO,1
   RTV,KMO,1
   SUW,KMO,1
   JLN,LNP,1
   KMO,LNP,1
   MOQ,LNP,1
   NPR,LNP,1
   OQS,LNP,1
   PRT,LNP,1
   QSU,LNP,1
   RTV,LNP,1
   SUW,LNP,1
   JLN,MOQ,1
   KMO,MOQ,1
   LNP,MOQ,1
   NPR,MOQ,1
   OQS,MOQ,1
   PRT,MOQ,1
   QSU,MOQ,1
   RTV,MOQ,1
   SUW,MOQ,1
   JLN,NPR,1
   KMO,NPR,1
   LNP,NPR,1
   MOQ,NPR,1
   OQS,NPR,1
   PRT,NPR,1
   QSU,NPR,1
   RTV,NPR,1
   SUW,NPR,1
   JLN,OQS,1
   KMO,OQS,1
   LNP,OQS,1
   MOQ,OQS,1
   NPR,OQS,1
   JLN,PRT,1
   KMO,PRT,1
   LNP,PRT,1
   MOQ,PRT,1
   NPR,PRT,1
   JLN,QSU,1
   KMO,QSU,1
   LNP,QSU,1
   MOQ,QSU,1
   NPR,QSU,1
   JLN,RTV,1
   KMO,RTV,1
   LNP,RTV,1
   MOQ,RTV,1
   NPR,RTV,1
   JLN,SUW,1
   KMO,SUW,1
   LNP,SUW,1
   MOQ,SUW,1
   NPR,SUW,1

Save this file as example_edge_list.csv. 

Loading this network and applying an algorithm are done by  

.. code-block:: python

   import networkx as nx
   import pandas as pd
   import cpalgorithm as cp
   
   df = pd.read_csv('example_edge_list.csv', sep=',') # Load a list of edges (comma-separated file)
   
   G = nx.from_pandas_edgelist(df) # NetworkX graph object
   
   algorithm = cp.KM_ER()
   algorithm.detect(G)
   c = algorithm.get_pair_id()
   x = algorithm.get_coreness()
 
   print('Name\tPairID\tCoreness')
   for key, value in sorted(c.items(), key=lambda x: x[1]):
       print('%s\t%d\t%f' %(key, c[key], x[key]))


which displays

.. code-block:: bash

   Name	PairID	Coreness
   OQS	0	0.000000
   LNP	0	1.000000
   QSU	0	0.000000
   SUW	0	0.000000
   KMO	0	1.000000
   MOQ	0	1.000000
   PRT	0	0.000000
   JLN	0	1.000000
   NPR	0	1.000000
   RTV	0	0.000000
   IKM	1	0.000000
   HJL	1	0.000000
   BDF	1	1.000000
   @BD	1	1.000000
   EGI	1	0.000000
   GIK	1	0.000000
   FHJ	1	0.000000
   CEG	1	1.000000
   ACE	1	1.000000
   DFH	1	1.000000

