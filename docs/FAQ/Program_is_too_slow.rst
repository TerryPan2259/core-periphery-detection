.. _program_is_too_slow:

.. role:: python(code)
    :language: python

############################
My program is too slow. Why? 
############################

cpalgorithm can be slow for some reasons and here is the workaround.


=====================
Combine multi-edges
=====================

If there are multiple edges between the same pair of nodes in your neworkx graph object :python:`G`, 
consider aggregating the multiple edges into one edge, that is to replace the multiple edges with a single edge having a weight, where the weight indicates the number of edges (or the sum of the weight over the multiple edges) between the pair of nodes. 

This can be done by 

.. code-block:: python

   Gnew = nx.Graph()
   for u,v,data in G.edges(data=True):
       w = data['weight'] if 'weight' in data else 1.0
       if Gnew.has_edge(u,v):
           Gnew[u][v]['weight'] += w
       else:
           Gnew.add_edge(u, v, weight=w)

where :python:`G` is the network containing multi-edges, and :python:`Gnew` is the network without multiple edges. 


=================
Shorten node name
=================

cpalgorithm internally converts nodes' names, which can be string or number, into unique integers. 
This process can be slow if the nodes' names are long.

One can check the nodes' names by  

.. code-block:: python

   G.nodes()

One can reduce the computational time by shortening the name of each node, e.g., integers.


======================
Use parallel computing
======================

Some algorithms provide different results on different runs. 
cpalgorithm runs such algorithms several times with different random seeds. 
Then, it chooses the one yielding the largest quality value.
If we run the algorithm in parallel, then computational time would be reduced. 

Here is an example of how to do this:
 
.. code-block:: python
   :linenos:

   import networkx as nx
   import cpalgorithm as cp
   import numpy as np
   from multiprocessing import Pool
   
   # A function for reducing computation time in algorithm using parallel computing  
   def par_detect_cp(num_of_cores, num_runs = 10):
   	pool = Pool(num_of_cores)
   	results = pool.map(_detect_cp, list(range(num_runs)))
   	algorithm = results[np.argmax(list(map(lambda x : sum(x.score()), results)))]
   	pool.close()
   	return algorithm
   
   # This is an internal function of detect_cp
   def _detect_cp(_rubbish):
   	alg = cp.KM_config(num_runs = 1)
   	alg.detect(G)
   	return alg
   
   
   # Construct Graph object
   global G # Declare the graph object as a global variable to save memory 
   G = nx.karate_club_graph()
   
   algorithm = par_detect_cp(num_of_cores = 10, num_runs=10)
   
   c = algorithm.get_pair_id()
   x = algorithm.get_coreness()
   
   print('Name\tPairID\tCoreness')
   for key, value in sorted(c.items(), key=lambda x: x[1]):
   	print('%s\t%d\t%f' %(key, c[key], x[key]))


In this example, the networkx graph object is constructed in lines 22 and 23. 
The algorithm is set in line 16.
In line 25, we specify the number of cores by :python:`num_cores=10`.
:python:`num_runs` is the number of times that cpalglrithm runs the algorithm.
Default is :python:`num_runs=10`. :python:`num_cores` should be equal or smaller than  :python:`num_runs`.
