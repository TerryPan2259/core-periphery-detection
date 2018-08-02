=======================
Tutorial
=======================

Preparing a graph
---------------------

cpalgorithm takes a graph object implemented in NetworkX as input 
For example, making an empty graph is done by 

.. code-block:: python

  import networkx as nx
  G = nx.Graph()

One can create a graph object from a file "example.csv": 

.. code-block:: python

  import networkx as nx
  G = nx.read_edgelist("example.csv")

.. note:: 

  The "example.csv" is a space-separated file consisting of two columns, where
  each row corresponds to a pair of adjacent nodes connected by an edge.  

See details in `NetworkX documentation <https://networkx.github.io/documentation/stable/>`_.

 
Core-periphery detection
------------------------
.. role:: python(code)
    :language: python


We'll demonstrate how to use the KM algorithm [1].


Create an object called KM_config:

.. code-block:: python
  
   import cpalgorithm as cp
   import networkx as nx

   algorithm = cp.KM_config()


Then, pass the network to :python:`detect()` method of the algorithm:

.. code-block:: python
  
   G = nx.karate_club_graph() # loading the karate club network
 
   algorithm.detect(G)

Finally, retrieve the results by


.. code-block:: python
  
   c = algorithm.get_pair_id()
   x = algorithm.get_coreness()
  
which gives two dictionaries :python:`c` and :python:`x`.
Dictionary :python:`c` consists of the names of nodes as keys and ID of the core-periphery pair to which the nodes belong as values, e.g.,  
 
.. code-block:: python

   c = {NodeA: 0, NodeB: 1, NodeC: 0, NodeD: 2 ..., 

In dictionary :python:`x`, keys are the node names, and values indicate coreness values, where coreness value 1 and 0 indicates a core or a peripheral node, respectively, e.g.,  

.. code-block:: python

   x = {NodeA: 1, NodeB: 1, NodeC: 0, NodeD: 1 ...,

Note that coreness values can be float values, indicating the extent to which the node belongs to the core. 

One can use other algorithms in the same way. 
For example, one can use the Borgatti-Everet algorithm as follows. 
 

.. code-block:: python
  
   import cpalgorithm as cp
   import networkx as nx

   algorithm = cp.BE()

   G = nx.karate_club_graph() # loading the karate club network
   algorithm.detect(G)
 
   c = algorithm.get_pair_id()
   x = algorithm.get_coreness()

The available algorithms are listed in :ref:`reference`. 


Statistical test
------------------------

Core or peripheral nodes may largely correspond to large-degree or small-degree nodes, respectively.
A question prompted by this observation is that does the detected core-periphery structure reveal something that cannot be explained by the degrees of nodes?
To examine this point, cpalgorithm provides a statistical test for individual core-periphery pair [3].
The statistical test examines to what extent can the detected core-periphery structure be explained by the degree. 
One can carry out the statistical test by writing a line of code: 

.. code-block:: python

   significant, p_values = cp.qstest(c, x, G, algorithm)

where :python:`significant` and :python:`p_values` are list objects.
List :python:`significant` is a boolean list, where :python:`significant[c]=True` or :python:`significant[c]=False` flag indicates that the cth core-periphery pair is significant (i.e., cannot be explained by nodes' degrees) or insignificant (i.e., can be explained by the nodes' degree), respectively, e.g.,  

.. code-block:: python

   significant = [True, False, False, True, ...,

List :python:`p_values` is a float list, where :python:`p_values[c]` is the p-value of the cth core-periphery pair under the configuration model, e.g.,  

.. code-block:: python

   p_values = [0.00001, 0.587, 0.443, 0.0001, ...,

.. note:: 

  The statistical test examines the significance of each core-periphery pair individually, which causes the multiple-comparisons problem. 
  To suppress the false positives, we adopt the e Šidák correction. 
  In other words, the :python:`p_values` is computed first. Then, we apply the Šidák correction. The result is  :python:`significant`.
  The default significance level is 0.05.
  

References
----------

- [1] S. Kojaku and N. Masuda, N. J. Phys. 20, 043012 (2018)
- [2] S. P. Borgatti and M. G. Everett, Soc. Netw. 21, 375 (2000) 
- [3] S. Kojaku and N. Masuda, Sci. Rep. 8, 7351 (2018)


