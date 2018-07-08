=======================
Tutorial
=======================

Preparing a graph
---------------------

cpalgorithm takes a graph object implemented NetworkX as input 
For example, making an empty graph is done by 

.. code-block:: python

  import networkx as nx
  G = nx.Graph()

One can create the graph object from a file: 

.. code-block:: python

  import networkx as nx
  G = nx.read_edgelist("example.csv")

.. note:: 

  The "example.csv" consists of two columns, where
  each row corresponds to a pair of adjacent nodes connected by an edge. 

See details in NetworkX documentation.

 
Core-periphery detection
------------------------
.. role:: python(code)
    :language: python

We'll demonstrate how to use the KM algorithm.


Create an object called KM_config:

.. code-block:: python
  
   import cpalgorithm as cp
   import networkx as nx

   algorithm = cp.KM_config()


Then, pass the network of interest to :python:`detect()` method of the algorithm:

.. code-block:: python
  
   G = nx.karate_club_graph() # loading the karate club network
 
   algorithm.detect(G)

Finally, retrieve the results by


.. code-block:: python
  
   c = algorithm.get_pair_id()
   x = algorithm.is_core()
  
which gives two dictionaries :python:`c` and :python:`x`.
Dictionary python:`c` consists of the names of nodes as keys and ID of the core-periphery pair to which the nodes belong as values, e.g.,  
 
.. code-block:: python

   c = {NodeA: 0, NodeB: 1, NodeC: 0, NodeD: 2 ..., 

Dictionary python:`x` consists of the names of nodes as keys and flags as values, where True or False flag indicates a core or a periphery, respectively, e.g.,  

.. code-block:: python

   x = {NodeA: True, NodeB: True, NodeC: False, NodeD: False ...,


Statistical test
------------------------

Detected core-periphery structure may be insignificant, e.g., the core-periphery structure may not be that strong. 
To examine this point, one needs to write a line of additional code. 


.. code-block:: python

   significant, p_values = cp.qstest(c, x, G, algorithm)

where :python:`significant` and :python:`p_values` are list onjects.



 
