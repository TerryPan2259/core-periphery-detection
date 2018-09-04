=======================
Tutorial
=======================

Preparing a graph
-----------------

The algorithms in cpalgorithm take NetworkX graph object as input. 
Making an empty graph is done by 

.. code-block:: python

  import networkx as nx
  G = nx.Graph()

One can create a graph object from a file "example.csv": 

.. code-block:: python

  import networkx as nx
  G = nx.read_edgelist("example.csv")

.. note:: 

  The "example.csv" is a space-separated file consisting of two columns, where
  each row corresponds to a pair of adjacent nodes connected by an edge. See :ref:`examples`.

See details in `NetworkX documentation <https://networkx.github.io/documentation/stable/>`_.

 
Core-periphery detection
------------------------
.. role:: python(code)
    :language: python

cpalgorithm contains several algorithms to find core-periphery structure in networks.
Here, we apply the KM-config algorithm to the karate club network.

Create an object called KM_config:

.. code-block:: python
  
   import cpalgorithm as cp
   import networkx as nx

   algorithm = cp.KM_config()

Load the karate club network:
 
.. code-block:: python

   G = nx.karate_club_graph() # loading the karate club network


Then, pass the graph to :python:`detect()` method:

.. code-block:: python
  
   algorithm.detect(G)


Retrieve the results by

.. code-block:: python
  
   c = algorithm.get_pair_id()
   x = algorithm.get_coreness()

  
:python:`c` and :python:`x` are python dict objects. 
Dictionary :python:`c` takes keys representing the node names, and integer values representing the IDs of the core-periphery pair to which to the node belongs.  
For example,
 
.. code-block:: python

   c = {NodeA: 0, NodeB: 1, NodeC: 0, NodeD: 2 ..., 

means that NodeA and NodeC belong to core-periphery pair 0, NoedB belongs to core-periphery pair 1 and NodeD belongs to core-periphery pair 2. 

Dictionary :python:`x` takes keys representing the node names, and float values representing the coreness values ranging between 0 and 1.
Coreness value 1 and 0 indicates a core or a peripheral node, respectively. For example, 

.. code-block:: python

   x = {NodeA: 1, NodeB: 1, NodeC: 0, NodeD: 1 ...,

means NodeA, NodeB NodeD are core nodes and NodeC is a peripheral node. 
Note that some algorithms set coreness values between 0 and 1, which indicates the extent to which the node belongs to the core. 


One can use other algorithms in the same way. 
For example, one needs to modify one line to use the Borgatti-Everet algorithm (e.g, cp.BE()). 

.. code-block:: python
  
   import cpalgorithm as cp
   import networkx as nx

   algorithm = cp.BE()
   #algorithm = cp.KM_config()

   G = nx.karate_club_graph() 
   algorithm.detect(G)
 
   c = algorithm.get_pair_id()
   x = algorithm.get_coreness()

The available algorithms are listed in :ref:`reference`. 


Statistical test
----------------

Likewise rich-club and various centrality measures, heterogeneous degree distribution alone may explain core-periphery structure.
cpalgorithm provides a statistical test to examine the significance of individual core-periphery pairs. 
The statistical test judges each detected core-periphery pair as insignificant if it can be explained largely by the degree (i.e., hub and non-hub nodes largely correspond to core and peripheral nodes, respectively). Otherwise, it judges a core-periphery pair as significant.  
One can carry out the statistical test by writing a line of code: 

.. code-block:: python

   sig_c, sig_x, significant, p_values = cp.qstest(c, x, G, algorithm)

where :python:`significant` and :python:`p_values` are list objects.
`sig_c` and `sig_x` are dict objects in which the insignificant core-periphery pairs are excluded. 
List :python:`significant` is a boolean list, where :python:`significant[c]=True` or :python:`significant[c]=False` flag indicates that the cth core-periphery pair is significant or insignificant, respectively, e.g., 

.. code-block:: python

   significant = [True, False, False, True, ...,

List :python:`p_values` is a float list, where :python:`p_values[c]` is the p-value of the cth core-periphery pair under the configuration model, e.g.,  

.. code-block:: python

   p_values = [0.00001, 0.587, 0.443, 0.0001, ...,

.. note:: 

  The statistical test examines the significance of each core-periphery pair individually, which causes the multiple-comparisons problem. 
  To suppress the false positives, we adopt the e Šidák correction. 
  The default significance level is 0.05.


