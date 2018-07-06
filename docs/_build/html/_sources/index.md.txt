=======================
Overview of cpalgorithm
=======================

cpalgorithm is a Python package for finding core-periphery structure in networks.

cpalgorithm provides:
 
* various algorithms for finding single or multiple core-periphery pairs in networks, and
* statistical tests for detected core-periphery structure.

Despite the importance of core-periphery structure, a few software implements core-periphery detection algorithms. 
Some authors provide the code for their algorithms, which is, however, often written in different programming languages and requires users to understand the usage of each. 
The aim of cpalgorithm is to provide an easy access to the algorithms.  

========================
Core-periphery structure
========================

Core-periphery structure is a mesoscale structure of networks, where core is a group of densely interconnected nodes, and periphery is a group of sparsely interconnected nodes. 
The core can be densely interconnected with the periphery or not.
Core-periphery structure has been found in various empirical networks such as social networks, biological networks and transportation networks.
For example, a political blog network consists of two core-periphery pairs, each of which consists of the blogs sharing the same political agenda.
Each core-periphery pair consists of a set of core blogs that are linking each other and a set of peripheral blogs that links mainly to the core blogs.  

[some figures go here]

In a worldwide airport network, there are 28 core-periphery pairs, each of which mostly consists of the airports in the same geographical region.
The core and peripheral ports largely correspond to the hubs and regional airports. 

[some figures go here]

Core-periphery structure may play a crucial role in complex systems.
For example, core-periphery structure is the most robust structure against random attacks to nodes. 
Transportation networks with core-periphery structure are economically efficient.


========================
Scope
========================

Among various algorithms, cpalgorithm focuses on the algorithms for finding density-based core-periphery structure in networks.
Other types of core-periphery structure such as transportation-based core-periphery structure is out of scope in the current version. 

 
.. toctree::
   :maxdepth: 2
   :caption: Contents:

   Installation 
   Examples 
   Reference 

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
