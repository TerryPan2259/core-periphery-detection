=======================
Overview of cpalgorithm
=======================

cpalgorithm is a Python package for finding core-periphery structure in networks.

Despite a growing interest to core-periphery structure in networks, there are a few libraries for the analysis of core-periphery structure. 
The aim of cpalgorithm is to provide core-periphery detection algorithms.

cpalgorithm provides:
 
* various algorithms for finding single or multiple core-periphery pairs in networks, and
* statistical tests for detected core-periphery structure.


=================================
What is core-periphery structure?
=================================

Core-periphery structure is a mesoscale structure of networks, where a core is a group of densely interconnected nodes, and a periphery is a group of sparsely interconnected nodes. 
Core-periphery structure has been found in various empirical networks such as social networks, biological networks and transportation networks [1, 2].
The core and peripheral nodes may play different roles in networks. 
An example is social network, where core and peripheal nodes often correspond to leaders and followers [1].
Another example is transportation network, where core and peripheral nodes largely correspond to hub and regional airports [3].  
 

Networks may contain a single core-periphery pair or multiple core-periphery pairs.
Figure 1 shows the multiple core-periphery pairs in a political blog networks [4] detected by KM-ER algorithm [3].   
There are two detected core-periphery pairs, each of which consists of the blogs sharing the same political agenda.

.. figure:: fig/poliblog.png
   :scale: 70 %
   :align: center 

**Figure 1**: Two core-periphery pairs in the political blog network detected by our algorithm. 
The filled or open circles indicate core nodes or peripheral nodes, respectively.  
The colour of each circle indicates the core-periphery pair to which the node belongs.
The open circles with a grey border indicate nodes that do not belong to any core-periphery pair. 


References
----------

- [1] P. Csermely, A. London, L.-Y. Wu, and B. Uzzi, `Journal of Complex Networks, 1, 93 (2013) <https://doi.org/10.1093/comnet/cnt016>`_
- [2] M. P. Rombach, M. A. Porter, J. H. Fowler, P. J. Mucha, `SIAM Review, 59, 619-646 (2017) <https://doi.org/10.1137/17M1130046>`_
- [3] S. Kojaku and N. Masuda, `Physical Review E, 96, 052313 (2017) <https://doi.org/10.1103/PhysRevE.96.052313>`_
- [4] L. A. Adamic and N. Glance, `in Proceedings of the 3rd International Workshop on Link Discovery, 36–43 (ACM, New York, USA, 2005) <https://doi.org/10.1145/1134271.1134277>`_

.. toctree::
   :maxdepth: 2 
   :caption: Contents:
   :glob:

   Installation 
   Tutorial 
   Reference 
   Examples
   FAQ/FAQ 
   Contact 

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
