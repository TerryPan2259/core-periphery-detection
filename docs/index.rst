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
Many studies have assumed that networks consist of one core and one periphery.  
However, networks may contain multiple core-periphery pairs.

We developed an algorithm to find multiple core-periphery pairs in networks [3].
In the political blog network [4] shown in Fig. 1, the algorithm found two core-periphery pairs, each of which consists of the blogs sharing the same political agenda.
We applied this algorithm or its variant to an inter-bank trade network [5], a maritime transportation network [6] and a short-text classification task [7].

We also showed that heterogeneous degree distributions alone explain single core-periphery structure [8]. In other words, when one says that a network is composed of a single core and a single periphery, hubs (i.e. nodes with large degrees, or the number of edges) are core nodes and nodes with small degrees are peripheral nodes. We proved that a core-periphery structure that is not merely explained by heterogeneous degree distributions necessarily involves at least three blocks of nodes. An example is a combination of a single core, a single periphery and a community, which yields three blocks. Another example is two core-periphery pairs, which yields four blocks in total. Based on this result, we extended our first algorithm [7] to find multiple core-periphery pairs in networks that are not merely explained by heterogeneous degree distributions [8].
The extended algorithm allows one to find small-degree core nodes and large-degree peripheral nodes.

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
- [5] S. Kojaku, G. Cimini, G. Caldarelli, N. Masuda, Journal of Network Theory in Finance, In press (2018)
- [6] X. Cui, S. Kojaku, N. Masuda, and D. Bollegala, `in Proceedings of the 7th Joint Conference on Lexical and Computational Semantics, 225-264 (ACL, New Orleans, USA, 2018) <https://doi.org/10.18653/v1/S18-2030>`_
- [7] S. Kojaku, M. Xu, H. Xia, N. Masuda, `Preprint arXiv:1808.04549 (2018) <https://arxiv.org/abs/1808.04549>`_
- [8] S. Kojaku and N. Masuda, `New Journal of Physics 20, 043012 (2018) <https://doi.org/10.1088/1367-2630/aab547>`_

.. toctree::
   :maxdepth: 2 
   :caption: Contents:
   :glob:

   Installation 
   Tutorial 
   FAQ/FAQ 
   Examples
   Reference 
   Contact 

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
