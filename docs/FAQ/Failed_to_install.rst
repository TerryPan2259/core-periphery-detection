.. _failed_to_install_cpalgorithm:

.. role:: python(code)
    :language: python

#############################
Failed to install cpalgorithm 
#############################

cpalgorithm contains c++ extension modules which need to be compiled.

First, please make sure that you have one of the following compilers:

 - Clang/LLVM 3.3 or newer (for Apple Xcode's clang, this is 5.0.0 or newer)
 - GCC 4.8 or newer
 - Microsoft Visual C++ Build Tools 2015 or newer
 - Intel C++ compiler 17 or newer
 - Cygwin/GCC (tested on 2.5.1)

Also, you need to install `pybind11 <https://github.com/pybind/pybind11>`_. 

========================================
cpalgorithm is not compatible with Conda
========================================

Installing packages via pip may destroy `Conda <https://conda.io/docs/index.html>`_.
Consider installing the package from source. See :ref:`install_from_source`. 
 

====================
Check python version
====================

cpalgorithm is not compatible with python 2.x. 
Use python version 3.4 or newer. 


=========
Change OS
=========

The library is frequently tested on Ubuntu 16.04 and CentOS 7, and occasionally on Windows 10 and MacOS.
So if you failed to install on Windows and MacOS, then consider using Ubuntu or CentOS. 

.. _install_from_source:

===================
Install from source
===================

If you failed to install via pip, another option is to build from source. 
Download the source code from GitHub, e.g., 

.. code-block:: bash

   git clone https://github.com/skojaku/core-periphery-detection


Then, move to the directory 

.. code-block:: bash

   cd core-periphery-detection

Type

.. code-block:: bash

   python -m pip install .

If you are luckily, you can install the package (Congratulations!). 
If you could not, then find out the error messages. 
In many cases, the compilers do not support std=c++11 or do not meet the requirements mentioned above.
In this case, upgrade the compilers.

 
