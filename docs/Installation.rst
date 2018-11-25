############
Installation
############

.. role:: python(code)
    :language: python

Requirements
------------
cpalgorithm is compatible with python 3.4 or newer.  

You need to have the latest :python:`pip` (a package manager for Python) and one of the following compilers: 

 - Clang/LLVM 3.3 or newer (for Apple Xcode's clang, this is 5.0.0 or newer)
 - GCC 4.8 or newer
 - Microsoft Visual C++ Build Tools 2015 or newer
 - Intel C++ compiler 17 or newer
 - Cygwin/GCC (tested on 2.5.1)

Currently, cpalgorithm is tested on Ubuntu and CentOS.
We confirmed that cpalgorithm can also run on MacOS. 
For Windows, cpalgorithm requires Visual C++ Build Tools 2015 or newer. To install, see `Windows Compilers <https://wiki.python.org/moin/WindowsCompilers/>`_.

cpalgorithm may not be compatible with conda.

Install
-------

Before installing cpalgorithm, install pybind11:

.. code-block:: bash

  $ pip3 install pybind11

Then, 

.. code-block:: bash

  $ pip3 install cpalgorithm

If you don't have a root access, try :python:`--user` flag:

.. code-block:: bash

  $ pip3 install --user cpalgorithm

You can upgrade a newer release by 
  
.. code-block:: bash

  $ pip3 install --upgrade cpalgorithm


