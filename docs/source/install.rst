Installation
============

Using pip
---------

.. code:: bash

    $ pip install opentiva


Build from source
-----------------

.. code:: bash

    $ git clone https://github.com/opentiva/opentiva
    $ cd opentiva/
    $ pip install .

Both pip and the above installation method use the precompiled opentiva.pkpd 
cython module.

To rebuild the opentiva.pkpd cython module during package installation install 
the cython package prior to running pip install or building from source. For 
example:

.. code:: bash

    $ git clone https://github.com/opentiva/opentiva
    $ cd opentiva/
    $ pip install cython
    $ pip install .
