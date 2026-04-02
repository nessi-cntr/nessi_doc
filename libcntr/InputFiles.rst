.. _PMan14:

Creating and Passing Input Files
================================

.. contents::
   :local:
   :depth: 2

Along with the ``libcntr`` library, we also provide an input file parser, which can be used in custom programs by including:

.. code-block:: cpp

   #include "cntr/utils/read_inputfile.hpp"

Input files allow for passing variables to ``NESSi``-based programs. We use the following format:

.. code-block:: sh

   __var=value

Note the double underscores in front of the variable name. Let us consider an example input file ``input.txt``:

.. code-block:: sh

   __nt=1000
   __ntau=400
   __h=0.01
   __beta=10.0

which is parsed by the code snippet:

.. code-block:: cpp

   int nt,ntau;
   double h,beta;

   find_param("input.txt", "__nt=", nt)
   find_param("input.txt", "__ntau=", ntau)
   find_param("input.txt", "__h=", h)
   find_param("input.txt", "__beta=", beta)

Furthermore, the input file parser can read vectors from file. A typical application is a time-dependent parameter such as an external electric field :math:`E(t)`. We assume the field to be represented on the same grid as in the calculation: :math:`E_j = E(j h)`, :math:`j=0,\dots,n_t`. We also include the value in thermal equilibrium :math:`E_{-1}` (hence, :math:`E(t)` is effectively represented by a ``cntr::function``). Preparing an input file ``field.txt`` like:

.. code-block:: sh

   0.0
   0.0
   0.1
   0.2
   ...

with ``nt+2`` entries, we can read the electric field by referring to ``field.txt`` in the input file ``input.txt``:

.. code-block:: sh

   __Efield=--field.txt

In the C++ program, the field can be read as follows:

.. code-block:: cpp

   std::vector<double> Efield;

   read_tvector("input.txt", "__Efield=", Efield, nt);

Resizing to the corresponding size is part of the function ``read_tvector``.

Multidimensional data can also be read. Extending the above example to a two-dimensional electric field, the input file ``field.txt`` would look like:

.. code-block:: sh

   0.0  0.0
   0.0  0.0
   0.1  0.1
   0.2  0.2
   ...

Reading this field in a program:

.. code-block:: cpp

   int dim=2;
   std::vector<dvector> Efield;

   read_tvector("input.txt", "__Efield=", Efield, dim, nt);

To ease the usage of input files, we provide the python module ``ReadCNTR.py`` (available for both python2 and python3), which allows to create input files from python dictionaries. For more details see :ref:`PMan13`.
