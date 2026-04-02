.. _PMan02:

Timeslices
==========

.. contents::
   :local:
   :depth: 2

.. _PMan02S01:

Overview
--------

.. list-table::
   :header-rows: 0

   * - class
     - ``cntr::herm_matrix_timestep<T>``

For a Green's function :math:`C` of type ``cntr::herm_matrix``, we define the timestep ``tstp``, with ``tstp`` :math:`\in\{` ``-1,0,1,...,nt`` :math:`\}`, as the following data set:

- If ``tstp=-1``, the timestep contains the Matsubara Green's function:

  - Matsubara component :math:`C^\mathrm{M}(j\Delta \tau)` for ``j=0,...,sntau``,

- If ``tstp>=0``, it contains the following real-time components:

  - :math:`C^\mathrm{R}({\tt tstp}\Delta t,j \Delta t)` for ``j=0,...,tstp``,
  - :math:`C^<(j \Delta t,{\tt tstp} \Delta t)` for ``j=0,...,tstp``,
  - :math:`C^\rceil({\tt tstp}\Delta t,j\Delta \tau)` for ``j=0,...,ntau``.

The data of the ``herm_matrix`` Green's function :math:`C` are therefore precisely a union of the timesteps ``-1,...,nt``. The class ``cntr::herm_matrix_timestep<T>`` is the container to store these data. It is characterized by the following parameters:

- ``T`` (template parameter): Precision, usually set to ``double``; we use the definition

  .. code-block:: cpp

     #define GREEN_TSTP cntr::herm_matrix_timestep<double>

- ``tstp`` (integer): The timestep variable, ``tstp>=-1``.
- ``ntau`` (integer): number of discretization points on the imaginary time axis, ``ntau>=0``
- ``size1`` (integer): orbital dimension. Each element :math:`C(t,t')` is a square matrix of dimension ``size1`` :math:`\times` ``size1``.
- ``sig`` (``FERMION`` or ``BOSON``).

For a summary of all member functions, see :ref:`PMan01S02Special`. The following paragraphs give detailed explanations of some functionalities.

.. _PMan02S02:

Constructors
------------

.. list-table::
   :header-rows: 0

   * - ``herm_matrix_timestep<T>()``
     - Default constructor, does not allocate memory and sets ``nt=-2``.
   * - ``herm_matrix_timestep<T>(int tstp,int ntau,int size1,int sig)``
     - Allocate memory, set all entries to ``0``. It requires ``tstp>=-1``, ``ntau>0``, ``size1>0``, and ``sig=FERMION`` or ``sig=BOSON``.

.. _PMan02S03:

Accessing individual matrix elements
-------------------------------------

Element access, as well as read-out of the density matrix follows the same syntax as for the Green's function ``herm_matrix<T>``. See Sections :ref:`PMan01S04` and :ref:`PMan01S03`.

**Example 1:**

.. code-block:: cpp

   int ntau=100;
   int size=2;
   int tstp=10;
   cdmatrix M(size1,size); // an Eigen matrix
   GREEN_TSTP tA(tstp,ntau,size,FERMION); // allocate a timestep
   // ... set entries of M ...
   tA.set_les(8,tstp,M);  // ok: A^<(8,10) set to M
   tA.set_les(tstp,8,M);  // error: time arguments do not match the domain of tA
   cntr::get_les(tstp,8,M,tA,tA); // ok: M set to A^<(8,10), assuming that A is hermitian

**Example 2:**

Simple print of all elements in case of a scalar Green's function, at tstp=-1 (Matsubara):

.. code-block:: cpp

   int ntau=100;
   int size=1;
   tstp=-1;
   GREEN_TSTP tA(tstp,ntau,size,FERMION);
   // ... do something with tA ...
   for(int j=0;j<=ntau;j++){
     cdouble x;
     tA.get_mat(j,x);
     cout << "Amat(j) at j=" << j << " : "  << x << endl;
   }

Remark: Some un-supported routines exist in which the dummy argument ``tstp`` can be omitted.

.. _PMan02S04:

Timestep-wise data exchange between Green's functions
------------------------------------------------------

In typical applications, data between different Green's functions are exchanged at once for an entire timestep. We provide two main member functions which allow the exchange of data of one entire timestep between variables of type ``cntr::herm_matrix`` and ``cntr::herm_matrix_timestep``:

.. list-table::
   :header-rows: 0

   * - ``A.set_timestep(tstp, B)``
     - Copies the data of timestep ``tstp`` from ``B`` to ``A``

- ``A`` and ``B`` are either ``cntr::herm_matrix<T>`` or ``cntr::herm_matrix_timestep<T>``; Types can be mixed to enable data exchange between the types in both directions
- ``tstp`` (integer): The timestep to be accessed
- For ``tstp==-1`` this copies the data from ``B`` to ``A`` as:

  - :math:`A^M(i\Delta \tau) = B^M(i\Delta \tau)` for ``i=0,...,ntau``

- For ``tstp>=0`` this copies the data from B to A as:

  - :math:`A^R({\tt tstp}\Delta t,j \Delta t) = B^R({\tt tstp}\Delta t,j \Delta t)` for ``j=0,...,tstp``
  - :math:`A^<(j \Delta t,{\tt tstp}\Delta t) = B^<(j \Delta t, {\tt tstp}\Delta t)` for ``j=0,...,tstp``
  - :math:`A^{\rceil}({\tt tstp}\Delta t,j\Delta\tau) = B^{\rceil}({\tt tstp}\Delta t,j\Delta\tau)` for ``j=0,...,ntau``

- Size consistency is required:

  - If ``A`` is of type ``herm_matrix_timestep``, ``tstp==A.tstp()`` is required, same for ``B``
  - If ``A`` is of type ``herm_matrix``, ``tstp<=A.nt()`` is required, same for ``B``
  - ``A.ntau()==B.ntau()`` and ``A.size1()==B.size1()`` is required

**Example:**

The following code copies the first 5 timesteps and the Matsubara timestep from a Green's function ``B`` to a Green's function ``A``:

.. code-block:: cpp

   int nt=10;
   int ntau=10
   int size1=2;
   GREEN A(nt,ntau,size1,FERMION);
   GREEN B(5,ntau,size1,FERMION);
   // ... do something with B ...
   for(int tstp=-1;tstp<=5;tstp++) A.set_timestep(tstp,B);

.. list-table::
   :header-rows: 0

   * - ``A.set_matrixelement(tstp, i1,i2,B,j1,j2)``
     - Read out only certain (orbital) matrix elements at timestep ``tstp`` from ``B`` and write them to ``A``

- ``A`` and ``B`` are either ``cntr::herm_matrix`` or ``cntr::herm_matrix_timestep``; The two types can be mixed.
- ``tstp`` (integer): The timestep to be accessed
- ``i1``, ``i2``, ``j1``, ``j2`` (integer): Orbital matrix indices
- For ``tstp==-1`` this copies the data from ``B`` to ``A`` as:

  - :math:`A^M(i\Delta \tau)_{i1,i2} = B^M(i\Delta \tau)_{j1,j2}` for ``i=0,...,ntau``

- For ``tstp>=0`` this copies the data from B to A as:

  - :math:`A^R({\tt tstp}\Delta t,i \Delta t)_{i1,i2} = B^R({\tt tstp}\Delta t,i \Delta t)_{j1,j2}` for ``i=0,...,tstp``
  - :math:`A^<(i \Delta t,{\tt tstp}\Delta t)_{i1,i2} = B^<(i \Delta t,{\tt tstp}\Delta t)_{j1,j2}` for ``i=0,...,tstp``
  - :math:`A^{\rceil}({\tt tstp}\Delta t,i\Delta\tau)_{i1,i2} = B^{\rceil}({\tt tstp}\Delta t,i\Delta\tau)_{j1,j2}` for ``i=0,...,ntau``

- Size consistency is required:

  - If ``A`` is of type ``herm_matrix_timestep``, ``tstp==A.tstp()`` is required, same for ``B``
  - If ``A`` is of type ``herm_matrix``, ``tstp<=A.nt()`` is required, same for ``B``
  - ``A.ntau()==B.ntau()`` is required
  - ``0 <= i1,i2 <= A.size1()``, ``0 <= j1,j2 <= B.size1()``

``set_matrixelement`` is useful in the context of constructing Feynman diagrams with multiorbital Green's functions.

**Example:**

The following (little useful) code permutes the ``(0,1)`` and ``(1,0)`` matrix elements of timestep ``7`` of a ``herm_matrix`` Green's function, using the ``herm_matrix_timestep`` as a temporary variable.

.. code-block:: cpp

   int nt=10;
   int ntau=10
   int size1=2;
   int sig=FERMION;
   GREEN G(nt,ntau,size1,sig);
   // ... do something with G
   {
       int tstp=7;
       GREEN tG(tstp,ntau,1,sig);              // allocate a timestep of orbital dimension 1 as temp. storage
       tG.set_matrixelement(tstp,0,0,G,0,1);   // save (0,1) matrix element
       G.set_matrixelement(tstp,0,1,G,1,0);    // copy (1,0) matrix element to (0,1)
       G.set_matrixelement(tstp,1,0,tG,0,0);   // copy (0,0) of tG to (1,0) matrix-element of G
   }
