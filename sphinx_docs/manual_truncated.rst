.. _Ptrunc:

Manual for truncated Kadanoff-Baym equations
=============================================

.. contents::
   :local:
   :depth: 2

.. _Ptrunc00:

Formalism: Memory truncated time propagation
--------------------------------------------

The Keldysh formalism for nonequilibrium Green's functions presented in :ref:`P3` is a powerful theoretical framework for the
description of the electronic structure, spectroscopy, and dynamics of strongly correlated systems. However,
the underlying Kadanoff-Baym equations (KBE) for the two-time Keldysh Green's functions involve a memory
kernel, which results in a high computational cost for long simulation times ``tmax``, with a cubic scaling of the
computation time with ``tmax``.

Truncation of the memory kernel can reduce the computational cost to linear scaling
with ``tmax``, but the required memory times will depend on the model and the diagrammatic approximation to the
self-energy. We explain here how the truncation algorithm of the memory kernel is implemented in NESSi and show in :ref:`Trunc`
how it can be used for a simple physical system, allowing the long time propagation of the system until the full thermalization is reached.

**Memory Constraint**

We explain the truncation for the Dyson equation. The Volterra integral equation (vie2) works entirely analogous. We aim to solve one of the two equivalent equations

.. math::
   :label: dyson1

   i\partial_tG(t,t^\prime)-\epsilon(t)G(t,t^\prime)-\int_{\mathcal{C}} d \bar{t}\Sigma(t,\bar{t})G(\bar{t},t^\prime)=\delta(t,t^\prime),

.. math::
   :label: dyson2

   -i\partial_{t^\prime}G(t,t^\prime)-G(t,t^\prime)\epsilon(t^\prime)-\int_{\mathcal{C}} d \bar{t}G(t,\bar{t})\Sigma(\bar{t},t^\prime)=\delta(t,t^\prime),

assuming that the self-energy is short-range in time:

.. math::
   :label: constraint

   \Sigma^R(t,t-s)=\Sigma^<(t,t-s)=0\quad \text{for} \quad |s|>t_c.

Here :math:`t_{c}` is some memory cutoff.

**Memory-truncated time slices and windows**

We explain the solution of these equations on an equidistant time mesh of timestep :math:`\Delta t`. In the following, we choose the cutoff :math:`t_c={\tt tc}\Delta t`, integer ``tc``. For any contour function :math:`X`, we define the truncated timeslice

.. math::
   :label: timeslice

   \mathcal{T}[X]_n^{\tt tc}\equiv \{X^R(n\Delta t, n\Delta t-m\Delta t),X^<(n\Delta t, n\Delta t-m\Delta t),0\leq m \leq {\tt tc}\},

We also define a moving window

.. math::
   :label: window

   \mathcal{M}[X]_n^{\tt tc}\equiv \bigcup_{m=0}^{\tt tc}\mathcal{T}[X]_{n-m}^{\tt tc},

and a triangular moving window

.. math::
   :label: tri_window

   \Delta[X]_n^{\tt tc}\equiv \bigcup_{m=0}^{\tt tc}\mathcal{T}[X]_{n-m}^{{\tt tc}-m},

The windows are shown in the figure in section :ref:`Ptrunc01S01`, showing the timeslice :math:`\mathcal{T}` as a red shaded area, the triangular window :math:`\Delta` by the solid red line and :math:`\mathcal{M}` by the dashed red contour.

**Timestepping**

The key fact underlying the truncated Kadanoff-Baym equations is the following: If the self-energy satisfies the constraint :eq:`constraint`, then the Dyson equation for the evaluation of :math:`\mathcal{T}[G]_n^{\tt tc}` is closed with the knowledge of :math:`\Delta[\Sigma]_{n}^{\tt tc}` and :math:`\Delta[G]_{n}^{\tt tc}` (where of course the dependence of :math:`\mathcal{T}[G]_{n}^{\tt tc}` on itself due to :math:`\mathcal{T}[G]_{n}^{\tt tc}\subset\Delta[\Sigma]_{n}^{\tt tc}`, just means in the end we solve an implicit equation for :math:`\mathcal{T}[G]_{n}^{\tt tc}`). The previous statement allows a time-propagation scheme for the object :math:`\Delta[G]_{n}^{\tt tc}` instead of a solution of the Kadanoff-Baym equation. The extreme limit is ``tc=0``: :math:`\Delta[G]_{n}^{0}` basically contains the one-particle density-matrix :math:`G^<(n\Delta t, n\Delta t)`. For a time-local equation such as Hartree-Fock (no memory integral), one indeed has a closed equation for the one-particle density matrix.

For the implementation, it is convenient to formulate a closed time-evolution based on the extended window :eq:`window` instead of the triangular :eq:`tri_window`. This is partly redundant in memory (:math:`\Delta[G]_{n}^{\tt tc} \subset \mathcal{M}[G]_{n}^{\tt tc}`) but it allows for an easier memory handling: We can define a time-stepping function, which calculates :math:`\mathcal{T}[G]_{n}^{\tt tc}` based on the knowledge of :math:`\mathcal{M}[G]_{n}^{\tt tc}` and :math:`\mathcal{M}[\Sigma]_{n}^{\tt tc}` with the exclusion of :math:`\mathcal{T}[G]_{n}^{\tt tc}`. After that, one can simply shift the moving window to :math:`\mathcal{T}[G]_{n+1}^{\tt tc}`, where now only :math:`\mathcal{T}[G]_{n+1}^{\tt tc}` is unknown, repeat the dyson timestep. If the self-energy depends on :math:`G`, this dependence is usually such that :math:`\mathcal{T}[\Sigma]_{n}^{\tt tc}` can be calculated from :math:`\mathcal{T}[G]_{n}^{\tt tc}`, so that the above time-stepping can be supplemented with an iterative calculation of :math:`\Sigma`.

The time-stepping of :math:`\mathcal{M}[G]_n^{\tt tc}` is implemented using the class ``herm_matrix_moving<T>`` discussed in detail in :ref:`Ptrunc01`, which contains the time window :math:`\mathcal{M}[X]_{n}^{\tt tc}`.

.. _Ptrunc01:

Truncated Green's function
--------------------------

.. _green_trunc_def:

.. _Ptrunc01S01:

Overview
~~~~~~~~

.. list-table::
   :header-rows: 0

   * - class
     - ``cntr::herm_matrix_moving<T>``

This class contains the data structures for representing two-time contour functions :math:`G(t,t')` on the equidistantly discretized time-contour :math:`\mathcal{C}` with points in a moving time window
:math:`\mathcal{M}[G]^{\tt tc}_{\tt t_{0}}`, see :eq:`window`. ``cntr::herm_matrix_moving<T>`` stores the following Keldysh components on the equidistantly discretized Keldysh contour with time discretization :math:`\Delta t` on the real time axis:

- retarded component :math:`G^\mathrm{R}(t_{0}-i\Delta t,t_{0}-i\Delta t -j\Delta t)` for ``i=0,...,tc``, ``j=0,...,tc``,
- lesser component :math:`G^<(t_{0}-i\Delta t,t_{0}-i\Delta t -j\Delta t)` for ``j=0,...,tc``, ``i=0,...,tc``.

.. note::

   - Here and in the manual below, we will consistently use the notation :math:`t_0` to indicate the leading physical time at which the moving window is attached.
   - The value :math:`t_0` itself is however not stored as a member of the class. The data of the ``cntr::herm_matrix_moving<T>`` object ``G`` are addressed with time indices **relative to the leading time step**, instead of the physical time, for example ``G.get_ret(i,j)`` :math:`\rightarrow G^R(t_0-i\Delta t,t_0-(i+j)\Delta t)`.
   - We will refer to the data with a given ``i`` as **timestep** or **timeslice** ``i`` of the moving window. Timestep ``0`` will also be referred to as the **leading timestep**. To store data at one timeslice, one can use the class ``cntr::herm_matrix_timestep_moving<T>``, see :ref:`Ptrunc02`.
   - Note that the values in the storage domain of a ``cntr::herm_matrix_moving<T>`` are sufficient to represent a Green's function which has **hermitian symmetry** (:math:`C=C^\ddagger`). If needed, in order to reconstruct a general two-time-function :math:`C` one can store the storage domain of :math:`C` and of :math:`C^\ddagger` in two separate ``cntr::herm_matrix_moving<T>`` variables.

The truncated Green's functions are defined by the following parameters:

- ``T`` (template parameter): Precision, usually set to ``double``; we use the global definition

  .. code-block:: cpp

     #define GTRUNC cntr::herm_matrix_moving<double>

- ``tc`` (integer): number of discretization points on the real time axis.
- ``size1`` (integer): orbital dimension. Each element :math:`C(t,t')` is a square matrix of dimension ``size1`` :math:`\times` ``size1``.
- ``sig`` (integer): Takes the values ``FERMION`` (defined as ``-1``) or ``BOSON`` (``+1``).

Size arguments of ``cntr::herm_matrix_moving<T>`` are returned by the following read-only member functions:

.. code-block:: cpp

   GTRUNC G;
   cout << "The object G of type herm_matrix_moving has the following dimensions:" << endl;
   cout << "tc= " << G.tc() << endl;
   cout << "size1= " << G.size1() << endl;
   cout << "sig= " << G.sig() << endl;

.. _Ptrunc01S02:

Constructors
~~~~~~~~~~~~

.. list-table::
   :header-rows: 0

   * - ``herm_matrix_moving<T>()``
     - Default constructor, does not allocate memory, sets ``tc=-1``
   * - ``herm_matrix_moving<T>(int tc,int size1,int sig)``
     - Allocate memory, and sets all entries to ``0``. It requires ``tc>0``, ``size1>0``, and ``sig`` :math:`\in\{` ``FERMION``, ``BOSON`` :math:`\}`.
   * - ``herm_matrix_moving<T>(herm_matrix_moving<T> &G)``
     - Copy assignment: Constructs a new instance with its own data.

.. _Ptrunc01S03:

Accessing individual elements
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following routines allow to read/write the elements of a contour function :math:`C(t_0-i\Delta t,t_0-i\Delta t-j\Delta t)` stored as ``cntr::herm_matrix_moving<T>`` at individual time arguments from/to another variable ``M``. All time arguments of the element access functions are understood relative to the leading timestep :math:`t_0` of the window.

The following member functions write components of a contour function :math:`C` from ``M``:

.. list-table::
   :header-rows: 0

   * - ``C.set_les(i,j,M)``
     - :math:`C^<(t_0-i\Delta t,t_0-i\Delta t-j\Delta t)` is set to ``M``
     - required: ``0 <= i,j <= C.tc()``
   * - ``C.set_ret(i,j,M)``
     - :math:`C^R(t_0-i\Delta t,t_0-i\Delta t-j\Delta t)` is set to ``M``
     - required: ``0 <= i,j <= C.tc()``
   * - ``C.get_les(i,j,M)``
     - ``M`` is set to :math:`C^<(t_0-i\Delta t,t_0-i\Delta t-j\Delta t)`
     - required: ``0 <= i,j <= C.tc()``
   * - ``C.get_ret(i,j,M)``
     - ``M`` is set to :math:`C^R(t_0-i\Delta t,t_0-i\Delta t-j\Delta t)`
     - required: ``0 <= i,j <= C.tc()``
   * - ``C.get_gtr(i,j,M)``
     - ``M`` is set to :math:`C^>(t_0-i\Delta t,t_0-i\Delta t-j\Delta t)`
     - required: ``0 <= i,j <= C.tc()``

- If ``C.size1()>1``, ``M`` must be a complex eigen matrix (``cdmatrix`` for double precision)
- If ``C.size1()==1``, ``M`` can be a scalar (``cdouble``) or a matrix
- In the ``set_[X]`` routines, matrices must have the correct dimension ``size1``; in the ``get_[X]`` routines, matrices are resized to dimension ``size1``.

.. _Ptrunc01S04:

Density matrix
~~~~~~~~~~~~~~

The following member function of a ``cntr::herm_matrix_moving<T>`` object ``C`` writes the Density matrix to an object ``M``:

.. list-table::
   :header-rows: 0

   * - ``C.density_matrix(int i,M)``
     - ``M`` is set to :math:`\text{DensityMatrix}[C]({t_0-i\Delta t})`
   * - ``C.density_matrix(int i)``
     - returns the ``(0,0)`` component of the :math:`\text{DensityMatrix}[C]({t_0-i\Delta t})` as a ``cplx`` number

- ``0 <=i <= C.tc()`` is required.
- If ``C.size1()>1``, ``M`` must be a complex eigen matrix (``cdmatrix`` for double precision)
- If ``C.size1()==1``, ``M`` can be a scalar (``cdouble``) or a matrix
- If ``M`` is a matrix, it is resized to a square matrix of dimension ``C.size1()``.

.. _Ptrunc05S04:

Forward move of the truncated time window
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The member function ``C.forward()`` of a ``cntr::herm_matrix_moving`` object is used to shift the time-window forward by one step, such that the new leading timestep corresponds to a physical time :math:`t_0+\Delta t`. This shift is not done by copying data, but by re-assigning pointers:

.. list-table::
   :header-rows: 0

   * - ``C.forward(void)``
     - Transfers the pointer to the timestep ``i`` to ``(i+1)`` for :math:`0\leq {\tt i} \leq {\tt tc}-1` and the pointer to timestep ``tc`` to ``0``. (All timesteps are understood relative to the leading timestep).

If ``C`` stores a window :math:`\mathcal{M}^{\tt tc}_{n}[C]` with a leading time :math:`t_0=n\Delta t`, then after ``C.forward()``, ``C`` stores the window :math:`\mathcal{M}^{\tt tc}_{n+1}[C]`, with the exception of the data at the leading timestep ``0`` of ``C`` (physical timestep :math:`(n+1)\Delta t`). After a move forward of the time window, the data at the new leading timestep will usually be recalculated (timestepping), hence the values after the shift at this time-step are not important.

**Examples:**

.. code-block:: cpp

   cntr::herm_matrix_moving<double> C(tc,size1,sig);

   // initialize C (with a leading physical time t0)
   // ...

   cdmatrix M1,M2;
   int j=3,i=4;  
   C.get_les(i,j,M1);   // M1 set to C^les(t0-i,t0-i-j)
   C.forward();
   C.get_les(i+1,j,M2);   // M2 now equals M1

   C.get_les(tc,j,M1);   // M2 set to C^les(t0-tc,t0-tc-j), a value at the last timeslice of C
   C.forward();
   C.get_les(0,j,M2);
   // M1 will now equal M2:  the value at the new leading timestep
   // is initialized with the value at the previous last step 

   // A typical operation after C.forward() is to replace the new leading
   // timestep 0 by an extrapolation from timesteps i=1,2,3,...:
   C.forward();
   int ExtrapolationOrder=MAX_SOLVE_ORDER; // (tc=>ExtrapolationOrder required).
   cntr::extrapolate_timestep(C,ExtrapolationOrder);

.. _Ptrunc05S04_B:

Initialization from ``cntr::herm_matrix``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following routines are used to transfer data from a ``cntr::herm_matrix``, ``cntr::herm_matrix_timestep`` or ``cntr::herm_matrix_timestep_view`` object to a ``cntr::herm_matrix_moving`` object ``C``.

.. list-table::
   :header-rows: 0

   * - ``C.set_timestep(int i,herm_matrix<T> &G,herm_matrix<T> &Gcc,int tstp)``
     - Set the data at timestep ``i`` of the moving window ``C`` from the physical timestep ``tstp`` of the Green's function ``G``; ``0<=tstp<=G.nt()`` and ``0<=i<=C.tc()`` required.
   * - ``C.set_timestep(int i,herm_matrix_timestep[_view]<T> &G,herm_matrix_timestep[_view]<T> &Gcc)``
     - Set the data at timestep ``i`` of the moving window from the Green's function timeslice ``G``

- The routine sets :math:`C^X(t_0-i\Delta t,t_0-(i+j)\Delta t)=G^X({\tt tstp}\Delta t,({\tt tstp}-j)\Delta t)`, for ``X=R,<`` and ``j=0,...,tc``, where ``tstp`` is either the given argument, or the member ``G.tstp_`` if ``G`` is of type ``herm_matrix_timestep[_view]``
- If ``tstp<tc``, only the values with ``tstp-j >=0`` are set, while the values with ``tstp-j<0`` are set to zero.
- ``G`` and ``Gcc`` store the hermitian domain of :math:`G` and :math:`G^\ddagger`. If :math:`G` is hermitian, call the routine with ``Gcc=G``.

The following routine is the most common way to initialize a moving window to start of a truncated time evolution:

.. list-table::
   :header-rows: 0

   * - ``C_trunc.set_from_G_backward(herm_matrix<T>& G,herm_matrix<T>& Gcc,int tstp)``
     - Equivalent to ``C.set_timestep(i,G,Gcc,tstp-i)`` for ``i=0,...,tc``; ``tstp <= G.nt()`` and ``tstp>=tc`` required

**Example:**

.. code-block:: cpp

   // Typical start of a time stepping with truncated window tc:
   int size1=1,tc=10;
   // first determine the full G up to nt=tc
   int nt=tc;
   int ntau=100;
   cntr::herm_matrix<double> G(nt,ntau,size1,FERMION);
   // do some initialization/calculation of G up to time nt ...

   // store the moving time window of size tc=nt with leading
   // physical time tstp=nt in a GTRUNC variable:
   int tstp=nt;
   cntr::herm_matrix_moving<double> G_trunc(tc,size1,FERMION);
   G_trunc.set_from_G_backward(G,G,tstp); // hermitian symmetry of G is assumed
   // equivalent to 
   G_trunc.clear(); //all elements set to 0
   for(int i=0;i<=tc;i++){
      for(int j=0;j<=tc;j++){
         cdmatrix M;      
         if(tstp-i-j>=0){
           cntr::get_les(tstp-i,tstp-i-j,M,G,G); // M = G^<(tstp-i,tstp-i-j)
           G_trunc.set_les(i,j,M);
           cntr::get_ret(tstp-i,tstp-i-j,M,G,G);
           G_trunc.set_ret(i,j,M);
         }
      }
   }

.. _Ptrunc01S05:

Text file I/O
~~~~~~~~~~~~~

The member functions ``print_to_file`` and ``read_to_file`` of a ``herm_matrix_moving`` implement text-file access.

.. list-table::
   :header-rows: 0

   * - ``C.print_to_file(const char *filename,int precision=16)``
     - Create a text ``filename`` and write content of ``C`` to the file
   * - ``C.read_from_file(const char *filename)``
     - Read content of ``C`` from a file previously written by ``print_to_file``. ``C`` is automatically resized to match parameters ``tc`` and ``size1``.

**Example:**

.. code-block:: cpp

   int tc=10;
   int size1=1;
   int sig=FERMION;
   GTRUNC A(tc,size1,sig);
   // ... set elements of  A ...

   //WRITE:
   // create a file filename.txt and store the data of A:
   A.print_to_file("filename.txt");

   // READ
   GTRUNC B;
   // if the file "filename.txt" has been written previously with print_to_file ,
   // the parameters (tc,size1,sig) and the data of B are modified
   // according to the information in the file:
   B.read_from_file("filename.txt");
   // now B is resized to B.tc()=10,B.size1()=1,B.sig()=FERMION, and the data of B match the data of A.

The file format is rather explicit (and storage intensive):

- First line: ``# tc size1 sig``
- The following lines list the retarded and lesser component in the format:
  ``XXX: i j Re_C^XXX(i,j)_{0,0} Im_C^XXX(i,j)_{0,0} ... Re_C^XXX(i,j)_{size1,size1} Im_C^XXX(i,j)_{size1,size1}``
  where ``XXX`` is ``ret`` for :math:`C^R`, ``les`` for :math:`C^<`, and ``i,j`` loop through the corresponding time arguments of moving time window.

.. _Ptrunc01S05HDF5:

HDF5 file I/O
~~~~~~~~~~~~~

Writing to text files is supposed as a quick way to generate human-readable data. For compressed storage, one should use the HDF5 format below. General info on the HDF5 format, as well as some scripts to interpret the HDF5 files are discussed in separate sections :ref:`Ptrunc13` and :ref:`Ptrunc14` below. HDF5 I/O is accessible if ``libcntr`` is compiled with the flag ``hdf5=ON``.

.. list-table::
   :header-rows: 0

   * - ``C.write_to_hdf5(const char *filename, const char *groupname)``
     - Creates an HDF5 file (or overwrites an existing one) with name ``filename``, creates a data group with name ``groupname``, and writes ``C`` into that group.
   * - ``C.write_to_hdf5(hid_t group_id, const char *groupname)``
     - Creates a data group with handle ``groupname`` inside the existing group with handle ``group_id``, and writes ``C`` into that group ``groupname``.
   * - ``C.write_to_hdf5(hid_t group_id)``
     - Writes ``C`` into a data group with handle ``group_id``.
   * - ``C.read_from_hdf5(const char *filename, const char *groupname)``
     - Open a file ``filename``, look for a data group with name ``groupname``, and read ``C`` from that group. ``C`` is automatically resized to match parameters ``tc`` and ``size1``.
   * - ``C.read_from_hdf5(hid_t group_id, const char *groupname)``
     - Look for a data group with name ``groupname`` in the data group with handle ``group_id``, and read ``C`` from that group. ``C`` is automatically resized to match parameters ``tc`` and ``size1``.
   * - ``C.read_from_hdf5(hid_t group_id)``
     - Read ``C`` from a data group with handle ``group_id``. ``C`` is automatically resized to match parameters ``tc`` and ``size1``.

**Examples:**

.. code-block:: cpp

   #include <cntr/cntr.hpp>
   ..
   // Create a truncated  contour Green's function
   int tc = 200, norb = 1;
   GTRUNC A(tc, norb, FERMION);

   // Open HDF5 file and write components of the Green's function A into a group g.
   std::string filename = "data.h5";
   A.write_to_hdf5(filename.c_str(), "g");

To understand the structure of the resulting HDF5 file ``data.h5`` we inspect it with the ``h5ls`` command line program:

.. code-block:: sh

   $ h5ls -r data.h5
   /                        Group
   /g                       Group
   /g/element_size          Dataset {1}
   /g/les                   Dataset {2601, 1, 1}
   /g/ret                   Dataset {2601, 1, 1}
   /g/sig                   Dataset {1}
   /g/size1                 Dataset {1}
   /g/size2                 Dataset {1}
   /g/tc                    Dataset {1}

where we see that apart from the contour components the Green's function group ``g`` contains additional information about the dimensions and the Fermi/Bose statistics (``sig`` :math:`= \mp 1`).
Since the cutoff time is :math:`t_c = 50`, the size of the truncated window of the ``/g/ret`` and ``/g/les`` components contains :math:`(t_c + 1)(t_c + 1) = 2601` elements.

If the file ``data.h5`` has been written previously with ``write_to_hdf5``, one can read it with the member function ``read_from_hdf5``:

.. code-block:: cpp

   // Open HDF5 file and read group g. The result is saved into the Green's function B
   GTRUNC B;
   B.read_from_hdf5(filename.c_str(), "g");
   // the parameters (tc,size1,sig) and the data of B are modified according to
   // the information in the file.

The following example illustrates how to write several Green's functions into one HDF5 file:

.. code-block:: cpp

   int tc = 200, norb = 1;
   GTRUNC A(tc, norb, FERMION);
   GTRUNC B(tc, norb, FERMION);

   // set data of A and B ...
   // create a hdf5 file
   std::string filename = "data.h5";
   hid_t file_id = open_hdf5_file(filename);
   A.write_to_hdf5(file_id, "a"); // file now has group a/ with a/ret/, a/les/, ...
   B.write_to_hdf5(file_id, "b"); // file now has group b/ with a/ret/, a/les/, ...
   close_hdf5_file(file_id);

   // Open HDF5 file for reading using the cntr hdf5 utils:
   file_id = read_hdf5_file(filename.c_str());
   B.read_from_hdf5(file_id, "b");
   A.read_from_hdf5(file_id, "a");
   close_hdf5_file(file_id);

To simplify postprocessing of contour Green's functions, NESSi also provides the python module ``ReadCNTRhdf5.py`` for reading the HDF5 format (using the python modules ``numpy`` and ``h5py``), producing python objects with the contour components as members.

.. _Ptrunc02:

Truncated Timeslices
--------------------

.. _Ptrunc02S01:

Overview
~~~~~~~~

.. list-table::
   :header-rows: 0

   * - class
     - ``cntr::herm_matrix_timestep_moving<T>``

We define the timestep with memory truncation ``tc`` at times ``t0`` as the data:

- :math:`C^\mathrm{R}(t_0,t_0 - j \Delta t)` for ``j=0,...,tc``,
- :math:`C^<(t_0,t_0 - j \Delta t)` for ``j=0,...,tc``.

The data of the ``herm_matrix_moving`` Green's function :math:`C` are therefore a union of the timesteps at times :math:`t_0, t_0-\Delta t, ..., t_0-t_c \Delta t`. The class ``cntr::herm_matrix_timestep_moving<T>`` is the container to store a truncated timeslice. It is characterized by the following parameters:

- ``T`` (template parameter): Precision, usually set to ``double``; we use the definition

  .. code-block:: cpp

     #define GTRUNC_TSTP cntr::herm_matrix_timestep_moving<double>

- ``tc`` (integer): number of discretization points on the real time axis.
- ``size1`` (integer): orbital dimension. Each element :math:`C(t,t')` is a square matrix of dimension ``size1`` :math:`\times` ``size1``.
- ``sig`` (``FERMION`` or ``BOSON``).

.. _Ptrunc02S02:

Constructors
~~~~~~~~~~~~

.. list-table::
   :header-rows: 0

   * - ``herm_matrix_timestep_moving<T>()``
     - Default constructor, does not allocate memory and sets ``tc=-1``.
   * - ``herm_matrix_timestep_moving<T>(int tc,int size1,int sig)``
     - Allocate memory, set all entries to ``0``. It requires ``tc>=-1``, ``size1>0``, and ``sig=FERMION`` or ``sig=BOSON``.
   * - ``herm_matrix_timestep_moving<T>(herm_matrix_timestep_moving<T> &g)``
     - Copy constructor. Create a new instance with new data, identical to ``g``.
   * - ``herm_matrix_timestep_moving<T>(herm_matrix_moving<T> &g,int i)``
     - Set ``tc``\ =\ ``g.tc()``, ``sig``\ =\ ``g.sig()``, ``size1``\ =\ ``g.size1()``. If ``tc=-1`` no data are transferred; if ``tc>=0`` the data from timeslice ``i`` of ``g`` are copied to the ``herm_matrix_timestep_moving``.

.. _Ptrunc02S03:

Accessing individual matrix elements
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Element access, as well as read-out of the density matrix follows a similar syntax as for the Green's function ``herm_matrix_moving<T>``. All time arguments of the element access functions are understood relative to the leading timestep :math:`t_0`.

.. list-table::
   :header-rows: 0

   * - ``C.set_les(j,M)``
     - :math:`C^<(t_0,t_0-j\Delta t)` is set to ``M``
     - required: ``0 <= j <= C.tc()``
   * - ``C.set_ret(j,M)``
     - :math:`C^R(t_0,t_0-j\Delta t)` is set to ``M``
     - required: ``0 <= j <= C.tc()``
   * - ``C.get_les(j,M)``
     - ``M`` is set to :math:`C^<(t_0,t_0-j\Delta t)`
     - required: ``0 <= j <= C.tc()``
   * - ``C.get_ret(j,M)``
     - ``M`` is set to :math:`C^R(t_0,t_0-j\Delta t)`
     - required: ``0 <= j <= C.tc()``
   * - ``C.get_gtr(j,M)``
     - ``M`` is set to :math:`C^>(t_0,t_0-j\Delta t)`
     - required: ``0 <= j <= C.tc()``

**Density matrix:**

.. list-table::
   :header-rows: 0

   * - ``C.density_matrix(M)``
     - ``M`` is set to :math:`\text{DensityMatrix}[C](t_0)`
   * - ``C.density_matrix(void)``
     - returns the (0,0) component of the :math:`\text{DensityMatrix}[C]({t_0-i\Delta t})` as a ``cplx`` number

**Example:**

.. code-block:: cpp

   int tc=10;
   int size=2;
   GTRUNC A(tc,size,FERMION); // allocate a full window
   // .. set data

   GTRUNC_TSTP tA(A,0); // has create a timestep variable which contains a copy of the data at the leading timestep 0 of A
   cdmatrix M2,M1; 
   tA.get_les(8,M1);  // ok: M1 = A^<(t0,t0-8), where t0 is the timeslice on which A has been created
   A.get_les(0,8,M2);  // now M1=M2

   //print of all elements of tA:
   for(int j=0;j<=tc;j++){
      cdouble x;
      double y;
      tA.get_ret(j,x);
      tA.get_les(j,y);
      cout << "Aret(t0,t0-j) at j=" << j << " : "  << x << endl;
      cout << "Ales(t0,t0-j) at j=" << j << " : "  << y << endl;
   }

.. _Ptrunc03:

Truncated Map
-------------

.. _Ptrunc03S01:

Overview
~~~~~~~~

.. list-table::
   :header-rows: 0

   * - class
     - ``cntr::herm_matrix_timestep_moving_view<T>``

The purpose of this class is to create a map to a pre-existing object of type ``cntr::herm_matrix_moving`` or ``cntr::herm_matrix_timestep_moving``. One can then use this class without making a physical copy of the original data. The structure is inherited from the existing object. It is characterized by the following parameters:

- ``T`` (template parameter): Precision, usually set to ``double``; we use the definition

  .. code-block:: cpp

     #define GTRUNC_TSTP_VIEW cntr::herm_matrix_timestep_view<double>

- ``tc`` (integer): number of discretization points on the real time axis.
- ``size1`` (integer): orbital dimension. Each element :math:`C(t,t')` is a square matrix of dimension ``size1`` :math:`\times` ``size1``.
- ``sig`` (``FERMION`` or ``BOSON``).

.. note::

   All classes use raw pointers. If a map ``G_view`` is created on an object ``G``, then the pointers in ``G_view`` become invalid if ``G`` moves out of scope. As for ``herm_matrix_timestep_view``, ``herm_matrix_timestep_view_moving`` should therefore be used with care to avoid dangling pointers.

.. _Ptrunc03S02:

Constructors
~~~~~~~~~~~~

.. list-table::
   :header-rows: 0

   * - ``herm_matrix_timestep_moving_view<T>()``
     - Default constructor. Data pointers remain invalid (zero).
   * - ``herm_matrix_timestep_moving_view(herm_matrix_moving<T> &g,int i)``
     - Creates a map to timeslice ``i`` of the ``herm_matrix_moving`` object.
   * - ``herm_matrix_timestep_moving_view(herm_matrix_moving<T> &g)``
     - Creates a map to the ``herm_matrix_timestep_moving`` object.

.. _Ptrunc03S03:

Accessing and manipulation
~~~~~~~~~~~~~~~~~~~~~~~~~~~

All functions described for ``herm_matrix_timestep_moving`` in :ref:`Ptrunc02` can be replaced by a ``herm_matrix_timestep_moving_view`` argument following the same syntax.

**Example:**

.. code-block:: cpp

   cdmatrix M,M1;
   herm_matrix_timestep_view<double> gtmp;
   {
     int tc=10;
     int size=1;
     GTRUNC G(tc,size,FERMION);
     // ... do something with G ...
     int i=5;
     gtmp=herm_matrix_timestep_view<double>(G,i);
     int j=3;
     G.get_les(i,j,M);  
     gtmp.get_les(j,M1); // now M=M1
   }
   // after G is out of scope, calls to the data of gtmp cause a memory error:
   gtmp.get_les(j,M1); // invalid

.. _Ptrunc05:

Truncated Contour Functions
----------------------------

.. _Ptrunc05S01:

Overview
~~~~~~~~

.. list-table::
   :header-rows: 0

   * - class
     - ``cntr::function_moving<T>``

A contour function is a matrix or scalar-valued function :math:`f(t)` which depends on physical time only. On the real-time branches, its value is the same for the same time :math:`t` on the upper and lower contour branch. A contour function is typically used to store time-dependent parameters of a system, or a time-dependent Hamiltonian.

The class ``cntr::function_moving<T>`` is the container to store these data. It is characterized by the following parameters:

- ``T`` (template parameter): Precision, usually set to ``double``; we use the definition

  .. code-block:: cpp

     #define CTRUNC cntr::function_moving<double>

- ``tc`` (integer): number of discretization points on the real time axis.
- ``size1`` (integer): orbital dimension. Each element :math:`f(t)` is a square matrix of dimension ``size1`` :math:`\times` ``size1``.

.. _Ptrunc05S02:

Constructors
~~~~~~~~~~~~

.. list-table::
   :header-rows: 0

   * - ``function_moving<T>()``
     - Default constructor, does not allocate memory and sets ``tc=-1``.
   * - ``function_moving<T>(int tc,int size1)``
     - Allocate memory, set all entries to ``0``. It requires ``tc>=-1``, ``size1>0``.

.. _Ptrunc05S03:

Accessing individual matrix elements
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following routines allow to read/write the elements of a contour function :math:`f(t)` stored as ``cntr::function_moving`` at individual time arguments from/to another variable ``M``. The latter can be either a scalar of type ``std::complex<T>``, or a complex square matrix defined in Eigen.

The following member functions of ``cntr::function_moving<T>`` set components of a contour function :math:`f(t)` from ``M``:

.. list-table::
   :header-rows: 0

   * - ``f.set_value(int i,M)``
     - For ``0<=i<=f.tc()``, :math:`f(t_0-i\Delta t)` is set to ``M``

- If ``f.size1()>1``, ``M`` must be a square matrix (``cdmatrix`` for ``double`` precision)
- If ``f.size1()==1``, ``M`` can be also a scalar (``cdouble``) or a square matrix
- If ``M`` is a matrix, it must be a square matrix of dimension ``f.size1()``

The following member functions of ``cntr::function_moving<T>`` read components of a contour function :math:`f` to ``M``:

.. list-table::
   :header-rows: 0

   * - ``f.get_value(i,M)``
     - For ``0<=i<=f.tc()``, ``M`` is set to :math:`f(t_0-i\Delta t)`

- If ``M`` is a matrix, it is resized to a square matrix of dimension ``f.size1()``
- If ``M`` is a scalar, only the (0,0) entry of :math:`f(t)` is read.

.. _Ptrunc05S04_fwd:

Manipulating the truncated time window
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following function shuffles the pointer to the individual timesteps:

.. list-table::
   :header-rows: 0

   * - ``f.forward(void)``
     - Transfers the pointer to the timestep ``t0-i`` to ``t0-(i+1)`` and ``t0-tc`` to ``t0``.

The following routine is used to transfer data from a ``cntr::function<T>`` to a ``cntr::function_moving<T>`` object:

.. list-table::
   :header-rows: 0

   * - ``f.set_from_function_backward(function<T>& f,int tstp)``
     - For ``0<=tstp<=f.nt()`` sets the data of ``function_moving`` to the data of ``function`` from ``tstp`` to ``tstp``-``tc``

**Example:**

.. code-block:: cpp

   int nt=100,size1=1;
   cntr::function<double> ft(nt,size1);
   // set data of ft
   int tc=10;
   cntr::function_moving<double> ft_trunc(tc,size1);
   int tstp=40;
   ft_trunc.set_from_function_backward(ft,tstp);
   // equivalent to:
   for(int i=0;i<=tc;i++){
     cdmatrix M;
     ft.get_value(tstp-i,M);
     ft_trunc.set_value(i,M);
   }

.. _Ptrunc02S04:

Timestep-wise data manipulation and access
-------------------------------------------

The following routines are used for transferring data between timeslices of different objects, or doing simple algebra on timeslices. They all work with a similar syntax for the types ``cntr::herm_matrix_moving``, ``cntr::herm_matrix_timestep_moving``, and ``cntr::herm_matrix_timestep_moving_view``, and all types can be mixed.

.. _Ptrunc02S04_A:

Copying data between timeslices
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In typical applications, data between different Green's functions are exchanged at once for an entire timestep. We provide the member functions which allow the exchange of data of one entire timestep between variables of type ``cntr::herm_matrix_moving``, ``cntr::herm_matrix_timestep_moving``, and ``cntr::herm_matrix_timestep_moving_view``, with a common syntax (all time variables of ``_moving`` objects are understood relative to the leading timestep of the object, which is denoted by ``t0`` below):

.. list-table::
   :header-rows: 0

   * - ``A.set_timestep( Timeslice_A , B ,  Timeslice_B )``
     - Copy the data of ``Timeslice_B`` of ``B`` to ``Timeslice_A`` of ``A``

- ``A`` and ``B`` can be ``herm_matrix_timestep_moving_view``, ``herm_matrix_timestep_moving``, or ``herm_matrix_moving``. If ``X``\ =\ ``A,B`` is of type ``herm_matrix_timestep_moving`` or ``herm_matrix_timestep_moving_view``, the argument ``Timeslice_X`` must be omitted.
- ``B.tc()==A.tc()`` is required
- ``0 <= timeslice_X <= X.tc()`` for ``X``\ =\ ``A,B`` required.

**Example:**

The following code copies the 5 leading timesteps from a Green's function ``B`` to a Green's function ``A``:

.. code-block:: cpp

   int tc=10;
   int size1=2;
   GTRUNC A(tc,size1,FERMION);
   GTRUNC B(tc,size1,FERMION);
   // ... do something with B ...
   for(int i=0;i<=5;i++) A.set_timestep(i,B,i);

   // equivalent to:
   for(int i=0;i<=5;i++){
     for(int j=0;j<=tc;j++){
       cdmatrix M;
       B.get_les(i,j,M);
       A.set_les(i,j,M);
       B.get_ret(i,j,M);
       A.set_ret(i,j,M);
     }
   }

**Set only selected matrix elements:**

The routines ``set_matrixelement`` are similar to ``set_timestep`` but do the data exchange only between selected matrix elements:

.. list-table::
   :header-rows: 0

   * - ``A.set_matrixelement( Timeslice_A, int i1,int i2, B ,  Timeslice_B,int j1,int j2)``
     - Set matrix-element ``(i1,i2)`` at ``Timeslice_A`` of ``A`` to matrix-element ``(j1,j2)`` of ``Timeslice_B`` of ``B``.

.. _Ptrunc06S01:

Scalar multiplication and incrementation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 0

   * - ``A.incr_timeslice(Timeslice_A, B , Timeslice_B,scalar alpha=1.0)``
     - ``A(..)`` at timeslice ``Timeslice_A`` incremented by ``alpha*B(..)`` at ``Timeslice B`` of B, with a scalar (complex or real) alpha. Default ``alpha=1.0``.
   * - ``A.set_timeslice_zero(Timeslice_A)``
     - ``A(..)`` at timeslice ``Timeslice_A`` set to zero
   * - ``A.smul(Timeslice_A,scalar alpha)``
     - ``A(..)`` at timeslice ``Timeslice_A`` multiplied with scalar (complex or real) alpha

**Example:**

.. code-block:: cpp

   int tc=10;
   int size1;
   GTRUNC A(tc,size1,FERMION);
   GTRUNC B(tc,size1,FERMION);
   // ... do something with B ...

   // increment leading timestep of B to leading timestep of A
   A.incr_timestep(0,B,0);
   // equivalent to:
   cntr::herm_matrix_timestep_moving_view B_view(B,0);
   A.incr_timestep(0,Bview);
   // also equivalent to:
   for(int j=0;j<=tc;j++){
     cdmatrix M,M1;
     B.get_les(0,j,M);  
     A.get_les(0,j,M1);
     M1+=M;
     A.set_les(0,j,M1);
     // ... analogous for ret
     B.get_ret(0,j,M);  
     A.get_ret(0,j,M1);
     M1+=M;
     A.set_ret(0,j,M1);
   }

**Example: Fourier transform**

The following operation considers a set of Green's functions :math:`G_k` indexed by momentum :math:`k`, and computes the Fourier transform at a given ``x`` at the leading timestep:

.. math::

   G_{j}(t,t') = \frac{1}{N_k}\sum_{k=0}^{N_k-1} e^{i2\pi k x /N_k} G_k(t,t').

.. code-block:: cpp

   int tc=10;
   int t0=10;
   int size1=1;
   int nk=10;
   std::vector<GTRUNC> Gk(nk); // a vector of nk=10 Green's functions
   for(int k=0;k<nk;k++) Gk[k]=GTRUNC(tc,size1,FERMION); // allocation
   // ...
   // ... do something with Gk
   // ...
   // do Fourier transform on all timesteps 
   // and store result in timestep variable tGav:
   GTRUNC_TSTP tGav(tc,size1,FERMION); // note: tGav is initialized with 0 already
   double nkinv=1.0/nk;
   double x=2.0;
   for(int k=0;k<nk;k++){
     cdouble weight=exp(II*2.0*M_PI*x*k*nkinv);
     tGav.incr_timestep(Gk[k],0,weight); // arg '0' indicates leading timestep
   }
   tGav.smul(nkinv);

.. _Ptrunc06S02:

Multiplication with one-time contour functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 0

   * - ``A.left_multiply(cntr::function_moving<T> &f)``
     - :math:`A(t,t')` set to :math:`f(t) A(t,t')` **at leading timestep of** ``A``
   * - ``A.right_multiply(cntr::function_moving<T> &f)``
     - :math:`A(t,t')` set to :math:`A(t,t') f(t')` **at leading timestep of** ``A``
   * - ``A.left_multiply_hermconj(cntr::function_moving<T> &f)``
     - :math:`A(t,t')` set to :math:`f(t)^\dagger A(t,t')` **at leading timestep of** ``A``
   * - ``A.right_multiply_hermconj(cntr::function_moving<T> &f)``
     - :math:`A(t,t')` set to :math:`A(t,t') f(t')^\dagger` **at leading timestep of** ``A``

- ``A`` is of type ``cntr::herm_matrix_moving<T>``, ``cntr::herm_matrix_timestep_moving<T>``, or ``cntr::herm_matrix_timestep_moving_view<T>``
- ``A.tc()<=f.tc()`` and ``f.size1()==A.size1()`` is required.

**Example:**

The following example computes a retarded interaction :math:`W(t,t') = U(t) \chi(t,t') U(t')` out of a given susceptibility :math:`\chi`, where :math:`U(t)` is a time-dependent bare interaction matrix element.

.. code-block:: cpp

   int tc=10;
   int size1=1;
   cntr::herm_matrix_moving<double> W(tc,size1,BOSON);
   cntr::herm_matrix_moving<double> chi(tc,size1,BOSON);
   cntr::function_moving<double> U(tc,size1);
   // ...
   // ... do something to set U and  chi: both are initialized with same leading timestep t0
   // ...
   // compute W on the leading timestep
   W.set_timestep(0,chi,0);
   W.left_multiply(U);
   W.right_multiply(U);

.. _Ptrunc11S01:

Comparing Green's functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To compare the data of two Green's functions on a given timestep, one can use:

.. list-table::
   :header-rows: 0

   * - ``T cntr::distance_norm2(GType &A, Timestep_A, GType &B, Timestep_A)``
     - Returns a difference measure :math:`\Delta[A,B]` between ``Timestep_A`` of ``A`` and ``Timestep_B`` of ``B``

The difference is defined as the :math:`L_2`-norm difference :math:`|| M ||` for :math:`M=A(t,t')-B(t,t')` of the individual elements, summed over all time-arguments of the timestep, as well over retarded and lesser components:

.. math::
   :label: eq:distnorm_trunc

   \Delta[A,B]
   =
   \sum_{j=0}^{\tt tc}
   \big(
   ||A^\mathrm{R}(t_0,t_0-j\Delta t)-B^\mathrm{R}(t_0,t_0-j\Delta t)||
   +
   ||A^<(t_0,t_0-j\Delta t)-B^<(t_0,t_0-j\Delta t)||
   \big)

**Example:**

.. code-block:: cpp

   GTRUNC G(tc,size1,sig);
   // some iterative procedure to determine G on timestep t0:
   {
       GTRUNC_TSTP tG(tc,size1,sig);   //temporary variable
       double convergence_error;
       int iter_max=100;
       for(int iter=0;iter<=iter_max;iter++){
           tG.set_timestep(G,0);          // store values of G before iteration
           // ... some code to update G on timestep tstp ...
           convergence_error = cntr::distance_norm2(G,0,tG);
           if( convergence_error < some_sufficiently_small_number) break;
       }
       if(iter>iter_max){
           cout << "no convergence!" << endl;
       }
   }

.. _Ptrunc11S02:

Timestep extrapolation
~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 0

   * - ``void cntr::extrapolate_timestep(GTRUNC &A,int ExtrapolationOrder)``
     - Extrapolate from timesteps ``t=t0-ExtrapolationOrder,...,t0-1`` to timestep ``t0``, using polynomial extrapolation.

- Size requirements: ``t0>=ExtrapolationOrder``
- The ``ExtrapolationOrder`` must be between ``1`` and ``MAX_SOLVE_ORDER`` (``=5``).

.. _Ptrunc07:

Truncated diagram utilities
----------------------------

The basic building blocks of Feynman diagrams are particle-hole bubbles (``Bubble1``) and particle-particle bubbles (``Bubble2``), as described in the Diagrams section of the manual for the full Green's functions. We provide two basic functions that allow to compute particle-particle and particle-hole Bubbles on a given timestep. The calls are essentially the same as for full Green's functions. **The operation is always performed on the leading timestep of the functions**.

.. list-table::
   :header-rows: 0

   * - ``cntr::Bubble1(GG &C,int c1,int c2,GG &A,GG &Acc,int a1,int a2,GG &B,GG &Bcc,int b1,int b2)``
     - :math:`C_{c1,c2}(t,t')=i A_{a1,a2}(t,t')  B_{b2,b1}(t',t)` is calculated on the leading timestep of ``C`` for given two-time objects ``A``, ``B`` with matrix indices ``a1``, ``a2`` and ``b1``, ``b2``, respectively.
   * - ``cntr::Bubble2(GG &C,int c1,int c2,GG &A,GG &Acc,int a1,int a2,GG &B,GG &Bcc,int b1,int b2)``
     - :math:`C_{c1,c2}(t,t')=i A_{a1,a2}(t,t')  B_{b1,b2}(t,t')` is calculated on the leading timestep of ``C`` for given two-time objects ``A``, ``B`` with matrix indices ``a1``, ``a2`` and ``b1``, ``b2``, respectively.

- ``C``, ``A``, ``Acc``, ``B``, ``Bcc`` can be of type ``cntr::herm_matrix_moving``, ``cntr::herm_matrix_timestep_moving`` or ``cntr::herm_matrix_timestep_moving_view``
- For ``X``\ =\ ``A`` or ``B``: ``Xcc`` contains the hermitian conjugate :math:`X^\ddagger` of :math:`X`. If the argument ``Xcc`` is omitted, it is assumed that ``X`` is hermitian, :math:`X=X^\ddagger`.
- If the index pairs ``(c1,c2)``, ``(a1,a2)``, or ``(b1,b2)`` are omitted, they are assumed to be ``(0,0)``

**Example:**

The following example calculates the second-order self-energy for a general time-dependent interaction :math:`U(t)` out of a scalar local Green's function, :math:`\Sigma(t,t') = U(t) G(t,t') G(t',t) U(t') G(t,t')`. The diagram is factorized in a particle-hole bubble :math:`\chi` and a particle-particle bubble:

.. math::

   \Sigma(t,t') = -i W(t,t') G(t,t'),
   \quad
   W(t,t')=U(t) \chi(t,t') U(t'),
   \quad
   \chi(t,t') = iG(t,t') G(t',t).

.. code-block:: cpp

   int tc=10;
   int size1=1;
   GTRUNC G(tc,size1,FERMION);
   GTRUNC Sigma(tc,size1,FERMION);
   CTRUNC U(tc,size1);
   // ...
   // ... do something to set G and U
   // ...
   // compute Sigma on the leading timestep:
   GTRUNC_TSTP tW(tc,size1,BOSON); // used as temporary variable; Note that a bubble of two Fermion GF is bosonic
   cntr::Bubble1(tW,G,G); // W(t,t') is set to chi=ii*G(t,t')G(t',t);
   // set tW  to W(t,t')=U(t)chi(t,t')U(t'):
   tW.left_multiply(U);
   tW.right_multiply(U); 
   cntr::Bubble2(Sigma,G,W); // Sigma(t,t') is set to Sigma=ii*G(t,t')W(t',t);
   Sigma.smul(0,-1.0); // the final -1 sign on the leading timestep 0.

.. _Ptrunc09:

Truncated Dyson Equation
-------------------------

The Dyson equation for the Green's function :math:`G(t,t')` can be written as

.. math::
   :label: trunc_dyson

   i\partial_t G(t,t^\prime) + \mu G(t,t^\prime)  - \epsilon(t) G(t,t^\prime) -
   \int_\mathcal{C} d\bar t\, \Sigma(t,\bar t) G(\bar t,t^\prime) = \delta_{\mathcal{C}}(t,t^\prime).

This equation is to be solved for :math:`G(t,t^\prime)` for given input :math:`\epsilon(t)` and :math:`\Sigma(t,t^\prime)`, and the KMS boundary conditions. All quantities :math:`\Sigma(t,t^\prime)`, :math:`G(t,t^\prime)`, and :math:`\epsilon(t)` are square matrices. The equation is an integro-differential form of the Dyson series

.. math::

   G=G_0+G_0\ast \Sigma\ast G,

where the free Green's function :math:`G_0` is determined by the differential equation :math:`i\partial_t G(t,t^\prime) + \mu G(t,t^\prime) - \epsilon(t) G(t,t^\prime) = \delta_{\mathcal{C}}(t,t^\prime)`.

It is assumed that :math:`\Sigma=\Sigma^\ddagger` is hermitian, and :math:`\epsilon(t) =\epsilon(t)^\dagger`, which
implies that also the solution :math:`G` possesses hermitian symmetry. Because of the hermitian symmetry, :math:`G` can also be determined from the equivalent conjugate equation

.. math::

   -i\partial_{t^\prime} G(t,t^\prime) + \mu G(t,t^\prime)-   G(t,t^\prime) \epsilon(t^\prime)-
   \int_\mathcal{C} d\bar t \,G(t,\bar t) \Sigma(\bar t,t^\prime) = \delta_{\mathcal{C}}(t,t^\prime).

Because of its causal nature, the Dyson equation can be solved in a time-stepping manner. The following routines are used to solve the Dyson equation:

.. list-table::
   :header-rows: 0

   * - ``void cntr::dyson_timestep(GG &G, GG &Sigma, cntr::function_moving<T> &H, T mu, int SolveOrder, T dt)``
     - Solve :eq:`trunc_dyson` for ``G`` on timestep ``G.t0()==Sigma.t0()`` with ``t0 > SolveOrder`` with given ``Sigma``, ``mu``, ``H``

- ``G`` and ``Sigma`` are Green's functions of type ``cntr::herm_matrix_moving<T>``, assumed to be hermitian
- ``H`` is a ``cntr::function_moving<T>``, the function :math:`\epsilon` in :eq:`trunc_dyson`.
- ``SolveOrder`` :math:`\in` ``1,...,MAX_SOLVE_ORDER``, the order of accuracy for the solution. Use ``SolveOrder=5`` (``=MAX_SOLVE_ORDER``) if there is no good reason against it.
- ``dt`` is the time-discretization step :math:`\Delta t` on the real-time branch
- Size requirements: ``G.size1()==Sigma.size1()==H.size1()``, ``G.t0()==Sigma.t0()``, ``G.tc()==Sigma.tc()``, ``G.t0() > SolveOrder``

**Input/Output relation and time-stepping:**

- ``dyson_timestep``: ``Sigma`` is read on timestep ``t=t0-tc,...,t0``, ``H`` is read on times ``t=t0-tc,...,t0``, ``G`` is read on timestep ``t=t0-tc,...t0-1``. ``G`` is written on timestep ``t=t0``.

Because of this causal structure, the Dyson equation can be solved by time-stepping: The equation is solved successively for timesteps ``tstp=SolveOrder+1,SolveOrder+2,...``, where the result at timestep ``tstp`` depends on all previous timesteps.

**Example:**

Typical solution of :eq:`trunc_dyson` with a non-self-consistent :math:`\Sigma`

.. code-block:: cpp

   int tc=10;
   int size1=1;
   int SolveOrder=5;
   double dt=0.01; // time-discretiation
   double mu=0;
   GTRUNC G(tc,size1,FERMION);
   GTRUNC S(tc,size1,FERMION);
   CTRUNC H(tc,size1);
   // ... initialize G, Sigma, H 

   for (int tstp=tc+1;tstp<=nt;tstp++){
       G.forward();
       S.forward();
       H.forward();
       H.set_value(0,...); // set H at leading timestep
       // compute Sigma at its new leading timestep.
       cntr::dyson_timestep(G,S,H,mu,SolveOrder,dt); // solve dyson
   }

See :ref:`Trunc` for an example of the complete time propagation scheme.

.. _Ptrunc10:

Truncated VIE2
--------------

The second important equation is an integral equation of the form

.. math::
   :label: trunc_vie2

   G(t,t^\prime) + \int_\mathcal{C} d\bar t\, F(t,\bar t) G(\bar t,t^\prime) = Q(t,t^\prime) \quad \Leftrightarrow\quad (1+F)*G=Q,

This linear equation is to be solved for :math:`G(t,t^\prime)` for a given input kernel :math:`F(t,t^\prime)`, its hermitian conjugate :math:`F^\ddagger(t,t^\prime)`, and a source term :math:`Q(t,t^\prime)`.
A typical physical application of :eq:`trunc_vie2` is the summation of a random phase approximation (RPA) series for a susceptibility :math:`\chi`,

.. _rpa_def:

.. math::
   :label: eq:rpa_vie2_trunc

   \chi
   = \chi_0 + \chi_0\ast V \ast \chi_0 +  \chi_0\ast V \ast \chi_0\ast V \ast \chi_0 + \cdots
   = \chi_0 + \chi_0\ast V \ast \chi.

In the solution of the equation, we assume that both :math:`Q` and :math:`G` are hermitian. In general, the hermitian symmetry would not hold for an arbitrary input :math:`F` and :math:`Q`. However, it does hold when :math:`F` and :math:`Q` satisfy the relation :math:`F\ast Q=Q\ast F^\ddagger` and :math:`Q=Q^\ddagger`, which is the case for the typical applications such as the RPA series. In this case, there is an equivalent conjugate equation

.. math::

   G(t,t^\prime) + \int_\mathcal{C} d\bar t \,G(t,\bar t) F^\ddagger(\bar t,t^\prime) = Q(t,t^\prime) \quad \Leftrightarrow\quad G*(1+F^\ddagger)=Q.

Because of its causal nature, the VIE2 equation can be solved in a time-stepping manner. The following routines are used to solve the VIE2 equation:

.. list-table::
   :header-rows: 0

   * - ``void cntr::vie2_timestep(GG &G, GG &F, GG &Fcc, GG &Q,int SolveOrder, T dt)``
     - Solve :eq:`trunc_vie2` for ``G`` on timestep ``G.t0()`` with given two-time objects ``F`` and ``Q``

- The type GG of ``G``, ``F``, ``Fcc``, and ``Q`` is ``cntr::herm_matrix_moving<T>``; ``F`` and ``Fcc`` store the hermitian domain of :math:`F` and :math:`F^\ddagger`, respectively.
- ``SolveOrder`` :math:`\in` ``1,...,MAX_SOLVE_ORDER``, the order of accuracy for the solution. Use ``SolveOrder=5`` (``=MAX_SOLVE_ORDER``) if there is no good reason against it.
- ``dt`` is the time-discretisation step :math:`\Delta t` on the real-time branch
- Size requirements: ``G``, ``F``, ``Fcc``, and ``Q`` must have the same ``size1`` and ``tc``; for ``vie2_timestep``: ``G``, ``F``, ``Fcc``, and ``Q`` must have ``tc >= kt*2+2``

**Input/Output relation and time-stepping:**

- ``vie2_timestep``: ``F`` and ``Q`` are read on timestep ``t=t0-tc,...,t0``, ``G`` is read on timestep ``t=t0-tc,...t0-1``. ``G`` is written on timestep ``t=t0``.

Because of this causal structure, the VIE2 equation can be solved by time-stepping, analogous to ``dyson``.

**Example:**

Solution of the RPA series :math:`\chi = \chi_0 + \chi_0\ast W\ast \chi` for :math:`\chi`, where :math:`\chi_0` is a known (hermitian) two-time function, and :math:`W` is a known one-time function (also hermitian, :math:`W(t)=W(t)^\dagger`). The equation is cast in the form :eq:`trunc_vie2` with :math:`F (t,t')= -\chi_0(t,t') W(t')`, :math:`F^\ddagger (t,t')= -W(t)\chi_0(t,t')`, and :math:`Q=\chi_0`.

.. code-block:: cpp

   int tc=10;
   int t0=20;
   int size1=2;
   int SolveOrder=5;
   double dt=0.01; // time-discretiation
   GTRUNC chi(tc,t0,size1,BOSON);
   GTRUNC chi0(tc,t0,size1,BOSON);
   CTRUNC W(tc,size1);
   GTRUNC F(tc,t0,size1,BOSON);
   GTRUNC Fcc(tc,t0,size1,BOSON);
   //
   // ... do something to determine chi0 and W on timestep -1
   // determine F, Fcc, and solve for chi on timestep -1:
   tstp=-1;
   F.set_timestep(tstp,chi0);
   F.right_multiply(tstp,W,-1.0);
   Fcc.set_timestep(tstp,chi0);
   Fcc.left_multiply(tstp,W,-1.0);
   cntr::vie2_mat(chi,F,Fcc,chi0,beta,SolveOrder,CNTR_MAT_FIXPOINT);
   // do something to determine chi0 and W  on timesteps t0-tc ... t0
   // get F, Fcc, and solve for chi on timestep t0:
   for(int i;i<=tc;i++){
     F.set_timestep(i,chi0,i);
     F.right_multiply(i,W,-1.0);
     Fcc.set_timestep(i,chi0,i);
     Fcc.left_multiply(i,W,-1.0);
   }
   cntr::vie2_timestep(chi,F,Fcc,chi0,SolveOrder,dt);

.. _Ptrunc12:

MPI parallelization
-------------------

In large-scale applications, one often encounters problems that can be
parallelized by distributing the calculation of different Green's
functions across different computing ranks. A typical example are lattice
simulations, where the solution of the Dyson equation for the Green's
function :math:`G_k` on each point :math:`k` on a discrete momentum grid
can be performed in parallel. The library provides some features for
distributed-memory parallelization via MPI:

- A class ``cntr::distributed_timestep_array_moving`` can be used to handle all-to-all communication of Green's functions on a given timestep.
- The class ``cntr::distributed_timestep_array_moving`` is derived from a class ``distributed_array``, which implements an array of data with distributed ownership.

.. _Ptrunc12S01:

Distributed timestep array
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The distributed timestep array allows an all-to-all communication of a set of Green's functions at one timestep. We provide a class ``cntr::distributed_timestep_array_moving`` which is customized for this application:

.. list-table::
   :header-rows: 0

   * - class
     - ``cntr::distributed_timestep_array_moving<T>``

- The ``cntr::distributed_timestep_array_moving<T>`` contains an array ``T_k`` of ``n`` ``herm_matrix_timestep_moving``s (each with the same ``tc``, ``t0``, ``size1`` and ``sig`` parameters), which are indexed by an index ``k`` :math:`\in\{` ``0,...,n-1`` :math:`\}`. ``k`` can be, e.g., a momentum label for lattice simulations.
- The full array is accessible at each MPI rank, but each timestep ``T_k`` is owned by precisely one rank (a vector ``tid_map[k]`` stores the id of the rank which owns timestep ``k``). Typically, each rank first performs local operations on the timesteps ``T_k`` which it owns. Then all timesteps are broadcasted from the owning rank to all other ranks, so that finally each timestep ``T_k`` is known to all ranks.

Important member functions of ``cntr::distributed_timestep_array_moving<T>``:

.. list-table::
   :header-rows: 0

   * - ``distributed_timestep_array_moving(int n,int tc,int t0,int size,int sig,bool mpi)``
     - Constructor. ``tc``, ``size`` and ``sig`` are parameters of the individual timesteps. ``mpi`` should be ``true`` (otherwise each rank owns each timestep, and MPI operations are ignored). The constructor automatically determines the ownership map ``tid_map``.
   * - ``int tc()``, ``int size()``, ``int t0()``, ``int sig()``, ``int n()``
     - Return the respective parameters.
   * - ``int tid()``, ``int ntasks()``
     - Returns rank id (``tid``) of the local process and number of MPI ranks (``ntasks``).
   * - ``void clear(void)``
     - All data set to ``0``.
   * - ``std::vector<cntr::herm_matrix_timestep_view<T> > G(void)``
     - Return a handle to the vector of timesteps.
   * - ``void mpi_bcast_block(int k)``
     - Broadcast timestep :math:`T_k` from its owner to all other ranks.
   * - ``void mpi_bcast_all(void)``
     - Broadcast all timesteps from their owner to all other ranks.
   * - ``std::vector<int> data().tid_map()``
     - Return the vector ``tid_map``, where ``tid_map[k]`` for ``0 <= k < n`` is the rank which owns timestep ``T_k``.

.. note::

   In the implementation, the ``distributed_timestep_array_moving`` is not an array of ``cntr::herm_matrix_timestep_moving<T>`` variables, but a sufficiently large contiguous data block and an array of shallow ``cntr::herm_matrix_timestep_moving_view<T>`` objects. This allows for an easier adjustment of the size of the timesteps.

**Example:**

.. code-block:: cpp

   // MPI initialized with ntasks=6 ranks, local rank has rank-ID tid
   int nk=6;
   // ... set tc,t0,size1,sig ...
   GTRUNC Gloc(tc,t0,size1,sig);
   std::vector<GTRUNC> Gk(nk);
   cntr::distributed_timestep_array_moving<double> TARRAY(nk,tc,t0,size1,FERMION,true);
   for(int k=0;k<nk;k++){
       if(tid_map[k]==k){
           cout << "rank " << tid << " owns k= " << k << endl;
           Gk[k]=GTRUNC(tc,t0,size1,sig); // allocate memory for full Green's function Gk only on ranks which own k
       }else{
           cout << "rank " << tid << " does not own k= " << k << endl;
       }
   }

   // typical simulations at a given timestep:
   TARRAY.set_t0(tstp);
   for(int k=0;k<nk;k++){
       if(tid_map[k]==k){
           // ... rank owns k do some heavy numerics on Gk[k]
           // then copy timestep to TARRAY:
           TARRAY.G(k).set_timestep(0,Gk[k]); // read in data from Gk[k] to the array
       }
   }
   TARRAY.mpi_bcast_all();
   // now on all ranks and for all k TARRAY.G(k) contains the data of Gk[k]
   // this can be used, e.g., to calculate a local Green's function:
   Gloc.clear_timestep(0);
   for(int k=0;k<nk;k++) Gloc.incr_timestep(0,TARRAY.G(k),1.0/nk);

.. _Ptrunc13:

HDF5 usage
----------

.. _Ptrunc13S1:

Overview
~~~~~~~~

``libcntr`` uses HDF5 to store the basic data types for contour functions to disk.
HDF5 is an open source library and file format for numerical data which is widely used in the field of scientific computing. The format has two building blocks:

- *data sets*: general multi-dimensional arrays of a single type
- *groups*: containers which can hold data sets and other groups.

Hence, by nesting groups, it is possible to store arbitrarily complicated structured data, and to create a file-system-like hierarchy where groups can be indexed using standard POSIX format, e.g. ``/path/to/data``.

The ``libcntr`` library comes with helper functions to store the basic contour response function data types in HDF5 with a predefined structure of groups and data sets, defined in the header ``cntr/hdf5/hdf5_interface.hpp``. For example a ``herm_matrix`` response function is stored as a group with a data set for each contour component ``mat`` (:math:`g^M(\tau)`), ``ret`` (:math:`g^R(t, t')`), ``les`` (:math:`g^<(t, t')`), and ``tv`` (:math:`g^\rceil(t, \tau)`), respectively. The retarded and lesser components are stored in upper and lower triangular contiguous time order respectively.

.. _Ptrunc13S2:

Reading/writing to HDF5 files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To store a contour Green's function ``G`` of type ``cntr::herm_matrix_moving`` or read it from an HDF5 file, one can use the member functions ``cntr::herm_matrix_moving::write_to_hdf5`` and ``cntr::herm_matrix_moving::read_from_hdf5``. This stores/reads the attributes ``tc``, ``sig``, ``size1``, ``size2``, ``element_size`` and the Green's function component's data sorted in groups ``ret``, ``les``.

**HDF5 writing functions:**

.. list-table::
   :header-rows: 0

   * - ``G.write_to_hdf5(hid_t group_id)``
     - Stores the ``cntr::herm_matrix_moving`` (attributes and data) to a given HDF5 group.
   * - ``G.write_to_hdf5(hid_t group_id, const char *groupname)``
     - Stores the ``cntr::herm_matrix_moving`` (attributes and data) to a given HDF5 group with given groupname.
   * - ``G.write_to_hdf5(const char *filename, const char *groupname)``
     - Write data and attributes from the ``cntr::herm_matrix_moving`` to an HDF5 file under the given group name.

**HDF5 reading functions:**

.. list-table::
   :header-rows: 0

   * - ``G.read_from_hdf5(hid_t group_id)``
     - Reads the ``cntr::herm_matrix_moving`` (attributes and data) from a given HDF5 group handle.
   * - ``G.read_from_hdf5(hid_t group_id, const char *groupname)``
     - Reads the ``cntr::herm_matrix_moving`` (attributes and data) from a given HDF5 group handle with given group name.
   * - ``G.read_from_hdf5(const char *filename, const char *groupname)``
     - Read all data and attributes from an HDF5 file from a given group into the ``cntr::herm_matrix_moving``.

**HDF5 writing: given timeslice:**

.. list-table::
   :header-rows: 0

   * - ``G.write_timestep_to_hdf5(int delti,hid_t group_id)``
     - Write data and attributes from the ``delti``-th time step of ``cntr::herm_matrix_moving`` object ``G`` to a given HDF5 group handle.
   * - ``G.write_timestep_to_hdf5(int delti,hid_t group_id, const char *groupname)``
     - Write data and attributes from the ``delti``-th time step of ``cntr::herm_matrix_moving`` object ``G`` to a given HDF5 group handle with given group name.
   * - ``G.write_timestep_to_hdf5(int delti,const char *filename, const char *groupname)``
     - Write data and attributes from the ``delti``-th time step of ``cntr::herm_matrix_moving`` object ``G`` to an HDF5 file under the given group name.

**Example:**

.. code-block:: cpp

   #include <cntr/cntr.hpp>
   ..
   // Create a contour Green's function
   int tc = 200, t0 = 400, norb = 1;
   GTRUNC A(tc, t0, norb, FERMION);

   // Open HDF5 file and write components of the Green's function A into a group g.
   std::string filename = "data.h5";
   A.write_to_hdf5(filename.c_str(), "g");

If the file ``data.h5`` has been written previously with ``write_to_hdf5``, one can read it with the member function ``read_from_hdf5``

.. code-block:: cpp

   // Open HDF5 file and read group g. The result is saved into the Green's function B
   GTRUNC B;
   B.read_from_hdf5(filename.c_str(), "g");

the parameters ``(tc,t0,size1,sig)`` and the data of ``B`` are modified according to the information in the file.

To simplify postprocessing of contour Green's functions, NESSi also provides the python module ``ReadCNTRhdf5.py`` for reading the HDF5 format (using the python modules ``numpy`` and ``h5py``), producing python objects with the contour components as members.

.. _Ptrunc14:

Python tools
------------

As a part of the NESSi package, we provide python tools, which include:

- scripts for pre-processing to assist the use of programs based on ``libcntr``
- scripts for reading and post-processing Green's functions from HDF5 format via the ``h5py`` python package.

The python tools can be found in ``libcntr/python`` (for python2) or ``libcntr/python3`` (for python3 compatibility).

.. _Ptrunc14S01:

Creation of an input file
~~~~~~~~~~~~~~~~~~~~~~~~~~

An input file for a custom program can be generated using the ``ReadCNTR.py`` python module:

**Example:**

.. code-block:: python

   from ReadCNTR import write_input_file

       inp = {
           'nt': nt,    # number of points on the real axis before the truncation
           'beta':beta, # inverse temperature f
           'ntau':ntau, # number of points on the Matsubara (imaginary) axis
           'h':h,       # timestep interval
           'tmax':tmax, # absolute number of timesteps 
           'tc':tc      # truncation time
           }

       write_input_file(inp, input_file)

.. note::

   To make the python module available to a python script, the module ``ReadCNTR.py`` needs to be in the python include path. For python2, this is achieved by setting ``export PYTHONPATH=<path>/nessi/libcntr/python``, while ``export PYTHONPATH=<path>/nessi/libcntr/python3`` makes the python3 version accessible.

.. _Ptrunc14S02:

Reading Green's functions from HDF5
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The python module unrolls the moving window of the ``ret`` and ``les`` components making it simple to work with time slices.

**Example 1:**

To store the imaginary part of the retarded Green's function :math:`\textrm{Im}[G_{k,k^\prime}^\mathrm{R}(t={\tt (t0-i) \Delta t}, t^\prime={\tt (t0-i-j)}\Delta t)]` at orbital indices :math:`k=k^\prime=0` as a function of :math:`j` into a text file at :math:`i=5` we may use the commands

.. code-block:: python

   import h5py
   from ReadCNTRhdf5 import read_group_trunc

   # read the components of the contour function from hdf5 file
   with h5py.File('data.h5', 'r') as fd:
       g = read_group_trunc(fd).g

   # choose a timeslice for the retarded Green's function
   indxt=5
   G_ret=g.ret[indxt,:,0,0].imag
   with open("data_output.txt",'w') as output:     # create a text file
       output.write("%s \t" % G_ret)               # write data into it

**Example 2:**

The same way one can save the lesser component of the Green's function :math:`\textrm{Im}[G^<(t={\tt (t0-i)}\Delta t, t^\prime={\tt (t0-i-j)\Delta t})]` as a function of ``j`` into a text file

.. code-block:: python

   import h5py
   from ReadCNTRhdf5 import read_group

   # read the components of the contour function from hdf5 file
   with h5py.File('data.h5', 'r') as fd:
       g = read_group(fd).g

   # choose a timeslice for the lesser Green's function
   indxt=5
   G_les=g.les[indxt,:,0,0].imag
   with open("data_output.txt",'w') as output:     # create a text file
       output.write("%s \t" % G_les)               # write data into it

There is also a ``read_imp_trunc_h5file`` routine providing the same functionality for truncated contour objects as ``read_imp_h5file`` for general contour objects.
