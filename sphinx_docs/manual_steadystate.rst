.. _PNess:

Manual for steady-state Kadanoff-Baym equations
================================================

.. contents::
   :local:
   :depth: 2

.. _PNess00:

Steady state Keldysh formalism
-------------------------------

Another framework for solving the KB equations is the Keldysh formalism for non-equilibrium steady states (NESSs). In comparison to two-time functions defined on the three-legged KB contour :math:`\mathcal C = \mathcal C_1 \cup \mathcal C_2 \cup \mathcal C_3`, the main assumption in the steady-state formalism is the absence of correlations with the initial noninteracting state on the imaginary time branch :math:`\mathcal C_3`, such that the mixed self-energies :math:`\Sigma^{\rceil}` and :math:`\Sigma^{\lceil}` vanish. The vertical branch can then be shifted to :math:`t = - \infty` and is eliminated from the equations, such that time arguments are restricted to the two-branch contour :math:`\mathcal C_K = \mathcal C_1 \cup \mathcal C_2`. Second, any contour-ordered two-time Green's function :math:`G(t,t') = -i \langle \mathcal{T}_{\mathcal C} c(t) c^{\dagger}(t') \rangle` in a NESS exhibits time-translational invariance :math:`G(t,t') = G(t-t')`, so that a numerical treatment is possible in Fourier representation.

**Green's functions**

To represent NESS Green's functions, we keep the retarded component :math:`G^R(t)` and lesser component :math:`G^<(t)` on an equidistant real-time grid as the two non-redundant components.
Frequency-dependent Green's functions are given by the spectral representation

.. _GFs:

.. math::
   :label: fft1

   G^R(t) = -i\theta(t) \int d\omega \, A(\omega) e^{-i\omega t},

.. math::
   :label: fft2

   G^R(\omega+i0) = \int_0^\infty dt \, G^R(t) e^{i(\omega+i0) t},

.. math::
   :label: fft3

   G^<(t) = \int \frac{d\omega}{2\pi} \, G^<(\omega) e^{-i\omega t},

.. math::
   :label: fft4

   G^<(\omega) = \int dt \, G^<(t) e^{i\omega t},

with the spectral function

.. math::

   A(\omega) =-\frac{1}{2\pi} \left(G^R(\omega+i0)-G^R(\omega+i0)^\dagger\right).

In equilibrium at temperature :math:`1/\beta`, the fluctuation dissipation function implies

.. math::
   :label: Feq

   G^<(\omega)= -\xi 2\pi i F_\mu(\omega) A(\omega),

where the sign :math:`\xi` is :math:`\xi=\pm1` for bosonic (fermionic) Green's functions, and the distribution function is :math:`F_\mu(\omega) =1/(e^{\beta(\omega-\mu)}-\xi)`.

**Dyson equation**

.. _Dyson_ness:

In order to solve the Dyson equation in the steady state, we rewrite the Kadanoff-Baym equations in the time-translationally invariant case. The equation for the retarded Green's function becomes

.. math::
   :label: nessdyson1

   G^R(\omega+i0^{+}) = \big[\omega+i0^{+} - \epsilon-\Sigma^R(\omega+i0^{+})\big]^{-1}

in Fourier representation and the lesser component is given by

.. math::
   :label: nessdyson2

   G^<(\omega) = G^R(\omega+i0) \Sigma^<(\omega) G^R(\omega+i0)^\dagger.

.. _PNess01:

Overview
--------

The current implementation of the steady state extension to NESSi has two versions, which are structured into two different namespaces ``ness`` and ``ness2``. We suggest to the user to start with the newer version ``ness2``, it covers and extends all functionalities of the older version. For the sake of completeness, the old version ``ness`` is also documented on this page.

.. list-table::
   :header-rows: 0

   * - ``ness2``
     - Namespace for current implementation of the NESS code, recommended version.
   * - ``ness``
     - Namespace for old implementation of the NESS code.

.. _PNess01S01:

ness2
~~~~~

This is a quick overview over the most important classes and functions in the ``ness2`` namespace:

.. list-table::
   :header-rows: 0

   * - ``ness2::fft_array``
     - Basic class containing the data structure for holding arrays of matrices for frequency and time domain, supports element access and algebraic manipulations.
     - Section :ref:`PNess02S01`
   * - ``ness2::herm_matrix_ness``
     - Class to represent steady state Green's function (lesser and retarded component), supports element access, Fourier transforms, equilibrium Green's function construction and reading and writing to file and two-time functions and other utilities.
     - Section :ref:`PNess02S02`
   * - ``ness2::dyson``
     - Solve the Dyson equation in the steady state.
     - Section :ref:`PNess05`
   * - ``ness2::green_equilibrium_ness``
     - Construct equilibrium Green's function from a given density of states.
     - Section :ref:`PNess03`
   * - ``ness2::Bubble1_ness``
     - Calculate particle-hole bubble diagram for two given Green's functions.
     - Section :ref:`PNess04`
   * - ``ness2::Bubble2_ness``
     - Calculate particle-particle bubble diagram for two given Green's functions.
     - Section :ref:`PNess04`
   * - ``ness2::cntr2ness``
     - Write from two-time function into a steady state Green's function, also vice versa is supported.
     - Section :ref:`PNess07S02`

.. _PNess01S02:

ness
~~~~

This is a quick overview over the most important classes and functions in the ``ness`` namespace:

.. list-table::
   :header-rows: 0

   * - ``ness::GF``
     - Class to represent steady state Green's function (lesser and retarded component), supports element access, Fourier transforms, equilibrium Green's function construction and reading and writing to file and two-time functions and other utilities (either frequency or time domain).
     - Section :ref:`PNess02S03`
   * - ``ness::GF_pair``
     - Class to represent steady state Green's function (lesser and retarded component), both in frequency and time domain, helper class ``fft_solver`` provides Fourier transform.
     - Section :ref:`PNess02S05`
   * - ``ness::fft_solver``
     - Helper class to provide Fourier forward and back transform.
     - Section :ref:`PNess02S04`
   * - ``ness::dyson_ness``
     - Solve the Dyson equation in the steady state.
     - Section :ref:`PNess05`
   * - ``ness::green_equilibrium_ness``
     - Construct equilibrium Green's function from a given density of states.
     - Section :ref:`PNess03`
   * - ``ness::Bubble1_ness``
     - Calculate particle-hole bubble diagram for two given Green's functions.
     - Section :ref:`PNess04`
   * - ``ness::Bubble2_ness``
     - Calculate particle-particle bubble diagram for two given Green's functions.
     - Section :ref:`PNess04`

.. _PNess02:

Green's functions
-----------------

.. _ness_green_def:

.. _PNess02S01:

``ness2::fft_array``
~~~~~~~~~~~~~~~~~~~~

.. _PNess02S01A:

Overview
^^^^^^^^^

.. _DFT:

.. list-table::
   :header-rows: 0

   * - class
     - ``ness2::fft_array``

This class is a basic data structure for an object ``C`` containing two arrays ``C.time_`` and ``C.freq_``, which hold ``Nft`` square matrices stored as consecutive complex numbers. The class provides plans to perform a (non-normalized) discrete Fourier transform

.. math::
   :label: dft1

   \mathcal{F}\left\{ C \right\}[l] = \sum_{j=0}^{N_{\rm ft}-1} e^{i 2 \pi j l / N_{\rm ft}} \, C[j],

.. math::
   :label: dft2

   \bar{\mathcal{F}}\left\{ C \right\}[l] = \sum_{j=0}^{N_{\rm ft}-1} e^{-i 2 \pi j l / N_{\rm ft}} \, C[j].

using the FFT algorithm. In addition, the class offers routines for standard element access and basic algebraic operations. The ``ness2::fft_array`` models an array that is periodic in both time and frequency dimensions, with a period of ``Nft``. This periodic wrapping is handled automatically when accessing array elements. For most efficient performance of the FFT, the size ``Nft`` should be a power of two. An object of the class is characterized by the following two arguments:

- ``Nft`` (integer): FFT array length.
- ``size`` (integer): matrix size (``size x size``).

.. _PNess02S01B:

Constructors
^^^^^^^^^^^^^

.. list-table::
   :header-rows: 0

   * - ``fft_array()``
     - Default constructor, with ``Nft = 4``, ``size = 1`` and ``FFTW_FLAG = FFTW_ESTIMATE``.
   * - ``fft_array(int Nft, int size, unsigned FFTW_FLAG = FFTW_ESTIMATE)``
     - Standard constructor.
   * - ``fft_array(const fft_array& other)``
     - Copy constructor given a ``ness2::fft_array`` object ``other``.

The constructor has an optional third argument ``FFTW_FLAG``, which chooses the plan for the construction of the plan in the FFTW library used to do the Fourier transforms (default is ``FFTW_ESTIMATE``).

.. _PNess02S01C:

Accessing individual elements
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**General element access:**

The following member functions read and set elements of the time (``ness2::fft_array::time_``) or frequency (``ness2::fft_array::freq_``) part of the ``ness2::fft_array``, specified by the ``domain`` argument, which can take the values ``ness2::fft_domain::freq`` or ``ness2::fft_domain::time``. The parameter ``M`` is a ``Matrix`` type (complex eigen matrices are supported). For ``[set|get]_element``, time or frequency arguments outside the interval ``[0, Nft-1]`` are mapped back into this interval by adding an integer multiple of ``Nft``.

.. list-table::
   :header-rows: 0

   * - ``C.set_element(i,M,domain)``
     - :math:`(C_{\rm time|freq})[i]` set to matrix ``M``.
   * - ``C.get_element(i,M,domain)``
     - Matrix ``M`` set to :math:`C_{\rm time|freq}[i]`, ``M`` is automatically resized as needed.

To set a matrix element of :math:`(C_{\rm time|freq})` one can use the following utility, which operates on all ``Nft`` entries:

.. list-table::
   :header-rows: 0

   * - ``C.set_matrixelement(i1,i2, fft array &B, j1,j2, domain)``
     - :math:`(C_{\rm time|freq})_{i_1,i_2}` set to :math:`(B_{\rm time|freq})_{j_1,j_2}`.

.. _PNess02S01D:

Simple Algebra
^^^^^^^^^^^^^^^

The following methods cover basic algebraic operations on the ``ness2::fft_array`` object. The ``domain`` argument is set to ``ness2::fft_domain::time`` or ``ness2::fft_domain::freq`` to specify whether the operation applies to ``C.time_`` or ``C.freq_``.
All routines operate on all ``Nft`` entries. The parameter ``M`` is a ``Matrix`` type (complex eigen matrices are supported).

.. list-table::
   :header-rows: 0

   * - ``C.incr(fft_array &B, cplx a, domain)``
     - :math:`C_{\rm time|freq}` set to :math:`C_{\rm time|freq}+aB_{\rm time|freq}`.
   * - ``C.smul(cplx a, domain)``
     - :math:`C_{\rm time|freq}` set to :math:`aC_{\rm time|freq}`.
   * - ``C.set_zero(domain)``
     - :math:`C_{\rm time|freq}` set to 0 for all entries.
   * - ``C.left_multiply(M, domain)``
     - :math:`C_{\rm time|freq}` set to :math:`MC_{\rm time|freq}`.
   * - ``C.right_multiply(M, domain)``
     - :math:`C_{\rm time|freq}` set to :math:`C_{\rm time|freq}M`.
   * - ``C.left_multiply_hermconj(M, domain)``
     - :math:`C_{\rm time|freq}` set to :math:`M^\dagger C_{\rm time|freq}`.
   * - ``C.right_multiply_hermconj(M, domain)``
     - :math:`C_{\rm time|freq}` set to :math:`C_{\rm time|freq}M^\dagger`.

.. _PNess02S01E:

Fourier transform
^^^^^^^^^^^^^^^^^^

The following member functions perform the discrete Fourier transforms :eq:`dft1` and :eq:`dft2` using FFT:

.. list-table::
   :header-rows: 0

   * - ``C.fft_to_time()``
     - Set :math:`C_{\rm time}` to :math:`\bar{\mathcal{F}}\{C_{\rm freq}\}` using :eq:`dft2` and FFT.
   * - ``C.fft_to_freq()``
     - Set :math:`C_{\rm freq}` to :math:`{\mathcal{F}}\{C_{\rm time}\}` using :eq:`dft1` and FFT.

**Example:**

.. code-block:: cpp

   int N_grid = 100;
   int size = 2;
   double h = 0.1;
   int ex_id = 25;

   // create arrays
   fft_array new_array(N_grid, size);
   fft_array new_array2(N_grid, size);

   // matrices for operations
   cdmatrix A(size,size);
   cdmatrix B(size,size);
   cdmatrix C(size,size);
   A(1,0) = 1;
   A(0,1) = II; // fill A with values

   // element access
   new_array.set_element(ex_id,A,ness2::fft_domain::time); // read into new_array from A
   new_array.get_element(ex_id,B,ness2::fft_domain::time); // read into B from new_array
   new_array2.set_matrixelement(1,1, new_array, 0, 1, ness2::fft_domain::time); // write one matrix element in new_array2 from new_array

   // example operations
   new_array.incr(new_array2, 0.5, ness2::fft_domain::time);
   new_array.smul(2.0, ness2::fft_domain::time);
   new_array.left_multiply(A, ness2::fft_domain::time);
   new_array.right_multiply(B, ness2::fft_domain::time);
   new_array2.set_zero(ness2::fft_domain::time);
   new_array.fft_to_freq(); // FFT to frequency
   new_array.get_element(ex_id,C,ness2::fft_domain::freq); // e.g. read into C

.. _PNess02S01F:

Grid info
^^^^^^^^^^

.. list-table::
   :header-rows: 0

   * - class
     - ``ness2::grid_info``

This helper class describes the time and frequency grids with steps ``h`` and ``dw`` and length ``Nft`` used in ``ness2::fft_array`` and ``ness2::herm_matrix_ness``.

.. list-table::
   :header-rows: 0

   * - ``grid_info(Nft, h)``
     - Constructor to initialize grid with time step ``h`` and length ``Nft``.
   * - ``grid.set_h(new_h)``
     - Update time step ``h`` to ``new_h`` and adjust frequency step ``dw`` accordingly.
   * - ``grid.set_dw(new_dw)``
     - Update frequency step ``dw`` to ``new_dw`` and adjust time step ``h`` accordingly.
   * - ``grid.time_at(n)``
     - Get time value :math:`t` at grid index ``n``.
   * - ``grid.freq_at(k)``
     - Get frequency value :math:`\omega` at grid index ``k``.
   * - ``grid.time_grid()``
     - Generate the full wrapped time grid (``std::vector<double>``).
   * - ``grid.freq_grid()``
     - Generate the full wrapped frequency grid (``std::vector<double>``).

**Example:**

.. code-block:: cpp

   int N_grid = 100;
   double h = 0.1;
   double new_h = 0.2;
   double new_dw = 0.1;
   int ex_id = 25;

   // manipulate and show grid info
   grid_info new_grid(N_grid, h); // create a grid
   cout << "Number of grid points: " << new_grid.Nft_ << endl;
   cout << "Time step: " << new_grid.h_ << endl;
   cout << "Frequency step: " << new_grid.dw_ << endl;

   new_grid.set_h(new_h); // update time step
   cout << "New time step: " << new_grid.h_ << endl;
   cout << "New frequency step: " << new_grid.dw_ << endl;

   new_grid.set_dw(new_dw); // update frequency step
   cout << "New time step: " << new_grid.h_ << endl;
   cout << "New frequency step: " << new_grid.dw_ << endl;
   cout << "Time at example index: " << new_grid.time_at(ex_id) << endl;
   cout << "Frequency at example index: " << new_grid.freq_at(ex_id) << endl;

   std::vector<double> full_time_grid = new_grid.time_grid(); // store full time grid in vector
   std::vector<double> full_freq_grid = new_grid.freq_grid(); // store full frequency grid in vector

.. _PNess02S02:

``ness2::herm_matrix_ness``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _PNess02S02A:

Overview
^^^^^^^^^

.. list-table::
   :header-rows: 0

   * - class
     - ``ness2::herm_matrix_ness``

The data type ``ness2::herm_matrix_ness`` is used to represent a steady-state function ``G`` with hermitian symmetry (non-hermitian Green's functions are not supported in the current extension of ``libcntr``).
An object of this class is initialized by:

- ``Nft`` (integer): FFT array length.
- ``size`` (integer): matrix size (``size x size``).

An object ``G`` of type ``ness2::herm_matrix_ness`` contains a pair of two members ``ness2::herm_matrix_ness::ret_`` and ``ness2::herm_matrix_ness::les_`` of type ``ness2::fft_array``, representing :math:`G^R` and :math:`G^<` in time and frequency. Time arguments thereby correspond to an equidistant grid

.. math::
   :label: timegrid

   t\in \{jh: j=-N_{\rm ft}/2,...,N_{\rm ft}/2-1\},

and the frequencies represent the dual FFT grid

.. math::
   :label: w-grid

   \omega \in  \{ \Delta_\omega j: j= -N_{\rm ft}/2,...,N_{\rm ft}/2-1,\Delta_\omega= \frac{2\pi}{hN_{\rm ft}}\}.

We require :math:`N_{\rm ft}` to be even, and use the conventional layout of the arrays where :math:`G^{<|R}(jh)` is represented by the :math:`j`-th element of ``G.[les|ret]_.time_``, if :math:`0\le j< N_{\rm ft}/2`, and by the :math:`(j+N_{\rm ft})`-th element if :math:`-N_{\rm ft}/2\le j< 0` (analogous for frequency). For convenience we provide a small helper class ``ness2::fft_grid``, see :ref:`PNess02S01F`.

.. _PNess02S02B:

Constructors
^^^^^^^^^^^^^

.. list-table::
   :header-rows: 0

   * - ``herm_matrix_ness()``
     - Default constructor, with ``Nft = 4``, ``size = 1``.
   * - ``herm_matrix_ness(Nft, size)``
     - Standard constructor.
   * - ``herm_matrix_ness(const herm_matrix_ness& other)``
     - Copy constructor given a ``herm_matrix_ness`` object ``other``.

The constructor again has an optional third argument ``FFTW_FLAG``, for the plan of the FFT construction on FFTW.

.. _PNess02S02C:

Accessing individual elements
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**General element access:**

The following member functions read and set elements of the time or frequency part (``domain`` argument set to ``ness2::fft_domain::freq`` or ``ness2::fft_domain::time``) of the retarded or lesser component of the ``ness2::herm_matrix_ness`` object. The parameter ``M`` is a ``Matrix`` type (complex eigen matrices are supported). As for ``ness2::fft_array`` ``[set|get]_element``, time or frequency (``i``) arguments outside the interval ``[0, Nft-1]`` are mapped back into this interval by adding an integer multiple of ``Nft``.

.. list-table::
   :header-rows: 0

   * - ``G.set_[les|ret](i,M,domain)``
     - :math:`(G_{\rm time|freq}^{<|R})[i]` set to matrix ``M``.
   * - ``G.get_[les|ret](i,M,domain)``
     - Matrix ``M`` set to :math:`G_{\rm time|freq}^{<|R}[i]`.
   * - ``G.retarded()``
     - Return reference to ``G.ret_``.
   * - ``G.lesser()``
     - Return reference to ``G.les_``.

.. _PNess02S02D:

Simple Algebra
^^^^^^^^^^^^^^^

The same functions as for ``ness2::fft_array`` are also defined for ``ness2::herm_matrix_ness``, the operation is performed for both retarded and lesser component.

.. list-table::
   :header-rows: 0

   * - ``G.incr(herm_matrix_ness &B, cplx a, domain)``
     - :math:`G^{R,<}_{\rm time|freq}` set to :math:`G^{R,<}_{\rm time|freq}+aB^{R,<}_{\rm time|freq}`.
   * - ``G.smul(cplx a, domain)``
     - :math:`G^{R,<}_{\rm time|freq}` set to :math:`a G^{R,<}_{\rm time|freq}`.
   * - ``G.set_zero(domain)``
     - :math:`G^{R,<}_{\rm time|freq}` set to 0 for all entries.
   * - ``G.left_multiply(M, domain)``
     - :math:`G^{R,<}_{\rm time|freq}` set to :math:`M G^{R,<}_{\rm time|freq}`.
   * - ``G.right_multiply(M, domain)``
     - :math:`G^{R,<}_{\rm time|freq}` set to :math:`G^{R,<}_{\rm time|freq}M`.
   * - ``G.left_multiply_hermconj(M, domain)``
     - :math:`G^{R,<}_{\rm time|freq}` set to :math:`M^\dagger G^{R,<}_{\rm time|freq}`.
   * - ``G.right_multiply_hermconj(M, domain)``
     - :math:`G^{R,<}_{\rm time|freq}` set to :math:`G^{R,<}_{\rm time|freq}M^\dagger`.

.. _PNess02S02E:

Integral transforms
^^^^^^^^^^^^^^^^^^^^

For the Fourier transform of Green's functions, we must distinguish between plain :ref:`DFT <DFT>`, which can be accessed by the call ``[ret|les]_.to_[time|freq]()`` to the members of ``ness2::herm_matrix_ness``, and approximations to the :ref:`GFs <GFs>` integral transforms.
The transformation ``G.transform_to_freq(h, METHOD)`` computes an approximation to the integral transforms :eq:`fft2`, :eq:`fft4` on the frequency grid :eq:`w-grid`, assuming a time grid with timestep ``h``. The integrals are computed with boundaries ``[0,T]`` (for :math:`G^R`) and ``[-T, T]`` (for :math:`G^<`), with :math:`T = (Nft/2 - 1)h`. If the argument ``METHOD`` is omitted or given by the keyword ``FFT_TRAPEZ`` (default), we compute the Fourier integral in the trapezoidal rule. The reverse transformation ``G.transform_to_time(h, METHOD)`` computes an approximation to the Fourier integrals :eq:`fft1`, :eq:`fft3` on the time grid :eq:`timegrid`. For both directions, we also provide an implementation ``METHOD=FFT_CUBIC``, where the integrals are approximated by the exact Fourier transform of a piecewise cubic interpolating function.

.. list-table::
   :header-rows: 0

   * - ``G.integral_transform_to_freq(h, METHOD)``
     - Compute the frequency-domain Fourier integrals :eq:`fft2`, :eq:`fft4` of the Green's function, assuming a timestep ``h``.
   * - ``G.integral_transform_to_time(h, METHOD)``
     - Compute the time-domain Fourier integrals :eq:`fft1`, :eq:`fft3` of the Green's function, assuming a timestep ``h``.

.. _PNess02S02F:

Density matrix
^^^^^^^^^^^^^^^

The density matrix in the steady state is directly evaluated by calculating the lesser component of the Green's function:

.. math::

   \rho_{A} = \xi i A^<(0),

with :math:`\xi=\pm1` for bosonic (fermionic) Green's functions. The equal-time convolution

.. math::

   \rho_{AB} = \xi_A \xi_B i (A\ast B)^<(0),

can be obtained from the Fourier transform of :math:`C=A\ast B`.
The corresponding function calls are ``ness2::density_matrix`` and ``ness2::convolution_density_matrix``:

.. _PNess02S02E_ref:

.. list-table::
   :header-rows: 0

   * - ``density_matrix(result, bosefermi, A)``
     - ``result`` (complex matrix) is set to :math:`\rho_{A}`, for a function ``A`` of type ``ness2::herm_matrix_ness``.
   * - ``convolution_density_matrix( result, bosefermi, A, B, double h)``
     - ``result`` (complex matrix) is set to the steady-state :math:`\rho_{AB}`, for the convolution of two functions ``A`` and ``B`` of type ``ness2::herm_matrix_ness``.

.. _PNess02S02G:

File I/O
^^^^^^^^^

File input and output from text files and also HDF5 files is supported.

**Reading and writing to text files:**

Text files start with an initial header line ``# Nft size precision`` containing apart from the size properties also the optional accuracy argument. The data is structured into blocks for the retarded and lesser component. Reading expects the same format as writing to text files.

.. list-table::
   :header-rows: 0

   * - ``G.print_to_file(filename, precision = 12)``
     - Print the ``ness2::herm_matrix_ness`` object ``G`` (``ret_``, ``les_`` and attributes ``Nft``, ``size``) to a text file. Precision argument is optional.
   * - ``G.read_from_file(filename, FFTW_FLAG = FFTW_ESTIMATE)``
     - Read ``ret_``, ``les_`` and attributes ``Nft``, ``size`` from a text file and store in ``ness2::herm_matrix_ness`` object ``G``, expects same format as ``print_to_file``. Second argument ``FFTW_FLAG`` optional to set the FFT plan of object ``G`` manually.

**Reading and writing to HDF5 files:**

To store a ``ness2::herm_matrix_ness`` or read it from a HDF5 file, one can use the member functions ``ness2::herm_matrix_ness::write_to_hdf5`` and ``ness2::herm_matrix_ness::read_from_hdf5``. This stores/reads the attributes ``Nft``, ``size1``, and the Green's function component's data sorted in groups ``ret``, ``les``.

.. list-table::
   :header-rows: 0

   * - ``G.write_to_hdf5(hid_t group_id)``
     - Stores the ``ness2::herm_matrix_ness`` (attributes and data) to a given HDF5 group.
   * - ``G.write_to_hdf5(hid_t group_id, const char *groupname)``
     - Stores the ``ness2::herm_matrix_ness`` (attributes and data) to a given HDF5 group with given groupname.
   * - ``G.write_to_hdf5(const char *filename, const char *groupname)``
     - Write data and attributes from the ``ness2::herm_matrix_ness`` to an HDF5 file under the given group name.
   * - ``G.read_from_hdf5(hid_t group_id)``
     - Reads the ``ness2::herm_matrix_ness`` (attributes and data) from a given HDF5 group handle.
   * - ``G.read_from_hdf5(hid_t group_id, const char *groupname)``
     - Reads the ``ness2::herm_matrix_ness`` (attributes and data) from a given HDF5 group handle with given group name.
   * - ``G.read_from_hdf5(const char *filename, const char *groupname)``
     - Read all data and attributes from an HDF5 file from a given group into the ``ness2::herm_matrix_ness``.

.. note::

   The HDF5 reading functions have a third optional argument ``FFTW_FLAG`` (default is ``FFTW_ESTIMATE``) to set the FFT plan of object ``G`` manually.

**Example:**

.. code-block:: cpp

   int N_grid = 100;
   int size = 2;
   double h = 0.1;
   int ex_id = 25;

   // Green's function object
   herm_matrix_ness G(N_grid, size);

   // complex matrices
   cdmatrix M1(size,size);
   cdmatrix M2(size,size);

   // Do something with G ...

   G.set_zero(ness2::fft_domain::freq);
   G.get_les(ex_id,M1,ness2::fft_domain::time); // read in from G into M1
   density_matrix(M2, -1, G); // get (fermionic) density matrix from G, save in M2
   G.integral_transform_to_freq(h, FFT_TRAPEZ); // evaluate frequency-domain Fourier integrals

   // save G in hdf5 file
   hid_t file_id = open_hdf5_file(flout); // char *flout points to the filename
   store_int_attribute_to_hid(file_id, "N_grid", N_grid); // write N_grid to group "/N_grid"
   // same for size, h ...
   G.write_to_hdf5(file_id, "G"); // new group "/G" in file
   close_hdf5_file(file_id);

.. _PNess02S03:

``ness::GF``
~~~~~~~~~~~~~

.. _PNess02S03A:

Overview
^^^^^^^^^

.. list-table::
   :header-rows: 0

   * - class
     - ``ness::GF``

This class contains the data structures for representing one-time contour Green's functions :math:`G(t)` on the equidistantly discretized real-time Keldysh contour :math:`\mathcal{C}_K` with points :math:`\{i h: i=0,1,2,...,{\tt nt-1}\}` or equivalently on an equidistant real-frequency grid :math:`\{ i d \omega: i=0,1,2,...,{\tt nfreq}-1\}` with spacing :math:`h` or :math:`d \omega` respectively. Note that the time argument of a Green's function object :math:`G(t)` is understood as a time difference in Keldysh steady state theory and the class ``ness::GF`` stores the positive time part. The Green's functions are initialized by the following parameters:

- ``dgrid`` (double): spacing of the time or frequency discretization.
- ``ngrid`` (integer): number of discretization points on the time or frequency axis.
- ``size1`` (integer): first orbital or any quantum number dimension.
- ``size2`` (integer): second orbital or any quantum number dimension. Each Green's function :math:`G(t)` is a matrix of dimension ``size1`` :math:`\times` ``size2``.
- ``sig`` (integer): takes the values ``-1`` (fermion) or ``+1`` (boson).
- ``gf_type`` (integer): documents the type of the Green's function and takes ``0`` (time type: ``time_gf``) or ``1`` (frequency type: ``freq_gf``).

.. _PNess02S03B:

Constructors
^^^^^^^^^^^^^

.. list-table::
   :header-rows: 0

   * - ``GF()``
     - Default constructor, does not allocate memory, sets all parameters to ``0``.
   * - ``GF(double dgrid, int ngrid, int size1, int size2, int sign, int gf_type)``
     - Allocates memory and sets all entries to ``0``.
   * - ``GF(double dgrid, int ngrid, int size1, int sign, int gf_type)``
     - Allocates memory and sets all entries to ``0`` with ``size1=size2``.
   * - ``GF(const GF& g)``
     - Constructs a ``GF`` object as a copy of a given ``GF`` object ``g``.

.. _PNess02S03C:

Accessing individual elements
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 0

   * - ``G.Retarded``
     - Return reference to retarded component.
   * - ``G.Lesser``
     - Return reference to lesser component.

The following member functions read and set components or submatrices of a Green's function ``G`` from another Green's function ``gf``:

.. list-table::
   :header-rows: 0

   * - ``G.set_element(int i1, int i2, const GF& gf, int s1, int s2)``
     - Set ``G.Retarded[n](i1,i2) = gf.Retarded[n](s1,s2)`` and ``G.Lesser[n](i1,i2) = gf.Lesser[n](s1,s2)`` for all ``n``.
   * - ``G.get_element(int i1, int i2, const GF& gf, int s1, int s2)``
     - Set ``gf.Retarded[n](s1,s2) = G.Retarded[n](i1,i2)`` and ``gf.Lesser[n](s1,s2) = G.Lesser[n](i1,i2)`` for all ``n``.
   * - ``G.set_submatrix(const std::vector<int> & subid1, const std::vector<int> & subid2, const GF& gf)``
     - Set ``G.Retarded[n](subid1[i1],subid2[i2]) = gf.Retarded[n](i1,i2)`` and corresponding lesser for all ``i1,i2`` and ``n``.
   * - ``G.get_submatrix(const std::vector<int> & subid1, const std::vector<int> & subid2, const GF& gf)``
     - Set ``gf.Retarded[n](i1,i2) = G.Retarded[n](subid1[i1],subid2[i2])`` and corresponding lesser for all ``i1,i2`` and ``n``.

.. _PNess02S03D:

Scalar multiplication and incrementation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 0

   * - ``G.smul(cplx weight)``
     - ``G`` set to ``weight`` :math:`\cdot` ``G`` for all times / frequencies.
   * - ``G.incr(GF &gf, cplx weight)``
     - ``G`` set to ``G`` + ``weight`` :math:`\cdot` ``gf`` for all times / frequencies.

- Size consistency ``G.ngrid_==gf.ngrid_`` and ``G.size1_`` :math:`\cdot` ``G.size2_==gf.size1_`` :math:`\cdot` ``gf.size2_`` as well as ``G.gf_type_==gf.gf_type_`` required.

.. _PNess02S03E:

Density matrix
^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 0

   * - ``G.Lesser[0]``
     - Density matrix :math:`\rho_{ab} = i \xi G^<_{a,b}(0)` (:math:`\xi` is the bosonic/fermionic sign).

.. _PNess02S03F:

File I/O
^^^^^^^^^

.. list-table::
   :header-rows: 0

   * - ``void ness::outputGF(GF &G, std::string filename)``
     - Writes the specified Green's function ``G`` into a text file ``filename``.
   * - ``void ness::readGF(GF &G, std::string filename)``
     - Reads the contents of a text file ``filename`` into a Green's function object ``G``.

The output or input text file is ordered: Frequency / time, :math:`\Re \{G^{\mathrm{R}}\}`, :math:`\Im \{G^{\mathrm{R}}\}`, :math:`\Im \{G^{<}\}`, :math:`\Re \{G^{<}\}`, ... (row-major order for orbital indices ``size1_``, ``size2_``).

.. _PNess02S04:

``ness::fft_solver``
~~~~~~~~~~~~~~~~~~~~~

.. _PNess02S04A:

Overview
^^^^^^^^^

.. list-table::
   :header-rows: 0

   * - class
     - ``ness::fft_solver``

This class contains the data structures for transforming a Green's function :math:`G(t)` on the equidistantly discretized real-time Keldysh contour :math:`\mathcal{C}_K` with points :math:`\{i h: i=0,1,2,...,{\tt nt-1}\}` to a Green's function :math:`G(\omega)` on an equidistant real-frequency grid :math:`\{ i d \omega: i=0,1,2,...,{\tt nfreq}-1\}` with spacing :math:`h` or :math:`d \omega` respectively. We use the fast Fourier transform (FFT) library FFTW3. The class is defined by two parameters:

- ``nt`` (integer): number of time points of the ``GF`` object.
- ``FFTW_FLAG`` (size_t): The FFTW3 flag, this argument is usually either ``FFTW_MEASURE`` or ``FFTW_ESTIMATE``.

.. _PNess02S04B:

Constructor
^^^^^^^^^^^^

.. list-table::
   :header-rows: 0

   * - ``fft_solver(int nt, size_t FFTW_FLAG)``
     - Default constructor, creates necessary FFTW3 structures (input and output array, FFT plans).

.. _PNess02S04C:

Member functions
^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 0

   * - ``fft.to_time(GF &outG, const GF &inG, int set_dt = 0)``
     - Transforms the contents of the Green's function ``inG`` of type ``freq_gf`` to time and writes them into the Green's function ``outG`` of type ``time_gf``. If ``set_dt`` is set to 1, ``outG.dgrid_`` is automatically set from ``inG.dgrid_``. It is important that ``inG.ngrid_ = fft.nfreq_`` and ``outG.ngrid_ = 3 fft.nfreq_ / 2 + 1``.
   * - ``fft.to_freq(GF &outG, const GF &inG, int set_dw = 0)``
     - Transforms the contents of the Green's function ``inG`` of type ``time_gf`` to frequency and writes them into the Green's function ``outG`` of type ``freq_gf``. If ``set_dw`` is set to 1, ``outG.dgrid_`` is automatically set from ``inG.dgrid_``. It is important that ``inG.ngrid_ = 3 fft.nfreq_ / 2 + 1`` and ``outG.ngrid_ = fft.nfreq_``.

.. _PNess02S05:

``ness::GF_pair``
~~~~~~~~~~~~~~~~~~

.. _PNess02S05A:

Overview
^^^^^^^^^

.. list-table::
   :header-rows: 0

   * - class
     - ``ness::GF_pair``

This class contains the data structures for simultaneously handling a time Green's function object and frequency Green's function object. The advantage of this class is that both objects (``freq`` and ``time``) are constructed with the necessary amount of data points for a Fourier transform. The Fourier transform can be called as a member function (``to_time|to_freq``) on the members of the class. The class is thus defined by all the parameters of a time ``GF`` object.

.. _PNess02S05B:

Constructors
^^^^^^^^^^^^^

.. list-table::
   :header-rows: 0

   * - ``GF_pair()``
     - Default constructor, does not allocate memory.
   * - ``GF_pair(double h, int nt, int size, int sign, fft_solver & fft, int long_freq = 0)``
     - Allocates memory and sets all entries to ``0`` with ``size1=size2=size``.
   * - ``GF_pair(const GF_pair & gp)``
     - Constructs a ``GF_pair`` object as a copy of a given ``GF_pair`` object ``gp``.

.. _PNess02S05C:

Accessing members
^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 0

   * - ``gp.freq``
     - Frequency ``GF`` object with correct spacing and number of points.
   * - ``gp.time``
     - Time ``GF`` object with correct spacing and number of points for the FFT.

.. _PNess02S05D:

Member functions
^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 0

   * - ``gp.to_time()``
     - Fourier transforms ``gp.freq`` to time and writes the result into ``gp.time``.
   * - ``gp.to_freq()``
     - Fourier transforms ``gp.time`` to frequency and writes the result into ``gp.freq``.

.. _PNess03:

Equilibrium Green's functions
------------------------------

Also in the steady state we provide tools for calculating the Green's function
from a given spectral function solving the integrals :eq:`fft1`, :eq:`fft3` or to force a given Green's function into equilibrium using relation :eq:`Feq`. The provided functions also work in the namespace ``ness``.

.. _PNess03S01:

Constructing a Green's function from a given DOS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

As defined in the manual, we can construct a Green's function from a given density of states :math:`A(\omega)` in the steady state:

.. math::

   G(t) = \int d \omega \,A(\omega) g_{\omega}(t) \ ,

where :math:`g_\omega(t)` is the single-orbital noninteracting Green's function with energy :math:`\omega`. The latter is defined as (:math:`F_\xi(x)=1/(e^{\beta x}-\xi )` is the Fermi/Bose distribution):

- :math:`g_\omega^R (t) = -i e^{-i t \omega}`
- :math:`g_\omega^< (t) = -i\xi  e^{-i t \omega} F_\xi(\omega-\mu)`

For the class ``ness2::herm_matrix_ness`` this is accomplished by the following call of ``ness2::green_equilibrium_ness``:

.. list-table::
   :header-rows: 0

   * - ``ness2::green_equilibrium_ness(sign, G, DOS &dos, beta, mu,h, METHOD, limit=100, nn=20)``
     - Set :math:`G^{<|R}(t)` to equilibrium Green's function for Bosons (``sign=+1``) or Fermions (``sign=-1``).

The user provides the inverse temperature ``beta``, the time step ``h`` and the chemical potential ``mu`` (set to 0.0 by default). ``DOS`` is a class representing :math:`A(\omega)`, which provides an operation ``dos(double omega)`` to return :math:`A(\omega)`, and the numbers ``dos.lo_`` and ``dos.hi_`` representing the lower and upper bound of the support of :math:`A(\omega)`. One may use keywords ``BOSON`` and ``FERMION`` for the ``sign`` :math:`\pm 1`. ``METHOD`` can be ``FFT_TRAPEZ`` or ``FFT_ADAPTIVE``, for ``FFT_TRAPEZ`` arguments ``nn`` and ``limit`` are ignored.

- If the argument ``METHOD`` is ``FFT_TRAPEZ`` the integral is computed as in ``transform_to_time``, after initializing the imaginary part of the frequency-dependent components ``ret_.freq_`` and ``les_.freq_`` using the ``dos`` function on the frequency grid :eq:`w-grid`. Note that the routine is implemented only for scalar :math:`A(\omega)`; if :math:`G` is matrix-valued, it is set to a diagonal matrix.
- Alternatively, if ``METHOD = FFT_ADAPTIVE``, integrals are computed using the adaptive Fourier integral of ``libcntr``, with a subdivision of the domain of :math:`A(\omega)` in at most ``limit`` intervals of ``nn`` points (default is ``limit=100`` and ``nn=20``). The adaptive Fourier integration is more accurate but slower, because it does not exploit the FFT algorithm.

**Example:**

.. code-block:: cpp

   int N_grid = 100;
   int size = 2;
   double h = 0.1;
   double beta= 0.1;
   double mu= 0.0;
   int sign= -1;
   herm_matrix_ness G(N_grid, size);
    
   cntr::bethedos dos;
   dos.V_ = 0.25;
   dos.lo_ = -0.5;
   dos.hi_ = 0.5;
    
   ness2::green_equilibrium_ness(-1, G, dos, beta, mu, h, FFT_TRAPEZ);

.. _PNess03S02:

Setting a Green's function to thermal equilibrium
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Also for this class, we provide a method ``ness2::herm_matrix_ness::force_equilibrium``, which does a Fourier transform of :math:`G^R(t)` to frequency (analogous to ``ness2::herm_matrix_ness::integral_transform_to_freq`` with ``METHOD = FFT_TRAPEZ``) to initialize :math:`G^<(\omega)` on the frequency grid :eq:`w-grid` according to the equilibrium dissipation fluctuation relation :eq:`Feq`, and transforms :math:`G^<(\omega)` to :math:`G^<(t)` analogous to ``ness2::herm_matrix_ness::integral_transform_to_time``.

.. list-table::
   :header-rows: 0

   * - ``G.force_equilibrium(sign, beta, mu, h)``
     - Set :math:`G^<(t)` according to :eq:`fft3` and the equilibrium dissipation fluctuation relation :eq:`Feq` for Bosons (``sign=+1``) or Fermions (``sign=-1``), with :math:`A(\omega)` from :math:`G^R(t)`.

.. _PNess04:

Diagram utilities
-----------------

.. _PNess04S01:

Bubble diagrams
~~~~~~~~~~~~~~~~

Apart from algebraic manipulations there is also a steady state implementation of the bubble products :math:`C(t)=iA_{a,a'}(t) B_{b',b}(-t)` and :math:`C(t)=iA_{a,a'}(t) B_{b,b'}(t)`, see the Diagrams section in the main code. Applying the Langreth rules we can provide the two basic functions, ``ness2::Bubble1_ness`` and ``ness2::Bubble2_ness``, as in the main code, which work also in the namespace ``ness``.

.. list-table::
   :header-rows: 0

   * - ``ness2::Bubble1_ness(herm_matrix_ness &C,int c1,int c2,herm_matrix_ness &A,int a1,int a2,herm_matrix_ness &B,int b1,int b2)``
     - :math:`C_{c1,c2}(t)=i A_{a1,a2}(t)  B_{b2,b1}(-t)` is calculated for given ``herm_matrix_ness`` objects ``A,B`` with matrix indices ``a1``, ``a2`` and ``b1``, ``b2``, respectively.
   * - ``ness2::Bubble1_ness(herm_matrix_ness &C,herm_matrix_ness &A,herm_matrix_ness &B)``
     - :math:`C(t)=i A(t)  B(-t)` is calculated for given scalar ``herm_matrix_ness`` objects ``A,B`` with hermitian symmetry (indices default to zero).
   * - ``ness2::Bubble2_ness(herm_matrix_ness &C,int c1,int c2,herm_matrix_ness &A,int a1,int a2,herm_matrix_ness &B,int b1,int b2)``
     - :math:`C_{c1,c2}(t)=i A_{a1,a2}(t)  B_{b1,b2}(t)` is calculated for given ``herm_matrix_ness`` objects ``A,B`` with matrix indices ``a1``, ``a2`` and ``b1``, ``b2``, respectively.
   * - ``ness2::Bubble2_ness(herm_matrix_ness &C,herm_matrix_ness &A,herm_matrix_ness &B)``
     - :math:`C(t)=i A(t)  B(t)` is calculated for given scalar ``herm_matrix_ness`` objects ``A,B`` with hermitian symmetry (indices default to zero).

**Example:**

The following example calculates the second-order self-energy for a general time-dependent interaction :math:`U(t)` out of a scalar local Green's function, :math:`\Sigma(t,t') = U(t) G(t,t') G(t',t) U(t') G(t,t')`. The diagram is factorized in a particle-hole bubble :math:`\chi` and a particle-particle bubble:

.. math::

   \Sigma(t,t') = -i W(t,t') G(t,t'),
   \quad
   W(t,t')=U(t) \chi(t,t') U(t'),
   \quad
   \chi(t,t') = iG(t,t') G(t',t).

.. code-block:: cpp

   int N_grid = 100;
   int size = 1;
   double h = 0.1;
   int sign= -1;
   cdmatrix U(size,size);
   herm_matrix_ness G(N_grid, size);
   herm_matrix_ness Sigma(N_grid, size);
   herm_matrix_ness W(N_grid, size);

   // set G and U ...

   // compute Sigma
   ness2::Bubble1_ness(W,G,G); // W(t) is set to chi=ii*G(t)G(-t);
   W.left_multiply(U, ness2::fft_domain::time); // set W  to W(t)=U chi(t) U:
   W.right_multiply(U, ness2::fft_domain::time);
   ness2::Bubble2_ness(Sigma,G,W); // Sigma(t) is set to Sigma=ii*G(t)W(t);
   Sigma.smul(-1.0, ness2::fft_domain::time); // the final -1 sign

.. _PNess05:

Dyson equation
--------------

At present, for the solution of the Dyson equation we provide only one method ``METHOD=FFT_TRAPEZ``, which is based on a straightforward discrete Fourier transform: The self energy is transformed to the frequency grid :eq:`w-grid` using the integral transform ``ness2::herm_matrix_ness::integral_transform_to_freq()``, corresponding to a trapezoidal evaluation of :eq:`fft2` and :eq:`fft4`. Then the frequency-dependent functions :math:`G^{<|R}(\omega)` are computed on the grid :eq:`w-grid` using :eq:`nessdyson1` and :eq:`nessdyson2`. For the back-transform to :eq:`fft1` and :eq:`fft3`, we again use the implementation ``ness2::herm_matrix_ness::integral_transform_to_time()``.

.. note::

   In order to reduce Fourier artifacts from the large frequency region we only use :math:`G^{<|R}(\omega)` at frequency grid points :math:`-(N_{\rm freq}/2-1),...,N_{\rm freq}/2-1`, with :math:`N_{\rm freq}=N_{\rm ft}/3`, and set :math:`G^{<|R}(\omega)` to zero outside this interval. With the stepsize :math:`\Delta_\omega`, the maximum frequency is therefore :math:`\omega_{\rm max} = \frac{2\pi}{6h}`. An accurate solution of the Dyson equation requires :math:`h` to be small enough such that spectral functions are sufficiently decayed outside :math:`[-\omega_{\rm max},\omega_{\rm max}]`. Moreover, :math:`N_{\rm ft}` must be sufficiently large such that the functions :math:`G^{R|<}(t)` decay at the largest time :math:`t_{\rm max} = \frac{hN_{\rm ft}}{2}`, and :math:`\Delta_\omega` can resolve the most narrow structures in frequency.

.. list-table::
   :header-rows: 0

   * - ``ness2::dyson(herm_matrix_ness &G, double mu, cdmatrix &epsilon, herm_matrix_ness &Sigma, double h, fft_integral_method method=FFT_TRAPEZ, [ETA])``
     - Solve the frequency Dyson equation with timestep ``h`` for ``G``, for a given self-energy ``Sigma``, local Hamiltonian ``epsilon``, and further parameters ``[ETA]``. For ``METHOD``, the only currently implemented option is ``FFT_TRAPEZ``.

- ``G`` and ``Sigma`` are the Green's functions of type ``ness2::herm_matrix_ness`` which must be defined on the same gridsize ``Nft`` and must have the same matrix dimension.
- ``epsilon`` is a complex Eigen matrix.
- ``[ETA]`` summarizes further optional parameters which introduce a long-time regularization of the Dyson equation. Possible options are ``REG_CONST``, ``REG_GAUSS``, ``REG_OHMIC``
- This routine is also implemented in the namespace ``ness`` using the function call ``ness::dyson_ness(GF & G, GF & Sigma, const cdmatrix & H, double eta)``.

**Example:**

.. code-block:: cpp

   int N_grid = 100;
   int size = 1;
   double h = 0.1;
   double mu=0;
   cdmatrix Ham(size,size);
   herm_matrix_ness G(N_grid, size);
   herm_matrix_ness Sigma(N_grid, size);

   // Initialize G, Sigma, Ham, ...

   ness2::dyson(G, mu, Ham, Sigma, h, FFT_TRAPEZ);

**VIE2**

In the old namespace ``ness``, there is also an explicit implementation for solving a VIE2, the involved Green's functions need to match in their ``size1_``, ``size2_`` and ``ngrid_`` properties and need to be hermitian.

.. list-table::
   :header-rows: 0

   * - ``ness::vie2_ness(GF & G, const GF & F, const GF & Q)``
     - :math:`(1 + F) \ast G = Q`, ``F``, ``Q`` given ``GF`` objects in frequency space, ``G`` is output ``GF`` in frequency space.
   * - ``ness::vie2_ness(GF_pair & G, const GF_pair & F, const GF_pair & Q)``
     - :math:`(1 + F) \ast G = Q`, ``F``, ``Q`` given ``GF_pair`` objects in time, ``G`` is output ``GF_pair`` in time.

.. _PNess06:

Parallelization with OMP
-------------------------

.. list-table::
   :header-rows: 0

   * - class
     - ``ness2::FFT_OMP_Manager``

The Fourier transforms can be parallelized using a shared-memory parallelization, using the built in functionalities of the FFTW library. In order to make this functionality available, we provide a small helper class ``ness2::FFT_OMP_Manager``. ``ness2::FFT_OMP_Manager::initialize`` must be called for initialization at the beginning of the program. Before creating an ``ness2::fft_array`` or ``ness2::herm_matrix_ness`` object (which contains an ``fft_plan`` of the FFTW library), one needs to call ``ness2::FFT_OMP_Manager::set_threads_for_new_plans``. All plans which are created after this call will be generated such that their execution spans over a specified number of ``nomp`` threads. Note that new plans are also created once a new ``ness2::fft_array`` or ``ness2::herm_matrix_ness`` is created via an assignment and or copy assignment. At the end of the program, call ``ness2::FFT_OMP_Manager::finalize``.

.. list-table::
   :header-rows: 0

   * - ``ness2::FFT_OMP_Manager::initialize()``
     - Initialize parallelization helper class, call at beginning of program.
   * - ``ness2::FFT_OMP_Manager::set_threads_for_new_plans(int nomp)``
     - Fix amount of threads ``nomp`` for all plans to be created (in ``ness2::fft_array`` or ``ness2::herm_matrix_ness``).
   * - ``ness2::FFT_OMP_Manager::finalize()``
     - Finalize, call at end of program.

Inside the code, you can access the number of threads currently used in the construction of plans from the variable ``FFT_OMP_Manager::current_threads``. To test the parallel setup, we provide a notebook ``utils/test_ness2_omp.ipynb``.

.. note::

   - This functionality is added automatically when the library is built with ``omp=ON`` and ``ness=ON``. In addition it requires the FFTW library to be compiled with ``--enable-threads``. If only the non-threaded FFTW library is available (or if ``omp=OFF``), the calls to the FFT OMP Manager have simply no effect, and all Fourier transforms are single threaded.
   - The generation of plans (i.e., the creation of any ``fft_herm_matrix_ness`` or ``fft_array``) is not thread-safe, and should be called from a serial region of the code.
   - In general, one should avoid nesting outer parallelization and parallelization of the Fourier transforms. Using outer parallelization in combination with FFT plans with more than one thread may lead to oversubscription. If outer parallelization is used, plans should therefore simply be created with ``nomp=1`` threads.

.. _PNess07:

Utilities
---------

.. _PNess07S01:

Comparing Green's functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To compare the data of two Green's functions, e.g. for convergence checks, one can analyze the difference using the :math:`L_2`-norm difference :math:`|| M ||` for :math:`M=A-B` of the individual elements. The difference is defined as the :math:`L_2`-norm difference, summed over all time or frequency points:

.. math::
   :label: eq:ness_distnorm

   \Delta[A,B]
   =
   \sum_{i=0}^{\tt ngrid-1}
   \big(
   ||A^\mathrm{R}(i \cdot \Delta x)-B^\mathrm{R}(i \cdot \Delta x)||
   +
   ||A^<(i\cdot \delta x)-B^<(i\cdot \delta x)||
   \big).

:math:`|| M ||` is the standard matrix :math:`L_2` norm, :math:`\Delta x = \Delta \omega` for ``domain=ness2::fft_domain::freq`` and :math:`\Delta x = h` for ``domain=ness2::fft_domain::time``. The corresponding function call is ``ness2::distance_norm2``:

.. list-table::
   :header-rows: 0

   * - ``ness2::distance_norm2(herm_matrix_ness &A,herm_matrix_ness &B,domain)``
     - Returns difference norm :math:`\|A-B\|_2` (double) for ``ness2::herm_matrix_ness`` or ``ness2::fft_array`` objects, on the time or frequency grid (``domain=ness2::fft_domain::[time|freq]``).

In the old namespace the function call is ``double ness::GF2norm(GF &g1, GF &g2)``, the type of ``g1`` and ``g2`` has to match and they must have equal ``size1_``, ``size2_`` and ``ngrid_``.

**Example:**

.. code-block:: cpp

   int N_grid = 100;
   int size = 1;
   double h = 0.1;
   herm_matrix_ness G(N_grid, size);

   // some iterative procedure to determine G
   {
       herm_matrix_ness tG(N_grid, size);   //temporary variable
       double convergence_error;
       int iter_max=100;
       for(int iter=0;iter<=iter_max;iter++){
           tG = G;        // store values of G before iteration
           // ... some code to update G ...
           convergence_error = ness2::distance_norm2(tG,G,ness2::fft_domain::time);
           if( convergence_error < some_sufficiently_small_number) break;
       }
       if(iter>iter_max){
           cout << "no convergence!" << endl;
       }
   }

.. _PNess07S02:

Data exchange with two-time functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Data exchange between two-time and steady state functions is enabled by the following functions:
The function ``ness2::cntr2ness`` allows to read the ``ness2::fft_domain::time`` data of a steady state function :math:`G_{\rm ness}` from a given timeslice :math:`t_0` of a two-time object :math:`G`, such that :math:`G_{\rm ness}^{R|<}(t)\leftarrow G_{\rm cntr}^{R|<}(t_0,t_0-t)` for all :math:`t`.
At negative times :math:`t`, :math:`\bar G` is set assuming hermitian symmetry; if ``tstp`` is smaller than the maximum time :math:`N_{\rm ft}/2-1` in :math:`G_{\rm ness}`, the remaining entries in :math:`G_{\rm ness}` are left zero.
The reverse function ``ness2::ness2cntr`` sets :math:`G_{\rm cntr}^{R|<}(t,t')\leftarrow G_{\rm ness}^{R|<}(t-t')` for all arguments :math:`(t,t')` where :math:`t-t'` is in the domain of :math:`G_{\rm ness}`, and zero otherwise. The two function calls are:

.. list-table::
   :header-rows: 0

   * - ``ness2::cntr2ness(herm_matrix_ness& Gness, cntr::herm_matrix<T>& Gcntr, int tstp = -1)``
     - Set values :math:`G^{<|R}` of a ``ness2::herm_matrix_ness`` object ``Gness`` from timeslice ``tstp`` of a two-time function ``Gcntr``.
   * - ``ness2::ness2cntr(cntr::herm_matrix<T>& Gcntr, const herm_matrix_ness& Gness)``
     - Set :math:`G^{<|R}` of a two-time function ``Gcntr`` from a steady state function ``Gness``, assuming time-translational invariance.

- ``Gcntr`` can be ``cntr::herm_matrix``, ``cntr::herm_matrix_timestep``, ``cntr::herm_matrix_moving`` or ``cntr::herm_matrix_timestep_moving``.
- For the latter three types, the argument ``tstp`` in ``cntr2ness`` is ignored and can be omitted, because always the leading physical timestep is addressed.
- If the argument ``tstp`` is omitted for ``Gcntr`` of type ``cntr::herm_matrix``, it defaults to the largest physical time (``Gcntr.nt``).

.. _PNess07S03:

Resampling
~~~~~~~~~~~

In the steady state, calculations on much finer time grids are possible than in two-time calculations. Therefore downsampling and upsampling routines are necessary to read and write Green's functions defined on different grids effectively between the two interfaces. The following routines resample steady state ``ness2::herm_matrix_ness`` objects, so that they can be brought into the same shape as their ``cntr`` counterpart and be processed with the ``ness2::cntr2ness`` and ``ness2::ness2cntr`` functions.

**Upsampling:**

Upsampling a ``ness2::herm_matrix_ness`` object to one with finer time grid is equivalent to increasing the maximum frequency on frequency domain. One can thus transform the ``ness2::herm_matrix_ness`` object of length ``Nft`` to frequency using ``integral_transform_to_freq`` and enlarge the frequency interval by adding zeros until the desired length ``Nft*factor`` is achieved. Then a back-transform to time using ``integral_transform_to_time`` on the new grid ``Nft*factor`` returns the upsampled object in time.

.. list-table::
   :header-rows: 0

   * - ``ness2::upsample(double h_in,herm_matrix_ness &in, int factor)``
     - Upsample a ``ness2::herm_matrix_ness`` object ``in`` with initial timestep ``h_in`` to a finer grid with a factor of ``factor`` more points and spacing ``h_in/factor``.

**Downsampling:**

Downsampling a ``ness2::herm_matrix_ness`` object to one with coarser time grid is equivalent to keeping only every ``factor``-th element in a Green's function of length ``Nft`` in time domain to reduce it to the length ``Nft/factor``. In frequency space this corresponds to dropping all high frequency elements outside the length ``Nft/factor`` of the resampled object. Returns ``ness2::herm_matrix_ness`` or ``ness2::fft_array``.

.. list-table::
   :header-rows: 0

   * - ``ness2::downsample(herm_matrix_ness &in, int factor)``
     - Downsample a ``ness2::herm_matrix_ness`` or ``ness2::fft_array`` object ``in`` to a coarser grid with a factor of ``factor`` fewer points. Returns ``ness2::herm_matrix_ness`` or ``ness2::fft_array``.

.. _PNess07S04:

HDF5 python utilities
~~~~~~~~~~~~~~~~~~~~~~

Also for the steady state code we provide python tools, which include scripts for reading and post-processing Green's functions from HDF5 format via the ``h5py`` python package.
The python tools can be found in ``libcntr/python3`` in the file ``ReadNESS.py``. An example usage is demonstrated in the example programs :ref:`NessEx01` and :ref:`NessEx02`.

.. note::

   To make the python module available to a python script, the module needs to be in the python include path. This is achieved by setting ``export PYTHONPATH=<path>/nessi/libcntr/python3``.
