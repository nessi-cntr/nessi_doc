.. _PMan12:

Free Green's Functions
======================

.. contents::
   :local:
   :depth: 2

We provide a number of tools for calculating Green's functions for non-interacting systems, either from a single-particle Hamiltonian or from a given spectral function.

.. _PMan12S01:

Equilibrium Green's functions for a given density of states
------------------------------------------------------------

.. _dos_def:

In many applications (for representing a fermionic bath, for example), one needs to construct Green's functions with a given density of states (DoS) :math:`A(\omega)`:

.. math::

   G(t,t') = \int d \omega \,A(\omega) g_{\omega}(t,t') \ ,

where :math:`g_\omega(t,t')` is the single-orbital noninteracting Green's function with energy :math:`\omega`. The latter is defined as (:math:`F_\xi(x)=1/(e^{\beta x}-\xi)` is the Fermi/Bose distribution):

- :math:`g_\omega^M (\tau) = - e^{-\tau(\omega-\mu)} F_\xi(-(\omega-\mu))`
- :math:`g_\omega^R (t,t') = -i e^{-i(t-t')\omega}`
- :math:`g_\omega^< (t,t') = -i\xi  e^{-i(t-t')\omega} F_\xi(\omega-\mu)`
- :math:`g_\omega^{\rceil} (t,\tau) = -i\xi e^{-it\omega} e^{\tau(\omega-\mu)} F_\xi(\omega-\mu)`

This task is accomplished by ``cntr::green_equilibrium``:

.. list-table::
   :header-rows: 0

   * - ``cntr::green_equilibrium(herm_matrix<T> &G, dos_function &dos, double beta, double h, double mu=0.0, int limit = 100, int nn = 20)``
     - Computes the free Green's function with given density of states ``dos``. **output:** object ``G`` of type ``cntr::herm_matrix``

The user provides the inverse temperature ``beta``, the time step ``h`` and the chemical potential ``mu`` (set to 0.0 by default). Parameters ``limit`` and ``nn`` determine the adaptive sampling method implemented in ``fourier``. Use the default settings unless there is a special reason not to. ``dos`` is an object representing a DoS class. This can represent an arbitrary density of states. As a useful example to show the minimal call properties, we provide the Bethe DoS :math:`A(\omega) = V^2/(2\pi) \sqrt{ 4 V^2 - \omega^2}`, represented by the class:

.. code-block:: cpp

   class bethedos {
     public:
       double hi_; // Higher edge of the density of states
       double lo_; // Lower edge of the density of states
       double V_; // 4V corresponds to the bandwidth
       bethedos() {
           V_ = 1;
           lo_ = -2;
           hi_ = 2;
       }
       double operator()(double x) {
           double arg = 4.0 * V_ * V_ - x * x;
           double num = V_ * V_ * 3.14159265358979323846 * 2;
           return (arg < 0 ? 0.0 : sqrt(arg) / num);
       }
   };

The member variables ``lo_``, ``hi_`` and the function ``double operator()(double x)`` are the basic building blocks for any DoS class. For constructing a Green's function for a Bethe DoS with band width of 1 centered around :math:`\omega=0`, we would use the following code:

.. code-block:: cpp

   // define parameters, nt, ntau, beta,  ...
   // define herm_matrix G

   cntr::bethedos dos;
   dos.V_ = 0.25;
   dos.lo_ = -0.5;
   dos.hi_ = 0.5;

   cntr::green_equilibrium(G, dos, beta, h);

We also provide the function ``cntr::green_equilibrium_bethe``, which constructs the Green's function for the Bethe DoS with ``V_=1.0`` without invoking the ``dos`` class directly.

Another example of a useful DoS is a smoothened box-like shape :math:`A(\omega) = C (f_\nu(\omega-b) - f_\nu(\omega-a))`, where :math:`f_\nu(\omega)` is a Fermi function with effective temperature :math:`\nu`, while :math:`[a, b]` defines the interval for the DoS. Further details can be found in the documentation of ``cntr::smooth_box``.

.. _PMan12S02:

Green's function for a constant or time-dependent Hamiltonian
--------------------------------------------------------------

Another often encountered example is the Green's function of a single-particle Hamiltonian :math:`\epsilon(t)`, defined by:

.. math::

   \label{tttsfwg}
   \big[ i \partial_t + \mu - \epsilon(t) \big] G(t,t')= \delta_{\mathcal{C}}(t,t')  \ .

We provide routines for both constant and time-dependent Hamiltonians:

.. list-table::
   :header-rows: 0

   * - ``cntr::green_from_H(herm_matrix<T> &G,T mu,cdmatrix &eps,T beta,T h)``
     - Computes the free Green's function for constant Hamiltonian ``eps``. **output:** object ``G`` of type ``cntr::herm_matrix``
   * - ``cntr::green_from_H(int tstp, GG &G,T mu,cdmatrix &eps,T beta,T h)``
     - Computes the free Green's function for constant Hamiltonian ``eps`` at time step ``tstp``. **output:** object ``G``
   * - ``cntr::green_from_H(herm_matrix<T> &G,T mu,cntr::function<T> &eps,T beta,T h,bool fixHam,int SolveOrder,int cf_order)``
     - Computes the free Green's function for the time-dependent Hamiltonian ``eps``. **output:** object ``G`` of type ``cntr::herm_matrix``
   * - ``cntr::green_from_H(herm_matrix<T> &G,T mu,cntr::function<T> &eps,T beta,T h,int SolveOrder,int cf_order)``
     - Computes the free Green's function for the time-dependent Hamiltonian ``eps``. **output:** object ``G`` of type ``cntr::herm_matrix``
   * - ``cntr::green_from_H(int tstp, GG &G,T mu,cntr::function<T> &eps,T beta,T h,bool fixHam,int SolveOrder,int cf_order)``
     - Computes the free Green's function for the time-dependent Hamiltonian ``eps`` at time step ``tstp``. **output:** object ``G``

- In the routines with a ``tstp``-argument, ``GG`` can be ``herm_matrix_timestep<T>`` or ``herm_matrix<T>``
- For a time-dependent Hamiltonian, ``eps`` stores the Hamiltonian at all times. The Green's function is obtained by numerically solving the equation of motion. The numerical parameters (``cf_order``, ``SolveOrder``, ``fixHam``) are explained in the second example below.

**Example 1:** The Green's function for a time-independent Hamiltonian can be constructed as follows:

.. code-block:: cpp

   int nt=100;
   int ntau=100;
   int size1=2;
   double mu=0.0;
   double h=0.05;
   double beta=10.0;

   GREEN G(nt,ntau,size1,FERMION);
   cdmatrix eps0(size1,size1);

   // define the Hamiltonian
   eps0(0,0) = -1.0;
   eps0(0,1) = 0.5 + 0.5*II;
   eps0(1,0) = std::conj(eps(0,1));
   eps0(1,1) = 1.0;

   // construct GF
   cntr:::green_from_H(G, mu, eps0, beta, h);

**Example 2:** For time-dependent Hamiltonians, we have implemented the commutator-free decomposition of the time-evolution operator in second (``cf_order=2``) and fourth (``cf_order=4``) order. By default, the fourth-order scheme is used. The method requires the evaluation of :math:`\epsilon(t)` at times between the interval :math:`[(n-1) h, n h]`; thus, we use polynomial interpolation of order ``SolveOrder``. By default, the maximum order is used. The following example shows how to use ``cntr::green_from_H`` for this case:

.. code-block:: cpp

   // define nt, ntau, h, beta, mu ...
   int size1=2;

   cdmatrix eps0(size1,size1);
   // defines eps0 as in the above example
   cdmatrix eps1(size1,size1);
   eps1(0,1) = 1.0;
   eps1(1,0) = 1.0;

   GREEN G(nt,ntau,size1,FERMION);
   CFUNC eps(nt, size1);

   eps.set_value(-1, eps0);

   cdmatrix eps_at_t(size1,size1);
   for(int tstp=0; tstp<=nt; tstp++){
       eps_at_t = eps0 + sin(tstp*h)*eps1;
       eps.set_value(tstp, eps_at_t);
   }

   cntr::green_from_H(G, mu, eps, beta, h, true);

A typical example where the Hamiltonian is not yet known at a certain time step is a time-dependent mean-field calculation, where the density matrix at :math:`t_n = nh` is needed to construct :math:`\epsilon(n h)`. For this purpose, we have implemented a polynomial extrapolation as a predictor for :math:`\epsilon(n h)`, triggered by the argument ``fixHam=False``. In contrast, ``fixHam=True`` indicates that :math:`\epsilon(n h)` is already known and an interpolation in :math:`[(n-1) h, n h]` can be employed. The following code snippet exemplifies a corresponding predictor-corrector scheme:

.. code-block:: cpp

   // define nt, ntau, h, beta, mu ...
   // define SolveOrder

   GREEN G(nt,ntau,size1,FERMION);
   CFUNC eps(nt, size1);

   // determine eps, G for time steps tstp=-1,...,SolveOrder

   for(int tstp=SolveOrder+1; tstp<=nt; tstp++){
       // predictor
       cntr::green_from_H(tstp, G, mu, eps, beta, h, false, SolveOrder);

       // update density from G at tstp, update eps at tstp

       // corrector
       cntr::green_from_H(tstp, G, mu, eps, beta, h, true, SolveOrder);
   }

.. _PMan12S03:

Bosonic Green's functions
--------------------------

All the above functions work for fermionic or bosonic Green's functions. However, for problems with electron-phonon coupling (or similar effects), we provide functions for calculating the free phonon propagator :math:`D_0(t,t') = -i\langle T_\mathcal{C} \hat{X}(t)\hat{X}(t')\rangle`, with :math:`\hat{X} = (\hat{a} + \hat{a}^\dagger)/\sqrt{2}`:

.. list-table::
   :header-rows: 0

   * - ``cntr::green_single_pole_XX_timestep(int tstp, herm_matrix_timestep<T> &D0, double w, double beta, double h)``
     - Computes the free phonon-type propagator with frequency ``w`` at time step ``tstp``. **output:** object ``D0`` of type ``cntr::herm_matrix_timestep``
   * - ``cntr::green_single_pole_XX_timestep(int tstp, herm_matrix<T> &D0, double w, double beta, double h)``
     - Computes the free phonon-type propagator with frequency ``w`` at time step ``tstp``. **output:** object ``D0`` of type ``cntr::herm_matrix``
   * - ``cntr::green_single_pole_XX_timestep(herm_matrix<T> &D0, double w, double beta, double h)``
     - Computes the free phonon-type propagator with frequency ``w`` for all time steps. **output:** object ``D0`` of type ``cntr::herm_matrix``
