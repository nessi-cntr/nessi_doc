.. _PMan07:

Dyson Equation
==============

.. contents::
   :local:
   :depth: 2

The Dyson equation for the Green's function :math:`G(t,t^\prime)` can be written as:

.. math::
   :label: dyson

   i\partial_t G(t,t^\prime) + \mu G(t,t^\prime)  - \epsilon(t) G(t,t^\prime) -
   \int_\mathcal{C} d\bar t\, \Sigma(t,\bar t) G(\bar t,t^\prime) = \delta_{\mathcal{C}}(t,t^\prime).

This equation is to be solved for :math:`G(t,t^\prime)` for given input :math:`\epsilon(t)` and :math:`\Sigma(t,t^\prime)`, and the KMS boundary conditions. All quantities :math:`\Sigma(t,t^\prime)`, :math:`G(t,t^\prime)`, and :math:`\epsilon(t)` are square matrices. The equation is an integro-differential form of the Dyson series:

.. math::

   G=G_0+G_0\ast \Sigma\ast G,

where the free Green's function :math:`G_0` is determined by the differential equation :math:`i\partial_t G(t,t^\prime) + \mu G(t,t^\prime)  - \epsilon(t) G(t,t^\prime)  = \delta_{\mathcal{C}}(t,t^\prime)` (cf Sec. :ref:`PMan12`).

It is assumed that :math:`\Sigma=\Sigma^\ddagger` is hermitian, and :math:`\epsilon(t) =\epsilon(t)^\dagger`, which implies that also the solution :math:`G` possesses hermitian symmetry. Because of the hermitian symmetry, :math:`G` can also be determined from the equivalent conjugate equation:

.. math::
   :label: dyson_cc

   -i\partial_{t^\prime} G(t,t^\prime) + \mu G(t,t^\prime)-   G(t,t^\prime) \epsilon(t^\prime)-
   \int_\mathcal{C} d\bar t \,G(t,\bar t) \Sigma(\bar t,t^\prime) = \delta_{\mathcal{C}}(t,t^\prime).

Because of its causal nature, the Dyson equation can be solved in a time-stepping manner (see explanation at input/output relation below). The following routines are used to solve the Dyson equation:

.. list-table::
   :header-rows: 0

   * - ``void cntr::dyson(GG &G, T mu, cntr::function<T> &H, GG &Sigma, T beta, T dt, int SolveOrder, MAT_METHOD)``
     - Solve :eq:`dyson` for ``G`` in the full two-time plane with given ``Sigma``, ``mu``, ``H``.
   * - ``void cntr::dyson_mat(GG &G, T mu, cntr::function<T> &H, GG &Sigma, T beta, int SolveOrder, MAT_METHOD)``
     - Solve :eq:`dyson` for ``G`` on timestep ``tstp=-1`` with given ``Sigma``, ``mu``, ``H``.
   * - ``void cntr::dyson_start(herm_matrix<T> &G, T mu, cntr::function<T> &H,GG &Sigma, T beta, T dt, int SolveOrder)``
     - Solve :eq:`dyson` for ``G`` of type ``cntr::herm_matrix`` on all timesteps ``0 <= tstp <= SolveOrder`` with given ``Sigma``, ``mu``, ``H``
   * - ``void cntr::dyson_timestep(int tstp,GG &G, T mu, cntr::function<T> &H,GG &Sigma, T beta, T dt, int SolveOrder)``
     - Solve :eq:`dyson` for ``G`` on a timestep ``tstp`` with ``tstp > SolveOrder`` with given ``Sigma``, ``mu``, ``H``
   * - ``void cntr::dyson_timestep_omp(int omp_num_threads, int tstp,GG &G, T mu, cntr::function<T> &H,GG &Sigma, T beta, T dt, int SolveOrder)``
     - Same as ``cntr::dyson_timestep``, using shared memory parallelization. **output:** object ``G``

- ``G`` and ``Sigma`` are Green's functions of type ``cntr::herm_matrix<T>``, assumed to be hermitian
- ``H`` is a ``cntr::function<T>``, the function :math:`\epsilon` in :eq:`dyson`.
- ``SolveOrder`` :math:`\in` ``1,...,MAX_SOLVE_ORDER``, the order of accuracy for the solution. Use ``SolveOrder=5`` (``=MAX_SOLVE_ORDER``) if there is no good reason against it.
- ``dt`` is the time-discretization step :math:`\Delta t` on the real-time branch
- ``beta`` is the length of the imaginary-time contour (inverse temperature)
- Size requirements:

  - ``G.size1()==Sigma.size1()==H.size1()``, ``G.ntau()==Sigma.ntau()``, ``G.ntau() > SolveOrder``
  - For ``dyson_start``: ``G.nt(),Sigma.nt() >= SolveOrder``
  - For ``dyson_timestep``: ``G.nt(),Sigma.nt() >= tstp``

- Special numerical parameters for ``dyson_mat``: ``MAT_METHOD`` = ``CNTR_MAT_FOURIER`` or ``CNTR_MAT_FIXPOINT``. Different methods to solve the equation on the Matsubara contour. ``CNTR_MAT_FIXPOINT`` is default and should be used.

**Input/Output relation and time-stepping:**

- ``dyson_mat``: ``Sigma`` is read on timestep ``tstp=-1``, ``H`` is read on times ``tstp=-1``; ``G`` is written on timestep ``tstp=-1``
- ``dyson_start``: ``Sigma`` is read on timestep ``tstp=-1,...,SolveOrder``, ``H`` is read on times ``tstp=-1,...,SolveOrder``, ``G`` is read on timestep ``tstp=-1``; ``G`` is written on timestep ``tstp=0,...,SolveOrder``
- ``dyson_timestep``: ``Sigma`` is read on timestep ``t=-1,...,tstp``, ``H`` is read on times ``t=-1,...,tstp``, ``G`` is read on timestep ``t=-1,...tstp-1``; ``G`` is written on timestep ``t=tstp``

Because of this causal structure, the Dyson equation can be solved by time-stepping: First ``G`` is determined on timestep ``tstp=-1`` using ``dyson_mat``. The result enters the determination of ``G`` on timesteps ``0,...,SolveOrder``, done with ``dyson_start``. The equation is solved successively for timesteps ``tstp=SolveOrder+1,SolveOrder+2,...``, where the result at timestep ``tstp`` depends on all previous timesteps.

**Example:**

Solution of :eq:`dyson` on all times for given self-energy:

.. code-block:: cpp

   int nt=10;
   int ntau=20;
   int size1=2;
   int SolveOrder=5;
   double dt=0.01; // time-discretization
   double beta=10.0; // inverse temperature
   double mu=1.433; // chemical potential
   GREEN G(nt,ntau,size1,FERMION);
   GREEN Sigma(nt,ntau,size1,FERMION);
   CFUNC H(nt,size1);
   //...
   // ... do something to set H and Sigma on timestep -1
   // solve for G on timestep -1:
   cntr::dyson_mat(G,mu,H,Sigma,beta,SolveOrder,CNTR_MAT_FIXPOINT);
   // ... do something to set H and Sigma on timesteps 0 ... SolveOrder=5
   // solve for G on timesteps 0 ... SolveOrder=5
   cntr::dyson_start(G,mu,H,Sigma,beta,dt,SolveOrder);
   //all other timesteps:
   for(int tstp=SolveOrder+1;tstp<=nt;tstp++){
       // ... do something to set H and Sigma on timestep tstp
       cntr::dyson_timestep(tstp,G,mu,H,Sigma,beta,dt,SolveOrder);
   }

**OpenMP parallelization:**

The variant ``dyson_timestep_omp`` uses shared memory ``openMP`` parallelization to distribute the evaluation of :math:`G(t,t')` for different time arguments of the timestep ``tstp`` over ``omp_num_threads`` tasks. Note that ``dyson_timestep_omp`` and ``dyson_timestep`` follow a slightly different implementation, thus the result differs by the numerical error.

The ``openMP`` parallelization is handled internally in ``dyson_timestep_omp``. A race condition cannot occur if ``omp_num_threads`` is set to the number of ``openMP`` threads in the current team.
