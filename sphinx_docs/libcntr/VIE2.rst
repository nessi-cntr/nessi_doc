.. _PMan08:

VIE2
====

.. contents::
   :local:
   :depth: 2

.. _PMan08S01:

Non-singular VIE2
-----------------

The second important equation is an integral equation of the form:

.. math::
   :label: vie2

   G(t,t^\prime) + \int_\mathcal{C} d\bar t\, F(t,\bar t) G(\bar t,t^\prime) = Q(t,t^\prime) \quad \Leftrightarrow\quad (1+F)*G=Q,

This linear equation is to be solved for :math:`G(t,t^\prime)` for a given input kernel :math:`F(t,t^\prime)`, its hermitian conjugate :math:`F^\ddagger(t,t^\prime)`, and a source term :math:`Q(t,t^\prime)`. A typical physical application of :eq:`vie2` is the summation of a random phase approximation (RPA) series for a susceptibility :math:`\chi`,

.. _rpa_def:

.. math::
   :label: eq:rpa_vie2

   \chi
   = \chi_0 + \chi_0\ast V \ast \chi_0 +  \chi_0\ast V \ast \chi_0\ast V \ast \chi_0 + \cdots
   = \chi_0 + \chi_0\ast V \ast \chi.

In the solution of the equation, we assume that both :math:`Q` and :math:`G` are hermitian. In general, the hermitian symmetry would not hold for an arbitrary input :math:`F` and :math:`Q`. However, it does hold when :math:`F` and :math:`Q` satisfy the relation :math:`F\ast Q=Q\ast F^\ddagger` and :math:`Q=Q^\ddagger`, which is the case for the typical applications such as the RPA series. In this case, there is an equivalent conjugate equation:

.. math::
   :label: vie2_cc

    G(t,t^\prime) + \int_\mathcal{C} d\bar t \,G(t,\bar t) F^\ddagger(\bar t,t^\prime) = Q(t,t^\prime) \quad \Leftrightarrow\quad G*(1+F^\ddagger)=Q.

In the implementation, the solution of Eq :eq:`vie2` is reduced to a **V**\olterra **I**\ntegral **E**\quation of **2**\nd kind,

.. _vie2_def:

which is why the corresponding routines are called ``vie2``.

Because of its causal nature, the VIE2 equation can be solved in a time-stepping manner (see explanation of the input/output relation below). The following routines are used to solve the VIE2 equation:

.. list-table::
   :header-rows: 0

   * - ``void cntr::vie2_mat(GG &G, GG &F, GG &Fcc, GG &Q, T beta, MAT_METHOD, int SolveOrder)``
     - Solve :eq:`vie2` for ``G`` on timestep ``tstp=-1`` with given two-time objects ``F`` and ``Q``
   * - ``void cntr::vie2_start(GG &G, GG &F, GG &Fcc, GG &Q,T beta,T dt,int SolveOrder)``
     - Solve :eq:`vie2` for ``G`` on all timesteps ``0 <= tstp <= SolveOrder`` with given two-time objects ``F`` and ``Q``
   * - ``void cntr::vie2_timestep(int tstp,GG &G, GG &F, GG &Fcc, GG &Q,T beta,T dt,int SolveOrder)``
     - Solve :eq:`vie2` for ``G`` on a timestep ``tstp`` with ``tstp > SolveOrder`` with given two-time objects ``F`` and ``Q``
   * - ``void cntr::vie2_timestep_omp(int omp_num_threads,int tstp,GG &G, GG &F, GG &Fcc, GG &Q,T beta,T dt,int SolveOrder)``
     - Same as ``cntr::vie2_timestep``, using shared memory parallelization. **output:** object ``G``

- The type ``GG`` of ``G``, ``F``, ``Fcc``, and ``Q`` is ``cntr::herm_matrix<T>``. ``F`` and ``Fcc`` store the hermitian domain of :math:`F` and :math:`F^\ddagger`, respectively.
- ``SolveOrder`` :math:`\in` ``1,...,MAX_SOLVE_ORDER``, the order of accuracy for the solution. Use ``SolveOrder=5`` (``=MAX_SOLVE_ORDER``) if there is no good reason against it.
- ``dt`` is the time-discretisation step :math:`\Delta t` on the real-time branch
- ``beta`` is the length of the imaginary-time contour (inverse temperature)
- Size requirements:

  - ``G``, ``F``, ``Fcc``, and ``Q`` must have the same ``size1`` and ``ntau``
  - For ``vie2_start``: ``G``, ``F``, ``Fcc``, and ``Q`` must have ``nt>= SolveOrder``
  - For ``vie2_timestep``: ``G``, ``F``, ``Fcc``, and ``Q`` must have ``nt >= tstp``

- Special numerical parameters for ``vie2_mat``: ``MAT_METHOD`` = ``CNTR_MAT_FOURIER`` or ``CNTR_MAT_FIXPOINT``. ``CNTR_MAT_FIXPOINT`` is default and should be used.

**Input/Output relation and time-stepping:**

- ``vie2_mat``: ``F`` and ``Q`` are read on timestep ``tstp=-1``; ``G`` is written on timestep ``tstp=-1``
- ``vie2_start``: ``F`` and ``Q`` are read on timestep ``tstp=-1,...,SolveOrder``, ``G`` is read on timestep ``tstp=-1``; ``G`` is written on timestep ``tstp=0,...,SolveOrder``
- ``vie2_timestep``: ``F`` and ``Q`` are read on timestep ``t=-1,...,tstp``, ``G`` is read on timestep ``t=-1,...tstp-1``; ``G`` is written on timestep ``t=tstp``

Because of this causal structure, the VIE2 equation can be solved by time-stepping, analogous to ``dyson``: First ``G`` is determined on timestep ``tstp=-1`` using ``vie2_mat``. The result enters the determination of ``G`` on timesteps ``0,...,SolveOrder`` with ``vie2_start``. The equation is solved successively for timesteps ``tstp=SolveOrder+1,SolveOrder+2,...``, where the result at timestep ``tstp`` depends on all previous timesteps.

**Example:**

Solution of the RPA series :math:`\chi = \chi_0 + \chi_0\ast W\ast \chi` for :math:`\chi`, where :math:`\chi_0` is a known (hermitian) two-time function, and :math:`W` is a known one-time function (also hermitian, :math:`W(t)=W(t)^\dagger`). The equation is cast in the form :eq:`vie2` with :math:`F (t,t')= -\chi_0(t,t') W(t')`, :math:`F^\ddagger (t,t')=  -W(t)\chi_0(t,t')`, and :math:`Q=\chi_0`.

.. code-block:: cpp

   int tstp;
   int nt=10;
   int ntau=20;
   int size1=2;
   int SolveOrder=5;
   double dt=0.01; // time-discretization
   double beta=10.0; // inverse temperature
   GREEN chi(nt,ntau,size1,BOSON);
   GREEN chi0(nt,ntau,size1,BOSON);
   CFUNC W(nt,size1);
   GREEN F(nt,ntau,size1,BOSON);
   GREEN Fcc(nt,ntau,size1,BOSON);
   //
   // ... do something to determine chi0 and W on timestep -1
   // determine F, Fcc, and solve for chi on timestep -1:
   tstp=-1;
   F.set_timestep(tstp,chi0);
   F.right_multiply(tstp,W,-1.0);
   Fcc.set_timestep(tstp,chi0);
   Fcc.left_multiply(tstp,W,-1.0);
   cntr::vie2_mat(chi,F,Fcc,chi0,beta,SolveOrder,CNTR_MAT_FIXPOINT);
   // .... do something to determine chi0 and W  on timesteps 0 ... SolveOrder=5
   // get F, Fcc, and solve for chi on timesteps 0 ... SolveOrder=5:
   for(int tstp=0;tstp<=SolveOrder;tstp++){
       F.set_timestep(tstp,chi0);
       F.right_multiply(tstp,W,-1.0);
       Fcc.set_timestep(tstp,chi0);
       Fcc.left_multiply(tstp,W,-1.0);
   }
   cntr::vie2_start(chi,F,Fcc,chi0,beta,dt,SolveOrder);
   //all other timesteps:
   for(int tstp=SolveOrder+1;tstp<=nt;tstp++){
       // .... do something to determine chi0 and W  on timesteps tstp
       F.set_timestep(tstp,chi0);
       F.right_multiply(tstp,W,-1.0);
       Fcc.set_timestep(tstp,chi0);
       Fcc.left_multiply(tstp,W,-1.0);
       cntr::vie2_timestep(tstp,chi,F,Fcc,chi0,beta,dt,SolveOrder);
   }

**OpenMP parallelization:**

The variant ``vie2_timestep_omp`` uses shared memory ``openMP`` parallelization to distribute the evaluation of :math:`G(t,t')` for different time arguments of the timestep ``tstp`` over ``omp_num_threads`` tasks. Note that ``vie2_timestep_omp`` and ``vie2_timestep`` follow a slightly different implementation, thus the result differs by the numerical error. Similar to ``dyson_timestep_omp``, the implementation prevents any race conditions if ``omp_num_threads`` is set to the number of threads in the current team.

.. _PMan08S02:

Singular VIE2
-------------

The singular ``vie2`` is a generalized solver of VIE2, where propagators have an instantaneous term, for example:

.. math::
   :label: Gsin

   G(t,t^\prime) = \delta_{\mathcal{C}}(t,t^\prime) G^\delta(t) + G^R(t,t^\prime),

where :math:`G^\delta` is the instantaneous part and :math:`G^R(t,t^\prime)` is the retarded part. The generalization of Eq. :eq:`vie2` to the following example leads to:

.. math::
   :label: vie2sin

   (1+F^\delta + F^R)*(G^\delta + G^R)=Q^\delta + Q^R

First, for a given timestep :math:`t`, we solve this equation for an instantaneous part leading to:

.. math::

   G^\delta(t)= (1+ F^\delta(t))^{-1} Q^\delta(t) .

The retarded part of the VIE2 can be rewritten in a form where a routine for the non-singular VIE2 can be employed, see :ref:`PMan08S01`:

.. math::

   (1+[F^\delta + F^R])*G^R= Q^R - F^R * G^\delta,

where :math:`G^\delta` is known from the solution of the instantaneous part.

A typical physical usage of the singular version of VIE2 is an evaluation of the retarded interaction within the GW approximation, see :ref:`SecGW` for an example:

.. math::

   [1-V(t) \chi_0(t,t^\prime)]* W(t,t^\prime)= V(t) \delta(t,t^\prime),

where :math:`V` is an instantaneous interaction, :math:`\chi_0` is a bare bubble and :math:`W` is a retarded interaction entering the GW self-energy.

The following routines are used to solve Eq. :eq:`vie2sin`:

.. list-table::
   :header-rows: 0

   * - ``void cntr::vie2_timestep_sin(int tstp,GG &G, cntr::function<T> &Gsin, GG &F, GG &Fcc, cntr::function<T> &Fsin,GG &Q, cntr::function<T> &Qsin,T beta,T dt,int SolveOrder)``
     - Solve :eq:`vie2sin` for ``G`` for a given timestep ``tstp``
   * - ``void cntr::vie2_timestep_sin_omp(int tstp,GG &G, cntr::function<T> &Gsin, GG &F, GG &Fcc, cntr::function<T> &Fsin,GG &Q, cntr::function<T> &Qsin,T beta,T dt,int SolveOrder)``
     - Same as ``cntr::vie2_timestep_sin``, using shared memory parallelization.

These routines use the same input as ``cntr::vie2_timestep``, see :ref:`PMan08S01`, with an addition of a singular part of the propagator ``Gsin`` presented by ``cntr::function<T>``.
