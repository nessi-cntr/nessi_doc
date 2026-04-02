.. _PMan09:

Utilities
=========

.. contents::
   :local:
   :depth: 2

.. _PMan09S01:

Comparing Green's functions
----------------------------

To compare the data of two Green's functions on a given timestep, one can use:

.. list-table::
   :header-rows: 0

   * - ``T cntr::distance_norm2(int tstp,GType &A,GType &B)``
     - Returns a difference measure :math:`\Delta[A,B]_{\tt tstp}` (see below) between ``A`` and ``B`` on timestep ``tstp``.

- The type GType of ``A`` and ``B`` is ``cntr::herm_matrix<T>`` or ``cntr::herm_matrix_timestep<T>``, the return type is ``T`` (template parameter)
- Size requirements:

  - ``A`` and ``B`` must have equal ``size1`` and ``ntau``
  - For ``X=A,B``: If ``X`` is ``herm_matrix``, then ``X.nt()>=tstp`` is required; if ``X`` is ``herm_matrix_timestep``, then ``X.tstp()==tstp`` is required.

The difference is defined as the :math:`L_2`-norm difference :math:`|| M ||` for :math:`M=A(t,t')-B(t,t')` of the individual elements, summed over all time-arguments of the timestep:

.. math::
   :label: eq:distnorm

   \text{For tstp}=-1:\,\,\,
   \Delta[A,B]
   =&
   \sum_{i=0}^{\tt ntau}
   ||A^\mathrm{M}(i\Delta\tau)-B^M(i\Delta\tau)||
   \\
   \text{For tstp}>=0:\,\,\,
   \Delta[A,B]
   =&
   \sum_{i=0}^{\tt tstp}
   \big(
   ||A^\mathrm{R}({\tt tstp}\,\Delta t,i\Delta t)-B^\mathrm{R}({\tt tstp}\,\Delta t,i\Delta t)||
   +
   ||A^<(i\Delta t,{\tt tstp}\,\Delta t)-B^<(i\Delta t,{\tt tstp}\,\Delta t)||
   \big)
   \nonumber
   \\
   &+
   \sum_{i=0}^{\tt ntau}
   ||A^{\rceil}({\tt tstp}\,\Delta t,i\Delta\tau)-B^{\rceil}({\tt tstp}\,\Delta t,i\Delta\tau)||

:math:`|| M ||` is the standard matrix :math:`L_2` norm.

It is also possible to compare the individual components entering :eq:`eq:distnorm` separately. For ``tstp>=0``, one obtains the differences of the retarded/lesser/left-mixing components by ``cntr::distance_norm2_ret``/``cntr::distance_norm2_les``/``cntr::distance_norm2_tv``. The interface is identical to ``cntr::distance_norm2``.

**Example:**

A typical application is to check for convergence in self-consistent simulations:

.. code-block:: cpp

   // int nt= ...
   // int ntau= ...
   // int sig=...
   GREEN G(nt,ntau,size1,sig);
   // int tstp= ...
   // some iterative procedure to determine G on timestep tstp:
   {
       GREEN_TSTP tG(tstp,ntau,size1,sig);   //temporary variable
       double convergence_error;
       int iter_max=100;
       for(int iter=0;iter<=iter_max;iter++){
           tG.set_timestep(tstp,G);          // store values of G before iteration
           // ... some code to update G on timestep tstp ...
           convergence_error = cntr::distance_norm2(tstp,G,tG);
           if( convergence_error < some_sufficiently_small_number) break;
       }
       if(iter>iter_max){
           cout << "no convergence!" << endl;
       }
   }

.. _PMan09S02:

Timestep extrapolation
-----------------------

**Extrapolation:**

.. list-table::
   :header-rows: 0

   * - ``void cntr::extrapolate_timestep(int tstp,Gtype &A,int ExtrapolationOrder)``
     - Extrapolate from timesteps ``t=tstp,tstp-1,...,tstp-ExtrapolationOrder`` to timestep ``tstp+1``, using polynomial extrapolation. **output:** extrapolated object ``A``

- The type GType of ``A`` is ``cntr::herm_matrix<T>`` or ``cntr::function<T>``
- If ``A`` is ``herm_matrix``, all entries :math:`A(t,t')` on timestep ``tstp+1`` are written; if ``A`` is ``function``, the value ``A(tstp+1)`` is written.
- Size requirements: ``tstp>=ExtrapolationOrder`` and ``A.nt()>=tstp+1`` required
- The ``ExtrapolationOrder`` must be between ``1`` and ``MAX_SOLVE_ORDER`` (``=5``).

**Continuation from Imaginary time to Real time 0:**

.. list-table::
   :header-rows: 0

   * - ``void cntr::set_t0_from_mat(herm_matrix<T> &A)``
     - Sets the values at timestep ``tstp=0`` from ``tstp=-1`` assuming that :math:`A(t,t')` is continuous.

- Size requirements: ``A.nt()>=0`` required

Continuity of :math:`A(t,t')` implies that the point ``0`` on the lower real-time branch and on the imaginary time branch of the contour are identical. This implies timestep ``0`` can be set from timestep ``-1`` using the relation (:math:`\xi` is the ``FERMION`` or ``BOSON`` sign):

- :math:`G^{\rceil}(0,j\Delta \tau) = \xi G^M(({\tt ntau} - j)\Delta \tau)` for ``j=0,...,ntau``
- :math:`G^{<}(0,0) = G^{\rceil}(0,0)`
- :math:`G^{R}(0,0) = i G^M(0) - G^{<}(0,0)`

.. note::

   A two-time function need not be continuous if it explicitly depends on a parameter which changes discontinuously from the imaginary time to real-time contour. For example, the self-energy depends on the interaction ``U``, which changes discontinuously for an interaction quench at time 0.

.. _PMan09S03:

Differentiation
---------------

In order to perform differentiation on the real part of the contour, one can use:

.. list-table::
   :header-rows: 0

   * - ``void deriv1_timestep(int tstp, Gtype &dA, Gtype &A, Gtype &Acc, T beta, T h, int SolverOrder)``
     - Calculates the left derivative ``dA`` at time step ``tstp``, i.e. :math:`dA(t,t') = id/dt A(t,t')`, when ``A`` and its complex conjugate ``Acc`` are given.
   * - ``void deriv2_timestep(int tstp, Gtype &dA, Gtype &A, Gtype &Acc, T beta, T h, int SolverOrder)``
     - Calculates the right derivative ``dA`` at time step ``tstp``, i.e. :math:`dA(t,t') = -id/dt' A(t,t')`, when ``A`` and its complex conjugate ``Acc`` are given.

- The type ``Gtype`` of ``A``, ``Acc``, and ``dA`` is ``cntr::herm_matrix<T>``. Currently, the routines are implemented only for objects of ``size1=1``
- ``Acc`` is the conjugate function to ``A`` (similar as in :ref:`PMan06`). If ``A`` is hermitian, simply provide ``Acc=A``.
- For the computation of timestep ``tstp``, :math:`A(t,t')` is addressed at :math:`t,t' \leq \max(n,k)`, where :math:`k` is the solver order. Note that ``tstp=0..k`` can be computed only if ``A`` is given for :math:`t,t'\leq k`.
- The possible contribution proportional to :math:`\delta(t,t')`, arising from the step discontinuity of a Green's function at equal times, is subtracted in these routines.
