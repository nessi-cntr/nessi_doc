.. _PMan06:

Convolution
===========

.. contents::
   :local:
   :depth: 2

.. _PMan06S01:

Full convolution
----------------

The convolution of two Green's functions ``A`` and ``B`` is defined by the integral:

.. math::
   :label: convolution001

   C(t,t^\prime) \equiv [A\ast B](t,t^\prime) = \int_{\mathcal{C}} d \bar{t}\, A(t,\bar{t}) B(\bar{t},t^\prime),

where the product of :math:`A(t,\bar{t})` and :math:`B(\bar{t},t^\prime)` is understood as a matrix product in orbital indices. The convolution is one of the most basic operations when dealing with nonequilibrium Green's functions. The evaluation of the integral uses Langreth rules, to derive separate integrals for each Keldysh component of ``C`` in terms of the Keldysh components of ``A`` and ``B``.

We provide a number of routines to calculate the generalized convolution, where a one-time function ``f(t)`` is sandwiched between ``A`` and ``B``:

.. math::
   :label: convolution002

   C(t,t^\prime) = \int_{\mathcal{C}} d \bar{t}\, A(t,\bar{t}) f(\bar{t}) B(\bar{t},t^\prime).

This integral is computed by the routine:

.. list-table::
   :header-rows: 0

   * - ``cntr::convolution(GG &C, GG &A, GG &Acc,cntr::function<T> &f, GG &B, GG &Bcc, T beta, T dt, int SolveOrder)``
     - Obtain :math:`C(t,t^\prime) = \int_{\mathcal{C}} d \bar{t}\, A(t,\bar{t}) f(\bar{t}) B(\bar{t},t^\prime)` in the full two-time plane with given two-time objects ``A`` and ``B``, and a function ``f``. **output:** object ``C``
   * - ``cntr::convolution_timestep(int tstp, GG &C, GG &A, GG &Acc,cntr::function<T> &f, GG &B, GG &Bcc, T beta, T dt, int SolveOrder)``
     - Obtain :math:`C(t,t^\prime) = \int_{\mathcal{C}} d \bar{t}\, A(t,\bar{t}) f(\bar{t}) B(\bar{t},t^\prime)` on timestep ``tstp`` of ``C`` with given two-time objects ``A`` and ``B``, and a function ``f``. **output:** object ``C``
   * - ``cntr::convolution_timestep_omp(int omp_num_threads,int tstp, GG &C, GG &A, GG &Acc,cntr::function<T> &f, GG &B, GG &Bcc, T beta, T dt, int SolveOrder)``
     - Same as ``cntr::convolution_timestep``, using shared memory parallelization. **output:** object ``C``

- ``C``, ``A``, ``Acc``, ``B``, ``Bcc`` are of type ``cntr::herm_matrix``
- ``f`` is of type ``cntr::function``. If the argument is omitted, :math:`f(t)=1` is assumed.
- ``dt`` is the time-discretisation step :math:`\Delta t` on the real-time branch
- ``beta`` is the length of the imaginary-time contour (inverse temperature)
- For ``X=A,B``, ``Xcc`` is the hermitian conjugate of ``X``. If the argument is omitted, hermitian symmetry ``Xcc=X`` is assumed.
- The arguments must satisfy the following size consistency:

  - ``X.nt()>=tstp`` for all ``C``, ``A``, ``Acc``, ``B``, ``Bcc``, ``f``
  - ``X.ntau()`` equal for all ``C``, ``A``, ``Acc``, ``B``, ``Bcc``
  - ``X.size1()`` equal for all ``C``, ``A``, ``Acc``, ``B``, ``Bcc``, ``f``

- ``0 <= SolveOrder <= MAX_SOLVE_ORDER`` is the integration order (currently, ``MAX_SOLVE_ORDER=5``). Use ``SolveOrder=5`` (``=MAX_SOLVE_ORDER``) if there is no good reason against it.
- To be able to perform the integration up to that order, the additional size requirements hold:

  - ``X.ntau()>SolveOrder`` required for all ``C``, ``A``, ``Acc``, ``B``, ``Bcc``
  - If ``tstp>-1``: ``X.nt()>=SolveOrder`` required for all ``C``, ``A``, ``Acc``, ``B``, ``Bcc``, ``f``

- The second variant ``convolution_timestep_omp`` uses shared memory ``openMP`` parallelization to distribute the evaluation of :math:`C(t,t')` for different time arguments on the timestep ``tstp`` over ``omp_num_threads`` tasks.

**Note on hermitian conjugate:**

For the hermitian conjugate one finds:

.. math::
   :label: convcc

   [C^\ddagger](t,t^\prime) =
   \int_{\mathcal{C}} d \bar{t}\, B^\ddagger(t,\bar{t}) f(\bar{t})^\dagger A^\ddagger(\bar{t},t^\prime).

In general, therefore, :math:`C` is not hermitian. ``cntr::convolution_timestep`` computes ``C`` only on the hermitian domain. To restore the full function ``C``, one needs :math:`C^\ddagger` on the hermitian domain. Following :eq:`convcc`, :math:`C^\ddagger` is computed by a second call to ``cntr::convolution_timestep``, storing the result to a separate ``herm_matrix`` variable ``Ccc``:

.. code-block:: cpp

   cntr::convolution_timestep(tstp,C,A,Acc,f,B,Bcc,beta,dt,SolveOrder);
   cntr::convolution_timestep(tstp,Ccc,Acc,A,fcc,Bcc,B,beta,dt,SolveOrder); // fcc is the function fcc(t)=f(t)^+

Remark: The library provides a number of un-supported routines which calculate only certain Keldysh components for the convolution, as well as variants with different argument lists.

**Note on causal time-dependence and time-stepping:**

The exact convolution is always causal. The numerical implementation is exactly causal for timesteps ``tstp>SolveOrder``, while the numerical error on timesteps ``0<=tstp<=SolveOrder`` can depend on the input (``A``, ``Acc``, ``B``, ``Bcc``, ``f``) on all timesteps ``-1<=tstp<=SolveOrder``. More precisely, the result for ``C`` on timestep ``tstp`` depends on the following input data:

- If ``tstp>SolveOrder`` or ``tstp==-1``: Results depend only on timesteps ``t<=tstp`` of the input ``A``, ``Acc``, ``B``, ``Bcc``, and on ``f(t)`` for ``t <=tstp``. This is called a causal time-dependence.
- If ``0<tstp<=SolveOrder``: Results depend on timesteps ``-1<=t<=SolveOrder`` of the input ``A``, ``Acc``, ``B``, ``Bcc``, and ``f``.

**Example:**

The example computes the convolution

.. math::

   C = G_0 \ast \Sigma \ast G_0

which is part of the Dyson series. The convolution is calculated in two steps:

.. math::

   A = G_0 \ast \Sigma,
   \,\,\,\,
   C = A \ast G_0.

Because ``A`` is not hermitian, three calls to ``cntr::convolution_timestep`` are needed, which calculate both :math:`A` and :math:`A^\ddagger` on the hermitian domain:

.. code-block:: cpp

   int nt=10;
   int ntau=20;
   int size1=2;
   int SolveOrder=5;
   double dt=0.01; // problem dependent
   double beta=10.0;
   GREEN G0(nt,ntau,size1,FERMION);
   GREEN Sigma(nt,ntau,size1,FERMION);
   GREEN C(nt,ntau,size1,FERMION);
   // the temporary A and its herm. conjugate Acc
   GREEN A(nt,ntau,size1,FERMION);
   GREEN Acc(nt,ntau,size1,FERMION);
   // do anything to fill G0 and Sigma
   // calculate A and Acc
   for(int tstp=-1;tstp<=nt;tstp++){
       // because G0 and Sigma are hermitian, some arguments are omitted
       cntr::convolution_timestep(tstp,A,G0,Sigma,beta,dt,SolveOrder);
       cntr::convolution_timestep(tstp,Acc,Sigma,G0,beta,dt,SolveOrder);
   }
   // now compute C
   for(int tstp=-1;tstp<=nt;tstp++){
       cntr::convolution_timestep(tstp,C,A,Acc,G0,G0,beta,dt,SolveOrder);
   }

.. note::

   As mentioned in the comment on causal time-dependency, the input ``A`` and ``Acc`` for ``C`` must be determined on all timesteps ``0<=tstp<=SolveOrder`` before ``C`` is computed on timesteps ``0<=tstp<=SolveOrder``.

.. code-block:: cpp

   // ...
   int SolveOrder=5;
   // INCORRECT:
   for(int tstp=-1;tstp<=nt;tstp++){
       cntr::convolution_timestep(tstp,A,G0,Sigma,beta,dt,SolveOrder);
       cntr::convolution_timestep(tstp,Acc,Sigma,G0,beta,dt,SolveOrder);
       // does not work for tstp=0...4, because for tstp=0...5 routine needs input A and Acc on all times 0...5
       cntr::convolution_timestep(tstp,C,A,Acc,G0,G0,beta,dt,SolveOrder);
   }

   // Alternative (CORRECT):
   // convolution on matsubara:
   tstp=-1;
   cntr::convolution_timestep(tstp,A,G0,Sigma,beta,dt,SolveOrder);
   cntr::convolution_timestep(tstp,Acc,Sigma,G0,beta,dt,SolveOrder);
   cntr::convolution_timestep(tstp,C,A,Acc,G0,G0,beta,dt,SolveOrder);
   // convolution on warm-up steps:
   for(int tstp=0;tstp<=SolveOrder;tstp++){
       cntr::convolution_timestep(tstp,A,G0,Sigma,beta,dt,SolveOrder);
       cntr::convolution_timestep(tstp,Acc,Sigma,G0,beta,dt,SolveOrder);
   }
   for(int tstp=0;tstp<=SolveOrder;tstp++){
       cntr::convolution_timestep(tstp,C,A,Acc,G0,G0,beta,dt,SolveOrder);
   }
   // convolution on all other steps:
   for(int tstp=SolveOrder+1;tstp<=nt;tstp++){
       cntr::convolution_timestep(tstp,A,G0,Sigma,beta,dt,SolveOrder);
       cntr::convolution_timestep(tstp,Acc,Sigma,G0,beta,dt,SolveOrder);
       cntr::convolution_timestep(tstp,C,A,Acc,G0,G0,beta,dt,SolveOrder);
   }

**OpenMP parallelization:**

The variant ``convolution_timestep_omp`` uses shared memory ``openMP`` parallelization to distribute the evaluation of :math:`C(t,t')` for different time arguments of the timestep ``tstp`` over ``omp_num_threads`` tasks. The implementation is free of any race condition if ``omp_num_threads`` is set to the number of threads in the current team.

.. _PMan06S02:

Convolution density matrix
---------------------------

In typical applications one is often only interested in the density matrix of the convolution integral :eq:`convolution001` or :eq:`convolution002` (see section :ref:`PMan01S04` for the definition of the density matrix). The following routine may be of use:

.. list-table::
   :header-rows: 0

   * - ``cntr::convolution_density_matrix(int tstp, Type &M, GG &C, GG &A, GG &Acc,cntr::function<T> &f, GG &B, GG &Bcc,T beta, T dt,int SolveOrder)``
     - Obtain density matrix for Green function :math:`C(t,t^\prime) = \int_{\mathcal{C}} d \bar{t}\, A(t,\bar{t}) f(\bar{t}) B(\bar{t},t^\prime)` with given two-time objects ``A,B`` at time ``tstp``, and write result to ``M``

- Same comments on input ``A``, ``Acc``, ``B``, ``Bcc`` as for ``cntr::convolution_timestep`` above
- If ``A.size1()==1``, ``M`` can be a scalar ``std::complex<T>`` or a complex eigenmatrix ``cdmatrix``; if ``A.size1()>1``, ``M`` must be complex eigenmatrix ``cdmatrix``
- If ``M`` is a complex eigenmatrix, it is automatically resized to ``size1`` :math:`\times` ``size1``
- ``f`` is of type ``cntr::function``. If the argument is omitted, :math:`f(t)=1` is assumed.

**Example:**

Calculate the correlation energy for a two-particle interaction from the self-energy:

.. math::

   E_{int} = \frac{1}{2} \text{tr} \Big[  \text{DensityMatrix}\big(\Sigma \ast G \big)(t)\Big]

.. code-block:: cpp

   int nt=10;
   int ntau=20;
   int size1=2;
   int SolveOrder=5;
   double dt=0.01; // time-discretization
   double beta=10.0;
   GREEN G(nt,ntau,size1,FERMION);
   GREEN Sigma(nt,ntau,size1,FERMION);
   // calculate G and Sigma on all timesteps
   for(int tstp=-1;tstp<=nt;tstp++){
       cdmatrix mm;
       cntr::convolution_density_matrix(tstp,mm,Sigma,G,beta,dt,SolveOrder);
       cout << "Eint at t= " << tstp << " : " << 0.5*mm.trace().real() << endl;
   }

.. note::

   The user can also use the function ``cntr::correlation_energy``, which implements the above integral directly.
