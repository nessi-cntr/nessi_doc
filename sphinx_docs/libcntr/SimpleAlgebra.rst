.. _PMan04:

Simple Algebra
==============

.. contents::
   :local:
   :depth: 2

Simple algebra operations include scalar multiplication, incrementation of a Green's function with another Green's function, and multiplication of Green's function with one-time functions. These operations are performed usually on one timestep of a Green's function. They are implemented as member functions of the ``cntr::herm_matrix`` and ``cntr::herm_matrix_timestep`` classes, with the same syntax for both objects.

.. _PMan04S01:

Scalar multiplication and incrementation
-----------------------------------------

.. list-table::
   :header-rows: 0

   * - ``A.smul(int tstp, std::complex<T> x)``
     - :math:`A(t,t')` set to :math:`xA(t,t')` on timestep ``tstp`` of ``A``
   * - ``A.set_timestep_zero(int tstp)``
     - :math:`A(t,t')` set to :math:`0` on timestep ``tstp`` of ``A``
   * - ``A.incr_timestep(int tstp,GG &B, std::complex<T> x)``
     - :math:`A(t,t')` set to :math:`A(t,t')+xB(t,t')` on timestep ``tstp`` of ``A``

- ``A`` and ``B`` can be of type ``herm_matrix<T>`` and ``herm_matrix_timestep<T>``
- If ``A`` is of type ``herm_matrix_timestep``, then ``A.tstp()==tstp`` is required (same for ``B``); if ``A`` is of type ``herm_matrix``, then ``A.nt()>=tstp`` is required (same for ``B``)
- Size consistency ``B.ntau()==A.ntau()`` and ``B.size1()==A.size1()`` required.

**Example:**

The following operation considers a set of Green's functions :math:`G_k` indexed by momentum :math:`k`, and computes the Fourier transform at a given ``j=2``:

.. math::

   G_{j}(t,t') = \frac{1}{N_k}\sum_{k=0}^{N_k-1} e^{i2\pi k j /N_k} G_k(t,t').

.. code-block:: cpp

   int nt=10;
   int ntau=100;
   int size1=1;
   int nk=10;
   GREEN_TSTP tGav;
   std::vector<GREEN> Gk(nk); // a vector of nk=10 Green's functions
   for(int k=0;k<nk;k++) Gk[k]=GREEN(nt,ntau,size1,FERMION); // allocation
   // ...
   // ... do something with Gk
   // ...
   // do Fourier transform on timestep -1 (Matsubara)
   // and store result in timestep variable tGav:
   {
       double nkinv=1.0/nk;
       int j=2;
       tstp=-1;
       tGav=GREEN_TSTP(tstp,ntau,size1,FERMION); // properly resize tGav
       tGav.set_timestep_zero(tstp); // in principle not needed, as tGav is initialized with 0 already
       for(int k=0;k<nk;k++){
           cdouble weight=exp(II*2.0*M_PI*j*k*nkinv);
           tGav.incr_timestep(tstp,Gk[k],weight);
       }
       tGav.smul(tstp,nkinv);
   }
   // the syntax is analogous if tGav is replaced by a herm_matrix type Green's function

.. _PMan04S02:

Multiplication with one-time contour functions
-----------------------------------------------

.. list-table::
   :header-rows: 0

   * - ``A.left_multiply(int tstp,cntr::function<T> &f, T x)``
     - :math:`A(t,t')` set to :math:`x f(t) A(t,t')` at timestep ``tstp`` of ``A``
   * - ``A.right_multiply(int tstp,cntr::function<T> &f, T x)``
     - :math:`A(t,t')` set to :math:`x A(t,t') f(t')` at timestep ``tstp`` of ``A``
   * - ``A.left_multiply_hermconj(int tstp,cntr::function<T> &f, T x)``
     - :math:`A(t,t')` set to :math:`x f(t)^\dagger A(t,t')` at timestep ``tstp`` of ``A``
   * - ``A.right_multiply_hermconj(int tstp,cntr::function<T> &f, T x)``
     - :math:`A(t,t')` set to :math:`x A(t,t') f(t')^\dagger` at timestep ``tstp`` of ``A``

- ``A`` is of type ``cntr::herm_matrix<T>`` or ``cntr::herm_matrix_timestep<T>``
- ``x`` is a scalar of type T, can be omitted (default ``x=1``)
- The following size consistency is required:

  - If ``A`` is ``herm_matrix_timestep``, then ``tstp==A.tstp()`` is required; if ``A`` is ``herm_matrix``, then ``tstp<=A.nt()`` is required
  - ``tstp<=f.nt()`` and ``f.size1()==A.size1()`` is required.

- The last two routines multiply with the adjoint matrix :math:`f(t)^\dagger` of :math:`f(t)`

**Example:**

The following example computes a retarded interaction :math:`W(t,t') = U(t) \chi(t,t') U(t')` out of a given susceptibility :math:`\chi`, where :math:`U(t)` is a time-dependent bare interaction matrix element.

.. code-block:: cpp

   int nt=10;
   int ntau=100;
   int size1=1;
   GREEN W(nt,ntau,size1,BOSON);
   GREEN chi(nt,ntau,size1,BOSON);
   CFUNC U(nt,size1);
   // ...
   // ... do something to set U and chi
   // ...
   // compute W on all timesteps
   for(int tstp=-1;tstp<=nt;tstp++){
       W.set_timestep(tstp,chi);
       W.left_multiply(tstp,U);
       W.right_multiply(tstp,U);
   }
