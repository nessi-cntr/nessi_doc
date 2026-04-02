.. _PMan00:

Preliminaries
=============

.. contents::
   :local:
   :depth: 2

.. _PMan00S01:

Namespaces
----------

.. list-table::
   :header-rows: 0

   * - ``cntr``
     - All relevant classes and functions of ``libcntr`` are defined here
   * - ``fourier``
     - Fourier transformation routines. Used almost exclusively internally.

Needed for developers only:

.. list-table::
   :header-rows: 0

   * - ``integration``
     - Weights for quadrature and differentiation routines.
   * - ``linalg``
     - Linear algebra routines (out-dated).

.. _PMan00S02:

Scalar and matrix types
-----------------------

**Complex numbers:**

The numerical routines in ``NESSi`` are built on standard complex-valued algebra.

For the most commonly-used double precision numbers, the following shortcuts are defined:

.. code-block:: cpp

   #include <complex>
   #define cdouble std::complex<double>
   #define II cdouble(0.0,1.0)

**Matrices and linear algebra:**

Green's functions :math:`C(t,t')` can be matrix-valued. We use the `eigen3 <http://eigen.tuxfamily.org/index.php?title=Main_Page>`_ library for handling matrices.

The following shortcuts are defined:

.. code-block:: cpp

   #include <eigen3/Eigen/Dense>
   // integer precision variable-size vector
   #define ivector VectorXi
   // double precision variable-size vector
   #define dvector VectorXd
   // complex-valued double precision variable-size vector
   #define cdvector VectorXcd
   // integer precision variable-size matrix
   #define imatrix MatrixXi
   // double precision variable-size matrix
   #define dmatrix MatrixXd
   // complex-valued double precision variable-size matrix
   #define cdmatrix MatrixXcd

**Example:**

.. code-block:: cpp

   cdmatrix mm(2,2);  // a complex-valued 2x2 matrix
   cdmatrix mm2;
   // set mm to sigma_x
   mm.setZero();
   mm(0,1)=1;
   mm(1,0)=1;
   mm2=mm*mm;  // mm^2
   // etc ...

(For a complete documentation, see `eigen3 <http://eigen.tuxfamily.org/index.php?title=Main_Page>`_)

.. _PMan00S03:

Some global #defines
--------------------

Listed for reference only. Explanation in the sections below.

.. code-block:: cpp

   #define CNTR_PI 3.14159265358979323846
   #define CNTR_MAT_FOURIER 0
   #define CNTR_MAT_CG 1
   #define CNTR_MAT_FIXPOINT 2
   #define FERMION -1
   #define BOSON 1
   #define GREEN cntr::herm_matrix<double>
   #define GREEN_TSTP cntr::herm_matrix_timestep<double>
   #define CFUNC cntr::herm_function<double>
