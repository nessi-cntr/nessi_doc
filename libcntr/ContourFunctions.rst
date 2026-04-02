.. _PMan03:

Contour Functions
=================

.. contents::
   :local:
   :depth: 2

.. _PMan03S01:

Overview
--------

.. list-table::
   :header-rows: 0

   * - class
     - ``cntr::function<T>``

A contour function is a matrix or scalar-valued function :math:`f(t)` which depends on physical time only:

- On the imaginary time branch, it takes one given constant value :math:`f_{-1}`
- On the real-time branches, its value is the same for the same time :math:`t` on the upper and lower contour branch.

A contour function is typically used to store time-dependent parameters of a system, or a time-dependent Hamiltonian.

The class ``cntr::function<T>`` is the container to store these data. It is characterized by the following parameters:

- ``T`` (template parameter): Precision, usually set to ``double``; we use the definition

  .. code-block:: cpp

     #define CFUNC cntr::herm_function<double>

- ``nt`` (integer): number of discretization points on the real time axis. For ``nt=-1``, only the (constant) value on the imaginary time axis is stored.
- ``size1`` (integer): orbital dimension. Each element :math:`f(t)` is a square matrix of dimension ``size1`` :math:`\times` ``size1``.

.. _PMan03S02:

Constructors
------------

.. list-table::
   :header-rows: 0

   * - ``function<T>()``
     - Default constructor, does not allocate memory and sets ``nt=-2``.
   * - ``function<T>(int nt,int size1)``
     - Allocate memory, set all entries to ``0``. It requires ``tstp>=-1``, ``size1>0``

.. _PMan03S03:

Accessing individual matrix elements
--------------------------------------

The following routines allow to read/write the elements of a contour function :math:`f(t)` stored as ``cntr::function`` at individual time arguments from/to another variable ``M``. The latter can be either a scalar of type ``std::complex<T>``, or a complex square matrix defined in Eigen (see :ref:`PMan00S02`).

The following member functions of ``cntr::function<T>`` set components of a contour function :math:`f(t)` from ``M``:

.. list-table::
   :header-rows: 0

   * - ``f.set_value(j,M)``
     - For ``j=-1``: :math:`f_{-1}` is set to ``M``; For ``0<=j<=f.nt()``: :math:`f(j\Delta t)` is set to ``M``

- If ``f.size1()>1``, ``M`` must be a square matrix (``cdmatrix`` for ``double`` precision)
- If ``f.size1()==1``, ``M`` can be also a scalar (``cdouble``) or a square matrix
- If ``M`` is a matrix, it must be a square matrix of dimension ``f.size1()``

The following member functions of ``cntr::function<T>`` read components of a contour function :math:`f` to ``M``:

.. list-table::
   :header-rows: 0

   * - ``f.get_value(j,M)``
     - For ``j=-1``: ``M`` is set to :math:`f_{-1}`; For ``0<=j<=f.nt()``: ``M`` is set to :math:`f(j\Delta t)`

- If ``M`` is a matrix, it is resized to a square matrix of dimension ``f.size1()``
- If ``M`` is a scalar, only the (0,0) entry of :math:`f(t)` is read.
