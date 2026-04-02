.. _PManView:

Map
===

.. contents::
   :local:
   :depth: 2

.. _PManViewS01:

Overview
--------

.. list-table::
   :header-rows: 0

   * - class
     - ``cntr::herm_matrix_timestep_view<T>``

The purpose of this class is to create a map to a pre-existing object of type ``cntr::herm_matrix`` or ``cntr::herm_matrix_timestep``. One can then use this class without making a physical copy of the original data. The structure is inherited from the existing object. The usage of this class is mainly reserved for active developers as it is employed to enhance a performance of low-lying routines and ``MPI`` communications. It is characterized by the following parameters:

- ``T`` (template parameter): Precision, usually set to ``double``; we use the definition

  .. code-block:: cpp

     #define GREEN_TSTP_VIEW cntr::herm_matrix_timestep_view<double>

- ``tstp`` (integer): The timestep variable, ``tstp>=-1``.
- ``ntau`` (integer): number of discretization points on the imaginary time axis, ``ntau>=0``
- ``size1`` (integer): orbital dimension. Each element :math:`C(t,t')` is a square matrix of dimension ``size1`` :math:`\times` ``size1``.
- ``sig`` (``FERMION`` or ``BOSON``).

.. _PManViewS02:

Constructors
------------

.. list-table::
   :header-rows: 0

   * - ``herm_matrix_timestep_view<T>()``
     - Default constructor, does not allocate memory and sets ``nt=-2``.
   * - ``herm_matrix_timestep_view<T>(int tstp,int ntau,int size1, int size2, int sig)``
     - Creates a map with all entries to ``0``. It requires ``tstp>=-1``, ``ntau>0``, ``size1>0``, ``size2>0``, and ``sig=FERMION`` or ``sig=BOSON``.
   * - ``herm_matrix_timestep_view(int tstp, herm_matrix<T> &g)``
     - Creates a map to the predefined ``herm_matrix`` object ``g`` at a given time ``tstp``. Similar constructor exists for ``herm_matrix_timestep``.

.. _PManViewS03:

Accessing and manipulation
--------------------------

In general, all functions described for ``herm_matrix_timestep`` argument in :ref:`PMan02S01` can be replaced by a ``herm_matrix_timestep_view`` argument following the same syntax.
