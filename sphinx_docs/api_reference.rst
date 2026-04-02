.. _api_reference:

C++ API Reference
=================

.. contents::
   :local:
   :depth: 1

This reference is generated automatically from the source code via Doxygen and Breathe.

libcntr — Green's function classes
------------------------------------

.. doxygenclass:: cntr::herm_matrix
   :members:

.. doxygenclass:: cntr::herm_matrix_timestep
   :members:

.. doxygenclass:: cntr::herm_matrix_timestep_view
   :members:

.. doxygenclass:: cntr::herm_matrix_moving
   :members:

.. doxygenclass:: cntr::herm_matrix_timestep_moving
   :members:

.. doxygenclass:: cntr::herm_matrix_timestep_moving_view
   :members:

.. doxygenclass:: cntr::herm_pseudo
   :members:

.. doxygenclass:: cntr::function
   :members:

.. doxygenclass:: cntr::function_moving
   :members:

libcntr — Utility classes
--------------------------

.. doxygenclass:: cntr::cyclic_timestep
   :members:

libcntr — MPI / distributed classes
--------------------------------------

.. doxygenclass:: cntr::distributed_array
   :members:

.. doxygenclass:: cntr::distributed_timestep_array
   :members:

.. doxygenclass:: cntr::distributed_timestep_array_moving
   :members:

Steady-state routines (namespace ``ness2``)
--------------------------------------------

.. doxygenclass:: ness2::herm_matrix_ness
   :members:

.. doxygenclass:: ness2::grid_info
   :members:

.. doxygenclass:: ness2::fft_array
   :members:

.. doxygenclass:: ness2::FFT_OMP_Manager
   :members:

.. doxygenclass:: ness2::bethedos
   :members:

.. doxygenclass:: ness2::gauss
   :members:

.. doxygenclass:: ness2::ohmic_sym
   :members:

Steady-state routines — legacy (namespace ``ness``)
-----------------------------------------------------

.. doxygenclass:: ness::GF
   :members:

.. doxygenclass:: ness::gauss
   :members:

.. doxygenclass:: ness::ohmic_sym
   :members:
