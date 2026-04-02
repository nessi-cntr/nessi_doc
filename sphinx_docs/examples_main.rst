.. _P4:

.. _PEx:

Example Programs
================

.. contents::
   :local:
   :depth: 2

.. _S4S1:

.. _PEx00:

Installation
------------

After downloading or cloning the NESSi repository `nessi <https://github.com/nessi-cntr/nessi>`_, the example programs are compiled in a similar fashion as the ``libcntr`` library. We have prepared a CMake build environment and it is most convenient to write a configuration shell script ``configure.sh`` of this type:

.. code-block:: sh

   CC=[C compiler] CXX=[C++ compiler]
   cmake \
        -DCMAKE_BUILD_TYPE=[Debug|Release] \
        -Domp=[ON|OFF] \
        -Dhdf5=[ON|OFF] \
        -Dmpi=[ON|OFF] \
        -DCMAKE_INCLUDE_PATH=[include directory] \
        -DCMAKE_LIBRARY_PATH=[library directory] \
        -DCMAKE_CXX_FLAGS="[compiling flags]" \
        ..

The dependencies for the example programs are the same as for ``libcntr``. The `eigen3 <http://eigen.tuxfamily.org/index.php?title=Main_Page>`_ library and, optionally (turned on by ``hdf5=ON``), the `HDF5 <https://www.hdfgroup.org/solutions/hdf5/>`_ library are required. Make sure the corresponding libraries can be found in the library path ``CMAKE_LIBRARY_PATH``, while the corresponding headers should be placed in ``CMAKE_INCLUDE_PATH``. Furthermore, the examples depend on the ``libcntr`` library. Therefore, the library path for ``libcntr`` should be provided by ``CMAKE_LIBRARY_PATH``, while the header ``cntr/cntr.hpp`` should be found in the include path.

Create a build directory (for instance, ``cbuild/``), navigate there and run the configure script:

.. code-block:: sh

   sh ../configure.sh

After successful configuration, compile via

.. code-block:: sh

   make

The executables are then found in the ``exe/`` directory.

The ``utils/`` directory contains useful python driver scripts which simplify the execution of the example programs. In order to run the python script, we need to make sure to set the python path to ``nessi/libcntr/python`` and/or ``nessi/libcntr/python3``. The scripts should be run from ``nessi/examples``. In the following table, we summarize the python scripts and the corresponding exe files and provide brief explanations of what is done in the python scripts.

.. list-table:: Summary of python scripts to run example codes in ``nessi``
   :header-rows: 1

   * - Script
     - Description
   * - ``test_equilibrium.py``
     - Runs the execute file ``test_equilibrium.x`` to show the scaling of accuracy of the Matsubara Dyson solvers with the specified order as an input. A figure for the scaling against :math:`N_{\tau}` is created.
   * - ``test_nonequilibrium.py``
     - Runs ``test_nonequilibrium.x`` to show the scaling of accuracy of the integro-differential (Dyson) and integral (VIE2) formulation with the specified order as an input. A figure for the scaling against :math:`N_{t}` is created.
   * - ``demo_hubbard_chain.py``
     - Runs ``hubbard_chain_**.x`` to simulate quench dynamics of the Hubbard chain. Here, ``**`` (= ``tt 2b, gw, tpp``) indicates different many body approximations, which can be specified in the python script. Figure for the time evolution of density and energies are created.
   * - ``demo_Holstein_impurity.py``
     - Runs ``Holstein_impurity_singlebath_**.x`` to simulate dynamics against modulation of system parameters in the Holstein-type impurity with a single bath site. Here, ``**`` (= ``Migdal, uMig``) indicates different approximate impurity solvers, which can be specified in the python script. The spectra of electrons and phonons and the time evolution of phonon displacement and energies are plotted.
   * - ``demo_Holstein.py``
     - Runs ``Holstein_bethe_**.x`` to simulate dynamics of the Holstein model against modulation of system parameters within DMFT. The rest is the same as ``demo_Holstein_impurity.py``.
   * - ``demo_Holstein_sc.py``
     - Runs ``Holstein_bethe_Nambu_**.x``, which is a generalized version of ``Holstein_bethe_**.x`` to treat s-wave superconductor (SC). In addition to figures for the spectra and evolution of densities and energies, evolution of the SC order parameter is plotted.
   * - ``demo_gw.py``
     - Runs ``gw.x`` to simulate the 1dim chain of the extended Hubbard model within the GW approximation using MPI parallelization. The script creates figures for the electric field and the change in the kinetic energy.
   * - ``demo_integration.py``
     - Runs ``integration.x`` to demonstrate the accuracy of the Gregory integration implemented in ``nessi``.


The truncated example can be found in :ref:`Trunc` and instructions for the steady-state examples in :ref:`NessEx`.

Example Pages
-------------

.. toctree::
   :maxdepth: 1

   examples/accuracy
   examples/hubbard_chain
   examples/holstein_impurity
   examples/holstein_dmft
   examples/hubbard_trans_inv
   examples/truncation
   examples/ness
