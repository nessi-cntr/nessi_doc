.. _P2:

Getting Started with NESSi
==========================

.. contents::
   :local:
   :depth: 2

Getting the Code
----------------

The **latest stable release** of the NESSi package is **NESSi-1.0.2** and can be
cloned from the GitHub repository:

https://github.com/nessi-cntr/nessi

.. code-block:: sh

   git clone https://github.com/nessi-cntr/nessi

Download the archive:

- `tar.gz <nessi-1.0.2.tar>`_
- `zip <nessi-1.0.2.zip>`_

.. note::

   This will be updated for v2.0.0.

Installation
------------

The NESSi repository (``nessi``) contains two folders:

- ``libcntr`` — the library and Python tools (see :ref:`PMan13`)
- ``examples`` — example programs (see :ref:`P4`)

Dependencies for ``libcntr``:

.. list-table::
   :header-rows: 0
   :align: left

   * - ``eigen3``
     - required
   * - ``HDF5``
     - optional, recommended
   * - ``FFTW3``
     - optional, needed for steady state code

``libcntr`` is compiled and installed using **CMake**. Several CMake variables
must be set (``-Dvar=value``). While CMake can be called directly, it is more
convenient to create a configuration script:

Example ``configure.sh``:

.. code-block:: sh

   CC=[C compiler] CXX=[C++ compiler]
   cmake \
     -DCMAKE_INSTALL_PREFIX=[install directory] \
     -DCMAKE_BUILD_TYPE=[Debug|Release] \
     -Domp=[ON|OFF] \
     -Dhdf5=[ON|OFF] \
     -Dmpi=[ON|OFF] \
     -DBUILD_DOC=[ON|OFF] \
     -Dness=[ON|OFF] \
     -DCMAKE_INCLUDE_PATH=[include directory] \
     -DCMAKE_LIBRARY_PATH=[library directory] \
     -DCMAKE_CXX_FLAGS="[compiling flags]" \
     ..

In the first line the C and C++ compilers are set (tested with GNU ``gcc/g++``
and Intel ``icc/icpc``). The installation directory (example: ``/home/opt``) is
provided via ``CMAKE_INSTALL_PREFIX``.

``CMAKE_BUILD_TYPE=Debug`` enables assertions and checks; ``Release`` is
recommended for production runs.

Optional features:

- ``omp=ON`` — enable OpenMP parallelization  
- ``mpi=ON`` — enable MPI routines  
- ``ness=ON`` — enable steady-state routines (requires OpenMP + FFTW with
  ``--enable-threads``)  
- ``hdf5=ON`` — use the ``hdf5`` library (recommended version ≥ 1.12)

Paths to required libraries are provided via ``CMAKE_INCLUDE_PATH`` and
``CMAKE_LIBRARY_PATH``.

You must include the following compilation flag:

.. code-block:: sh

   -std=c++11

**Compiling and installing**

Create a build directory (e.g. ``cbuild``), then run:

.. code-block:: sh

   sh ../configure.sh

Compile:

.. code-block:: sh

   make

Install:

.. code-block:: sh

   make install

Supported Platforms
-------------------

We have tested installation on:

- macOS (MacPorts)
- macOS (Homebrew)
- Linux (Ubuntu, Debian, CentOS, Arch Linux)

Detailed instructions follow.

macOS with MacPorts
~~~~~~~~~~~~~~~~~~~

Install dependencies:

.. code-block:: sh

   sudo port install gcc9 eigen3-devel hdf5
   sudo port install openmpi-devel-gcc9

Install Doxygen:

.. code-block:: sh

   sudo port install doxygen graphviz

Python packages (``numpy``, ``scipy``, ``matplotlib``, ``h5py``) can be
installed via ``pip`` or MacPorts.

Example ``configure.sh``:

.. code-block:: sh

   CC=mpicc CXX=mpicxx
   cmake \
     -DCMAKE_INSTALL_PREFIX=$HOME/opt \
     -DCMAKE_INSTALL_NAME_DIR=$HOME/opt/lib \
     -DCMAKE_BUILD_TYPE=Release \
     -Domp=ON \
     -Dhdf5=ON \
     -Dmpi=ON \
     -DBUILD_DOC=ON \
     -Dness=ON \
     -DCMAKE_INCLUDE_PATH=/opt/local/include \
     -DCMAKE_LIBRARY_PATH=/opt/local/lib \
     -DCMAKE_CXX_FLAGS="-std=c++11 -O3" \
     ..

macOS with Homebrew
~~~~~~~~~~~~~~~~~~~

Unlink legacy Python:

.. code-block:: sh

   brew unlink python@2

Install dependencies:

.. code-block:: sh

   brew install eigen hdf5 open-mpi
   brew install doxygen graphviz

Example ``configure.sh``:

.. code-block:: sh

   CC=mpicc CXX=mpicxx
   cmake \
     -DCMAKE_INSTALL_PREFIX=$HOME/opt \
     -DCMAKE_INSTALL_NAME_DIR=$HOME/opt/lib \
     -DCMAKE_BUILD_TYPE=Release \
     -Domp=OFF \
     -Dhdf5=ON \
     -Dmpi=ON \
     -DBUILD_DOC=ON \
     -Dness=ON \
     -DCMAKE_INCLUDE_PATH=/usr/local/include \
     -DCMAKE_LIBRARY_PATH=/usr/local/lib \
     -DCMAKE_CXX_FLAGS="-std=c++11 -O3" \
     ..

Note: ``omp=OFF`` because Apple ``clang`` does not support OpenMP directly.
Workarounds exist (see `here <https://stackoverflow.com/questions/46414660/macos-cmake-and-openmp>`_).

Linux
~~~~~

On Debian/Ubuntu:

.. code-block:: sh

   apt-get install -y --allow-unauthenticated libhdf5-serial-dev libopenmpi-dev libeigen3-dev doxygen graphviz cmake

On Arch Linux:

.. code-block:: sh

   pacman -Sy --noconfirm hdf5 gcc openmpi eigen doxygen graphviz make cmake

Example ``configure.sh`` (paths usually not needed):

.. code-block:: sh

   CC=mpicc CXX=mpicxx
   cmake \
     -DCMAKE_INSTALL_PREFIX=$HOME/opt \
     -DCMAKE_BUILD_TYPE=Release \
     -Domp=ON \
     -Dhdf5=ON \
     -Dmpi=ON \
     -DBUILD_DOC=ON \
     -Dness=ON \
     -DCMAKE_INCLUDE_PATH="" \
     -DCMAKE_LIBRARY_PATH="" \
     -DCMAKE_CXX_FLAGS="-std=c++11 -O3" \
     ..

If ``$HOME/opt/lib`` is in ``LD_LIBRARY_PATH``, no ``CMAKE_INSTALL_NAME_DIR`` is needed.

.. _P2S2:

Creating a Custom Program
-------------------------

Minimal working example:

.. code-block:: cpp

   #include <sys/stat.h>
   #include <iostream>
   #include "cntr/cntr.hpp"

   int main(int argc, char *argv[]) {
       GREEN G;
       int Nt   = 100; // # of timesteps
       int Ntau = 100; // # of imaginary-time steps
       int size = 3;   // GF size

       G = GREEN(Nt, Ntau, size, FERMION);

       std::cout << "Number of timesteps: " << G.nt() << std::endl;
       std::cout << "Number of Matsubara points: " << G.ntau() << std::endl;
       std::cout << "Size of Greens function: " << G.size1() << std::endl;

       return 0;
   }

Test Suite and Documentation
----------------------------

To run the test suite (Catch framework):

.. code-block:: sh

   make test

After completion:

``All tests passed``

To test MPI routines:

.. code-block:: sh

   make test_mpi

Documentation is generated automatically when ``BUILD_DOC=ON`` is set in the
configure script. After building, the HTML docs are at:

``doc/html/index.html``

Links
-----

- `eigen3 <http://eigen.tuxfamily.org/index.php?title=Main_Page>`_
- `hdf5 <https://www.hdfgroup.org/solutions/hdf5/>`_
- `doxygen <http://www.doxygen.nl>`_
- `nessi <https://github.com/nessi-cntr/nessi>`_
- `MacPorts <https://www.macports.org>`_
- `homebrew <https://brew.sh>`_
- `here <https://stackoverflow.com/questions/46414660/macos-cmake-and-openmp>`_
- `FFTW3 <https://www.fftw.org>`_