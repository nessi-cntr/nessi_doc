.. _PMan11:

HDF5 Usage
==========

.. contents::
   :local:
   :depth: 2

.. _PMan11S1:

Overview
--------

``libcntr`` uses the `HDF5 <https://www.hdfgroup.org/solutions/hdf5/>`_ format to store the basic data types for contour functions to disk. HDF5 is an open source library and file format for numerical data which is widely used in the field of scientific computing. The format has two building blocks:

- *data sets*: general multi-dimensional arrays of a single type
- *groups*: containers which can hold data sets and other groups.

Hence, by nesting groups, it is possible to store arbitrarily complicated structured data, and to create a file-system-like hierarchy where groups can be indexed using standard `POSIX <https://pubs.opengroup.org/onlinepubs/9699919799/>`_ format, e.g. ``/path/to/data``.

The ``libcntr`` library comes with helper functions to store the basic contour response function data types in HDF5 with a predefined structure of groups and data sets, defined in the header ``cntr/hdf5/hdf5_interface.hpp``. For example a ``herm_matrix`` response function is stored as a group with a data set for each contour component ``mat`` (:math:`g^M(\tau)`), ``ret`` (:math:`g^R(t, t')`), ``les`` (:math:`g^<(t, t')`), and ``tv`` (:math:`g^\rceil(t, \tau)`), respectively, see Section :ref:`PMan01`. The retarded and lesser components are stored in upper and lower triangular contiguous time order respectively. In the ``libcntr`` HDF5 format each component is stored as a rank 3 array where the first index is time, imaginary time, or triangular contiguous two-time, and the remaining two indices are orbital indices.

.. _PMan11S2:

Reading/writing to hdf5 files
------------------------------

To store a contour Green's function ``G`` of type ``cntr::herm_matrix`` or read it from a file, one can use the member functions ``cntr::herm_matrix::write_to_hdf5`` and ``cntr::herm_matrix::read_from_hdf5``. This stores/reads the attributes ``nt``, ``ntau``, ``sig``, ``size1``, ``size2``, ``element_size`` and the Green's function component's data sorted in groups ``ret``, ``les``, ``tv``, ``mat``.

**Implementation:**

The reading and writing functions are overloaded in a hierarchical fashion, so that hdf5 data can be accessed using filename, groupname and group id.

- The top-level writing functions, which take the file and group name as arguments, by default internally use the standard hdf5 C API to create a new file via ``H5Fcreate`` using the truncated write mode ``H5F_ACC_TRUNC``. If the file does not exist, a new file is created in truncated writing mode and if the file exists it will be overwritten. The subordinate function then creates/opens a hdf5 group with given handle and name using hdf5 API function ``H5Gcreate``. Then the final subordinate functions stores the attributes and data to the group handle, using hdf5 API functions ``H5Dcreate`` and ``H5Tcreate``, ``H5Tinsert`` to create a complex compound hdf5 type. Files and groups are closed after writing, using hdf5 API (``H5Fclose``, ``H5Gclose``, ``H5Dclose``, ...).

- The top-level reading function, which take the file and group name as arguments, internally uses the standard hdf5 C API to open the file via ``H5Fopen`` in read only mode (``H5F_ACC_RDONLY``). The subordinate function then opens a hdf5 group with given handle and name using hdf5 API function ``H5Gopen``. Then the final subordinate functions reads the data from the group handle, using hdf5 API functions ``H5Dopen`` and ``H5Dread`` to read the data. Files and groups are closed after reading, using hdf5 API (``H5Fclose``, ``H5Gclose``, ``H5Dclose``, ...).

.. note::

   The functions also work for ``cntr::herm_matrix_timestep`` objects and store the attribute ``tstp`` instead of ``nt``.

**hdf5 writing functions:**

.. list-table::
   :header-rows: 0

   * - ``G.write_to_hdf5(hid_t group_id)``
     - Stores the ``cntr::herm_matrix`` (attributes and data) to a given hdf5 group.
   * - ``G.write_to_hdf5(hid_t group_id, const char *groupname)``
     - Stores the ``cntr::herm_matrix`` (attributes and data) to a given hdf5 group with given groupname.
   * - ``G.write_to_hdf5(const char *filename, const char *groupname)``
     - Write data and attributes from the ``cntr::herm_matrix`` to a hdf5 file under the given group name.

**hdf5 reading functions:**

.. list-table::
   :header-rows: 0

   * - ``G.read_from_hdf5(hid_t group_id)``
     - Reads the ``cntr::herm_matrix`` (attributes and data) from a given hdf5 group handle.
   * - ``G.read_from_hdf5(hid_t group_id, const char *groupname)``
     - Reads the ``cntr::herm_matrix`` (attributes and data) from a given hdf5 group handle with given group name.
   * - ``G.read_from_hdf5(const char *filename, const char *groupname)``
     - Read all data and attributes from a hdf5 file from a given group into the ``cntr::herm_matrix``.

**hdf5 writing: timeslices:**

.. list-table::
   :header-rows: 0

   * - ``G.write_to_hdf5_slices(hid_t group_id, int dt)``
     - Write data and attributes from every ``dt``-th time step of ``cntr::herm_matrix`` object ``G`` to a given hdf5 group handle.
   * - ``G.write_to_hdf5_slices(hid_t group_id, const char *groupname, int dt)``
     - Write data and attributes from every ``dt``-th time step of ``cntr::herm_matrix`` object ``G`` to a given hdf5 group handle with given group name.
   * - ``G.write_to_hdf5_slices(const char *filename, const char *groupname, int dt)``
     - Write data and attributes from every ``dt``-th time step of ``cntr::herm_matrix`` object ``G`` to a hdf5 file under the given group name.

**hdf5 writing: Wigner representation:**

.. list-table::
   :header-rows: 0

   * - ``G.write_to_hdf5_tavtrel(hid_t group_id, int dt)``
     - Stores greater and lesser components of a ``cntr::herm_matrix`` object ``G`` in Wigner time representation (average and relative time) to a given hdf5 group handle.
   * - ``G.write_to_hdf5_tavtrel(hid_t group_id, const char *groupname, int dt)``
     - Stores greater and lesser components of a ``cntr::herm_matrix`` object ``G`` in Wigner time representation (average and relative time) to a given hdf5 group handle with given group name.
   * - ``G.write_to_hdf5_tavtrel(const char *filename, const char *groupname, int dt)``
     - Stores greater and lesser components of a ``cntr::herm_matrix`` object ``G`` in Wigner time representation (average and relative time) to a given HDF5 file under a specified group name.

**hdf5 reading: up to given time step:**

.. list-table::
   :header-rows: 0

   * - ``G.read_from_hdf5(int nt1, hid_t group_id)``
     - Reads the ``cntr::herm_matrix`` (attributes and data) from a given hdf5 group up to a given number of time steps ``nt1``.
   * - ``G.read_from_hdf5(int nt1, hid_t group_id, const char *groupname)``
     - Reads the ``cntr::herm_matrix`` (attributes and data) from a given hdf5 group with given group name up to a given number of time steps ``nt1``.
   * - ``G.read_from_hdf5(int nt1, const char *filename, const char *groupname)``
     - Reads the ``cntr::herm_matrix`` (attributes and data) from a given hdf5 file and given group name up to a given number of time steps ``nt1``.

**Example:**

In C++ this takes the form:

.. code-block:: cpp

   #include <cntr/cntr.hpp>
   ..
   // Create a contour Green's function
   int nt = 200, ntau = 400, norb = 1;
   GREEN A(nt, ntau, norb, FERMION);

   // Open HDF5 file and write components of the Green's function A into a group g.
   std::string filename = "data.h5";
   A.write_to_hdf5(filename.c_str(), "g");

If the file ``data.h5`` has been written previously with ``write_to_hdf5``, one can read it with the member function ``read_from_hdf5``:

.. code-block:: cpp

   // Open HDF5 file and read group g. The result is saved into the Green's function B
   GREEN B;
   B.read_from_hdf5(filename.c_str(), "g");

The parameters ``(nt,ntau,size1,sig)`` and the data of ``B`` are modified according to the information in the file (similar to reading/writing to text files discussed in :ref:`PMan01S05`).

To understand the structure of the resulting HDF5 file ``data.h5`` we inspect it with the ``h5ls`` command line program:

.. code-block:: sh

   $ h5ls -r data.h5
   ...
   /g                       Group
   /g/element_size          Dataset {1}
   /g/les                   Dataset {20301, 1, 1}
   /g/mat                   Dataset {401, 1, 1}
   /g/nt                    Dataset {1}
   /g/ntau                  Dataset {1}
   /g/ret                   Dataset {20301, 1, 1}
   /g/sig                   Dataset {1}
   /g/size1                 Dataset {1}
   /g/size2                 Dataset {1}
   /g/tv                    Dataset {80601, 1, 1}

Apart from the contour components the Green's function group ``g`` contains additional information about the dimensions and the Fermi/Bose statistics (``sig`` :math:`= \mp 1`). To understand the dimensions of the contour components we can look at the number of imaginary time steps ``ntau`` and number of real time steps ``nt`` using the ``h5dump`` command line utility:

.. code-block:: sh

   $ h5dump -d /g/ntau data.h5
   HDF5 "data.h5" {
   DATASET "/g/ntau" {
      DATATYPE  H5T_STD_I32LE
      DATASPACE  SIMPLE { ( 1 ) / ( 1 ) }
      DATA {
      (0): 400
      }
   }
   }
   $ h5dump -d /g/nt data.h5
   HDF5 "data.h5" {
   DATASET "/g/nt" {
      DATATYPE  H5T_STD_I32LE
      DATASPACE  SIMPLE { ( 1 ) / ( 1 ) }
      DATA {
      (0): 200
      }
   }
   }

This shows that the dimensions are :math:`n_\tau = 400` and :math:`n_t=200`. The size of the ``/g/mat`` component reveals that this corresponds to :math:`n_\tau + 1 = 401` imaginary time points. The mixed ``/g/tv`` component has a slow time index and a fast imaginary time index and is of size :math:`(n_t + 1)(n_\tau + 1) = 80601` while the two time triangular storage of the ``/g/ret`` and ``/g/les`` components contains :math:`(n_t + 1)(n_t + 2)/2 = 20301` elements.

To simplify postprocessing of contour Green's functions ``NESSi`` also provides the python module ``ReadCNTRhdf5.py`` for reading the HDF5 format (using the python modules ``numpy`` and ``h5py``) producing python objects with the contour components as members. For details see :ref:`PMan13S02`.
