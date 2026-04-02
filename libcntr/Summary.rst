.. _PMan01S02Special:

Summary: Member Functions of herm_matrix and herm_matrix_timestep[_view]
=========================================================================

.. contents::
   :local:
   :depth: 2

The following table gives a brief summary of member functions of a ``herm_matrix`` and ``herm_matrix_timestep`` object ``A``, with reference to the documentation section:

.. list-table::
   :header-rows: 1

   * - Constructors and file I/O of ``herm_matrix``
     - Description
     - Reference
   * - ``herm_matrix<T>()``
     - Default constructor
     - :ref:`PMan01S02`
   * - ``herm_matrix<T>(int nt,int ntau,int size1,int sig)``
     - Constructor
     - :ref:`PMan01S02`
   * - ``void print_to_file(const char* filename)``
     - Create a textfile ``filename`` and write human readable data (large files!)
     - :ref:`PMan01S05`
   * - ``void read_from_file(const char* filename)``
     - Read a textfile ``filename`` previously created with ``print_to_file``, and resize/set ``A`` to content of file
     - :ref:`PMan01S05`
   * - ``void write_to_hdf5(const char *filename, const char *groupname)``
     - Write a data group ``groupname`` into a HDF5 file ``filename``
     - :ref:`PMan11`
   * - ``void read_from_hdf5(const char *filename, const char *groupname)``
     - Read a group ``groupname`` from HDF5 file ``filename`` and resize/set ``A`` to content of file
     - :ref:`PMan11`

.. list-table::
   :header-rows: 1

   * - Constructors of ``herm_matrix_timestep``
     - Description
     - Reference
   * - ``herm_matrix_timestep()``
     - Default constructor
     - :ref:`PMan02S02`
   * - ``herm_matrix_timestep(int nt,int ntau,int size1,int sig)``
     - Constructor
     - :ref:`PMan02S02`

.. list-table::
   :header-rows: 1

   * - Constructors of ``herm_matrix_timestep_view``
     - Description
     - Reference
   * - ``herm_matrix_timestep_view()``
     - Default constructor
     - :ref:`PManViewS02`
   * - ``herm_matrix_timestep_view(int nt,int ntau,int size1,int size2,int sig)``
     - Constructor
     - :ref:`PManViewS02`
   * - ``herm_matrix_timestep_view(int tstp, GG &g)``
     - Constructor, GG = ``herm_matrix<T>`` or ``herm_matrix_timestep<T>``
     - :ref:`PManViewS02`

.. list-table::
   :header-rows: 1

   * - Member functions of ``herm_matrix``, ``herm_matrix_timestep``, ``herm_matrix_timestep_view`` with the same syntax
     - Description
     - Reference
   * - ``void set_[XXX](time_arguments,matrix &M)``
     - Set component ``XXX=mat,ret,les,tv`` at given time arguments of the hermitian domain from a matrix or scalar M
     - :ref:`PMan01S03`
   * - ``void get_[XXX](time_arguments,matrix &M)``
     - Read component ``XXX=mat,ret,les,tv`` at given time arguments of the hermitian domain to a matrix or scalar M
     - :ref:`PMan01S03`
   * - ``void density_matrix(int tstp,matrix &M)``
     - Read density matrix at time ``tstp`` to a matrix or scalar M
     - :ref:`PMan01S04`
   * - ``void set_timestep(int tstp,GG &B)``
     - Set data at timestep ``tstp`` to data of ``B``, where ``B`` is ``herm_matrix<T>`` or ``herm_matrix_timestep<T>``
     - :ref:`PMan02S04`
   * - ``void set_matrixelement(int tstp,int i1,int i2,GG &B,int j1,int j2)``
     - Set matrixelement ``(i1,i2)`` at timestep ``tstp`` to matrixelement ``(j1,j2)`` of ``B``, where ``B`` is ``herm_matrix<T>`` or ``herm_matrix_timestep<T>``
     - :ref:`PMan02S04`
   * - ``void set_timestep_zero(int tstp)``
     - Set elements at timestep ``tstp`` to ``0``
     - :ref:`PMan04`
   * - ``void smul(int tstp, std::complex<T> x)``
     - :math:`A(t,t')` set to :math:`xA(t,t')` on timestep ``tstp`` of ``A``
     - :ref:`PMan04`
   * - ``void incr_timestep(int tstp,GG &B, std::complex<T> x)``
     - :math:`A(t,t')` set to :math:`A(t,t')+xB(t,t')` on timestep ``tstp`` of ``A``, where ``B`` is ``herm_matrix<T>`` or ``herm_matrix_timestep<T>``
     - :ref:`PMan04`
   * - ``void left_multiply(int tstp,cntr::function<T> &f, T x)``
     - :math:`A(t,t')` set to :math:`x f(t) A(t,t')` at timestep ``tstp`` of ``A``
     - :ref:`PMan04S02`
   * - ``void right_multiply(int tstp,cntr::function<T> &f, T x)``
     - :math:`A(t,t')` set to :math:`x A(t,t') f(t')` at timestep ``tstp`` of ``A``
     - :ref:`PMan04S02`
   * - ``void left_multiply_hermconj(int tstp,cntr::function<T> &f, T x)``
     - :math:`A(t,t')` set to :math:`x f(t)^\dagger A(t,t')` at timestep ``tstp`` of ``A``
     - :ref:`PMan04S02`
   * - ``void right_multiply_hermconj(int tstp,cntr::function<T> &f, T x)``
     - :math:`A(t,t')` set to :math:`x A(t,t') f(t')^\dagger` at timestep ``tstp`` of ``A``
     - :ref:`PMan04S02`
   * - ``void Bcast_timestep(int tstp, int root)``
     - Broadcast the data of :math:`A(t,t')` on timestep ``tstp`` from rank ``root`` to all other ranks.
     - :ref:`PMan10S01`
   * - ``void Reduce_timestep(int tstp, int root)``
     - In place reduction, replacing :math:`A(t,t')` on rank ``root`` by :math:`\sum_{\text{rank}\,j} A_{\text{at rank}\,j}(t,t')` for one timestep ``tstp``.
     - :ref:`PMan10S01`
   * - ``void Send_timestep(int tstp, int dest, int tag)``
     - Send the data of :math:`A(t,t')` on timestep ``tstp`` to rank ``dest`` with a message tag ``tag``.
     - :ref:`PMan10S01`
   * - ``void Recv_timestep(int tstp, int root, int tag)``
     - Receive a message with tag ``tag`` from rank ``root`` which contains the data of :math:`A(t,t')` on timestep.
     - :ref:`PMan10S01`
