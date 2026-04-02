.. _PMan13:

Python Tools
============

.. contents::
   :local:
   :depth: 2

As a part of the ``NESSi`` package, we provide python tools, which include:

- Scripts for pre-processing to assist the use of programs based on ``libcntr``
- Scripts for reading and post-processing Green's functions from HDF5 format via the ``h5py`` python package.

The python tools can be found in ``libcntr/python`` (for python2) or ``libcntr/python3`` (for python3 compatibility).

.. _PMan13S01:

Creation of an input file
--------------------------

An input file for a custom program (see details in :ref:`PMan14`) can be generated using ``ReadCNTR.py`` python module:

**Example:**

.. code-block:: python

   from ReadCNTR import write_input_file

   inp = {
     'nt': 1000,                 # number of points on the real axis
     'ntau': 500,                # number of points on the Matsubara (imaginary) axis
     'h': 0.01,                  # timestep interval
     'beta': 10.0,               # inverse temperature
     'Efield': '--field.txt'     # electrical field read from a file
   }

   write_input_file(inp, 'input.txt')      # creates an input file

.. note::

   To make the python module available to a python script, the module ``ReadCNTR.py`` needs to be in the python include path. For python2, this is achieved by setting ``export PYTHONPATH=<path>/nessi/libcntr/python``, while ``export PYTHONPATH=<path>/nessi/libcntr/python3`` makes the python3 version accessible.

.. _PMan13S02:

Reading Green's functions from hdf5
-------------------------------------

The python module unrolls the triangular storage of the ``ret`` and ``les`` components making it simple to work with time slices.

**Example 1:**

To store the imaginary part of the retarded Green's function :math:`\textrm{Im}[G_{i,j}^\mathrm{R}(t, t^\prime=5 h)]` at orbital indices :math:`i=j=0` as a function of :math:`t` into a text file we may use the commands:

.. code-block:: python

   import h5py
   from ReadCNTRhdf5 import read_group

   # read the components of the contour function from hdf5 file
   with h5py.File('data.h5', 'r') as fd:
       g = read_group(fd).g

   # choose a timeslice for the retarded Green's function
   indxt=5
   G_ret=g.ret[:,indxt,0,0].imag
   with open("data_output.txt",'w') as output:     # create a text file
       output.write("%s \t" % G_ret)               # write data into it

**Example 2:**

The same way one can save the lesser component of the Green's function :math:`\textrm{Im}[G^<(t=5 h, t^\prime)]` as a function of :math:`t^\prime` into a text file:

.. code-block:: python

   import h5py
   from ReadCNTRhdf5 import read_group

   # read the components of the contour function from hdf5 file
   with h5py.File('data.h5', 'r') as fd:
       g = read_group(fd).g

   # choose a timeslice for the lesser Green's function
   indxt=5
   G_les=g.les[indxt,:,0,0].imag
   with open("data_output.txt",'w') as output:     # create a text file
       output.write("%s \t" % G_les)               # write data into it

.. _PMan13S03:

Postprocessing analysis
------------------------

The ``ReadCNTRhdf5`` python module contains also functions for reading parameters and observables from an output file. The following minimal example script is used to interpret the output files in the examples provided in section :ref:`PEx`.

**Example:**

.. code-block:: python

   import h5py
   from ReadCNTRhdf5 import read_imp_h5file
   # Plot

   ax0=plt.subplot(1,1,1)
   imp_filename = '{}/data.h5'.format(output_file)
   data = read_imp_h5file(imp_filename,['obs','parm'])

   Nt=data.parm.Nt                                                       # read out the number of timesteps on the real axis from the output file
   dt=data.parm.h                                                        # read out the timestep interval from the output file
   tpts = np.linspace(0.0,Nt*dt,Nt+1)                                    # define the real time axis
   ekin=np.real(data.obs.Ekin.data[1:,0,0]-data.obs.Ekin.data[0,0,0])    # read out the value of an observable (Ekin) from the output file and calculate the difference between the equilibrium and nonequilibrium values

   plt.plot(tpts,ekin,label="Ekin")                                      # plot data

   ax0.set_xlabel(r't ')                                                 # set x label
   ax0.set_ylabel(r"$\mathrm{E}_\mathrm{kin}$")                          # set y label

   ax0.legend(loc="best",bbox_to_anchor=(0.45,0.3),frameon=False,fancybox=False,handlelength=1)
   plt.savefig('Ekin.pdf')
   plt.show()
