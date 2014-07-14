mVMC-mini
=========

* version: 1.0.0 (based on mVMC 2013.11.08 version)
* date: 2014/07/17
* contact: miniapp@riken.jp


About mVMC-mini and mVMC
-------------------------

The mVMC-mini package is based on mVMC simulation program.
The primary purpose of mVMC-mini is to enable the performance study of mVMC
on various platforms with less effort using compact installation procedure
while representing a typical computational work load of mVMC.
It contains the mVMC source program and the test data as described in the
following sections.

Regarding the license condition to use this mVMC-mini package, there are
several related files.
Please refer to the "LICENSE" file included in the top directory of
this package.
Please also refer to the two "LICENSE" files each included in
src/pfapack/ and src/sfmt/ subdirectories which are used as part of mVMC.


Installation
------------

The mVMC-mini package "mVMC-mini.v1.tar.gz" can be downloaded from the 
repository.  
https://github.com/fiber-miniapp/mVMC-mini  
The following software suite is necessary for installing the package.

+ Prerequisite software:
  - C compiler and Fortran compiler with OpenMP support
  - BLACS library.
  - Scalapack library.
  - MPI library.

####example - K computer
As a quick start, the installation example on K computer is shown first.
The installation should be quite easy as below.

    $ tar -zxf mVMC-mini.v1.tar.gz
    $ cd src
    $ make Kei

When the installation completes successfully, the executable file is created
as "vmc.out" in the src/ directory as below.

	$ ls -go vmc.out 
	-rwxr-xr-x 1 463739 Jun 12 17:28 vmc.out

####example - Intel compiler and libraries
The installation example on Intel software environment is shown next.

    $ tar -zxf mVMC-mini.v1.tar.gz
    $ cd src
    $ make intel


####General installation steps on other Linux platform

The installation of mVMC-mini on other Linux based platforms should be fairly
simple as well.  A typical installation step is explained below.

#####step 1.

Obtain the mVMC-mini package "mVMC-mini.v1.tar.gz" from the repository.
https://github.com/fiber-miniapp/mVMC-mini  
Extract its contents using tar command.
This readme.txt (readme.md) file is included in the package.

	$ tar -zxf mVMC-mini.v1.tar.gz
	$ ls
	LICENSE    job_middle  makeDef           result
	README.md  job_tiny    readme.asis.utf8  src

The directories contain the following files.  

+	job_middle/	# medium size test job definition directory
+	job_tiny/	# tiny size test job definition directory
+	makeDef/	# contains Python script to produce the job definition file
+	result/		# contains computed results on a reference platform (FX10)
+	src/		# source directory, with following additional subdirectories
 -	pfapack/ # subdir. containing Pfaffian computation library
 -	sfmt/ # subdir. containing SIMD-oriented Fast Mersenne Twister


In src/ directory, there are "Makefile_*" for several platforms.

+ Predefined platform : Makefile_${platform}
 -	Makefile_intel : for Intel compiler+MPI on Intel Xeon Linux cluster
 -	Makefile_fx10 : for Fujitsu compiler+MPI on FX10 and K computer
 -	Makefile_pgi  : for PGI compiler + OpenMPI/pgi on Intel Xeon Linux
 -	Makefile_gnu  : for GNU compiler + OpenMPI/gcc on Intel Xeon Linux

If your testing platform is covered by one of these, set the value of platform,
and you can move to the step 3, without taking step 2.  
For example of intel, set the value as below, and you can move to step 3.

	$ platform=intel

If your testing platform is NOT covered by above, create your
"Makefile_${platform}", by following the step 2.


#####step 2.

Create "Makefile_${platform}" _file_.

If you are installing mVMC on a different platform, then create a file
"Makefile_${platform}" accordingly.
The value for ${platform} can be any, such as akb.
It may be easier to copy Makefile_skelton and edit it.
Edit the file and give some commonly used compiler options via CFLAGS.
In addition to the library software included in the package, mVMC-mini
also requires the **BLACS** library and **Scalapack** library,
whose location is system dependent.
The location of these BLACS and Scalapack library should be specified via
LIB setting in the "Makefile_${platform}" file.

	$ cd src
	$ platform=akb
	$ cp Makefile_intel Makefile_${platform}
	$ vi Makefile_${platform}		# edit CFLAGS and LIB, etc.

#####step 3.

Run make command in the src/ directory. 
After make command finishes successfully, there should be an executable file 
named "vmc.out" in the directory.

	$ make -f Makefile_${platform}
	$ ls -go vmc.out 
	-rwxr-xr-x 1 463739 Jun 12 17:28 vmc.out


How to run an example job
----------------------

#### Test job directory and job script
A couple of test directories are included in the package.

+ job_middle/
	Medium size test job definition directory.
	It takes 4 minutes on K computer 128 nodes (128MPIx8OpenMP)
+ job_tiny/
	Tiny size test job definition directory.
	It takes 12 seconds on Intel E5-2420. Mostly used for a quick debug.

There is a shell script file "job.sh" in each directory.
The script can be run in foreground.
It can be run as a batch job with appropriate batch directives added.

After the successful job execution, the output files named "zvo_*.dat"
are saved under a directory. The name of the directory is defined in the
input data file "multiDir.def".
In the stdout, there will be an error message printed as >
"Error: opt.init does not exist."
which can be safely ignored.

For medium size test job, example batch job script files for
K computer, FX10 and Intel cluster are provided as:

- "job-K.sh"
- "job-FX10.sh"
- "job-Intel.sh"

On K computer, the input/output files must be staged.
Please note that the stgin/stgout-basedir and stgin/stgout directives
in "job-K.sh" must be modified according to the installed path.
The jobscript lines starting with:

    #PJM --stgin-basedir
    #PJM --stgout-basedir

will have to be changed to reflect the base directories and files accordingly.

The simulated model is 2D square Kondo lattice model (J/t=1.0, half-filling)
and 20 steps of parameter optimization is carried out.
The number of parameter optimization steps is defined by
NSROptItrStep in the "zmodpara.def" file.
The computing time should be proportional to this value, so it can be used
to control the job elapse time.

#### Verifying the computed results

The previously computed result files on a reference platform (FX10)
are also saved under result/ directory for comparison.
Verifying the numerical results can be best done by checking the values
in "zvo_out_000.dat" file.
Its first column shows the value of the energy expectation.
For jobs with the same number of MPI processes, the computed values
should roughly match.


Testing at scale
-----------------------

mVMC adopts two different parallel processing approaches for its
computing phases.
For the parameter optimization phase, it calls ScaLAPACK using all the
MPI processes.
For the Monte Carlo computation phase, it divides the MPI processes into
NSplitSize groups specified in zmodpara.def file, and each group
produces NVMCSample Monte Carlo samples and computes physical variables.
The computed values from all the processes are gathered to calculate
the expectation.
If a job is configured to run with Nmpi processes, the total number of
the Monte Carlo samples will be (Nmpi/NSplitSize)*NVMCSample .
Nmpi is the number of MPIs, which is not included in the zmodpara.def file,
but is declared as mpirun command option.

#### weak scaling test
The default setting will effectively define weak scaling test variations.
Copy the directory job_middle to another one, and change the number of MPIs
to run weak scaling test.

The parameters in the default job_middle zmodpara.def includes
NSplitSize=4 and NVMCSample=192, so if a job is run using 4 MPIs,
the number of Monte Carlo sampling will be 192, and if a job is run using
8 MPIs, the number will be 384, etc.

#### strong scaling test
For strong scaling tests, the number of Monte Carlo samples should be kept
the same regardless of Nmpi, which can done by changing the number of groups
NSplitSize according to Nmpi.
The default job_middle parameters, however, is not suitable for
scaling tests up to large scale.

A separate Python script named "makeDef_large.py" is provided
under makeDef/ directory, which will produce the input definition files
and the jobscript.
To produce the input definition files and the jobscript for Nmpi process job,
run the following command on a system which supports Python.

    ./makeDef_large.py Nmpi

After this script is run, there will be a jobscript "job_mpi${Nmpi}" and
a set of input definition files named "*.def".
The number of groups are set as 64, so jobs with the multiple of 64 MPIs
are good choices.  The jobs up to 4096 MPIs can be created.
A job with 1024 MPIs takes about one hour.
Jobs with small number of MPIs will run many hours, accordingly.
When running strong scaling tests, it is generally recommended to
reduce the number of steps, i.e. NSROptItrStep.

#### OpenMP threads
The number of OpenMP threads can be set independent of the input difinition
files. If the number of threads exceeds (NSPGaussLeg*NMPTrans)/NSplitSize,
then the thread load imbalance may occur.
By default, the job.sh file created by above makeDef_large.py sets the
threads value as:

    export OMP_NUM_THREADS=8



What mVMC-mini does - brief explanation
---------------------------------------

mVMC-mini is a subset of mVMC full application, and essentially retains
the same feature as mVMC.
The original mVMC has the following features.
Using the multi-variable variational parameters, mVMC analyzes the
physical characteristics of the strongly correlated electron systems,
by configuring the variational wave functions closely representing
the ground state of such electron systems.  It executes the Monte Carlo
sampling of the real space configuration of the electrons.

The variational wave function in mVMC is composed of three parts, the
singlet paring wave function, correlation factors, and quantum-number
projections. mVMC computes physical properties of the variational wave
function in the strongly correlated electron systems, such as energy,
magnetism and superconductivity. To carry out such computation, it
computes Pfaffian of the skew-symmetric matrix to obtain the inner
product between the paring wave function and the real space
electron configuration.

In optimizing the variational parameters, it applies the stochastic
reconfiguration method, which is more stable than the steepest descent
method.  The Monte Carlo sampling and the matrix computation process
in the optimization process are both parallelized.

For detail explanation of mVMC, refer to the paper:
Tahara D and Imada M, J. Phys. Soc. Jpn. 77, 114701 (2008)
"Variational Monte Carlo method combined with quantum-number
projection and multi-variable optimization."

Contact point of the original mVMC:
Dr. Satoshi Morita <morita@issp.u-tokyo.ac.jp>


Target exa-scale problem setting
--------------------------------

The milestone performance of the 10,000 atom system simulation using mVMC
is expected to be 8 hours on exa-scale platform. Such milestone is
under review.


