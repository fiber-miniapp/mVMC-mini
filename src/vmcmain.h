/*-------------------------------------------------------------
 * Variational Monte Carlo
 * main program header
 *-------------------------------------------------------------
 * by Satoshi Morita and Ryui Kaneko
 *-------------------------------------------------------------*/

#ifndef _VMC_INCLUDE_FILES
#define _VMC_INCLUDE_FILES

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <errno.h>
#include <limits.h>
#include <string.h>

#ifdef _mpi_use
  #include <mpi.h>
#else
typedef int MPI_Comm;
MPI_Comm MPI_COMM_WORLD=0;
inline void MPI_Init(int argc, char* argv[]) {return;}
inline void MPI_Finalize() {return;}
inline void MPI_Abort(MPI_Comm comm, int errorcode) {exit(errorcode); return;}
inline void MPI_Barrier(MPI_Comm comm) {return;}
inline void MPI_Comm_size(MPI_Comm comm, int *size) {*size = 1; return;}
inline void MPI_Comm_rank(MPI_Comm comm, int *rank) {*rank = 0; return;}
#endif /* _mpi_use */

extern int omp_get_max_threads(void);
extern int omp_get_thread_num(void);

#include "sfmt/SFMT.h"
#include "global.h"

#include "safempi.c"
#include "time.c"
#include "workspace.c"

#ifdef _lapack
 #include "stcopt_dposv.c"
#else
 #include "stcopt_pdposv.c"
#endif

#include "gauleg.c"
#include "legendrepoly.c"
#include "splitloop.c"
#include "avevar.c"
#include "average.c"

#include "parameter.c"
#include "projection.c"
#include "slater.c"
#include "qp.c"
#include "matrix.c"
#include "pfupdate.c"
#include "pfupdate_two.c"
#include "locgrn.c"
#include "calham.c"
#include "calgrn.c"

#include "setmemory.c"
#include "readdef.c"
#include "initfile.c"

#include "vmcmake.c"
#include "vmccal.c"

#endif /* _VMC_INCLUDE_FILES */
