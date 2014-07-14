/*-------------------------------------------------------------
 * Variational Monte Carlo
 * calculate weighted averages of physical quantities
 *-------------------------------------------------------------
 * by Satoshi Morita
 *-------------------------------------------------------------*/
void WeightAverageWE(MPI_Comm comm);
void WeightAverageSROpt(MPI_Comm comm);
void WeightAverageGreenFunc(MPI_Comm comm);

void weightAverageReduce(int n, double *vec, MPI_Comm comm);

/* calculate average of Wc, Etot and Etot2 */
/* All processes will have the result */
void WeightAverageWE(MPI_Comm comm) {
  const int n=3;
  double invW;
  int rank,size;
  double send[n], recv[n];
  MPI_Comm_rank(comm,&rank);
  MPI_Comm_size(comm,&size);

  /* Wc, Etot and Etot2 */
  if(size>1) {
    send[0] = Wc;
    send[1] = Etot;
    send[2] = Etot2;

    SafeMpiAllReduce(send,recv,n,comm);

    Wc    = recv[0];
    invW  = 1.0/Wc;
    Etot  = recv[1]*invW;
    Etot2 = recv[2]*invW;
  } else {
    invW  = 1.0/Wc;
    Etot  *= invW;
    Etot2 *= invW;
  }

  return;
}

/* calculate average of SROptOO and SROptHO */
/* All processes will have the result */
void WeightAverageSROpt(MPI_Comm comm) {
  int i,n;
  double invW = 1.0/Wc;
  double *vec,*buf;
  int rank,size;
  MPI_Comm_rank(comm,&rank);
  MPI_Comm_size(comm,&size);

  /* SROptOO and SROptHO */
  n = SROptSize*(SROptSize+1);
  vec = SROptOO;
  if(size>1) {
    RequestWorkSpaceDouble(n);
    buf = GetWorkSpaceDouble(n);

    SafeMpiAllReduce(vec,buf,n,comm);

    #pragma omp parallel for default(shared) private(i)
    #pragma loop noalias
#pragma ivdep
    for(i=0;i<n;i++) vec[i] = buf[i] * invW;

    ReleaseWorkSpaceDouble();
 } else {
    #pragma omp parallel for default(shared) private(i)
    #pragma loop noalias
#pragma ivdep
    for(i=0;i<n;i++) vec[i] *= invW;
  }

  return;
}

/* calculate average of Green functions */
/* Only rank=0 process will have the result */
void WeightAverageGreenFunc(MPI_Comm comm) {
  int n;
  double *vec;

  /* Green functions */
  /* CisAjs and CisAjsCktAlt */
  n = NCisAjs+NCisAjsCktAlt;
  vec = PhysCisAjs;
  weightAverageReduce(n,vec,comm);
  
  return;
}

void weightAverageReduce(int n, double *vec, MPI_Comm comm) {
  int i;
  const double invW = 1.0/Wc;
  double *buf;
  int rank,size;
  MPI_Comm_rank(comm,&rank);
  MPI_Comm_size(comm,&size);

  if(size>1) {
    RequestWorkSpaceDouble(n);
    buf = GetWorkSpaceDouble(n);

    SafeMpiReduce(vec,buf,n,comm);
    if(rank==0) {
      #pragma omp parallel for default(shared) private(i)
      #pragma loop noalias
#pragma ivdep
      for(i=0;i<n;i++) vec[i] = buf[i] * invW;
    }

    ReleaseWorkSpaceDouble();
  } else {
    #pragma omp parallel for default(shared) private(i)
    #pragma loop noalias
#pragma ivdep
    for(i=0;i<n;i++) vec[i] *= invW;
  }

  return;
}
