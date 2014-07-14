/*-------------------------------------------------------------
 * Variational Monte Carlo
 * Variational parameters
 *-------------------------------------------------------------
 * by Satoshi Morita
 *-------------------------------------------------------------*/
#define D_AmpMax             4.0

void InitParameter();
int ReadInitParameter(char *initFile);
void SyncModifiedParameter(MPI_Comm comm);

void shiftGJ();

/* initialize variational parameters */
void InitParameter() {
  int i;

  #pragma omp parallel for default(shared) private(i)
  for(i=0;i<NProj;i++) Proj[i] = 0.0;

  for(i=0;i<NSlater;i++){
    if(OptFlag[i+NProj] > 0){
      Slater[i] = genrand_real2(); /* uniform distribution [0,1) */
    } else {
      Slater[i] = 0.0;
    }
  }

  return;
}

/* read initial vaules of variational parameters from initFile */
int ReadInitParameter(char *initFile) {
  FILE *fp;
  int i,xi;
  double xtmp;

  fp = fopen(initFile, "r");
  if(fp!=NULL){
    while(fscanf(fp, "%lf ", &xtmp)!=EOF){
      for(i=1;i<4;i++) fscanf(fp, "%lf ", &xtmp);
      for(xi=0;xi<NProj;xi++) {
        fscanf(fp, "%lf %lf ", &(Proj[xi]), &xtmp);
      }
      for(xi=0;xi<NSlater;xi++) {
        fscanf(fp, "%lf %lf ", &(Slater[xi]), &xtmp);
      }
    }
    fclose(fp);
  } else { fprintf(stderr, "Error: %s does not exist.\n",initFile); }
  
  return 0;
}

/* sync and modify variational parameters */
void SyncModifiedParameter(MPI_Comm comm) {
  double xmax,ratio;
  int i;

#ifdef _mpi_use
  int size;
  MPI_Comm_size(comm, &size);
  if(size>1) MPI_Bcast(Para, NPara, MPI_DOUBLE, 0, comm);
#endif /* _mpi_use */

  /***** shift correlation factors *****/
  shiftGJ();
 
  /***** rescale Slater *****/
  xmax = fabs(Slater[0]);
  for(i=1;i<NSlater;i++){
    if(xmax < fabs(Slater[i])) xmax = fabs(Slater[i]);
  }
  ratio = D_AmpMax/xmax;
  #pragma omp parallel for default(shared) private(i)
  for(i=0;i<NSlater;i++) Slater[i] *= ratio;

  return;
}

/* shift Gutzwiller-Jastrow factor */
void shiftGJ() {
  double shift=0.0;
  const int n = NGutzwillerIdx+NJastrowIdx;
  int i;

  if(NGutzwillerIdx==0||NJastrowIdx==0) return;

  for(i=0;i<n;i++) {
    shift += Proj[i];
  }
  shift /= (double)n;

  for(i=0;i<n;i++) {
    Proj[i] -= shift;
  }

  return;
}
