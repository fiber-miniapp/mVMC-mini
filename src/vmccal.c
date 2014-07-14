/*-------------------------------------------------------------
 * Variational Monte Carlo
 * calculate physical quantities
 *-------------------------------------------------------------
 * by Satoshi Morita and Ryui Kaneko
 *-------------------------------------------------------------*/

void VMCMainCal(MPI_Comm comm);
void clearPhysQuantity();
void calculateOO(double *srOptOO, double *srOptHO, const double *srOptO,
                 const double w, const double e, const int srOptSize);

void VMCMainCal(MPI_Comm comm) {
  int *eleIdx,*eleCfg,*eleNum,*eleProjCnt;
  double e,x,w,ip;

  const int qpStart=0;
  const int qpEnd=NQPFull;
  int sample,sampleStart,sampleEnd;
  int i,info;

  /* optimazation for Kei */
  const int nProj=NProj;
  double *srOptO = SROptO;

  int rank,size;
  MPI_Comm_size(comm,&size);
  MPI_Comm_rank(comm,&rank);

  SplitLoop(&sampleStart,&sampleEnd,NVMCSample,rank,size);

  /* initialization */
  clearPhysQuantity();

  for(sample=sampleStart;sample<sampleEnd;sample++) {
    eleIdx = EleIdx + sample*Nsize;
    eleCfg = EleCfg + sample*Nsite2;
    eleNum = EleNum + sample*Nsite2;
    eleProjCnt = EleProjCnt + sample*NProj;

    StartTimer(40);
    info = CalculateMAll(eleIdx,qpStart,qpEnd);
    StopTimer(40);

    if(info!=0) {
      fprintf(stderr,"waring: VMCMainCal rank:%d sample:%d info:%d (CalculateMAll)\n",rank,sample,info);
      continue;
    }

    ip = CalculateIP(PfM,qpStart,qpEnd,MPI_COMM_SELF);
    x = LogProjVal(eleProjCnt);
    /* calculate reweight */
    w = exp(2.0*(log(fabs(ip))+x) - logSqPfFullSlater[sample]);
    if( !isfinite(w) ) {
      fprintf(stderr,"waring: VMCMainCal rank:%d sample:%d w=%e\n",rank,sample,w);
      continue;
    }

    StartTimer(41);
    /* calculate energy */
    e = CalculateHamiltonian(ip,eleIdx,eleCfg,eleNum,eleProjCnt);
    StopTimer(41);
    if( !isfinite(e) ) {
      fprintf(stderr,"waring: VMCMainCal rank:%d sample:%d e=%e\n",rank,sample,e);
      continue;
    }

    Wc += w;
    Etot  += w * e;
    Etot2 += w * e * e;

    if(NVMCCalMode==0) {
      /* Calculate O for correlation fauctors */
      SROptO[0] = 1.0;
      #pragma loop noalias
      for(i=0;i<nProj;i++) srOptO[i+1] = (double)(eleProjCnt[i]);

      StartTimer(42);
      /* SlaterElmDiff */
      SlaterElmDiff(SROptO+NProj+1,ip,eleIdx);
      StopTimer(42);
      
      StartTimer(43);
      /* Calculate OO and HO */
      calculateOO(SROptOO,SROptHO,SROptO,w,e,SROptSize);

      StopTimer(43);

    } else if(NVMCCalMode==1) {
      StartTimer(42);
      /* Calculate Green Function */
      CalculateGreenFunc(w,ip,eleIdx,eleCfg,eleNum,eleProjCnt);
      StopTimer(42);
    }
  } /* end of for(sample) */

  return;
}

void clearPhysQuantity(){
  int i,n;
  double *vec;
  Wc = Etot = Etot2 = 0.0;
  if(NVMCCalMode==0) {
    /* SROptOO, SROptHO, SROptO */
    n = SROptSize*(SROptSize+2);
    vec = SROptOO;
    for(i=0;i<n;i++) vec[i] = 0.0;
  } else if(NVMCCalMode==1) {
    /* CisAjs, CisAjsCktAlt */
    n = NCisAjs+NCisAjsCktAlt;
    vec = PhysCisAjs;
    for(i=0;i<n;i++) vec[i] = 0.0;
  }
  return;
}

void calculateOO(double *srOptOO, double *srOptHO, const double *srOptO,
                 const double w, const double e, const int srOptSize) {
  double we=w*e;

  #define M_DAXPY daxpy_
  #define M_DGER dger_

  extern int M_DAXPY(const int *n, const double *alpha, const double *x, const int *incx,
                     double *y, const int *incy);
  extern int M_DGER(const int *m, const int *n, const double *alpha,
                    const double *x, const int *incx, const double *y, const int *incy, 
                    double *a, const int *lda);
  int m,n,incx,incy,lda;
  m=n=lda=srOptSize;
  incx=incy=1;

  /* OO[i][j] += w*O[i]*O[j] */
  M_DGER(&m, &n, &w, srOptO, &incx, srOptO, &incy, srOptOO, &lda);

  /* HO[i] += w*e*O[i] */
  M_DAXPY(&n, &we, srOptO, &incx, srOptHO, &incy);

  return;
}

