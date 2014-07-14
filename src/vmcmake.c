/*-------------------------------------------------------------
 * Variational Monte Carlo
 * make sample
 *-------------------------------------------------------------
 * by Satoshi Morita
 *-------------------------------------------------------------*/

void VMCMakeSample(MPI_Comm comm);
int makeInitialSample(int *eleIdx, int *eleCfg, int *eleNum, int *eleProjCnt,
                      const int qpStart, const int qpEnd, MPI_Comm comm);
void copyFromBurnSample(int *eleIdx, int *eleCfg, int *eleNum, int *eleProjCnt);
void copyToBurnSample(const int *eleIdx, const int *eleCfg, const int *eleNum, const int *eleProjCnt);
void saveEleConfig(const int sample, const double logIp,
                   const int *eleIdx, const int *eleCfg, const int *eleNum, const int *eleProjCnt);
void sortEleConfig(int *eleIdx, int *eleCfg, const int *eleNum);

void ReduceCounter(MPI_Comm comm);

void VMCMakeSample(MPI_Comm comm) {
  int outStep,nOutStep;
  int inStep,nInStep;
  int updateType;
  int mi,mj,ri,rj,s,t,i;
  int nAccept=0;
  int sample;

  double logIpOld,logIpNew; /* logarithm of inner product <phi|L|x> */
  int projCntNew[NProj];
  double pfMNew[NQPFull];
  double x,w;

  int qpStart,qpEnd;
  int rank,size;
  MPI_Comm_size(comm,&size);
  MPI_Comm_rank(comm,&rank);

  SplitLoop(&qpStart,&qpEnd,NQPFull,rank,size);

  StartTimer(30);
  if(BurnFlag==0) {
    makeInitialSample(TmpEleIdx,TmpEleCfg,TmpEleNum,TmpEleProjCnt,
                      qpStart,qpEnd,comm);
  } else {
    copyFromBurnSample(TmpEleIdx,TmpEleCfg,TmpEleNum,TmpEleProjCnt);
  }
  
  CalculateMAll(TmpEleIdx,qpStart,qpEnd);
  logIpOld = CalculateLogIP(PfM,qpStart,qpEnd,comm);
  if( !isfinite(logIpOld) ) {
    if(rank==0) fprintf(stderr,"waring: VMCMakeSample remakeSample logIpOld=%e\n",logIpOld);
    makeInitialSample(TmpEleIdx,TmpEleCfg,TmpEleNum,TmpEleProjCnt,
                      qpStart,qpEnd,comm);
    CalculateMAll(TmpEleIdx,qpStart,qpEnd);
    logIpOld = CalculateLogIP(PfM,qpStart,qpEnd,comm);
    BurnFlag = 0;
  }
  StopTimer(30);

  nOutStep = (BurnFlag==0) ? NVMCWarmUp+NVMCSample : NVMCSample+1;
  nInStep = NVMCIniterval * Nsite;

  for(i=0;i<4;i++) Counter[i]=0;  /* reset counter */

  for(outStep=0;outStep<nOutStep;outStep++) {
    for(inStep=0;inStep<nInStep;inStep++) {

      StartTimer(31);
      /* generate the candidate configuration */
      updateType = -1;
      while(updateType<0) {
        updateType=0; /* hopping */
        if(NExUpdatePath==1) {
          if(genrand_real2()<0.5) updateType = 1; /* exchange */
        }
        mi = gen_rand32()%Ne;
        s = (genrand_real2()<0.5) ? 0 : 1;
        ri = TmpEleIdx[mi+s*Ne];
        do {
          rj = gen_rand32()%Nsite;
        } while (TmpEleCfg[rj+s*Nsite] != -1);

        /* check */
        if(updateType==0) {
          /* The mi-th electron with spin s hops to site rj */
          if(LocSpn[ri]==0 || LocSpn[rj]==0) updateType=-1;
        } else {
          /* The mi-th electron with spin s exchanges with the electron on site rj with spin 1-s */
          if(TmpEleCfg[rj+(1-s)*Nsite] == -1 || TmpEleCfg[ri+(1-s)*Nsite] != -1) updateType=-1;
        }
      }
      StopTimer(31);
      
      if(updateType==0) { /* hopping */
        StartTimer(32);
        /* The mi-th electron with spin s hops to site rj */
        TmpEleIdx[mi+s*Ne] = rj;
        TmpEleCfg[ri+s*Nsite] = -1;
        TmpEleCfg[rj+s*Nsite] = mi;
        TmpEleNum[ri+s*Nsite] = 0;
        TmpEleNum[rj+s*Nsite] = 1;

          StartTimer(60);
        UpdateProjCnt(ri,rj,s,projCntNew,TmpEleProjCnt,TmpEleNum);
          StopTimer(60);
          StartTimer(61);
        CalculateNewPfM2(mi,s,pfMNew,TmpEleIdx,qpStart,qpEnd);
          StopTimer(61);

          StartTimer(62);
        /* calculate inner product <phi|L|x> */
        logIpNew = CalculateLogIP(pfMNew,qpStart,qpEnd,comm);
          StopTimer(62);

        /* Metroplis */
        x = LogProjRatio(projCntNew,TmpEleProjCnt);
        w = exp(2.0*(x+logIpNew-logIpOld));
        if( !isfinite(w) ) w = -1.0; /* should be rejected */

        if(w > genrand_real2()) { /* accept */
            StartTimer(63);
            UpdateMAll(mi,s,TmpEleIdx,qpStart,qpEnd);
            StopTimer(63);

          for(i=0;i<NProj;i++) TmpEleProjCnt[i] = projCntNew[i];
          logIpOld = logIpNew;
          nAccept++;
          Counter[1]++;
        } else { /* reject */
          TmpEleIdx[mi+s*Ne] = ri;
          TmpEleCfg[ri+s*Nsite] = mi;
          TmpEleCfg[rj+s*Nsite] = -1;
          TmpEleNum[ri+s*Nsite] = 1;
          TmpEleNum[rj+s*Nsite] = 0;
        }
        Counter[0]++;
        StopTimer(32);
      } else { /* exchange */
        StartTimer(33);
        StartTimer(65);
        /* The mi-th electron with spin s exchanges with the electron on site rj with spin 1-s */
        t = 1-s;
        mj = TmpEleCfg[rj+t*Nsite];
        /* The mi-th electron with spin s hops to rj */
        TmpEleIdx[mi+s*Ne] = rj;
        TmpEleCfg[ri+s*Nsite] = -1;
        TmpEleCfg[rj+s*Nsite] = mi;
        TmpEleNum[ri+s*Nsite] = 0;
        TmpEleNum[rj+s*Nsite] = 1;
        UpdateProjCnt(ri,rj,s,projCntNew,TmpEleProjCnt,TmpEleNum);
        /* The mj-th electron with spin t hops to ri */
        TmpEleIdx[mj+t*Ne] = ri;
        TmpEleCfg[rj+t*Nsite] = -1;
        TmpEleCfg[ri+t*Nsite] = mj;
        TmpEleNum[rj+t*Nsite] = 0;
        TmpEleNum[ri+t*Nsite] = 1;
        UpdateProjCnt(rj,ri,t,projCntNew,projCntNew,TmpEleNum);

        StopTimer(65);
        StartTimer(66);

        CalculateNewPfMTwo2(mi, s, mj, t, pfMNew, TmpEleIdx, qpStart, qpEnd);

        StopTimer(66);
        StartTimer(67);

        /* calculate inner product <phi|L|x> */
        logIpNew = CalculateLogIP(pfMNew,qpStart,qpEnd,comm);
        
        StopTimer(67);

        /* Metroplis */
        x = LogProjRatio(projCntNew,TmpEleProjCnt);
        w = exp(2.0*(x+logIpNew-logIpOld));
        if( !isfinite(w) ) w = -1.0; /* should be rejected */
        
        if(w > genrand_real2()) { /* accept */
          StartTimer(68);
          UpdateMAllTwo(mi, s, mj, t, ri, rj, TmpEleIdx,qpStart,qpEnd);
          StopTimer(68);

/*           printf("%d %lf %lf\n",inStep,InvM[(mi+s*Ne)*Nsize+mj+t*Ne],InvM[(mi+s*Ne)*Nsize+mj+t*Ne+1]); */
/*           CalculateMAll(TmpEleIdx,qpStart,qpEnd); */
/*           printf("%d %lf %lf\n",inStep,InvM[(mi+s*Ne)*Nsize+mj+t*Ne],InvM[(mi+s*Ne)*Nsize+mj+t*Ne+1]); */

          for(i=0;i<NProj;i++) TmpEleProjCnt[i] = projCntNew[i];
          logIpOld = logIpNew;
          nAccept++;
          Counter[3]++;
        } else { /* reject */
          TmpEleIdx[mj+t*Ne] = rj;
          TmpEleCfg[rj+t*Nsite] = mj;
          TmpEleCfg[ri+t*Nsite] = -1;
          TmpEleNum[rj+t*Nsite] = 1;
          TmpEleNum[ri+t*Nsite] = 0;

          TmpEleIdx[mi+s*Ne] = ri;
          TmpEleCfg[ri+s*Nsite] = mi;
          TmpEleCfg[rj+s*Nsite] = -1;
          TmpEleNum[ri+s*Nsite] = 1;
          TmpEleNum[rj+s*Nsite] = 0;
        }
        Counter[2]++;
        StopTimer(33);
      }

      if(nAccept>Nsite) {
        StartTimer(34);
        /* recal PfM and InvM */
        CalculateMAll(TmpEleIdx,qpStart,qpEnd);
        logIpOld = CalculateLogIP(PfM,qpStart,qpEnd,comm);
        StopTimer(34);
        nAccept=0;
      }
    } /* end of instep */

    StartTimer(35);
    /* save Electron Configuration */
    if(outStep >= nOutStep-NVMCSample) {
      sample = outStep-(nOutStep-NVMCSample);
      saveEleConfig(sample,logIpOld,TmpEleIdx,TmpEleCfg,TmpEleNum,TmpEleProjCnt);
    }
    StopTimer(35);

  } /* end of outstep */

  copyToBurnSample(TmpEleIdx,TmpEleCfg,TmpEleNum,TmpEleProjCnt);
  BurnFlag=1;
  return;
}

int makeInitialSample(int *eleIdx, int *eleCfg, int *eleNum, int *eleProjCnt,
                      const int qpStart, const int qpEnd, MPI_Comm comm) {
  const int nsize = Nsize;
  const int nsite2 = Nsite2;
  int flag=1,flagRdc,loop=0;
  int ri,mi,si,msi,rsi;
  int rank,size;
  MPI_Comm_size(comm,&size);
  MPI_Comm_rank(comm,&rank);

  do {
    /* initialize */
    #pragma omp parallel for default(shared) private(msi)
    for(msi=0;msi<nsize;msi++) eleIdx[msi] = -1;
    #pragma omp parallel for default(shared) private(rsi)
    for(rsi=0;rsi<nsite2;rsi++) eleCfg[rsi] = -1;
    
    /* local spin */
    for(ri=0;ri<Nsite;ri++) {
      if(LocSpn[ri]==0) {
        do {
          mi = gen_rand32()%Ne;
          si = (genrand_real2()<0.5) ? 0 : 1;
        } while(eleIdx[mi+si*Ne]!=-1);
        eleCfg[ri+si*Nsite] = mi;
        eleIdx[mi+si*Ne] = ri;
      }
    }
    
    /* itinerant electron */
    for(si=0;si<2;si++) {
      for(mi=0;mi<Ne;mi++) {
        if(eleIdx[mi+si*Ne]== -1) {
          do {
            ri = gen_rand32()%Nsite;
          } while (eleCfg[ri+si*Nsite]!= -1 || LocSpn[ri]==0);
          eleCfg[ri+si*Nsite] = mi;
          eleIdx[mi+si*Ne] = ri;
        }
      }
    }
    
    /* EleNum */
    #pragma omp parallel for default(shared) private(rsi)
    #pragma loop noalias
    for(rsi=0;rsi<nsite2;rsi++) {
      eleNum[rsi] = (eleCfg[rsi] < 0) ? 0 : 1;
    }
    
    MakeProjCnt(eleProjCnt,eleNum);

    flag = CalculateMAll(eleIdx,qpStart,qpEnd);
    if(size>1) {
      MPI_Allreduce(&flag,&flagRdc,1,MPI_INT,MPI_MAX,comm);
      flag = flagRdc;
    }

    loop++;
    if(loop>100) {
      if(rank==0) fprintf(stderr, "error: makeInitialSample: Too many loops\n");
      MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
    }
  } while (flag>0);
 
  return 0;
}

void copyFromBurnSample(int *eleIdx, int *eleCfg, int *eleNum, int *eleProjCnt) {
  int i,n;
  const int *burnEleIdx = BurnEleIdx;
  n = Nsize + 2*Nsite + 2*Nsite + NProj;
  #pragma loop noalias
  for(i=0;i<n;i++) eleIdx[i] = burnEleIdx[i];
  return;
}

void copyToBurnSample(const int *eleIdx, const int *eleCfg, const int *eleNum, const int *eleProjCnt) {
  int i,n;
  int *burnEleIdx = BurnEleIdx;
  n = Nsize + 2*Nsite + 2*Nsite + NProj;
  #pragma loop noalias
  for(i=0;i<n;i++) burnEleIdx[i] = eleIdx[i];
  return;
}

void saveEleConfig(const int sample, const double logIp,
                   const int *eleIdx, const int *eleCfg, const int *eleNum, const int *eleProjCnt) {
  int i,offset;
  double x;
  const int nsize=Nsize;
  const int nsite2 = Nsite2;
  const int nProj = NProj;

  offset = sample*nsize;
  #pragma loop noalias
  for(i=0;i<nsize;i++) EleIdx[offset+i] = eleIdx[i];
  offset = sample*nsite2;
  #pragma loop noalias
  for(i=0;i<nsite2;i++) EleCfg[offset+i] = eleCfg[i];
  #pragma loop noalias
  for(i=0;i<nsite2;i++) EleNum[offset+i] = eleNum[i];
  offset = sample*nProj;
  #pragma loop noalias
  for(i=0;i<nProj;i++) EleProjCnt[offset+i] = eleProjCnt[i];
  
  x = LogProjVal(eleProjCnt);
  logSqPfFullSlater[sample] = 2.0*(x+logIp);
  
  return;
}

void sortEleConfig(int *eleIdx, int *eleCfg, const int *eleNum) {
/*   int ri,mi=0; */
/*   for(ri=0;ri<Nsite;ri++) { */
/*     if(eleNum[ri]>0) { */
/*       eleCfg[ri]=mi; */
/*       eleIdx[mi]=ri; */
/*       mi++; */
/*     } else { */
/*       eleCfg[ri]=-1; */
/*     } */
/*   } */
/*   mi=0; */
/*   for(ri=0;ri<Nsite;ri++) { */
/*     if(eleNum[ri+Nsite]>0) { */
/*       eleCfg[ri+Nsite]=mi; */
/*       eleIdx[mi+Ne]=ri; */
/*       mi++; */
/*     } else { */
/*       eleCfg[ri+Nsite]=-1; */
/*     } */
/*   } */

  return;
}

void ReduceCounter(MPI_Comm comm) {
  #ifdef _mpi_use
  int n=4;
  int recv[n];
  int i;
  int rank,size;
  MPI_Comm_size(comm,&size);
  MPI_Comm_rank(comm,&rank);

  MPI_Allreduce(Counter,recv,n,MPI_INT,MPI_SUM,comm);
  if(rank==0) {
    for(i=0;i<n;i++) Counter[i] = recv[i];
  }
  #endif
  return;
}

