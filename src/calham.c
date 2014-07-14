/*-------------------------------------------------------------
 * Variational Monte Carlo
 * calculate Hamiltonian
 *-------------------------------------------------------------
 * by Satoshi Morita
 *-------------------------------------------------------------*/

double CalculateHamiltonian(const double ip, int *eleIdx, const int *eleCfg,
                             int *eleNum, const int *eleProjCnt);

double CalculateHamiltonian(const double ip, int *eleIdx, const int *eleCfg,
                             int *eleNum, const int *eleProjCnt) {
  double e=0.0;
  int idx;
  int ri,rj,s,rk,rl,t;
  int *myEleIdx, *myEleNum, *myProjCntNew;
  double *myBuffer;
  double myEnergy;

  RequestWorkSpaceThreadInt(Nsize+Nsite2+NProj);
  RequestWorkSpaceThreadDouble(NQPFull+2*Nsize);
  /* GreenFunc1: NQPFull, GreenFunc2: NQPFull+2*Nsize */

#pragma omp parallel default(shared)\
  private(myEleIdx,myEleNum,myProjCntNew,myBuffer,myEnergy)\
  reduction(+:e)
  {
    myEleIdx = GetWorkSpaceThreadInt(Nsize);
    myEleNum = GetWorkSpaceThreadInt(Nsite2);
    myProjCntNew = GetWorkSpaceThreadInt(NProj);
    myBuffer = GetWorkSpaceThreadDouble(NQPFull+2*Nsize);

    #pragma loop noalias
    for(idx=0;idx<Nsize;idx++) myEleIdx[idx] = eleIdx[idx];
    #pragma loop noalias
    for(idx=0;idx<Nsite2;idx++) myEleNum[idx] = eleNum[idx];
    #pragma omp barrier
    
    myEnergy = 0.0;

    #pragma omp master
    {StartTimer(70);}

    /* Coulomb */
    #pragma omp for private(idx,ri,rj) nowait
    for(idx=0;idx<NCoulomb;idx++) {
      ri = Coulomb[idx][0];
      rj = Coulomb[idx][1];
      myEnergy += ParaCoulomb[idx] * eleNum[ri] * eleNum[rj];
    }

    #pragma omp master
    {StopTimer(70);StartTimer(71);}

    /* Transfer */
    #pragma omp for private(idx,ri,rj,s) schedule(dynamic) nowait
    for(idx=0;idx<NTransfer;idx++) {
      ri = Transfer[idx][0];
      rj = Transfer[idx][1];
      s  = Transfer[idx][2];
      
      myEnergy -= ParaTransfer[idx]
        * GreenFunc1(ri,rj,s,ip,myEleIdx,eleCfg,myEleNum,eleProjCnt,myProjCntNew,myBuffer);
      /* Caution: negative sign */
    }

    #pragma omp master
    {StopTimer(71);StartTimer(72);}

    /* Interaction */
    #pragma omp for private(idx,ri,rj,s,rk,rl,t) schedule(dynamic) nowait
    for(idx=0;idx<NInteraction;idx++) {
      ri = Interaction[idx][0];
      rj = Interaction[idx][1];
      s  = Interaction[idx][2];
      rk = Interaction[idx][3];
      rl = Interaction[idx][4];
      t  = Interaction[idx][5];
      
      myEnergy += ParaInteraction[idx]
        * GreenFunc2(ri,rj,rk,rl,s,t,ip,myEleIdx,eleCfg,myEleNum,eleProjCnt,myProjCntNew,myBuffer);
    }

    #pragma omp master
    {StopTimer(72);}

    e += myEnergy;
  }

  ReleaseWorkSpaceThreadInt();
  ReleaseWorkSpaceThreadDouble();
  return e;
}
