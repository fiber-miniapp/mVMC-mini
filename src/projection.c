/*-------------------------------------------------------------
 * Variational Monte Carlo
 * Gutzwiller-Jastrow Projection and DH correlation factors
 *-------------------------------------------------------------
 * by Satoshi Morita
 *-------------------------------------------------------------*/

inline double LogProjVal(const int *projCnt);
inline double LogProjRatio(const int *projCntNew, const int *projCntOld);
inline double ProjRatio(const int *projCntNew, const int *projCntOld);
void MakeProjCnt(int *projCnt, const int *eleNum);
void UpdateProjCnt(const int ri, const int rj, const int s,
                   int *projCntNew, const int *projCntOld,
                   const int *eleNum);

inline double LogProjVal(const int *projCnt) {
  int idx;
  double z=0;
  for(idx=0;idx<NProj;idx++) {
    z += Proj[idx] * (double)(projCnt[idx]);
  }
  return z;
}

inline double LogProjRatio(const int *projCntNew, const int *projCntOld) {
  int idx;
  double z=0;
  for(idx=0;idx<NProj;idx++) {
    z += Proj[idx] * (double)(projCntNew[idx]-projCntOld[idx]);
  }
  return z;
}

inline double ProjRatio(const int *projCntNew, const int *projCntOld) {
  int idx;
  double z=0;
  for(idx=0;idx<NProj;idx++) {
    z += Proj[idx] * (double)(projCntNew[idx]-projCntOld[idx]);
  }
  return exp(z);
}

void MakeProjCnt(int *projCnt, const int *eleNum) {
  const int *n0=eleNum;
  const int *n1=eleNum+Nsite;
  int idx,offset;
  int ri,rj;
  int xi,xj;
  /* optimization for Kei */
  const int nProj=NProj;
  const int nSite=Nsite;

  /* initialization */
  for(idx=0;idx<nProj;idx++) projCnt[idx] = 0;

  /* Gutzwiller factor */
  if(NGutzwillerIdx>0) {
    for(ri=0;ri<nSite;ri++) {
      projCnt[ GutzwillerIdx[ri] ] += n0[ri]*n1[ri];
    }
  }

  /* Jastrow factor exp(sum {v_ij * (ni-1) * (nj-1)}) */
  if(NJastrowIdx>0) {
    offset = NGutzwillerIdx;
    for(ri=0;ri<nSite;ri++) {
      xi = n0[ri]+n1[ri]-1;
      if(xi==0) continue;
      for(rj=ri+1;rj<nSite;rj++) {
        xj = n0[rj]+n1[rj]-1;
        idx = offset + JastrowIdx[ri][rj];
        projCnt[idx] += xi*xj;
      }
    }
  }

  return;
}

/* An electron with spin s hops from ri to rj. */
void UpdateProjCnt(const int ri, const int rj, const int s,
                   int *projCntNew, const int *projCntOld,
                   const int *eleNum) {
  const int *n0=eleNum;
  const int *n1=eleNum+Nsite;
  int idx,offset;
  int rk;
  /* optimization for Kei */
  const int nProj=NProj;
  const int nSite=Nsite;

  if(projCntNew!=projCntOld) {
    for(idx=0;idx<nProj;idx++) projCntNew[idx] = projCntOld[idx];
  }
  if(ri==rj) return;

  if(NGutzwillerIdx>0){
    idx = GutzwillerIdx[ri];
    projCntNew[idx] -= n0[ri]+n1[ri];
    idx = GutzwillerIdx[rj];
    projCntNew[idx] += n0[rj]*n1[rj];
  }

  if(NJastrowIdx>0){
    offset = NGutzwillerIdx; 
    /* update [ri][rj] */
    if(ri<rj) idx = offset + JastrowIdx[ri][rj];
    else  idx = offset + JastrowIdx[rj][ri];
    projCntNew[idx] += n0[ri]+n1[ri]-n0[rj]-n1[rj]+1;
    /* update [ri][rk] (rk != ri, rj) */
    for(rk=0;rk<nSite;rk++) {
      if(rk==rj) continue;
      if(rk==ri) continue;
      if(rk>ri) idx = offset + JastrowIdx[ri][rk];
      else idx = offset + JastrowIdx[rk][ri]; 
      projCntNew[idx] -= n0[rk]+n1[rk]-1;
    }
    /* update [rj][rk] (rk != ri, rj) */
    for(rk=0;rk<nSite;rk++) {
      if(rk==ri) continue;
      if(rk==rj) continue;
      if(rk>rj) idx = offset + JastrowIdx[rj][rk];
      else idx = offset + JastrowIdx[rk][rj]; 
      projCntNew[idx] += n0[rk]+n1[rk]-1;
    }
  }
  
  return;

}
