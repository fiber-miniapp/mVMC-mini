/*-------------------------------------------------------------
 * Variational Monte Carlo
 * slater elements
 *-------------------------------------------------------------
 * by Satoshi Morita
 *-------------------------------------------------------------*/

void UpdateSlaterElm();
void SlaterElmDiff(double *srOptO, const double ip, int *eleIdx);

void UpdateSlaterElm() {
  int ri,tri,rsi0,rsi1;
  int rj,trj,rsj0,rsj1;
  int qpidx,mpidx,spidx;
  double cs,cc,ss;
  double slt_ij,slt_ji;
  int *xqp;
  double *sltE,*sltE_i0,*sltE_i1;

  #pragma omp parallel for default(shared)                      \
    private(qpidx,mpidx,spidx,                                  \
            xqp,cs,cc,ss,sltE,                                  \
            ri,tri,rsi0,rsi1,sltE_i0,sltE_i1,                   \
            rj,trj,rsj0,rsj1,slt_ij,slt_ji)
  #pragma loop noalias
#pragma ivdep
  for(qpidx=0;qpidx<NQPFull;qpidx++) {
    mpidx = qpidx / NSPGaussLeg;
    spidx = qpidx % NSPGaussLeg;

    xqp = QPTrans[mpidx];
    cs = SPGLCosSin[spidx];
    cc = SPGLCosCos[spidx];
    ss = SPGLSinSin[spidx];
    
    sltE = SlaterElm + qpidx*Nsite2*Nsite2;
    
    for(ri=0;ri<Nsite;ri++) {
      tri = xqp[ri];
      rsi0 = ri;
      rsi1 = ri+Nsite;
      sltE_i0 = sltE + rsi0*Nsite2;
      sltE_i1 = sltE + rsi1*Nsite2;
      
      for(rj=0;rj<Nsite;rj++) {
        trj = xqp[rj];
        rsj0 = rj;
        rsj1 = rj+Nsite;
        
        slt_ij = Slater[ OrbitalIdx[tri][trj] ];
        slt_ji = Slater[ OrbitalIdx[trj][tri] ];
        
        sltE_i0[rsj0] = -(slt_ij - slt_ji)*cs;
        sltE_i0[rsj1] = slt_ij*cc + slt_ji*ss;
        sltE_i1[rsj0] = -slt_ij*ss - slt_ji*cc;
        sltE_i1[rsj1] = (slt_ij - slt_ji)*cs;
      }
    }
  }

  return;
}

void SlaterElmDiff(double *srOptO, const double ip, int *eleIdx) {
  const int nBuf=NSlater*NQPFull;
  const int nsize = Nsize;
  const int ne = Ne;
  const int nQPFull = NQPFull;
  const int nMPTrans = NMPTrans;
  const int nSlater = NSlater;

  const double invIP = 1.0/ip;
  int msi,msj,ri,rj,tri,trj;
  int mpidx,spidx,orbidx,qpidx,i;
  double cs,cc,ss;
  int *xqp;
  double *invM,*invM_i;

  int *orbitalIdx_i;
  int *transOrbIdx;
  int *tOrbIdx,*tOrbIdx_i;
  double *buf, *buffer;
  double tmp;

  RequestWorkSpaceInt(nMPTrans*Nsize*Nsize);
  RequestWorkSpaceDouble(NQPFull*NSlater);

  transOrbIdx = GetWorkSpaceInt(nMPTrans*Nsize*Nsize); /* transOrbIdx[mpidx][msi][msj] */
  buffer = GetWorkSpaceDouble(NQPFull*NSlater);

  for(i=0;i<nBuf;i++) buffer[i]=0.0;

  #pragma omp parallel for default(shared)  \
    private(qpidx,mpidx,msi,msj,xqp,        \
            ri,tri,rj,trj,                  \
            tOrbIdx,tOrbIdx_i,orbitalIdx_i)
  #pragma loop noalias
  for(mpidx=0;mpidx<nMPTrans;mpidx++) {
    xqp = QPTrans[mpidx];
    tOrbIdx = transOrbIdx + mpidx*nsize*nsize;
    for(msi=0;msi<nsize;msi++) {
      ri = eleIdx[msi];
      tri = xqp[ri];
      tOrbIdx_i = tOrbIdx + msi*nsize;
      orbitalIdx_i = OrbitalIdx[tri];
      for(msj=0;msj<nsize;msj++) {
        rj = eleIdx[msj];
        trj = xqp[rj];
        tOrbIdx_i[msj] = orbitalIdx_i[trj];
      }
    }
  }

  #pragma omp parallel for default(shared)        \
    private(qpidx,mpidx,spidx,cs,cc,ss,           \
            tOrbIdx,invM,buf,msi,msj,             \
            tOrbIdx_i,invM_i,orbidx)
  #pragma loop noalias
  for(qpidx=0;qpidx<nQPFull;qpidx++) {
    mpidx = qpidx / NSPGaussLeg;
    spidx = qpidx % NSPGaussLeg;

    cs = PfM[qpidx] * SPGLCosSin[spidx];
    cc = PfM[qpidx] * SPGLCosCos[spidx];
    ss = PfM[qpidx] * SPGLSinSin[spidx];

    tOrbIdx = transOrbIdx + mpidx*nsize*nsize;
    invM = InvM + qpidx*Nsize*Nsize;
    buf = buffer + qpidx*NSlater;

    #pragma loop norecurrence
    for(msi=0;msi<ne;msi++) {
      tOrbIdx_i = tOrbIdx + msi*nsize;
      invM_i = invM + msi*nsize;
      for(msj=0;msj<ne;msj++) {
        /* si=0 sj=0*/
        orbidx = tOrbIdx_i[msj];
        buf[orbidx] += invM_i[msj]*cs;
      }
      for(msj=ne;msj<nsize;msj++) {
        /* si=0 sj=1*/
        orbidx = tOrbIdx_i[msj];
        buf[orbidx] -= invM_i[msj]*cc;
      }
    }
    #pragma loop norecurrence
    for(msi=ne;msi<nsize;msi++) {
      tOrbIdx_i = tOrbIdx + msi*nsize;
      invM_i = invM + msi*nsize;
      for(msj=0;msj<ne;msj++) {
        /* si=1 sj=0*/
        orbidx = tOrbIdx_i[msj];
        buf[orbidx] += invM_i[msj]*ss;
      }
      for(msj=ne;msj<nsize;msj++) {
        /* si=1 sj=1*/
        orbidx = tOrbIdx_i[msj];
        buf[orbidx] -= invM_i[msj]*cs;
      }
    }
  }

  /* store SROptO[] */
  for(orbidx=0;orbidx<nSlater;orbidx++) {
    srOptO[orbidx] = 0.0;
  }
  #pragma loop noalias
  for(qpidx=0;qpidx<nQPFull;qpidx++) {
    tmp = QPFullWeight[qpidx];
    buf = buffer + qpidx*nSlater;
    for(orbidx=0;orbidx<nSlater;orbidx++) {
      srOptO[orbidx] += tmp * buf[orbidx];
    }
  }
  for(orbidx=0;orbidx<nSlater;orbidx++) {
    srOptO[orbidx] *= invIP;
  }

  ReleaseWorkSpaceInt();
  ReleaseWorkSpaceDouble();
  return;
}
