/*-------------------------------------------------------------
 * Variational Monte Carlo
 * Allocate and free memory for global array
 *-------------------------------------------------------------
 * by Satoshi Morita and Ryui Kaneko
 *-------------------------------------------------------------*/

void SetMemoryDef();
void FreeMemoryDef();
void SetMemory();
void FreeMemory();

void SetMemoryDef() {  
  int i;
  int *pInt;
  double *pDouble;

  /* Int */
  LocSpn = (int*)malloc(sizeof(int)*NTotalDefInt);
  pInt = LocSpn + Nsite;

  Transfer = (int**)malloc(sizeof(int*)*NTransfer);
  for(i=0;i<NTransfer;i++) {
    Transfer[i] = pInt;
    pInt += 3;
  }

  Coulomb = (int**)malloc(sizeof(int*)*NCoulomb);
  for(i=0;i<NCoulomb;i++) {
    Coulomb[i] = pInt;
    pInt += 2;
  }

  Interaction = (int**)malloc(sizeof(int*)*NInteraction);
  for(i=0;i<NInteraction;i++) {
    Interaction[i] = pInt;
    pInt += 6;
  }

  GutzwillerIdx = pInt;
  pInt += Nsite;

  JastrowIdx = (int**)malloc(sizeof(int*)*Nsite);
  for(i=0;i<Nsite;i++) {
    JastrowIdx[i] = pInt;
    pInt += Nsite;
  }

  OrbitalIdx = (int**)malloc(sizeof(int*)*Nsite);
  for(i=0;i<Nsite;i++) {
    OrbitalIdx[i] = pInt;
    pInt += Nsite;
  }

  QPTrans = (int**)malloc(sizeof(int*)*NQPTrans);
  for(i=0;i<NQPTrans;i++) {
    QPTrans[i] = pInt;
    pInt += Nsite;
  }

  CisAjsIdx = (int**)malloc(sizeof(int*)*NCisAjs);
  for(i=0;i<NCisAjs;i++) {
    CisAjsIdx[i] = pInt;
    pInt += 3;
  }

  CisAjsCktAltIdx = (int**)malloc(sizeof(int*)*NCisAjsCktAlt);
  for(i=0;i<NCisAjsCktAlt;i++) {
    CisAjsCktAltIdx[i] = pInt;
    pInt += 6;
  }

  OptFlag = pInt;

  /* Double */
  ParaTransfer = (double*)malloc(sizeof(double)*(NTotalDefDouble));
  pDouble = ParaTransfer + NTransfer;

  ParaCoulomb = pDouble;
  pDouble += NCoulomb;

  ParaInteraction = pDouble;
  pDouble +=  NInteraction;
  
  ParaQPTrans = pDouble;
  pDouble +=  NQPTrans;

  return;
}

void FreeMemoryDef() {
  free(ParaTransfer);

  free(CisAjsCktAltIdx);
  free(CisAjsIdx);
  free(QPTrans);
  free(OrbitalIdx);
  free(JastrowIdx);
  free(Interaction);
  free(Coulomb);
  free(Transfer);
  free(LocSpn);

  return;
}

void SetMemory() {

  /***** Variational Parameters *****/
  Para = (double*)malloc(sizeof(double)*(NPara));
  Proj = Para;
  Slater = Para + NProj;

  /***** Electron Configuration ******/
  EleIdx = (int*)malloc(sizeof(int)*( NVMCSample*2*Ne ));
  EleCfg = (int*)malloc(sizeof(int)*( NVMCSample*2*Nsite ));
  EleNum = (int*)malloc(sizeof(int)*( NVMCSample*2*Nsite ));
  EleProjCnt = (int*)malloc(sizeof(int)*( NVMCSample*NProj ));
  logSqPfFullSlater = (double*)malloc(sizeof(double)*(NVMCSample));

  TmpEleIdx = (int*)malloc(sizeof(int)*(2*Ne+2*Nsite+2*Nsite+NProj));
  TmpEleCfg = TmpEleIdx + 2*Ne;
  TmpEleNum = TmpEleCfg + 2*Nsite;
  TmpEleProjCnt = TmpEleNum + 2*Nsite;

  BurnEleIdx = (int*)malloc(sizeof(int)*(2*Ne+2*Nsite+2*Nsite+NProj));
  BurnEleCfg = BurnEleIdx + 2*Ne;
  BurnEleNum = BurnEleCfg + 2*Nsite;
  BurnEleProjCnt = BurnEleNum + 2*Nsite;

  /***** Slater Elements ******/
  SlaterElm = (double*)malloc( sizeof(double)*(NQPFull*(2*Nsite)*(2*Nsite)) );

  InvM = (double*)malloc( sizeof(double)*(NQPFull*(Nsize*Nsize+1)) );
  PfM = InvM + NQPFull*Nsize*Nsize;

  /***** Quantum Projection *****/
  QPFullWeight = (double*)malloc(sizeof(double)*(NQPFull+5*NSPGaussLeg));
  SPGLCos    = QPFullWeight + NQPFull;
  SPGLSin    = SPGLCos + NSPGaussLeg;
  SPGLCosSin = SPGLCos + 2*NSPGaussLeg;
  SPGLCosCos = SPGLCos + 3*NSPGaussLeg;
  SPGLSinSin = SPGLCos + 4*NSPGaussLeg;

  /***** Stocastic Reconfiguration *****/
  if(NVMCCalMode==0){
    SROptOO = (double*)malloc( sizeof(double)*(SROptSize*(SROptSize+2)) );
    SROptHO = SROptOO + SROptSize*SROptSize;
    SROptO  = SROptHO + SROptSize;

    SROptData = (double*)malloc( sizeof(double)*(NSROptItrSmp*(2+NPara)) );
  }

  /***** Physical Quantity *****/
  if(NVMCCalMode==1){
    PhysCisAjs  = (double*)malloc(sizeof(double)*(NCisAjs+NCisAjsCktAlt));
    PhysCisAjsCktAlt = PhysCisAjs + NCisAjs;
  }

  initializeWorkSpaceAll();
  return;
}

void FreeMemory() {
  FreeWorkSpaceAll();

  if(NVMCCalMode==1){
    free(PhysCisAjs);
  }

  if(NVMCCalMode==0){
    free(SROptData);
    free(SROptOO);
  }

  free(QPFullWeight);

  free(InvM);
  free(SlaterElm);

  free(BurnEleIdx);
  free(TmpEleIdx);
  free(logSqPfFullSlater);
  free(EleProjCnt);
  free(EleIdx);
  free(EleCfg);

  free(Para);

  return;
}
