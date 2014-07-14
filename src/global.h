/*-------------------------------------------------------------
 * Variational Monte Carlo
 * global variables
 *-------------------------------------------------------------
 * by Satoshi Morita
 *-------------------------------------------------------------*/

#ifndef _INCLUDE_GLOBAL
#define _INCLUDE_GLOBAL

#define D_FileNameMax 256

/***** definition *****/
char CDataFileHead[D_FileNameMax]; /* prefix of output files */
char CParaFileHead[D_FileNameMax]; /* prefix for optimized variational parameters */

int NVMCCalMode; /* calculation mode
                    0: optimization of variational paraneters,
                    1: calculation of expectation values */

int NDataIdxStart; /* starting value of the file index */
int NDataQtySmp; /* the number of output files */

int Nsite; /* the number of sites */
int Ne;    /* the number of electrons with up spin */
int Nsize; /* the number of electrons = 2*Ne */
int Nsite2; /* 2*Nsite */

int NSPGaussLeg; /* the number of points for the Gauss-Legendre quadrature */
int NSPStot; /* S of Spin projection */
int NMPTrans; /* the number of quantum projection for translation and point group symmetry */
int NQPFull; /* the total number of quantum projection = NSPGaussLeg*NMPTrans*NQPTransOpt */

int NSROptItrStep; /* the number of SR method steps */
int NSROptItrSmp; /* the number of SR method steps for calculation of average value */
int NSROptFixSmp; /* the number of SR method steps with fixed samples (1 is recommended) */

double DSROptRedCut; /* SR stabilizing factor for truncation of redundant directions */
double DSROptStaDel; /* SR stabiliaing factor for diagonal element modification */
double DSROptStepDt; /* step width of the SR method */

int NVMCWarmUp; /* Monte Carlo steps for warming up */
int NVMCIniterval; /* sampling interval [MCS] */ 
int NVMCSample; /* the number of samples */
int NExUpdatePath; /* update by exchange hopping  0: off, 1: on */

int RndSeed; /* seed for pseudorandom number generator */
int NSplitSize; /* the number of inner MPI processes */
 
/* total length of def array */
int NTotalDefInt, NTotalDefDouble;

/* zlocspin.def */
int NLocSpn; /* the number of local spin */
int *LocSpn; /* [Nsite] */
/* local spin flag  0: local spin, 1: itinerant electron */

/* for Hamiltonian */
int NTransfer;
int **Transfer; /* [NTransfer][3] */
double *ParaTransfer;

int NCoulomb;
int **Coulomb; /* [NCoulomb][2] */
double *ParaCoulomb;

int NInteraction;
int **Interaction; /* [NInteraction][6] */
double *ParaInteraction;

/* for variational parameters */
int NGutzwillerIdx, *GutzwillerIdx; /* [Nsite] */
int NJastrowIdx, **JastrowIdx; /* [Nsite][Nsite] */
int NOrbitalIdx, **OrbitalIdx; /* [Nsite][Nsite] */

/* zqptransidx.def */
int NQPTrans, **QPTrans; /* [NQPTrans][Nsite] */
double *ParaQPTrans;

/* for Green functions */
int NCisAjs,         **CisAjsIdx;         /* [NCisAjs][3] */
int NCisAjsCktAlt, **CisAjsCktAltIdx; /* [NCisAjsCktAlt][6] */

/* Optimization flag */
int *OptFlag; /* [NPara]  1: optimized, 0 or 2: fixed */

/***** Variational Parameters *****/
int NPara;    /* the total number of variational prameters = NSlater + NProj */
int NProj;    /* the number of correlation factor */
int NSlater;  /* the number of pair orbital (f_ij) = NOrbitalIdx */
double *Para;   /* variatonal parameters */
double *Proj;   /* correlation factor (Proj  =Para) */
double *Slater; /* pair orbital       (Slater=Para+NProj) */

/***** Electron Configuration ******/
int *EleIdx; /* EleIdx[sample][mi+si*Ne] */
int *EleCfg; /* EleCfg[sample][ri+si*Nsite] */
int *EleNum; /* EleIdx[sample][ri+si*Nsite] */
int *EleProjCnt; /* EleProjCnt[sample][proj] */
double *logSqPfFullSlater; /* logSqPfFullSlater[sample] */

int *TmpEleIdx;
int *TmpEleCfg;
int *TmpEleNum;
int *TmpEleProjCnt;

int *BurnEleIdx;
int *BurnEleCfg;
int *BurnEleNum;
int *BurnEleProjCnt;
int BurnFlag=0; /* 0: off, 1: on */

/***** Slater Elements ******/
double *SlaterElm; /* SlaterElm[QPidx][ri+si*Nsite][rj+sj*Nsite] */

double *InvM; /* InvM[QPidx][mi+si*Ne][mj+sj*Ne] */
double *PfM; /* PfM[QPidx] */

/***** Quantum Projection *****/
double *QPFullWeight; /* QPFullWeight[NQPFull] */
double *SPGLCos, *SPGLSin; /* [NSPGaussLeg]  cos(beta/2) and sin(beta/2) */
double *SPGLCosSin, *SPGLCosCos, *SPGLSinSin; /* [NSPGaussLeg] */

/***** Stocastic Reconfiguration *****/
int SROptSize; /* 1+NPara */
double *SROptOO; /* [SROptSize*SROptSize] <O^\dagger O> */
double *SROptHO; /* [SROptSize]            < HO > */
double *SROptO;  /* [SROptSize] calculation buffar */

double *SROptData; /* [2+NPara] storage for energy and variational parameters */

/***** Physical Quantity *****/
double Wc; /* Weight for correlation sampling = <psi|x> */
double Etot; /* <H> */
double Etot2; /* <H^2> */

double *PhysCisAjs; /* [NCisAjs] */
double *PhysCisAjsCktAlt; /* [NCisAjsCktAlt] */

/***** Output File *****/
/* FILE *FileCfg; */
FILE *FileOut;
FILE *FileVar;
FILE *FileTime;
FILE *FileSRinfo; /* zvo_SRinfo.dat */
FILE *FileCisAjs;
FILE *FileCisAjsCktAlt;

/* FILE *FileTimerList; */
/* FILE *FileOpt;    /\* zqp_opt *\/ */

/***** HitachiTimer *****/
const int NTimer=100;
double Timer[100], TimerStart[100];

/***** openMP *****/
int NThread;

/***** for DGETRI and DSKPFA in CalculateMAll *****/
int LapackLWork;

/***** counter for vmcMake *****/
int Counter[4] = {0,0,0,0};
/* 0: hopping, 1: hopping accept, 2: exchange try, 3: exchange accept */

#endif /*  _INCLUDE_GLOBAL */
