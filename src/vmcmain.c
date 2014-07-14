/*-------------------------------------------------------------
 * Variational Monte Carlo
 * main program
 *-------------------------------------------------------------
 * by Satoshi Morita
 *-------------------------------------------------------------*/
/* #include "fjcoll.h" */
#include "vmcmain.h"

int VMCParaOpt(MPI_Comm comm_parent, MPI_Comm comm_child1, MPI_Comm comm_child2);
int VMCPhysCal(MPI_Comm comm_parent, MPI_Comm comm_child1, MPI_Comm comm_child2);
void outputData();
void printUsageError();
void initMultiDef(char *fileMultiDef, char *fileDefList, char *fileInitPara, 
                  MPI_Comm comm_parent, MPI_Comm *comm_child1);

/*main program*/
int main(int argc, char* argv[])
{
  /* input file name */
  char fileDefList[D_FileNameMax];
  char fileInitPara[D_FileNameMax];
  
  int info=0;

  /* for MPI */
  int rank0=0,size0=1;
  int group1=0,group2=0,rank1=0,rank2=0,size1=1,size2=1;
  MPI_Comm comm0,comm1,comm2;

  /* check the number of arguments */
  if(argc<2) {
    printUsageError();
    exit(EXIT_FAILURE);
  }

  MPI_Init(&argc, &argv);
  NThread = omp_get_max_threads();

  InitTimer();
  StartTimer(0);
  StartTimer(1);
  StartTimer(10);

  initMultiDef(argv[1],fileDefList,fileInitPara,MPI_COMM_WORLD,&comm0);

  MPI_Comm_rank(comm0, &rank0);
  MPI_Comm_size(comm0, &size0);

  StopTimer(10);

  StartTimer(11);
  ReadDefFileNInt(fileDefList, comm0);
  StopTimer(11);
  
  StartTimer(12);
  SetMemoryDef();
  StopTimer(12);
  
  StartTimer(11);
  ReadDefFileIdxPara(fileDefList, comm0);
  StopTimer(11);
  
  StartTimer(12);
  SetMemory();
  StopTimer(12);
  
  /* split MPI coummunicator */
#ifdef _mpi_use
  StartTimer(10);
  group1 = rank0/NSplitSize;
  MPI_Comm_split(comm0,group1,rank0,&comm1);
  MPI_Comm_size(comm1,&size1);
  MPI_Comm_rank(comm1,&rank1);
  group2 = rank1;
  MPI_Comm_split(comm0,group2,rank0,&comm2);
  MPI_Comm_size(comm2,&size2);
  MPI_Comm_rank(comm2,&rank2);

  if(size0%NSplitSize!=0 && rank0==0) {
    fprintf(stderr,"warning: load imbalance. MPI_size0=%d NSplitSize=%d\n",size0,NSplitSize);
  }
  /*   printf("rank=%d group1=%d rank1=%d rank2=%d size1=%d size2=%d\n", */
  /*      rank,group1,rank1,rank2,size1,size2); */
  StopTimer(10);
#endif

  /* initialize Mersenne Twister */
  init_gen_rand(RndSeed+group1);
  /* get the size of work space for LAPACK and PFAPACK */
  LapackLWork = getLWork();

  StartTimer(13);
  /* initialize variational parameters */
  InitParameter(); /* Run parallelly for synchronization of random generator */
  if(rank0==0 && strlen(fileInitPara)>0) ReadInitParameter(fileInitPara);
  SyncModifiedParameter(comm0);
  StopTimer(13);

  /* initialize variables for quantum projection */
  InitQPWeight();
  /* initialize output files */
  if(rank0==0) InitFile(fileDefList, rank0);

  StopTimer(1);

  if(NVMCCalMode==0) {
    StartTimer(2);
    /*-- VMC Parameter Optimization --*/
    VMCParaOpt(comm0, comm1, comm2);
    StopTimer(2);
  } else if(NVMCCalMode==1) {
    StartTimer(2);
    /*-- VMC Physical Quantity Calculation --*/
    VMCPhysCal(comm0, comm1, comm2);
    StopTimer(2);
  } else {
    info=1;
    if(rank0==0) fprintf(stderr,"error: NVMCCalMode must be 0 or 1.\n");
  }

  StopTimer(0);
  if(rank0==0) {
    if(NVMCCalMode==0) {
      OutputTimerParaOpt();
    } else if(NVMCCalMode==1) {
      OutputTimerPhysCal();
    } 
  }

  /* close output files */
  if(rank0==0) CloseFile(rank0);

  FreeMemory();
  FreeMemoryDef();
  MPI_Finalize();

  return info;
}

/*-- VMC Parameter Optimization --*/
int VMCParaOpt(MPI_Comm comm_parent, MPI_Comm comm_child1, MPI_Comm comm_child2) {
  int step;
  int info;
  int rank;
  MPI_Comm_rank(comm_parent, &rank);

  for(step=0;step<NSROptItrStep;step++) {
    if(rank==0) OutputTime(step);
      StartTimer(20);
    UpdateSlaterElm();
      StopTimer(20);
      StartTimer(3);
    VMCMakeSample(comm_child1);
      StopTimer(3);
      StartTimer(4);
    VMCMainCal(comm_child1);
      StopTimer(4);
      StartTimer(21);
    WeightAverageWE(comm_parent);
    WeightAverageSROpt(comm_parent);
    ReduceCounter(comm_child2);
      StopTimer(21);
      StartTimer(22);
    /* output zvo_out and zvo_var */
    if(rank==0) outputData();
      StopTimer(22);
      StartTimer(5);
    info = StochasticOpt(comm_parent);
      StopTimer(5);

    if(info!=0) {
      if(rank==0) fprintf(stderr, "Error: StcOpt info=%d step=%d\n",info,step);
      return info;
    }

      StartTimer(23);
    SyncModifiedParameter(comm_parent);
      StopTimer(23);

    if(step >= NSROptItrStep-NSROptItrSmp) {
      StoreOptData(step-(NSROptItrStep-NSROptItrSmp));
    }
  }

  if(rank==0) OutputTime(NSROptItrStep);

  /* output zqp_opt */
  if(rank==0) OutputOptData();

  return 0;
}

/*-- VMC Physical Quantity Calculation --*/
int VMCPhysCal(MPI_Comm comm_parent, MPI_Comm comm_child1, MPI_Comm comm_child2) {
  int ismp;
  int rank;
  MPI_Comm_rank(comm_parent, &rank);

  StartTimer(20);
  UpdateSlaterElm();
  StopTimer(20);

  for(ismp=0;ismp<NDataQtySmp;ismp++) {
    if(rank==0) OutputTime(ismp);

    InitFilePhysCal(ismp, rank);
    
    StartTimer(3);

    VMCMakeSample(comm_child1);

    StopTimer(3);
    StartTimer(4);

    VMCMainCal(comm_child1);

    StopTimer(4);
    StartTimer(21);

    WeightAverageWE(comm_parent);
    WeightAverageGreenFunc(comm_parent);
    ReduceCounter(comm_child2);

    StopTimer(21);
    StartTimer(22);
    /* output zvo_out and green functions */
    if(rank==0) outputData();
    CloseFilePhysCal(rank);

    StopTimer(22);
    StopTimer(5);
  }

  if(rank==0) OutputTime(NDataQtySmp);

  return 0;
}

void outputData() {
  int i;

  /* zvo_out.dat */
  fprintf(FileOut, "% .8e % .8e % .8e \n", Etot, Etot2, (Etot2 - Etot*Etot)/(Etot*Etot));
  
  /* zvo_var.dat */
  fwrite(Para,sizeof(double),NPara,FileVar);

  if(NVMCCalMode==1) {
    for(i=0;i<NCisAjs;i++) fprintf(FileCisAjs, "% .18e  ", PhysCisAjs[i]);
    fprintf(FileCisAjs, "\n");
    
    for(i=0;i<NCisAjsCktAlt;i++) fprintf(FileCisAjsCktAlt, "% .18e  ", PhysCisAjsCktAlt[i]);
    fprintf(FileCisAjsCktAlt, "\n");
  }

  return;
}

void printUsageError() {
  fprintf(stderr,"Usage: vmc.out MultiDefFile\n");
  return;
}

/* This function splits MPI communicator,
   reads dirName, fileDefList, fileInitPara from fileMultiDef.
   This function sets fileDefList, fileInitPara and changes current working directory
   only in processes with rank1=0 (i.e. group2=0). */
void initMultiDef(char *fileMultiDef, char *fileDefList, char *fileInitPara, 
                  MPI_Comm comm_parent, MPI_Comm *comm_child1) {
  const int lengthTmpString=3*D_FileNameMax;
  char tmpString[lengthTmpString];
  char dirName[D_FileNameMax];
  char *dirNameList;
  FILE *fp;
  int i;
  int nMultiDef=0;

  int rank, size;
  int group1, group2, rank1;
  int div, mod, threshold;
  MPI_Comm comm_child2;

  MPI_Comm_rank(comm_parent, &rank);
  MPI_Comm_size(comm_parent, &size);

  if(rank==0) {
    if( (fp=fopen(fileMultiDef, "r")) != NULL ) {
      if( fgets(tmpString, lengthTmpString, fp) != NULL) {
        sscanf(tmpString, "%d\n", &nMultiDef);
      } else {
        fprintf(stderr,"error: %s is incomplete. (line=1)\n",fileMultiDef);
        MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
      }
    } else {
      fprintf(stderr,"error: %s does not exist.\n",fileMultiDef);
      MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
    }

    /* check MPI size and nMultiDef */
    if(nMultiDef<1) {
      fprintf(stderr,"error: nMultiDef (%d) sets should be larger than 1.\n",nMultiDef);
      MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
    } else if(size<nMultiDef) {
      fprintf(stderr,"error: nMultiDef (%d) sets should be smaller than MPI size (%d).\n",
              nMultiDef,size);
      MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
    } else if(size%nMultiDef!=0) {
      fprintf(stderr,"warning: load imbalance. MPI_size=%d nMultiDef=%d\n",size,nMultiDef);
    }
  }

  MPI_Bcast(&nMultiDef, 1, MPI_INT, 0, comm_parent);

  /* split MPI communicator */
  div = size / nMultiDef;
  mod = size % nMultiDef;
  threshold = (div+1)*mod;
  if(rank < threshold) {
    group1 = rank / (div+1);
  } else {
    group1 = mod + (rank-threshold)/div;
  }
  MPI_Comm_split(comm_parent,group1,rank,comm_child1);
  MPI_Comm_rank((*comm_child1), &rank1);
  group2 = rank1;
  MPI_Comm_split(comm_parent,group2,rank,&comm_child2);

  /* read fileDirList (only at rank=0 process) */
  if(rank==0) {
    dirNameList = (char*)malloc(nMultiDef*lengthTmpString*sizeof(char));
    for(i=0;i<nMultiDef;i++) {
      /* read nMultiDef lines */
      if( fgets(dirNameList + i*lengthTmpString, lengthTmpString, fp) == NULL) {
        fprintf(stderr,"error: %s is incomplete. (line=%d)\n",fileMultiDef,i+2);
        MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
      }
    }
    fclose(fp);
  }    

  if(group2==0) { /* rank1==0 */
    /* scatter and broadcast dirName */
    MPI_Scatter(dirNameList,lengthTmpString,MPI_CHAR,tmpString,lengthTmpString,MPI_CHAR,0,comm_child2);
    i = sscanf(tmpString, "%s %s %s\n", dirName, fileDefList, fileInitPara);
    if(i<2) {
      fprintf(stderr,"error: %s is incomplete. (line:%d)\n",fileMultiDef,group1+2);
      MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
    } else if (i==2) {
      fprintf(stderr,"warning: fileInitPara is not set. (group=%d)\n",group1);
      fileInitPara[0]='\0';
    }

    /* change current working directory */
    if( chdir(dirName) != 0) {
      /* error handle */
      fprintf(stderr,"error: chdir(): %s: ",dirName);
      perror("");
      MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
    }
  }

  MPI_Comm_free(&comm_child2);
  return;
}
