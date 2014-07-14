/*-------------------------------------------------------------
 * Variational Monte Carlo
 * timer program
 * "-lrt" option is needed for clock_gettime().
 *-------------------------------------------------------------
 * by Satoshi Morita
 *-------------------------------------------------------------*/

void OutputTime(int step);
void InitTimer();
inline void StartTimer(int n);
inline void StopTimer(int n);
void OutputTimerParaOpt();
void OutputTimerPhysCal();

void OutputTime(int step) {
  time_t tx;
  double pHop,pEx;

  tx = time(NULL);
  if(step==0) {
    fprintf(FileTime, "%05d  acc_hop acc_ex  n_hop    n_ex     : %s", step, ctime(&tx));
  } else {
    pHop = (Counter[0] == 0) ? 0.0 : (double)Counter[1] / (double)Counter[0];
    pEx  = (Counter[2] == 0) ? 0.0 : (double)Counter[3] / (double)Counter[2];
    fprintf(FileTime, "%05d  %.5lf %.5lf %-8d %-8d : %s", step, pHop,pEx,
            Counter[0], Counter[2], ctime(&tx));
  }
}

void InitTimer() {
  int i;
  for(i=0;i<NTimer;i++) Timer[i]=0.0;
  for(i=0;i<NTimer;i++) TimerStart[i]=0.0;
  return;
}

void StartTimer(int n) {
#ifdef _mpi_use
  TimerStart[n]=MPI_Wtime();
#else
  struct timespec ts;
  clock_gettime(CLOCK_REALTIME,&ts);
  TimerStart[n]=ts.tv_sec + ts.tv_nsec*1.0e-9;
#endif
  return;
}

void StopTimer(int n) {
#ifdef _mpi_use
  Timer[n] += MPI_Wtime() - TimerStart[n];
#else
  struct timespec ts;
  clock_gettime(CLOCK_REALTIME,&ts);
  Timer[n] += ts.tv_sec + ts.tv_nsec*1.0e-9 - TimerStart[n];
#endif
  return;
}

void OutputTimerParaOpt() {
  char fileName[D_FileNameMax];
  FILE *fp;
  sprintf(fileName, "%s_HitachiTimer.dat", CDataFileHead); 
  fp = fopen(fileName, "w");

  fprintf(fp,"All                         [0] %12.5lf\n",Timer[0]);
  fprintf(fp,"Initialization              [1] %12.5lf\n",Timer[1]);
  fprintf(fp,"  read options             [10] %12.5lf\n",Timer[10]);
  fprintf(fp,"  ReadDefFile              [11] %12.5lf\n",Timer[11]);
  fprintf(fp,"  SetMemory                [12] %12.5lf\n",Timer[12]);
  fprintf(fp,"  InitParameter            [13] %12.5lf\n",Timer[13]);
  fprintf(fp,"VMCParaOpt                  [2] %12.5lf\n",Timer[2]);
  fprintf(fp,"  VMCMakeSample             [3] %12.5lf\n",Timer[3]);
  fprintf(fp,"    makeInitialSample      [30] %12.5lf\n",Timer[30]);
  fprintf(fp,"    make candidate         [31] %12.5lf\n",Timer[31]);
  fprintf(fp,"    hopping update         [32] %12.5lf\n",Timer[32]);
  fprintf(fp,"      UpdateProjCnt        [60] %12.5lf\n",Timer[60]);
  fprintf(fp,"      CalculateNewPfM2     [61] %12.5lf\n",Timer[61]);
  fprintf(fp,"      CalculateLogIP       [62] %12.5lf\n",Timer[62]);
  fprintf(fp,"      UpdateMAll           [63] %12.5lf\n",Timer[63]);
  fprintf(fp,"    exchange update        [33] %12.5lf\n",Timer[33]);
  fprintf(fp,"      UpdateProjCnt        [65] %12.5lf\n",Timer[65]);
  fprintf(fp,"      CalculateNewPfMTwo2  [66] %12.5lf\n",Timer[66]);
  fprintf(fp,"      CalculateLogIP       [67] %12.5lf\n",Timer[67]);
  fprintf(fp,"      UpdateMAllTwo        [68] %12.5lf\n",Timer[68]);
  fprintf(fp,"    recal PfM and InvM     [34] %12.5lf\n",Timer[34]);
  fprintf(fp,"    save electron config   [35] %12.5lf\n",Timer[35]);
  fprintf(fp,"  VMCMainCal                [4] %12.5lf\n",Timer[4]);
  fprintf(fp,"    CalculateMAll          [40] %12.5lf\n",Timer[40]);
  fprintf(fp,"    LocEnergyCal           [41] %12.5lf\n",Timer[41]);
  fprintf(fp,"      CalHamiltonian0      [70] %12.5lf\n",Timer[70]);
  fprintf(fp,"      CalHamiltonian1      [71] %12.5lf\n",Timer[71]);
  fprintf(fp,"      CalHamiltonian2      [72] %12.5lf\n",Timer[72]);
  fprintf(fp,"    ReturnSlaterElmDiff    [42] %12.5lf\n",Timer[42]);
  fprintf(fp,"    calculate OO and HO    [43] %12.5lf\n",Timer[43]);
  fprintf(fp,"  StochasticOpt             [5] %12.5lf\n",Timer[5]);
  fprintf(fp,"    preprocess             [50] %12.5lf\n",Timer[50]);
  fprintf(fp,"    stcOptMain             [51] %12.5lf\n",Timer[51]);
  fprintf(fp,"      initBLACS            [55] %12.5lf\n",Timer[55]);
  fprintf(fp,"      calculate S and g    [56] %12.5lf\n",Timer[56]);
  fprintf(fp,"      DPOSV                [57] %12.5lf\n",Timer[57]);
  fprintf(fp,"      gatherParaChange     [58] %12.5lf\n",Timer[58]);
  fprintf(fp,"    postprocess            [52] %12.5lf\n",Timer[52]);
  fprintf(fp,"  UpdateSlaterElm          [20] %12.5lf\n",Timer[20]);
  fprintf(fp,"  WeightAverage            [21] %12.5lf\n",Timer[21]);
  fprintf(fp,"  outputData               [22] %12.5lf\n",Timer[22]);
  fprintf(fp,"  SyncModifiedParameter    [23] %12.5lf\n",Timer[23]);

  fclose(fp);
}

void OutputTimerPhysCal() {
  char fileName[D_FileNameMax];
  FILE *fp;
  sprintf(fileName, "%s_HitachiTimer.dat", CDataFileHead); 
  fp = fopen(fileName, "w");

  fprintf(fp,"All                         [0] %12.5lf\n",Timer[0]);
  fprintf(fp,"Initialization              [1] %12.5lf\n",Timer[1]);
  fprintf(fp,"  read options             [10] %12.5lf\n",Timer[10]);
  fprintf(fp,"  ReadDefFile              [11] %12.5lf\n",Timer[11]);
  fprintf(fp,"  SetMemory                [12] %12.5lf\n",Timer[12]);
  fprintf(fp,"  InitParameter            [13] %12.5lf\n",Timer[13]);
  fprintf(fp,"VMCPhysCal                  [2] %12.5lf\n",Timer[2]);
  fprintf(fp,"  VMCMakeSample             [3] %12.5lf\n",Timer[3]);
  fprintf(fp,"    makeInitialSample      [30] %12.5lf\n",Timer[30]);
  fprintf(fp,"    make candidate         [31] %12.5lf\n",Timer[31]);
  fprintf(fp,"    hopping update         [32] %12.5lf\n",Timer[32]);
  fprintf(fp,"      UpdateProjCnt        [60] %12.5lf\n",Timer[60]);
  fprintf(fp,"      CalculateNewPfM2     [61] %12.5lf\n",Timer[61]);
  fprintf(fp,"      CalculateLogIP       [62] %12.5lf\n",Timer[62]);
  fprintf(fp,"      UpdateMAll           [63] %12.5lf\n",Timer[63]);
  fprintf(fp,"    exchange update        [33] %12.5lf\n",Timer[33]);
  fprintf(fp,"      UpdateProjCnt        [65] %12.5lf\n",Timer[65]);
  fprintf(fp,"      CalculateNewPfMTwo2  [66] %12.5lf\n",Timer[66]);
  fprintf(fp,"      CalculateLogIP       [67] %12.5lf\n",Timer[67]);
  fprintf(fp,"      UpdateMAllTwo        [68] %12.5lf\n",Timer[68]);
  fprintf(fp,"    recal PfM and InvM     [34] %12.5lf\n",Timer[34]);
  fprintf(fp,"    save electron config   [35] %12.5lf\n",Timer[35]);
  fprintf(fp,"  VMCMainCal                [4] %12.5lf\n",Timer[4]);
  fprintf(fp,"    CalculateMAll          [40] %12.5lf\n",Timer[40]);
  fprintf(fp,"    LocEnergyCal           [41] %12.5lf\n",Timer[41]);
  fprintf(fp,"      CalHamiltonian0      [70] %12.5lf\n",Timer[70]);
  fprintf(fp,"      CalHamiltonian1      [71] %12.5lf\n",Timer[71]);
  fprintf(fp,"      CalHamiltonian2      [72] %12.5lf\n",Timer[72]);
  fprintf(fp,"    CalculateGreenFunc     [42] %12.5lf\n",Timer[42]);
  fprintf(fp,"      GreenFunc1           [50] %12.5lf\n",Timer[50]);
  fprintf(fp,"      GreenFunc2           [51] %12.5lf\n",Timer[51]);
  fprintf(fp,"  UpdateSlaterElm          [20] %12.5lf\n",Timer[20]);
  fprintf(fp,"  WeightAverage            [21] %12.5lf\n",Timer[21]);
  fprintf(fp,"  outputData               [22] %12.5lf\n",Timer[22]);

  fclose(fp);
}
