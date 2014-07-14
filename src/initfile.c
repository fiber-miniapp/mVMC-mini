/*-------------------------------------------------------------
 * Variational Monte Carlo
 * initialization of files
 *-------------------------------------------------------------
 * by Satoshi Morita and Ryui Kaneko
 *-------------------------------------------------------------*/

void InitFile(char *xNameListFile, int rank);
void InitFilePhysCal(int i, int rank);
void CloseFile(int rank);
void CloseFilePhysCal(int rank);
void writeConfig(char *xNameFile, char *fileName);
int fileCopyAdd(char *inputFileName, FILE *outputFile);

void InitFile(char *xNameListFile, int rank) {
  char fileName[D_FileNameMax];

  if(rank!=0) return;

  sprintf(fileName, "%s_cfg_%03d.dat", CDataFileHead, NDataIdxStart);
  writeConfig(xNameListFile, fileName);

  sprintf(fileName, "%s_time_%03d.dat", CDataFileHead, NDataIdxStart);
  FileTime = fopen(fileName, "w");

  if(NVMCCalMode==0) {
    sprintf(fileName, "%s_SRinfo.dat", CDataFileHead);
    FileSRinfo = fopen(fileName, "w");
    fprintf(FileSRinfo,
            "#Npara Msize optCut diagCut sDiagMax  sDiagMin    absRmax       imax\n");

    sprintf(fileName, "%s_out_%03d.dat", CDataFileHead, NDataIdxStart);
    FileOut = fopen(fileName, "w");

    sprintf(fileName, "%s_varbin_%03d.dat", CDataFileHead, NDataIdxStart);
    FileVar = fopen(fileName, "wb");
    fwrite(&NPara,sizeof(int),1,FileVar);
    fwrite(&NSROptItrStep,sizeof(int),1,FileVar);
  }

  return;
}

void InitFilePhysCal(int i, int rank) {
  char fileName[D_FileNameMax];
  int idx = i+NDataIdxStart;
  int one = 1;

  if(rank!=0) return;

  sprintf(fileName, "%s_out_%03d.dat", CDataFileHead, idx);
  FileOut = fopen(fileName, "w");
 
  sprintf(fileName, "%s_varbin_%03d.dat", CDataFileHead, idx);
  FileVar = fopen(fileName, "wb");
  fwrite(&NPara,sizeof(int),1,FileVar);
  fwrite(&one,sizeof(int),1,FileVar);

  /* Green function */
  sprintf(fileName, "%s_cisajs_%03d.dat", CDataFileHead, idx);
  FileCisAjs = fopen(fileName, "w");

  sprintf(fileName, "%s_cisajscktalt_%03d.dat", CDataFileHead, idx);
  FileCisAjsCktAlt = fopen(fileName, "w");

  return;
}

void CloseFile(int rank) {
  if(rank!=0) return;

  fclose(FileTime);

  if(NVMCCalMode==0) {
    fclose(FileSRinfo);
    fclose(FileOut);
    fclose(FileVar);
  }

  return;
}

void CloseFilePhysCal(int rank) {
  if(rank!=0) return;

  fclose(FileOut);
  fclose(FileVar);
  fclose(FileCisAjs);
  fclose(FileCisAjsCktAlt);

  return;
}

void writeConfig(char *xnamefile, char *fileName) {
  FILE *ofp,*fplist;
  char defname[D_FileNameMax];

  /* clear configfile */
  ofp = fopen(fileName, "w");

  fileCopyAdd(xnamefile,ofp);

  fplist = fopen(xnamefile,"r");
  if(fplist==NULL){
    fprintf(stderr, "Error writecfg.c");
  } else {
    while(fscanf(fplist, "%s\n", defname)!=EOF){
      fileCopyAdd(defname,ofp);
    }
    fclose(fplist);
  }

  fclose(ofp);
  return;
}

int fileCopyAdd(char *inputfileName, FILE *ofp){
  int i;
  FILE *ifp;

  ifp = fopen(inputfileName, "r");
  if(ifp == NULL) return 1;

  fprintf(ofp,"################\n#%s\n################\n",inputfileName);
  while((i=getc(ifp))!=EOF) putc(i,ofp);
  fclose(ifp);

  return 0;
}
