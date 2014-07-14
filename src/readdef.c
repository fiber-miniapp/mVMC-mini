/*-------------------------------------------------------------
 * Variational Monte Carlo
 * Read Definition Files
 *-------------------------------------------------------------
 * by Satoshi Morita
 *-------------------------------------------------------------*/

int ReadDefFileError(char *defname);
int ReadDefFileNInt(char *xNameListFile, MPI_Comm comm);
int ReadDefFileIdxPara(char *xNameListFile, MPI_Comm comm);

int ReadDefFileError(char *defname){
  fprintf(stderr, "error: %s (Broken file or Not exist)\n", defname);
  return 1;
}

int ReadDefFileNInt(char *xNameListFile, MPI_Comm comm){
  FILE *fp, *fplist;
  char defname[D_FileNameMax];
  char ctmp[D_FileNameMax];
  int itmp;

  int rank, info=0;
  const int nBufInt=36;
  const int nBufDouble=3;
  const int nBufChar=D_FileNameMax;
  int bufInt[nBufInt];
  double bufDouble[nBufDouble];

  MPI_Comm_rank(comm, &rank);

  if(rank==0) {
    fplist = fopen(xNameListFile, "r");
    if(fplist!=NULL) {
      /* zmodpara.def */
      if(info==0) {
        if(fscanf(fplist, "%s\n", defname)!=EOF) {
          fp = fopen(defname, "r");
          if(fp!=NULL) {
            fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
            fscanf(fp,"%s %d\n", ctmp, &itmp);
            fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
            fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
            fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
            fscanf(fp,"%s %s\n", ctmp, CDataFileHead);
            fscanf(fp,"%s %s\n", ctmp, CParaFileHead);
            fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
            fscanf(fp,"%s %d\n", ctmp, &(bufInt[ 0])); /* NVMCCalMode */
            fscanf(fp,"%s %d\n", ctmp, &(bufInt[ 1])); /* NLanczosMode */
            fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
            fscanf(fp,"%s %d\n", ctmp, &(bufInt[ 2])); /* NDataIdxStart */
            fscanf(fp,"%s %d\n", ctmp, &(bufInt[ 3])); /* NDataQtySmp */
            fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
            fscanf(fp,"%s %d\n", ctmp, &(bufInt[ 4])); /* Nsite */
            fscanf(fp,"%s %d\n", ctmp, &(bufInt[ 5])); /* Ne */
            fscanf(fp,"%s %d\n", ctmp, &(bufInt[ 6])); /* NSPGaussLeg */
            fscanf(fp,"%s %d\n", ctmp, &(bufInt[ 7])); /* NSPStot */
            fscanf(fp,"%s %d\n", ctmp, &(bufInt[ 8])); /* NMPTrans */
            fscanf(fp,"%s %d\n", ctmp, &(bufInt[ 9])); /* NSROptItrStep */
            fscanf(fp,"%s %d\n", ctmp, &(bufInt[10])); /* NSROptItrSmp */
            fscanf(fp,"%s %d\n", ctmp, &(bufInt[11])); /* NSROptFixSmp */
            fscanf(fp,"%s %lf\n",ctmp, &(bufDouble[0])); /* DSROptRedCut */
            fscanf(fp,"%s %lf\n",ctmp, &(bufDouble[1])); /* DSROptStaDel */
            fscanf(fp,"%s %lf\n",ctmp, &(bufDouble[2])); /* DSROptStepDt */
            fscanf(fp,"%s %d\n", ctmp, &(bufInt[12])); /* NVMCWarmUp */
            fscanf(fp,"%s %d\n", ctmp, &(bufInt[13])); /* NVMCIniterval */
            fscanf(fp,"%s %d\n", ctmp, &(bufInt[14])); /* NVMCSample */
            fscanf(fp,"%s %d\n", ctmp, &(bufInt[15])); /* NExUpdatePath */
            fscanf(fp,"%s %d\n", ctmp, &(bufInt[16])); /* RndSeed */
            if(bufInt[16]<0) {
              bufInt[16] = (int)time(NULL);
              fprintf(stderr, "remark: Seed = %d\n", bufInt[16]);
            }
            fscanf(fp,"%s %d\n", ctmp, &(bufInt[17])); /* NSplitSize */
            fclose(fp);
          } else { info = ReadDefFileError(defname); }
        } else { info = ReadDefFileError(xNameListFile); }
      }
      /*locspn.def----------------------------------------*/
      if(info==0) {
        if(fscanf(fplist, "%s\n", defname)!=EOF) {
          fp = fopen(defname, "r");
          if(fp!=NULL) {
            fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
            fscanf(fp,"%s %d\n", ctmp, &(bufInt[18])); /* NLocSpn */
            fclose(fp);
          } else { info =  ReadDefFileError(defname); }
        } else { info = ReadDefFileError(xNameListFile); }
      }
      /*transfer.def--------------------------------------*/
      if(info==0) {
        if(fscanf(fplist, "%s\n", defname)!=EOF) {
          fp = fopen(defname, "r");
          if(fp!=NULL) {
            fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
            fscanf(fp,"%s %d\n", ctmp, &(bufInt[19])); /* NTransfer */
            fclose(fp);
          } else { info= ReadDefFileError(defname);}
        } else {info = ReadDefFileError(xNameListFile);}
      }
      /*coulomb.def----------------------------------*/
      if(info==0) {
        if(fscanf(fplist, "%s\n", defname)!=EOF) {
          fp = fopen(defname, "r");
          if(fp!=NULL) {
            fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
            fscanf(fp,"%s %d\n", ctmp, &(bufInt[20])); /* NCoulomb */
            fclose(fp);
          } else { info = ReadDefFileError(defname);}
        } else { info = ReadDefFileError(xNameListFile);}
      }
      /*interaction.def----------------------------------*/
      if(info==0) {
        if(fscanf(fplist, "%s\n", defname)!=EOF) {
          fp = fopen(defname, "r");
          if(fp!=NULL) {
            fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
            fscanf(fp,"%s %d\n", ctmp, &(bufInt[21])); /* NInteraction */
            fclose(fp);
          } else { info = ReadDefFileError(defname);}
        } else { info = ReadDefFileError(xNameListFile);}
      }
      /*gutzwilleridx.def---------------------------------*/
      if(info==0) {
        if(fscanf(fplist, "%s\n", defname)!=EOF) {
          fp = fopen(defname, "r");
          if(fp!=NULL) {
            fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
            fscanf(fp,"%s %d\n", ctmp, &(bufInt[22])); /* NGutzwillerIdx */
            fclose(fp);
          } else { info = ReadDefFileError(defname);}
        } else { info = ReadDefFileError(xNameListFile);}
      }
      /*jastrowidx.def------------------------------------*/
      if(info==0) {
        if(fscanf(fplist, "%s\n", defname)!=EOF) {
          fp = fopen(defname, "r");
          if(fp!=NULL) {
            fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
            fscanf(fp,"%s %d\n", ctmp, &(bufInt[23])); /* NJastrowIdx */
            fclose(fp);
          } else { info = ReadDefFileError(defname);}
        } else { info = ReadDefFileError(xNameListFile);}
      }
      /*orbitalidx.def------------------------------------*/
      if(info==0) {
        if(fscanf(fplist, "%s\n", defname)!=EOF) {
          fp = fopen(defname, "r");
          if(fp!=NULL) {
            fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
            fscanf(fp,"%s %d\n", ctmp, &(bufInt[24])); /* NOrbitalIdx */
            fclose(fp);
          } else { info = ReadDefFileError(defname);}
        } else { info = ReadDefFileError(xNameListFile);}
      }
      /*qptransidx.def------------------------------------*/
      if(info==0) {
        if(fscanf(fplist, "%s\n", defname)!=EOF) {
          fp = fopen(defname, "r");
          if(fp!=NULL) {
            fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
            fscanf(fp,"%s %d\n", ctmp, &(bufInt[25])); /* NQPTrans */
            fclose(fp);
          } else { info = ReadDefFileError(defname);}
        } else { info = ReadDefFileError(xNameListFile);}
      }
      /*cisajs.def----------------------------------------*/
      if(info==0) {
        if(fscanf(fplist, "%s\n", defname)!=EOF) {
          fp = fopen(defname, "r");
          if(fp!=NULL) {
            fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
            fscanf(fp,"%s %d\n", ctmp, &(bufInt[26])); /* NCisAjs */
            fclose(fp);
          } else { info = ReadDefFileError(defname);}
        } else { info = ReadDefFileError(xNameListFile);}
      }
      /*cisajscktalt.def----------------------------------*/
      if(info==0) {
        if(fscanf(fplist, "%s\n", defname)!=EOF) {
          fp = fopen(defname, "r");
          if(fp!=NULL) {
            fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
            fscanf(fp,"%s %d\n", ctmp, &(bufInt[27])); /* NCisAjsCktAlt */
            fclose(fp);
          } else { info = ReadDefFileError(defname);}
        } else { info = ReadDefFileError(xNameListFile);}
      }

      fclose(fplist);
    } else { info=ReadDefFileError(xNameListFile);}
  } /* if(rank==0) */

  if(info!=0) {
    if(rank==0) {
      fprintf(stderr, "error: Definition files(*.def) are incomplete.\n");
    }
    MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
  }

#ifdef _mpi_use
  MPI_Bcast(bufInt, nBufInt, MPI_INT, 0, comm);
  MPI_Bcast(bufDouble, nBufDouble, MPI_DOUBLE, 0, comm);
  MPI_Bcast(CDataFileHead, nBufChar, MPI_CHAR, 0, comm);
  MPI_Bcast(CParaFileHead, nBufChar, MPI_CHAR, 0, comm);
#endif /* _mpi_use */

  NVMCCalMode            =  bufInt[ 0];
  /* NLanczosMode           =  bufInt[ 1]; */
  NDataIdxStart          =  bufInt[ 2];
  NDataQtySmp            =  bufInt[ 3];
  Nsite                  =  bufInt[ 4];
  Ne                     =  bufInt[ 5];
  NSPGaussLeg            =  bufInt[ 6];
  NSPStot                =  bufInt[ 7];
  NMPTrans               =  bufInt[ 8];
  NSROptItrStep          =  bufInt[ 9];
  NSROptItrSmp           =  bufInt[10];
  NSROptFixSmp           =  bufInt[11];
  NVMCWarmUp             =  bufInt[12];
  NVMCIniterval          =  bufInt[13];
  NVMCSample             =  bufInt[14];
  NExUpdatePath          =  bufInt[15];
  RndSeed                =  bufInt[16];
  NSplitSize             =  bufInt[17];
  NLocSpn                =  bufInt[18];
  NTransfer              =  bufInt[19];
  NCoulomb               =  bufInt[20];
  NInteraction           =  bufInt[21];
  NGutzwillerIdx         =  bufInt[22];
  NJastrowIdx            =  bufInt[23];
  NOrbitalIdx            =  bufInt[24];
  NQPTrans               =  bufInt[25];
  NCisAjs                =  bufInt[26];
  NCisAjsCktAlt          =  bufInt[27];

  DSROptRedCut = bufDouble[0];
  DSROptStaDel = bufDouble[1];
  DSROptStepDt = bufDouble[2];

  Nsize   = 2*Ne;
  Nsite2  = 2*Nsite;
  NSlater = NOrbitalIdx;
  NProj   = NGutzwillerIdx + NJastrowIdx;
  NPara   = NSlater + NProj;
  NQPFull = NSPGaussLeg * NMPTrans;
  SROptSize = NPara+1;

  NTotalDefInt = Nsite /* LocSpn */
    + 3*NTransfer /* Transfer */
    + 2*NCoulomb /* Coulomb */
    + 6*NInteraction /* Interaction */
    + Nsite /* GutzwillerIdx */
    + Nsite*Nsite /* JastrowIdx */
    + Nsite*Nsite /* OrbitalIdx */
    + Nsite*NQPTrans /* QPTrans */
    + 3*NCisAjs /* CisAjs */
    + 6*NCisAjsCktAlt /* CisAjsCktAlt */
    + NPara; /* OptFlag */
  NTotalDefDouble = NTransfer /* ParaTransfer */
    + NCoulomb /* ParaCoulomb */
    + NInteraction /* ParaInteraction */
    + NQPTrans; /* ParaQPTrans */

  return 0;
}

int ReadDefFileIdxPara(char *xNameListFile, MPI_Comm comm){
  FILE *fp, *fplist;
  char defname[D_FileNameMax];
  char ctmp[D_FileNameMax];
  int itmp;

  int i,j,idx,idx0,idx1,info=0;
  int fidx=0; /* index for OptFlag */
  int x0,x1,x2,x3,x4,x5;
  int rank;

  MPI_Comm_rank(comm, &rank);

  if(rank==0) {
    fplist = fopen(xNameListFile, "r");
    if(fplist!=NULL) {
      /*modpara.def---------------------------------------*/
      fscanf(fplist, "%s\n", defname);

      /*locspn.def----------------------------------------*/
      if(info==0) {
        if(fscanf(fplist, "%s\n", defname)!=EOF) {
          fp = fopen(defname, "r");
          if(fp!=NULL) {
            for(i=0;i<5;i++) fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
            idx = 0;
            while( fscanf(fp, "%d %d\n", &(x0), &(x1) )!=EOF){
              LocSpn[x0] = x1;
              idx++;
            }
            if(NLocSpn>2*Ne){
              fprintf(stderr, "Error: 2*Ne must be (2*Ne >= NLocalSpin).\n");
              info=1;
            }
            if(NLocSpn>0 && NExUpdatePath==0){
              fprintf(stderr, "Error: NExUpdatePath (in modpara.def) must be 1.\n");
              info=1;
            }
            if(idx!=Nsite) info=ReadDefFileError(defname);
            fclose(fp);
          } else { info = ReadDefFileError(defname);}
        } else { info = ReadDefFileError(xNameListFile);}
      }

      /*transfer.def--------------------------------------*/
      if(info==0) {
        if(fscanf(fplist, "%s\n", defname)!=EOF) {
          if(NTransfer>0){
            fp = fopen(defname, "r");
            if(fp!=NULL) {
              for(i=0;i<5;i++) fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
              idx = 0;
              while( fscanf(fp, "%d %d %d %lf\n", 
                            &(Transfer[idx][0]),
                            &(Transfer[idx][1]),
                            &(Transfer[idx][2]),
                            &(ParaTransfer[idx]))!=EOF){
                idx++;
              }
              if(idx!=NTransfer) info = ReadDefFileError(defname);
              fclose(fp);
            } else { info = ReadDefFileError(defname);}
          }
        } else { info = ReadDefFileError(xNameListFile);}
      }

      /*coulomb.def----------------------------------*/
      if(info==0) {
        if(fscanf(fplist, "%s\n", defname)!=EOF) {
          if(NCoulomb>0){
            fp = fopen(defname, "r");
            if(fp!=NULL) {
              for(i=0;i<5;i++) fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
              idx = 0;
              while( fscanf(fp, "%d %d %d %d %lf\n", 
                            &(x0),&(x1),&(x2),&(x3),
                            &(ParaCoulomb[idx]) )!=EOF){
                Coulomb[idx][0] = x0+x1*Nsite;
                Coulomb[idx][1] = x2+x3*Nsite;
                idx++;
              }
              if(idx!=NCoulomb) info=ReadDefFileError(defname);
              fclose(fp);
            } else { info=ReadDefFileError(defname);}
          }
        } else { info=ReadDefFileError(xNameListFile);}
      }

      /*interaction.def---------------------------------------*/
      if(info==0) {
        if(fscanf(fplist, "%s\n", defname)!=EOF) {
          if(NInteraction>0){
            fp = fopen(defname, "r");
            if(fp!=NULL) {
              for(i=0;i<5;i++) fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
              idx = 0;
              while( fscanf(fp, "%d %d %d %d %d %d %lf\n", 
                            &(Interaction[idx][0]),
                            &(Interaction[idx][1]),
                            &(Interaction[idx][2]),
                            &(Interaction[idx][3]),
                            &(Interaction[idx][4]),
                            &(Interaction[idx][5]),
                            &(ParaInteraction[idx]) )!=EOF ){
                idx++;
              }
              if(idx!=NInteraction) info=ReadDefFileError(defname);
              fclose(fp);
            } else { info=ReadDefFileError(defname);}
          }
        } else {
          info=ReadDefFileError(xNameListFile);
        }
      }

      /*gutzwilleridx.def---------------------------------*/
      if(info==0) {
        if(fscanf(fplist, "%s\n", defname)!=EOF) {
          if(NGutzwillerIdx>0){
            fp = fopen(defname, "r");
            if(fp!=NULL) {
              for(i=0;i<5;i++) fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
              idx0 = idx1 = 0;
              while( fscanf(fp, "%d ", &i) != EOF){
                fscanf(fp, "%d\n", &(GutzwillerIdx[i]));
                idx0++;
                if(idx0==Nsite) break;
              }
              while( fscanf(fp, "%d ", &i) != EOF){
                fscanf(fp, "%d\n", &(OptFlag[fidx]));
                fidx++;
                idx1++;
              }
              if(idx0!=Nsite || idx1!=NGutzwillerIdx) {
                info=ReadDefFileError(defname);
              }
              fclose(fp);
            } else { info=ReadDefFileError(defname);}
          }
        } else { info=ReadDefFileError(xNameListFile);}
      }

      /*jastrowidx.def------------------------------------*/
      if(info==0) {
        if(fscanf(fplist, "%s\n", defname)!=EOF) {
          if(NJastrowIdx>0){
            fp = fopen(defname, "r");
            if(fp!=NULL) {
              for(i=0;i<5;i++) fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
              idx0 = idx1 = 0;
              while( fscanf(fp, "%d %d ", &i, &j) != EOF){
                if(i==j){
                  fprintf(stderr, "Error in %s: [Condition] i neq j\n", defname);
                  info=1;
                  break;
                }
                fscanf(fp, "%d\n", &(JastrowIdx[i][j]));
                JastrowIdx[i][i] = -1; // This case is Gutzwiller.
                idx0++;
                if(idx0==Nsite*(Nsite-1)) break;
              }
              while( fscanf(fp, "%d ", &i) != EOF){
                fscanf(fp, "%d\n", &(OptFlag[fidx]));
                fidx++;
                idx1++;
              }
              if(idx0!=Nsite*(Nsite-1) || idx1!=NJastrowIdx) {
                info=ReadDefFileError(defname);
              }
              fclose(fp);
            } else { info=ReadDefFileError(defname);}
          }
        } else { info=ReadDefFileError(xNameListFile);}
      }

      /*orbitalidx.def------------------------------------*/
      if(info==0) {
        if(fscanf(fplist, "%s\n", defname)!=EOF) {
          if(NOrbitalIdx>0){
            fp = fopen(defname, "r");
            if(fp!=NULL) {
              for(i=0;i<5;i++) fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
              idx0 = idx1 = 0;
              while( fscanf(fp, "%d %d ", &i, &j) != EOF){
                fscanf(fp, "%d\n", &(OrbitalIdx[i][j]));
                idx0++;
                if(idx0==Nsite*Nsite) break;
              }
              while( fscanf(fp, "%d ", &i) != EOF){
                fscanf(fp, "%d\n", &(OptFlag[fidx]));
                fidx += 1;
                idx1++;
              }
              if(idx0!=Nsite*Nsite || idx1!=NOrbitalIdx) {
                info=ReadDefFileError(defname);
              }
              fclose(fp);
            } else { info=ReadDefFileError(defname);}
          }
        } else { info=ReadDefFileError(xNameListFile);}
      }

      /*qptransidx.def------------------------------------*/
      if(info==0) {
        if(fscanf(fplist, "%s\n", defname)!=EOF) {
          if(NQPTrans>0){
            fp = fopen(defname, "r");
            if(fp!=NULL) {
              for(i=0;i<5;i++) fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
              for(i=0;i<NQPTrans;i++){
                fscanf(fp, "%d ", &itmp);
                fscanf(fp, "%lf\n", &(ParaQPTrans[itmp]));
              }
              idx = 0;
              while( fscanf(fp, "%d %d ", &i, &j) != EOF){
                fscanf(fp, "%d\n", &(QPTrans[i][j]));
                idx++;
              }
              if(idx!=NQPTrans*Nsite) info=ReadDefFileError(defname);
              fclose(fp);
            } else { info=ReadDefFileError(defname);}
          }
        } else { info=ReadDefFileError(xNameListFile);}
      }

      /*cisajs.def----------------------------------------*/
      if(info==0){
        if(fscanf(fplist, "%s\n", defname)!=EOF) {
          if(NCisAjs>0){
            fp = fopen(defname, "r");
            if(fp!=NULL) {
              for(i=0;i<5;i++) fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
              idx = 0;
              while( fscanf(fp, "%d %d %d %d\n",
                            &(x0), &(x1), &(x2), &(x3)) != EOF){
                CisAjsIdx[x0][0] = x1;
                CisAjsIdx[x0][1] = x2;
                CisAjsIdx[x0][2] = x3;
                idx++;
              }
              if(idx!=NCisAjs) info=ReadDefFileError(defname);
              fclose(fp);
            } else { info=ReadDefFileError(defname);}
          }
        } else { info=ReadDefFileError(xNameListFile);}
      }

      /*cisajscktalt.def--------------------------------*/
      if(info==0) {
        if(fscanf(fplist, "%s\n", defname)!=EOF) {
          if(NCisAjsCktAlt>0){
            fp = fopen(defname, "r");
            if(fp!=NULL) {
              for(i=0;i<5;i++) fgets(ctmp, sizeof(ctmp)/sizeof(char), fp);
              idx = 0;
              while( fscanf(fp, "%d %d %d %d %d %d\n", 
                            &(x0), &(x1), &(x2), &(x3), &(x4), &(x5) ) != EOF ){
                CisAjsCktAltIdx[idx][0] = x0;
                CisAjsCktAltIdx[idx][1] = x1;
                CisAjsCktAltIdx[idx][2] = x2;
                CisAjsCktAltIdx[idx][3] = x3;
                CisAjsCktAltIdx[idx][4] = x4;
                CisAjsCktAltIdx[idx][5] = x5;
                idx++;
              }
              if(idx!=NCisAjsCktAlt) info=ReadDefFileError(defname);
              fclose(fp);
            } else { info=ReadDefFileError(defname);}
          }
        } else { info=ReadDefFileError(xNameListFile);}
      }

      fclose(fplist);

    } else { info = ReadDefFileError(xNameListFile);}

    if(fidx!=NPara){
      fprintf(stderr, "error: OptFlag is incomplete.\n");
      info=1;
    }
  } /* if(rank==0) */

  if(info!=0) {
    if(rank==0) {
      fprintf(stderr, "error: Indices and Parameters of Definition files(*.def) are incomplete.\n");
    }
    MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
  }

#ifdef _mpi_use
  SafeMpiBcastInt(LocSpn, NTotalDefInt, comm);
  SafeMpiBcast(ParaTransfer, NTotalDefDouble, comm);
#endif /* _mpi_use */

  return 0;
}
