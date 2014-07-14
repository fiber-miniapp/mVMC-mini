#!/usr/bin/env python
import sys,time,os,random
from math import *
import numpy as np

T={(1,0): 1.0, (-1,0): 1.0,
   (0,1): 1.0, (0,-1): 1.0}
#J=1.45
J=1.0

def position(r):
    x=(r[0]+Lx)%Lx
    y=(r[1]+Ly)%Ly
    return (x,y)

def indexToPosition(i):
    x=i%Lx
    y=i/Lx
    return (x,y)

def positionToIndex(r):
    x=(r[0]+Lx)%Lx
    y=(r[1]+Ly)%Ly
    return x+Lx*y

def neighborIndex(i,dr):
    bi,ri = i/Nlattice, i%Nlattice
    r=indexToPosition(ri)
    x=r[0]+dr[0]
    y=r[1]+dr[1]
    return positionToIndex([x,y])+bi*Nlattice

def direction(i,j):
    ri=indexToPosition(i)
    rj=indexToPosition(j)
    dx=(rj[0]-ri[0]+Lx)%Lx
    dy=(rj[1]-ri[1]+Ly)%Ly
    return (dx,dy)

def locgrnIdx(i,j,s):
    return (Nsite*i+j)+s*Nsite*Nsite

def subIndex(i):
    r=indexToPosition(i)
    sx=r[0]%Sx
    sy=r[1]%Sy
    return sx+Sx*sy

def rotateIndex(i, deg):
    bi,ri = i/Nlattice, i%Nlattice
    rx,ry = indexToPosition(ri)
    if deg == 1:   x,y = -ry,rx
    elif deg == 2: x,y = -rx,-ry
    elif deg == 3: x,y =  ry,-rx
    else: x,y = rx,ry
    return positionToIndex([x,y])+bi*Nlattice

if len(sys.argv)<2:
    print "./makeJob.py num_mpi_process"
    sys.exit()

Lx = Ly = 12
J = 1.0
Nsite = 2*Lx*Ly
Nlattice = Lx*Ly
Nelectron = Nsite/2 # half-filled

Sx = 4
Sy = 4
Nsub = Sx*Sy
NQPTrans = 4*Nsub

NProcess = int(sys.argv[1])
NTotalSample = 12288
NOuterMPI = 64
NSplitSize = NProcess/NOuterMPI
NVMCSample = NTotalSample/NOuterMPI

if NProcess<NOuterMPI:
    NSplitSize=1
    NVMCSample = NTotalSample/NProcess

# seed = int(time.time())
seed = 3141592
random.seed()

dir='job_mpi{0}'.format(NProcess)
if os.access(dir, os.F_OK):
    print "error: '{0}' exists".format(dir)
    sys.exit()
else:
    print dir
    os.mkdir(dir)
os.chdir(dir)

dir='Lx{0}Ly{1}_J{2}'.format(Lx,Ly,J)
if os.access(dir, os.F_OK):
    print "error: '{0}' exists".format(dir)
    sys.exit()
else:
    print dir
    os.mkdir(dir)

### job.sh ###
f = open("job.sh",'w')
f.write(
    "#!/bin/bash\n\n"+
    "export OMP_NUM_THREADS=8\n\n"+
    "mpiexec -np {0} ./vmc.out multiDir.def\n".format(NProcess)
    )
f.close()
os.chmod("job.sh",0755)

### multiDir.def ###
f = open("multiDir.def",'w')
f.write(
    "1\n"+
    "./{0} xnamelist.def opt.init\n".format(dir)
    )
f.close()

os.chdir(dir)

separator = '--------------------\n'

pre='z'
fileName = []
fileName.append('xnamelist.def')
fileName.append(pre+'modpara.def')
fileName.append(pre+'locspn.def')
fileName.append(pre+'transfer.def')
fileName.append(pre+'coulomb.def')
fileName.append(pre+'interaction.def')
fileName.append(pre+'gutzwilleridx.def')
fileName.append(pre+'jastrowidx.def')
fileName.append(pre+'orbitalidx.def')
fileName.append(pre+'qptransidx.def')
fileName.append(pre+'cisajs.def')
fileName.append(pre+'cisajscktalt.def')

### xnamelist.def ###
f = open(fileName[0],'w')
for x in fileName[1:]:
    f.write(x+"\n")
f.close()

### modpara.def ###
f = open(fileName[1],'w')
f.write(
    separator+
    "Model_Parameters  0\n"+
    separator+
    "VMC_Cal_Parameters\n"+
    separator+
    "CDataFileHead  zvo\n"+
    "CParaFileHead  zqp\n"+
    separator+
    "NVMCCalMode    0\n"+
    "NLanczosMode   0\n"+
    separator+
    "NDataIdxStart  0\n"+
    "NDataQtySmp    10\n"+
    separator+
    "Nsite          {0}\n".format(Nsite)+
    "Nelectron      {0}\n".format(Nelectron)+
    "NSPGaussLeg    16\n"+
    "NSPStot        0\n"+
    "NMPTrans       {0}\n".format(NQPTrans)+
    "NSROptItrStep  20\n"+
    "NSROptItrSmp   5\n"+
    "NSROptFixSmp   1\n"+
    "DSROptRedCut   1e-10\n"+
    "DSROptStaDel   0.01\n"+
    "DSROptStepDt   0.01\n"+
    "NVMCWarmUp     10\n"+
    "NVMCIniterval  1\n"+
    "NVMCSample     {0}\n".format(NVMCSample)+
    "NExUpdatePath  1\n"+
    "RndSeed        {0}\n".format(seed)+
    "NSplitSize     {0}\n".format(NSplitSize))
f.close()

### locspn.def ###
f = open(fileName[2],'w')
f.write(separator+
        "NLocalSpin\t{0}\n".format(Nlattice)+
        separator+
        "i_0LocSpn_1IteElc\n"+
        separator)
for i in range(Nsite):
    if i<Nlattice: f.write("{0}\t1\n".format(i)) # conduction electrons
    else: f.write("{0}\t0\n".format(i)) # localized electrons
f.close()

### transfer.def ###
paraList = []
for ri in range(Nlattice):
    for dr,t in T.iteritems():
        rj = neighborIndex(ri,dr)
        paraList.append([ri,rj,t])
paraList.sort()

f = open(fileName[3],'w')
f.write(separator+
        "NTransfer\t"+str(2*len(paraList))+"\n"+
        separator+
        "i_j_s_tijs\n"+
        separator)
for s in [0,1]:
    for para in paraList:
        f.write("{0}\t{1}\t{2}\t{3}\n".format(para[0],para[1], s, para[2]))
f.close()

## J ##
paraList = []
for i in range(Nlattice):
    paraList.append([i,i+Nlattice,J])
paraList.sort()

### coulomb.def ###
f = open(fileName[4],'w')
f.write(separator+
        "NCoulomb\t"+ str(4*len(paraList)) +"\n"+
        separator+
        "i_s_j_t_Visjt\n"+
        separator)
for para in paraList:
    f.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(para[0],0,para[1],0,0.25*para[2]))
    f.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(para[0],0,para[1],1,-0.25*para[2]))
    f.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(para[0],1,para[1],0,-0.25*para[2]))
    f.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(para[0],1,para[1],1,0.25*para[2]))
f.close()

### interaction.def ###
f = open(fileName[5],'w')
f.write(separator+
        "NInteraction\t"+ str(2*len(paraList)) +"\n"+
        separator+
        "i_j_s_k_l_t_Jijsklt\n"+
        separator)
for para in paraList:
    f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(para[0],para[1],0,para[1],para[0],1,-0.5*para[2]))
    f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(para[0],para[1],1,para[1],para[0],0,-0.5*para[2]))
f.close()

### gutzwilleridx.def ###
NGutzwiller = 2
f = open(fileName[6],'w')
f.write(separator+
        "NGutzwillerIdx\t"+ str(NGutzwiller) +"\n"+
        separator+
        "i_GutzwillerIdx\n"+
        separator)
for i in range(Nsite):
    if i<Nlattice: f.write("{0}\t{1}\n".format(i,0))
    else: f.write("{0}\t{1}\n".format(i,1))
for idx in range(NGutzwiller):
    f.write("{0}\t{1}\n".format(idx,1)) # optimized
f.close()

### jastrow.def ###
jastrow = {}
idx=-1
for i in range(Nlattice):
    dr = indexToPosition(i)
    dr_rev = (Lx-dr[0], Ly-dr[1])
    if dr_rev in jastrow:
        jastrow[dr] = jastrow[dr_rev]
    else:
        if dr[0]==0 and dr[1] > Ly/2: 
            jastrow[dr] = jastrow[(dr[0], Ly-dr[1])]
            continue
        elif dr[1]==0 and dr[0] > Lx/2:
            jastrow[dr] = jastrow[(Lx-dr[0], dr[1])]
            continue
        else:
            jastrow[dr] = idx
            idx += 1
NJastrow = idx+1 # Ns/2+1 (Lx,Ly:even), Ns/2 (Lx:even, Ly:odd)
f = open(fileName[7],'w')
f.write(separator+
        "NJastrowIdx\t"+ str(NJastrow) +"\n"+
        separator+
        "i_j_JastrowIdx\n"+
        separator)
for i in range(Nsite):
    for j in range(Nsite):
        if i==j: continue
        elif i<Nlattice and j<Nlattice:
            dr=direction(i,j)
            idx = jastrow[dr]
            f.write("{0}\t{1}\t{2}\n".format(i,j,idx))
        else:
            f.write("{0}\t{1}\t{2}\n".format(i,j,NJastrow-1))
for idx in range(NJastrow):
    f.write("{0}\t{1}\n".format(idx,1)) #optimized
f.close()

### orbitalidx.def ###
NOrbital = 4*Nsub*Nlattice
f = open(fileName[8],'w')
f.write(separator+
        "NOrbitalIdx\t"+ str(NOrbital) +"\n"+
        separator+
        "i_j_OrbitalIdx\n"+
        separator)
orbital=[[dj+isub*Nlattice for dj in range(Nlattice)] for isub in range(Nsub)]
for i in range(Nsite):
    for j in range(Nsite):
        bi,ri = i/Nlattice, i%Nlattice
        bj,rj = j/Nlattice, j%Nlattice
        isub = subIndex(ri)
        dj = positionToIndex(direction(ri,rj))
        idx = orbital[isub][dj] + (2*bi+bj)*Nsub*Nlattice
        f.write("{0}\t{1}\t{2}\n".format(i,j,idx))
for idx in range(NOrbital):
    f.write("{0}\t{1}\n".format(idx,1))
f.close()

### qptransidx.def ###
f = open(fileName[9],'w')
f.write(separator+
        "NQPTrans\t"+ str(NQPTrans) +"\n"+
        separator+
        "TrIdx_TrWeight_and_TrIdx_i_xi\n"+
        separator)
for idx in range(NQPTrans):
    f.write("{0}\t{1}\n".format(idx, 1.0))
for idx in range(NQPTrans):
    si, deg = idx%Nsub, idx/Nsub
    dr = (si%Sx, si/Sx)
    for i in range(Nsite):
        j = rotateIndex(i,deg)
        j = neighborIndex(j,dr)
        f.write("{0}\t{1}\t{2}\n".format(idx, i, j))
f.close()

### cisajs.def ###
f = open(fileName[10],'w')
f.write(separator+
        "NCisAjs\t"+ str(Nsite*Nsite*2) +"\n"+
        separator+
        "idx_i_j_s\n"+
        separator)
for i in range(Nsite):
    for j in range(Nsite):
        f.write("{0}\t{1}\t{2}\t{3}\n".format(locgrnIdx(i,j,0), i, j, 0))
for i in range(Nsite):
    for j in range(Nsite):
        f.write("{0}\t{1}\t{2}\t{3}\n".format(locgrnIdx(i,j,1), i, j, 1))
f.close()

### cisajscktaltdc.def ###
f = open(fileName[11],'w')
f.write(separator+
        "NCisAjsCktAltDC\t"+ str(Nlattice) +"\n"+
        separator+
        "i_j_s_k_l_t\n"+
        separator)
for i in range(Nlattice):
    f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(i,i,0,i,i,1))
f.close()

########## Phys Cal ##########
fileName = []
fileName.append('xnamelist_aft.def')
fileName.append(pre+'modpara_aft.def')
fileName.append(pre+'locspn.def')
fileName.append(pre+'transfer.def')
fileName.append(pre+'coulomb.def')
fileName.append(pre+'interaction.def')
fileName.append(pre+'gutzwilleridx.def')
fileName.append(pre+'jastrowidx.def')
fileName.append(pre+'orbitalidx.def')
fileName.append(pre+'qptransidx.def')
fileName.append(pre+'cisajs.def')
fileName.append(pre+'cisajscktalt.def')

### xnamelist_aft.def ###
f = open(fileName[0],'w')
for x in fileName[1:]:
    f.write(x+"\n")
f.close()

### modpara_aft.def ###
f = open(fileName[1],'w')
f.write(
    separator+
    "Model_Parameters  0\n"+
    separator+
    "VMC_Cal_Parameters\n"+
    separator+
    "CDataFileHead  zvo_aft\n"+
    "CParaFileHead  zqp_aft\n"+
    separator+
    "NVMCCalMode    1\n"+
    "NLanczosMode   0\n"+
    separator+
    "NDataIdxStart  10\n"+
    "NDataQtySmp    10\n"+
    separator+
    "Nsite          {0}\n".format(Nsite)+
    "Nelectron      {0}\n".format(Nelectron)+
    "NSPGaussLeg    16\n"+
    "NSPStot        0\n"+
    "NMPTrans       {0}\n".format(NQPTrans)+
    "NSROptItrStep  2000\n"+
    "NSROptItrSmp   500\n"+
    "NSROptFixSmp   1\n"+
    "DSROptRedCut   1e-10\n"+
    "DSROptStaDel   0.01\n"+
    "DSROptStepDt   0.01\n"+
    "NVMCWarmUp     1000\n"+
    "NVMCIniterval  10\n"+
    "NVMCSample     {0}\n".format(NVMCSample)+
    "NExUpdatePath  1\n"+
    "RndSeed        {0}\n".format(seed+1000)+
    "NSplitSize     {0}\n".format(NSplitSize))
f.close()
