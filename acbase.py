#!/usr/bin/python

import csv
import sys, os
import traceback
import math
import cmath
from myutils import *
from log import *
import data
import numpy
import scipy
from acn import *
#from knitro import *
#from knitroNumPy import *
from amplpy import *
#from line_profiler import *
import time
#import evaluation2
import acn_evaluation2
import multiprocessing
from acn_solution import *
from acnsolverAMPL import acnexposevars, acnsolveit, resolve_fixed, write_base_files


def acbase(acn, procid):

 log = acn.alldata['log']
 log.joint("\n***\n solving base\n")

 # Get info on cores and set parallelism
 numcores = multiprocessing.cpu_count()
 print('numcores=',numcores)

 # use this to suppress stdout in parallel contingency loop
 nooutput = 0
 if nooutput:
     sys.stdout =  open(os.devnull, 'w')
     sys.stderr =  open(os.devnull, 'w')

 # Initialize
 division = int(acn.division)
 if (division == 2 or division == 4 or (division == 1 and acn.numbuses < 10000) or (division == 3 and acn.numbuses < 10000)):
     if acn.numnodes == 1:
        max_passes = 1 # 2
     else:
        max_passes = 1
 else:
     max_passes = 1
 passnum = 1
 stop = 0
 solvenum = 0
 maxsolves = 6

 # Create some numpy 2-d arrays for storing solutions for each solve
 # The # of rows in each arrays is maxsolves
 acn.basesol_busVmags = np.zeros((maxsolves,acn.numbuses), dtype=float)
 acn.basesol_busVangles = np.zeros((maxsolves,acn.numbuses), dtype=float)
 acn.basesol_clearedloads = np.zeros((maxsolves,acn.numloads+1), dtype=float)
 acn.basesol_ploads = np.zeros((maxsolves,acn.numloads+1), dtype=float)
 acn.basesol_genP = np.zeros((maxsolves,acn.numgens), dtype=float)
 acn.basesol_genQ = np.zeros((maxsolves,acn.numgens), dtype=float)
 acn.basesol_genON = np.zeros((maxsolves,acn.numgens), dtype=int)
 acn.basesol_gensu = np.zeros((maxsolves,acn.numgens), dtype=int)
 acn.basesol_gensd = np.zeros((maxsolves,acn.numgens), dtype=int)
 acn.basesol_ntrON = np.zeros((maxsolves,acn.nontranscount_original), dtype=int)
 acn.basesol_trON = np.zeros((maxsolves,acn.transcount_original), dtype=int)
 acn.basesol_taps = np.zeros((maxsolves,acn.transcount_original), dtype=int)
 # Get size of basesw so we can create corresponding numpy array
 numxha0 = 0
 for h in range(1,1 + acn.numswitchedshunts):
   oursh = acn.ourswshunts[h]
   if oursh.status:
     for a in range(1,1+oursh.len):
       if oursh.n[a] > 0:
         numxha0 += 1
 #print('numxha0=',numxha0)
 #breakexit('numxha0')
 acn.numxha0 = numxha0
 acn.basesol_sw = np.zeros((maxsolves,numxha0), dtype=int)

 # BASE SOLVE LOOP BEGIN
 useparallel = 0
 while stop == 0:
   if useparallel == 0:
     log.joint('---> sequential base solves\n')
     # Sequential implementation of base solve loop
     for h in range(maxsolves):
         solvenum = h
         acbasesub(acn, solvenum, passnum)
         #breakexit('passnum')
         stop = 0
         passnum += 1
         elapsedtime = time.time()-acn.timebeg
         remaintime = acn.maxtime - elapsedtime
         if remaintime < 30:
           stop = 1
         if passnum > max_passes:
           stop = 1
         #stop = 0 # use to force stop
   else:
     # Parallel implementation of base solve loop using forks
     # Currently fork as many child processes as threads/cores available
     log.joint('---> parallel base solves\n')
     #TBD
     for i in range(numcores):
         try:
             pid = os.fork()
         except:
             #print("Could not create a child process")
             continue
         if pid == 0:
             for h in range(maxsolves):
                 if 1: # may impose condition here later
                     elapsedtime = time.time()-acn.timebeg
                     remaintime = acn.maxtime - elapsedtime
                     needtosolve = 1
                     modulo = (h+1)%numcores
                     #if passnum > 2 and numnegscores > 0 and acn.scores[h] >= 0:
                     #    needtosolve = 0 #if negative scores, focus on improving only these for now
                     limit = 60.0
                     #print('i=',i,' modulo
                     if i==modulo and needtosolve == 1 and remaintime >= limit:
                         try:
                             log.joint("Solve base "+str(h)+" with procid "+str(procid)+" core "+str(i)+" remain_time="+str(remaintime)+"\n")
                             acbasesub(acn, h, passnum)
                         except:
                             log.joint("base solve "+str(h)+" encountered a problem inside acbasesub!\n")
                             #exit()
                             pass
             os._exit(os.EX_OK)
         else:
             pass
             #print("In the parent process after forking the child {}".format(pid))
         #breakexit("loop")
     #breakexit('parallel contingency loop')
     for i in range(numcores):
         finished = os.waitpid(0, 0)
     stop=1

 # BASE SOLVE LOOP FINISHED

 if nooutput:
     sys.stdout = sys.__stdout__
     sys.stderr = sys.__stderr__

 # Now evaluate all base solutions and return the best one
 for solvenum in range(maxsolves):
     basefilename = "base_" + str(solvenum) + ".txt"
     try:
         pobj = -1e20
         log.joint("\nEvaluating base_*.txt...\n")
         base.obj, base.infeas, exist, base.summary = acn_evaluation2.base_from_file(acn, basefilename)
         log.joint("\nEvaluating solution_BASECASE.txt...\n")
         pobj, pinfeas, pexist, base.psummary = acn_evaluation2.base_from_file(acn, "solution_BASECASE.txt")
         log.joint("\n------------------------------------------\n")
         if (pexist == True):
             log.joint('current solution: obj=' + str(pobj)+'\n')
             log.joint('current solution: infeas=' + str(pinfeas)+'\n')
         log.joint('base exists=' +  str(exist)+'\n')
         log.joint('solution_BASECASE exists=' +  str(pexist)+'\n')
         if (exist == True):
             log.joint('base.obj=' + str(base.obj)+'\n')
             log.joint('base.infeas=' + str(base.infeas)+'\n')
             if base.infeas == 0 and (pexist == False or base.obj > pobj):
                 log.joint('We found a base solution ' +str(solvenum)+' better than current solution.\n')
                 # We've computed a base case solution that is better than the
                 # current best solution, so overwrite with this one
                 # TBD: BETTER WAY TO HANDLE THIS COPY!
                 if solvenum==0:
                   os.popen('\cp base_0.txt solution_BASECASE.txt')
                   os.popen('\cp base_0.bin solution_BASECASE.bin')
                 elif solvenum==1:
                   os.popen('\cp base_1.txt solution_BASECASE.txt')
                   os.popen('\cp base_1.bin solution_BASECASE.bin')
                 elif solvenum==2:
                   os.popen('\cp base_2.txt solution_BASECASE.txt')
                   os.popen('\cp base_2.bin solution_BASECASE.bin')
                 elif solvenum==3:
                   os.popen('\cp base_3.txt solution_BASECASE.txt')
                   os.popen('\cp base_3.bin solution_BASECASE.bin')
                 elif solvenum==4:
                   os.popen('\cp base_4.txt solution_BASECASE.txt')
                   os.popen('\cp base_4.bin solution_BASECASE.bin')
                 elif solvenum==5:
                   os.popen('\cp base_5.txt solution_BASECASE.txt')
                   os.popen('\cp base_5.bin solution_BASECASE.bin')
                 # Now write out individual base solution components
                 # needed for contingency modeling to files so that they can be read in by
                 # MyPython2.py->acnsolverAMPL2.py
                 acn.basebusVmags = acn.basesol_busVmags[solvenum,:]
                 acn.basebusVangles = acn.basesol_busVangles[solvenum,:]
                 acn.clearedloads = acn.basesol_clearedloads[solvenum,:]
                 acn.ploads = acn.basesol_ploads[solvenum,:]
                 acn.basegenP = acn.basesol_genP[solvenum,:]
                 acn.basegenQ = acn.basesol_genQ[solvenum,:]
                 acn.basegenON = acn.basesol_genON[solvenum,:]
                 acn.basegensu = acn.basesol_gensu[solvenum,:]
                 acn.basegensd = acn.basesol_gensd[solvenum,:]
                 acn.basentrON = acn.basesol_ntrON[solvenum,:]
                 acn.basetrON = acn.basesol_trON[solvenum,:]
                 acn.basetaps = acn.basesol_taps[solvenum,:]
                 #print('basesw0=',acn.basesw0)
                 j = 0
                 acn.basesw={}
                 for (h,a) in acn.basesw0:
                   acn.basesw[h,a] = acn.basesol_sw[solvenum,j]
                   j=j+1
                 #print('basesw=',basesw)
                 #breakexit('zzz')
                 write_base_files(acn,acn.basebusVmags,acn.basebusVangles,acn.ploads,acn.basegenP,acn.basegenON,acn.basegensu,acn.basegensd,acn.basentrON,acn.basetrON,acn.basesw)
                 acn.improvedbase = 1
     except:
         log.joint("Error evaluating the base solution!\n")
         var = traceback.format_exc()
         traceback.print_exc()
         acn_evaluation2.print_alert(var, raise_exception=False)
         pass

def acbasesub(acn, solvenum, passnum):

 log  = acn.alldata['log']
 log.joint("\n***\n running base solution heuristic 1\n")

 division = int(acn.division)

 print("*** solvenum: ", solvenum, " pass: ",passnum, " division=",division)
 #breakexit("pass")

 heatmap = acn.heatmap
 totviolvector = acn.totviolvector


 numbuses = acn.numbuses
 actionable1 = np.zeros(1+numbuses, dtype = 'I')
 sortedorder = np.argsort(-totviolvector)

 N = 100
 if N > numbuses:
     N = numbuses

 showthresh = 1e-04
 totactionable = 0
 log.joint(' top %d violation buses with viol > %g\n' %(N,showthresh))
 for j in range(N):
     k = sortedorder[j]
     if totviolvector[k] > showthresh:
      log.joint('%d. count1 %d i %d viol %.4e heat %g pv %g qv %g\n' %(j,k+1,acn.ourbuses[k+1].i, totviolvector[k], heatmap[k], acn.pimb[k], acn.qimb[k]))
      actionable1[k+1] = 1
      totactionable += 1

 log.joint(' total actionable: ' + str(totactionable) + '\n')
 priorsol = acn.options['priorsol']

 alldata = acn.alldata
 raw = alldata['RAW']
 con = alldata['CON']
 sup = alldata['JSON']

 pi = math.pi

 numpy.set_printoptions(precision=16)
 elapsedtime = time.time()-acn.timebeg
 remaintime = acn.maxtime - elapsedtime



 # Create AMPL object
 ampl = AMPL()

 # Set Knitro as the solver to use
 ampl.setOption('solver', 'knitroampl')
 #ampl.setOption('solver', '/scratch/knitroampl')
 #ampl.setOption('knitroampl_auxfiles', 'rc') #print variable/constraint names inside knitroampl
 #ampl.setOption('TMPDIR','/scratch') # write AMPL temporary files to "/scratch"

 # Set AMPL presolve tolerance and option
 #ampl.setOption('presolve_eps', 1.0e-8)
 ampl.setOption('presolve_eps', 1.0e-0)
 #ampl.setOption('presolve', 0)
 #ampl.setOption('solver_msg', 0)

 #breakexit('premodel')
 # Read the model file                                                  #modelDirectory = argv[2] if argc == 3 else os.path.join('..', 'models')
 usealternate = 0
 useexact = 0
 userelax = 0
 if usealternate == 0:
     if useexact == 0:
         if userelax:
            ampl.read('./constrainedNLP_heur1.mod')
         else:
            ampl.read('./base_norelax.mod')
     else:
         ampl.read('./constrainedNLPexact_heur1.mod')
 else:
     ampl.read('./alternate_heur1.mod')
 #ampl.read('./constrainedNLP.mod.dano')

 Iampl = ampl.getParameter('I')
 Iampl.set(acn.numbuses)
 Gampl = ampl.getParameter('G')
 Gampl.set(acn.numgens)
 Lampl = ampl.getParameter('L')
 Lampl.set(acn.numloads)



 Deltaampl = ampl.getParameter('Delta')
 Deltaampl.set(acn.options['delta'])
 deltarampl = ampl.getParameter('deltar')
 deltarampl.set(acn.options['deltar'])
 NTRcntampl = ampl.getParameter('NTRcnt')
 NTRcntampl.set(acn.nontranscount)
 NTRcntampl_original = ampl.getParameter('NTRcnt_original')
 NTRcntampl_original.set(acn.nontranscount_original)

 TRcntampl = ampl.getParameter('TRcnt')
 TRcntampl_original = ampl.getParameter('TRcnt_original')
 TRcntampl_original.set(acn.transcount_original)
 TRcntampl.set(acn.transcount)
 '''
 '''

 SWSHcntampl = ampl.getParameter('SWSHcnt')
 SWSHcntampl.set(acn.numswitchedshunts)

 numcblocks = ampl.getParameter('numcblocks')
 numcb = {}
 pgcostcblock = ampl.getParameter('pgcostcblock')
 pgcost = {}
 pgmaxcblock = ampl.getParameter('pgmaxcblock')
 pgmax = {}

 Oncost = ampl.getParameter('Oncost')
 oncost = {}
 Sucost = ampl.getParameter('Sucost')
 sucost = {}
 Sdcost = ampl.getParameter('Sdcost')
 sdcost = {}
 maxpgAMPL = ampl.getParameter('maxpg')
 minpgAMPL = ampl.getParameter('minpg')
 maxpgVals = {}
 minpgVals = {}
 maxqgAMPL = ampl.getParameter('maxqg')
 minqgAMPL = ampl.getParameter('minqg')
 maxqgVals = {}
 minqgVals = {}
 genonexistAMPL = ampl.getParameter('genonexist')
 genonexistVals = {}
 genpgexistAMPL = ampl.getParameter('genpgexist')
 genpgexistVals = {}
 genprdAMPL = ampl.getParameter('genprd')
 genpruAMPL = ampl.getParameter('genpru')
 genprdVals = {}
 genpruVals = {}
 gensdqualAMPL = ampl.getParameter('gensdqual')
 gensuqualAMPL = ampl.getParameter('gensuqual')
 gensdqualVals = {}
 gensuqualVals = {}

 #breakexit('base1break1')

 for gencount in range(1, acn.numgens+1):
  ourgen = acn.ourgens[gencount]
  buscount = acn.busnumbertocount[ourgen.i]

  ourbus = acn.ourbuses[buscount]
  numcb[gencount] = ourgen.numcblocks
  oncost[gencount] = ourgen.oncost #well, could have done this to begin with

  sucost[gencount] = ourgen.sucost
  sdcost[gencount] = ourgen.sdcost
  maxpgVals[gencount] = ourgen.maxpg
  minpgVals[gencount] = ourgen.minpg
  maxqgVals[gencount] = ourgen.maxqg
  minqgVals[gencount] = ourgen.minqg
  genonexistVals[gencount] = ourgen.status
  genpgexistVals[gencount] = ourgen.pg
  genprdVals[gencount] = ourgen.prdmax
  genpruVals[gencount] = ourgen.prumax

  gensdqualVals[gencount] = ourgen.sdqual
  gensuqualVals[gencount] = ourgen.suqual
  #print('doing %i\n' %gencount)

  for j in range(1, numcb[gencount]+1):
   pgcost[gencount,j] = ourgen.pgcostcblock[j]
   pgmax[gencount,j] = ourgen.pgmaxcblock[j]

  #costblocks[count] = cblocks['c']

 numcblocks.setValues(numcb)
 #print('numcb',numcb)
 pgcostcblock.setValues(pgcost)
 pgmaxcblock.setValues(pgmax)
 Oncost.setValues(oncost)
 Sucost.setValues(sucost)
 Sdcost.setValues(sdcost)
 #print(sdcost)

 maxpgAMPL.setValues(maxpgVals)
 minpgAMPL.setValues(minpgVals)
 maxqgAMPL.setValues(maxqgVals)
 minqgAMPL.setValues(minqgVals)
 genonexistAMPL.setValues(genonexistVals)
 genpgexistAMPL.setValues(genpgexistVals)
 genprdAMPL.setValues(genprdVals)
 genpruAMPL.setValues(genpruVals)
 gensdqualAMPL.setValues(gensdqualVals)
 gensuqualAMPL.setValues(gensuqualVals)

 # Set bounds on voltage magnitude for AMPL model
 nvloAMPL = ampl.getParameter('nvlo')
 nvhiAMPL = ampl.getParameter('nvhi')
 nvloVals = {}
 nvhiVals = {}
 nvaloAMPL = ampl.getParameter('nvalo')
 nvahiAMPL = ampl.getParameter('nvahi')
 nvaloVals = {}
 nvahiVals = {}
 for i in range(1,1+acn.numbuses):
  ourbus = acn.ourbuses[i]
  nvloVals[i] = ourbus.nvlo
  nvhiVals[i] = ourbus.nvhi
  #nvaloVals[i] = ourbus.nvalo
  #nvahiVals[i] = ourbus.nvahi
  nvaloVals[i] = ourbus.nvalo*math.pi/180
  nvahiVals[i] = ourbus.nvahi*math.pi/180
 nvloAMPL.setValues(nvloVals)
 nvhiAMPL.setValues(nvhiVals)
 nvaloAMPL.setValues(nvaloVals)
 nvahiAMPL.setValues(nvahiVals)

 PcostcblockAMPL = ampl.getParameter('Pcostcblock')
 Pcost = {}
 PmaxcblockAMPL = ampl.getParameter('Pmaxcblock')
 Pmax = {}
 QcostcblockAMPL = ampl.getParameter('Qcostcblock')
 Qcost = {}
 QmaxcblockAMPL = ampl.getParameter('Qmaxcblock')
 Qmax = {}
 if useexact:
     ScostcblockAMPL = ampl.getParameter('Scostcblock')
     Scost = {}
     SmaxcblockAMPL = ampl.getParameter('Smaxcblock')
     Smax = {}
     scblocksAMPL = ampl.getParameter('Snumcblocks')
     scblocksAMPL.set(acn.numscblocks)
 pcblocksAMPL = ampl.getParameter('Pnumcblocks')
 qcblocksAMPL = ampl.getParameter('Qnumcblocks')
 pcblocksAMPL.set(acn.numpcblocks)
 qcblocksAMPL.set(acn.numqcblocks)
 for n in range(1, acn.numpcblocks+1):
  Pcost[n] = acn.pcostcblock[n]
  Pmax[n] = acn.pmaxcblock[n]
  #print(n,': Pcost=',Pcost[n],' Pmax=',Pmax[n],' sbase=',acn.options['sbase'])
 PcostcblockAMPL.setValues(Pcost)
 PmaxcblockAMPL.setValues(Pmax)
 for n in range(1, acn.numqcblocks+1):
  Qcost[n] = acn.qcostcblock[n]
  Qmax[n] = acn.qmaxcblock[n]
  #print(n,': Qcost=',Qcost[n],' Qmax=',Qmax[n],' sbase=',acn.options['sbase'])
 QcostcblockAMPL.setValues(Qcost)
 QmaxcblockAMPL.setValues(Qmax)
 if useexact:
     for n in range(1, acn.numscblocks+1):
         Scost[n] = acn.scostcblock[n]
         Smax[n] = acn.smaxcblock[n]
         #print(n,': Scost=',Scost[n],' Smax=',Smax[n],' sbase=',acn.options['sbase'])
     ScostcblockAMPL.setValues(Scost)
     SmaxcblockAMPL.setValues(Smax)
 #sys.exit('buses')

 nontransatbus = ampl.getSet('nontransatbus')
 nontransatbusVals = {}
 transatbus = ampl.getSet('transatbus')
 transatbusVals = {}
 for i in range(acn.numbuses):
  ourbus = acn.ourbuses[i+1]
  nontrans = []
  for ournontratbus in ourbus.starnontrans.values():
   cnt1 = ournontratbus.count1
   nontrans.append(cnt1)
  nontransatbusVals[i+1] = nontrans
  #print(i,"*",nontrans,ourbus.nontrdegree)
  trans = []
  for ourtransatbus in ourbus.startrans.values():
   cnt1 = ourtransatbus.count1
   trans.append(cnt1)
  transatbusVals[i+1] = trans

 for i in range(1,acn.numbuses+1):
  '''
  log.joint('bus ' + str(i) + " lines: ")
  for j in range(len(nontransatbusVals[i])): log.joint(' ' + str(nontransatbusVals[i][j]))
  log.joint('\n')

  '''
  nontransatbus[i].setValues(nontransatbusVals[i])
  transatbus[i].setValues(transatbusVals[i])

 #breakexit('')


 GffntrAMPL = ampl.getParameter('Gffntr')
 GftntrAMPL = ampl.getParameter('Gftntr')
 BffntrAMPL = ampl.getParameter('Bffntr')
 BftntrAMPL = ampl.getParameter('Bftntr')
 ntrbusiAMPL = ampl.getParameter('ntrbusi')
 ntrbusjAMPL = ampl.getParameter('ntrbusj')
 ntrpAMPL = ampl.getParameter('ntrp')
 ntrratingAMPL = ampl.getParameter('ntrrating0')

 GffntrVals = {}
 GftntrVals = {}
 BffntrVals = {}
 BftntrVals = {}
 ntrbusiVals = {}
 ntrbusjVals = {}
 ntrpVals={}
 ntrratingVals={}

 for i in range(acn.nontranscount):
  GffntrVals[i+1] = acn.nontransATbuses[i].Gff
  GftntrVals[i+1] = acn.nontransATbuses[i].Gft
  BffntrVals[i+1] = acn.nontransATbuses[i].Bff
  BftntrVals[i+1] = acn.nontransATbuses[i].Bft
  ntrbusiVals[i+1] = acn.nontransATbuses[i].busi0+1  #in acn.py busi0 (and busj0) refer to parent nontrans
  ntrbusjVals[i+1] = acn.nontransATbuses[i].busj0+1
  ntrpVals[i+1] = acn.nontransATbuses[i].parentcount1
  ntrratingVals[i+1] = acn.nontransATbuses[i].ournon.rating_a #use rating_c for ctgy
  #print(i+1,ntrpVals[i+1])


 GffntrAMPL.setValues(GffntrVals)
 GftntrAMPL.setValues(GftntrVals)
 BffntrAMPL.setValues(BffntrVals)
 BftntrAMPL.setValues(BftntrVals)
 ntrbusiAMPL.setValues(ntrbusiVals)
 ntrbusjAMPL.setValues(ntrbusjVals)
 ntrpAMPL.setValues(ntrpVals)
 ntrratingAMPL.setValues(ntrratingVals)

 ntrqualswVals={}
 ntrcostswVals={}
 ntrstatVals={}
 ntrqualswAMPL = ampl.getParameter('ntrqualsw')
 ntrcostswAMPL = ampl.getParameter('ntrcostsw')
 ntrstatAMPL = ampl.getParameter('ntrstat')
 for i in range(1,1+acn.nontranscount_original):
  #if (division == 3 and acn.numbuses < 10000 and passnum == 1) or (division == 4 and passnum == 1):
  if (division == 3 and solvenum > 2) or (division == 4 and solvenum > 2):
    ntrqualswVals[i] = acn.ournontrans[i].swqual
  else:
    ntrqualswVals[i] = 0
  ntrcostswVals[i] = acn.ournontrans[i].csw
  ntrstatVals[i] = acn.ournontrans[i].status

 ntrqualswAMPL.setValues(ntrqualswVals)
 ntrcostswAMPL.setValues(ntrcostswVals)
 ntrstatAMPL.setValues(ntrstatVals)

 log.joint(' setting transformer values in AMPL\n')


 trqualswVals={}
 trqualswAMPL = ampl.getParameter('trqualsw')
 trcostswVals={}
 trcostswAMPL = ampl.getParameter('trcostsw')
 trstatVals={}
 trstatAMPL = ampl.getParameter('trstat')
 if useexact:
     trratingorigVals={}
     trratingorigAMPL = ampl.getParameter('trrating_original')
 for i in range(1,1+acn.transcount_original):
  #if (division == 3 and acn.numbuses < 10000 and passnum == 1) or (division == 4 and passnum == 1):
  if (division == 3 and solvenum > 2) or (division == 4 and solvenum > 2):
    trqualswVals[i] = acn.ourtrans[i].swqual
  else:
    trqualswVals[i] = 0
  trcostswVals[i] = acn.ourtrans[i].csw
  trstatVals[i] = acn.ourtrans[i].stat
  if useexact:
      trratingorigVals[i] = acn.ourtrans[i].rating_a #use rating_c for ctgy
 trqualswAMPL.setValues(trqualswVals)
 trcostswAMPL.setValues(trcostswVals)
 trstatAMPL.setValues(trstatVals)
 if useexact:
     trratingorigAMPL.setValues(trratingorigVals)

 trbusiAMPL = ampl.getParameter('trbusi')
 trbusjAMPL = ampl.getParameter('trbusj')
 tr_o_or_rAMPL = ampl.getParameter('tr_o_or_r')
 trpAMPL = ampl.getParameter('trp')
 trratingAMPL = ampl.getParameter('trrating0')

 trbusiVals = {}
 trbusjVals = {}
 tr_o_or_rVals = {}
 trpVals={}
 trratingVals={}

 for i in range(acn.transcount):
  trbusiVals[i+1] = acn.transATbuses[i].busi0+1  #in acn.py, this is busi for the parent transformer
  trbusjVals[i+1] = acn.transATbuses[i].busj0+1
  tr_o_or_rVals[i+1] = 2*(acn.transATbuses[i].original_or_reversed) - 1 # so, = 1 if original, = -1 if reversed
  trpVals[i+1] = acn.transATbuses[i].parentcount1
  trratingVals[i+1] = acn.transATbuses[i].ournon.rating_a #use rating_c for ctgy
  #print(i+1,ntrpVals[i+1])

 trbusiAMPL.setValues(trbusiVals)
 trbusjAMPL.setValues(trbusjVals)
 tr_o_or_rAMPL.setValues(tr_o_or_rVals)
 trpAMPL.setValues(trpVals)
 trratingAMPL.setValues(trratingVals)

 gmfVals={}
 gmfAMPL = ampl.getParameter('gmf')
 bmfVals={}
 bmfAMPL = ampl.getParameter('bmf')

 barxstfVals={}
 barxstfAMPL = ampl.getParameter('barxstf')
 maxxstfVals={}
 maxxstfAMPL = ampl.getParameter('maxxstf')
 minxstfVals={}
 minxstfAMPL = ampl.getParameter('minxstf')

 heurmaxxstfVals={}
 heurmaxxstfAMPL = ampl.getParameter('heurmaxxstf')
 heurminxstfVals={}
 heurminxstfAMPL = ampl.getParameter('heurminxstf')

 brevetaufVals={}
 brevetaufAMPL = ampl.getParameter('brevetauf')
 brevethetafVals={}
 brevethetafAMPL = ampl.getParameter('brevethetaf')
 taustfVals={}
 taustfAMPL = ampl.getParameter('taustf')
 tau0_fixedVals={}
 tau0_fixedAMPL = ampl.getParameter('tau0_fixed')
 thetastfVals={}
 thetastfAMPL = ampl.getParameter('thetastf')
 theta0fVals={}
 theta0fAMPL = ampl.getParameter('theta0f')

 FtauvalidVals={}
 FtauvalidAMPL = ampl.getParameter('Ftauvalid')
 FthetavalidVals={}
 FthetavalidAMPL = ampl.getParameter('Fthetavalid')
 FnuvalidVals={}
 FnuvalidAMPL = ampl.getParameter('Fnuvalid')

 numnuxVals={}
 numnuxAMPL = ampl.getParameter('numnux')

 for cnt1 in range(1,1+acn.transcount_original):
  thistrans = acn.ourtrans[cnt1]
  gmfVals[cnt1] = acn.ourtrans[cnt1].gmf
  bmfVals[cnt1] = acn.ourtrans[cnt1].bmf
  barxstfVals[cnt1] = acn.ourtrans[cnt1].barxstf

  maxxstfVals[cnt1] = acn.ourtrans[cnt1].maxxstf
  minxstfVals[cnt1] = acn.ourtrans[cnt1].minxstf



  brevetaufVals[cnt1] = acn.ourtrans[cnt1].brevetauf
  brevethetafVals[cnt1] = acn.ourtrans[cnt1].brevethetaf
  taustfVals[cnt1] = acn.ourtrans[cnt1].taustf
  tau0_fixedVals[cnt1] = acn.ourtrans[cnt1].tau0f
  thetastfVals[cnt1] = acn.ourtrans[cnt1].thetastf
  theta0fVals[cnt1] = acn.ourtrans[cnt1].theta0f
  FtauvalidVals[cnt1] = acn.ourtrans[cnt1].Ftauvalid
  FthetavalidVals[cnt1] = acn.ourtrans[cnt1].Fthetavalid
  FnuvalidVals[cnt1] = acn.ourtrans[cnt1].Fnuvalid
  numnuxVals[cnt1] = acn.ourtrans[cnt1].numnux

 barxstfAMPL.setValues(barxstfVals)
 gmfAMPL.setValues(gmfVals)
 bmfAMPL.setValues(bmfVals)
 maxxstfAMPL.setValues(maxxstfVals)
 minxstfAMPL.setValues(minxstfVals)
 brevetaufAMPL.setValues(brevetaufVals)
 brevethetafAMPL.setValues(brevethetafVals)
 tau0_fixedAMPL.setValues(tau0_fixedVals)
 taustfAMPL.setValues(taustfVals)
 theta0fAMPL.setValues(theta0fVals)
 thetastfAMPL.setValues(thetastfVals)
 FtauvalidAMPL.setValues(FtauvalidVals)
 FthetavalidAMPL.setValues(FthetavalidVals)
 FnuvalidAMPL.setValues(FnuvalidVals)
 numnuxAMPL.setValues(numnuxVals)



 log.joint(' setting transformer page 54 values in AMPL\n')

 truegAMPL = ampl.getParameter('trueg')
 truegVals = {}
 truebAMPL = ampl.getParameter('trueb')
 truebVals = {}
 for i in range(1,1+acn.transcount_original):

  if acn.ourtrans[i].impedancecorrected==1:
      log.joint('      trans %d numnux %d impedancecorrected %d\n' %(i, acn.ourtrans[i].numnux, acn.ourtrans[i].impedancecorrected))
      if usealternate:
          truegVals[i, -1 + acn.ourtrans[i].minxstf] = 0
          truebVals[i, -1 + acn.ourtrans[i].minxstf] = 0

      for xposition in range(1,1+acn.ourtrans[i].numnux):
          xstf = xposition - 1 + acn.ourtrans[i].minxstf
          truegVals[i,xstf] = acn.ourtrans[i].trueg[xposition]
          truebVals[i,xstf] = acn.ourtrans[i].trueb[xposition]


          #log.joint('trans %d xposition %d xstf %d trueg %g trueb %g\n' %(i, xposition, xstf, truegVals[i,xstf], truebVals[i,xstf]))
  else:
      #log.joint('      trans %d not impedancecorrected numnux %d\n' %(i, acn.ourtrans[i].numnux))
      if usealternate:
          truegVals[i,-1] = acn.ourtrans[i].trueg[1]
          truebVals[i,-1] = acn.ourtrans[i].trueb[1]

      truegVals[i,0] = acn.ourtrans[i].trueg[1]
      truebVals[i,0] = acn.ourtrans[i].trueb[1]

 truegAMPL.setValues(truegVals)
 truebAMPL.setValues(truebVals)

 #breakexit('')


 '''

 '''
 log.joint(' setting generator data\n')
 gencountAMPL = ampl.getParameter('gencount')
 gencountat0AMPL = ampl.getParameter('gencountat0')
 gensatbusAMPL = ampl.getSet('gensatbus')
 gensatbusat0AMPL = ampl.getSet('gensatbusat0')
 gencountVals = {}
 gencountat0Vals = {}
 gensatbusVals = {}
 gensatbusat0Vals = {}
 for count in range(acn.numbuses):
  ourbus = acn.ourbuses[count+1]
  gencountVals[count+1] = ourbus.gencount
  gencountat0Vals[count+1] = ourbus.gencountat0
  gensatthisbus = []
  gensatthisbusat0 = []
  for gc in range(ourbus.gencount):
   gen = ourbus.gens1[gc+1]
   gensatthisbus.append(gen.count)
   if gen.status:
    gensatthisbusat0.append(gen.count)
  gensatbusVals[count+1] = gensatthisbus
  gensatbusat0Vals[count+1] = gensatthisbusat0
 gencountAMPL.setValues(gencountVals)
 gencountat0AMPL.setValues(gencountat0Vals)

 for count in range(1,1+acn.numbuses):
  if gencountVals[count] > 0:
   gensatbusAMPL[count].setValues(gensatbusVals[count])
  if gencountat0Vals[count] > 0:
   gensatbusat0AMPL[count].setValues(gensatbusat0Vals[count])

 loadcountAMPL = ampl.getParameter('loadcount')
 loadisactiveAMPL = ampl.getParameter('loadisactive')
 activeloadcountat0AMPL = ampl.getParameter('activeloadcountat0')
 loadsatbusAMPL = ampl.getSet('loadsatbus')
 activeloadsatbusat0AMPL = ampl.getSet('activeloadsatbusat0')
 loadcountVals = {}
 loadisactiveVals = {}
 activeloadcountat0Vals = {}
 loadsatbusVals = {}
 activeloadsatbusat0Vals = {}
 for count in range(acn.numbuses):
  ourbus = acn.ourbuses[count+1]
  loadcountVals[count+1] = ourbus.loadcount
  activeloadcountat0Vals[count+1] = ourbus.activeloadcountat0
  loadsatthisbus = []
  activeloadsatthisbusat0 = []
  for lc in range(ourbus.loadcount):
   load = ourbus.loads[lc]
   loadsatthisbus.append(load.count)
   loadisactiveVals[load.count] = load.status
   if load.status:
    activeloadsatthisbusat0.append(load.count)
  loadsatbusVals[count+1] = loadsatthisbus
  activeloadsatbusat0Vals[count+1] = activeloadsatthisbusat0
 loadcountAMPL.setValues(loadcountVals)
 loadisactiveAMPL.setValues(loadisactiveVals)
 activeloadcountat0AMPL.setValues(activeloadcountat0Vals)

 for count in range(1,1+acn.numbuses):
  if loadcountVals[count] > 0:
   loadsatbusAMPL[count].setValues(loadsatbusVals[count])
  if activeloadcountat0Vals[count] > 0:
   activeloadsatbusat0AMPL[count].setValues(activeloadsatbusat0Vals[count])

 tminAMPL = ampl.getParameter('tmin')
 tmaxAMPL = ampl.getParameter('tmax')
 pexistAMPL = ampl.getParameter('pexist')
 qexistAMPL = ampl.getParameter('qexist')
 prdAMPL = ampl.getParameter('prd')
 pruAMPL = ampl.getParameter('pru')
 tminVals={}
 tmaxVals={}
 pexistVals = {}
 qexistVals = {}
 prdVals = {}
 pruVals = {}

 loadnumcblocks = ampl.getParameter('loadnumcblocks')
 loadnumcb = {}
 loadcostcblock = ampl.getParameter('loadcostcblock')
 loadcost = {}
 loadmaxcblock = ampl.getParameter('loadmaxcblock')
 loadmax = {}

 for count in range(1,1+acn.numloads):
  load = acn.ourloads[count]
  loadnumcb[count] = load.numcblocks
  tminVals[count] =load.tmin*load.status
  tmaxVals[count] =load.tmax*load.status
  pexistVals[count]= load.plpu
  qexistVals[count]= load.qlpu
  prdVals[count]= load.prd
  pruVals[count]= load.pru

  for j in range(1, loadnumcb[count]+1):
   loadcost[count,j] = load.costcblock[j]
   loadmax[count,j] = load.maxcblock[j]

 loadnumcblocks.setValues(loadnumcb)
 loadcostcblock.setValues(loadcost)
 loadmaxcblock.setValues(loadmax)
 tminAMPL.setValues(tminVals)
 tmaxAMPL.setValues(tmaxVals)
 pexistAMPL.setValues(pexistVals)
 qexistAMPL.setValues(qexistVals)
 prdAMPL.setValues(prdVals)
 pruAMPL.setValues(pruVals)

 fsglpuAMPL = ampl.getParameter('fsglpu')
 fsglpuVals = {}
 fsblpuAMPL = ampl.getParameter('fsblpu')
 fsblpuVals = {}
 vmag = ampl.getVariable('busVmags')
 vmagInit = {}
 vangle = ampl.getVariable('busVangles')
 vangleInit = {}
 for count in range(acn.numbuses):
  ourbus = acn.ourbuses[count+1]
  fsglpuVals[count+1] = ourbus.fsglpu
  fsblpuVals[count+1] = ourbus.fsblpu
  vmagInit[count+1] = 1.0 #flat start or prior value?
  vangleInit[count+1] = 0.0 #flat start or prior value?
 fsglpuAMPL.setValues(fsglpuVals)
 fsblpuAMPL.setValues(fsblpuVals)
 vmag.setValues(vmagInit)
 vangle.setValues(vangleInit)

 activeswshatbusAMPL = ampl.getSet('activeswshatbus')
 activeswshatbusVals = {}
 numactiveswshatbusAMPL = ampl.getParameter('numactiveswshatbus')
 numactiveswshatbusVals = {}
 for count in range(1,1+acn.numbuses):
  ourbus = acn.ourbuses[count]
  activeswatthisbus = []
  numactiveswshatbusVals[count] = ourbus.numactiveswshunts
  for sc in range(ourbus.numactiveswshunts):
   thissw = ourbus.activeswshunts[sc]
   #print('--->',thissw.count)
   activeswatthisbus.append(thissw.count)
  activeswshatbusVals[count] = activeswatthisbus
  #print(count, ourbus.numactiveswshunts, ":::",activeswshatbusVals[count])
  if ourbus.numactiveswshunts > 0:
   activeswshatbusAMPL[count].setValues(activeswshatbusVals[count])
  #breakexit('poooo')

 numactiveswshatbusAMPL.setValues(numactiveswshatbusVals)

 swshactiveAMPL = ampl.getParameter('swshactive')
 swshactiveVals = {}
 for h in range(1,1+acn.numswitchedshunts):
  swshactiveVals[h] = acn.ourswshunts[h].status
 swshactiveAMPL.setValues(swshactiveVals)

 bhasizeAMPL = ampl.getParameter('bhasize')
 bhasizeVals = {}
 for h in range(1,1 + acn.numswitchedshunts):
  oursh = acn.ourswshunts[h]
  bhasizeVals[h] = oursh.len
 bhasizeAMPL.setValues(bhasizeVals)

 bhaAMPL = ampl.getParameter('bha')
 bhaVals = {}
 nhaAMPL = ampl.getParameter('nha')
 nhaVals = {}
 #breakexit('switched')
 for h in range(1,1 + acn.numswitchedshunts):
  oursh = acn.ourswshunts[h]
  #log.joint(' switched shunt ' + str(h) + ' stat ' + str(oursh.status) + '\n')
  for a in range(1,1+oursh.len):
   bhaVals[h,a] = oursh.b[a]
   nhaVals[h,a] = oursh.n[a]
   #log.joint('    h = %d a = %d b = %f n = %d\n' %(h,a,bhaVals[h,a],nhaVals[h,a]))
 bhaAMPL.setValues(bhaVals)
 nhaAMPL.setValues(nhaVals)
 #breakexit('switched')

 forcefixing = 0
 log.joint(' now fixing some integer/binary variables; forcefixing = ' + str(forcefixing) +'\n')

 #breakexit('fixing')
 log.joint(' fixing nontransformer switching variables\n')
 UPntr0Vals ={}
 LOntr0Vals ={}

 fixednontranscount = 0
 for cnt1 in range(1,1+acn.nontranscount_original):
  thisnontrans = acn.ournontrans[cnt1]
  UPntr0Vals[cnt1] = 1.0
  LOntr0Vals[cnt1] = 0.0
  doit = actionable1[thisnontrans.ourfrombus.count] + actionable1[thisnontrans.ourtobus.count]
  fixit = (doit == 0) + forcefixing

  if fixit > 0:

   value = alldata['priorntrON'][cnt1] # acn.ournontrans[cnt1].status

   #log.joint(' nontrans %d being fixed to %d\n' %(cnt1, value))
   UPntr0Vals[cnt1] = value
   LOntr0Vals[cnt1] = value

   fixednontranscount += 1
  else:
   UPntr0Vals[cnt1] = 1
   LOntr0Vals[cnt1] = 0

   log.joint(' nontrans %d loose\n' %(cnt1))
   #breakexit('fixed')

 log.joint(' fixed %d nontrans vars\n' %(fixednontranscount))

 UPntr0AMPL = ampl.getParameter('UPntr0')
 LOntr0AMPL = ampl.getParameter('LOntr0')

 UPntr0AMPL.setValues(UPntr0Vals)
 LOntr0AMPL.setValues(LOntr0Vals)

 log.joint(' fixing transformer switching variables\n')
 UPtr0Vals ={}
 LOtr0Vals ={}

 fixedtranscount = 0
 for cnt1 in range(1,1+acn.transcount_original):
  thistrans = acn.ourtrans[cnt1]

  UPtr0Vals[cnt1] = 1.0
  LOtr0Vals[cnt1] = 0.0

  doit = actionable1[thistrans.ourfrombus.count] + actionable1[thistrans.ourtobus.count]
  fixit = (doit == 0) + forcefixing

  if fixit > 0:
   value =  alldata['priortrON'][cnt1] # acn.ournontrans[cnt1].statusthistrans.stat
   #log.joint(' trans %d being fixed to %d\n' %(cnt1, value))
   UPtr0Vals[cnt1] = value
   LOtr0Vals[cnt1] = value

   fixedtranscount += 1
  else:

   log.joint(' trans %d loose\n' %(cnt1))
   #breakexit('fixed')
 log.joint(' fixed %d trans vars\n' %(fixedtranscount))
 UPtr0AMPL = ampl.getParameter('UPtr0')
 LOtr0AMPL = ampl.getParameter('LOtr0')

 UPtr0AMPL.setValues(UPtr0Vals)
 LOtr0AMPL.setValues(LOtr0Vals)

 #breakexit('f')

 log.joint(' fixing transformer tap variables\n')

 fixed0bin0 = ampl.getSet('fixed0binxstf0')
 fixed0bin0Vals = {}
 numfixedtap0 = 0
 numloosetap0 = 0
 for cnt1 in range(1,1+acn.transcount_original):
  thistrans = acn.ourtrans[cnt1]
  fixedtap = []
  #fix xst0f if neither end is actionable1
  keeploose = actionable1[thistrans.ourfrombus.count] + actionable1[thistrans.ourtobus.count]

  if forcefixing:
   keeploose = 0

  if solvenum > 2 and division <= 2:
   keeploose = 1
  #print("solvenum=",solvenum," division=",division)
  #breakexit("hhhhh")

  if keeploose > 0:
      heurmaxxstfVals[cnt1] = acn.ourtrans[cnt1].maxxstf
      heurminxstfVals[cnt1] = acn.ourtrans[cnt1].minxstf
      #print(' <<<<< ', cnt1, thistrans.ourfrombus.count, thistrans.ourtobus.count, doit, heurmaxxstfVals[cnt1], heurminxstfVals[cnt1])
  else:
      perturbdelta = 0
      heurmaxxstfVals[cnt1] = priorsol.xst0fdict[cnt1] + perturbdelta
      if heurmaxxstfVals[cnt1] > acn.ourtrans[cnt1].maxxstf:
       heurmaxxstfVals[cnt1] = acn.ourtrans[cnt1].maxxstf
      heurminxstfVals[cnt1] = priorsol.xst0fdict[cnt1] - perturbdelta
      if heurminxstfVals[cnt1] < acn.ourtrans[cnt1].minxstf:
       heurminxstfVals[cnt1] = acn.ourtrans[cnt1].minxstf
      if perturbdelta < 2:
       intmax = int(round(heurmaxxstfVals[cnt1]))
       intmin = int(round(heurminxstfVals[cnt1]))
       #log.joint(' fixing for trans %d (%s,%s,%s) heurmin %d heurmax %d truemin %g truemax %g\n' %(cnt1, thistrans.i, thistrans.j, thistrans.ckt, intmin, intmax, acn.ourtrans[cnt1].minxstf, acn.ourtrans[cnt1].maxxstf))

       thisfixedtap = 0
       for m in range(minxstfVals[cnt1], maxxstfVals[cnt1]+1):
        if m < intmin or m > intmax:
         fixedtap.append(m)
         thisfixedtap += 1
        #else:
         #log.joint(' skipped ' + str(m) + '\n')
       #breakexit('fixedtapi')
       numloosetap0 += intmax - intmin
       #if thisfixedtap > 0:
       # log.joint(' fixedtap ' + str(thisfixedtap) + ' bin variables\n')
  fixed0bin0Vals[cnt1] = fixedtap
  numfixedtap0 += len(fixedtap)
  #transatbus[cnt1].setValues(transatbusVals[cnt1])
  fixed0bin0[cnt1].setValues(fixed0bin0Vals[cnt1])
 log.joint(' number of binxstf0 variables fixed at 0: ' + str(numfixedtap0) + '; loose: ' + str(numloosetap0) + '\n'\
)
 heurmaxxstfAMPL.setValues(heurmaxxstfVals)
 heurminxstfAMPL.setValues(heurminxstfVals)

 UPgen0Vals = {}
 LOgen0Vals = {}
 UPgenstartup0Vals = {}
 LOgenstartup0Vals = {}
 UPgenshutdown0Vals = {}
 LOgenshutdown0Vals = {}
 fgen = 0
 f0 = f1 = 0
 log.joint('fixing generator ON variables\n')
 for gencount in range(1, acn.numgens+1):

  thisgen = acn.ourgens[gencount]
  doit = actionable1[thisgen.ourbus.count]
  fixit = (doit == 0) + forcefixing
  UPgen0Vals[gencount] = 1
  LOgen0Vals[gencount] = 0
  UPgenstartup0Vals[gencount] = 1
  LOgenstartup0Vals[gencount] = 0

  UPgenshutdown0Vals[gencount] = 1
  LOgenshutdown0Vals[gencount] = 0

  if fixit > 0:
   fgen += 1
   thevalue = alldata['priorgenON'][gencount]
   #if thevalue == -1:
   # breakexit('wowie ' + str(gencount))
   UPgen0Vals[gencount] = thevalue
   LOgen0Vals[gencount] = thevalue

   UPgenstartup0Vals[gencount] = 0
   LOgenstartup0Vals[gencount] = 0

   UPgenshutdown0Vals[gencount] = 0
   LOgenshutdown0Vals[gencount] = 0


   #if thevalue == 1:
   # log.joint('gencount %d fixed at 1\n' %(gencount))


   #log.joint('fixed gen count %d to %g\n' %(gencount,value)
   f0 += (thevalue == 0)
   f1 += (thevalue == 1)

 UPgen0AMPL = ampl.getParameter('UPgen0')
 LOgen0AMPL = ampl.getParameter('LOgen0')

 UPgen0AMPL.setValues(UPgen0Vals)
 LOgen0AMPL.setValues(LOgen0Vals)

 UPgenstartup0AMPL = ampl.getParameter('UPgenstartup0')
 LOgenstartup0AMPL = ampl.getParameter('LOgenstartup0')

 UPgenstartup0AMPL.setValues(UPgenstartup0Vals)
 LOgenstartup0AMPL.setValues(LOgenstartup0Vals)

 UPgenshutdown0AMPL = ampl.getParameter('UPgenshutdown0')
 LOgenshutdown0AMPL = ampl.getParameter('LOgenshutdown0')


 UPgenshutdown0AMPL.setValues(UPgenshutdown0Vals)
 LOgenshutdown0AMPL.setValues(LOgenshutdown0Vals)

 log.joint('fixed %d generators (0: %d, 1: %d)\n' %(fgen, f0, f1))


 log.joint(' fixing switchable shunt variables\n')
 numfixedxha0 = 0
 numxha0 = 0
 UPxha0Vals = {}
 LOxha0Vals = {}
 for h in range(1,1 + acn.numswitchedshunts):
  oursh = acn.ourswshunts[h]
  doit = actionable1[oursh.ourbus.count]
  fixit = (doit == 0) + forcefixing
  if oursh.status:

    #print(oursh.n)

    if oursh.sizeAh > 0:
     solution = alldata['priorswitchsol'][h]

    for a in range(1,1+oursh.len):
        if oursh.n[a] > 0:
            numxha0 += 1
            UPxha0Vals[oursh.count, a] = oursh.n[a]
            LOxha0Vals[oursh.count, a] = 0

            if fixit>0:
             UPxha0Vals[oursh.count, a] = solution[a]
             LOxha0Vals[oursh.count, a] = solution[a]
             numfixedxha0 += 1
        else:
            break

 UPxha0AMPL = ampl.getParameter('UPxha0')
 LOxha0AMPL = ampl.getParameter('LOxha0')

 UPxha0AMPL.setValues(UPxha0Vals)
 LOxha0AMPL.setValues(LOxha0Vals)

 log.joint('fixed %d switchable shunt variables\n' %(numfixedxha0))
 #breakexit('done fixing')

 #amplstate = 'expand; display {j in 1.._nvars} (_varname[j],_var[j].lb0,_var[j].\
#lb,_var[j].ub,_var[j].ub0);'     #shows full model
 #filename = 'heur1model.out'
 #modelout = ampl.getOutput(amplstate)
 #outfile = open(filename,"w")
 #outfile.write("model = " + str(modelout) + "\n")
 #outfile.close()

 # Set Knitro options
 elapsedtime = time.time()-acn.timebeg
 remaintime = acn.maxtime - elapsedtime
 print('Before solve 1: remain time = ',remaintime)
 #timelimit = min(remaintime,acn.maxtime-75)
 timelimit = min(remaintime-35,acn.maxtime-75)
 timelimit = max(timelimit, 0.001)
 if solvenum <= 2 or division > 2 or 1:
     muruleopt = ' bar_murule=0'  #pred-corr
 else:
     muruleopt = ' bar_murule=1'  #monotone
 logname = ' outname='+'"'+'knitro'+str(solvenum)+'.log'+'"'
 maxtime = ' maxtime_real=' + str(timelimit) + ' mip_maxtime_real=' + str(timelimit) + ' ms_maxtime_real=' + str(timelimit) # just in case we use multi-start
 #maxtime = ' maxtime_real=600 '
 base_options = 'outlev=4 outmode=2 debug=0 presolve_tol=0.5 linsolver=6 feastol=1e-5 ftol=1e-3 scale=0 honorbnds=0 cg_maxit=1 bar_refinement=1 bar_switchrule=0 bar_feasible=1 restarts=3 restarts_maxit=1000 maxit=30000'  # use bar_initpt=2 bar_directinterval=0 ?
 just_feasible = 1 #ignoring optimality works quite well
 if just_feasible:
     base_options += ' opttol=1.0e30 opttol_abs=1e30'
 else:
     base_options += ' opttol=1.0e-3'
 if (acn.numbuses < 10000):
     base_options += ' bar_maxcrossit=1'
 mip_options = ' mip_multistart=1 mip_terminate=0 mip_outinterval=1 mip_outlevel=1 mip_debug=1 mip_outsub=2 mip_nodealg=1 mip_intvar_strategy=0'
 if (1):
     intvar_strategy = ' mip_intvar_strategy=0'  # handle integer vars directly
 else:
     intvar_strategy = ' mip_intvar_strategy=2'  # convert binary vars to complementarity constraints
 #if (acn.numbuses < 100 or passnum == 2):
 if (0):
   relax = ' relax=0'
   log.joint(' using relax = 0 ...\n')
 else:
   relax = ' relax=1'
   log.joint(' using relax = 1 ...\n')
 options = base_options + muruleopt + mip_options + maxtime + relax + logname + intvar_strategy
 #print "options=",options
 ampl.setOption('knitro_options', options)

 # Solve

 usescript = 1
 if usescript:

  try:
      if solvenum == 0 or solvenum == 3:
         if passnum == 1:
           log.joint(' !!! solvenum=0/3, fixed solve (solve_script)\n')
           ampl.eval('include solve_script;')
         else: #passnum=2
           log.joint(' !!! solvenum=0/3, unfixed solve (solve_script12)\n')
           ampl.eval('include solve_script12;')
      elif solvenum == 1 or solvenum == 4:
         log.joint(' !!! solvenum=1/4, fixed swsh solve (solve_script11)\n')
         ampl.eval('include solve_script11;')
      elif solvenum == 2 or solvenum == 5:
         log.joint(' !!! solvenum=2/5, unfixed solve (solve_script12)\n')
         ampl.eval('include solve_script12;')
      else:
         ampl.eval('include solve_script;')
  except:
      log.joint("Error calling ampl solve_script!\n")
      var = traceback.format_exc()
      traceback.print_exc()
      acn_evaluation2.print_alert(var, raise_exception=False)
      pass

 else: #old way

  try:

   totalcost, maxviol, timeneeded = acnsolveit(acn, log, ampl, 'ourobjective', 'first')
   '''
   time1 = time.time()
   ampl.solve()

   time2 = time.time()
   totalcost = ampl.getObjective('ourobjective')
   log.joint(" >>> base1 objective value: %.6e in time %g\n" %(totalcost.value(), time2-time1))
   '''
   acnexposevars(acn, log, ampl, truegVals, truebVals, 'first')

   time2 = time.time()
  except:
   log.joint("Error calling ampl.solve!\n")
   var = traceback.format_exc()
   traceback.print_exc()
   acn_evaluation2.print_alert(var, raise_exception=False)
   pass

 #breakexit('base1break2')
 #amplstate = 'display binxstf0, xstf0, gf0, bf0;'
 #filename = 'structure1.out'
 #modelout = ampl.getOutput(amplstate)
 #outfile = open(filename,"w")
 #outfile.write("objective terms: \n" + str(modelout) + "\n")
 #outfile.close()

 # Re-solve with integer vars fixed
 elapsedtime = time.time()-acn.timebeg
 remaintime = acn.maxtime - elapsedtime
 print('Before solve 2: remain time = ',remaintime, ' elapsed time= ',elapsedtime)
 log.joint('Before solve 2: remain time = ' + str(remaintime) + ' elapsed time= ' + str(elapsedtime) +'\n')
 if remaintime > 30:
     ampl.setOption('presolve_eps', 1.0e-6)
     resolve_fixed(acn, ampl, usescript, solvenum, userelax)

 #amplstate = 'display binxstf0, xstf0, gf0, bf0;'
 #filename = 'finalstructure.out'
 #modelout = ampl.getOutput(amplstate)
 #outfile = open(filename,"w")
 #outfile.write("objective terms: \n" + str(modelout) + "\n")
 #outfile.close()

 elapsedtime = time.time()-acn.timebeg
 remaintime = acn.maxtime - elapsedtime
 print('After solve 2: remain time = ',remaintime, ' elapsed time= ',elapsedtime)
 log.joint('After solve 2: remain time = ' + str(remaintime) + ' elapsed time= ' + str(elapsedtime) +'\n')

 # Write some AMPL solution files to file
 ##amplstate = 'display gencost0, genon0, genstartup0, genshutdown0, genP0, genQ0, minpg, maxpg, minqg, maxqg;'
 #amplstate = 'display sum{g in 1..G} Delta*(-gencost0[g] - Oncost[g]*genon0[g]) - sum{g in 1..G} (Sucost[g]*genstartup0[g] + Sdcost[g]*genshutdown0[g]);'
 ##amplstate += 'display -sum{i in 1..I} Delta*buscost0[i];'
 #amplstate += 'display -sum{i in 1..I} Delta*(sum{n in 1..Pnumcblocks} Pcostcblock[n]*(Pplusblock0[i,n] + Pminusblock0[i,n]));'
 #amplstate += 'display -sum{i in 1..I} Delta*(sum{n in 1..Qnumcblocks} Qcostcblock[n]*(Qplusblock0[i,n] + Qminusblock0[i,n]));'
 #amplstate += 'display sum{l in 1..L: loadisactive[l] > 0} Delta*loadcost0[l];'
 #amplstate += 'display -sum{e in 1..NTRcnt_original: ntrqualsw[e]>0 and ntrstat[e]=0} ntrcostsw[e]*ntrON0[e];'
 #amplstate += 'display -sum{e in 1..NTRcnt_original: ntrqualsw[e]>0 and ntrstat[e]>0} ntrcostsw[e]*(ntrstat[e]-ntrON0[e]);'
 #amplstate += 'display -sum{f in 1..TRcnt_original: trqualsw[f]>0 and trstat[f]=0} trcostsw[f]*trON0[f];'
 #amplstate += 'display -sum{f in 1..TRcnt_original: trqualsw[f]>0 and trstat[f]>0} trcostsw[f]*(trstat[f]-trON0[f]);'
 #filename = 'amplsolheur.out'
 #modelout = ampl.getOutput(amplstate)
 #outfile = open(filename,"w")
 #outfile.write("objective terms: \n" + str(modelout) + "\n")
 #outfile.close()


 # Extract base solution from AMPL (move to function?)
 sbase = acn.options['sbase'] #acn.raw.case_identification.sbase #needed?
 basebusVmags = numpy.zeros(acn.numbuses)
 v = ampl.getVariable('busVmags')
 df = v.getValues()
 vals = df.toDict()
 for i in range(acn.numbuses):
    basebusVmags[i] = vals[i+1]

 basebusVangles = numpy.zeros(acn.numbuses)
 v = ampl.getVariable('busVangles')
 df = v.getValues()
 vals = df.toDict()
 for i in range(acn.numbuses):
    basebusVangles[i] = vals[i+1]

 clearedloads = numpy.zeros(1+acn.numloads)
 v = ampl.getVariable('t0')
 dv = v.getValues()
 vals = dv.toDict()

 for i in range(acn.numloads):
     clearedloads[i+1] = vals[i+1]

 for count in range(acn.numbuses):
  ourbus = acn.ourbuses[count+1]
  for lc in range(ourbus.loadcount):
   load = ourbus.loads[lc]
   if load.status:
       i = load.count
       #log.joint(' load %d t0 %g\n' %(i, clearedloads[i]))

 ploads = numpy.zeros(acn.numloads+1)
 qloads = numpy.zeros(acn.numloads+1)
 u = ampl.getVariable('loadP0')
 du = u.getValues()
 uals = du.toDict()
 w = ampl.getVariable('loadQ0')
 dw = w.getValues()
 wals = dw.toDict()
 for i in range(acn.numloads):
    ploads[i+1] = uals[i+1]
    qloads[i+1] = wals[i+1]
    count += 1

 '''

 #report the long way
 for count in range(acn.numbuses):
  ourbus = acn.ourbuses[count+1]
  for lc in range(ourbus.loadcount):
   load = ourbus.loads[lc]
   if load.status:
       i = load.count
       log.joint(' load %d t0 %g p %g q %g\n' %(i, clearedloads[i], ploads[i], qloads[i]))
 '''

 basegenP = numpy.zeros(acn.numgens)
 v = ampl.getVariable('genP0')
 df = v.getValues()
 vals = df.toDict()
 for i in range(acn.numgens):
    basegenP[i] = vals[i+1]

 basegencost = numpy.zeros(acn.numgens)
 v = ampl.getVariable('gencost0')
 df = v.getValues()
 vals = df.toDict()
 #log.joint(' in acnsolverAMPL basegencosts\n')
 sumgenbase = 0
 for i in range(acn.numgens):
    basegencost[i] = vals[i+1]
    sumgenbase += vals[i+1]
    #log.joint(' gen %d cost %g sumgenbase %.8e\n' %(i+1, vals[i+1], sumgenbase))

 basegenQ = numpy.zeros(acn.numgens)
 v = ampl.getVariable('genQ0')
 df = v.getValues()
 vals = df.toDict()
 for i in range(acn.numgens):
    basegenQ[i] = vals[i+1]

 basegenON = numpy.zeros(acn.numgens,dtype=int)
 basegensu = numpy.zeros(acn.numgens,dtype=int)
 basegensd = numpy.zeros(acn.numgens,dtype=int)
 v = ampl.getVariable('genon0')
 df = v.getValues()
 valson = df.toDict()

 #log.joint(' in acnsolverAMPL basegenONcosts\n')
 sumgenon = 0
 for i in range(acn.numgens):
    basegenON[i] = valson[i+1]
    thiscost = valson[i+1]*oncost[i+1]
    sumgenon += thiscost
    #log.joint(' gen %d on %g cost %g sumgenon %.8e\n' %(i+1, valson[i+1], thiscost, sumgenon))

 #log.joint(' in acnsolverAMPL basegensucosts\n')
 vsu = ampl.getVariable('genstartup0')
 df = vsu.getValues()
 valssu = df.toDict()
 sumgensu = 0
 for i in range(acn.numgens):
    basegensu[i] = valssu[i+1]
    thiscost = valssu[i+1]*sucost[i+1]
    sumgensu += thiscost
    #log.joint(' gen %d on %g cost %g sumgensu %.8e\n' %(i+1, valssu[i+1], thiscost, sumgensu))

 vsd = ampl.getVariable('genshutdown0')
 df = vsd.getValues()
 valssd = df.toDict()
 for i in range(acn.numgens):
    basegensd[i] = valssd[i+1]

 sumgensd = 0
 for i in range(acn.numgens):
    basegensd[i] = valssd[i+1]
    thiscost = valssd[i+1]*sdcost[i+1]
    sumgensd += thiscost
    #log.joint(' gen %d on %g cost %g sumgensd %.8e\n' %(i+1, valssd[i+1], thiscost, sumgensd))

 '''
 log.joint(' gencost: %.8e\n' %(sumgenbase))
 log.joint(' gencost + oncost: %.8e\n' %(sumgenbase + sumgenon))
 log.joint(' gencost + oncost + sucost: %.8e\n' %(sumgenbase + sumgenon + sumgensu))
 log.joint(' gencost + oncost + sucost + sdcost: %.8e\n' %(sumgenbase + sumgenon + sumgensu + sumgensd))
 '''

 #relaxPplus0 = numpy.zeros(acn.numbuses)
 #v = ampl.getVariable('relaxPplus0')
 #df = v.getValues()
 #vals = df.toDict()
 #for i in range(acn.numbuses):
 #   relaxPplus0[i] = vals[i+1]
 #print(relaxPplus0)
 #relaxPminus0 = numpy.zeros(acn.numbuses)
 #v = ampl.getVariable('relaxPminus0')
 #df = v.getValues()
 #vals = df.toDict()
 #for i in range(acn.numbuses):
 #   relaxPminus0[i] = vals[i+1]
 #print(relaxPminus0)

 '''
 log.joint('active power imbalances:\n')
 for i in range(acn.numbuses):
    log.joint(' bus %d Pplus %g Pminus %g sum %g\n' %(i+1, relaxPplus0[i], relaxPminus0[i], relaxPplus0[i] + relaxPminus0[i]))
    #log.joint(' bus %d Pplus %g ' Pminus %g total %g\n' %(i+1, relaxPplus0[i], relaxPminus0[i], relaxPlus0[i] + relaxPminus0[i]))
    #log.joint(' bus ' + str(i+1)+ ' Pplus ' + str(relaxPplus0[i]) + ' Pminus ' + str(relaxPminus0[i]) + ' total ' + str(relaxPplus0[i] + relaxPminus0[i]) + '\n')
 '''

 basentrP = numpy.zeros(acn.nontranscount)
 v = ampl.getVariable('ntrP0')
 df = v.getValues()
 vals = df.toDict()
 for i in range(acn.nontranscount):
    basentrP[i] = vals[i+1]

 '''
 for i in range(1,1+acn.numbuses):
  ourbus = acn.ourbuses[i]
  log.joint('bus ' + str(i) + " p ntrps:\n")
  for j in range(len(nontransatbusVals[i])):
      k = nontransatbusVals[i][j]
      pow = basentrP[k-1]
      log.joint(' ' + str(k) + ' p ' + str(pow) + '\n')
 '''

 basentrQ = numpy.zeros(acn.nontranscount)
 v = ampl.getVariable('ntrQ0')
 df = v.getValues()
 vals = df.toDict()
 for i in range(acn.nontranscount):
    basentrQ[i] = vals[i+1]

 '''
 for i in range(1,1+acn.numbuses):
  ourbus = acn.ourbuses[i]
  log.joint('bus ' + str(i) + " q ntrps:\n")
  for j in range(len(nontransatbusVals[i])):
      k = nontransatbusVals[i][j]
      qow = basentrQ[k-1]
      log.joint(' ' + str(k) + ' q ' + str(qow) + '\n')
 '''



 #breakexit('got ntr flows')
 basentrON = numpy.zeros(acn.nontranscount_original,dtype=int)
 v = ampl.getVariable('ntrON0')
 df = v.getValues()
 vals = df.toDict()
 for i in range(acn.nontranscount_original):
    basentrON[i] = int(vals[i+1])

 basetrON = numpy.zeros(acn.transcount_original,dtype=int)
 basetaps = numpy.zeros(acn.transcount_original,dtype=int)
 v = ampl.getVariable('trON0')
 df = v.getValues()
 vals = df.toDict()
 for i in range(acn.transcount_original):
     basetrON[i] = int(vals[i+1])
 alldata['basetrON'] = basetrON

 v = ampl.getVariable('xstf0')
 df = v.getValues()
 vals = df.toDict()
 for i in range(acn.transcount_original):
     basetaps[i] = int(vals[i+1])

 basetrP = numpy.zeros(acn.transcount)
 v = ampl.getVariable('trP0')
 df = v.getValues()
 vals = df.toDict()
 for i in range(acn.transcount):
    basetrP[i] = vals[i+1]

 mykeybus = alldata['mykeybus']
 mykeytrdual = alldata['mykeytrdual']

 acnexposevars(acn, log, ampl, truegVals, truebVals, 'last')


 basetrQ = numpy.zeros(acn.transcount)
 v = ampl.getVariable('trQ0')
 df = v.getValues()
 vals = df.toDict()
 for i in range(acn.transcount):
    basetrQ[i] = vals[i+1]


 log.joint('Retrieve xha0: numxha0 = %d (%d)\n' %(numxha0, acn.numswitchedshunts))
 basesw={}
 j=0
 if numxha0 > 0:
    v = ampl.getVariable('xha0')
    df = v.getValues()
    vals = df.toDict()
    for (h,a), step in vals.items():
      #print(h,a,step)
      basesw[h,a] = int(step)
      acn.basesol_sw[solvenum,j] = int(step)
      j=j+1
    #lsw = len(basesw)
    #print('lsw=',lsw,' numxha0=',numxha0)
    #print('basesw=',basesw)
    #breakexit("sw")

 elapsedtime = time.time()-acn.timebeg
 remaintime = acn.maxtime - elapsedtime
 print('Before sol eval: remain time = ',remaintime)

 # Create base case solution data structure and file
 heurfilename = "base_" + str(solvenum) + ".txt"
 heurfilename_noext = "base_" + str(solvenum)
 log.joint("Write base case solution\n")
 s1 = acnSolver()
 s1.data.read(raw, sup, con)
 s1.write_sol_base_only('.',heurfilename_noext,acn,basebusVmags,basebusVangles,clearedloads,basegenP,basegenQ,basegenON,basentrON,basetrON,basetaps,basesw)
 acn.basesol_busVmags[solvenum,:] = basebusVmags[:]
 acn.basesol_busVangles[solvenum,:] = basebusVangles[:]
 acn.basesol_clearedloads[solvenum,:] = clearedloads[:]
 acn.basesol_ploads[solvenum,:] = ploads[:]
 acn.basesol_genP[solvenum,:] = basegenP[:]
 acn.basesol_genQ[solvenum,:] = basegenQ[:]
 acn.basesol_genON[solvenum,:] = basegenON[:]
 acn.basesol_gensu[solvenum,:] = basegensu[:]
 acn.basesol_gensd[solvenum,:] = basegensd[:]
 acn.basesol_ntrON[solvenum,:] = basentrON[:]
 acn.basesol_trON[solvenum,:] = basetrON[:]
 acn.basesol_taps[solvenum,:] = basetaps[:]
 #acn.basesol_sw filled in above; store the initial one
 if solvenum == 0:
   acn.basesw0 = basesw


#############

def acnexposevars_heur(acn, log, ampl, truegVals, truebVals, name):

 alldata = acn.alldata

 mykeybus = alldata['mykeybus']
 mykeytrdual = alldata['mykeytrdual']
 log.joint('exposing vars, ' + str(name)+ ', mykeybus = %d, mykeytrdual = %d\n' %(mykeybus, mykeytrdual))
 if mykeybus <= 0 and mykeytrdual <= 0:
     log.joint('nothing to expose\n')
     return


 Gf0a = ampl.getVariable('Gf0')
 Gf0f = Gf0a.getValues()
 Gf0vals = Gf0f.toDict()

 gf0a = ampl.getVariable('gf0')
 gf0f = gf0a.getValues()
 gf0vals = gf0f.toDict()

 bf0a = ampl.getVariable('bf0')
 bf0f = bf0a.getValues()
 bf0vals = bf0f.toDict()

 phasea = ampl.getVariable('thetaf0')
 phasef = phasea.getValues()
 phasevals = phasef.toDict()

 gf0a = ampl.getVariable('gf0')
 gf0f = gf0a.getValues()
 gf0vals = gf0f.toDict()

 bina = ampl.getVariable('binxstf0')
 binf = bina.getValues()
 binvals = binf.toDict()

 # Transformer tap settings
 vx = ampl.getVariable('xstf0')
 fx = vx.getValues()
 xvals = fx.toDict()

 basetrP = numpy.zeros(acn.transcount)
 v = ampl.getVariable('trP0')
 df = v.getValues()
 vals = df.toDict()
 for i in range(acn.transcount):
    basetrP[i] = vals[i+1]

 basetrQ = numpy.zeros(acn.transcount)
 v = ampl.getVariable('trQ0')
 df = v.getValues()
 vals = df.toDict()
 for i in range(acn.transcount):
    basetrQ[i] = vals[i+1]

 basebusVmags = numpy.zeros(acn.numbuses)
 v = ampl.getVariable('busVmags')
 df = v.getValues()
 vals = df.toDict()
 for i in range(acn.numbuses):
    basebusVmags[i] = vals[i+1]

 basebusVangles = numpy.zeros(acn.numbuses)
 v = ampl.getVariable('busVangles')
 df = v.getValues()
 vals = df.toDict()
 for i in range(acn.numbuses):
    basebusVangles[i] = vals[i+1]


 basetrON = numpy.zeros(acn.transcount_original,dtype=int)
 v = ampl.getVariable('trON0')
 df = v.getValues()
 vals = df.toDict()
 for i in range(acn.transcount_original):
     basetrON[i] = int(vals[i+1])

 relQmv = ampl.getVariable('relaxQminus0')
 relQmf = relQmv.getValues()
 relQmvals = relQmf.toDict()

 Qswshv = ampl.getVariable('Qswsh')
 Qswshf = Qswshv.getValues()
 Qswshvals = Qswshf.toDict()

 if mykeybus >= 1 and mykeybus <= acn.numbuses:
    thisbus = acn.ourbuses[mykeybus]
    log.joint(' Qwsh at bus1 ' + str(mykeybus) + ' ' + str(Qswshvals[mykeybus]) + '\n')
    for trATbus in thisbus.startrans.values():
        cnt1 = trATbus.count1
        #print(cnt1, acn.transcount)
        log.joint(' !!! ' + str(cnt1) + ' pf ' + str(basetrP[cnt1-1]) + ' qf ' + str(basetrQ[cnt1-1]) + ' parent ' + str(trATbus.parentcount1) + '\n')
        #breakexit('1')
        trbusi = trATbus.busi0+1
        #breakexit('2')
        trbusj = trATbus.busj0+1
        #breakexit('3')
        log.joint('    ' + ' busi1 ' + str(trbusi) + ' busj1 ' + str(trbusj) + ' tron ' + str(basetrON[trATbus.parentcount1-1])+ ' o_or_r ' + str(trATbus.original_or_reversed) + '\n')
        #breakexit('')
        if cnt1 == mykeytrdual:
            mykeytr1 = trATbus.parentcount1
            vGf0f = Gf0vals[cnt1]
            vgf0f = gf0vals[mykeytr1]
            vbf0f = bf0vals[mykeytr1]
            i = trbusi
            j = trbusj
            log.joint('   angle at ' + str(i) + ' ' + str(basebusVangles[i-1]) )
            log.joint(' angle at ' + str(j) + ' ' + str(basebusVangles[j-1]) )
            log.joint(' phase ' + str(phasevals[mykeytr1]) + '\n')
            log.joint('   tap ' + str(xvals[mykeytr1]) + '\n')
            minx = acn.ourtrans[mykeytr1].minxstf #minxstfVals[mykeytr1]
            maxx = acn.ourtrans[mykeytr1].maxxstf #maxxstfVals[mykeytr1]
            log.joint('   min ' + str(minx) )
            log.joint('   max ' + str(maxx) + '\n')

            for m in range(acn.ourtrans[mykeytr1].minxstf, 1 + acn.ourtrans[mykeytr1].maxxstf):

                if binvals[mykeytr1,m] > 1e-5:
                    trueg = truegVals[mykeytr1,m]
                    trueb = truebVals[mykeytr1,m]
                    log.joint('   m  %d bin %g trueg %g trueb %g\n' %(m, binvals[mykeytr1,m], trueg, trueb))
            log.joint('    G %g g %g b %g\n' %(vGf0f, vgf0f, vbf0f))
            #breakexit('minmax')


 log.joint('exposed vars\n')
