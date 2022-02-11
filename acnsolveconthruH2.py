#!/usr/bin/python

import csv
import sys, os, shutil
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
from acnextra import *
#from acnsolverAMPL2.py import resolve_fixed2

def trace(frame, event, arg):
    print("%s, %s:%d" % (event, frame.f_code.co_filename, frame.f_lineno))
    return trace


# the bulk of the code for setting up and solving individual contingencies is in acnsolverAMPLcon
def acnsolveconthruH2(alldata,acn, log, passnum, ctgcount, ctgidx, conlabel, ctg, basebusVmags, basebusVangles, basentrON, basetrON, baseloadP, basegenP, basegenON, basegenSU, basegenSD, basesw, basexstf, genout_index, ntrout_index, trout_index, genoutIndex, genoutID, lineout_i, lineout_j, lineout_ckt, exceptions, numexceptions):

 log.joint("\n***\n Entering acnsolveconthruH2 for contingency: %s with %d exceptions\n" %(conlabel,numexceptions))

 numbuses = acn.numbuses
 #put exceptions in actionable
 division = int(acn.division)

 inexceptions1 = acn.actionable1
 inexceptions1.fill(0)

 for j in range(1,1+numbuses):
     inexceptions1[j] = acn.ourbuses[j].count in exceptions.values()

 #Whether or not to use base solution as initial point
 initwithbase = 1

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

 # Read the model file
 #modelDirectory = argv[2] if argc == 3 else os.path.join('..', 'models')
 #ampl.read('./constrainedNLP2.mod')
 #ampl.read('./constrainedNLP2exact.mod')
 ampl.read('./constrainedNLP2exact_heur2.mod')

 genoutAmpl = ampl.getParameter('genout_index')
 genoutAmpl.set(genout_index)
 #genoutAmpl.set(-1)
 ntroutAmpl = ampl.getParameter('ntrout_index')
 ntroutAmpl.set(ntrout_index)
 #ntroutAmpl.set(-1)
 troutAmpl = ampl.getParameter('trout_index')
 troutAmpl.set(trout_index)
 #troutAmpl.set(-1)
 #breakexit('genout_index')

 Iampl = ampl.getParameter('I')
 Iampl.set(acn.numbuses)
 Gampl = ampl.getParameter('G')
 Gampl.set(acn.numgens)
 Lampl = ampl.getParameter('L')
 Lampl.set(acn.numloads)


 Deltactampl = ampl.getParameter('Deltact')
 Deltactampl.set(acn.options['deltactg'])
 deltarctgampl = ampl.getParameter('deltarctg')
 deltarctgampl.set(acn.options['deltarctg'])
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
 #genonexistVals = {}
 genpgexistAMPL = ampl.getParameter('genpgexist')
 #genpgexistVals = {}
 genprdAMPL = ampl.getParameter('genprdctg')
 genpruAMPL = ampl.getParameter('genpructg')
 genprdVals = {}
 genpruVals = {}
 gensdqualAMPL = ampl.getParameter('gensdqualctg')
 gensuqualAMPL = ampl.getParameter('gensuqualctg')
 gensdqualVals = {}
 gensuqualVals = {}
 gensd0AMPL = ampl.getParameter('genshutdown0')
 gensu0AMPL = ampl.getParameter('genstartup0')
 if initwithbase:
     genonAMPL = ampl.getVariable('genon')
     genonAMPL.setValues(basegenON)
     genPAMPL = ampl.getVariable('genP')
     genPAMPL.setValues(basegenP)

 for gencount in range(1, acn.numgens+1):
  ourgen = acn.ourgens[gencount]
  buscount = acn.busnumbertocount[ourgen.i]

  ourbus = acn.ourbuses[buscount]
  numcb[gencount] = ourgen.numcblocks
  oncost[gencount] = ourgen.oncost #well, could have done this to begin with

  sucost[gencount] = ourgen.sucost
  sdcost[gencount] = ourgen.sdcost
  if ourgen.i == genoutIndex and ourgen.id == genoutID:
      # Generator is out.  How to force genon=0? ---> below
      maxpgVals[gencount] = 0
      minpgVals[gencount] = 0
      maxqgVals[gencount] = 0
      minqgVals[gencount] = 0
  else:
      maxpgVals[gencount] = ourgen.maxpg
      minpgVals[gencount] = ourgen.minpg
      maxqgVals[gencount] = ourgen.maxqg
      minqgVals[gencount] = ourgen.minqg
  #genonexistVals[gencount] = ourgen.status
  #genpgexistVals[gencount] = ourgen.pg
  genprdVals[gencount] = ourgen.prdmaxctg
  genpruVals[gencount] = ourgen.prumaxctg

  gensdqualVals[gencount] = ourgen.sdqualctg
  gensuqualVals[gencount] = ourgen.suqualctg
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
 #genonexistAMPL.setValues(genonexistVals)
 genonexistAMPL.setValues(basegenON)
 #genpgexistAMPL.setValues(genpgexistVals)
 genpgexistAMPL.setValues(basegenP)
 genprdAMPL.setValues(genprdVals)
 genpruAMPL.setValues(genpruVals)
 gensdqualAMPL.setValues(gensdqualVals)
 gensuqualAMPL.setValues(gensuqualVals)
 gensd0AMPL.setValues(basegenSD)
 gensu0AMPL.setValues(basegenSU)

 # Set bounds on voltage magnitude for AMPL model
 evloAMPL = ampl.getParameter('evlo')
 evhiAMPL = ampl.getParameter('evhi')
 evloVals = {}
 evhiVals = {}
 nvaloAMPL = ampl.getParameter('nvalo')
 nvahiAMPL = ampl.getParameter('nvahi')
 nvaloVals = {}
 nvahiVals = {}
 for i in range(1,1+acn.numbuses):
  ourbus = acn.ourbuses[i]
  evloVals[i] = ourbus.evlo
  evhiVals[i] = ourbus.evhi
  #nvaloVals[i] = ourbus.nvalo
  #nvahiVals[i] = ourbus.nvahi
  nvaloVals[i] = ourbus.nvalo*math.pi/180
  nvahiVals[i] = ourbus.nvahi*math.pi/180
 evloAMPL.setValues(evloVals)
 evhiAMPL.setValues(evhiVals)
 nvaloAMPL.setValues(nvaloVals)
 nvahiAMPL.setValues(nvahiVals)

 # Set cost blocks
 PcostcblockAMPL = ampl.getParameter('Pcostcblock')
 Pcost = {}
 PmaxcblockAMPL = ampl.getParameter('Pmaxcblock')
 Pmax = {}
 QcostcblockAMPL = ampl.getParameter('Qcostcblock')
 Qcost = {}
 QmaxcblockAMPL = ampl.getParameter('Qmaxcblock')
 Qmax = {}
 ScostcblockAMPL = ampl.getParameter('Scostcblock')
 Scost = {}
 SmaxcblockAMPL = ampl.getParameter('Smaxcblock')
 Smax = {}
 pcblocksAMPL = ampl.getParameter('Pnumcblocks')
 qcblocksAMPL = ampl.getParameter('Qnumcblocks')
 scblocksAMPL = ampl.getParameter('Snumcblocks')
 pcblocksAMPL.set(acn.numpcblocks)
 qcblocksAMPL.set(acn.numqcblocks)
 scblocksAMPL.set(acn.numscblocks)
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
 for n in range(1, acn.numscblocks+1):
  Scost[n] = acn.scostcblock[n]
  Smax[n] = acn.smaxcblock[n]
  #print(n,': Scost=',Scost[n],' Smax=',Smax[n],' sbase=',acn.options['sbase'])
 ScostcblockAMPL.setValues(Scost)
 SmaxcblockAMPL.setValues(Smax)


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

 GffntrVals = {}
 GftntrVals = {}
 BffntrVals = {}
 BftntrVals = {}
 ntrbusiVals = {}
 ntrbusjVals = {}
 ntrpVals={}

 for i in range(acn.nontranscount):
  GffntrVals[i+1] = acn.nontransATbuses[i].Gff
  GftntrVals[i+1] = acn.nontransATbuses[i].Gft
  BffntrVals[i+1] = acn.nontransATbuses[i].Bff
  BftntrVals[i+1] = acn.nontransATbuses[i].Bft
  ntrbusiVals[i+1] = acn.nontransATbuses[i].busi0+1 #in acn.py busi0 (and busj0) refer to parent nontrans
  ntrbusjVals[i+1] = acn.nontransATbuses[i].busj0+1
  ntrpVals[i+1] = acn.nontransATbuses[i].parentcount1
  #print(i+1,ntrpVals[i+1])


 GffntrAMPL.setValues(GffntrVals)
 GftntrAMPL.setValues(GftntrVals)
 BffntrAMPL.setValues(BffntrVals)
 BftntrAMPL.setValues(BftntrVals)
 ntrbusiAMPL.setValues(ntrbusiVals)
 ntrbusjAMPL.setValues(ntrbusjVals)
 ntrpAMPL.setValues(ntrpVals)

 ntrqualswVals={}
 ntrcostswVals={}
 #ntrstatVals={}
 ntrratingVals={}
 ntrqualswAMPL = ampl.getParameter('ntrqualsw')
 ntrcostswAMPL = ampl.getParameter('ntrcostsw')
 ntrstatAMPL = ampl.getParameter('ntrstat')
 ntrratingAMPL = ampl.getParameter('ntrratingct')
 for i in range(1,1+acn.nontranscount_original):
  if passnum == 2 or passnum == 4 or passnum == 6:
     ntrqualswVals[i] = acn.ournontrans[i].swqual
  else:
     ntrqualswVals[i] = 0
  ntrcostswVals[i] = acn.ournontrans[i].csw
  #ntrstatVals[i] = acn.ournontrans[i].status
  ntrratingVals[i] = acn.ournontrans[i].rating_c #use rating_c for ctgy
 ntrqualswAMPL.setValues(ntrqualswVals)
 ntrcostswAMPL.setValues(ntrcostswVals)
 #ntrstatAMPL.setValues(ntrstatVals)
 ntrstatAMPL.setValues(basentrON)
 ntrratingAMPL.setValues(ntrratingVals)
 if initwithbase:
     ntrONAMPL = ampl.getVariable('ntrON')
     ntrONAMPL.setValues(basentrON)

 log.joint(' setting transformer values in AMPL\n')

 trqualswVals={}
 trqualswAMPL = ampl.getParameter('trqualsw')
 trcostswVals={}
 trcostswAMPL = ampl.getParameter('trcostsw')
 #trstatVals={}
 trstatAMPL = ampl.getParameter('trstat')
 trratingVals={}
 trratingAMPL = ampl.getParameter('trratingct')
 for i in range(1,1+acn.transcount_original):
  if passnum == 2 or passnum == 4 or passnum == 6:
     trqualswVals[i] = acn.ourtrans[i].swqual
  else:
     trqualswVals[i] = 0
  trcostswVals[i] = acn.ourtrans[i].csw
  #trstatVals[i] = acn.ourtrans[i].stat
  trratingVals[i] = acn.ourtrans[i].rating_c #use rating_c for ctgy
 trqualswAMPL.setValues(trqualswVals)
 trcostswAMPL.setValues(trcostswVals)
 #trstatAMPL.setValues(trstatVals)
 trstatAMPL.setValues(basetrON)
 trratingAMPL.setValues(trratingVals)
 if initwithbase:
     trONAMPL = ampl.getVariable('trON')
     trONAMPL.setValues(basetrON)

 trbusiAMPL = ampl.getParameter('trbusi')
 trbusjAMPL = ampl.getParameter('trbusj')
 tr_o_or_rAMPL = ampl.getParameter('tr_o_or_r')
 trpAMPL = ampl.getParameter('trp')

 trbusiVals = {}
 trbusjVals = {}
 tr_o_or_rVals = {}
 trpVals={}

 for i in range(acn.transcount):
  trbusiVals[i+1] = acn.transATbuses[i].busi0+1 #in acn.py, this is busi for the parent transformer
  trbusjVals[i+1] = acn.transATbuses[i].busj0+1
  tr_o_or_rVals[i+1] = 2*(acn.transATbuses[i].original_or_reversed) - 1 # so, = 1 if original, = -1 if reversed
  trpVals[i+1] = acn.transATbuses[i].parentcount1
  #print(i+1,ntrpVals[i+1])

 trbusiAMPL.setValues(trbusiVals)
 trbusjAMPL.setValues(trbusjVals)
 tr_o_or_rAMPL.setValues(tr_o_or_rVals)
 trpAMPL.setValues(trpVals)

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
 for i in range(1,1+acn.transcount_original):
  gmfVals[i] = acn.ourtrans[i].gmf
  bmfVals[i] = acn.ourtrans[i].bmf
  barxstfVals[i] = acn.ourtrans[i].barxstf
  maxxstfVals[i] = acn.ourtrans[i].maxxstf
  minxstfVals[i] = acn.ourtrans[i].minxstf
  brevetaufVals[i] = acn.ourtrans[i].brevetauf
  brevethetafVals[i] = acn.ourtrans[i].brevethetaf
  taustfVals[i] = acn.ourtrans[i].taustf
  tau0_fixedVals[i] = acn.ourtrans[i].tau0f
  thetastfVals[i] = acn.ourtrans[i].thetastf
  theta0fVals[i] = acn.ourtrans[i].theta0f
  FtauvalidVals[i] = acn.ourtrans[i].Ftauvalid
  FthetavalidVals[i] = acn.ourtrans[i].Fthetavalid
  FnuvalidVals[i] = acn.ourtrans[i].Fnuvalid
  numnuxVals[i] = acn.ourtrans[i].numnux


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
  #log.joint('      trans %d numnux %d impedancecorrected %d\n' %(i, acn.ourtrans[i].numnux, acn.ourtrans[i].impedancecorrected))
  if acn.ourtrans[i].impedancecorrected==1:
      for xposition in range(1,1+acn.ourtrans[i].numnux):
          xstf = xposition - 1 + acn.ourtrans[i].minxstf
          truegVals[i,xstf] = acn.ourtrans[i].trueg[xposition]
          truebVals[i,xstf] = acn.ourtrans[i].trueb[xposition]

      #print(i, xposition, xstf)
     #log.joint('trans %d m %d trueg %g trueb %g\n' %(i, xposition, acn.ourtrans[i].trueg[j], acn.ourtrans[i].trueb[j]))
  else:
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
 loadP0AMPL = ampl.getParameter('loadP0')
 prdAMPL = ampl.getParameter('loadprdctg')
 pruAMPL = ampl.getParameter('loadpructg')
 tminVals={}
 tmaxVals={}
 pexistVals = {}
 qexistVals = {}
 prdVals = {}
 pruVals = {}
 #if initwithbase:
 #    loadPAMPL = ampl.getVariable('loadP')
 #    loadPAMPL.setValues(baseloadP)

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
  prdVals[count]= load.prdctg
  pruVals[count]= load.pructg

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
 loadP0AMPL.setValues(baseloadP)
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
  if initwithbase:
      vmagInit[count+1] = basebusVmags[count+1] #flat start or prior value?
      vangleInit[count+1] = basebusVangles[count+1] #flat start or prior value?
      #vmagInit[count+1] = 1.0 #flat start or prior value?
      #vangleInit[count+1] = 0.0 #flat start or prior value?
  else:
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
  #log.joint(' switched shunt ' + str(h) + '\n')

  oursh = acn.ourswshunts[h]
  for a in range(1,1+oursh.len):
   bhaVals[h,a] = oursh.b[a]
   nhaVals[h,a] = oursh.n[a]
   #log.joint('    h = %d a = %d b = %f n = %d\n' %(h,a,bhaVals[h,a],nhaVals[h,a]))
 bhaAMPL.setValues(bhaVals)
 nhaAMPL.setValues(nhaVals)
 #breakexit('done switched')
 forcefixing = 0
 log.joint(' now fixing some integer/binary variables; forcefixing = ' + str(forcefixing) +'\n')

 log.joint(' fixing nontransformer switching variables\n')
 UPntrVals ={}
 LOntrVals ={}

 fixednontranscount = 0
 for cnt1 in range(1,1+acn.nontranscount_original):
  thisnontrans = acn.ournontrans[cnt1]
  UPntrVals[cnt1] = 1.0
  LOntrVals[cnt1] = 0.0
  doit = inexceptions1[thisnontrans.ourfrombus.count] + inexceptions1[thisnontrans.ourtobus.count]
  fixit = (doit == 0) + forcefixing

  if cnt1 == ntrout_index:
      #breakexit('le match')
      fixit = 0


  #fixit = 0

  if fixit > 0:

   value = basentrON[cnt1] # acn.ournontrans[cnt1].status

   #log.joint(' nontrans %d being fixed to %d\n' %(cnt1, value))
   UPntrVals[cnt1] = value
   LOntrVals[cnt1] = value

   fixednontranscount += 1
  else:
   UPntrVals[cnt1] = 1
   LOntrVals[cnt1] = 0

   log.joint(' nontrans %d loose\n' %(cnt1))
   #breakexit('fixed')

 log.joint(' fixed %d nontrans vars\n' %(fixednontranscount))

 UPntrAMPL = ampl.getParameter('UPntr')
 LOntrAMPL = ampl.getParameter('LOntr')

 UPntrAMPL.setValues(UPntrVals)
 LOntrAMPL.setValues(LOntrVals)

 log.joint(' fixing transformer switching variables\n')
 UPtrVals ={}
 LOtrVals ={}

 fixedtranscount = 0
 for cnt1 in range(1,1+acn.transcount_original):
  thistrans = acn.ourtrans[cnt1]

  UPtrVals[cnt1] = 1.0
  LOtrVals[cnt1] = 0.0

  doit = inexceptions1[thistrans.ourfrombus.count] + inexceptions1[thistrans.ourtobus.count]
  fixit = (doit == 0) + forcefixing

  if cnt1 == trout_index:
      #breakexit('le match')
      fixit = 0

  #fixit = 0

  if fixit > 0:
   value =  basetrON[cnt1] # acn.ournontrans[cnt1].statusthistrans.stat
   #log.joint(' trans %d being fixed to %d\n' %(cnt1, value))
   UPtrVals[cnt1] = value
   LOtrVals[cnt1] = value

   fixedtranscount += 1
  else:

   log.joint(' trans %d loose\n' %(cnt1))
   #breakexit('fixed')
 log.joint(' fixed %d trans vars\n' %(fixedtranscount))
 UPtrAMPL = ampl.getParameter('UPtr')
 LOtrAMPL = ampl.getParameter('LOtr')

 UPtrAMPL.setValues(UPtrVals)
 LOtrAMPL.setValues(LOtrVals)

 log.joint(' fixing transformer tap variables\n')

 fixedbin = ampl.getSet('fixedbinxstf')
 fixedbinVals = {}
 numfixedtap = 0
 numloosetap = 0
 for cnt1 in range(1,1+acn.transcount_original):
  thistrans = acn.ourtrans[cnt1]
  fixedtap = []
  #fix xstf if either end is in the exceptions set
  doit = inexceptions1[thistrans.ourfrombus.count] + inexceptions1[thistrans.ourtobus.count]

  fixit = (doit == 0) + forcefixing

  if cnt1 == trout_index:
      #breakexit('le match')
      fixit = 0


  #fixit = 0

  if fixit == 0:
      heurmaxxstfVals[cnt1] = acn.ourtrans[cnt1].maxxstf
      heurminxstfVals[cnt1] = acn.ourtrans[cnt1].minxstf
      #print(' <<<<< ', cnt1, thistrans.ourfrombus.count, thistrans.ourtobus.count, doit, heurmaxxstfVals[cnt1], heurminxstfVals[cnt1])
  else:
      perturbdelta = 0
      heurmaxxstfVals[cnt1] = basexstf[cnt1] + perturbdelta
      if heurmaxxstfVals[cnt1] > acn.ourtrans[cnt1].maxxstf:
       heurmaxxstfVals[cnt1] = acn.ourtrans[cnt1].maxxstf
      heurminxstfVals[cnt1] = basexstf[cnt1] - perturbdelta
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
       numloosetap += intmax - intmin
       #if thisfixedtap > 0:
       # log.joint(' fixedtap ' + str(thisfixedtap) + ' bin variables\n')
  fixedbinVals[cnt1] = fixedtap
  numfixedtap += len(fixedtap)
  #transatbus[cnt1].setValues(transatbusVals[cnt1])
  fixedbin[cnt1].setValues(fixedbinVals[cnt1])
 log.joint(' number of binxstf variables fixed at 0: ' + str(numfixedtap) + '; loose: ' + str(numloosetap) + '\n')

 heurmaxxstfAMPL.setValues(heurmaxxstfVals)
 heurminxstfAMPL.setValues(heurminxstfVals)

 log.joint('fixing generator binary variables\n')
 UPgenVals = {}
 LOgenVals = {}
 UPgenstartupVals = {}
 LOgenstartupVals = {}
 UPgenshutdownVals = {}
 LOgenshutdownVals = {}
 fgen = 0
 f0 = f1 = 0

 for gencount in range(1, acn.numgens+1):

  thisgen = acn.ourgens[gencount]
  dofixit = inexceptions1[thisgen.ourbus.count] == 0 # and (gencount != genoutIndex or thisgen.id != genoutID)
  #dofixit = 1 when bus is not in the list of exceptions
  fixit = dofixit + (forcefixing and (gencount != genoutIndex or thisgen.id != genoutID))
  #second term is positive if we are forcing fixing and the generator is not the generator out

  UPgenVals[gencount] = 1
  LOgenVals[gencount] = 0
  UPgenstartupVals[gencount] = 1
  LOgenstartupVals[gencount] = 0

  UPgenshutdownVals[gencount] = 1
  LOgenshutdownVals[gencount] = 0

  #print(gencount, basegenON[gencount])
  if fixit > 0:
   fgen += 1
   thevalue = basegenON[gencount] #basegenOn set in acnsolverAMPL indexed from 0 (written to file around line 51)
                                  #but it is read in acnsolverAMPL2 indexed from 1 (around line 345)
   #if thevalue == -1:
   # breakexit('wowie ' + str(gencount))
   UPgenVals[gencount] = thevalue
   LOgenVals[gencount] = thevalue

   UPgenstartupVals[gencount] = 0
   LOgenstartupVals[gencount] = 0

   UPgenshutdownVals[gencount] = 0
   LOgenshutdownVals[gencount] = 0



   #log.joint('gencount %d fixed at %d\n' %(gencount,thevalue))


   #log.joint('fixed gen count %d to %g\n' %(gencount,value)
   f0 += (thevalue == 0)
   f1 += (thevalue == 1)

 UPgenAMPL = ampl.getParameter('UPgen')
 LOgenAMPL = ampl.getParameter('LOgen')

 UPgenAMPL.setValues(UPgenVals)
 LOgenAMPL.setValues(LOgenVals)

 UPgenstartupAMPL = ampl.getParameter('UPgenstartup')
 LOgenstartupAMPL = ampl.getParameter('LOgenstartup')

 UPgenstartupAMPL.setValues(UPgenstartupVals)
 LOgenstartupAMPL.setValues(LOgenstartupVals)

 UPgenshutdownAMPL = ampl.getParameter('UPgenshutdown')
 LOgenshutdownAMPL = ampl.getParameter('LOgenshutdown')


 UPgenshutdownAMPL.setValues(UPgenshutdownVals)
 LOgenshutdownAMPL.setValues(LOgenshutdownVals)

 log.joint('fixed %d generators (0: %d, 1: %d)\n' %(fgen, f0, f1))

 log.joint(' fixing switchable shunt variables\n')

 numxha = 0
 numfixedxha = 0
 UPxhaVals = {}
 LOxhaVals = {}

 for h in range(1,1 + acn.numswitchedshunts):
  oursh = acn.ourswshunts[h]
  #dofixit = (inexceptions1[oursh.ourbus.count] is False)
  dofixit = inexceptions1[oursh.ourbus.count] == 0
  #print('h',h,'buscount',oursh.ourbus.count,'shcount',oursh.count,'dofixit',dofixit)

  if oursh.status:
    #print('>>>', oursh.n)

    for a in range(1,1+oursh.len):
        #print('-->h',h,'a',a,'n',oursh.n[a])
        if oursh.n[a] > 0:
            numxha += 1

        if dofixit and (h,a) in basesw:
            UPxhaVals[h, a] = basesw[h,a]
            LOxhaVals[h, a] = basesw[h,a]
            numfixedxha += 1
        else:
            if oursh.n[a] > 0:
                UPxhaVals[h, a] = oursh.n[a]
                LOxhaVals[h, a] = 0
        #print('there')


 log.joint('fixed %d switchable shunt variables\n' %(numfixedxha))
 #breakexit('foooo')


 UPxhaAMPL = ampl.getParameter('UPxha')
 LOxhaAMPL = ampl.getParameter('LOxha')

 #print('UPxhaVals', UPxhaVals)
 UPxhaAMPL.setValues(UPxhaVals)
 LOxhaAMPL.setValues(LOxhaVals)

 if initwithbase and numxha>0:
     xhaAMPL = ampl.getVariable('xha')
     xhaAMPL.setValues(basesw)

 write_ampl_mod = 0
 if write_ampl_mod:
     amplstate = 'expand; display {j in 1.._nvars} (_varname[j],_var[j].lb0,_var[j].\
lb,_var[j].ub,_var[j].ub0);'     #shows full model
     filename = 'ctgmodel_' + str(conlabel) + '.out'
     modelout = ampl.getOutput(amplstate)
     outfile = open(filename,"w")
     outfile.write("model = " + str(modelout) + "\n")
     outfile.close()

 # Set Knitro options
 acn.maxtime = 900
 maxtime = ' maxtime_real=' + str(acn.maxtime) + ' mip_maxtime_real=' + str(acn.maxtime) + ' ms_maxtime_real=' + str(acn.maxtime) # just in case we use multi-start
 #maxtime = ' maxtime_real=600 '
 ctg_options = 'outlev=4 outmode=0 debug=0 linsolver=6 outname="knitro.log" feastol=1e-5 ftol=1e-3 scale=0 honorbnds=1 cg_maxit=1 bar_murule=0 bar_refinement=1 bar_switchrule=0 bar_feasible=1 restarts=3 maxit=3000 presolve_tol=0.5'  # use bar_initpt=2 bar_directinterval=0 ?
 if acn.useparallel:
     ctg_options += ' par_numthreads=1'
 just_feasible = 1 #ignoring optimality works quite well
 if just_feasible:
     ctg_options += ' opttol=1.0e30 opttol_abs=1e30'
 else:
     ctg_options += ' opttol=1.0e-3'
 if (acn.numbuses < 10000):
     ctg_options += ' bar_maxcrossit=1'
 mip_options = ' mip_terminate=1 mip_outinterval=1 mip_outlevel=3 mip_debug=0 mip_outsub=0 mip_nodealg=1'
 if (1):
     intvar_strategy = ' mip_intvar_strategy=0'  # handle integer vars directly
 else:
     intvar_strategy = ' mip_intvar_strategy=2'  # convert binary vars to complementarity constraints
 if (0):
   relax = ' relax=0'
   log.joint(' using relax = 0 ...\n')
 else:
   relax = ' relax=1'
   log.joint(' using relax = 1 ...\n')
 options = ctg_options + mip_options + maxtime + relax + intvar_strategy
 #print "options=",options
 ampl.setOption('knitro_options', options)

 # Set relaxation variables in balance constraints to make them feasible?
 initfeasslacks = 0
 if initfeasslacks == 1:
    print("Initializing relaxation variables...")
    maxviol = 0.0
    # Get current value of balance constraints and use these to
    # initialize relaxation variables
    pbalanceAMPL = ampl.getConstraint('PBalance')
    relaxPplusA = ampl.getVariable('relaxPplus')
    relaxPminusA = ampl.getVariable('relaxPminus')
    suffixes = ['body']
    df = pbalanceAMPL.getValues(suffixes)
    bodies = df.toDict()
    suffixes = ['lb']
    df = pbalanceAMPL.getValues(suffixes)
    lbs = df.toDict()
    plus = {}
    minus = {}
    for k in range(1, numbuses+1):
       #body = pbalanceAMPL[k].body()
       #lb = pbalanceAMPL[k].lb()
       val = bodies[k] - lbs[k]  #compute initial error/violation in the P balance constraint
       #print("val=", val)
       if (val > 0):
         plus[k] = val
         minus[k] = 0
         #relaxPplusA[k].setValue(val)
         maxviol = max(maxviol, val)
       else:
         plus[k] = 0
         minus[k] = -val
         maxviol = max(maxviol, -val)
         #relaxPminusA[k].setValue(-val)
    relaxPplusA.setValues(plus)
    relaxPminusA.setValues(minus)
    print("PBalance: maxviol=", maxviol)

    maxviol=0.0
    qbalanceAMPL = ampl.getConstraint('QBalance')
    relaxQplusA = ampl.getVariable('relaxQplus')
    relaxQminusA = ampl.getVariable('relaxQminus')
    suffixes = ['body']
    df = qbalanceAMPL.getValues(suffixes)
    bodies = df.toDict()
    suffixes = ['lb']
    df = qbalanceAMPL.getValues(suffixes)
    lbs = df.toDict()
    for k in range(1, numbuses+1):
       #body = qbalanceAMPL[k].body()
       #lb = qbalanceAMPL[k].lb()
       val = bodies[k] - lbs[k]  #compute initial error/violation in the Q balance constraint
       #print("val=", val)
       if (val > 0):
         plus[k] = val
         minus[k] = 0
         maxviol = max(maxviol, val)
         #relaxQplusA[k].setValue(val)
       else:
         plus[k] = 0
         minus[k] = -val
         maxviol = max(maxviol, -val)
         #relaxQminusA[k].setValue(-val)
    relaxQplusA.setValues(plus)
    relaxQminusA.setValues(minus)
    print("QBalance: maxviol=", maxviol)

    maxviol=0.0
    # Get current value of line current constraints and use these to
    # initialize relaxation variables
    #linecurrentT = ampl.getConstraint('LineCurrentTR')
    relaxlineT = ampl.getVariable('relaxStr')
    trPAMPL = ampl.getVariable('trP')
    dfP = trPAMPL.getValues()
    valsP = dfP.toDict()
    trQAMPL = ampl.getVariable('trQ')
    dfQ = trQAMPL.getValues()
    valsQ = dfQ.toDict()
    trpAMPL = ampl.getParameter('trp')
    trrateAMPL = ampl.getParameter('trratingct')
    viol = {}
    for j in range(1,1+acn.transcount_original):
        viol[j] = 0.0
    for i in range(1,1+acn.transcount):
        j = int(trpAMPL[i])
        #print("trP=",valsP[i]," trQ=",valsQ[i], " trrating=", trrateAMPL[j])
        val = numpy.sqrt(valsP[i]*valsP[i] + valsQ[i]*valsQ[i]) - trrateAMPL[j]
        if (val > 0):
            viol[j] = max(viol[j], val)
         #  relaxlineT[i].setValue(val)
            maxviol = max(maxviol,val)
    relaxlineT.setValues(viol)
    print("LineCurrentTR: maxviol=", maxviol)

    #breakexit("initfeas")

 # Solve
 usescript = 1
 time1=time2 = 0
 try:
  time1 = time.time()
  if usescript:
      ampl.eval('include solve_script2;')
      #solveout = ampl.getOutput('include solve_script2;')
  else:
      #ampl.solve()
      #ampl.eval('solve > junk.txt;')  # use this to re-direct ampl output to junk file
      #ampl.eval('solve > /dev/null;') # use this to re-direct ampl output to null device
      solveout = ampl.getOutput('solve;')
  time2 = time.time()
 except:
  log.joint("Error calling ampl.solve!")
  var = traceback.format_exc()
  #traceback.print_exc()
  #acn_evaluation2.print_alert(var, raise_exception=False)
  pass

 log.joint(' elapsed time in solve_script2: ' +str(time2 - time1) + '\n')
 # Re-solve with integer vars fixed

 log.joint(' now calling resolve_fixed2H2\n')
 time1 = time.time()
 ampl.setOption('presolve_eps', 1.0e-6)
 resolve_fixed2H2(acn, ampl, usescript, basegenON)
 time2 = time.time()
 log.joint(' came back from resolve_fixed2H2, time: %g\n' %(time2 - time1))

 # Write some AMPL solution files to file
 write_ampl_sol = 0
 if write_ampl_sol:
     #amplstate = 'display gencost0, genon, genstartup, genshutdown, genP, genQ, minpg, maxpg, minqg, maxqg;'
     amplstate = 'display sum{g in 1..G: g != genout_index} Deltact*(-gencost0[g] - Oncost[g]*genon[g]) - sum{g in 1..G: g != genout_index} (Sucost[g]*genstartup[g] + Sdcost[g]*genshutdown[g]);'
     amplstate += 'display -sum{i in 1..I} Deltact*buscost0[i];'
     #amplstate += 'display -sum{i in 1..I} Deltact*(sum{n in 1..Pnumcblocks} Pcostcblock[n]*(Pplusblock0[i,n] + Pminusblock0[i,n]));'
     #amplstate += 'display -sum{i in 1..I} Deltact*(sum{n in 1..Qnumcblocks} Qcostcblock[n]*(Qplusblock0[i,n] + Qminusblock0[i,n]));'
     amplstate += 'display sum{l in 1..L: loadisactive[l] > 0} Deltact*loadcost0[l];'
     amplstate += 'display -sum{e in 1..NTRcnt_original: ntrqualsw[e]>0 and ntrstat[e]=0 and e != ntrout_index} ntrcostsw[e]*ntrON[e];'
     amplstate += 'display -sum{e in 1..NTRcnt_original: ntrqualsw[e]>0 and ntrstat[e]>0 and e != ntrout_index} ntrcostsw[e]*(ntrstat[e]-ntrON[e]);'
     amplstate += 'display -sum{f in 1..TRcnt_original: trqualsw[f]>0 and trstat[f]=0 and f != trout_index} trcostsw[f]*trON[f];'
     amplstate += 'display -sum{f in 1..TRcnt_original: trqualsw[f]>0 and trstat[f]>0 and f != trout_index} trcostsw[f]*(trstat[f]-trON[f]);'
     amplstate += 'display gencommit;'
     amplstate += 'display PBalance;'
     amplstate += 'display QBalance;'
     amplstate += 'display LineCurrentNTR;'
     amplstate += 'display LineCurrentTR;'
     amplstate += 'display ramping;'
     filename = 'amplsol_'+ str(conlabel) +'.out'
     modelout = ampl.getOutput(amplstate)
     outfile = open(filename,"w")
     outfile.write("objective and constraint values: \n" + str(modelout) + "\n")
     outfile.close()


 # Write the contingency solution to file
 ctgytmpfile = "solution2_" + str(conlabel) + ".txt"
 write_ctg_solutionH2(log, acn, ampl, genoutIndex, genoutID, lineout_i, lineout_j, lineout_ckt, numxha, ctgytmpfile)
 #breakexit('done write_ctg')

 log.joint('wrote contingency solution to ' +  ctgytmpfile + '\n')

 # Evaluate contingency solution
 score = -1e20
 try:
  ctgsol_obj, ctgsol_infeas, exist, ctgsol_summary = acn_evaluation2.ctgy_from_file(acn, conlabel, ctgytmpfile)
  #breakexit("evalctg")
  if (ctgsol_infeas):
      log.joint(str(conlabel) + ':optimized solution obtained through H2 infeasible!\n')
  else:
      log.joint(str(conlabel) + ':optimized solution obtained through H2 feasible!\n')
  ctgyfile = "solution_" + str(conlabel) + ".txt"

  log.joint('\n===\nbegin evaluation of prior infeasible solution in file ' + ctgyfile +'\n')
  #print('ctgyfile=',ctgyfile)
  pobj, pinfeas, pexist, ctgsol_psummary = acn_evaluation2.ctgy_from_file(acn, conlabel, ctgyfile)
  if (pinfeas) or (pexist == False):
      log.joint(str(conlabel) + ':prior solution infeasible!\n')
      score = -1e20
  else:
      score = pobj
  log.joint('completed eval of prior infeasible solution\n===\n')
  log.joint('comparison: prior infeasible solution: obj=' + str(pobj)+'\n')
  log.joint('comparison: prior infeasible solution: infeas=' + str(pinfeas)+'\n')
  log.joint('contingency solution exists=' +  str(exist)+'\n')
  if (exist == True):
    log.joint('ctgsol_obj=' + str(ctgsol_obj)+'\n')
    log.joint('ctgsol_infeas=' + str(ctgsol_infeas)+'\n')
    if ctgsol_infeas == 0 and (pexist == False or ctgsol_obj > pobj):
        # We've computed a contingency case solution that is better than the
        # 'infeasible_solution', so overwrite with this one
        #os.popen('\cp solution2.txt solution_???.txt')
        log.joint(' solution computed through H2 is better than infeasible solution, so overwriting\n')
        shutil.copyfile(ctgytmpfile, ctgyfile)
        score = ctgsol_obj
        #breakexit("overwrite")
    else:
        log.joint('using infeasible solution!\n')
  #if (ctgsol_infeas or pinfeas):
  #    breakexit("infeasible")
  #if (ctgsol_obj < 0.0):
  #    log.joint(str(conlabel) + ': bad solution!\n')
  #    breakexit("2: bad ctg solution")
 except:
  #breakexit("acnsolverAMPLcon: Error evaluating the contingency solution!")
  print("acnsolveconthruH2: Error evaluating the contingency solution!")
  pass
  '''
  var = traceback.format_exc()
  traceback.print_exc()
  acn_evaluation2.print_alert(var, raise_exception=False)
  pass
  '''

 #Write the score for this contingency to a file
 scorefile = "score1_" + str(conlabel) + ".txt"
 score1 = open(scorefile,"w")
 score1.write(str(score)+"\n")
 score1.close()

##### Called after initial solve to fix all integer variables and re-solve
#Maybe put this in a different file?


def resolve_fixed2H2(acn, ampl, usescript, basegenON):

  log = acn.log
  log.joint('running resolve_fixed2H2\n')

  # if usescript = 1, then we already fixed everything inside the AMPL script after
  # the initial solve, so we can skip the fixing code below

  if usescript == 0:

      # Round and fix all integer variables in the model

      # Transformer tap settings
      v = ampl.getVariable('xstf')
      df = v.getValues()
      vals = df.toDict()
      fixedVals = {}
      for i in range(1,1+acn.transcount_original):
          tmpval = min(vals[i],acn.ourtrans[i].maxxstf)
          tmpval = max(tmpval, acn.ourtrans[i].minxstf)
          tmpval = int(round(tmpval))
          fixedVals[i] = tmpval
      v.setValues(fixedVals)
      v.fix()

      # Switched shunt susceptance
      v = ampl.getVariable('xha')
      df = v.getValues()
      vals = df.toDict()
      fixedVals = {}
      for h in range(1,1 + acn.numswitchedshunts):
          if acn.ourswshunts[h].status > 0:  #switch shunt active
              oursh = acn.ourswshunts[h]
              for a in range(1,1+oursh.len):
                  upbnd = int(oursh.n[a])
                  if upbnd > 0:
                      xhaval = int(vals[h,a])
                      tmpval = min(xhaval,upbnd)
                      tmpval = max(tmpval, 0)
                      fixedVals[h,a] = tmpval
                      #print('h=',h,' a=',a,' xhaval=',xhaval,' upbnd=',upbnd,' fixed=',tmpval)
      v.setValues(fixedVals)
      v.fix()

      # Line switching
      v = ampl.getVariable('ntrON')
      df = v.getValues()
      vals = df.toDict()
      fixedVals = {}
      for i in range(1,1+acn.nontranscount_original):
          if vals[i] < 0.5:
              fixedVals[i] = 0
          else:
              fixedVals[i] = 1
      v.setValues(fixedVals)
      v.fix()

      # Transformer switching
      v = ampl.getVariable('trON')
      df = v.getValues()
      vals = df.toDict()
      fixedVals = {}
      for i in range(1,1+acn.transcount_original):
          if vals[i] < 0.5:
              fixedVals[i] = 0
          else:
              fixedVals[i] = 1
      v.setValues(fixedVals)
      v.fix()

      # Generator commitment decision variables
      genon = ampl.getVariable('genon')
      genondf = genon.getValues()
      genonvals = genondf.toDict()
      gensu = ampl.getVariable('genstartup')
      gensudf = gensu.getValues()
      gensuvals = gensudf.toDict()
      gensd = ampl.getVariable('genshutdown')
      gensddf = gensd.getValues()
      gensdvals = gensddf.toDict()
      fixedonVals = {}
      fixedsuVals = {}
      fixedsdVals = {}
      for i in range(1, acn.numgens+1):
          ourgen = acn.ourgens[i]
          # First fix generator on or off and then fix start up/shutdown values
          if genonvals[i] < 0.5:
              fixedonVals[i] = 0
              if basegenON[i] > 0:  # generator currently on; shut it off
                  fixedsuVals[i] = 0
                  fixedsdVals[i] = 1
              else:
                  fixedsuVals[i] = 0
                  fixedsdVals[i] = 0
          else:
              fixedonVals[i] = 1
              if basegenON[i] > 0:
                  fixedsuVals[i] = 0
                  fixedsdVals[i] = 0
              else:                # generator currently off; turn it on
                  fixedsuVals[i] = 1
                  fixedsdVals[i] = 0
      genon.setValues(fixedonVals)
      genon.fix()
      gensu.setValues(fixedsuVals)
      gensu.fix()
      gensd.setValues(fixedsdVals)
      gensd.fix()
      #print("genonFIXED=",fixedonVals)
      #print("genstartupFIXED=",fixedsuVals)
      #print("genshutdownFIXED=",fixedsdVals)

  # Now set options and solve with integer variables fixed

  try:
    acn.maxtime = 900
    maxtime = ' maxtime_real=' + str(acn.maxtime) + ' mip_maxtime_real=' + str(acn.maxtime) + ' ms_maxtime_real=' + str(acn.maxtime) # just in case we use multi-start
    resolve_options = 'outlev=4 outmode=0 outname="knitro-fixed.log" debug=0 feastol=1e-5 feastol_abs=9e-5 ftol=1e-6 scale=0 honorbnds=0 cg_maxit=50 bar_murule=0 bar_feasible=0 bar_refinement=1 bar_initpi_mpec=0.0 maxit=3000 alg=1 strat_warm_start=1 bar_initpt=2 bar_initmu=1e-6 bar_slackboundpush=1e-6 infeastol=1e-5 restarts=1 presolve_initpt=1 presolve_tol=0.5'
    #resolve_options = 'outlev=4 outmode=2 debug=1 feastol=1e-5 opttol=1e-3 ftol=1e-6 scale=0 honorbnds=0 cg_maxit=50 bar_murule=1 strat_warm_start=1 bar_refinement=1 bar_initpi_mpec=0.0 maxit=100 alg=3'
    mip_options = ' mip_terminate=1 mip_outinterval=1 mip_outlevel=1 mip_debug=0 mip_outsub=2 mip_nodealg=1 mip_intvar_strategy=0 mip_numthreads=4'
    resolve_options += mip_options
    if acn.useparallel:
      resolve_options += ' par_numthreads=1'
    just_feasible = 0
    if just_feasible:
      resolve_options += ' opttol=1.0e30 opttol_abs=1e30'
    else:
      resolve_options += ' opttol=1.0e-3'
    if (acn.numbuses < 10000):
      resolve_options += ' bar_maxcrossit=1'
    ampl.setOption('knitro_options', resolve_options)
    ampl.eval('include resolve_script2;')
    #solveout = ampl.getOutput('include resolve_script2;')
    #ampl.solve()
    #ampl.eval('solve > junk.txt')  # use this to re-direct ampl output to junk file
    #ampl.eval('solve > /dev/null') # use this to re-direct ampl output to null device
    #breakexit("resolve")
  except:
    log.joint("Error in fixed re-solve!\n")
    pass
  #breakexit("bbb")


def write_ctg_solutionH2(log, acn, ampl, genoutIndex, genoutID, lineout_i, lineout_j, lineout_ckt, numxha, ctgytmpfile):

    # Extract and write contingency solution (move to function? create solution class?)
    sbase = acn.options['sbase'] #acn.raw.case_identification.sbase #needed?
    busVmags = numpy.zeros(acn.numbuses)
    v = ampl.getVariable('busVmags')
    df = v.getValues()
    vals = df.toDict()
    for i in range(acn.numbuses):
        busVmags[i] = vals[i+1]

    busVangles = numpy.zeros(acn.numbuses)
    v = ampl.getVariable('busVangles')
    df = v.getValues()
    vals = df.toDict()
    for i in range(acn.numbuses):
        busVangles[i] = vals[i+1]

    clearedloads = numpy.zeros(1+acn.numloads)
    v = ampl.getVariable('t0')
    dv = v.getValues()
    vals = dv.toDict()
    for i in range(acn.numloads):
        clearedloads[i+1] = vals[i+1]

    #ploads = numpy.zeros(acn.numloads+1)
    #qloads = numpy.zeros(acn.numloads+1)
    #u = ampl.getVariable('loadP')
    #du = u.getValues()
    #uals = du.toDict()
    #w = ampl.getVariable('loadQ')
    #dw = w.getValues()
    #wals = dw.toDict()
    #for i in range(acn.numloads):
    #    ploads[i+1] = uals[i+1]
    #    qloads[i+1] = wals[i+1]

    genP = numpy.zeros(acn.numgens)
    v = ampl.getVariable('genP')
    df = v.getValues()
    vals = df.toDict()
    for i in range(acn.numgens):
        genP[i] = vals[i+1]

    gencost = numpy.zeros(acn.numgens)
    v = ampl.getVariable('gencost0')
    df = v.getValues()
    vals = df.toDict()
    #log.joint(' in acnsolverAMPL gencosts\n')
    sumgen = 0
    for i in range(acn.numgens):
        gencost[i] = vals[i+1]
        sumgen += vals[i+1]
        #log.joint(' gen %d cost %g sumgen %.8e\n' %(i+1, vals[i+1], sumgen))

    genQ = numpy.zeros(acn.numgens)
    v = ampl.getVariable('genQ')
    df = v.getValues()
    vals = df.toDict()
    for i in range(acn.numgens):
        genQ[i] = vals[i+1]

    genON = numpy.zeros(acn.numgens,dtype=int)
    gensu = numpy.zeros(acn.numgens,dtype=int)
    gensd = numpy.zeros(acn.numgens,dtype=int)
    v = ampl.getVariable('genon')
    df = v.getValues()
    valson = df.toDict()

    #log.joint(' in acnsolverAMPL genONcosts\n')
    #sumgenon = 0
    for i in range(acn.numgens):
        genON[i] = valson[i+1]
        #thiscost = valson[i+1]*oncost[i+1]
        #sumgenon += thiscost
        #log.joint(' gen %d on %g cost %g sumgenon %.8e\n' %(i+1, valson[i+1], thiscost, sumgenon))


    #log.joint(' in acnsolverAMPL gensucosts\n')
    vsu = ampl.getVariable('genstartup')
    df = vsu.getValues()
    valssu = df.toDict()
    #sumgensu = 0
    for i in range(acn.numgens):
        gensu[i] = valssu[i+1]
        #thiscost = valssu[i+1]*sucost[i+1]
        #sumgensu += thiscost
        #log.joint(' gen %d on %g cost %g sumgensu %.8e\n' %(i+1, valssu[i+1], thiscost, sumgensu))

    vsd = ampl.getVariable('genshutdown')
    df = vsd.getValues()
    valssd = df.toDict()
    for i in range(acn.numgens):
        gensd[i] = valssd[i+1]

    #sumgensd = 0
    for i in range(acn.numgens):
        gensd[i] = valssd[i+1]
        #thiscost = valssd[i+1]*sdcost[i+1]
        #sumgensd += thiscost
        #log.joint(' gen %d on %g cost %g sumgensd %.8e\n' %(i+1, valssd[i+1], thiscost, sumgensd))

    relaxPplus = numpy.zeros(acn.numbuses)
    v = ampl.getVariable('relaxPplus')
    df = v.getValues()
    vals = df.toDict()
    for i in range(acn.numbuses):
        relaxPplus[i] = vals[i+1]
    #print(relaxPplus)

    relaxPminus = numpy.zeros(acn.numbuses)
    v = ampl.getVariable('relaxPminus')
    df = v.getValues()
    vals = df.toDict()
    for i in range(acn.numbuses):
        relaxPminus[i] = vals[i+1]
    #print(relaxPminus)

    ntrP = numpy.zeros(acn.nontranscount)
    v = ampl.getVariable('ntrP')
    df = v.getValues()
    vals = df.toDict()
    for i in range(acn.nontranscount):
        ntrP[i] = vals[i+1]

    ntrQ = numpy.zeros(acn.nontranscount)
    v = ampl.getVariable('ntrQ')
    df = v.getValues()
    vals = df.toDict()
    for i in range(acn.nontranscount):
        ntrQ[i] = vals[i+1]

    #breakexit('got ntr flows')
    ntrON = numpy.zeros(acn.nontranscount_original,dtype=int)
    v = ampl.getVariable('ntrON')
    df = v.getValues()
    vals = df.toDict()
    for i in range(acn.nontranscount_original):
        ntrON[i] = int(vals[i+1])

    trON = numpy.zeros(acn.transcount_original,dtype=int)
    taps = numpy.zeros(acn.transcount_original,dtype=int)
    v = ampl.getVariable('trON')
    df = v.getValues()
    vals = df.toDict()
    for i in range(acn.transcount_original):
        trON[i] = int(vals[i+1])
    v = ampl.getVariable('xstf')
    df = v.getValues()
    vals = df.toDict()
    for i in range(acn.transcount_original):
        taps[i] = int(vals[i+1])

    trP = numpy.zeros(acn.transcount)
    v = ampl.getVariable('trP')
    df = v.getValues()
    vals = df.toDict()
    for i in range(acn.transcount):
        trP[i] = vals[i+1]

    trQ = numpy.zeros(acn.transcount)
    v = ampl.getVariable('trQ')
    df = v.getValues()
    vals = df.toDict()
    for i in range(acn.transcount):
        trQ[i] = vals[i+1]

    #sw = numpy.zeros(acn.numswitchedshunts)
    sw={}
    if numxha > 0:
       v = ampl.getVariable('xha')
       df = v.getValues()
       vals = df.toDict()
       for (h,a), step in vals.items():
           #print(h,a,step)
           sw[h,a] = int(step)

    log.joint("Write contingency solution\n")
    sol2 = open(ctgytmpfile,"w")

    sol2.write("--bus section\n")
    sol2.write("i,v,theta\n")
    for j in range(acn.numbuses):
        #sol2.write(str(acn.ourbuses[j+1].theirbus.i) + "," + str(busVmags[j]) + "," + str(busVangles[j]*180/pi) + "\n")
        sol2.write(str(acn.ourbuses[j+1].theirbus.i) + "," + str(busVmags[j]) + "," + str(busVangles[j]) + "\n")

    sol2.write("--load section\n")
    sol2.write("i,id,t\n")
    countactive = 0
    for count in range(1,1+acn.numloads):
        load = acn.ourloads[count]
        if load.status:
            #sol2.write(str(acn.ourloads[count].theirload.i) + "," + str(acn.ourloads[count].theirload.id) + "," + str(clearedloads[countactive]) + "\n")
            #sol2.write(str(load.theirload.i) + "," + str(load.theirload.id) + "," + str(clearedloads[countactive]) + "\n")
            sol2.write(str(load.theirload.i) + "," + str(load.theirload.id) + "," + str(clearedloads[load.count]) + "\n")
            countactive += 1

    sol2.write("--generator section\n")
    sol2.write("i,id,p,q,x\n")
    count = 0
    for gencount in range(1, acn.numgens+1):
        gen = acn.ourgens[gencount]
        pgen = qgen = 0.0
        genon = genON[count]
        # Do not write solution info for removed generators!
        if gen.i != genoutIndex or gen.id != genoutID:
            if genon:
                pgen = genP[count]
                qgen = genQ[count]
            #sol2.write(str(gen.i) + "," + str(gen.id) + "," + str(s*pgen) + "," + str(s*qgen) +  "," + str(genon) + "\n")
            sol2.write(str(gen.i) + "," + str(gen.id) + "," + str(pgen) + "," + str(qgen) +  "," + str(genon) + "\n")
        count += 1

    sol2.write("--line section\n")
    sol2.write("iorig,idest,id,x\n")
    for j in range(acn.nontranscount_original):
        count = j+1
        ntr = acn.ournontrans[count]
        # Do not write solution info for removed line!
        if ntr.i != lineout_i or ntr.j != lineout_j or ntr.ckt != lineout_ckt:
            sol2.write(str(ntr.i) + "," + str(ntr.j) + "," + str(ntr.ckt) + "," + str(ntrON[j]) +  "\n")

    sol2.write("--transformer section\n")
    sol2.write("iorig,idest,id,x,xst\n")
    for j in range(acn.transcount_original):
        count = j+1
        tr = acn.ourtrans[count]
        # Do not write solution info for removed transformer!
        if tr.i != lineout_i or tr.j != lineout_j or tr.ckt != lineout_ckt:
            sol2.write(str(tr.i) + "," + str(tr.j) + "," + str(tr.ckt) + "," + str(trON[j]) + "," + str(taps[j]) + "\n")

    sol2.write("--switched shunt section\n")
    sol2.write("i,xst1,xst2,xst3,xst4,xst5,xst6,xst7,xst8\n")
    for h in range(1,1 + acn.numswitchedshunts):
        oursh = acn.ourswshunts[h]
        if oursh.status:
            sol2.write(str(oursh.i))
            for a in range(1,1+oursh.len):
                if (oursh.n[a]):
                    sol2.write("," + str(sw[h,a]))
            sol2.write("\n")

    sol2.close()


#############
