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
#import acn_evaluation2
from acn_solution import *
from acnsolveconthruH2 import *

def trace(frame, event, arg):
    print("%s, %s:%d" % (event, frame.f_code.co_filename, frame.f_lineno))
    return trace


def record_xstf(acn):
 log = acn.log
 log.joint(' recording xstf\n')
 xstfdict = {}
 xpositiondict = {}
 bfdict = {}
 gfdict = {}



 for count in range(1,1+acn.transcount_original):

     tr = acn.ourtrans[count]

     trueb = tr.trueb
     trueg = tr.trueg
     xstf = 0

     if tr.Ftauvalid:
         xstf = int( round((tr.tau0f - tr.brevetauf)/tr.taustf))

     if tr.Fthetavalid:
         xstf = int( round((tr.theta0f - tr.brevethetaf)/tr.thetastf))
         #print(tr.theta0f, tr.brevethetaf, tr.thetastf, tr.thetastf, tr.minxstf, tr.maxxstf)
         #breakexit('')
     if xstf > tr.maxxstf:
         xstf = tr.maxxstf
     if xstf < tr.minxstf:
         xstf = tr.minxstf
     #if count == 224:
     #    log.joint(' xagain ' + str(xstf) + ' max ' + str(tr.maxxstf) + '\n')
     #    breakexit('heeee')

     xposition = xstf + 1 - tr.minxstf

     xstfdict[count] = xstf

     xpositiondict[count] = xposition

     #print('count',count,'xstf',xstf,'numnux',tr.numnux, 'minxstf',tr.minxstf)
     #print(trueb)
     if tr.impedancecorrected:
         bfdict[count] = trueb[xposition]
         gfdict[count] = trueg[xposition]
     else:
         bfdict[count] = trueb[1]
         gfdict[count] = trueg[1]



 #breakexit('3')
 f = open('basexstf.txt', 'w')
 for count in range(1,1+acn.transcount_original):
     #print(count, xstfdict[count])
     f.write('%d\n' %(xstfdict[count]))
 f.close()

 #breakexit('2')

def write_base_files(acn, basebusVmags,basebusVangles,ploads,basegenP,basegenON,basegensu,basegensd,basentrON,basetrON, basesw):

    f = open('basebusVmags.txt', 'w')
    for i in range(acn.numbuses):
        f.write("%s\n" % basebusVmags[i])
    f.close()

    f = open('basebusVangles.txt', 'w')
    for i in range(acn.numbuses):
        f.write("%s\n" % basebusVangles[i])
    f.close()

    f = open('baseloadP.txt', 'w')
    for i in range(acn.numloads):
        f.write("%s\n" % ploads[i+1])
    f.close()

    f = open('basegenP.txt', 'w')
    for i in range(acn.numgens):
        f.write("%s\n" % basegenP[i])
    f.close()

    f = open('basegenON.txt', 'w')
    for i in range(acn.numgens):
        f.write("%s\n" % basegenON[i])
    f.close()

    f = open('basegenSU.txt', 'w')
    for i in range(acn.numgens):
        f.write("%s\n" % basegensu[i])
    f.close()

    f = open('basegenSD.txt', 'w')
    for i in range(acn.numgens):
        f.write("%s\n" % basegensd[i])
    f.close()

    f = open('basentrON.txt', 'w')
    for i in range(acn.nontranscount_original):
        f.write("%s\n" % basentrON[i])
    f.close()

    f = open('basetrON.txt', 'w')
    for i in range(acn.transcount_original):
        f.write("%s\n" % basetrON[i])
    f.close()

    f = open('basesw.txt', 'w')
    for (h,a) in basesw:
         f.write("%d %d %d\n" %(h,a,basesw[h,a]))
    f.close()

    record_xstf(acn)

def acnsolverAMPL(alldata):

 log  = alldata['log']
 log.joint("\n***\n solving AMPL\n")

 #sys.settrace(trace)

 acn = alldata['acn']
 raw = alldata['RAW']
 con = alldata['CON']
 sup = alldata['JSON']

 pi = math.pi

 numpy.set_printoptions(precision=16)
 elapsedtime = time.time()-acn.timebeg
 remaintime = acn.maxtime - elapsedtime


 mykeybus = alldata['mykeybus']
 mykeynontrdual = alldata['mykeynontrdual']
 mykeytrdual = alldata['mykeytrdual']

 # Create AMPL object
 ampl = AMPL()

 # Set Knitro as the solver to use
 ampl.setOption('solver', 'knitroampl')
 #ampl.setOption('solver', '/scratch/knitroampl')
 #ampl.setOption('knitroampl_auxfiles', 'rc') #print variable/constraint names inside knitroampl
 #ampl.setOption('TMPDIR','/scratch') # write AMPL temporary files to "/scratch"

 # Set AMPL presolve tolerance and option
 ampl.setOption('presolve_eps', 1.0e-8)
 #ampl.setOption('presolve', 0)
 #ampl.setOption('solver_msg', 0)

 # Read the model file                                                  #modelDirectory = argv[2] if argc == 3 else os.path.join('..', 'models')
 usealternate = 0
 useexact = 1
 if usealternate == 0:
     if useexact == 0:
         modelfile = 'constrainedNLP.mod'
         ampl.read('./constrainedNLP.mod')
     else:
         modelfile = 'constrainedNLPexact.mod'
         ampl.read('./constrainedNLPexact.mod')
 else:
     modelfile = 'alternate.mod'
     ampl.read('./alternate.mod')

 log.joint(' ====>>>> AMPL model: ' + modelfile + '\n')

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
  ntrqualswVals[i] = acn.ournontrans[i].swqual
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
  trqualswVals[i] = acn.ourtrans[i].swqual
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

  if acn.ourtrans[i].impedancecorrected==1:
      log.joint('      trans %d numnux %d impedancecorrected %d\n' %(i, acn.ourtrans[i].numnux, acn.ourtrans[i].impedancecorrected))
      if usealternate:
          truegVals[i, -1 + acn.ourtrans[i].minxstf] = 0
          truebVals[i, -1 + acn.ourtrans[i].minxstf] = 0

      for xposition in range(1,1+acn.ourtrans[i].numnux):
          xstf = xposition - 1 + acn.ourtrans[i].minxstf
          truegVals[i,xstf] = acn.ourtrans[i].trueg[xposition]
          truebVals[i,xstf] = acn.ourtrans[i].trueb[xposition]


          log.joint('trans %d xposition %d xstf %d trueg %g trueb %g\n' %(i, xposition, xstf, truegVals[i,xstf], truebVals[i,xstf]))
  else:
      log.joint('      trans %d not impedancecorrected numnux %d\n' %(i, acn.ourtrans[i].numnux))
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
  #log.joint(' switched shunt ' + str(h) + '\n')

  oursh = acn.ourswshunts[h]
  for a in range(1,1+oursh.len):
   bhaVals[h,a] = oursh.b[a]
   nhaVals[h,a] = oursh.n[a]
   #log.joint('    h = %d a = %d b = %f n = %d\n' %(h,a,bhaVals[h,a],nhaVals[h,a]))
 bhaAMPL.setValues(bhaVals)
 nhaAMPL.setValues(nhaVals)
 #breakexit('done switched')
 #amplstate = 'expand; display {j in 1.._nvars} (_varname[j],_var[j].lb0,_var[j].\
 #lb,_var[j].ub,_var[j].ub0);'     #shows full model
 #filename = 'basemodel.out'
 #modelout = ampl.getOutput(amplstate)
 #outfile = open(filename,"w")
 #outfile.write("model = " + str(modelout) + "\n")
 #outfile.close()

 # Set Knitro options
 maxtime = ' maxtime_real=' + str(acn.maxtime) + ' mip_maxtime_real=' + str(acn.maxtime) + ' ms_maxtime_real=' + str(acn.maxtime) # just in case we use multi-start
 #maxtime = ' maxtime_real=600 '
 base_options = 'outlev=4 outmode=0 debug=0 linsolver=6 feastol=1e-5 ftol=1e-3 scale=0 honorbnds=1 cg_maxit=1 bar_murule=1 bar_refinement=1 bar_switchrule=0 bar_feasible=1 restarts=3 restarts_maxit=1000 maxit=30000'  # use bar_initpt=2 bar_directinterval=0 ?
 if acn.useparallel:
     base_options += ' par_numthreads=1'
 just_feasible = 1 #ignoring optimality works quite well
 if just_feasible:
     base_options += ' opttol=1.0e30 opttol_abs=1e30'
 else:
     base_options += ' opttol=1.0e-3'
 if (acn.numbuses < 10000):
     base_options += ' bar_maxcrossit=1'
 mip_options = ' mip_heuristic=0 mip_terminate=1 mip_outinterval=1 mip_outlevel=0 mip_debug=1 mip_outsub=2 mip_nodealg=1 mip_intvar_strategy=0'
 if (acn.numbuses < 100):
     intvar_strategy = ' mip_intvar_strategy=0'  # handle integer vars directly
 else:
     intvar_strategy = ' mip_intvar_strategy=2'  # convert binary vars to complementarity constraints
 if (acn.numbuses < 100):
   relax = ' relax=1'
   log.joint(' using relax = 0 because buses < 100\n')
 else:
   relax = ' relax=1'
   log.joint(' using relax = 1 because buses >= 100\n')
 options = base_options + mip_options + maxtime + relax + intvar_strategy
 #print "options=",options
 ampl.setOption('knitro_options', options)

 # Solve

 usescript = 1
 if usescript:
     acnsolve_throughscript(acn, log, ampl, 'ourobjective', 'initial')
 else:

  try:
      totalcost, feaserror, timeneeded = acnsolveit(acn, log, ampl, 'ourobjective', 'initial')
  except:
      log.joint("Error calling ampl.solve!\n")
      var = traceback.format_exc()
      traceback.print_exc()
      acn_evaluation2.print_alert(var, raise_exception=False)
      pass

 amplstate = 'display binxstf0, xstf0, gf0, bf0;'
 filename = 'structure1.out'
 modelout = ampl.getOutput(amplstate)
 outfile = open(filename,"w")
 outfile.write("objective terms: \n" + str(modelout) + "\n")
 outfile.close()

 # Re-solve with integer vars fixed

 resolve_fixed(acn, ampl, usescript, 1)

 # Finished solve phase

 amplstate = 'display binxstf0, xstf0, gf0, bf0;'
 filename = 'finalstructure.out'
 modelout = ampl.getOutput(amplstate)
 outfile = open(filename,"w")
 outfile.write("objective terms: \n" + str(modelout) + "\n")
 outfile.close()

 acnexposevars(acn, log, ampl, truegVals, truebVals, 'first')
 # Write some AMPL solution files to file
 #amplstate = 'display gencost0, genon0, genstartup0, genshutdown0, genP0, genQ0, minpg, maxpg, minqg, maxqg;'
 amplstate = 'display sum{g in 1..G} Delta*(-gencost0[g] - Oncost[g]*genon0[g]) - sum{g in 1..G} (Sucost[g]*genstartup0[g] + Sdcost[g]*genshutdown0[g]);'
 #amplstate += 'display -sum{i in 1..I} Delta*buscost0[i];'
 amplstate += 'display -sum{i in 1..I} Delta*(sum{n in 1..Pnumcblocks} Pcostcblock[n]*(Pplusblock0[i,n] + Pminusblock0[i,n]));'
 amplstate += 'display -sum{i in 1..I} Delta*(sum{n in 1..Qnumcblocks} Qcostcblock[n]*(Qplusblock0[i,n] + Qminusblock0[i,n]));'
 if mykeybus >= 1:
     amplstate += 'display busQcost0['+str(mykeybus)+'];'
     amplstate += 'display Delta*(sum{n in 1..Qnumcblocks} Qcostcblock[n]*(Qplusblock0['+str(mykeybus)+',n] + Qminusblock0['+str(mykeybus)+',n]));'
 if mykeynontrdual >= 1:
     amplstate += 'display ntrP0['+str(mykeynontrdual)+'];'
 amplstate += 'display sum{l in 1..L: loadisactive[l] > 0} Delta*loadcost0[l];'
 amplstate += 'display -sum{e in 1..NTRcnt_original: ntrqualsw[e]>0 and ntrstat[e]=0} ntrcostsw[e]*ntrON0[e];'
 amplstate += 'display -sum{e in 1..NTRcnt_original: ntrqualsw[e]>0 and ntrstat[e]>0} ntrcostsw[e]*(ntrstat[e]-ntrON0[e]);'
 amplstate += 'display -sum{f in 1..TRcnt_original: trqualsw[f]>0 and trstat[f]=0} trcostsw[f]*trON0[f];'
 amplstate += 'display -sum{f in 1..TRcnt_original: trqualsw[f]>0 and trstat[f]>0} trcostsw[f]*(trstat[f]-trON0[f]);'
 filename = 'amplsol.out'
 modelout = ampl.getOutput(amplstate)
 outfile = open(filename,"w")
 outfile.write("objective terms: \n" + str(modelout) + "\n")
 outfile.close()

 basentrP = numpy.zeros(acn.nontranscount)
 v = ampl.getVariable('ntrP0')
 df = v.getValues()
 vals = df.toDict()
 for i in range(acn.nontranscount):
    basentrP[i] = vals[i+1]


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

 relaxPplus0 = numpy.zeros(acn.numbuses)
 v = ampl.getVariable('relaxPplus0')
 df = v.getValues()
 vals = df.toDict()
 for i in range(acn.numbuses):
    relaxPplus0[i] = vals[i+1]
 #print(relaxPplus0)
 relaxPminus0 = numpy.zeros(acn.numbuses)
 v = ampl.getVariable('relaxPminus0')
 df = v.getValues()
 vals = df.toDict()
 for i in range(acn.numbuses):
    relaxPminus0[i] = vals[i+1]
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


 acnexposevars(acn, log, ampl, truegVals, truebVals, 'last')
 #breakexit('heyyy')
 '''
 for i in range(1,1+acn.numbuses):
  ourbus = acn.ourbuses[i]
  if len(transatbusVals[i]) > 0:
   log.joint('bus ' + str(i) + " p trps:\n")
   for j in range(len(transatbusVals[i])):
      k = transatbusVals[i][j]
      pow = basetrP[k-1]
      log.joint(' ' + str(k) + ' p ' + str(pow) + '\n')
 '''


 basetrQ = numpy.zeros(acn.transcount)
 v = ampl.getVariable('trQ0')
 df = v.getValues()
 vals = df.toDict()
 for i in range(acn.transcount):
    basetrQ[i] = vals[i+1]

 '''
 for i in range(1,1+acn.numbuses):
  ourbus = acn.ourbuses[i]
  if len(transatbusVals[i]) > 0:
   log.joint('bus ' + str(i) + " q trps:\n")
   for j in range(len(transatbusVals[i])):
      k = transatbusVals[i][j]
      qow = basetrQ[k-1]
      log.joint(' ' + str(k) + ' q ' + str(qow) + '\n')
 '''

 #basesw = numpy.zeros(acn.numswitchedshunts)
 basesw={}
 v = ampl.getVariable('xha0')
 df = v.getValues()
 vals = df.toDict()
 for (h,a), step in vals.items():
   #print(h,a,step)
   basesw[h,a] = int(step)

 # Create base case solution data structure and file
 log.joint("Write base case solution\n")
 s1 = acnSolver()
 s1.data.read(raw, sup, con)
 s1.write_sol_base_only('.','solution1',acn,basebusVmags,basebusVangles,clearedloads,basegenP,basegenQ,basegenON,basentrON,basetrON,basetaps,basesw)


 # Evaluate base solution
 try:
  pobj = -1e20
  log.joint("\nEvaluating solution1.txt...\n")
  base.obj, base.infeas, exist, base.summary = acn_evaluation2.base_from_file(acn, "solution1.txt")
  log.joint("\nEvaluating solution_BASECASE.txt...\n")
  pobj, pinfeas, pexist, base.psummary = acn_evaluation2.base_from_file(acn, "solution_BASECASE.txt")
  log.joint("\n------------------------------------------\n")
  if (pexist == True):
    log.joint('infeasible solution: obj=' + str(pobj)+'\n')
    log.joint('infeasible solution: infeas=' + str(pinfeas)+'\n')
  log.joint('solution1 exists=' +  str(exist)+'\n')
  log.joint('solution_BASECASE exists=' +  str(pexist)+'\n')
  if (exist == True):
    log.joint('base.obj=' + str(base.obj)+'\n')
    log.joint('base.infeas=' + str(base.infeas)+'\n')
    if base.infeas == 0 and (pexist == False or base.obj > pobj):
        # We've computed a base case solution that is better than the
        # 'infeasible_solution', so overwrite with this one
        os.popen('\cp solution1.txt solution_BASECASE.txt')
        os.popen('\cp solution1.bin solution_BASECASE.bin')
        # Now write out individual base solution components
        # needed for contingency modeling to files so that they can be read in by
        # MyPython2.py->acnsolverAMPL2.py
        write_base_files(acn, basebusVmags,basebusVangles,ploads,basegenP,basegenON,basegensu,basegensd,basentrON,basetrON, basesw)

 except:
  log.joint("Error evaluating the base solution!\n")
  var = traceback.format_exc()
  traceback.print_exc()
  acn_evaluation2.print_alert(var, raise_exception=False)
  pass

##### Called after initial solve to fix all integer variables and re-solve
#Maybe put this in a different file?

def resolve_fixed(acn, ampl, usescript, solvenum, userelax):

  log = acn.log
  log.joint('running resolve_fixed: usescript=%d\n' %(usescript))

  # if usescript = 1, then we already fixed everything inside the AMPL script after
  # the initial solve, so we can skip the fixing code below

  if usescript == 0:

      # Round and fix all integer variables in the model
      timef1 = time.time()
      log.joint(' fixing tap settings\n')
      numfixed = 0
      # Transformer tap settings
      v = ampl.getVariable('xstf0')
      df = v.getValues()
      vals = df.toDict()
      fixedVals = {}
      for i in range(1,1+acn.transcount_original):
          tmpval = min(vals[i],acn.ourtrans[i].maxxstf)
          tmpval = max(tmpval, acn.ourtrans[i].minxstf)
          tmpval = int(round(tmpval))
          fixedVals[i] = tmpval
          numfixed += 1
      v.setValues(fixedVals)
      v.fix()
      log.joint(' %d fixed\n' %(numfixed))

      log.joint(' fixing shunt susceptances\n')
      numfixed = 0

      # Switched shunt susceptance
      v = ampl.getVariable('xha0')
      df = v.getValues()
      vals = df.toDict()
      fixedVals = {}
      for h in range(1,1 + acn.numswitchedshunts):
          if acn.ourswshunts[h].status > 0:  #switch shunt active
              oursh = acn.ourswshunts[h]
              for a in range(1,1+oursh.len):
                  upbnd = int(oursh.n[a])
                  if upbnd > 0:
                      xha0val = int(vals[h,a])
                      tmpval = min(xha0val,upbnd)
                      tmpval = max(tmpval, 0)
                      fixedVals[h,a] = tmpval
                      #print('h=',h,' a=',a,' xha0val=',xha0val,' upbnd=',upbnd,' fixed=',tmpval)
                      numfixed += 1
      v.setValues(fixedVals)
      v.fix()
      log.joint(' %d fixed\n' %(numfixed))

      log.joint(' fixing line switching variables\n')
      numfixed = 0
      # Line switching
      v = ampl.getVariable('ntrON0')
      df = v.getValues()
      vals = df.toDict()
      fixedVals = {}
      for i in range(1,1+acn.nontranscount_original):
          if vals[i] < 0.5:
              fixedVals[i] = 0
          else:
              fixedVals[i] = 1
          numfixed += 1
      v.setValues(fixedVals)
      v.fix()

      log.joint(' %d fixed\n' %(numfixed))
      log.joint(' fixing transformer switching variables\n')
      numfixed = 0
      # Transformer switching
      v = ampl.getVariable('trON0')
      df = v.getValues()
      vals = df.toDict()
      fixedVals = {}
      for i in range(1,1+acn.transcount_original):
          if vals[i] < 0.5:
              fixedVals[i] = 0
          else:
              fixedVals[i] = 1
          numfixed += 1
      v.setValues(fixedVals)
      v.fix()

      log.joint(' %d fixed\n' %(numfixed))
      log.joint(' fixing generator commitment variables\n')
      numfixed = 0
      # Generator commitment decision variables
      genon = ampl.getVariable('genon0')
      genondf = genon.getValues()
      genonvals = genondf.toDict()
      gensu = ampl.getVariable('genstartup0')
      gensudf = gensu.getValues()
      gensuvals = gensudf.toDict()
      gensd = ampl.getVariable('genshutdown0')
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
              if ourgen.status > 0:  # generator currently on; shut it off
                  fixedsuVals[i] = 0
                  fixedsdVals[i] = 1
              else:
                  fixedsuVals[i] = 0
                  fixedsdVals[i] = 0
              numfixed += 2
          else:
              fixedonVals[i] = 1
              if ourgen.status > 0:
                  fixedsuVals[i] = 0
                  fixedsdVals[i] = 0
              else:                # generator currently off; turn it on
                  fixedsuVals[i] = 1
                  fixedsdVals[i] = 0
              numfixed += 2
      genon.setValues(fixedonVals)
      genon.fix()
      gensu.setValues(fixedsuVals)
      gensu.fix()
      gensd.setValues(fixedsdVals)
      gensd.fix()
      timef2 = time.time()
      log.joint(' %d fixed\n' %(numfixed))
      log.joint(' total time spent fixing %.3e\n' %(timef2 - timef1))

  # Now solve with integer variables fixed

  try:
    elapsedtime = time.time()-acn.timebeg
    remaintime = acn.maxtime - elapsedtime
    remaintime = 600
    logname = ' outname='+'"'+'knitro-fixed'+str(solvenum)+'.log'+'"'
    #print(logname)
    maxtime = ' maxtime_real=' + str(remaintime) + ' mip_maxtime_real=' + str(remaintime) + ' ms_maxtime_real=' + str(remaintime) # just in case we use multi-start
    resolve_options = 'outlev=4 outmode=2 debug=0 feastol=1e-5 feastol_abs=9e-5 ftol=1e-6 scale=0 honorbnds=0 cg_maxit=50 bar_murule=0 bar_feasible=0 bar_refinement=1 bar_initpi_mpec=0.0 maxit=3000 alg=1 strat_warm_start=1 bar_initpt=2 bar_initmu=1e-6 bar_slackboundpush=1e-6 infeastol=1e-5 presolve=1 presolve_tol=0.5 presolve_initpt=1'
    mip_options = ' mip_heuristic=0 mip_terminate=1 mip_outinterval=1 mip_outlevel=0 mip_debug=1 mip_outsub=2 mip_nodealg=1 mip_intvar_strategy=0'
    resolve_options += mip_options + maxtime
    just_feasible = 0
    if just_feasible:
      resolve_options += ' opttol=1.0e30 opttol_abs=1e30'
    else:
      resolve_options += ' opttol=1.0e-3'
    if (acn.numbuses < 10000):
      resolve_options += ' bar_maxcrossit=1'
    resolve_options += logname
    ampl.setOption('knitro_options', resolve_options)

    log.joint(' about to start resolve\n')
    '''
    timea = time.time()
    #ampl.solve()
    #ampl.eval('solve > /dev/null;') # use this to re-direct ampl output to null device
    solveout = ampl.getOutput('solve;')
    log.joint("Finish ampl.solve\n output=%s" %(solveout))
    timeb = time.time()
    totalcost = ampl.getObjective('ourobjective')
    log.joint(" >>> resolve objective value: %g in time %g\n" %(totalcost.value(), timeb-timea))

    '''

    resolveobj, resolvemaxviol, timediff = acnsolveit(acn, log, ampl, 'ourobjective', 'resolve')

    totalgencostvar = ampl.getVariable('totalgencost')
    totalgencost = totalgencostvar.value()
    if userelax:
        totalbuscostvar = ampl.getVariable('totalbuscost')
        totalbuscost = totalbuscostvar.value()
        log.joint('     totalgencost %.4e, totalbuscost %.4e\n' %(totalgencost, totalbuscost))
        totalbusPcostvar = ampl.getVariable('totalbusPcost0')
        totalbusPcost = totalbusPcostvar.value()
        totalbusQcostvar = ampl.getVariable('totalbusQcost0')
        totalbusQcost = totalbusQcostvar.value()
        log.joint('     totalbusPcost %.4e\n' %(totalbusPcost))
        log.joint('     totalbusQcost %.4e\n' %(totalbusQcost))
  except:
    log.joint("Error in fixed re-solve!\n")
    pass
  #breakexit("bbb")

#############

def acnexposevars(acn, log, ampl, truegVals, truebVals, name):

 alldata = acn.alldata

 mykeybus = alldata['mykeybus']
 mykeynontrdual = alldata['mykeynontrdual']
 mykeytrdual = alldata['mykeytrdual']
 log.joint('exposing vars, ' + str(name)+ ', mykeybus = %d, mykeynontrdual = %d mykeytrdual = %d\n' %(mykeybus, mykeynontrdual, mykeytrdual))
 if mykeybus <= 0 and mykeytrdual <= 0:
     log.joint('nothing to expose\n')
     return

 basegenQ = numpy.zeros(acn.numgens)
 v = ampl.getVariable('genQ0')
 df = v.getValues()
 vals = df.toDict()
 for i in range(acn.numgens):
    basegenQ[i] = vals[i+1]


 genQatbus = numpy.zeros(acn.numbuses)
 gqv = ampl.getVariable('genQatbus')
 gqdf = gqv.getValues()
 gqvals = gqdf.toDict()
 for i in range(acn.numbuses):
    genQatbus[i] = gqvals[i+1]

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

 basenontrP = numpy.zeros(acn.nontranscount)
 vP = ampl.getVariable('ntrP0')
 dPf = vP.getValues()
 Pvals = dPf.toDict()
 for i in range(acn.nontranscount):
    basenontrP[i] = Pvals[i+1]

 basenontrQ = numpy.zeros(acn.nontranscount)
 v = ampl.getVariable('ntrQ0')
 df = v.getValues()
 vals = df.toDict()
 for i in range(acn.nontranscount):
    basenontrQ[i] = vals[i+1]


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

 Qfshv = ampl.getVariable('Qfsh')
 Qfshf = Qfshv.getValues()
 Qfshvals = Qfshf.toDict()

 Qswshv = ampl.getVariable('Qswsh')
 Qswshf = Qswshv.getValues()
 Qswshvals = Qswshf.toDict()

 busQcost0v = ampl.getVariable('busQcost0')
 busQcost0f = busQcost0v.getValues()
 busQcost0vals = busQcost0f.toDict()

 relaxQplus0v = ampl.getVariable('relaxQplus0')
 relaxQplus0f = relaxQplus0v.getValues()
 relaxQplus0vals = relaxQplus0f.toDict()
 Qplusblock0v = ampl.getVariable('Qplusblock0')
 Qplusblock0f = Qplusblock0v.getValues()
 Qplusblock0vals = Qplusblock0f.toDict()
 relaxQminus0v = ampl.getVariable('relaxQminus0')
 relaxQminus0f = relaxQminus0v.getValues()
 relaxQminus0vals = relaxQminus0f.toDict()
 Qminusblock0v = ampl.getVariable('Qminusblock0')
 Qminusblock0f = Qminusblock0v.getValues()
 Qminusblock0vals = Qminusblock0f.toDict()

 Qswshv = ampl.getVariable('Qswsh')
 Qswshf = Qswshv.getValues()
 Qswshvals = Qswshf.toDict()
 loadQatbusv = ampl.getVariable('loadQatbus')
 loadQatbusf = loadQatbusv.getValues()
 loadQatbusvals = loadQatbusf.toDict()

 if mykeybus >= 1 and mykeybus <= acn.numbuses:
    thisbus = acn.ourbuses[mykeybus]


    #for n in range(1, acn.numqcblocks+1):
    #    print(n, Qminusblock0vals[mykeybus,n])


    #breakexit('exppp')
    sumnontp = 0
    sumnontq = 0
    for nontrATbus in thisbus.starnontrans.values():
        cnt1 = nontrATbus.count1
        log.joint(' NNN ' + str(cnt1) + ' pf ' + str(basenontrP[cnt1-1]) + ' qf ' + str(basenontrQ[cnt1-1]) + ' parent ' + str(nontrATbus.parentcount1) + '\n')
        sumnontp += basenontrP[cnt1-1]
        sumnontq += basenontrQ[cnt1-1]

    sumtp = 0
    sumtq = 0
    for trATbus in thisbus.startrans.values():
        cnt1 = trATbus.count1
        #print(cnt1, acn.transcount)
        log.joint(' TTT ' + str(cnt1) + ' pf ' + str(basetrP[cnt1-1]) + ' qf ' + str(basetrQ[cnt1-1]) + ' parent ' + str(trATbus.parentcount1) + '\n')

        sumtp += basetrP[cnt1-1]
        sumtq += basetrQ[cnt1-1]
        #breakexit('1')
        trbusi = trATbus.busi0+1
        #breakexit('2')
        trbusj = trATbus.busj0+1
        #breakexit('3')

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
    log.joint(' at %d, sum nontransP  %.4e, sum nontransQ %.4e\n' %(mykeybus, sumnontp, sumnontq))
    log.joint(' at %d, sum transP  %.4e, sum transQ %.4e\n\n' %(mykeybus, sumtp, sumtq))

    #let's again compute genQ
    thisbus = acn.ourbuses[mykeybus]
    gq = 0
    for gc in range(thisbus.gencount):
        gen = thisbus.gens1[gc+1]
        gq += basegenQ[gen.count-1]



    log.joint(' genQ at bus1 ' + str(mykeybus) + ' ' + str(genQatbus[mykeybus]) + ' or ' + str(gq) +'\n')
    log.joint(' Qfsh at bus1 ' + str(mykeybus) + ' ' + str(Qfshvals[mykeybus]) + '\n')
    log.joint(' at %d, sum nontransQ %.4e\n' %(mykeybus, sumnontq))
    log.joint(' at %d, sum transQ %.4e\n' %(mykeybus, sumtq))
    log.joint(' Qswsh at bus1 ' + str(mykeybus) + ' ' + str(Qswshvals[mykeybus]) + '\n')
    log.joint(' loadQ at bus1 ' + str(mykeybus) + ' ' + str(loadQatbusvals[mykeybus]) + '\n')
    net = genQatbus[mykeybus] - sumnontq - sumtq - Qfshvals[mykeybus] - Qswshvals[mykeybus] - loadQatbusvals[mykeybus]

    log.joint(' net at bus1 %d = %.4e\n' %(mykeybus, net))
    log.joint(' relaxQplus0 at bus1 %d = %.4e\n' %(mykeybus, relaxQplus0vals[mykeybus]))
    log.joint(' relaxQminus0 at bus1 %d = %.4e\n' %(mykeybus, relaxQminus0vals[mykeybus]))

    log.joint(' busQcost0 at bus1 ' + str(mykeybus) + ' ' + str(busQcost0vals[mykeybus]) + '\n')

 log.joint('exposed vars\n')

def acnsolve_throughscript(acn, log, ampl, objname, taskname):
  totalobj = 1e20
  feaserror = 1e20
  timediff = 1e20
  try:
    time1 = time.time()
    #ampl.eval('include solve_script;')
    solveout = ampl.getOutput('include solve_script;')
    time2 = time.time()
    timediff = time2-time1
  except:
      log.joint("Error calling ampl solve_script!\n")
      var = traceback.format_exc()
      traceback.print_exc()
      acn_evaluation2.print_alert(var, raise_exception=False)
      pass
  totalobjfun = ampl.getObjective(objname)
  totalobj = totalobjfun.value()
  amplout = 'display '+objname+'.feaserror;'
  feaserrorStr = ampl.getOutput(amplout)
  feaserrorsplit = feaserrorStr.split()[-1]
  feaserror = float(feaserrorsplit)

  log.joint(" >>> solved " + taskname +", objective %.6e feaserror %.3e time %g\n" %(totalobj, feaserror, timediff))


  return totalobj, feaserror, timediff


def acnsolveit(acn, log, ampl, objname, taskname):
     totalobj = 1e20
     feaserror = 1e20
     timediff = 1e20
     division = int(acn.division)
     print("division=",division)
     try:
         time1 = time.time()
         if (division == 2 or division == 4):
            ampl.eval('include resolve_script_div24;')
         else:
            ampl.eval('include resolve_script;')
         #ampl.solve()
         time2 = time.time()
         timediff = time2-time1

     except:
         log.joint("In acnsolveit, error calling ampl.solve!\n")
         var = traceback.format_exc()
         traceback.print_exc()
         breakexit('error')
         pass
     totalobjfun = ampl.getObjective(objname)
     totalobj = totalobjfun.value()
     amplout = 'display '+objname+'.feaserror;'
     feaserrorStr = ampl.getOutput(amplout)
     feaserrorsplit = feaserrorStr.split()[-1]
     #print("feaserror = ",float(feaserrorsplit))
     feaserror = float(feaserrorsplit)

     log.joint(" >>> solved " + taskname +", objective %.6e feaserror %.3e time %g\n" %(totalobj, feaserror, timediff))

     return totalobj, feaserror, timediff
