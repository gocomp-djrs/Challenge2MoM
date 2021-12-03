#!/usr/bin/python

import csv
import sys, os
import traceback
import math
import cmath
from myutils import *
from acnheurbase1 import *
from log import *
import data
import numpy as np
import scipy
from acn import *
#from knitro import *
#from knitroNumPy import *
from amplpy import *
#from line_profiler import *
import time
#import evaluation2
import acn_evaluation2
from acnsolverAMPL import record_xstf
from acnxtra import  acnprintgraph


def trace(frame, event, arg):
    print("%s, %s:%d" % (event, frame.f_code.co_filename, frame.f_lineno))
    return trace

def acncomputepriorbase(alldata):

 log = alldata['log']
 acn = alldata['acn']

 mykeybus = alldata['mykeybus']
 mykeytrdual = alldata['mykeytrdual']

 acn.options['priorsol'] = ACNpriorsol('prior1')
 priorsol = acn.options['priorsol']


 # Get info on nodes/procid
 try:
     #procidenv = os.getenv('SLURM_PROCID',-1)
     procidenv = os.environ.get('SLURM_PROCID')
     if procidenv is None:
         print("No environment variable")
         procid = 0
     else:
         procid = int(procidenv)
 except:
     procid = 0
     if procid < 0:
         procid = 0 # just to be sure
     #log.joint("procid = " +str(procid) + "\n")
     procid_type = type(procid)
     log.joint("procid = " +str(procid) + " type = ", str(procid_type) + "\n")


 log.joint("procid = " +str(procid) +  "\n")
 #breakexit("env")

 priorfilename = "priorbase_" + str(procid) + ".txt"
 log.joint(' computing prior base solution and writing it to ' + priorfilename + '\n')

 prior = open(priorfilename,"w")
 #sys.settrace(trace)



 vmagdict = {}
 vadict = {}
 prior.write("--bus section\n")
 ourbuses = acn.ourbuses
 prior.write("i,v,theta\n")
 if procid == 0:
     fvmag = open('basebusVmags.txt', 'w')
     fvang = open('basebusVangles.txt', 'w')
 for j in range(1,1+acn.numbuses):
   thisbus = ourbuses[j]
   vmag = thisbus.vm
   if vmag > thisbus.nvhi:
      vmag = thisbus.nvhi
   if vmag < thisbus.nvlo:
      vmag = thisbus.nvlo
   va = thisbus.va
   vmagdict[j] = vmag
   vadict[j] = va
   prior.write(str(thisbus.i) + "," + str(vmag) + "," + str(va) + "\n")
   if procid == 0:
       fvmag.write("%s\n" % vmag)
       fvang.write("%s\n" % va)
 if procid == 0:
     fvmag.close()
     fvang.close()

 switchsol = handleswitchedshunts(acn)
 alldata['priorswitchsol'] = switchsol
 #print(switchsol)
 #breakexit('switchsol')

 pgendict = {}
 qgendict = {}
 genondict = {}
 for gencount in range(1, acn.numgens+1):
     gen = acn.ourgens[gencount]
     pgen = qgen = 0.0
     genon = gen.status
     pgen = gen.pg
     qgen = gen.qg
     if pgen > gen.maxpg:
         pgen = gen.maxpg
     if pgen < gen.minpg:
         pgen = gen.minpg
     if qgen > gen.maxqg:
         qgen = gen.maxqg
     if qgen < gen.minqg:
         qgen = gen.minqg

     pgendict[gencount] = pgen
     qgendict[gencount] = qgen
     genondict[gencount] = genon

 pfntrdict = {}
 qfntrdict = {}
 for ntrcount0 in range(acn.nontranscount): #0-based
     thisnontransATbus = acn.nontransATbuses[ntrcount0]
     Gff = thisnontransATbus.Gff
     Gft = thisnontransATbus.Gft
     Bff = thisnontransATbus.Bff
     Bft = thisnontransATbus.Bft
     stat = thisnontransATbus.ournon.status
     parent1 = thisnontransATbus.parentcount1
     busi1 = thisnontransATbus.busi0 + 1
     busj1 = thisnontransATbus.busj0 + 1
     vmi = vmagdict[busi1]
     vmj = vmagdict[busj1]
     vai = vadict[busi1]
     vaj = vadict[busj1]
     vaij = vai - vaj

     #pfntr =  Gff*vmi**2 + (Gft*math.cos(vaij) + Bft*math.sin(vaij))*vmi*vmj )
     pfntr = stat*( Gff*vmi**2 + (Gft*math.cos(vaij) + Bft*math.sin(vaij))*vmi*vmj )
     pfntrdict[ntrcount0+1] = pfntr

     qfntr = stat*( -Bff*vmi**2 + (-Bft*math.cos(vaij) + Gft*math.sin(vaij))*vmi*vmj )
     qfntrdict[ntrcount0+1] = qfntr

     if busi1 == mykeybus:
         log.joint(' ntr1 %d parent1 %d pfntr %.4e\n' %(ntrcount0+1, parent1 ,pfntr))

 xst0fdict = {}
 xpositiondict = {}
 bf0dict = {}
 gf0dict = {}

 #breakexit('2')

 for count in range(1,1+acn.transcount_original):

     tr = acn.ourtrans[count]

     trueb = tr.trueb
     trueg = tr.trueg
     xst0f = 0

     if tr.Ftauvalid:
         xst0f = int( round((tr.tau0f - tr.brevetauf)/tr.taustf))

     if tr.Fthetavalid:
         xst0f = int( round((tr.theta0f - tr.brevethetaf)/tr.thetastf))
         #print(tr.theta0f, tr.brevethetaf, tr.thetastf, tr.thetastf, tr.minxstf, tr.maxxstf)
         #breakexit('')
     if xst0f > tr.maxxstf:
         xst0f = tr.maxxstf
     if xst0f < tr.minxstf:
         xst0f = tr.minxstf
     #if count == 224:
     #    log.joint(' xagain ' + str(xst0f) + ' max ' + str(tr.maxxstf) + '\n')
     #    breakexit('heeee')

     xposition = xst0f + 1 - tr.minxstf

     xst0fdict[count] = xst0f

     xpositiondict[count] = xposition

     #print('count',count,'xst0f',xst0f,'numnux',tr.numnux, 'minxstf',tr.minxstf)
     #print(trueb)
     if tr.impedancecorrected:
         bf0dict[count] = trueb[xposition]
         gf0dict[count] = trueg[xposition]
     else:
         bf0dict[count] = trueb[1]
         gf0dict[count] = trueg[1]

 priorsol.xst0fdict = xst0fdict
 Gf0dict = {}
 Bf0dict = {}
 trpcdict = {}
 angledict = {}

 for trcount0 in range(acn.transcount):

     trATbus = acn.transATbuses[trcount0]
     tr_o_or_r = 2*acn.transATbuses[trcount0].original_or_reversed - 1
     trparent = trATbus.parent
     trpc = trATbus.parentcount1

     bf0 = bf0dict[trpc]
     gf0 = gf0dict[trpc]


     if tr_o_or_r == 1:
         Bf0 = bf0/(trparent.tau0f)**2 + trparent.bmf
         Gf0 = gf0/(trparent.tau0f)**2 + trparent.gmf
     else:
         Bf0 = bf0
         Gf0 = gf0

     busi1 = trATbus.busi0 + 1
     busj1 = trATbus.busj0 + 1

     if busi1 == mykeybus:
         log.joint(' > trc1 ' + str(trcount0+1) + ' trpc ' + str(trpc) + ' o_d ' + str(tr_o_or_r) + '\n')

     vmi = vmagdict[busi1]
     vmj = vmagdict[busj1]
     vai = vadict[busi1]
     vaj = vadict[busj1]
     diff = vai - vaj
     if tr_o_or_r == 1:
         angle = diff - trATbus.parent.theta0f
     else:
         angle = diff + trATbus.parent.theta0f
     angledict[trcount0+1] = angle
     tau = trATbus.parent.tau0f

     if trcount0+1 == mykeytrdual:
      log.joint('    tr ' + str(trcount0+1) + ' busi1 ' + str(busi1) + ' busj1 ' + str(busj1) + ' tr_o_or_r ' + str(tr_o_or_r) + '\n')
      log.joint('    vmi ' + str(vmi) + ' vmj ' + str(vmj) + ' (at ' + str(busi1) +')\n')
      log.joint('    vai ' + str(vai) + ' vaj ' + str(vaj) + '\n')

      log.joint('    angdiff ' + str(diff) + ' theta ' + str(trATbus.parent.theta0f) + ' angle ' + str(angle) + '\n')
      log.joint('    gf0 ' + str(gf0) + ' bf0 ' + str(bf0) + ' tau ' + str(tau) + ' gmf ' + str(trparent.gmf) + '\n')
      log.joint('    Gf0 ' + str(Gf0) + '\n')



     Bf0dict[trcount0+1] = Bf0
     Gf0dict[trcount0+1] = Gf0
     trpcdict[trcount0+1] = trpc

 pftrdict = {}
 qftrdict = {}
 for trcount0 in range(acn.transcount): #0-based
     trATbus = acn.transATbuses[trcount0]

     Gf0 = Gf0dict[trcount0+1]
     Bf0 = Bf0dict[trcount0+1]
     trpc = trpcdict[trcount0+1]
     tr_o_or_r = 2*acn.transATbuses[trcount0].original_or_reversed - 1
     bf0 = bf0dict[trpc]
     gf0 = gf0dict[trpc]

     stat = trATbus.stat

     busi1 = trATbus.busi0 + 1
     busj1 = trATbus.busj0 + 1

     vmi = vmagdict[busi1]
     vmj = vmagdict[busj1]

     angle = angledict[trcount0+1]
     tau = trATbus.parent.tau0f

     pftr = stat*( Gf0*vmi**2 - (gf0*math.cos(angle) + bf0*math.sin(angle))*vmi*vmj/tau )
     qftr = stat*( -Bf0*vmi**2 + (bf0*math.cos(angle) - gf0*math.sin(angle))*vmi*vmj/tau )

     pftrdict[trcount0+1] = pftr
     qftrdict[trcount0+1] = qftr

     if trcount0+1 == mykeytrdual:

         log.joint('  Pterm1 ' + str(Gf0*vmi**2) + '\n')
         log.joint('  Pterm2 ' + str(-gf0*math.cos(angle)) + '\n')
         log.joint('  Pterm3 ' + str(-bf0*math.sin(angle)) + '\n')

         log.joint('  P at ' + str(mykeytrdual) +' = ' + str(pftr) + '\n\n')
         '''
         log.joint('  Qterm1 ' + str(-Bf0*vmi**2) + '\n')
         log.joint('  Qterm2 ' + str(bf0*math.cos(angle)) + '\n')
         log.joint('  Qterm3 ' + str(-gf0*math.sin(angle)) + '\n')
         '''
         log.joint('  Q at ' + str(mykeytrdual) +' = ' + str(qftr) + '\n\n')

 priortloaddict = {}
 prior.write("--load section\n")
 prior.write("i,id,t\n")
 countactive = 0
 if procid == 0:
     f = open('baseloadP.txt', 'w')
 for count in range(1,1+acn.numloads):
     load = acn.ourloads[count]
     pload = 0;
     if load.status:
       priort = 1
       if load.tmin > priort:
           priort = load.tmin
       if load.tmax < priort:
           priort = load.tmax
       pload = load.plpu*priort
       prior.write(str(load.theirload.i) + "," + str(load.theirload.id) + "," + str(priort) + "\n")
       priortloaddict[load.count] = priort
       #if load.i == mykeybus:
       #    breakexit('hhheeeee')
       countactive += 1
     if procid == 0:
         f.write("%s\n" % pload)
 if procid == 0:
     f.close()


 acn_balances(log, acn, pgendict, qgendict, genondict, vmagdict, vadict, pfntrdict, pftrdict,qfntrdict, qftrdict, priortloaddict)

 if procid == 0:
     fgenp = open('basegenP.txt', 'w')
     fgenon = open('basegenON.txt', 'w')
     fgensu = open('basegenSU.txt', 'w')
     fgensd = open('basegenSD.txt', 'w')

 prior.write("--generator section\n")
 prior.write("i,id,p,q,x\n")
 alldata['priorgenON'] = {}
 for gencount in range(1, acn.numgens+1):
     gen = acn.ourgens[gencount]
     pgen = pgendict[gencount]
     qgen = qgendict[gencount]
     genon = genondict[gencount]

     pgen *= genon
     qgen *= genon

     alldata['priorgenON'][gencount] = genon

     #prior.write(str(gen.i) + "," + str(gen.id) + "," + str(sbase*pgen) + "," + str(sbase*qgen) +  "," + str(genon) + "\n")
     prior.write(str(gen.i) + "," + str(gen.id) + "," + str(pgen) + "," + str(qgen) +  "," + str(genon) + "\n")
     startup = 0
     shutdown = 0
     if genon:
         if gen.status == 0:
             startup = 1
     else:
         if gen.status:
             shutdown = 1

     if procid == 0:
         fgenp.write("%s\n" % pgen)
         fgenon.write("%s\n" % genon)
         fgensu.write("%s\n" % startup)
         fgensd.write("%s\n" % shutdown)

 if procid == 0:
     fgenp.close()
     fgenon.close()
     fgensu.close()
     fgensd.close()

 alldata['priorntrON'] = {}

 prior.write("--line section\n")
 prior.write("iorig,idest,id,x\n")
 if procid == 0:
     f = open('basentrON.txt', 'w')
 for j in range(acn.nontranscount_original):
     count = j+1
     ntr = acn.ournontrans[count]
     prior.write(str(ntr.i) + "," + str(ntr.j) + "," + str(ntr.ckt) + "," + str(ntr.status) +  "\n")
     alldata['priorntrON'][count] = ntr.status
     if procid == 0:
         f.write("%s\n" % ntr.status)
 if procid == 0:
     f.close()


 prior.write("--transformer section\n")
 prior.write("iorig,idest,id,x,xst\n")
 if procid == 0:
     f = open('basetrON.txt', 'w')
 alldata['priortrON'] = {}
 for j in range(acn.transcount_original):
     count = j+1
     tr = acn.ourtrans[count]
     xst0f = xst0fdict[count]

     prior.write(str(tr.i) + "," + str(tr.j) + "," + str(tr.ckt) + "," + str(tr.stat) + "," + str(xst0f) + "\n")
     alldata['priortrON'][count] = tr.stat
     if procid == 0:
         f.write("%s\n" % tr.stat)
 if procid == 0:
     f.close()
     record_xstf(acn)  #write backup basexstf.txt file

 #switched shunts
 if procid == 0:
     f = open('basesw.txt', 'w')
 prior.write("--switched shunt section\n")
 prior.write("i,xst1,xst2,xst3,xst4,xst5,xst6,xst7,xst8\n")
 for h in range(1,1 + acn.numswitchedshunts):
  oursh = acn.ourswshunts[h]

  if oursh.status:

    #print(oursh.n)

    prior.write(str(oursh.i))

    if oursh.sizeAh > 0:
        solution = switchsol[h]

    for a in range(1,1+oursh.len):
        if oursh.n[a] > 0:
            if procid == 0:
                f.write("%d %d %d\n" %(h,a,solution[a]))
            prior.write("," + str(solution[a]))
        else:
            break

    prior.write("\n")
 if procid == 0:
     f.close()
 prior.close()

 procedure = alldata['PROCEDURE']
 if procedure == 'COMPUTE_PRIOR_BASE_AND_EVAL' or procedure == 'HEURBASE1':

     base.obj, base.infeas, base.exist, base.summary = acn_evaluation2.base_from_file(acn, priorfilename)
     log.joint(' base.obj = %.8e base.infeas = %.8e \n' %(base.obj, base.infeas))

     #print(acn.pimb_abs)
     #print(acn.ntrmagviol)
     #print(acn.trmagviol)

     #breakexit('preexpand')

     if procedure == 'HEURBASE1':
         acn_expand_prior(acn) #heatmap, totviolvector =

         #breakexit('foooooo')
         #acnprintgraph(acn, 1)

         acnheurbase1(acn, procid)

 return 0

def acnsolverHeurBase1(alldata):

    acn = alldata['acn']
    acn_expand_prior(acn) #heatmap, totviolvector =
    procid = 0
    acnheurbase1(acn, procid)

def acn_expand_prior(acn):
    log = acn.log

    thresh = 1e-3
    log.joint(' \n***expand_prior thresh = ' +str(thresh) + '\n')
    heatmap = acn.heatmap
    totviolvector = acn.totviolvector
    heatmap.fill(0)
    totviolvector.fill(0)

    #heatmap = numpy.zeros(acn.numbuses)
    #totviolvector = numpy.zeros(acn.numbuses)


    threshvector = numpy.full(acn.numbuses, thresh)

    np.add( heatmap, np.greater_equal(acn.pimb_abs, threshvector, out = heatmap))
    np.add( totviolvector, acn.pimb_abs, out = totviolvector )


    np.add( heatmap, np.greater_equal(acn.qimb_abs, threshvector, out = heatmap))
    np.add( totviolvector, acn.qimb_abs, out = totviolvector )


    sumntbviol = 0
    for cnt1 in range(1,1+acn.nontranscount_original):
        thisnontrans = acn.ournontrans[cnt1]
        i0 = thisnontrans.ourfrombus.count-1
        j0 = thisnontrans.ourtobus.count-1
        heatmap[i0] += np.absolute( acn.ntrmagviol[cnt1-1] ) >= thresh
        heatmap[j0] += np.absolute( acn.ntrmagviol[cnt1-1] ) >= thresh

        totviolvector[i0] += np.absolute( acn.ntrmagviol[cnt1-1] )
        totviolvector[j0] += np.absolute( acn.ntrmagviol[cnt1-1] )
        sumntbviol += acn.ntrmagviol[cnt1-1]


    #print(heatmap)
    sumhotviol = 0
    sumpabsviol = 0
    sumqabsviol = 0
    numhotviol = 0
    for j in range(1, 1 + acn.numbuses):
        if heatmap[j-1] > 0:
            #log.joint(' bus1 %d totv %.4e pviol %.4e qviol %.4e\n' %(j-1,  totviolvector[j-1], acn.pimb[j-1], acn.qimb[j-1]))
            sumhotviol += totviolvector[j-1]
            numhotviol += 1
        sumpabsviol += acn.pimb_abs[j-1]
        sumqabsviol += acn.qimb_abs[j-1]
    log.joint(' number of hot violations %d\n' %(numhotviol))
    log.joint(' sum of hot violations %.4e\n' %(sumhotviol))
    log.joint(' sum of ntb violations %.4e\n' %(sumntbviol))
    log.joint(' sum of p abs violations %.4e, as fraction of total abs p load %.4e\n' %(sumpabsviol, sumpabsviol/acn.totalPloadabspu))
    log.joint(' sum of q abs violations %.4e, as fraction of total abs q load %.4e\n' %(sumqabsviol, sumqabsviol/acn.totalQloadabspu))



    #return heatmap, totviolvector

def handleswitchedshunts(acn):
  log = acn.log
  log.joint(' handling switched shunts\n')

  timestart = time.time()

  maxlen = 1
  for h in range(1,1 + acn.numswitchedshunts):
    oursh = acn.ourswshunts[h]
    if oursh.status == 1:
        if maxlen < oursh.sizeAh:
            maxlen = oursh.sizeAh

  log.joint(' max block size ' + str(maxlen) + '\n')

  bvector = np.zeros(maxlen)
  nvector = np.zeros(maxlen, dtype=int)

  # Create AMPL object for twoknap solve
  ampl = AMPL()

  # Set Knitro as the solver to use
  ampl.setOption('solver', 'knitroampl')
  #ampl.setOption('solver', '/scratch/knitroampl')
  #ampl.setOption('knitroampl_auxfiles', 'rc') #print variable/constraint names inside knitroampl
  #ampl.setOption('TMPDIR','/scratch') # write AMPL temporary files to "/scratch"

  # Set AMPL presolve tolerance and option
  ampl.setOption('presolve_eps', 1.0e-8)
  #ampl.setOption('presolve', 0)

  switchsol = {}
  optcount = 0
  cumerror = 0
  for swshindex in range(1,1 + acn.numswitchedshunts):
    oursh = acn.ourswshunts[swshindex]
    if oursh.status == 1:

     #log.joint('  switched shunt ' + str(swshindex) + ' binitpu ' + str(oursh.binitpu) +' stat ' + str(oursh.status) + ' size ' + str(oursh.sizeAh) + ' len ' + str(oursh.len) +'\n')

     if oursh.sizeAh > 0:
         #print(' n',oursh.n)
         #  bvector = np.zeros(oursh.sizeAh)
         #nvector = np.zeros(oursh.sizeAh, dtype=int)
         for a in range(1,1+oursh.len):
             if(oursh.n[a] > 0):
                 #log.joint('    h = %d a = %d b = %f n = %d\n' %(swshindex,a,oursh.b[a],oursh.n[a]))
                 bvector[a-1] = oursh.b[a]
                 nvector[a-1] = oursh.n[a]
         xval, error = twoknap(log, ampl, oursh.sizeAh, oursh.binitpu, bvector, nvector, swshindex)
         cumerror += error
         optcount += 1
         switchsol[swshindex] = xval
  timeend = time.time()
  log.joint(' total time spent on %d 2knaps:  %f; cum error: %.4e\n' %(optcount, timeend - timestart, cumerror))
  #breakexit(' st')
  return switchsol

def twoknap(log, ampl, Asize, binitpu, bvector, nvector, problemindex):
    #log.joint('   2k\n')

    #print('   ',binitpu, bvector, nvector, Asize)

    np.set_printoptions(precision=16)

    ampl.eval('reset;')

    # Create AMPL object
    #ampl = AMPL()

    # Set Knitro as the solver to use
    #ampl.setOption('solver', 'knitroampl')
    #ampl.setOption('knitroampl_auxfiles', 'rc') #print variable/constraint names inside knitroampl

    # Set AMPL presolve tolerance and option
    #ampl.setOption('presolve_eps', 1.0e-8)
    #ampl.setOption('presolve', 0)

    # The stuff above was pulled outside this routine so it is only
    # done once, instead of each time through the loop

    try:
        ampl.read('./2knap.mod')
    except:
        log.joint('cannot find 2knap.mod\n')
        #breakexit('2knapfailure')

    AsizeAMPL = ampl.getParameter('Asize')
    AsizeAMPL.set(Asize)

    binitAMPL = ampl.getParameter('binit')
    binitAMPL.set(binitpu)


    nAMPL = ampl.getParameter('n')
    bAMPL = ampl.getParameter('b')
    n = {}
    b = {}
    for j in range(Asize):
        n[j+1] = nvector[j]
        b[j+1] = bvector[j]

    #print(n,b)
    nAMPL.setValues(n)
    bAMPL.setValues(b)
    #amplstate = 'expand; display {j in 1.._nvars} (_varname[j],_var[j].lb0,_var[j].\
    #lb,_var[j].ub,_var[j].ub0);'     #shows full model
    #filename = '2kmodel.out'
    #modelout = ampl.getOutput(amplstate)
    #outfile = open(filename,"w")
    #outfile.write("model = " + str(modelout) + "\n")
    #outfile.close()

    maxtime = ' maxtime_real=600 '
    base_options = 'outlev=0 outmode=0 debug=0 linsolver=6 feastol=1e-5 ftol=1e-3 scale=0 honorbnds=1 cg_maxit=1 bar_murule=1 bar_refinement=1 bar_maxcrossit=1 bar_switchrule=0 bar_feasible=1 restarts=3 maxit=3000'
    #mip_options = ' mip_heuristic=0 mip_terminate=1 mip_outinterval=1 mip_outlevel=0 mip_debug=0 mip_outsub=0 mip_nodealg=1 mip_intvar_strategy=0'
    mip_options = ' mip_heuristic=0 mip_integral_gap_rel=1e-4 mip_outinterval=1 mip_outlevel=0 mip_debug=0 mip_outsub=0 mip_intvar_strategy=0'

    options = base_options + mip_options + maxtime
    ampl.setOption('knitro_options', options)


    # Solve
    try:
        time1 = time.time()
        #sys.stdout = open('stdout.log', 'a')

        sys.stdout =  open(os.devnull, 'w')
        ampl.solve()
        sys.stdout = sys.__stdout__

        #ampl.eval('solve > junk.txt;')  # use this to re-direct ampl output to junk file
        #ampl.eval('solve > /dev/null;') # use this to re-direct ampl output to null device
        time2 = time.time()
    except:
        log.joint("Error calling ampl.solve!\n")
        var = traceback.format_exc()
        traceback.print_exc()
        pass

    #breakexit('solve')
    #print('   ', Asize,'>',binitpu, bvector, nvector)
    sval = ampl.getVariable('s')
    df = sval.getValues()
    vals = df.toDict()
    totalerror = vals[1] + vals[2]
    if totalerror > 1e-02 or 1 == problemindex%100:
        log.joint(' total error in 2k ' + str(problemindex) +': ' + str(totalerror) + ' (size: ' + str(Asize) + ')\n')

    xval = ampl.getVariable('x')
    df = xval.getValues()
    xvals = df.toDict()

    if totalerror > 9e-4:#1e-4:
        log.joint(' total error %.3e\n' %(totalerror))
        #print(xvals)
        log.joint(' weak knapsack solution for problem ' + str(problemindex) + '\n')
        log.joint(' (rhs) binit ' + str(binitpu) +'\n')
        log.joint(' s[1] = ' +str(vals[1])+'\n')
        log.joint(' s[2] = ' +str(vals[2])+'\n')
        for j in range(1,1+Asize):
            log.joint(' x['+str(j)+'] = '+str(xvals[j])+'\n')

        print('upper bounds on integer variables, n = ', n)
        print('knapsack coefficients, b = ', b)
        breakexit('huh????')


    '''
    for j in range(Asize):
        log.joint('x['+str(j)+'] = ' + str(xvals[j+1]) +'\n')

    '''

    return xvals, totalerror
    #breakexit(' ran 2k')

def acn_balances(log, acn, pgendict, qgendict, genondict, vmagdict, vadict, pfntrdict, pftrdict,qfntrdict, qftrdict, priortloaddict):

    alldata=acn.alldata
    mykeybus = alldata['mykeybus']
    mykeytrdual = alldata['mykeytrdual']

    ourbuses = acn.ourbuses
    pgenoutput = np.zeros(1+acn.numbuses)
    pfxsh = np.zeros(1+acn.numbuses)
    pntrinj = np.zeros(1+acn.numbuses)
    ptrinj = np.zeros(1+acn.numbuses)
    pnetnotload = np.zeros(1+acn.numbuses)

    qgenoutput = np.zeros(1+acn.numbuses)
    qfxsh = np.zeros(1+acn.numbuses)
    qntrinj = np.zeros(1+acn.numbuses)
    qtrinj = np.zeros(1+acn.numbuses)
    qnetnotload = np.zeros(1+acn.numbuses)

    ploadbus = np.zeros(1+acn.numbuses)
    qloadbus = np.zeros(1+acn.numbuses)

    for j in range(1,1+acn.numbuses):
        thisbus = ourbuses[j]
        for k in range(1,1+thisbus.gencount):
            thisgen = thisbus.gens1[k]
            thisgencount = thisgen.count
            #if j == mykeybus:
            #    print(j,k,thisgencount,pgendict[thisgencount])

            pgenoutput[j] += pgendict[thisgencount]*genondict[thisgencount]
            qgenoutput[j] += qgendict[thisgencount]*genondict[thisgencount]
        pnetnotload[j] += pgenoutput[j]
        qnetnotload[j] += qgenoutput[j]

    if mykeybus >= 1:
        log.joint('vmag at %d = %.3e\n' %(mykeybus, acn.ourbuses[mykeybus].vm))
        log.joint('pgen at %d = %.4e\n' %(mykeybus, pgenoutput[mykeybus]))
        log.joint('qgen at %d = %.4e\n' %(mykeybus, qgenoutput[mykeybus]))


    for j in range(1,1+acn.numbuses):
        thisbus = ourbuses[j]
        pfxshuntval = thisbus.fsglpu*vmagdict[j]**2
        qfxshuntval = thisbus.fsblpu*vmagdict[j]**2
        #if j == mykeybus:
        #    print(j,pfxshuntval)
        pfxsh[j] += pfxshuntval
        pnetnotload[j] -= pfxsh[j]
        qfxsh[j] += qfxshuntval
        qnetnotload[j] -= qfxsh[j]

    if mykeybus >= 1:
        log.joint('pfxsh at %d = %.4e\n' %(mykeybus, pfxsh[mykeybus]))
        log.joint('qfxsh at %d = %.4e\n' %(mykeybus, qfxsh[mykeybus]))



    for j in range(1,1+acn.numbuses):
        thisbus = ourbuses[j]

        for ournontratbus in thisbus.starnontrans.values():
            cnt1 = ournontratbus.count1
            pntrinj[j] += pfntrdict[cnt1]
            qntrinj[j] += qfntrdict[cnt1]
        pnetnotload[j] -= pntrinj[j]
        qnetnotload[j] -= qntrinj[j]
            #if j == mykeybus:
            #   print('nnn',j, cnt1, pfntrdict[cnt1])
        for trATbus in thisbus.startrans.values():
            cnt1 = trATbus.count1
            ptrinj[j] += pftrdict[cnt1]
            qtrinj[j] += qftrdict[cnt1]
        pnetnotload[j] -= ptrinj[j]
        qnetnotload[j] -= qtrinj[j]
        for k in range(thisbus.loadcount):
            thisload = thisbus.loads[k]
            if thisload.status:
                ploadbus[j] += thisload.plpu*priortloaddict[thisload.count]
                qloadbus[j] += thisload.qlpu*priortloaddict[thisload.count]
        pnetnotload[j] -= ploadbus[j]
        qnetnotload[j] -= qloadbus[j]

    if mykeybus >= 1:
        log.joint('pntrinj at %d = %.4e\n' %(mykeybus, pntrinj[mykeybus]))
        log.joint('ptrinj at %d = %.4e\n' %(mykeybus, ptrinj[mykeybus]))
        log.joint('pload at %d = %.4e\n' %(mykeybus, ploadbus[mykeybus]))
        log.joint('pnetnotload at %d = %.4e\n' %(mykeybus, pnetnotload[mykeybus]))

        log.joint('qntrinj at %d = %.4e\n' %(mykeybus, qntrinj[mykeybus]))
        log.joint('qtrinj at %d = %.4e\n' %(mykeybus, qtrinj[mykeybus]))
        log.joint('qload at %d = %.4e\n' %(mykeybus, qloadbus[mykeybus]))
        log.joint('qnetnotload at %d = %.4e\n' %(mykeybus, qnetnotload[mykeybus]))


        thisbus = acn.ourbuses[mykeybus]
        for k in range(thisbus.loadcount):
            thisload = thisbus.loads[k]
            log.joint('load # ' + str(k) + ' (at this bus)')
            thisload.showload(log)
