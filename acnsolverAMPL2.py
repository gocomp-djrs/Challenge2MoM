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
from acnsolveconthruH2 import *
from acn_solution import *

def trace(frame, event, arg):
    print("%s, %s:%d" % (event, frame.f_code.co_filename, frame.f_lineno))
    return trace
def getsw(filename):
    sw = {}
    with open(filename, 'r') as f:
        fcontents = f.readlines()
        for line in fcontents:
            data = line.split()
            sw[int(data[0]), int(data[1])] = int(data[2])
    f.close()
    return sw

def file_to_dict(file, isint):
    x = {}
    i = 1
    with open(file, 'r') as f:
        fcontents = f.readlines()
        for line in fcontents:
            # remove linebreak which is the last character of the string
            if isint:
                x[i] = int(line[:-1])
            else:
                x[i] = float(line[:-1]) #does this give good precision?
            i += 1
    f.close()
    #print('x=',x)
    return x

def write_ctg_solution(log, acn, ampl, genoutIndex, genoutID, lineout_i, lineout_j, lineout_ckt, ctgytmpfile):

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



def acnsolverAMPL2(alldata):

 log = alldata['log']
 log.joint("\n***\n solving contingencies\n")

 #sys.settrace(trace)

 acn = alldata['acn']
 rawname = alldata['RAW']
 conname = alldata['CON']
 supname = alldata['JSON']

 pi = math.pi
 #acn.log = log

 numpy.set_printoptions(precision=16)
 elapsedtime = time.time()-acn.timebeg
 remaintime = acn.maxtime - elapsedtime

 # Get some info on contingencies
 contingencies = acn.con.contingencies
 numctgys = len(contingencies)
 log.joint("\n***\n probing %i contingencies\n\n" %(numctgys))

 # Get info on cores and processors
 numcores = multiprocessing.cpu_count()
 print('numcores=',numcores)
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
     #breakexit("env")
     cores = open("cores"+str(procid)+".txt","w",0)
     cores.write("num cores = "	+ str(numcores) + "\n")
     cores.write("procid = "	+ str(procid) + "\n")
     cores.close()

 basebusVmags = file_to_dict('basebusVmags.txt', 0)
 basebusVangles = file_to_dict('basebusVangles.txt', 0)
 acn.basebusVmags = numpy.zeros(acn.numbuses)
 acn.basebusVangles = numpy.zeros(acn.numbuses)
 for i in range(acn.numbuses):
     acn.basebusVmags[i] = basebusVmags[i+1]
     acn.basebusVangles[i] = basebusVangles[i+1]

 # Check if we have finished writing fallback contingency files
 if procid == 0:
     pname = "write_ctgs.log"
     try:
         tmp = open(pname,"r")
         tmp.close()
         # files written
     except:
         # Need to write fallback contingency files
         s = acnSolver()
         s.data.read(rawname, supname, conname)
         write_ctgs = True
         ctimebeg = time.time()
         if write_ctgs:
            print("Writing contingencies...")
            #s.write_sol('./')
            #print("WRITING NEW FALLBACK CONTINGENCIES!")
            s.write_sol_ctgs_only_new('.','solution_BASECASE',acn,acn.basebusVmags,acn.basebusVangles)
            #s.write_sol_ctgs_only('.','solution_BASECASE')
         ctime = time.time()-ctimebeg
         print("code2: write contingencies time = ",ctime)

 # Read in and store arrays needed from base solution
 basentrON = file_to_dict('basentrON.txt', 1)
 basetrON = file_to_dict('basetrON.txt', 1)
 baseloadP = file_to_dict('baseloadP.txt',0)
 basegenP = file_to_dict('basegenP.txt',0)
 basegenON = file_to_dict('basegenON.txt', 1)
 basegenSU = file_to_dict('basegenSU.txt', 1)
 basegenSD = file_to_dict('basegenSD.txt', 1)

 basexstf = file_to_dict('basexstf.txt', 1)
 #print(basexstf)


 #switched shunts stored differently
 basesw = getsw('basesw.txt')

 #print(basesw)



 # TBD: sort contingencies
 conlist = []
 sortedind = []
 for ctgidx, (conlabel, ctg) in enumerate(contingencies.items(), start=1):
     conlist.append(conlabel)
 conlabelsorted = []
 contingencysortedorder = [] #tbd: determine this order
 for h in range(numctgys):
     contingencysortedorder.append(h+1)
     j = contingencysortedorder[h]
     conlabelsorted.append(conlist[j-1])
     sortedind.append(h)
     #log.joint('%i %i %s\n' %(sortedind[h],j,conlabelsorted[h]))
 #breakexit('sort')

 # Run contingency loop
 ctgnum = 0
 first = 1
 numnodes = acn.numnodes
 useparallel = 0
 acn.useparallel = useparallel

 # use this to suppress stdout in parallel contingency loop
 nooutput = 0
 if nooutput:
     sys.stdout =  open(os.devnull, 'w')
     sys.stderr =  open(os.devnull, 'w')

 # Initialize
 max_passes = 1
 passnum = 1
 numnegscores = 0
 stop = 0
 acn.scores = np.zeros(numctgys)

 while stop == 0:
   if useparallel == 0:
     log.joint('---> sequential solves\n')
     # Sequential implementation of contingency loop
     #for conlabel, ctg in contingencies.items():
     for h in range(numctgys):
     #for ctgidx, (conlabel, ctg) in enumerate(contingencies.items(), start=1):
         conlabel = conlabelsorted[h]
         ctg = contingencies[conlabel]
         ctgidx = contingencysortedorder[h]
         if 1 or conlabel == 'CTG_000551': #pick out specific contingency for debugging
             #print("conlabel=",conlabel," h=",h," numctgys=",numctgys," passnum=",passnum)
             #print("ctgidx=",ctgidx)
             #breakexit("before ctgsolve")
             ctgsolve(alldata, passnum, h+1, ctgidx, conlabel, ctg, contingencies, basebusVmags, basebusVangles, basentrON, basetrON, baseloadP, basegenP, basegenON, basegenSU, basegenSD, basesw, basexstf)
             #breakexit("after ctgsolve")
             #sys.exit()
   else:
     # Parallel implementation of contingency loop using forks
     # Currently fork as many child processes as threads/cores available
     log.joint('---> parallel solves\n')
     numleft = numctgys - first + 1
     numpercore = int(math.ceil(float(numleft)/float(numcores)))
     numpernode = int(math.ceil(float(numleft)/float(numnodes)))
     #numpercore = numctgys/numcores
     # determine range of contingencies to be handled by this procid
     proclo = procid*numpernode
     prochi = (procid+1)*numpernode
     numconlist2 = prochi-proclo
     for i in range(numcores):
         try:
             pid = os.fork()
         except:
             #print("Could not create a child process")
             continue
         if pid == 0:
             #for ctgidx, (conlabel, ctg) in enumerate(contingencies.items(), start=1):
             for h in range(numctgys):
                 #if passnum > 1:
                 #    print("!!!!!!CURRENT SCORE: ",acn.scores[h])
                 if sortedind[h] >= proclo and sortedind[h] < prochi:
                     conlabel = conlabelsorted[h]
                     ctg = contingencies[conlabel]
                     ctgidx = contingencysortedorder[h]
                     modulo = (h+1)%numcores
                     elapsedtime = time.time()-acn.timebeg
                     remaintime = acn.maxtime - elapsedtime
                     needtosolve = 1
                     if passnum > 2 and numnegscores > 0 and acn.scores[h] >= 0:
                         needtosolve = 0 #if negative scores, focus on improving only these for now
                     limit = max(180,0.1*numctgys)
                     #if i==modulo and ctgidx >= first and needtosolve == 1:  # ignore time limit
                     if i==modulo and ctgidx >= first and needtosolve == 1 and remaintime >= limit:
                     #if i==modulo and ctgidx >= procid*numpernode + first and ctgidx < (procid+1)*numpernode + first and remaintime >= limit:
                         try:
                             log.joint("passnum="+str(passnum)+" conlabel="+str(conlabel)+"\n");
                             print("Solve contingency ",ctgidx, " / ",conlabel, " with procid ", procid, " core ", i, " remain_time=", remaintime)
                             log.joint("Solve contingency "+str(ctgidx)+" with procid "+str(procid)+" core "+str(i)+" remain_time="+str(remaintime)+"\n")
                             ctgsolve(alldata, passnum, h+1, ctgidx, conlabel, ctg, contingencies, basebusVmags, basebusVangles, basentrON, basetrON, baseloadP, basegenP, basegenON, basegenSU, basegenSD, basesw, basexstf)
                         except:
                             log.joint("Contingency "+str(ctgidx)+" encountered a problem inside ctgsolve!\n")
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

   # SORT AND THEN MAKE ANOTHER PASS
   conlist2 = conlabelsorted
   acn.scores = np.zeros(numctgys)
   numnegscores = 0
   for h in range(numctgys):
     conlabel = conlabelsorted[h]
     pname = "score1_"+str(conlabel)+".txt"
     try:
        with open(pname) as f:
            line = f.readline()
            acn.scores[h] = float(line[:-1])
        f.close()
     except:
        acn.scores[h] = -1e20
     if acn.scores[h] < 0:
        numnegscores += 1

     #score = acn.scores[h]
     #print("ctg label: ",conlabel, "score=",acn.scores[h])

   #print("scores=",acn.scores)
   sortedind = numpy.argsort(acn.scores)
   #print("sorted indices:",sortedind)
   conlist2sorted = []
   for h in range(numctgys):
     ind = sortedind[h]
     conlist2sorted.append(conlist2[ind])
   #print("sorted labels 1:",conlist2sorted)
   conlabelsorted = conlist2sorted
   #print("sorted labels 1:",conlabelsorted)
   for h in range(numctgys):
     contingencysortedorder[h] = sortedind[h]+1
   #print("order:",contingencysortedorder)
   #breakexit("scores")

   # Decide whether to make another pass
   passnum += 1
   stop = 0
   elapsedtime = time.time()-acn.timebeg
   remaintime = acn.maxtime - elapsedtime
   if remaintime < 30:
     stop = 1
   if passnum > max_passes:
     stop = 1
   #stop = 0 # use to force stop

 # CONTINGENCY SOLVE LOOP FINISHED

 if nooutput:
     sys.stdout = sys.__stdout__
     sys.stderr = sys.__stderr__
 #breakexit('sol')

def ctgsolve(alldata, passnum, ctgcount, ctgidx, conlabel, ctg, contingencies, basebusVmags, basebusVangles, basentrON, basetrON, baseloadP, basegenP, basegenON, basegenSU, basegenSD, basesw, basexstf):
#def ctgsolve(acn, passnum, ctgcount, ctgidx, conlabel, ctg, contingencies, basebusVmags, basebusVangles, basebusCshunts, basegenP, basegenQ, baselineP, baselineQ, busbaseP, busbaseQ):

    log = alldata['log']
    acn = alldata['acn']
    numctgys = len(contingencies)

    log.joint("\n***\n Enter ctgsolve to solve contingency: %s\n" %(conlabel))
    log.joint("*******************************************\n")
    log.joint("     Solving contingency %i of %i\n" %(ctgidx, numctgys))
    log.joint("*******************************************\n")

    # passnum used to indicate the pass number through the contingency loop

    skip = 0
    boe = ctg.branch_out_events
    goe = ctg.generator_out_events
    if boe:
        b = boe[0] # only one event
        i = b.i
        j = b.j
        ckt = b.ckt
        #log.joint(" branch out event -> i = %i j = %i ckt = %s\n" %(i,j,ckt))
    elif goe:
        g = goe[0] # only one event
        genoutIndex = g.i
        genoutID = g.id
        #log.joint(" generator out event -> i = %i ID = %s\n" %(g.i, g.id))
    else:
        log.joint(' ----> contingency is neither branch-out nor gen-out\n')
        # just use base solution in this case
        skip = 1
        #return

    # Compute score using base solution
    # If this score looks good we may skip the contingency solve

    # For large models, skip the solve if base solution seems reasonably good to save time

    # If less than x minutes left, abandon contingency solves and just copy base
    # solution before we run out of time

    if skip == 1:
        #tbd
        log.joint("Skipping contingency \n")
    else:
        #log.joint("Solving contingency \n")
        if boe:
            b = boe[0] # only one event
            i = b.i
            j = b.j
            ckt = b.ckt
            acnsolver_lineout(alldata, passnum, ctgcount, ctgidx, conlabel, ctg, contingencies, basebusVmags, basebusVangles, basentrON, basetrON, baseloadP, basegenP, basegenON, basegenSU, basegenSD, basesw, basexstf)
        elif goe:
            g = goe[0] # only one event
            genoutIndex = g.i
            genoutID = g.id
            acnsolver_genout(alldata, passnum, ctgcount, ctgidx, conlabel, ctg, contingencies, basebusVmags, basebusVangles, basentrON, basetrON, baseloadP, basegenP, basegenON, basegenSU, basegenSD, basesw, basexstf)
        else:
            log.joint(' ----> contingency is neither branch-out nor gen-out\n')
            # just use base solution in this case (should be caught already)
            skip = 1
            #return

        elapsedtime = time.time()-acn.timebeg
        remaintime = acn.maxtime - elapsedtime
        #if remaintime < 30.0:
        #   log.joint("Skipping writing contingency because running out of time: time left = "+ str(remaintime) + "\n")
        #   skip = 1

        #if skip < 1:
           # Now compute score for contingency solution and check feasibility
           # tbd

def acnsolver_lineout(alldata, passnum, ctgcount, ctgidx, conlabel, ctg, contingencies, basebusVmags, basebusVangles, basentrON, basetrON, baseloadP, basegenP, basegenON, basegenSU, basegenSD, basesw, basexstf):

    log = alldata['log']
    acn = alldata['acn']
    #log.joint("\n***\n Enter acn_lineout for contingency: %s\n" %(conlabel))

    acnsolverAMPLcon(alldata, passnum, ctgcount, ctgidx, conlabel, ctg, contingencies, basebusVmags, basebusVangles, basentrON, basetrON, baseloadP, basegenP, basegenON, basegenSU, basegenSD, basesw, basexstf)
    #return

def acnsolver_genout(alldata, passnum, ctgcount, ctgidx, conlabel, ctg, contingencies, basebusVmags, basebusVangles, basentrON, basetrON, baseloadP, basegenP, basegenON, basegenSU, basegenSD, basesw, basexstf):

    log = alldata['log']
    acn = alldata['acn']
    #log.joint("\n***\n Enter acn_genout for contingency: %s\n" %(conlabel))

    acnsolverAMPLcon(alldata, passnum, ctgcount, ctgidx, conlabel, ctg, contingencies, basebusVmags, basebusVangles, basentrON, basetrON, baseloadP, basegenP, basegenON, basegenSU, basegenSD, basesw, basexstf)
    #return


# the bulk of the code for setting up and solving individual contingencies is in acnsolverAMPLcon
def acnsolverAMPLcon(alldata, passnum, ctgcount, ctgidx, conlabel, ctg, contingencies, basebusVmags, basebusVangles, basentrON, basetrON, baseloadP, basegenP, basegenON, basegenSU, basegenSD, basesw, basexstf):

 log  = alldata['log']
 #log.joint("\n***\n Enter acnsolverAMPLcon for contingency: %s\n" %(conlabel))


 #sys.settrace(trace)

 acn = alldata['acn']
 rawname = alldata['RAW']
 conname = alldata['CON']
 supname = alldata['JSON']
 division = int(acn.division)

 pi = math.pi

 numpy.set_printoptions(precision=16)
 elapsedtime = time.time()-acn.timebeg
 remaintime = acn.maxtime - elapsedtime

 # Mark info for contingency
 genoutIndex = -1
 genoutID = -1
 lineout_i = -1
 lineout_j = -1
 lineout_ckt = -1
 boe = ctg.branch_out_events
 goe = ctg.generator_out_events
 if boe:
     b = boe[0] # only one event
     lineout_i = b.i
     lineout_j = b.j
     lineout_ckt = b.ckt
     #log.joint(" branch out event -> i = %i j = %i ckt = %s\n" %(b.i,b.j,b.ckt))
 elif goe:
     g = goe[0] # only one event
     genoutIndex = g.i
     genoutID = g.id
     #log.joint(" generator out event -> i = %i ID = %s\n" %(g.i, g.id))
 else:
     log.joint(' ----> contingency is neither branch-out nor gen-out\n')


 genout_index = -1
 for gencount in range(1, acn.numgens+1):
     gen = acn.ourgens[gencount]
     if gen.i == genoutIndex and gen.id == genoutID:
         genout_index = gencount
 ntrout_index = -1
 for ntrcount in range(1, acn.nontranscount_original+1):
     ntr = acn.ournontrans[ntrcount]
     if ntr.i == lineout_i and ntr.j == lineout_j and ntr.ckt == lineout_ckt:
         ntrout_index = ntrcount
 trout_index = -1
 for trcount in range(1, acn.transcount_original+1):
     tr = acn.ourtrans[trcount]
     if tr.i == lineout_i and tr.j == lineout_j and tr.ckt == lineout_ckt:
         trout_index = trcount

 if alldata['FLAVOR'] == 'H2':

     log.joint(' FLAVOR = H2')
     count1_of_affected_bus = count1_of_affected_from_bus = count1_of_affected_to_bus = -1
     if genout_index > 0:
         count1_of_affected_bus = acn.ourgens[genout_index].ourbus.count
         log.joint('  generator count1 %d at bus count1 %d is out\n' %(genout_index, count1_of_affected_bus))
     if ntrout_index > 0:
         count1_of_affected_from_bus = acn.ournontrans[ntrout_index].ourfrombus.count
         count1_of_affected_to_bus = acn.ournontrans[ntrout_index].ourtobus.count
         log.joint('  nontrans count1 %d at ( %d , %d ) is out\n' %(ntrout_index,count1_of_affected_from_bus, count1_of_affected_to_bus))
     if trout_index > 0:
         count1_of_affected_from_bus = acn.ourtrans[trout_index].ourfrombus.count
         count1_of_affected_to_bus = acn.ourtrans[trout_index].ourtobus.count
         log.joint('  trans count1 %d at ( %d , %d ) is out\n' %(trout_index,count1_of_affected_from_bus, count1_of_affected_to_bus))

     #breakexit('H2,1')

     exceptions = {}
     numexceptions = 0

     use_exceptions = 1
     alsohitgens = 1

     if acn.numbuses > 5000:
         maxexceptions = 2
     elif acn.numbuses > 1000:
         maxexceptions = 5
     else:
         maxexceptions = 5 #10

     if division < 3:
         maxexceptions = maxexceptions * passnum
     elif passnum > 1:
         maxexceptions = maxexceptions * (passnum-1)

     freeze = fixall = 1
     #print("division=",acn.division," passnum=",passnum)
     #breakexit("passnum")

     if genout_index > 0:
         #one dijkstra
         numexceptions = acn_dijkstra(acn, ctg.label, fixall, ntrout_index - 1, trout_index - 1, count1_of_affected_bus, exceptions, numexceptions, maxexceptions, alsohitgens)

     if ntrout_index > 0 or trout_index > 0:
         #two dijkstras
         numexceptions = acn_dijkstra(acn, ctg.label, fixall, ntrout_index - 1, trout_index - 1, count1_of_affected_from_bus, exceptions, numexceptions, maxexceptions, alsohitgens)


         numexceptions = acn_dijkstra(acn, ctg.label, fixall, ntrout_index - 1, trout_index - 1, count1_of_affected_to_bus, exceptions, numexceptions, maxexceptions, alsohitgens)


     for j in range(1,1+acn.numbuses):
      thisbus = acn.ourbuses[j]
      if thisbus.count in exceptions.values():
          log.joint(' bus count1 %d is excepted\n' %(j))


     acnsolveconthruH2(alldata,acn, log, passnum, ctgcount, ctgidx, conlabel, ctg, basebusVmags, basebusVangles, basentrON, basetrON, baseloadP, basegenP, basegenON, basegenSU, basegenSD, basesw, basexstf, genout_index, ntrout_index, trout_index, genoutIndex, genoutID, lineout_i, lineout_j, lineout_ckt, exceptions, numexceptions)
     return


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

 # Read the model file
 #modelDirectory = argv[2] if argc == 3 else os.path.join('..', 'models')
 #ampl.read('./constrainedNLP2.mod')
 ampl.read('./constrainedNLP2exact.mod')

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


 for gencount in range(1, acn.numgens+1):
  ourgen = acn.ourgens[gencount]
  buscount = acn.busnumbertocount[ourgen.i]

  ourbus = acn.ourbuses[buscount]
  numcb[gencount] = ourgen.numcblocks
  oncost[gencount] = ourgen.oncost #well, could have done this to begin with

  sucost[gencount] = ourgen.sucost
  sdcost[gencount] = ourgen.sdcost
  if ourgen.i == genoutIndex and ourgen.id == genoutID:
      # Generator is out.  How to force genon=0?
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
  ntrqualswVals[i] = acn.ournontrans[i].swqual
  ntrcostswVals[i] = acn.ournontrans[i].csw
  #ntrstatVals[i] = acn.ournontrans[i].status
  ntrratingVals[i] = acn.ournontrans[i].rating_c #use rating_c for ctgy
 ntrqualswAMPL.setValues(ntrqualswVals)
 ntrcostswAMPL.setValues(ntrcostswVals)
 #ntrstatAMPL.setValues(ntrstatVals)
 ntrstatAMPL.setValues(basentrON)
 ntrratingAMPL.setValues(ntrratingVals)

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
  trqualswVals[i] = acn.ourtrans[i].swqual
  trcostswVals[i] = acn.ourtrans[i].csw
  #trstatVals[i] = acn.ourtrans[i].stat
  trratingVals[i] = acn.ourtrans[i].rating_c #use rating_c for ctgy
 trqualswAMPL.setValues(trqualswVals)
 trcostswAMPL.setValues(trcostswVals)
 #trstatAMPL.setValues(trstatVals)
 trstatAMPL.setValues(basetrON)
 trratingAMPL.setValues(trratingVals)

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
 maxtime = ' maxtime_real=' + str(acn.maxtime) + ' mip_maxtime_real=' + str(acn.maxtime) + ' ms_maxtime_real=' + str(acn.maxtime) # just in case we use multi-start
 #maxtime = ' maxtime_real=600 '
 ctg_options = 'outlev=0 outmode=0 debug=0 linsolver=6 feastol=1e-5 ftol=1e-3 scale=0 honorbnds=1 cg_maxit=1 bar_murule=1 bar_refinement=1 bar_switchrule=0 bar_feasible=1 restarts=3 maxit=3000'  # use bar_initpt=2 bar_directinterval=0 ?
 if acn.useparallel:
     ctg_options += ' par_numthreads=1'
 just_feasible = 1 #ignoring optimality works quite well
 if just_feasible:
     ctg_options += ' opttol=1.0e30 opttol_abs=1e30'
 else:
     ctg_options += ' opttol=1.0e-3'
 if (acn.numbuses < 10000):
     ctg_options += ' bar_maxcrossit=1'
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
 options = ctg_options + mip_options + maxtime + relax + intvar_strategy
 #print "options=",options
 ampl.setOption('knitro_options', options)

 # Solve
 usescript = 1
 try:
  time1 = time.time()
  if usescript:
      #ampl.eval('include solve_script2;')
      solveout = ampl.getOutput('include solve_script2;')
  else:
      #ampl.solve()
      #ampl.eval('solve > junk.txt;')  # use this to re-direct ampl output to junk file
      #ampl.eval('solve > /dev/null;') # use this to re-direct ampl output to null device
      solveout = ampl.getOutput('solve;')
      #log.joint("Finish ampl.solve\n output=%s" %(solveout))
  time2 = time.time()
 except:
  log.joint("Error calling ampl.solve!")
  var = traceback.format_exc()
  #traceback.print_exc()
  #acn_evaluation2.print_alert(var, raise_exception=False)
  pass

 # Re-solve with integer vars fixed
 resolve_fixed2(acn, ampl, usescript, basegenON)

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
     amplstate += 'display -sum{f in 1..TRcnt_original: f != trout_index} (Deltact*trexceed[f]);'
     amplstate += 'display gencommit;'
     amplstate += 'display PBalance;'
     amplstate += 'display QBalance;'
     amplstate += 'display LineCurrentNTR;'
     amplstate += 'display LineCurrentTR;'
     amplstate += 'display Deltact;'
     amplstate += 'display trexceed, relaxStr;'
     amplstate += 'display Pcostcblock;'
     amplstate += 'display Qcostcblock;'
     amplstate += 'display Scostcblock;'
     amplstate += 'display Strblock;'
     amplstate += 'display minxstf, xstf, maxxstf;'
     filename = 'amplsol_'+ str(conlabel) +'.out'
     modelout = ampl.getOutput(amplstate)
     outfile = open(filename,"w")
     outfile.write("objective and constraint values: \n" + str(modelout) + "\n")
     outfile.close()

 # Write the contingency solution to file
 ctgytmpfile = "solution2_" + str(conlabel) + ".txt"
 write_ctg_solution(log, acn, ampl, genoutIndex, genoutID, lineout_i, lineout_j, lineout_ckt, ctgytmpfile)
 #breakexit('done write_ctg')

 # Evaluate contingency solution
 try:
  ctgsol_obj, ctgsol_infeas, exist, ctgsol_summary = acn_evaluation2.ctgy_from_file(acn, conlabel, ctgytmpfile)
  #breakexit("evalctg")
  if (ctgsol_infeas):
      log.joint(str(conlabel) + ':optimized solution infeasible!\n')
  ctgyfile = "solution_" + str(conlabel) + ".txt"
  #print('ctgyfile=',ctgyfile)
  pobj, pinfeas, pexist, ctgsol_psummary = acn_evaluation2.ctgy_from_file(acn, conlabel, ctgyfile)
  if (pinfeas):
      log.joint(str(conlabel) + ':prior solution infeasible!\n')
  log.joint('infeasible solution: obj=' + str(pobj)+'\n')
  log.joint('infeasible solution: infeas=' + str(pinfeas)+'\n')
  log.joint('contingency solution exists=' +  str(exist)+'\n')
  if (exist == True):
    log.joint('ctgsol_obj=' + str(ctgsol_obj)+'\n')
    log.joint('ctgsol_infeas=' + str(ctgsol_infeas)+'\n')
    if ctgsol_infeas == 0 and (pexist == False or ctgsol_obj > pobj):
        # We've computed a contingency case solution that is better than the
        # 'infeasible_solution', so overwrite with this one
        #os.popen('\cp solution2.txt solution_???.txt')
        shutil.copyfile(ctgytmpfile, ctgyfile)
        #breakexit("overwrite")
    else:
        log.joint('using infeasible solution!\n')
  if (ctgsol_infeas or ctgsol_obj<0):
      log.joint(str(conlabel) + ': bad solution!\n')
  #    breakexit("1: bad ctg solution")
  #    sys.exit()
 except:
  log.joint("Error evaluating the contingency solution!\n")
  #breakexit("acnsolverAMPLcon: Error evaluating the contingency solution!")
  '''
  var = traceback.format_exc()
  traceback.print_exc()
  acn_evaluation2.print_alert(var, raise_exception=False)
  pass
  '''


##### Called after initial solve to fix all integer variables and re-solve
#Maybe put this in a different file?

def resolve_fixed2(acn, ampl, usescript, basegenON):

  log = acn.log
  log.joint('running resolve_fixed2\n')

  # if usescript = 1, then we already fixed everything inside the AMPL script after
  # the initial solve, so we can skip the fixing code below

  if usescript == 0:

      log.joint(' usescript = 0\n')

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
    resolve_options = 'outlev=0 outmode=0 outname="knitro-fixed.log" debug=0 feastol=1e-5 feastol_abs=9e-5 ftol=1e-6 scale=0 honorbnds=0 cg_maxit=50 bar_murule=0 bar_feasible=0 bar_refinement=1 bar_initpi_mpec=0.0 maxit=3000 alg=1 strat_warm_start=1 bar_initpt=2 bar_initmu=1e-6 bar_slackboundpush=1e-6 infeastol=1e-5 restarts=1'
    mip_options = ' mip_heuristic=0 mip_terminate=1 mip_outinterval=1 mip_outlevel=1 mip_debug=0 mip_outsub=2 mip_nodealg=1 mip_intvar_strategy=0 mip_maxtime_real=60'
    resolve_options += mip_options
    just_feasible = 0
    if just_feasible:
      resolve_options += ' opttol=1.0e30 opttol_abs=1e30'
    else:
      resolve_options += ' opttol=1.0e-3'
    if (acn.numbuses < 10000):
      resolve_options += ' bar_maxcrossit=1'
    ampl.setOption('knitro_options', resolve_options)
    #ampl.solve()
    #ampl.eval('solve > junk.txt')  # use this to re-direct ampl output to junk file
    #ampl.eval('solve > /dev/null') # use this to re-direct ampl output to null device
    solveout = ampl.getOutput('include resolve_script2;')
    #solveout = ampl.getOutput('solve;')
    #log.joint("Finish ampl.solve\n output=%s" %(solveout))
    #breakexit("resolve")
  except:
    log.joint("Error in fixed re-solve!\n")
    pass
  #breakexit("bbb")

#############
