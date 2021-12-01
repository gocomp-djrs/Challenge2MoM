'''
MyPython1.py

Entry to Python infeasibility solver, code1.
Syntax:

$ python MyPython1.py CON SUP RAW RESERVED TIMELIMIT DIVISION NETWORK
'''

# python3 MyPython1.py case.con case.json case.raw 0 600 1 C2S1N00014

# built in imports
import sys, os
import traceback

# modules for this code
# better way to make this visible?
sys.path.append(os.path.normpath('.'))
sys.path.append(os.path.normpath('./py'))
#from infeasibility_solution import Solver
from infeasibility_solution import *
from acnprior import *
import evaluation2
import acn_evaluation2

from acnconvert import *
from myutils import *
from log import *
from versioner import stateversion
from acnsolverAMPL import *
from acn_solution import *
import multiprocessing

def trace(frame, event, arg):
    print("%s, %s:%d" % (event, frame.f_code.co_filename, frame.f_lineno))
    return trace

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

def startalldata(log, con_name, sup_name, raw_name):

    log.joint("Run startalldata\n")

    alldata = {}

    rawfile = confile = jsonfile = 'NONE'

    rawfile = raw_name
    confile = con_name
    jsonfile = sup_name


    for x in [('RAW',rawfile), ('CON',confile), ('JSON',jsonfile)]:
      if x[1] == 'NONE':
        log.stateandquit(' no ' + x[0] + ' input'+'\n')
      if x[1] != 'NONE': log.joint(x[0] + " : " + str(x[1])+'\n')
      alldata[x[0]] = x[1]

    return alldata

def run():

    timebeg = time.time()

    # To enable debug tracing
    #sys.settrace(trace)

    args = sys.argv

    con_name = args[1]
    sup_name = args[2]
    raw_name = args[3]
    reserved = args[4]
    time_limit = args[5]
    division = args[6]
    network = args[7]
    #sol1_name = 'solution1.txt'
    #sol2_name = 'solution2.txt'
    #write_ctgs_in_code1 = False

    print('\nPython code 1')
    print('syntax:')
    print('$ python MyPython1.py CON SUP RAW RESERVED TIMELIMIT DIVISION NETWORK')
    print('args:')
    print(args)


    # Get processor ID
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


    logfile = "acn_"+str(procid)+".log"
    if len(sys.argv) == 3:
        logfile = sys.argv[2]
    log = danoLogger(logfile)
    #log.both_off()   #turn off logging

    stateversion(log)

    #try:
    #  cmd = 'gunzip -N knitroampl.gz'
    #  os.system(cmd)
    #  log.joint("Unzipped knitroampl.gz\n")
    #except:
    #  log.joint("Failed to unzip knitroampl.gz\n")

    #alldata = read_config(log, sys.argv[1])
    alldata = startalldata(log, con_name, sup_name, raw_name)
    alldata['log'] = log
    alldata['mykeybus'] = -1
    alldata['mykeytrdual'] = -1
    alldata['mykeynontrdual'] = -1
    alldata['DIVISION'] = int(division)  #for consistency

    s = acnSolver()
    s.data.read(raw_name, sup_name, con_name)

    if procid == 0:
      # Writing infeasibility solution (from infeasibility_solution.py)
      # We'll fall back on these if necessary
      log.joint("Writing infeasible base solution\n")
      write_ctgs_in_code1 = False
      save_base_case_sol = True
      saved_base_case_sol_file_name = 'solution_BASECASE.bin'
      try:
         if write_ctgs_in_code1:
           s.write_sol('./')
         else:
           if save_base_case_sol:
               s.write_sol1('./', saved_sol_file_name=saved_base_case_sol_file_name)
           else:
               s.write_sol1('.')
         log.joint("Finished writing infeasible base solution\n")
      except:
         log.joint("Error writing infeasible base solution (MyPython1)!\n")
         pass


    alldata['arpae'] = s

    try:
        convertcode = acnconvert(alldata)
        acn = alldata['acn']
    except:
        log.joint("Error: MyPython1 calling acnconvert!\n")
        var = traceback.format_exc()
        traceback.print_exc()
        evaluation2.print_alert(var, raise_exception=False)
        pass

    # Call ACN Solver for base case solution
    acn.timebeg = timebeg
    acn.maxtime = 0.90*float(time_limit)
    acn.numnodes = 1 # specify number of nodes to use
    acn.division = division
    acn.improvedbase = 0
    acn.procid = procid
    alldata['PROCEDURE'] = 'SOLVEAMPL1'
    alldata['PROCEDURE'] = 'HEURBASE1'
    #alldata['PROCEDURE'] = 'COMPUTE_PRIOR_BASE'
    try:
        solveit = 1
        #alldata['PROCEDURE'] = 'COMPUTE_PRIOR_BASE'
        #priorcode = acncomputepriorbase(alldata)
        #evaluate infeasible/prior solution (is this needed?)
        # This evaluation changes acn.pimb*/acn.qimb, which affects computation
        # of acn.totviolvector, which is used in acnheurbase1 fixing heuristic threshold!
        #pobj, pinfeas, pexist, base.psummary = acn_evaluation2.base_from_file(acn, "priorbase.txt")
        #alldata['PROCEDURE'] = 'HEURBASE1'
        #if acn.numbuses > 20000 and 0:
        #    if (pexist == True) and (pinfeas==0) and (pobj > 0):
        #        solveit = 0
        #        print("Using infeasible/prior solution for base: numbuses=",acn.numbuses," obj=",pobj)
        if solveit:
          if alldata['PROCEDURE'] == 'SOLVEAMPL1':
            acnsolverAMPL(alldata)
          elif alldata['PROCEDURE'] == 'HEURBASE1':
            priorcode = acncomputepriorbase(alldata)
            #priorcode = acnsolverHeurBase1(alldata)
          elif alldata['PROCEDURE'] == 'COMPUTE_PRIOR_BASE':
            acncomputepriorbase(alldata)
            raw = alldata['RAW']
            con = alldata['CON']
            sup = alldata['JSON']
            pobj, pinfeas, pexist, base.psummary = acn_evaluation2.base_from_file(acn, "priorbase.txt")
            if (pexist == True):
                print('prior base solution: obj=' + str(pobj))
                print('prior base solution: infeas=' + str(pinfeas))
    except:
        log.joint("Error: MyPython1 calling acnsolverAMPL!\n")
        var = traceback.format_exc()
        traceback.print_exc()
        evaluation2.print_alert(var, raise_exception=False)
        pass

    # Mark that procid=0 is done
    if procid == 0:
          pname = "code1_proc0.log"
          proc = open(pname,"w")
          elapsedtime = time.time()-acn.timebeg
          proc.write("Done: elapsed time="+str(elapsedtime)+"\n")
          proc.close()

    # Check if procid==1 solution better than current
    if procid == 1:
          ready = 0
          while ready == 0:
             pname = "code1_proc0.log"
             try:
                 tmp = open(pname,"r")
                 tmp.close()
                 ready = 1
             except:
                 ready = 0 # keep waiting
          # Now compare procid=1 solution to procid=0 solution and
          # update base solution if procid=1 is better
          log.joint("\nEvaluating heur1_1.txt...\n")
          base.obj, base.infeas, exist, base.summary = acn_evaluation2.base_from_file(acn, "heur1_1.txt")
          log.joint("\nEvaluating solution_BASECASE.txt...\n")
          pobj, pinfeas, pexist, base.psummary = acn_evaluation2.base_from_file(acn, "solution_BASECASE.txt")
          log.joint("\n------------------------------------------\n")
          if (pexist == True):
              log.joint('current solution: obj=' + str(pobj)+'\n')
              log.joint('current solution: infeas=' + str(pinfeas)+'\n')
          log.joint('heur1 exists=' +  str(exist)+'\n')
          log.joint('solution_BASECASE exists=' +  str(pexist)+'\n')
          if (exist == True):
              log.joint('base.obj=' + str(base.obj)+'\n')
              log.joint('base.infeas=' + str(base.infeas)+'\n')
              if base.infeas == 0 and (pexist == False or base.obj > pobj):
                  log.joint('heur1_1 solution better than current solution.')
                  # We've computed a solution that is better than the
                  # current solution, so overwrite with this one
                  os.popen('\cp heur1_1.txt solution_BASECASE.txt')
                  os.popen('\cp heur1_1.bin solution_BASECASE.bin')
                  time.sleep(1.0) #allow a little to make sure copy above is finishes
                  # Now write out individual base solution components
                  # needed for contingency modeling to files so that they can be read in by
                  # MyPython2.py->acnsolverAMPL2.py
                  write_base_files(acn, acn.basebusVmags1,acn.basebusVangles1,acn.baseploads1,acn.basegenP1,acn.basegenON1,acn.basegensu1,acn.basegensd1,acn.basentrON1,acn.basetrON1,acn.basesw1)
                  acn.improvedbase = 1 #???
                  acn.basebusVmags = acn.basebusVmags1
                  acn.basebusVangles = acn.basebusVangles1
                  acn.clearedloads = acn.clearedloads1
                  acn.basegenP = acn.basegenP1
                  acn.basegenQ = acn.basegenQ1
                  acn.basegenON = acn.basegenON1
                  acn.basentrON = acn.basentrON1
                  acn.basetrON = acn.basetrON1
                  acn.basetaps = acn.basetaps1
                  acn.basesw = acn.basesw1
                  pname = "code1_proc1_improved.log"
                  proc = open(pname,"w")
                  proc.write("Improved: obj="+str(base.obj)+"\n")
                  proc.close()
          # Mark procid=1 is done
          pname = "code1_proc1.log"
          proc = open(pname,"w")
          elapsedtime = time.time()-acn.timebeg
          proc.write("Done: elapsed time="+str(elapsedtime)+"\n")
          proc.close()

    # Check if procid==2 solution better than current
    if procid == 2:
          ready = 0
          while ready == 0:
             pname = "code1_proc1.log"
             try:
                 tmp = open(pname,"r")
                 tmp.close()
                 ready = 1
             except:
                 ready = 0 # keep waiting
          # Now compare procid=2 solution to current solution and
          # update base solution if procid=2 is better
          log.joint("\nEvaluating heur1_2.txt...\n")
          base.obj, base.infeas, exist, base.summary = acn_evaluation2.base_from_file(acn, "heur1_2.txt")
          log.joint("\nEvaluating solution_BASECASE.txt...\n")
          pobj, pinfeas, pexist, base.psummary = acn_evaluation2.base_from_file(acn, "solution_BASECASE.txt")
          log.joint("\n------------------------------------------\n")
          if (pexist == True):
              log.joint('current solution: obj=' + str(pobj)+'\n')
              log.joint('current solution: infeas=' + str(pinfeas)+'\n')
          log.joint('heur1 exists=' +  str(exist)+'\n')
          log.joint('solution_BASECASE exists=' +  str(pexist)+'\n')
          if (exist == True):
              log.joint('base.obj=' + str(base.obj)+'\n')
              log.joint('base.infeas=' + str(base.infeas)+'\n')
              if base.infeas == 0 and (pexist == False or base.obj > pobj):
                  log.joint('heur1_2 solution better than current solution.')
                  # We've computed a solution that is better than the
                  # current solution, so overwrite with this one
                  os.popen('\cp heur1_2.txt solution_BASECASE.txt')
                  os.popen('\cp heur1_2.bin solution_BASECASE.bin')
                  time.sleep(1.0) #allow a little to make sure copy above is finishes
                  # Now write out individual base solution components
                  # needed for contingency modeling to files so that they can be read in by
                  # MyPython2.py->acnsolverAMPL2.py
                  write_base_files(acn, acn.basebusVmags2,acn.basebusVangles2,acn.baseploads2,acn.basegenP2,acn.basegenON2,acn.basegensu2,acn.basegensd2,acn.basentrON2,acn.basetrON2,acn.basesw2)
                  acn.improvedbase = 1 #???
                  acn.basebusVmags = acn.basebusVmags2
                  acn.basebusVangles = acn.basebusVangles2
                  acn.clearedloads = acn.clearedloads2
                  acn.basegenP = acn.basegenP2
                  acn.basegenQ = acn.basegenQ2
                  acn.basegenON = acn.basegenON2
                  acn.basentrON = acn.basentrON2
                  acn.basetrON = acn.basetrON2
                  acn.basetaps = acn.basetaps2
                  acn.basesw = acn.basesw2
                  pname = "code1_proc2_improved.log"
                  proc = open(pname,"w")
                  proc.write("Improved: obj="+str(base.obj)+"\n")
                  proc.close()
          # Mark procid=2 is done
          pname = "code1_proc2.log"
          proc = open(pname,"w")
          elapsedtime = time.time()-acn.timebeg
          proc.write("Done: elapsed time="+str(elapsedtime)+"\n")
          proc.close()

    # Check if procid==3 solution better than current
    if procid == 3:
          ready = 0
          while ready == 0:
             pname = "code1_proc2.log"
             try:
                 tmp = open(pname,"r")
                 tmp.close()
                 ready = 1
             except:
                 ready = 0 # keep waiting
          # Now compare procid=3 solution to current solution and
          # update base solution if procid=3 is better
          log.joint("\nEvaluating heur1_3.txt...\n")
          base.obj, base.infeas, exist, base.summary = acn_evaluation2.base_from_file(acn, "heur1_3.txt")
          log.joint("\nEvaluating solution_BASECASE.txt...\n")
          pobj, pinfeas, pexist, base.psummary = acn_evaluation2.base_from_file(acn, "solution_BASECASE.txt")
          log.joint("\n------------------------------------------\n")
          if (pexist == True):
              log.joint('current solution: obj=' + str(pobj)+'\n')
              log.joint('current solution: infeas=' + str(pinfeas)+'\n')
          log.joint('heur1 exists=' +  str(exist)+'\n')
          log.joint('solution_BASECASE exists=' +  str(pexist)+'\n')
          if (exist == True):
              log.joint('base.obj=' + str(base.obj)+'\n')
              log.joint('base.infeas=' + str(base.infeas)+'\n')
              if base.infeas == 0 and (pexist == False or base.obj > pobj):
                  log.joint('heur1_3 solution better than current solution.')
                  # We've computed a solution that is better than the
                  # current solution, so overwrite with this one
                  os.popen('\cp heur1_3.txt solution_BASECASE.txt')
                  os.popen('\cp heur1_3.bin solution_BASECASE.bin')
                  time.sleep(1.0) #allow a little to make sure copy above is finishes
                  # Now write out individual base solution components
                  # needed for contingency modeling to files so that they can be read in by
                  # MyPython2.py->acnsolverAMPL2.py
                  write_base_files(acn, acn.basebusVmags3,acn.basebusVangles3,acn.baseploads3,acn.basegenP3,acn.basegenON3,acn.basegensu3,acn.basegensd3,acn.basentrON3,acn.basetrON3,acn.basesw3)
                  acn.improvedbase = 1 #???
                  acn.basebusVmags = acn.basebusVmags3
                  acn.basebusVangles = acn.basebusVangles3
                  acn.clearedloads = acn.clearedloads3
                  acn.basegenP = acn.basegenP3
                  acn.basegenQ = acn.basegenQ3
                  acn.basegenON = acn.basegenON3
                  acn.basentrON = acn.basentrON3
                  acn.basetrON = acn.basetrON3
                  acn.basetaps = acn.basetaps3
                  acn.basesw = acn.basesw3
                  pname = "code1_proc3_improved.log"
                  proc = open(pname,"w")
                  proc.write("Improved: obj="+str(base.obj)+"\n")
                  proc.close()
          # Mark procid=3 is done
          pname = "code1_proc3.log"
          proc = open(pname,"w")
          elapsedtime = time.time()-acn.timebeg
          proc.write("Done: elapsed time="+str(elapsedtime)+"\n")
          proc.close()

    # Check if procid==4 solution better than current
    if procid == 4:
          ready = 0
          while ready == 0:
             pname = "code1_proc3.log"
             try:
                 tmp = open(pname,"r")
                 tmp.close()
                 ready = 1
             except:
                 ready = 0 # keep waiting
          # Now compare procid=4 solution to current solution and
          # update base solution if procid=4 is better
          log.joint("\nEvaluating heur1_4.txt...\n")
          base.obj, base.infeas, exist, base.summary = acn_evaluation2.base_from_file(acn, "heur1_4.txt")
          log.joint("\nEvaluating solution_BASECASE.txt...\n")
          pobj, pinfeas, pexist, base.psummary = acn_evaluation2.base_from_file(acn, "solution_BASECASE.txt")
          log.joint("\n------------------------------------------\n")
          if (pexist == True):
              log.joint('current solution: obj=' + str(pobj)+'\n')
              log.joint('current solution: infeas=' + str(pinfeas)+'\n')
          log.joint('heur1 exists=' +  str(exist)+'\n')
          log.joint('solution_BASECASE exists=' +  str(pexist)+'\n')
          if (exist == True):
              log.joint('base.obj=' + str(base.obj)+'\n')
              log.joint('base.infeas=' + str(base.infeas)+'\n')
              if base.infeas == 0 and (pexist == False or base.obj > pobj):
                  log.joint('heur1_4 solution better than current solution.')
                  # We've computed a solution that is better than the
                  # current solution, so overwrite with this one
                  os.popen('\cp heur1_4.txt solution_BASECASE.txt')
                  os.popen('\cp heur1_4.bin solution_BASECASE.bin')
                  time.sleep(1.0) #allow a little to make sure copy above is finishes
                  # Now write out individual base solution components
                  # needed for contingency modeling to files so that they can be read in by
                  # MyPython2.py->acnsolverAMPL2.py
                  write_base_files(acn, acn.basebusVmags4,acn.basebusVangles4,acn.baseploads4,acn.basegenP4,acn.basegenON4,acn.basegensu4,acn.basegensd4,acn.basentrON4,acn.basetrON4,acn.basesw4)
                  acn.improvedbase = 1 #???
                  acn.basebusVmags = acn.basebusVmags4
                  acn.basebusVangles = acn.basebusVangles4
                  acn.clearedloads = acn.clearedloads4
                  acn.basegenP = acn.basegenP4
                  acn.basegenQ = acn.basegenQ4
                  acn.basegenON = acn.basegenON4
                  acn.basentrON = acn.basentrON4
                  acn.basetrON = acn.basetrON4
                  acn.basetaps = acn.basetaps4
                  acn.basesw = acn.basesw4
                  pname = "code1_proc4_improved.log"
                  proc = open(pname,"w")
                  proc.write("Improved: obj="+str(base.obj)+"\n")
                  proc.close()
          # Mark procid=4 is done
          pname = "code1_proc4.log"
          proc = open(pname,"w")
          elapsedtime = time.time()-acn.timebeg
          proc.write("Done: elapsed time="+str(elapsedtime)+"\n")
          proc.close()

    # Check if procid==5 solution better than current
    if procid == 5:
          ready = 0
          while ready == 0:
             pname = "code1_proc4.log"
             try:
                 tmp = open(pname,"r")
                 tmp.close()
                 ready = 1
             except:
                 ready = 0 # keep waiting
          # Now compare procid=5 solution to current solution and
          # update base solution if procid=5 is better
          log.joint("\nEvaluating heur1_5.txt...\n")
          base.obj, base.infeas, exist, base.summary = acn_evaluation2.base_from_file(acn, "heur1_5.txt")
          log.joint("\nEvaluating solution_BASECASE.txt...\n")
          pobj, pinfeas, pexist, base.psummary = acn_evaluation2.base_from_file(acn, "solution_BASECASE.txt")
          log.joint("\n------------------------------------------\n")
          if (pexist == True):
              log.joint('current solution: obj=' + str(pobj)+'\n')
              log.joint('current solution: infeas=' + str(pinfeas)+'\n')
          log.joint('heur1 exists=' +  str(exist)+'\n')
          log.joint('solution_BASECASE exists=' +  str(pexist)+'\n')
          if (exist == True):
              log.joint('base.obj=' + str(base.obj)+'\n')
              log.joint('base.infeas=' + str(base.infeas)+'\n')
              if base.infeas == 0 and (pexist == False or base.obj > pobj):
                  log.joint('heur1_5 solution better than current solution.')
                  # We've computed a solution that is better than the
                  # current solution, so overwrite with this one
                  os.popen('\cp heur1_5.txt solution_BASECASE.txt')
                  os.popen('\cp heur1_5.bin solution_BASECASE.bin')
                  # Now write out individual base solution components
                  # needed for contingency modeling to files so that they can be read in by
                  # MyPython2.py->acnsolverAMPL2.py
                  write_base_files(acn, acn.basebusVmags5,acn.basebusVangles5,acn.baseploads5,acn.basegenP5,acn.basegenON5,acn.basegensu5,acn.basegensd5,acn.basentrON5,acn.basetrON5,acn.basesw5)
                  acn.improvedbase = 1 #???
                  acn.basebusVmags = acn.basebusVmags5
                  acn.basebusVangles = acn.basebusVangles5
                  acn.clearedloads = acn.clearedloads5
                  acn.basegenP = acn.basegenP5
                  acn.basegenQ = acn.basegenQ5
                  acn.basegenON = acn.basegenON5
                  acn.basentrON = acn.basentrON5
                  acn.basetrON = acn.basetrON5
                  acn.basetaps = acn.basetaps5
                  acn.basesw = acn.basesw5
                  pname = "code1_proc5_improved.log"
                  proc = open(pname,"w")
                  proc.write("Improved: obj="+str(base.obj)+"\n")
                  proc.close()
          # Mark procid=5 is done
          pname = "code1_proc5.log"
          proc = open(pname,"w")
          elapsedtime = time.time()-acn.timebeg
          proc.write("Done: elapsed time="+str(elapsedtime)+"\n")
          proc.close()

    # Temporary hack to make sure procid=0 finishes last.  Have each process
    # write a temp file when finished and procid=0 will not finish until it
    # sees that all other temp files for procid>0 have been written.

    #if procid > 0:  # now done above
    #      pname = "code1_proc"+str(procid)+".log"
    #      proc = open(pname,"w")
    #      elapsedtime = time.time()-acn.timebeg
    #      proc.write("Done: elapsed time="+str(elapsedtime)+"\n")
    #      proc.close()

    #procid=0 should exit last;
    proccnt = 1
    if procid == 0:
          #elapsedtime = time.time()-acn.timebeg
          #remaintime = acn.maxtime - elapsedtime
          while proccnt < acn.numnodes:
              pname = "code1_proc"+str(proccnt)+".log"
              try:
                  tmp = open(pname,"r")
                  tmp.close()
                  proccnt = proccnt + 1
              except:
                  proccnt = 1 # start over

    print("Done: procid=",procid," time=",time.time()-timebeg)
    show_final_base = 0
    if procid == 0 and show_final_base:
        pobj, pinfeas, pexist, base.psummary = acn_evaluation2.base_from_file(acn, "solution_BASECASE.txt")
        if (pexist == True):
            print('\n\nFINAL BASE solution: obj=' + str(pobj))
            print('FINAL BASE solution: infeas=' + str(pinfeas))


    # Writing infeasibility solution (from infeasibility_solution.py) for contingencies
    # We'll fall back on these if necessary
    ctimebeg = time.time()
    if procid == 0:
        #sys.stdout =  open(os.devnull, 'w')
        # Check if base solution was improved by solve on another node
        if acn.improvedbase == 0 and acn.numnodes>1:
            for a in range(1,acn.numnodes):
                pname = "code1_proc"+str(a)+"_improved.log"
                try:
                  tmp = open(pname,"r")
                  tmp.close()
                  acn.improvedbase = 1
                except:
                  pass
        write_ctgs = True
        if write_ctgs:
            print("Writing contingencies...improved base=",acn.improvedbase)
            #s.write_sol('./')
            if solveit and acn.improvedbase:
                basebusVmags = file_to_dict('basebusVmags.txt', 0)
                basebusVangles = file_to_dict('basebusVangles.txt', 0)
                basebusVmags_best = numpy.zeros(acn.numbuses)
                basebusVangles_best = numpy.zeros(acn.numbuses)
                for i in range(acn.numbuses):
                    basebusVmags_best[i] = basebusVmags[i+1]
                    basebusVangles_best[i] = basebusVangles[i+1]
                s.write_sol_ctgs_only_new('.','solution_BASECASE',acn,basebusVmags_best,basebusVangles_best)
            else:
                s.write_sol_ctgs_only('.','solution_BASECASE')
        #sys.stdout = sys.__stdout__
        #sys.exit()
        pname = "write_ctgs.log"
        wrctg = open(pname,"w")
        elapsedtime = time.time()-acn.timebeg
        wrctg.write("Finished writing fallback contingencies: elapsed time="+str(elapsedtime)+"\n")
        wrctg.close()
        ctime = time.time()-ctimebeg
        print("code1: write contingencies time = ",ctime)

    total_time = time.time()-timebeg
    print("code1: total time = ",total_time)
    log.joint('done\n')
    log.closelog()


if __name__ == '__main__':
    run()
