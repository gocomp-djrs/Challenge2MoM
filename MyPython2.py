'''
MyPython2.py

Entry to Python infeasibility solver, code1.
Syntax:

$ python MyPython2.py CON SUP RAW RESERVED TIMELIMIT DIVISION NETWORK
'''

# python3 MyPython2.py case.con case.json case.raw 0 600 1 C2S1N00014

# built in imports
import sys, os
import traceback

# modules for this code
# better way to make this visible?
sys.path.append(os.path.normpath('.'))
sys.path.append(os.path.normpath('./py'))
#from infeasibility_solution import Solver
from infeasibility_solution import *
from acn_solution import *
import evaluation2
import acn_evaluation2

from acnconvert import *
from myutils import *
from log import *
from versioner import stateversion
from acnsolverAMPL2 import *
from acn_solution import *
import multiprocessing

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

    print('\nPython infeasibility solver, code 2')
    print('syntax:')
    print('$ python MyPython2.py CON SUP RAW RESERVED TIMELIMIT DIVISION NETWORK')
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
    print("procid=",procid)

    #logfile = 'acn2.log'
    #logfile = os.devnull #suppress generating this log file
    logfile = "acn2_"+str(procid)+".log"
    #if len(sys.argv) == 3:
    #    logfile = sys.argv[2]
    log = danoLogger(logfile)
    #log.both_off()   #turn off logging

    stateversion(log)


    #alldata = read_config(log, sys.argv[1])
    alldata = startalldata(log, con_name, sup_name, raw_name)
    alldata['log'] = log
    alldata['mykeybus'] = -1
    alldata['mykeytrdual'] = -1
    alldata['mykeynontrdual'] = -1
    alldata['DIVISION'] = int(division)  #for consistency

    s = acnSolver()
    s.data.read(raw_name, sup_name, con_name)

    # Writing infeasibility solution (from infeasibility_solution.py) for contingencies
    # We'll fall back on these if necessary
    # use this to suppress stdout
    #if procid == 0:
    #    sys.stdout =  open(os.devnull, 'w')
    #    write_ctgs = True
    #    if write_ctgs:
    #        #s.write_sol('./')
    #        s.write_sol_ctgs_only('.','solution_BASECASE')
    #    sys.stdout = sys.__stdout__
    #    #sys.exit()

    alldata['arpae'] = s

    try:
        convertcode = acnconvert(alldata)
        acn = alldata['acn']
    except:
        log.joint("Error: MyPython2 calling acnconvert!\n")
        var = traceback.format_exc()
        traceback.print_exc()
        evaluation2.print_alert(var, raise_exception=False)
        pass

    # Call ACN Solver for contingencies
    acn.timebeg = timebeg
    acn.maxtime = 0.95*float(time_limit)
    acn.numnodes = 1 # specify number of nodes to use
    acn.division = division
    alldata['PROCEDURE'] = 'SOLVEAMPL1'
    alldata['FLAVOR'] = 'H2'
    try:
        if alldata['PROCEDURE'] == 'SOLVEAMPL1':
            acnsolverAMPL2(alldata)
            #acnsolverAMPL2par(alldata)
    except:
        log.joint("Error: MyPython2 calling acnsolverAMPL2!\n")
        var = traceback.format_exc()
        traceback.print_exc()
        evaluation2.print_alert(var, raise_exception=False)
        pass

    # Temporary hack to make sure procid=0 finishes last.  Have each process
    # write a temp file when finished and procid=0 will not finish until it
    # sees that all other temp files for procid>0 have been written.
    if procid > 0:
      pname = "proc"+str(procid)+".log"
      proc = open(pname,"w")
      elapsedtime = time.time()-acn.timebeg
      proc.write("Done: elapsed time="+str(elapsedtime)+"\n")
      proc.close()

    #procid=0 should exit last;
    proccnt = 1
    if procid == 0:
      #elapsedtime = time.time()-acn.timebeg
      #remaintime = acn.maxtime - elapsedtime
      while proccnt < acn.numnodes:
          pname = "proc"+str(proccnt)+".log"
          try:
              tmp = open(pname,"r")
              tmp.close()
              proccnt = proccnt + 1
          except:
              proccnt = 1 # start over

    # Evaluate score
    printsummary = 1
    if procid == 0 and printsummary == 1:
        line_switching_allowed = True if division == '3' or division == '4' else None
        xfmr_switching_allowed = True if division == '3' or division == '4' else None
        check_contingencies = True # set to False to check just the base case
        try:
            # use this to suppress stdout
            sys.stdout =  open(os.devnull, 'w')
            obj, infeas, exist, summary = evaluation2.run(raw_name, con_name, sup_name, './', "", "summary.csv", "detail.csv", line_switching_allowed, xfmr_switching_allowed, check_contingencies)
            #obj, infeas, exist, summary = acn_evaluation2.run(acn, './', "", "summary.csv", "detail.csv", line_switching_allowed, xfmr_switching_allowed, check_contingencies)
            sys.stdout = sys.__stdout__
            timeend = time.time()
            total_time = timeend-timebeg
            print("================FINAL SCORE==============")
            print("obj = ", obj)
            print("infeas = ", infeas)
            print("code2: total time = ", total_time)
            log.joint("================FINAL SCORE==============\n")
            log.joint("obj = " + str(obj) + "\n")
            log.joint("infeas = " + str(infeas)  + "\n")
            log.joint("code2: total time = " + str(total_time)  + "\n")
            if exist == False:
                print("WARNING: Final Solution does not exist!?")
                log.joint("WARNING: Final Solution does not exist!?\n")
            summary = open("summary.txt","a+")
            summary.write(str(obj) + "," + str(infeas) + "," + str(total_time) + "\n")
            summary.close()
        except:
            print("Error: MyPython2 evaluating the solution!")
            log.joint("Error: MyPython2 evaluating the solution!\n")
            var = traceback.format_exc()
            traceback.print_exc()
            evaluation2.print_alert(var, raise_exception=False)
            pass

    log.joint('done\n')
    log.closelog()

    print("Done: procid=",procid," time=",time.time()-timebeg)
    #sys.exit(0)
    os._exit(0)


if __name__ == '__main__':
    run()
