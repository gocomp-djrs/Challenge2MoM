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

    logfile = 'acn2.log'
    if len(sys.argv) == 3:
        logfile = sys.argv[2]
    log = danoLogger(logfile)

    stateversion(log)


    #alldata = read_config(log, sys.argv[1])
    alldata = startalldata(log, con_name, sup_name, raw_name)
    alldata['log'] = log
    alldata['mykeybus'] = -1
    alldata['mykeytrdual'] = -1
    alldata['mykeynontrdual'] = -1
    alldata['DIVISION'] = division  #for consistency


    s = acnSolver()
    s.data.read(raw_name, sup_name, con_name)

    # Writing infeasibility solution (from infeasibility_solution.py) for contingencies
    # We'll fall back on these if necessary
    write_ctgs = True
    if write_ctgs:
        #s.write_sol('./')
        s.write_sol_ctgs_only('.','solution_BASECASE')
    #sys.exit()

    alldata['arpae'] = s

    try:
        convertcode = acnconvert(alldata)
        acn = alldata['acn']
    except:
        log.joint("exception: MyPython2 calling acnconvert!\n")
        var = traceback.format_exc()
        traceback.print_exc()
        evaluation2.print_alert(var, raise_exception=False)
        pass

    # Call ACN Solver for contingencies
    acn.timebeg = timebeg
    acn.maxtime = 0.95*float(time_limit)
    acn.division = division
    alldata['PROCEDURE'] = 'SOLVEAMPL1'
    alldata['FLAVOR'] = 'NOTH2'
    try:
        if alldata['PROCEDURE'] == 'SOLVEAMPL1':
            acnsolverAMPL2(alldata)
    except:
        log.joint("exception: MyPython2 calling acnsolverAMPL2!\n")
        var = traceback.format_exc()
        traceback.print_exc()
        evaluation2.print_alert(var, raise_exception=False)
        pass

    # Evaluate score
    printsummary = 1
    if printsummary == 1:
        line_switching_allowed = True if division == '3' or division == '4' else None
        xfmr_switching_allowed = True if division == '3' or division == '4' else None
        check_contingencies = True # set to False to check just the base case
        try:
            obj, infeas, exist, summary = evaluation2.run(raw_name, con_name, sup_name, './', "", "summary.csv", "detail.csv", line_switching_allowed, xfmr_switching_allowed, check_contingencies)
            #obj, infeas, exist, summary = acn_evaluation2.run(acn, './', "", "summary.csv", "detail.csv", line_switching_allowed, xfmr_switching_allowed, check_contingencies)
            timeend = time.time()
            total_time = timeend-timebeg
            log.joint("================FINAL SCORE==============\n")
            log.joint("obj = " + str(obj) + "\n")
            log.joint("infeas = " + str(infeas)  + "\n")
            log.joint("code2: total time = " + str(total_time)  + "\n")
            if exist == False:
                log.joint("WARNING: Final Solution does not exist!?\n")
            summary = open("summary.txt","a+")
            summary.write(str(obj) + "," + str(infeas) + "," + str(total_time) + "\n")
            summary.close()
        except:
            log.joint("exception: MyPython2 evaluating the solution!\n")
            var = traceback.format_exc()
            traceback.print_exc()
            evaluation2.print_alert(var, raise_exception=False)
            pass

    log.joint('done\n')
    log.closelog()


if __name__ == '__main__':
    run()
