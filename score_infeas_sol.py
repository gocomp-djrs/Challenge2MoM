'''
score_infeas_sol.py

Entry to Python infeasibility solver.
Syntax:

$ python score_infeas_sol.py CON SUP RAW
'''

# python3 score_infeas_sol.py case.con case.json case.raw

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
#import acn_evaluation2

from acnconvert import *
from myutils import *
from log import *
from versioner import stateversion
from acnsolverAMPL import *


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
    reserved = 0
    time_limit = 600
    division = 1
    network = 'C2S6N00594'
    write_ctgs_in_code1 = True
    save_base_case_sol = True
    saved_base_case_sol_file_name = 'saved_solution_BASECASE.bin'

    print('\nCode to score infeasible solution')
    print('syntax:')
    print('$ python score_infeas_sol.py CON SUP RAW')
    print('args:')
    print(args)

    logfile = 'acn.log'
    if len(sys.argv) == 3:
        logfile = sys.argv[2]
    log = danoLogger(logfile)

    stateversion(log)

    #alldata = read_config(log, sys.argv[1])
    alldata = startalldata(log, con_name, sup_name, raw_name)
    alldata['log'] = log
    alldata['mykeybus'] = -1
    alldata['mykeytrdual'] = -1

    s = Solver()
    s.data.read(raw_name, sup_name, con_name)

    # Writing infeasibility solution (from infeasibility_solution.py)
    if write_ctgs_in_code1:
        s.write_sol('./')
    else:
        if save_base_case_sol:
            s.write_sol1('./', saved_sol_file_name=saved_base_case_sol_file_name)
        else:
            s.write_sol1('.')

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
            print("obj = ",obj)
            print("infeas = ",infeas)
            print("code2: total time = ",total_time)
            summary = open("summary.txt","a+")
            summary.write(str(obj) + "," + str(infeas) + "," + str(total_time) + "\n")
            summary.close()
        except:
            log.joint("Error: evaluating the solution!\n")
            var = traceback.format_exc()
            traceback.print_exc()
            evaluation2.print_alert(var, raise_exception=False)
            pass

    log.joint('done\n')
    log.closelog()


if __name__ == '__main__':
    run()
