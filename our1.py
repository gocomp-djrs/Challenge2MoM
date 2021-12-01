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
#import evaluation2
import acn_evaluation2

from acnconvert import *
from myutils import *
from log import *
from versioner import stateversion
from acnsolverAMPL import *

def trace(frame, event, arg):
    print("%s, %s:%d" % (event, frame.f_code.co_filename, frame.f_lineno))
    return trace

def readconfig1(filename):
    print("reading config file " + filename + "\n")
    try:
        f = open(filename, "r")
        lines = f.readlines()
        f.close()
    except:
        sys.exit("cannot open file " + filename + "\n")
    if len(lines) < 7:
        sys.exit("too short\n")

    linenum = 0
    out = {}
    while linenum < 7:
     thisline = lines[linenum].split()
     out[linenum+1] = thisline[0]
     linenum += 1

    out[5] = float(out[5])
    out[6] = int(out[6])
    
    return out[1],out[2],out[3],out[4],out[5],out[6],out[7]


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
    if len(args) > 2:
        con_name = args[1]
        sup_name = args[2]
        raw_name = args[3]
        reserved = args[4]
        time_limit = args[5]
        division = args[6]
        network = args[7]
    else:
        if len(args) == 2:
            filename = args[1]
        else:
            filename = "one.conf"
        con_name, sup_name, raw_name, reserved, time_limit, division, network = readconfig1(filename)


    
    #sol1_name = 'solution1.txt'
    #sol2_name = 'solution2.txt'
    #write_ctgs_in_code1 = False

    print('\nPython code 1')
    print('syntax:')
    print('$ python MyPython1.py CON SUP RAW RESERVED TIMELIMIT DIVISION NETWORK')
    print('args:')
    print(args)

    logfile = 'acn.log'
    if len(sys.argv) == 3:
        logfile = sys.argv[2]
    log = danoLogger(logfile)

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
    alldata['DIVISION'] = division  #for consistency

    s = Solver()
    s.data.read(raw_name, sup_name, con_name)

    # Writing infeasibility solution (from infeasibility_solution.py)
    # We'll fall back on these if necessary
    write_ctgs_in_code1 = False
    save_base_case_sol = True
    saved_base_case_sol_file_name = 'solution_BASECASE.bin'
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

    # Call ACN Solver for base case solution
    acn.timebeg = timebeg
    acn.maxtime = 0.95*float(time_limit)
    acn.division = division #not really needed any more
    alldata['PROCEDURE'] = 'SOLVEAMPL1'
    alldata['PROCEDURE'] = 'HEURBASE1'
    #alldata['PROCEDURE'] = 'COMPUTE_PRIOR_BASE'
    try:
        if alldata['PROCEDURE'] == 'SOLVEAMPL1':
            acnsolverAMPL(alldata)
        elif alldata['PROCEDURE'] == 'HEURBASE1':
            priorcode = acncomputepriorbase(alldata)
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
        '''
        var = traceback.format_exc()
        traceback.print_exc()
        evaluation2.print_alert(var, raise_exception=False)
        pass
        '''

    show_final_base = 1
    if show_final_base:
        pobj, pinfeas, pexist, base.psummary = acn_evaluation2.base_from_file(acn, "solution_BASECASE.txt")
        print("pexist=",pexist)
        if (pexist == True):
            print('\n\nFINAL BASE solution: obj=' + str(pobj))
            print('FINAL BASE solution: infeas=' + str(pinfeas))

    total_time = time.time()-timebeg
    print("code1: total time = ",total_time)
    log.joint('done\n')
    log.closelog()


if __name__ == '__main__':
    run()
