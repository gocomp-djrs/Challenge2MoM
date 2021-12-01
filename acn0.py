from myutils import *
from log import *
from versioner import stateversion
import sys, os
import os.path
from os import path
from acnsolverAMPL import *
from acnprior import *
from infeasibility_solution import *
from acnconvert import *
import acn_evaluation2
from acnxtra import  acnprintgraph
def read_config(log, filename):

    log.joint("reading config file " + filename + "\n")

    try:
        f = open(filename, "r")
        lines = f.readlines()
        f.close()
    except:
        log.stateandquit("cannot open file " + filename + "\n")

    alldata = {}

    rawfile = confile = jsonfile = procedure = 'NONE'
    mykeybus = mykeynontrdual = mykeytrdual = -1
    division = 1    
    evalfile = 'NONE'

    linenum = 0
    while linenum < len(lines):
     thisline = lines[linenum].split()

     if len(thisline) > 0:
       if thisline[0] == 'RAW':
         rawfile = thisline[1]
       elif thisline[0] == 'CON':
         confile = thisline[1]
       elif thisline[0] == 'JSON':
         jsonfile = thisline[1]
       elif thisline[0] == 'DIVISION':
         division = int(thisline[1])
       elif thisline[0] == 'PROCEDURE':
         procedure = thisline[1]
         if procedure == 'EVAL_BASE':
             evalfile = thisline[2]
       elif thisline[0] == 'mykeybus':
           mykeybus = int(thisline[1])
       elif thisline[0] == 'mykeytrdual':
           mykeytrdual = int(thisline[1])
       elif thisline[0] == 'mykeynontrdual':
           mykeynontrdual = int(thisline[1])
       elif thisline[0] == 'END':
           break
       else:
         log.stateandquit("illegal input " + thisline[0] + "\n")
     linenum += 1

    
    for x in [('RAW',rawfile), ('CON',confile), ('JSON',jsonfile), ('PROCEDURE',procedure), ('DIVISION',division), ('mykeybus',mykeybus), ('mykeynontrdual',mykeynontrdual), ('mykeytrdual',mykeytrdual)]:
      if x[1] == 'NONE':
        log.stateandquit(' no ' + x[0] + ' input'+'\n')
      if x[1] != 'NONE': log.joint(x[0] + " : " + str(x[1])+'\n')
      if x[1] == 'EVAL_BASE':
          log.joint("  file " + str(evalfile)+'\n')
          alldata['EVALFILE'] = evalfile
      alldata[x[0]] = x[1]

    #breakexit('')
    return alldata


def run():

    logfile = 'acn.log'
    if len(sys.argv) == 3:
        logfile = sys.argv[2]
    log = danoLogger(logfile)

    stateversion(log)

    if len(sys.argv) == 2:
        configfilename = sys.argv[1]
    else:
        configfilename = '../runs/acn.conf'

    alldata = read_config(log, configfilename)
    alldata['log'] = log
    
    arpae = Solver()
        
    arpae.data.read(alldata['RAW'],alldata['JSON'],alldata['CON'])

    alldata['arpae'] = arpae


    convertcode = acnconvert(alldata)
    acn = alldata['acn']

    if path.exists('solution_BASECASE.txt'):
        os.remove('solution_BASECASE.txt')
    if path.exists('solution1.txt'):
        os.remove('solution1.txt')

    #acnprintgraph(acn, 0)  # 0 means all nodes are black
    
    if alldata['PROCEDURE'] == 'SOLVEAMPL1':
     acnsolverAMPL(alldata)
    elif alldata['PROCEDURE'] == 'COMPUTE_PRIOR_BASE':
     priorcode = acncomputepriorbase(alldata)
    elif alldata['PROCEDURE'] == 'COMPUTE_PRIOR_BASE_AND_EVAL':
     priorcode = acncomputepriorbase(alldata)
    elif alldata['PROCEDURE'] == 'HEURBASE1':
     priorcode = acncomputepriorbase(alldata)
    elif alldata['PROCEDURE'] == 'EVAL_BASE':
     base.obj, base.infeas, base.exist, base.summary = acn_evaluation2.base_from_file(alldata['acn'], alldata['EVALFILE'])
     log.joint(' base.obj = %.8e base.infeas = %.8e \n' %(base.obj, base.infeas))
    else:
     log.joint('\n**** unhandled: PROCEDURE ' + alldata['PROCEDURE'] + '\n')
    log.joint('\ndone\n')
    log.closelog()

if __name__ == '__main__':
    run()
