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
import acn_evaluation2

def trace(frame, event, arg):
    print("%s, %s:%d" % (event, frame.f_code.co_filename, frame.f_lineno))
    return trace


 sol1 = open("solution1.txt","w")

 sol1.write("--bus section\n")
 sol1.write("i,v,theta\n")
 for j in range(acn.numbuses):
     #sol1.write(str(acn.ourbuses[j+1].theirbus.i) + "," + str(basebusVmags[j]) + "," + str(basebusVangles[j]*180/pi) + "\n")
     sol1.write(str(acn.ourbuses[j+1].theirbus.i) + "," + str(basebusVmags[j]) + "," + str(basebusVangles[j]) + "\n")

 sol1.write("--load section\n")
 sol1.write("i,id,t\n")
 countactive = 0
 for count in range(1,1+acn.numloads):
     load = acn.ourloads[count]
     if load.status:
       #sol1.write(str(acn.ourloads[count].theirload.i) + "," + str(acn.ourloads[count].theirload.id) + "," + str(clearedloads[countactive]) + "\n")
       #sol1.write(str(load.theirload.i) + "," + str(load.theirload.id) + "," + str(clearedloads[countactive]) + "\n")
       sol1.write(str(load.theirload.i) + "," + str(load.theirload.id) + "," + str(clearedloads[load.count]) + "\n")
       countactive += 1

 sol1.write("--generator section\n")
 sol1.write("i,id,p,q,x\n")
 count = 0
 for gencount in range(1, acn.numgens+1):
     gen = acn.ourgens[gencount]
     pgen = qgen = 0.0
     genon = basegenON[count]
     if genon:
       pgen = basegenP[count]
       qgen = basegenQ[count]
     count += 1
     #sol1.write(str(gen.i) + "," + str(gen.id) + "," + str(sbase*pgen) + "," + str(sbase*qgen) +  "," + str(genon) + "\n")
     sol1.write(str(gen.i) + "," + str(gen.id) + "," + str(pgen) + "," + str(qgen) +  "," + str(genon) + "\n")


 sol1.write("--line section\n")
 sol1.write("iorig,idest,id,x\n")
 for j in range(acn.nontranscount_original):
     count = j+1
     ntr = acn.ournontrans[count]
     sol1.write(str(ntr.i) + "," + str(ntr.j) + "," + str(ntr.ckt) + "," + str(basentrON[j]) +  "\n")


 sol1.write("--transformer section\n")
 sol1.write("iorig,idest,id,x,xst\n")
 for j in range(acn.transcount_original):
     count = j+1
     tr = acn.ourtrans[count]
     sol1.write(str(tr.i) + "," + str(tr.j) + "," + str(tr.ckt) + "," + str(basetrON[j]) + "," + str(basetaps[j]) + "\n")

 sol1.write("--switched shunt section\n")
 sol1.write("i,xst1,xst2,xst3,xst4,xst5,xst6,xst7,xst8\n")
 for h in range(1,1 + acn.numswitchedshunts):
  oursh = acn.ourswshunts[h]
  if oursh.status:
    sol1.write(str(oursh.i))
    for a in range(1,1+oursh.len):
        if (oursh.n[a]):
            sol1.write("," + str(basesw[h,a]))
    sol1.write("\n")

 sol1.close()

 # Evaluate base solution
 try:
  base.obj, base.infeas, exist, base.summary = acn_evaluation2.base_from_file(acn, "solution1.txt")
  pobj, pinfeas, pexist, base.psummary = acn_evaluation2.base_from_file(acn, "solution_BASECASE.txt")
  print('infeasible solution: obj=' + str(pobj))
  print('infeasible solution: infeas=' + str(pinfeas))
  print('base solution exists=' +  str(exist))
  if (exist == True):
    print('base.obj=' + str(base.obj))
    print('base.infeas=' + str(base.infeas))
  print('pexist = ', pexist)
 except:
  log.joint("Error evaluating the base solution!\n")
  var = traceback.format_exc()
  traceback.print_exc()
  acn_evaluation2.print_alert(var, raise_exception=False)
  pass


##### Called after initial solve to fix all integer variables and re-solve
#Maybe put this in a different file?

def resolve_fixed(acn, ampl):

  log = acn.log
  log.joint('running resolve_fixed\n')

  #pnnl = acn.pnnl
  #raw = pnnl.raw
  #rop = pnnl.rop
  #con = pnnl.con

  #solverdata = acn.solverdata
  #numbuses = acn.numbuses
  #numgens = acn.numgens
  fixtol = 1e-3

  # Round and fix all integer variables in the model

  # Transformer tap settings
  v = ampl.getVariable('xstf0')
  df = v.getValues()
  vals = df.toDict()
  for i in range(1,1+acn.transcount_original):
      tmpval = min(vals[i],acn.ourtrans[i].maxxstf)
      tmpval = max(tmpval, acn.ourtrans[i].minxstf)
      tmpval = int(round(tmpval))
      v[i].fix(tmpval)

  # Switched shunt susceptance
  v = ampl.getVariable('xha0')
  df = v.getValues()
  vals = df.toDict()
  for h in range(1,1 + acn.numswitchedshunts):
      if acn.ourswshunts[h].status > 0:  #switch shunt active
          oursh = acn.ourswshunts[h]
          for a in range(1,1+oursh.len):
              upbnd = int(oursh.n[a])
              if upbnd > 0:
                  xha0val = int(vals[h,a])
                  tmpval = min(xha0val,upbnd)
                  tmpval = max(tmpval, 0)
                  v[h,a].fix(tmpval)
                  #print('h=',h,' a=',a,' xha0val=',xha0val,' upbnd=',upbnd,' fixed=',tmpval)
  #sys.exit('bye')

  # Line switching
  v = ampl.getVariable('ntrON0')
  df = v.getValues()
  vals = df.toDict()
  for i in range(1,1+acn.nontranscount_original):
      if vals[i] < 0.5:
          v[i].fix(0)
      else:
          v[i].fix(1)

  # Transformer switching
  v = ampl.getVariable('trON0')
  df = v.getValues()
  vals = df.toDict()
  for i in range(1,1+acn.transcount_original):
      if vals[i] < 0.5:
          v[i].fix(0)
      else:
          v[i].fix(1)

  # Generator commitment decision variables
  genon = ampl.getVariable('genon0')
  genondf = genon.getValues()
  genonvals = genondf.toDict()
  gensu = ampl.getVariable('genstartup0')
  gensudf = gensu.getValues()
  gesuvals = gensudf.toDict()
  gensd = ampl.getVariable('genshutdown0')
  gensddf = gensd.getValues()
  gesdvals = gensddf.toDict()
  for i in range(1, acn.numgens+1):
      ourgen = acn.ourgens[i]
      # First fix generator on or off and then fix start up/shutdown values
      if genonvals[i] < 0.5:
          genon[i].fix(0)
          if ourgen.status > 0:  # generator currently on; shut it off
              gensu[i].fix(0)
              gensd[i].fix(1)
          else:
              gensu[i].fix(0)
              gensd[i].fix(0)
      else:
          genon[i].fix(1)
          if ourgen.status > 0:
              gensu[i].fix(0)
              gensd[i].fix(0)
          else:                # generator currently off; turn it on
              gensu[i].fix(1)
              gensd[i].fix(0)
  try:
    resolve_options = 'outlev=4 outmode=2 outname="knitro-fixed.log" debug=1 feastol=1e-5 feastol_abs=9e-5 ftol=1e-6 scale=0 honorbnds=1 cg_maxit=50 bar_murule=0 bar_feasible=1 bar_refinement=1 bar_maxcrossit=1 bar_initpi_mpec=0.0 maxit=3000 alg=1 strat_warm_start=1 bar_initmu=1e-6 infeastol=1e-5'
    #resolve_options = 'outlev=4 outmode=2 debug=1 feastol=1e-5 opttol=1e-3 ftol=1e-6 scale=0 honorbnds=0 cg_maxit=50 bar_murule=1 strat_warm_start=1 bar_refinement=1 bar_maxcrossit=1 bar_initpi_mpec=0.0 maxit=100 alg=3'
    mip_options = ' mip_heuristic=0 mip_terminate=1 mip_outinterval=1 mip_outlevel=3 mip_debug=1 mip_outsub=2 mip_nodealg=1 mip_intvar_strategy=0'
    resolve_options += mip_options
    just_feasible = 0
    if just_feasible:
      resolve_options += ' opttol=1.0e30 opttol_abs=1e30'
    else:
      resolve_options += ' opttol=1.0e-3'
    ampl.setOption('knitro_options', resolve_options)
    ampl.solve()
    #breakexit("resolve")
  except:
    log.joint("Error in fixed re-solve!\n")
    pass
  #breakexit("bbb")

#############
