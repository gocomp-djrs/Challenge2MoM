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
import time

def acnprintgraph(acn, donodecolors):
 log  = acn.alldata['log']
 rawfile = acn.alldata['RAW']
 log.joint("\n***\n printing graph from rawfile " +  rawfile + "\n")

 totviolvector = acn.totviolvector
 thresh = 1e-2

 numbuses = acn.numbuses
 busnumbertocount = acn.busnumbertocount
 f = open('graph.gv','w');
 f.write('graph {\n')
 f.write(' graph[bgcolor=white,\n   overlap=prism,\n   overlap_scaling=-30,\n   size="10!"];\n')
 f.write(' node[color=black,\n   height=0.5,\n   label="\\N",   shape=point,   width=0.5];\n')
 for buscnt1 in range(1, numbuses+1):
    ourbus = acn.ourbuses[buscnt1]
    if donodecolors and totviolvector[buscnt1-1] > thresh:
     log.joint('hey at bus ' + str(buscnt1) + ' totviol = ' + str(totviolvector[buscnt1-1]) + '\n')
     f.write(' ' + str(buscnt1) + ' [color=red,   height=1.0,   width = 1.0];\n')
     #breakexit('heyyyy')
 uniquepair = {}
 uniquecount = 0

 for cnt1 in range(1, 1 + acn.nontranscount_original):
  thisnontrans = acn.ournontrans[cnt1]
  fromcnt = thisnontrans.ourfrombus.count
  tocnt = thisnontrans.ourtobus.count
  if (fromcnt, tocnt) not in uniquepair.values():
   uniquecount += 1
   uniquepair[uniquecount] = (fromcnt, tocnt)
   #log.joint('hit ' + str(cnt1) + ' ' + str(uniquecount) + ' (' + str(fromcnt) + ',' + str(tocnt) + ')' + '\n')
  else:
   #theindex = list(uniquepair.keys()).index((fromcnt, tocnt))
   theindex = list(uniquepair.values()).index((fromcnt,tocnt)) + 1
   #log.joint('oooh ' + str(cnt1)+ ' (' + str(fromcnt) + ',' + str(tocnt) + ') also index ' + str(theindex) + '\n')
   #breakexit('ooh')

 log.joint(' unique count among nontransformers: ' + str(uniquecount) + '\n')

 for cnt1 in range(1, 1 + acn.transcount_original):
  thistrans = acn.ourtrans[cnt1]
  fromcnt = thistrans.ourfrombus.count
  tocnt = thistrans.ourtobus.count
  if ((fromcnt, tocnt) not in uniquepair.values()) and ((tocnt, fromcnt) not in uniquepair.values()):
   uniquecount += 1
   uniquepair[uniquecount] = (fromcnt, tocnt)
   #log.joint('hit ' + str(cnt1) + ' ' + str(uniquecount) + ' (' + str(fromcnt) + ',' + str(tocnt) + ')\n')
  else:
   #theindex = list(uniquepair.keys()).index((fromcnt, tocnt))
   theindex = list(uniquepair.values()).index((fromcnt,tocnt)) + 1
   #log.joint('oooh ' + str(cnt1)+ ' (' + str(fromcnt) + ',' + str(tocnt) + ') also index ' + str(theindex) + '\n')
   #breakexit('ooh')

 log.joint(' unique count also including transformers: ' + str(uniquecount) + '\n')
 for hit in range(1, uniquecount + 1):
  #log.joint(' ' + str(uniquepair[hit][0]) + ' -- ' + str(uniquepair[hit][1]) + ';\n')
  f.write(' ' + str(uniquepair[hit][0]) + ' -- ' + str(uniquepair[hit][1]) + ';\n')
 f.write('}\n')
 #acn.nontranscount_original
 #acn.nontranscount_original

 f.write('\n\n# rawfile: ' + rawfile + '\n')
 f.close()
 #breakexit('printed')
