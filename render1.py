#!/usr/bin/python
import os
import sys
import numpy as np
from log import Logger
from myconstants import setconstants
from myutils import *
from reader import *
from actions import *
from versioner import stateversion

def read_config(log, filename):

    log.joint("reading config file " + filename + "\n")

    try:
        f = open(filename, "r")
        lines = f.readlines()
        f.close()
    except IOError as (errno, strerror):
        log.stateandquit("cannot open file " + filename + "\n")

        

    bag = {}
    casefile = 'NONE'
    resultsfilename = 'NONE'
    nodesfilename = 'NONE'
    linenum = 0
    while linenum < len(lines):
        thisline = lines[linenum].split()
        if len(thisline) > 0:
            if thisline[0] == 'CASEFILE':
                casefile = thisline[1]
                log.joint("casefile = " + casefile+"\n")
            elif thisline[0] == 'RESULTSFILE':
                resultsfilename = thisline[1]
                log.joint("resultsfile = " + resultsfilename+"\n")
            elif thisline[0] == 'NODESFILE':
                nodesfilename = thisline[1]
                log.joint("nodesfile = " + nodesfilename+"\n")
            elif thisline[0] == 'END':
                break
            else:
                log.stateandquit("illegal term " + thisline[0] + "\n")
                    
        linenum += 1

    if casefile == 'NONE' or resultsfilename == 'NONE' or nodesfilename == 'NONE':
      log.stateandquit("missing file\n")
    
    bag['casefile'] = casefile
    bag['nodesfilename'] = nodesfilename
    bag['resultsfilename'] = resultsfilename    
    return bag

def process(log, bag):
    buses = bag['buses']
    nodeattr = bag['nodeattr']
    injection = bag['injection'] 
    generation = bag['generation'] 
    negdemand = bag['negdemand']
    maxgeneration = bag['maxgeneration']
    zerogen = bag['zerogen']


    newones = {}
    for bus in nodeattr:
        count = bus.count
        thisattr = nodeattr[bus]

        if generation[count] > 0 or negdemand[count] < 0:
            if generation[count] > 0:
                ratio = generation[count]/maxgeneration
                if ratio >= 0.5:
                    thiscolor = "red"
                elif ratio >= 0.1:
                    thiscolor = "orange"
                else:
                    thiscolor = "yellow"
                lenattr = len(thisattr)
                newstring = thisattr[0:lenattr-2] + " color="+thiscolor+", height=1.0, width=1.0];"
                newones[bus] = newstring
            elif negdemand[count] < 0:
                thiscolor = "green"
                lenattr = len(thisattr)
                newstring = thisattr[0:lenattr-2] + " color="+thiscolor+", height=1.0, width=1.0];"
                newones[bus] = newstring
            else: newones[bus] = thisattr
        else:
            print count, thisattr

            if zerogen[count] == 1:
                lenattr = len(thisattr)
                thiscolor = "white"
                newstring = thisattr[0:lenattr-2] + " color="+thiscolor+", height=1.0, width=1.0];"
                newones[bus] = newstring
            else: newones[bus] = thisattr

        log.joint("   " + str(bus.count) + " " + newones[bus] + "\n")


                
                
        breakexit
    

def getnodes(log, bag):
    filename = bag['nodesfilename']
    log.joint("reading nodes file " + filename + "\n")
    buses = bag['buses']
    try:
        f = open(filename, "r")
        lines = f.readlines()
        f.close()
    except IOError as (errno, strerror):
        log.stateandquit("cannot open file " + filename + "\n")
    linenum = 0
    nodeattr = {}
    while linenum < len(lines):
        thisline = lines[linenum].split()
        if len(thisline) > 0:
            if(thisline[0] == "END"):
                break
            else:
                nodecount = int(thisline[0])
                if nodecount < 1 or nodecount > len(buses):
                    breakexit("illegal node count " + str(nodecount))
                nodeattr[buses[nodecount]] = thisline[1]
        linenum += 1
    log.joint("done reading nodes\n")

    
    bag['nodeattr'] = nodeattr
def getinjections(log, bag):
    filename = bag['resultsfilename']
    log.joint("reading results file " + filename + "\n")
    buses = bag['buses']
    baseMVA = bag['baseMVA']

    injection = [0 for j in xrange(len(buses)+1)]
    generation = [0 for j in xrange(len(buses)+1)]
    negdemand = [0 for j in xrange(len(buses)+1)]
    zerogen = [0 for j in xrange(len(buses)+1)]

    try:
        f = open(filename, "r")
        lines = f.readlines()
        f.close()
    except IOError as (errno, strerror):
        log.stateandquit("cannot open file " + filename + "\n")

    
    for bus in buses.values():
        thisdemand = bus.Pd*baseMVA
        if bus.Pd < 0:
 #           log.joint("bus " + str(bus.count) + " has neg demand " + str(thisdemand) + "\n")
            negdemand[bus.count] = thisdemand
    linenum = 0
    while linenum < len(lines):
        thisline = lines[linenum].split()
        if len(thisline) > 0:
            if(thisline[0][0:4] == "gen"):
                print thisline
                gennumber = int(thisline[2])
                busnumber = int(thisline[4])
                thisgeneration = float(thisline[6])
                generation[busnumber] += thisgeneration
                if thisgeneration == 0:
#                    breakexit("ha")
                    zerogen[busnumber] = 1
#                log.joint("g " + str(gennumber) + " " + str(busnumber) + " "  + str(thisgeneration) + "\n")
#                if negdemand[busnumber] < 0:
#                    breakexit("hoooo")
        linenum += 1

    maxnegdemand = maxgeneration = 0.0
    indmaxnegdemand = indmaxgeneration = -1
    for j in xrange(len(buses)+1):
        if generation[j] > maxgeneration:
            maxgeneration = generation[j]
            indmaxgeneration = j
        if negdemand[j] < maxnegdemand:
            maxnegdemand = negdemand[j]
            indmaxnegdemand = j

    log.joint(" max generation " + str(maxgeneration) + " at " + str(indmaxgeneration) + "\n")
    log.joint(" max negdemand " + str(maxnegdemand) + " at " + str(indmaxnegdemand) + "\n")

    bag['injection'] = injection
    bag['generation'] = generation
    bag['negdemand'] = negdemand
    bag['maxgeneration'] = maxgeneration
    bag['maxnegdemand'] = maxnegdemand
    bag['zerogen'] = zerogen

###### main

if len(sys.argv) != 2 and len(sys.argv) !=3 :
  print 'Usage: driver.py configfile [logfile]'
  exit(0)

if len(sys.argv) ==3 :
  logfile = sys.argv[2]
else:
  logfile = "render1.log"

log = Logger(logfile)

stateversion(log)

bag = read_config(log, sys.argv[1])
readcase(log, bag, bag['casefile'])
setconstants(log,bag)

#for branch in bag['branches'].values():
#    log.joint(" " + str(branch.id_f) + " -- " + str(branch.id_t) +";\n")
#breakexit("")

getinjections(log, bag)

getnodes(log, bag)

process(log, bag)

log.closelog()



