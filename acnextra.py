import csv
import sys, os, shutil
import traceback
import math
import cmath
from myutils import *
from log import *
import data
import numpy as np
import scipy
from acn import *
import time
import acn_evaluation2
from acnextra import *


def acn_dijkstra(acn, ctglabel, appliestoall, ntrout0, trout0, sourcebuscount1, exceptions, numexceptions, maxexceptions, alsohitgens):
    log = acn.log
    log.joint(' dijkstra for ntrout0 %d trout0 %d buscount1 %d\n' %(ntrout0, trout0, sourcebuscount1))
    log.joint('   maxexceptions %d appliestoall %i alsohitgens %i\n' %(maxexceptions,appliestoall,alsohitgens))
    time1 = time.time()
    if not acn.Dijkstrasetup: #Dijkstrasetup

        log.joint('Dijkstra not yet set up\n')

        arccount = acn.nontranscount + acn.transcount

        row_ind = np.empty(arccount, dtype = int)
        col_ind = np.empty(arccount, dtype = int)
        w = np.empty(arccount)
        pairs = {}

        arccount = 0

        for h in range(acn.nontranscount):
            ournontrans = acn.nontransATbuses[h]
            parent = ournontrans.parent
            xvalue = parent.x
            if xvalue < 0:
                xvalue = 0
            rvalue = 1
            u = ournontrans.ourfrombus.count-1
            v = ournontrans.ourtobus.count-1
            #print(arccount,'N',h, 'f', u, 't', v, xvalue, rvalue)

            row_ind[arccount] = u
            col_ind[arccount] = v
            w[arccount] = xvalue/rvalue


            arccount += 1

        arccount_pretrans = arccount
        acn.arccount_pretrans = arccount_pretrans

        for h in range(acn.transcount):
            ourtrans = acn.transATbuses[h]
            parent = ourtrans.parent
            if parent.Bft > 1e-5:
                xvalue = np.absolute(1/parent.Bft)
            else:
                xvalue = 0
            rvalue = 1
            u = ourtrans.ourfrombus.count-1
            v = ourtrans.ourtobus.count-1
            #print(arccount,'t',h, 'f', u, 't', v, xvalue, rvalue)

            row_ind[arccount] = u
            col_ind[arccount] = v
            w[arccount] = xvalue/rvalue
            arccount += 1

        acn.row_ind = row_ind
        acn.col_ind = col_ind
        acn.w = w
        acn.arccount = arccount
        acn.Dijkstrasetup = True

        time2 = time.time()
        log.joint(' Dijkstra setup time %g\n' %(time2 - time1))
    else:
        log.joint('Dijkstra already set up\n')
        row_ind = acn.row_ind
        col_ind = acn.col_ind
        w = acn.w
        arccount = acn.arccount
        arccount_pretrans = acn.arccount_pretrans


    #in case of lineout or transout, modify its weight

    branchout0 = -1
    if ntrout0 >= 0:
        branchout0 = ntrout0
        partnercount0 = acn.nontransATbuses[ntrout0].partner.count1 - 1        
    if trout0 >= 0:
        branchout0 = arccount_pretrans + branchout0
        partnercount0 = arccount_pretrans + acn.transATbuses[trout0].partner.count1 - 1        
    if branchout0 >= 0:
        infinity = 1e20
        oldweight = w[branchout0]
        w[branchout0] = infinity
        oldweightpartner = w[partnercount0]
        w[partnercount0] = infinity
        log.joint(' set weight of branch %d to %g and partner %d to %g\n' %(branchout0, w[branchout0], partnercount0, w[partnercount0]))


    '''
    for j in range(arccount): 
        log.joint('%d ( %d , %d ) -> %g\n' %(j, row_ind[j], col_ind[j], w[j]))
    '''


    time2 = time.time()
    time2 = time.time()
    log.joint(' map time: ' + str(time2 - time1) + '\n')

    time1 = time.time()

    mat_coo = scipy.sparse.coo_matrix((w, (row_ind, col_ind)))

    #print(mat_coo)
    #breakexit('dijk')

    time2 = time.time()
    log.joint(' coo time: ' + str(time2 - time1) + '\n')

    #breakexit('coo')

    time1 = time.time()
    distances, pred = scipy.sparse.csgraph.dijkstra(mat_coo, indices = sourcebuscount1-1, return_predecessors = True)
    time2 = time.time()

    #print( "distances", distances)
    #print( "number of distances", len(distances))
    #print predecessors

    log.joint(' Dijkstra ran in time %g\n' %(time2-time1))

     #restore w!!!
    if branchout0 >= 0:
     w[branchout0] = oldweight
     w[partnercount0] = oldweightpartner
     log.joint(' restored weight of line %d to %g and partner %d to %g\n' %(branchout0, w[branchout0], partnercount0, w[partnercount0]))


    sortedorder = numpy.argsort(distances)
    numhit = 0
    numgenhit = 0

    neededgens = 1*alsohitgens

    numbuses = acn.numbuses


    numexceptions = 0
    
    for h in range(numbuses):
     j = sortedorder[h]
     ourbus = acn.ourbuses[j+1]

     if appliestoall or ourbus.gencount > 0:
      notinex = ourbus.count not in exceptions.values()
      if notinex:
        log.joint(' numhits %d numgenhit %d numexceptions %d sorder %d j %d bus.i %d bus.count %d distance %g pred %d ' %(numhit, numgenhit, numexceptions, h, j, ourbus.i, ourbus.count, distances[j], pred[j]))
        log.joint(" <<\n")
        exceptions[numexceptions] = ourbus.count
        numexceptions += 1

        numhit += 1

        if ourbus.gencount > 1 or ourbus.count != sourcebuscount1:
          numgenhit += 1
          log.joint(" hit bus count %i\n" %(ourbus.count))

      
        if numhit >= maxexceptions and numgenhit >= neededgens:
          log.joint(" ---> hit limit of %d exceptions (%d gens hit)\n" %(numhit,numgenhit))
          h = numbuses+1
          break


    #breakexit('finallydijk')

    return numexceptions
