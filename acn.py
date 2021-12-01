# Main data structure and utilities

import sys
import math
import cmath
import numpy

from myutils import *

from log import *

class ACN:
 def __init__(self, alldata):
  log = self.log = alldata['log']
  self.arpae = alldata['arpae']
  raw = self.raw = self.arpae.data.raw
  con = self.con = self.arpae.data.con
  sup = self.sup = self.arpae.data.sup
  self.alldata = alldata
  self.division = alldata['DIVISION']
  self.Dijkstrasetup = False
  self.row_ind = self.col_ind = self.w = self.arccount = None

  self.options = {}
  self.options['sbase'] = sbase = raw.case_identification.sbase
  self.options['epsilon_3'] = 1.0e-3
  self.options['epsilon_5'] = 1.0e-5
  self.options['epsilon_6'] = 1.0e-6
  self.options['delta'] = sup.sup_jsonobj["systemparameters"]["delta"]
  self.options['deltactg'] = sup.sup_jsonobj["systemparameters"]["deltactg"]
  self.options['deltar'] = sup.sup_jsonobj["systemparameters"]["deltar"]
  self.options['deltarctg'] = sup.sup_jsonobj["systemparameters"]["deltarctg"]
  log.joint("\nsbase = " + str(sbase) + "\n")

  log.joint('delta = %g\n' %(self.options['delta']))

  self.numcontingencies = len(con.contingencies)
  log.joint('number of contingencies = %d\n' %(self.numcontingencies))


  self.maxtime = 1e6
  self.timebeg = 0
  self.numbuses = len(self.raw.buses)
  self.numgens = len(self.raw.generators)
  self.numloads = len(self.raw.loads)
  self.totalPloadpu = self.totalQloadpu = 0
  self.totalPloadabspu = self.totalQloadabspu = 0
  self.numactiveloadsat0 = 0
  self.numfixedshunts = len(self.raw.fixed_shunts)
  self.numswitchedshunts = len(self.raw.switched_shunts)
  self.numactiveloads = self.raw.num_loads_active
  self.numvariabletransformers = 0
  self.numimpedancecorrectedtransformers = 0
  self.numtautransformers = 0
  self.numthetatransformers = 0

  self.pimb = numpy.zeros(self.numbuses)
  self.qimb = numpy.zeros(self.numbuses)

  self.pimb_abs = numpy.zeros(self.numbuses)
  self.qimb_abs = numpy.zeros(self.numbuses)
  self.pimb_penalty = numpy.zeros(self.numbuses)
  self.qimb_penalty = numpy.zeros(self.numbuses)


  self.totviolvector = numpy.zeros(self.numbuses)
  self.heatmap = numpy.zeros(self.numbuses)

  self.actionable1 = numpy.zeros(1+self.numbuses, dtype = int)


  #we have two types: nontrans's and nontransAtbus's (the latter is a child of the
  #former
  self.nontransATbuses = {}
  self.nontranscount = 0 #2x count of non-trans branches (once for each end)
  self.nontranscount_original = 0 # count of non-trans branches

  self.transATbuses = {}
  self.transcount = 0 #2x count of trans branches (once for each end)
  self.transcount_original = 0 # count of trans branches



  busnumbertocount = self.busnumbertocount = {}
  buscounttonumber = self.buscounttonumber = {}

  log.joint(' buses %i loads %i (%i active) generators %i fixedshunts %i switched %i\n' %(self.numbuses, self.numloads, self.numactiveloads, self.numgens, self.numfixedshunts, self.numswitchedshunts))

  log.joint(' mapping buses\n')
  ourbuses = self.ourbuses = {}
  count = 1
  for theirbus in raw.buses. values():
   ourbuses[count] = ACNbus(log, self, theirbus, sbase, count) #counting from 1
   busnumbertocount[theirbus.i] = count
   buscounttonumber[count] = theirbus.i
   #ourbuses[count].showbus_slim_nocr(log)
   count += 1
  log.joint(' mapped %i buses\n' %(count-1))

  log.joint(' mapping loads\n')

  ourloads = self.ourloads = {}
  count = 1
  for theirload in raw.get_loads():
   ourloads[count] = ACNload(log, sup, theirload, sbase, ourbuses[busnumbertocount[theirload.i]], count)
   count += 1
  log.joint(' mapped %i loads %i active\n' %(count-1, self.numactiveloadsat0))

  log.joint(' total Ploadpu: %g, Qloadplu: %g\n' %(self.totalPloadpu,self.totalQloadpu))
  log.joint(' total Ploadabspu: %g, Qloadabsplu: %g\n' %(self.totalPloadabspu,self.totalQloadabspu))

  #print(' num sup loads ', len(sup.loads))
  #for theirload in raw.get_loads():

  #print(' theirloadkey ', key, ' lencblocks', len(sup.loads[key]['cblocks']))

  #breakexit("")

  log.joint(' mapping %d fixed shunts\n' %(self.numfixedshunts))
  ourfxshunts = self.ourfxshunts = {}
  count = 1
  for fs in  raw.fixed_shunts.values():
   buscount = busnumbertocount[fs.i]
   ourbus = ourbuses[buscount]
   if fs.status:
    ourfxshunts[count] = ACNfxshunt(log, fs, ourbus, sbase, count)
    ourbus.fsglpu += fs.gl/sbase
    ourbus.fsblpu += fs.bl/sbase
    #print(fs.i, buscount, fs.bl)
    count += 1
  log.joint(' mapped %i active fixed shunts\n' %(count-1))


  log.joint(' mapping generators\n')
  generators = raw.generators
  ourgens = self.ourgens = {}
  gencounttonumberid = self.gencounttonumberid = {}
  gennumberidtocount = self.gennumberidtocount = {}
  gencounttogen = self.gencounttogen = {}
  count = 1
  for theirgen in generators.values():
   buscount = busnumbertocount[theirgen.i]
   cblocks = sup.generators[(theirgen.i, theirgen.id)]['cblocks']
   ourgen = ourgens[count] = ACNgen(log, sup, theirgen, sbase, ourbuses[buscount], count)
   ourbuses[buscount].addgenerator(log,cblocks,ourgen)
   gennumberidtocount[ourgen.i,ourgen.id] = count
   gencounttonumberid[count] = (ourgen.i, ourgen.id)
   gencounttogen[count] = ourgen

   #print(self.ourbuses[581].gencount)
   count += 1
  log.joint(' mapped %i generators\n' %(count-1))


  numnontrans = self.numnontrans = len(raw.nontransformer_branches)
  numtrans = self.numtrans = len(raw.transformers)

  self.ntrmagviol = numpy.zeros(self.numnontrans)
  self.trmagviol = numpy.zeros(self.numtrans)


  log.joint(' nontransformerbranches %i transformers %i\n' %(self.numnontrans, self.numtrans))

  log.joint(' mapping nontransformerbranches\n')
  nontransid = self.nontransid = {}
  nontransfromid = self.nontransfromid = {}
  ournontrans = self.ournontrans = {}
  count1 = 1
  for theirnontrans in raw.nontransformer_branches.values():
   fromcount = busnumbertocount[theirnontrans.i]
   tocount = busnumbertocount[theirnontrans.j]
   nontransid[theirnontrans] = count1
   nontransfromid[count1] = theirnontrans

   ournon = ournontrans[count1] = ACNnontrans(log, self, sup, theirnontrans, ourbuses[fromcount], ourbuses[tocount], sbase, count1)
   ourbuses[fromcount].addfromnontrans(log, ournon)
   newfromAT = self.nontransATbuses[self.nontranscount-1]
   ourbuses[tocount].addtonontrans(log, ournon)
   newtoAT = self.nontransATbuses[self.nontranscount-1]

   newfromAT.partner = newtoAT
   newtoAT.partner = newfromAT

   count1 += 1


  #breakexit('pooooooo')
  log.joint(' mapped %i nontransbranches\n' %(count-1))
  log.joint(' number of nontransATbuses: %i\n' %(self.nontranscount))
  log.joint(' number of nontransATbuses (original): %i\n' %(self.nontranscount_original))

  self.ticts={}
  count = 1
  foo = raw.get_transformer_impedance_correction_tables()
  log.joint(' mapping %d TICT records\n' %(len(foo)))
  for r in foo:
   #log.joint('  TICT %i with count %i \n' %(r.i, r.tict_point_count))
   self.ticts[count] = ACNtict(log, r, count)
   count += 1
  self.numticts = count -1
  log.joint(' numticts: %d\n' %(self.numticts))

  log.joint(' mapping transformers\n')
  self.numtrans = len(raw.transformers)
  transid = self.transid = {}
  transfromid = self.transfromid = {}
  ourtrans = self.ourtrans = {}
  count1 = 1
  for theirtrans in raw.transformers.values():
      fromcount = busnumbertocount[theirtrans.i]
      tocount = busnumbertocount[theirtrans.j]
      transid[theirtrans] = count1
      transfromid[count1] = theirtrans
      ours = ourtrans[count1] = ACNtrans(log, self, sup, theirtrans, ourbuses[fromcount], ourbuses[tocount], sbase, count1)

      ourbuses[fromcount].addfromtrans(log,ours)

      newfromAT = self.transATbuses[self.transcount-1]
      ourbuses[tocount].addtotrans(log,ours)
      newtoAT = self.transATbuses[self.transcount-1]

      newfromAT.partner = newtoAT
      newtoAT.partner = newfromAT

      #print(count1,newfromAT, newtoAT)

      #ourtrans.showtrans_slim(log)
      count1 += 1

  log.joint(' mapped %i transformers, %d variable, %d impedancecorr %d tau %d theta\n' %(count-1, self.numvariabletransformers, self.numimpedancecorrectedtransformers, self.numtautransformers, self.numthetatransformers))

  log.joint(' mapping %d switched shunts\n' %(self.numswitchedshunts))

  count = 1
  ourswshunts = self.ourswshunts = {}
  for ss in  raw.switched_shunts.values():
    buscount = busnumbertocount[ss.i]
    ourswshunts[count] = ACNswshunt(log, ss, ourbuses[buscount], sbase, count)
    count += 1
  log.joint(' mapped %i switched shunts\n' %(count-1))


  log.joint(' mapping blocks\n')
  numpcblocks = self.numpcblocks = len(sup.sup_jsonobj['pcblocks'])
  numqcblocks = self.numqcblocks = len(sup.sup_jsonobj['qcblocks'])
  numscblocks = self.numscblocks = len(sup.sup_jsonobj['scblocks'])

  log.joint(' num pcblocks %d, qcblocks %d, scblocks %d\n' %(numpcblocks, numqcblocks, numscblocks))

  self.pcostcblock = numpy.zeros(1+self.numpcblocks)
  self.pmaxcblock = numpy.zeros(1+self.numpcblocks)
  self.pcblock1 = {}
  self.pcblock1['c'] = [0 for j in range(1 + numpcblocks)]
  self.pcblock1['pmax'] = [0 for j in range(1 + numpcblocks)]
  j = 1
  for pcblock in sup.sup_jsonobj['pcblocks']:
   for x in ('c', 'pmax'):
    self.pcblock1[x][j] = pcblock[x]
   self.pcostcblock[j] = self.pcblock1['c'][j] * sbase
   self.pmaxcblock[j] = self.pcblock1['pmax'][j] / sbase
   log.joint('   pcblock # %d c %f pmax %f\n' %(j, self.pcblock1['c'][j], self.pcblock1['pmax'][j]))
   j += 1

  self.qcostcblock = numpy.zeros(1+self.numqcblocks)
  self.qmaxcblock = numpy.zeros(1+self.numqcblocks)
  self.qcblock1 = {}
  self.qcblock1['c'] = [0 for j in range(1 + numqcblocks)]
  self.qcblock1['qmax'] = [0 for j in range(1 + numqcblocks)]
  j = 1
  for qcblock in sup.sup_jsonobj['qcblocks']:
   for x in ('c', 'qmax'):
    self.qcblock1[x][j] = qcblock[x]
   self.qcostcblock[j] = self.qcblock1['c'][j] * sbase
   self.qmaxcblock[j] = self.qcblock1['qmax'][j] / sbase
   log.joint('   qcblock # %d c %f qmax %f\n' %(j, self.qcblock1['c'][j], self.qcblock1['qmax'][j]))
   j += 1

  self.scostcblock = numpy.zeros(1+self.numscblocks)
  self.smaxcblock = numpy.zeros(1+self.numscblocks)
  self.scblock1 = {}
  self.scblock1['c'] = [0 for j in range(1 + numscblocks)]
  self.scblock1['tmax'] = [0 for j in range(1 + numscblocks)]
  j = 1
  for scblock in sup.sup_jsonobj['scblocks']:
   for x in ('c', 'tmax'):
    self.scblock1[x][j] = scblock[x]
   self.scostcblock[j] = self.scblock1['c'][j] * sbase
   self.smaxcblock[j] = self.scblock1['tmax'][j]
   log.joint('   scblock # %d c %f tmax %f\n' %(j, self.scblock1['c'][j], self.scblock1['tmax'][j]))
   j += 1


class ACNbus:
     def __init__(self, log, acn, theirbus, sbase, count):
      self.theirbus = theirbus
      # fields such as nvhi, nvhlo, etc accessible through 'theirbus'
      self.count = count
      self.i = theirbus.i
      self.sbase = sbase
      self.gens1 = {}
      self.gencostblock1 = {}
      self.numgencostblock1= {}
      self.gencount = 0
      self.gencountat0 = 0
      self.acn = acn

      self.fromnontrans = {}
      self.tonontrans = {}
      self.starnontrans = {}
      self.outnontransdegree = self.innontransdegree = self.nontransdegree = 0

      self.fromtrans = {}
      self.totrans = {}
      self.startrans = {}
      self.outtransdegree = self.intransdegree = self.transdegree = 0

      self.loads = {}
      self.loadcount = 0
      self.activeloadcountat0 = 0
      self.Ploadpu = self.Qloadpu = 0
      self.fsglpu = self.fsblpu = 0
      self.vm = theirbus.vm
      self.va = theirbus.va*math.pi/180.0
      self.nvhi = theirbus.nvhi
      self.nvlo = theirbus.nvlo
      self.evhi = theirbus.evhi
      self.evlo = theirbus.evlo

      self.nvahi = self.va + 100*math.pi
      self.nvalo = self.va - 100*math.pi  #crude bounds

      self.fxshunts = {}
      self.numfxshunts = 0
      self.swshunts = {}
      self.activeswshunts = {}
      self.numswshunts = 0
      self.numactiveswshunts = 0

      #log.joint(' added bus count %i i %i vm %g va %g\n' %(self.count,self.i, self.vm, self.va))

     def addload(self, log, ourload):
      self.loads[self.loadcount] = ourload
      ourload.inbuscount0 = self.loadcount
      self.loadcount += 1
      self.activeloadcountat0 += ourload.status
      self.acn.numactiveloadsat0 += ourload.status
      self.Ploadpu += ourload.plpu
      self.Qloadpu += ourload.qlpu

      self.acn.totalPloadpu += ourload.plpu
      self.acn.totalQloadpu += ourload.qlpu

      self.acn.totalPloadabspu += numpy.absolute(ourload.plpu)
      self.acn.totalQloadabspu += numpy.absolute(ourload.qlpu)

      #log.joint(' adding load > %i %s < count %d ' %(ourload.i, ourload.id, ourload.count))
      #log.joint('inbuscount0 %d busnumber %i\n' %(ourload.inbuscount0, self.i))

     def addfxshunt(self, log, ourfxshunt):
      self.fxshunts[self.numfxshunts] = ourfxshunt
      ourfxshunt.inbuscount0 = self.numfxshunts
      self.numfxshunts += 1

      #log.joint(' adding fxshunt > %i %s < count %d ' %(ourfxshunt.i, ourfxshunt.id, ourfxshunt.count))
      #log.joint('inbuscount0 %d busnumber %i\n' %(ourfxshunt.inbuscount0, self.i))

     def addswshunt(self, log, ourswshunt):
      self.swshunts[self.numswshunts] = ourswshunt
      ourswshunt.inbuscount0 = self.numswshunts
      self.numswshunts += 1
      if ourswshunt.status:
       self.activeswshunts[self.numactiveswshunts] = ourswshunt
       self.numactiveswshunts += 1

      #log.joint(' adding fxshunt > %i %s < count %d ' %(ourfxshunt.i, ourfxshunt.id, ourfxshunt.count))
      #log.joint('inbuscount0 %d busnumber %i\n' %(ourfxshunt.inbuscount0, self.i))



     def addfromnontrans(self, log, nontrans):
      self.fromnontrans[self.outnontransdegree] = nontrans

      self.outnontransdegree += 1
      self.starnontrans[self.nontransdegree] = newfrom = ACNnontransATbus(log,nontrans,1,self.acn.nontranscount)
      #log.joint(" nontrans " + str(newfrom.count1) + " added_from at bus " + str(self.count) +"\n")
      self.nontransdegree += 1
      self.acn.nontransATbuses[self.acn.nontranscount] = newfrom
      #print('added from non trans %d\n' %(self.acn.nontranscount))
      self.acn.nontranscount += 1
      self.acn.nontranscount_original += newfrom.original_or_reversed


     def addtonontrans(self, log, nontrans):
      self.tonontrans[self.innontransdegree] = nontrans
      self.innontransdegree += 1
      self.starnontrans[self.nontransdegree] = newto = ACNnontransATbus(log,nontrans,0,self.acn.nontranscount)
      #log.joint(" nontrans " + str(newto.count1) + " added_to at   bus " + str(self.count) +"\n")
      self.nontransdegree += 1
      self.acn.nontransATbuses[self.acn.nontranscount] = newto
      self.acn.nontranscount += 1
      self.acn.nontranscount_original += newto.original_or_reversed


     def addfromtrans(self, log, trans):
      self.fromtrans[self.outtransdegree] = trans
      self.outtransdegree += 1
      self.startrans[self.transdegree] = newfrom = ACNtransATbus(log,trans,1,self.acn.transcount)
      self.transdegree += 1
      self.acn.transATbuses[self.acn.transcount] = newfrom
      self.acn.transcount += 1
      #log.joint(" trans " + str(newfrom.count1) + " added_from at bus " + str(self.count) +"\n")
      #log.joint('    transcount now ' + str(self.acn.transcount) + '\n')
      self.acn.transcount_original += newfrom.original_or_reversed


     def addtotrans(self, log, trans):
      self.totrans[self.intransdegree] = trans
      self.intransdegree += 1
      self.startrans[self.transdegree] = newto = ACNtransATbus(log,trans,0,self.acn.transcount)
      self.transdegree += 1
      self.acn.transATbuses[self.acn.transcount] = newto
      self.acn.transcount += 1
      #log.joint(" trans " + str(newto.count1) + " added_to  at bus " + str(self.count) +"\n")
      #log.joint('    transcount now ' + str(self.acn.transcount) + '\n')
      self.acn.transcount_original += newto.original_or_reversed


     def addgenerator(self, log, cblocks, ourgen):
      self.gencount += 1
      self.gencountat0 += ourgen.status
      self.gens1[self.gencount] = ourgen
      self.gencostblock1[self.gencount] = {}
      numcblocks = self.numgencostblock1[self.gencount] = len(cblocks)
      #log.joint('  at bus count %i i %i adding gen count %i with i %i id %s num cblocks %d now gencount is %d\n' %(self.count, self.i, ourgen.count, ourgen.i, ourgen.id, numcblocks, self.gencount))


      self.gencostblock1[self.gencount]['c'] = [0 for j in range(1+numcblocks)]
      self.gencostblock1[self.gencount]['pmax'] = [0 for j in range(1+ numcblocks)]

      j = 1
      foo = self.gencostblock1[self.gencount]
      for cblock in cblocks:
       foo['c'][j] = cblock['c'] #self.gencostblock1[self.gencount]['c'][j] = cblock['c']
       foo['pmax'][j] = cblock['pmax'] #self.gencostblock1[self.gencount]['pmax'][j] = cblock['pmax']
       #log.joint('   gen cblock # %d c %f pmax %f\n' %(j, self.gencostblock1[self.gencount]['c'][j], self.gencostblock1[self.gencount]['pmax'][j]))
       j += 1


     def showbus_slim_nocr(self, log):
      log.joint("(count1 %i busi %i vm %g va %g)" %(self.count, self.theirbus.i, self.vm, self.va))

class ACNload:
 def __init__(self, log, sup, theirload, sbase, ourbus, count):
  self.theirload = theirload
  self.i = theirload.i
  self.id = theirload.id
  self.plpu = theirload.pl/sbase
  self.qlpu = theirload.ql/sbase

  self.status = theirload.status
  self.count = count
  #if theirload.status == 0:
  # log.joint('load count %d not active\n' %(count))


  self.costblock1 = {}
  self.theirkey = theirkey = (theirload.i, theirload.id)
  self.tmin = sup.loads[theirkey]['tmin']
  self.tmax = sup.loads[theirkey]['tmax']
  self.pru = sup.loads[theirkey]['prumax']/sbase
  self.prd = sup.loads[theirkey]['prdmax']/sbase
  self.pructg = sup.loads[theirkey]['prumaxctg']/sbase
  self.prdctg = sup.loads[theirkey]['prdmaxctg']/sbase
  self.numcblocks = len(sup.loads[theirkey]['cblocks'])

  #log.joint(' creating ACN load with i %i id %s status %i num_cblocks %d\n' %(self.i, self.id, self.status, self.numcblocks))

  #self.costblock1['c'] = [0 for j in range(1+self.numcblocks)]
  #self.costblock1['pmax'] = [0 for j in range(1+self.numcblocks)]
  self.costcblock = numpy.zeros(1+self.numcblocks)
  self.maxcblock = numpy.zeros(1+self.numcblocks)
  j = 1
  for cblock in sup.loads[theirkey]['cblocks']:
   self.costcblock[j] = cblock['c']*sbase
   self.maxcblock[j] = cblock['pmax']/sbase
   #log.joint('  cblock # %d c %f pmax %f\n' %(j, self.costcblock[j], self.maxcblock[j]))
   j += 1

  ourbus.addload(log, self)
  #self.showload(log)

 def showload(self, log):
  log.joint(' load count1 ' + str(self.count) + ' key ' + str(self.theirkey) + ' stat ' + str(self.status))
  log.joint(' tmin ' +str(self.tmin) + ' tmax ' + str(self.tmax))
  log.joint(' plpu ' + str(self.plpu) + ' qlpu ' + str(self.qlpu) + '\n')


class ACNfxshunt:
 def __init__(self, log, theirfxshunt, ourbus, sbase, count):
  self.theirfxshunt = theirfxshunt
  self.i = theirfxshunt.i
  self.id = theirfxshunt.id
  self.glpu = theirfxshunt.gl/sbase
  self.blpu = theirfxshunt.bl/sbase

  self.status = theirfxshunt.status
  self.count = count

  #log.joint('  creating fxshunt count %d with i %i id %s status %d\n' %(count, self.i, self.id, self.status))

  ourbus.addfxshunt(log, self)

class ACNswshunt:
 def __init__(self, log, theirswshunt, ourbus, sbase, count):
  self.theirswshunt = theirswshunt
  self.i = theirswshunt.i
  self.status = theirswshunt.stat
  self.ourbus = ourbus
  self.count = count
  self.len = 8
  #log.joint('  creating swshunt count %d with i %i len %d\n' %(count, theirswshunt.i, self.len))

  self.binitpu = theirswshunt.binit/sbase

  self.b = {}
  self.b[1] = theirswshunt.b1/sbase
  self.b[2] = theirswshunt.b2/sbase
  self.b[3] = theirswshunt.b3/sbase
  self.b[4] = theirswshunt.b4/sbase
  self.b[5] = theirswshunt.b5/sbase
  self.b[6] = theirswshunt.b6/sbase
  self.b[7] = theirswshunt.b7/sbase
  self.b[8] = theirswshunt.b8/sbase

  self.n = {}
  self.n[1] = theirswshunt.n1
  self.n[2] = theirswshunt.n2
  self.n[3] = theirswshunt.n3
  self.n[4] = theirswshunt.n4
  self.n[5] = theirswshunt.n5
  self.n[6] = theirswshunt.n6
  self.n[7] = theirswshunt.n7
  self.n[8] = theirswshunt.n8

  for j in range(1,9):
   if self.n[j] <= 0:
    break

  if j == 8:
   self.sizeAh = 8
  else:
   self.sizeAh = j-1 #last positive

  ourbus.addswshunt(log, self)

  #if count == xxx:
  # breakexit('did xxx of size ' + str(self.sizeAh))
class ACNgen:
 def __init__(self, log, sup, theirgen, sbase, ourbus, count):
  self.i = theirgen.i
  self.id = theirgen.id
  self.ourbus = ourbus
  self.theirgen = theirgen
  self.status = theirgen.stat #should be on or off in original solution
  self.count = count
  self.pg = theirgen.pg/sbase
  self.qg = theirgen.qg/sbase
  self.maxpg = theirgen.pt/sbase
  self.minpg = theirgen.pb/sbase
  self.maxqg = theirgen.qt/sbase
  self.minqg = theirgen.qb/sbase
  self.area = ourbus.theirbus.area
  self.sup = sup
  self.theirjsongen = theirjsongen = sup.generators[(theirgen.i, theirgen.id)]

  self.prumax = theirjsongen['prumax']/sbase
  self.prdmax = theirjsongen['prdmax']/sbase
  self.prumaxctg = theirjsongen['prumaxctg']/sbase
  self.prdmaxctg = theirjsongen['prdmaxctg']/sbase
  self.oncost = theirjsongen['oncost']

  self.sucost = theirjsongen['sucost']
  self.sdcost = theirjsongen['sdcost']
  self.suqual = theirjsongen['suqual']
  self.sdqual = theirjsongen['sdqual']
  self.suqualctg = theirjsongen['suqualctg']
  self.sdqualctg = theirjsongen['sdqualctg']
  self.cblocks = theirjsongen['cblocks']
  self.numcblocks = len(self.cblocks)

  self.pgcostcblock = numpy.zeros(1+self.numcblocks)
  self.pgmaxcblock = numpy.zeros(1+self.numcblocks)
  j = 1
  for cblock in self.cblocks:
   self.pgcostcblock[j] = cblock['c']*sbase
   self.pgmaxcblock[j] = cblock['pmax']/sbase
   j += 1




  #log.joint(' creating ACN gen count %d with i %i id %s status %i pg %g qg %g\n' %(self.count, self.i, self.id, self.status, self.pg, self.qg))

class ACNnontrans:
     def __init__(self, log, acn, sup, theirnontrans, ourfrombus, ourtobus, sbase, count1):
         self.acn = acn
         self.ournontrans = theirnontrans
         self.ourfrombus = ourfrombus
         self.ourtobus = ourtobus
         self.ckt = theirnontrans.ckt
         self.count1 = count1
         self.i = theirnontrans.i
         self.j = theirnontrans.j
         self.st = theirnontrans.st
         self.key = (self.i,self.j,self.ckt)
         self.swqual = sup.lines[self.key]['swqual']
         if acn.division < 3:
          self.swqual = 0
         #print('nontrans div qual',acn.division, self.swqual)
         self.csw = sup.lines[self.key]['csw']
         self.status = theirnontrans.st
         #self.mag_max = theirnontrans.line_mag_max

         self.r = r = theirnontrans.r
         self.rbase = r/sbase
         self.x = x = theirnontrans.x
         self.b = theirnontrans.b
         self.ge = r/(r**2 + x**2)
         self.be = -x/(r**2 + x**2)
         self.z = z = r + x*1j
         self.y = y = 1/z
         self.Gff = self.ge
         self.Gtt = self.ge

         self.Gft = -self.ge
         self.Gtf = -self.ge

         self.Bft = -self.be
         self.Btf = -self.be

         self.Bff = (self.be + self.b/2.0)
         self.Btt = (self.be + self.b/2.0)

         self.rating_a = theirnontrans.ratea/sbase
         self.rating_c = theirnontrans.ratec/sbase

     def shownontsb(self, log):
         log.joint(' nontsb count1 %d from i %d to i %d (count1 %d %d) swqual %d status %d rating_a %g\n' %(self.count1, self.i, self.j, self.acn.busnumbertocount[self.i], self.acn.busnumbertocount[self.j], self.swqual, self.status, self.rating_a))

class ACNnontransATbus:
     def __init__(self, log, ournontrans, original_or_reversed, ntcount0):
      self.ournon = ournontrans
      self.original_or_reversed = original_or_reversed #1 if original
      self.count1 = ntcount0 + 1
      self.parentcount1 = ournontrans.count1
      self.swqual = ournontrans.swqual

      self.status = ournontrans.status
      self.parent = ournontrans
      self.partner = None

      if self.original_or_reversed == 1:
       self.Gff = ournontrans.Gff
       self.Gft = ournontrans.Gft
       self.Bff = ournontrans.Bff
       self.Bft = ournontrans.Bft
       self.ourfrombus = ournontrans.ourfrombus
       self.ourtobus = ournontrans.ourtobus
       self.busi0 = ournontrans.ourfrombus.count - 1
       self.busj0 = ournontrans.ourtobus.count - 1
       #log.joint(" added original ntsb count %i with busi0 %i busj0 %i\n" % (self.count1, self.busi0, self.busj0))
      else:
       self.Gff = ournontrans.Gtt
       self.Gft = ournontrans.Gtf
       self.Bff = ournontrans.Btt
       self.Bft = ournontrans.Btf
       self.ourfrombus = ournontrans.ourtobus
       self.ourtobus = ournontrans.ourfrombus
       self.busi0 = ournontrans.ourtobus.count - 1
       self.busj0 = ournontrans.ourfrombus.count - 1

class ACNtrans:
 def __init__(self, log, acn, sup, theirtrans, ourfrombus, ourtobus, sbase, count1):
  self.acn = acn
  self.theirtrans = theirtrans
  self.count1 = count1
  self.ourfrombus = ourfrombus
  self.ourtobus = ourtobus
  self.ckt = theirtrans.ckt
  self.i = theirtrans.i
  self.j = theirtrans.j
  self.key = (self.i,self.j,self.ckt)
  self.swqual = sup.transformers[self.key]['swqual']
  if acn.division < 3:
   self.swqual = 0
  #print('trans div qual',acn.division, self.swqual)

  self.csw = sup.transformers[self.key]['csw']
  self.stat = theirtrans.stat
  self.rating_a = theirtrans.rata1/sbase
  self.rating_c = theirtrans.ratc1/sbase
  self.r12 = r12 = theirtrans.r12
  self.x12 = x12 = theirtrans.x12
  self.mag1 = theirtrans.mag1
  self.mag2 = theirtrans.mag2
  self.windv1 = theirtrans.windv1
  self.angle1 = theirtrans.ang1
  self.cod1 = theirtrans.cod1
  self.rma1 = theirtrans.rma1
  self.rmi1 = theirtrans.rmi1
  self.ntp1 = theirtrans.ntp1
  self.tab1 = theirtrans.tab1
  self.windv2 = theirtrans.windv2

  self.gmf = self.mag1
  self.bmf = self.mag2
  self.gf = r12/(r12**2.0 + x12**2.0)
  self.bf = -x12/(r12**2.0 + x12**2.0)
  self.tau0f = theirtrans.windv1/theirtrans.windv2
  pi = math.pi
  self.theta0f = theirtrans.ang1*pi/180.0
  if self.ntp1 > 1:
   self.barxstf = int(round(0.5 * (self.ntp1 - 1)))
  else:
   self.barxstf = 0
  self.maxxstf = self.barxstf
  self.minxstf = -self.barxstf  #defaults as per eq (61)
  self.Ftauvalid = 0
  if self.cod1 == 1 or self.cod1 == -1:
   self.overtauf = self.rma1
   self.undertauf = self.rmi1
   self.Ftauvalid = 1
  else:
   self.undertauf = self.overtauf = self.tau0f
  self.acn.numtautransformers += self.Ftauvalid
  '''
  print('>>> i', self.i,' j', self.j, 'mag1', self.mag1, 'stat',self.stat,'r12',self.r12)
  print(' count',self.count,'ntp1',self.ntp1,'tab1',self.tab1,self.barxstf)

  breakexit('')
  '''

  self.taustf = 0
  self.brevetauf = (self.overtauf + self.undertauf)/2
  if self.ntp1 > 1:
   self.taustf = (self.overtauf - self.undertauf)/(2*self.barxstf)

  self.Fthetavalid = 0
  if self.cod1 == 3 or self.cod1 == -3:
   self.overthetaf = self.rma1*pi/180.0
   self.underthetaf = self.rmi1*pi/180.0
   self.Fthetavalid = 1
  else:
   self.underthetaf = self.overthetaf = self.theta0f
  self.thetastf = 0
  if self.ntp1 > 1:
   self.thetastf = (self.overthetaf - self.underthetaf)/(2*self.barxstf)
  self.brevethetaf = (self.overthetaf + self.underthetaf)/2
  self.barsf = self.rating_a/sbase
  self.barsctf = self.rating_c/sbase
  self.acn.numthetatransformers += self.Fthetavalid
  '''
  log.joint('  trans %d i %d j %d ckt %s (cnt1 %d cnt1 %d) stat %d r12 %g x12 %g g %g b %g cod1 %d tab1 %d ntp1 %d taustf %g thetastf %g,\n             tauvalid %d thetavalid %d\n' %(self.count1, self.i, self.j, self.ckt, self.ourfrombus.count, self.ourtobus.count, self.stat, self.r12, self.x12, self.gf, self.bf, self.cod1, self.tab1, self.ntp1, self.taustf, self.thetastf, self.Ftauvalid, self.Fthetavalid))
  log.joint('    brevetau %g tau0f %g\n    brevetheta %g theta0f %g\n' %(self.brevetauf, self.tau0f, self.brevethetaf, self.theta0f))

  '''
  self.variableratio = 0
  if self.cod1 == -1 or self.cod1 == 1:
   self.variableratio = 1
   #log.joint('   variable ratio\n')
  self.variablephase = 0
  if self.cod1 == -3 or self.cod1 == 3:
   self.variablephase = 1
   #log.joint('   variable phase\n')

  self.impedancecorrected = 0
  self.numm = 0
  self.Fnuvalid = 0

  self.impedancecorr = {}
  self.listnux = {}

  if self.tab1 != 0 and (self.cod1 == -3 or self.cod1 == -1 or self.cod1 == 1 or self.cod1 == 3):
   log.joint('   impedancecorrected\n')
   self.Fnuvalid = 1

   self.impedancecorrected = 1
   self.acn.numimpedancecorrectedtransformers += 1

   for tict in acn.ticts.values():
    if tict.i == self.tab1:
     log.joint('   trans %d stat %d g %g b %g cod1 %d tab1 %d ntp1 %d taustf %g thetastf %g,\n             tauvalid %d thetavalid %d\n' %(self.count1, self.stat, self.gf, self.bf, self.cod1, self.tab1, self.ntp1, self.taustf, self.thetastf, self.Ftauvalid, self.Fthetavalid))

     if self.variableratio:
      log.joint('   variable ratio\n')
     if self.variablephase:
      log.joint('   variable phase\n')

     self.acn.numvariabletransformers += 1


     log.joint('   tict# %d i %d point_count %d\n' %(tict.count, tict.i, tict.point_count))
     log.joint('               holy shit i ' + str(self.i) + ' j ' + str(self.j) + ' ckt ' + str(self.ckt) + '\n')
     self.numm = tict.point_count
     self.Mf = [j for j in range(1,self.numm+1)]
     self.nuf = numpy.zeros(1+self.numm+1)



     for m in range(1,self.numm+1):
      self.nuf[m] = tict.corrtable[2*(m-1)+1]
      log.joint('    nuf[' + str(m) + '] = ' + str(self.nuf[m]) + '\n')

     self.tauf = {}
     self.thetaf = {}


     if self.Ftauvalid: #self.cod1 == -1 or self.cod1 == 1:
      denominator = self.taustf
      base = self.brevetauf
      tict_base_smallest = 100000
      tict_base_largest = -100000

      log.joint('     brevetau %.9e taustf %.9e\n' %(self.brevetauf, self.taustf))
      for m in range(1,self.numm+1):
       self.tauf[m] = tict.corrtable[2*(m-1)]
       log.joint('    tauf[' + str(m) + '] = '+ str(self.tauf[m]) + '\n')
       if self.tauf[m] > tict_base_largest:
        tict_base_largest = self.tauf[m]
       if self.tauf[m] < tict_base_smallest:
         tict_base_smallest = self.tauf[m]

     if self.Fthetavalid: #self.cod1 == -3 or self.cod1 == 3:
      denominator = self.thetastf
      base = self.brevethetaf
      log.joint('     theta0 %g brevetheta %.9e thetastf %.9e\n' %(self.theta0f, self.brevethetaf, self.thetastf))
      log.joint('     compatible x: ' + str((self.theta0f - self.brevethetaf)/self.thetastf ) + '\n')
      tict_base_smallest = 100000
      tict_base_largest = -100000
      for m in range(1,self.numm+1):
       self.thetaf[m] = tict.corrtable[2*(m-1)]*pi/180.0
       if self.thetaf[m] > tict_base_largest:
        tict_base_largest = self.thetaf[m]
       if self.thetaf[m] < tict_base_smallest:
         tict_base_smallest = self.thetaf[m]

       log.joint('    thetaf[' + str(m) + '] = ' + str(self.thetaf[m]) + '\n')


     #print('corrtable',tict.corrtable)
     log.joint('     updating min, max xstf values\n')
     log.joint('     tict_base_smallest %g tict_base_largest %g denominator %g\n' %(tict_base_smallest, tict_base_largest, denominator))
     log.joint('     pre: minx %d maxx %d\n' %(self.minxstf, self.maxxstf))


     implied_largest = int(round( (tict_base_largest - base)/denominator ))
     if implied_largest < self.maxxstf:
      self.maxxstf = int(round(implied_largest))

     log.joint('     tict_base_largest %g implied_largest %d so max %d\n' %(tict_base_largest,implied_largest, self.maxxstf))

     implied_smallest = int(round( (tict_base_smallest - base)/denominator ))
     if implied_smallest > self.minxstf:
      self.minxstf = int(round(implied_smallest))

     log.joint('     tict_base_smallest %g implied_smallest %d so min %d\n' %(tict_base_smallest,implied_smallest, self.minxstf))

     #breakexit('ticted')
     log.joint('   now constructing evaluation table\n')
     segment = 1
     xvalue = math.floor(self.minxstf)
     self.impedancecorr = {}
     countnu = 1

     if self.Ftauvalid:
      while xvalue < self.maxxstf + 1:
       tauvalue = self.brevetauf + self.taustf*xvalue
       #log.joint('   doing xvalue %d tau %g prev nu %g next nu %g\n' %(xvalue, tauvalue, self.nuf[segment], self.nuf[segment + 1]))
       while tauvalue > self.tauf[segment + 1]:
        segment += 1
        log.joint('    > moved to segment %d\n' %(segment))
       base = tauvalue - self.tauf[segment]
       slope = (self.nuf[segment+1] - self.nuf[segment])/(self.tauf[segment+1] - self.tauf[segment])
       computednu = slope*base + self.nuf[segment]
       #log.joint('     slope %g base %g -> computed %g\n' %(slope, base, computednu))
       self.impedancecorr[countnu] = computednu
       self.listnux[countnu] = xvalue
       countnu += 1
       xvalue += 1


     if self.Fthetavalid:
      log.joint('  beginning computednu loop with xvalue = ' + str(xvalue))
      while xvalue < self.maxxstf + 1:
       thetavalue = self.brevethetaf + self.thetastf*xvalue
       #log.joint('   doing xvalue %d theta %g prevtheta %g nexttheta %g prevnu %g nextnu %g\n' %(xvalue, thetavalue, self.thetaf[segment], self.thetaf[segment + 1], self.nuf[segment], self.nuf[segment + 1]))

       while thetavalue > self.thetaf[segment + 1]:
        segment += 1
        #log.joint('    > moved to segment %d\n' %(segment))

       base = thetavalue - self.thetaf[segment]
       slope = (self.nuf[segment+1] - self.nuf[segment])/(self.thetaf[segment+1] - self.thetaf[segment])
       computednu = slope*base + self.nuf[segment]
       #log.joint('     >slope %g base %g -> computed %g\n' %(slope, base, computednu))
       self.impedancecorr[countnu] = computednu
       self.listnux[countnu] = xvalue
       countnu += 1
       xvalue += 1




   #print(self.impedancecorr)
   #breakexit('')
   #log.joint('       impedance correction has length %d\n' %(len(self.impedancecorr)))
   #print(' - > ', self.maxxstf)
   #breakexit('hoo')
  else:
   self.listnux[1] = 0
   #self.maxxstf = self.minxstf = 0
   self.impedancecorr[1] = 1
   #log.joint('     not impedancecorrected\n')

  if len(self.impedancecorr) == 0:
   log.joint(' no impedance correction\n')
   self.impedancecorr[1] = 1.0
   self.listnux[1] = 0

  #now get the true g values
  self.trueg = {}
  self.trueb = {}
  for countnu in range(1,1+len(self.listnux)):
   self.trueg[countnu] = self.gf/self.impedancecorr[countnu]
   self.trueb[countnu] = self.bf/self.impedancecorr[countnu]

  self.numnux = len(self.listnux)


  #log.joint('      numnux %d minxstf %d maxxstf %d -> num xstf %d\n' %(self.numnux,self.minxstf, self.maxxstf, 1 + self.maxxstf - self.minxstf))


  #if self.Ftauvalid and self.numnux != 1 + self.maxxstf - self.minxstf:
  # log.joint('      numnux %d minxstf %d maxxstf %d -> num xstf %d\n' %(self.numnux,self.minxstf, self.maxxstf, 1 + self.maxxstf - self.minxstf))
  # breakexit(" ALERT numnux (tau) inconsistent")


  #if self.Fthetavalid and self.numnux != 1 + self.maxxstf - self.minxstf:
  # log.joint('      numnux %d minxstf %d maxxstf %d -> num xstf %d\n' %(self.numnux,self.minxstf, self.maxxstf, 1 + self.maxxstf - self.minxstf))
  # breakexit(" ALERT numnux (theta) inconsistent")

  #breakexit('haa')
  #the rest is carried over but we may toss it

  ratio = theirtrans.windv1/theirtrans.windv2
  self.invratio2 = invratio2 = 1/ratio**2
  self.multtf = multtf = 1/(ratio*cmath.exp(1j*self.theta0f))
  self.multft = multft = 1/(ratio*cmath.exp(-1j*self.theta0f))
  self.z = z = r12 + x12*1j
  self.y = y = 1/z
  self.Yff = Yff = y*invratio2 + self.mag1 + self.mag2*1j
  self.Yft = Yft = -y*multft
  self.Ytf = Ytf = -y*multtf
  self.Ytt = Ytt = y
  #bypass the complex arithmetic
  self.Gff = self.gf/(ratio**2) + self.mag1
  self.Gft = -self.gf/ratio
  self.Bft = -self.bf/ratio
  self.Gtt = self.gf
  self.Gtf = -self.gf/ratio
  self.Btf = -self.bf/ratio
  self.Bff = (self.bf/(ratio**2) + self.mag2)
  self.Btt = self.bf
 def showtrans_slim(self, log):
  log.joint(" trans %i " %self.count)
  self.ourfrombus.showbus_slim_nocr(log)
  self.ourtobus.showbus_slim_nocr(log)
  log.joint(" r12 = %.8e x12 = %.8e gf = %.8e bf = %.8e mag1 = %.8e mag2 = %.8e tauf = %.8e\n theta0f_rad = %.16e\n" %(self.r12, self.x12, self.gf, self.bf, self.mag1, self.mag2, self.tauf, self.theta0f_rad))
  log.joint("y " + str(self.y) + "\n")
  log.joint("Yff " + str(self.Yff) + " , Yft " + str(self.Yft) + "\n")
  log.joint("Ytf " + str(self.Ytf) + " , Ytt " + str(self.Ytt) + "\n")

  log.joint("            Gff = %.16e Gft = %.16e Bft = %.16e\n" %(self.Gff, self.Gft, self.Bft))
  log.joint("            Gtt = %.16e Gtf = %.16e Btf = %.16e\n" %(self.Gtt, self.Gtf, self.Btf))
  log.joint("            Btt = %.16e Bff = %.16e\n" %(self.Btt, self.Bff))
  log.joint("   v1 %f v2 %f ratio %f\n" %(self.theirtrans.windv1, self.theirtrans.windv2, self.tauf))

class ACNtransATbus:
     def __init__(self, log, ourtrans, original_or_reversed,tcount0):
      self.ournon = ourtrans
      self.original_or_reversed = original_or_reversed #1 if original
      self.count1 = tcount0 + 1
      self.parentcount1 = ourtrans.count1
      self.swqual = ourtrans.swqual
      self.stat = ourtrans.stat
      self.parent = ourtrans
      self.partner = None

      if self.original_or_reversed == 1:
       self.Gff = ourtrans.Gff
       self.Gft = ourtrans.Gft
       self.Bff = ourtrans.Bff
       self.Bft = ourtrans.Bft
       self.ourfrombus = ourtrans.ourfrombus
       self.ourtobus = ourtrans.ourtobus
       self.busi0 = ourtrans.ourfrombus.count - 1
       self.busj0 = ourtrans.ourtobus.count - 1
       #log.joint("   added original tsb count %i with busi0 %i busj0 %i parent1 %d\n" % (self.count1, self.busi0, self.busj0, self.parentcount1))
      else:
       self.Gff = ourtrans.Gtt
       self.Gft = ourtrans.Gtf
       self.Bff = ourtrans.Btt
       self.Bft = ourtrans.Btf
       self.ourfrombus = ourtrans.ourtobus
       self.ourtobus = ourtrans.ourfrombus
       self.busi0 = ourtrans.ourtobus.count - 1
       self.busj0 = ourtrans.ourfrombus.count - 1
       #log.joint("   added reverse tsb count %i with busi0 %i busj0 %i parent1 %d\n" % (self.count1, self.busi0, self.busj0, self.parentcount1))


class ACNtict:
 def __init__(self, log, r, count):
  self.i = r.i
  self.count = count
  self.their_r = r
  self.point_count = r.tict_point_count
  self.corrtable = [r.t1, r.f1, r.t2, r.f2, r.t3, r.f3, r.t4, r.f4, r.t5, r.f5,
             r.t6, r.f6, r.t7, r.f7, r.t8, r.f8, r.t9, r.f9, r.t10, r.f10,
             r.t11, r.f11][0:(2*r.tict_point_count)]
  #to do: create a local data structure that has suggestive names

class ACNpriorsol:
 def __init__(self, name):
  self.name = name
  self.xst0fdict = {}
