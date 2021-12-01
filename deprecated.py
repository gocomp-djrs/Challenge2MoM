 bhaAMPL = ampl.getSet('bha')
 bhsVals = {}
 for h in range(1,1 + acn.numswitchedshunts):
  oursh = acn.ourswshunts[h]
  thisbh = []
  if oursh.status:
   for a in range(1,1+oursh.len):
    thisbh.append(oursh.b[a])
  bhsVals[h] = thisbh
  if oursh.status:
   print (h, bhsVals[h])
   breakexit('poo')
   bhaAMPL[h].setValues(bhsVals[h])

------------------
from acnsolverAMPL.py Jan 12 2021


 Gf0a = ampl.getVariable('Gf0')
 Gf0f = Gf0a.getValues()
 Gf0vals = Gf0f.toDict()
 
 gf0a = ampl.getVariable('gf0')
 gf0f = gf0a.getValues()
 gf0vals = gf0f.toDict()

 bf0a = ampl.getVariable('bf0')
 bf0f = bf0a.getValues()
 bf0vals = bf0f.toDict()

 phasea = ampl.getVariable('thetaf0')
 phasef = phasea.getValues()
 phasevals = phasef.toDict()

 gf0a = ampl.getVariable('gf0')
 gf0f = gf0a.getValues()
 gf0vals = gf0f.toDict()

 bina = ampl.getVariable('binxstf0')
 binf = bina.getValues()
 binvals = binf.toDict()

 # Transformer tap settings
 vx = ampl.getVariable('xstf0')
 fx = vx.getValues()
 xvals = fx.toDict()

 if mykeybus >= 1 and mykeybus <= acn.numbuses:
    thisbus = acn.ourbuses[mykeybus]
    for trATbus in thisbus.startrans.values():
        cnt1 = trATbus.count1
        #print(cnt1, acn.transcount)
        log.joint(' !!! ' + str(cnt1) + ' pf ' + str(basetrP[cnt1-1]) + ' parent ' + str(trATbus.parentcount1) + '\n')
        log.joint('    ' + ' busi ' + str(trbusiVals[cnt1]) + ' busj ' + str(trbusjVals[cnt1]) + ' tron ' + str(basetrON[trATbus.parentcount1-1])+ '\n')
        #breakexit('')
        if cnt1 == mykeytrdual:
            mykeytr1 = trATbus.parentcount1
            vGf0f = Gf0vals[cnt1]
            vgf0f = gf0vals[mykeytr1]
            vbf0f = bf0vals[mykeytr1]
            i = trbusiVals[cnt1]
            j = trbusjVals[cnt1]
            log.joint('   angle at ' + str(i) + ' ' + str(basebusVangles[i-1]) )
            log.joint(' angle at ' + str(j) + ' ' + str(basebusVangles[j-1]) )
            log.joint(' phase ' + str(phasevals[mykeytr1]) + '\n')
            log.joint('   tap ' + str(xvals[mykeytr1]) + '\n')
            minx = minxstfVals[mykeytr1]
            maxx = maxxstfVals[mykeytr1]
            log.joint('   min ' + str(minx) )
            log.joint('   max ' + str(maxx) + '\n')
            for m in range(acn.ourtrans[mykeytr1].minxstf, 1 + acn.ourtrans[mykeytr1].maxxstf):
                if binvals[mykeytr1,m] > 1e-5:
                    log.joint('   m  %d bin %g trueg %g trueb %g\n' %(m, binvals[mykeytr1,m], truegVals[mykeytr1,m], truebVals[mykeytr1,m]))
            log.joint('    G %g g %g b %g\n' %(vGf0f, vgf0f, vbf0f))
            breakexit('minmax')
-------------
