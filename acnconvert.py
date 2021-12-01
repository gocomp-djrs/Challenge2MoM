from myutils import *
from log import *
import sys, os
import numpy
from acn import *

def acnconvert(alldata):
    code = 0

    log = alldata['log']

    arpae = alldata['arpae']
    raw = arpae.data.raw


    log.joint('\nacnconvert: extracting data from arpae object\n')
    try:
      alldata['acn'] = ACN(alldata)
    except Exception as e:
      log.joint('Failed to set alldata: '+ str(e))


    log.joint('acnconvert returns with code ' + str(code) + '\n')
    return code
