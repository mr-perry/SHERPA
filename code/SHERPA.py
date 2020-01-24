from SHERPA_funcs import *
from datetime import datetime

def main():
  """
  SHARAD EDR Processing Algorithm: SHERPA

  Input:
    Full Path to SHARAD PDS Label File
      For the current version please ensure the SCIENCE and AUXILIARY data
      files are within the same directory as the LBL file

  Written By: Matthew R. Perry
  Last Updated: 20 January 2020
  """
  st = datetime.now()
  prog = "SHARAD EDR Processing Algorithm"
  vers = "0.1"
  iFiles, oFiles, TransID, OSTLine, OperMode, verb = parseargs(prog, vers)
  #
  # Parse Auxiliary File
  # 
  parseAuxFile(iFiles['AUX'], oFiles['AUX'], dic=False, csv=True)
  cnt = sepSAdata(iFiles['SCIENCE'], oFiles['SCIENCE'], oFiles['ANCILLARY'], OperMode['nrec'], oFiles['_log'])
  tt = datetime.now() - st
  writeLog(oFiles['_log'], 'Total time:\t{}'.format(tt), verb=verb)
  oFiles['_log'].close()
  return

if __name__ == '__main__':
  main()
