from SHERPA_funcs import *
from datetime import datetime

def sherpa(lblFile, outDir, roi=[None,None,None,None]):
  verb=False
  #
  # Get start time
  #
  st = datetime.now()
  #
  # Program information
  #
  prog = "SHARAD EDR Processing Algorithm"
  vers = "0.1"
  #
  # Find science and auxiliary files
  #
  iFiles = findFiles(lblFile)
  #
  # Form outDir
  #
  oDirs, oFiles = formOut(lblFile, outDir)
  #
  # Make output directories
  #
  makeOut(outDir, oDirs) 
  #
  # Parse the PDS file name for important information
  #
  TransID, OSTLine, OperMode = parseFileName(iFiles['SCIENCE'])
  #
  # Open Log File
  #
  oFiles['_log'] = open(oFiles['LOG'], 'w')
  writeLog(oFiles['_log'], '----- {} {} Processing Log -----'.format(prog, vers), verb=verb)
  writeLog(oFiles['_log'], 'PDS Label file:\t{}'.format(iFiles['LABEL']), verb=verb)
  writeLog(oFiles['_log'], 'Auxiliary file:\t{}'.format(iFiles['AUX']), verb=verb)
  writeLog(oFiles['_log'], 'Science file:\t{}'.format(iFiles['SCIENCE']), verb=verb)
  writeLog(oFiles['_log'], 'Output directory:\t{}'.format(outDir), verb=verb)
  writeLog(oFiles['_log'], 'TransID:\t{}'.format(TransID), verb=verb)
  writeLog(oFiles['_log'], 'OSTLine:\t{}'.format(OSTLine), verb=verb)
  writeLog(oFiles['_log'], 'OperMode:\t{}'.format(OperMode['Mode']), verb=verb)
  writeLog(oFiles['_log'], '\tOn-board presumming:\t{}'.format(OperMode['Presum']), verb=verb)
  writeLog(oFiles['_log'], '\tBits per Sample:\t{}'.format(OperMode['BitsPerSample']), verb=verb)
  writeLog(oFiles['_log'], 'PRF:\t{}'.format(OperMode['PRF']), verb=verb)
  writeLog(oFiles['_log'], 'Record Length:\t{}'.format(OperMode['recLen']), verb=verb)
  writeLog(oFiles['_log'], 'Number of Records:\t{}'.format(OperMode['nrec']), verb=verb)
  writeLog(oFiles['_log'], '', verb=verb)
  ##########################################################
  #
  # Start of actual processing...
  #
  ##########################################################
  #
  # Parse Auxiliary File
  #
  _, idx = parseAuxFile(iFiles['AUX'], oFiles['AUX'], dic=False, csv=True, roi=roi) 
  #
  # Separate the SCIENCE and ANCILLARY data
  #
  cnt = sepSAdata(iFiles['SCIENCE'], oFiles['SCIENCE'], oFiles['ANCILLARY'], OperMode['recLen'], OperMode['BitsPerSample'], OperMode['Presum'])
  #
  # Get the total time (tt)
  #
  tt = datetime.now() - st
  #
  # Report time and close log file
  #
  writeLog(oFiles['_log'], 'Total time:\t{}'.format(tt), verb=verb)
  oFiles['_log'].close()
  return
  
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
  #
  # Get start time
  #
  st = datetime.now()
  #
  # Program information
  #
  prog = "SHARAD EDR Processing Algorithm"
  vers = "0.1"
  #
  # Parse the arguments
  #
  iFiles, oFiles, TransID, OSTLine, OperMode, roi, verb = parseargs(prog, vers)
  #
  # Parse Auxiliary File
  # 
  _, idx = parseAuxFile(iFiles['AUX'], oFiles['AUX'], dic=False, csv=True, roi=roi)
  #
  # Separate the SCIENCE and ANCILLARY data
  #
  cnt = sepSAdata(iFiles['SCIENCE'], oFiles['SCIENCE'], oFiles['ANCILLARY'], OperMode['recLen'], OperMode['BitsPerSample'], OperMode['Presum'], idx=idx)
  #
  # Get the total time (tt)
  #
  tt = datetime.now() - st
  #
  # Report time and close log file
  # 
  writeLog(oFiles['_log'], 'Total time:\t{}'.format(tt), verb=verb)
  oFiles['_log'].close()
  return

if __name__ == '__main__':
  main()
