import glob, sys, os, struct, bitstring, argparse
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


#########################################################################
#
# Functions
#
#########################################################################
#
# Log Functions
#
###############################################
def writeLog(log, string):
  print(string)
  if log is not None:
    log.write(string + '\n')
###############################################
#
# General Functions
#
###############################################
def parseargs(prog, vers):
  #
  # Set the default values
  #
  runName = ['processEDR']
  chirp = ['ref']
  beta = [0]
  verb = [False]
  diag = [False]
  fil_type = ['match']
  presum_proc = [None]
  path = './input/'
  #
  # Initiate the parser
  #
  parser = argparse.ArgumentParser(description=str(prog + ' ' + vers))
  #
  # Ordered Arguments
  #
  parser.add_argument('ObsID', type=str, nargs=1,
                      help=str('7 character observation ID from PDS e.g., 0402801'))
  #
  # Optional Arguments
  #
  parser.add_argument('-a', '--path', nargs=1, type=str, default=path,
                     help=str('Path to data directory where the label, auxiliary, and science dataare found'))
  parser.add_argument('-r', '--runname', nargs=1, default=runName, type=str,
                     help=str('Desired output name'))
  parser.add_argument('-c', '--chirp', nargs=1, default=chirp, type=str,
                     help=str('[UPB, Ideal, ref, Vibro]'))
  parser.add_argument('-f', '--filtype', nargs=1, default=fil_type, type=str,
                     help=str('Match or Inverse'))
  parser.add_argument('-p', '--presumproc', nargs=1, default=presum_proc, type=int,
                     help=str('Processing Presum Value'))
  parser.add_argument('-b', '--beta', nargs=1, default=beta, type=float,
                     help=str('Kaiser window smoothing factor: 0 - uniform'))
  parser.add_argument('-v', '--verbose', action="store_true", default=verb,
                     help=str('Print to screen as well as log file'))
  parser.add_argument('-d', '--diagnostic', action="store_true", default=diag,
                     help=str('Show diagnostic plots. Will exit after 1 presum window'))
  #
  # Print help
  #
  if len(sys.argv[1:]) == 0:
    parser.print_help()
    parser.exit()
  #
  # Parse arguments
  #
  args = parser.parse_args()
  ObsID = args.ObsID[0]
  runName = args.runname[0]
  path = args.path[0]
  if args.chirp[0].lower() not in ['upb', 'ideal', 'ref', 'vibro']:
    print('Selection of chirp type not understood.')
    print('You selected: {}'.format(args.chirp[0]))
    print('Allowed values: [upb, ideal, ref, vibro]')
    print('Bye-bye for now...')
    parser.print_help()
    parser.exit()
  else:
    chirp = args.chirp[0].lower()
  if args.filtype[0] not in ['match', 'inverse']:
    print('Selection of filter type not understood.')
    print('You selected: {}'.format(args.filtype[0]))
    print('Allowed values: [match, inverse]')
    print('Bye-bye for now...')
    parser.print_help()
    parser.exit()
  else:
    fil_type = args.filtype[0].lower()
  beta = args.beta[0]
  try:
    verb = args.verbose[0]
  except:
    verb = args.verbose
  if args.diagnostic[0] == True:
    print("You've selected to produce diagnostic plots. Please be advised that after 1"\
          "presum processing window, the processor will exit and your plots will be saved")
  try:
    diag = args.diagnostic[0]
  except:
    diag = args.diagnostic
  presum_proc = args.presumproc[0]
  return ObsID, path, runName, chirp, fil_type, beta, presum_proc, verb, diag


def getData(path, ObsID):
  logs = []
  lblFile = []
  auxFile = []
  sciFile = []
  if path[0:4].lower() == 'http':
    logs.append(_log, 'This feature is not available at this time')
    logs.append(_log, 'Bye-bye for now...')
    return lblFile, auxFile, sciFile, logs
  else:
    path = path + '/' if path[-1] != '/' else path
    fileList = glob.glob(path + 'e_' + ObsID + '*')
    if len(fileList) > 3:
      writeLog(_log, 'Too many files found associated with ObsID {}'.format(ObsID))
      for _i in fileList:
        logs.append(str(fileList[_i]))
        logs.append(str('Bye-bye for now...'))
        return lblFile, auxFile, sciFile, logs
    else:
      logs.append('{} files found for ObsID {}'.format(str(len(fileList)), ObsID))
      try:
        lblFile = glob.glob(path + 'e_' + ObsID + '*.lbl')[0]
        auxFile = glob.glob(path + 'e_' + ObsID + '*a.dat')[0]
        sciFile = glob.glob(path + 'e_' + ObsID + '*s.dat')[0]
        logs.append('PDS Label File:\t{}'.format(lblFile))
        logs.append('Auxiliary Data File:\t{}'.format(auxFile))
        logs.append('Science Data File:\t{}'.format(sciFile))
        return lblFile, auxFile, sciFile, logs
      except:
        logs.append('Files associated with ObsID {} not found. Please check your path and try again.'.format(ObsID))
        return lblFile, auxFile, sciFile, logs



def main():
  prog = 'SHARAD EDR PROCESSING ALGORITHM (SHERPA)'
  vers = '0.1'
  ObsID, path, runName, chirp, fil_type, beta, presum_proc, verb, diag = parseargs(prog, vers)
  #
  # Open log file
  #
  logFile = '../output/' + str(runName) + '.log'
  _log = open(logFile, 'w')
  writeLog(_log, '--- {} {} ---'.format(prog, vers))
  #
  # Print summary
  #
  writeLog(_log, 'PDS Observation ID:\t{}'.format(ObsID))
  writeLog(_log, 'Run Name:\t{}'.format(runName))
  writeLog(_log, 'Type of chirp for range compression:\t{}'.format(chirp))
  writeLog(_log, 'Type of filtering for range compression:\t{}'.format(fil_type))
  writeLog(_log, 'Kaiser Window Smoothing Coefficient:\t{}'.format(str(beta)))
  writeLog(_log, 'Verbose Output:\t{}'.format(verb))
  writeLog(_log, 'Diagnostic Plots:\t{}'.format(diag))
  #
  # Find required data files
  #
  lblFile, auxFile, sciFile, logs = getData(path, ObsID)
  for log in logs:
    writeLog(_log, log)
  #
  # Parse file name for information
  #
  ObsID, OSTLine, OperMode, PRF, Version = parseFileName(sciFile)
  InstrPresum = OperMode['Presum']
  if presum_proc is None:
    presum_proc = OperMode['Presum']
  presum_fac = int(presum_proc / InstrPresum)
  if presum_fac < 1:
    presum_fac = 1
    print('WARNING: Processing presum less than onboard presumming....keeping native presum value')
  instrMode = OperMode['Mode']
  BitsPerSample = OperMode['BitsPerSample']
  #
  # Determine to Bits per Sample
  #
  if BitsPerSample == 4:
    recLen = 1986
  elif BitsPerSample == 6:
    recLen = 2886
  elif BitsPerSample == 8:
    recLen = 3786
  #
  # Get the number of records from the file size of the science data
  # divided by the record length
  #
  nrec = int(os.path.getsize(sciFile) / recLen)
  writeLog(_log, 'Reading Auxilliary file...')
  AuxDF = parseAuxFile(auxFile, df=True)
  writeLog(_log, 'Finished reading Auxilliary file...')
  if verb:
    writeLog(_log, 'Number of Records:\t{}'.format(nrec))
    writeLog(_log, '----- Presum Information -----')
    writeLog(_log, 'InstrPresum:\t{}'.format(InstrPresum))
    writeLog(_log, 'Processing Presum:\t{}'.format(presum_proc))
    if presum_fac == 1:
      writeLog(_log, 'No additional presumming performed')
    else:
      writeLog(_log, 'Presum Factor:\t{}'.format(presum_fac))
      writeLog(_log, '----- End Presum Information -----')
    writeLog(_log, 'Instrument Mode:\t{}'.format(instrMode))
    writeLog(_log, 'Bits Per Sample:\t{}'.format(BitsPerSample))
    writeLog(_log, 'Record Length:\t{}'.format(recLen))
    writeLog(_log, 'Longitude Range of Observation:\t{}, {}'.format(np.amin(AuxDF['SUB_SC_EAST_LONGITUDE']), np.amax(AuxDF['SUB_SC_EAST_LONGITUDE'])))
    writeLog(_log, 'Planetocentric Latitude Range of Observation:\t{}, {}'.format(np.amin(AuxDF['SUB_SC_PLANETOCENTRIC_LATITUDE']), np.amax(AuxDF['SUB_SC_PLANETOCENTRIC_LATITUDE'])))
    writeLog(_log, 'Planetographic Latitude Range of Observation:\t{}, {}'.format(np.amin(AuxDF['SUB_SC_PLANETOGRAPHIC_LATITUDE']), np.amax(AuxDF['SUB_SC_PLANETOGRAPHIC_LATITUDE'])))
    writeLog(_log, 'Solar Longitude Range of Observation:\t{}, {}'.format(np.amin(AuxDF['SOLAR_LONGITUDE']), np.amax(AuxDF['SOLAR_LONGITUDE'])))
  #
  # Construct window; this is windowing the time series raw data
  #
  if chirp == 'ideal' or chirp == 'upb' or chirp == 'UPB':
    window, win_str = makeWindow(beta, length=3600)
  else:
    if chirp == 'ref' or chirp == 'vibro':
      window, win_str = makeWindow(beta, length=2048)
  ######################################################################
  #
  # Begin Processing
  #
  ######################################################################
  if chirp == 'ideal' or chirp == 'upb':
    if chirp == 'ideal':
      print('Using ideal chirp')
    else:
      print('Using cal_filter chirp from UPB')
    EDRData = np.zeros([3600, int(np.ceil(nrec/presum_fac))], complex)
    presum_rec = np.zeros([3600, presum_fac], complex)
  else:
    EDRData = np.zeros([4096, int(np.ceil(nrec/presum_fac))], complex)
    presum_rec = np.zeros([4096, presum_fac], complex)
  writeLog(_log, 'Opening EDR File:\t{}'.format(sciFile))
  _file1 = open(sciFile, 'rb')
  if verb:
    writeLog(_log, 'Begin decompression at:\t{}'.format(datetime.now()))
  pl = -1
  plot_data = []
  for _i in range(0, nrec, presum_fac):
    #
    # Get the time since the start of the run for this pulse group
    #
    tcen = AuxDF['ELAPSED_TIME'][_i]
    for _k in range(0, presum_fac):
      #
      # Read in single record
      #
      _frame = _i + _k
      sci, ancil = readEDRrecord(_file1, _frame, recLen, BitsPerSample)
      #
      # Parse Ancilliary Datq
      #
      ancil = parseAncilliary(ancil)
      #
      # Decompress science data
      #
      sci = decompressSciData(sci,
              ancil['OST_LINE']['COMPRESSION_SELECTION'],
              InstrPresum,
              BitsPerSample,
              ancil['SDI_BIT_FIELD'])
      #
      # Determine calibrated chirp
      #
      calChirp = detChirpFiles(AuxDF['TX_TEMP'][_frame], AuxDF['RX_TEMP'][_frame], chirp=chirp)
      #
      # Perform chirp compression
      #
      presum_rec[:,_k] = rangeCompression(sci, calChirp, window, chirp=chirp, fil_type="Match", diag=diag)
#      amp = presum_rec[:,_k]
#      pwr = np.real(amp)**2
#      dB = 10 * np.log10(pwr / np.max(pwr))
#      plt.plot(dB)
#      plt.show()
#      exit()
      #
      # Perform on-ground calibration
      #  This step will have to wait. Apparently the angles given in the Auxilliary file are
      #  incorrect for some dates and change after a certain date. Fabrizio is checking on this
      #
#      presum_rec[:, _k] = calibrateData(presum_rec[:, _k])
    EDRData[:,int(_i/presum_fac)] = np.sum(presum_rec, axis=1)
    amp = EDRData[:,int(_i/presum_fac)]
    np.save(runName+'.npy', amp)
    pwr = np.abs(amp)**2
    pwr = pwr / np.max(pwr)
    dB = 10 * np.log10(pwr)
    plt.plot(dB)
    plt.savefig(runName+'.png')
    plt.show()
    exit()
  if verb:
    writeLog(_log, 'Decompression finished at:\t{}'.format(datetime.now()))
  fname = '../output/' + str(runName) + '.npy'
  np.save(fname, EDRData)
  return

if __name__ == '__main__':
  main()
