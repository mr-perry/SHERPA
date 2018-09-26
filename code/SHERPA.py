#
# Import necessary libraries
#
from SHERPA_func import *

def main():
  prog = 'processEDR'
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
  # Construct window
  #
  if chirp == 'ideal' or chirp == 'UPB':
    window, win_str = makeWindow(beta, length=3600)
  else:
    window, win_str = makeWindow(beta, length=2048)
  writeLog(_log, win_str)
  ######################################################################
  #
  # Begin Processing
  #
  ######################################################################
  if chirp == 'ideal' or chirp == 'UPB':
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
  for _i in range(0, nrec, presum_fac):
    #
    # Get the time since the start of the run for this pulse group
    #
    tcen = AuxDF['ELAPSED_TIME'][_i]
    for _k in range(0, presum_fac):
      #
      # Read in single record
      #
      sci, ancil = readEDRrecord(_file1, _i, recLen, BitsPerSample)
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
      calChirp = detChirpFiles(AuxDF['TX_TEMP'][_i], AuxDF['RX_TEMP'][_i], chirp=chirp)
      #
      # Perform chirp compression
      #
      presum_rec[:,_k] = rangeCompression(sci, calChirp, window, chirp=chirp, fil_type="Match", diag=diag)
      #
      # Perform on-ground calibration
      #  This step will have to wait. Apparently the angles given in the Auxilliary file are
      #  incorrect for some dates and change after a certain date. Fabrizio is checking on this
      #
#      presum_rec[:, _k] = calibrateData(presum_rec[:, _k])
    EDRData[:,int(_i/presum_fac)] = np.sum(presum_rec, axis=1)
  if verb:
    writeLog(_log, 'Decompression finished at:\t{}'.format(datetime.now()))
  fname = '../output/' + str(runName) + '.npy'
  np.save(fname, EDRData)
  return

if __name__ == '__main__':
  main()
