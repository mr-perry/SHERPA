import glob, sys, os, struct, bitstring, argparse
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

###############################################
#
# Log Functions
#
###############################################
def writeLog(log, string):
  print(string)
  if log is not None:
    log.write(string + '\n')

################################################
#
# Functions!!!
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

def makeWindow(beta, length=2048):
  """
  This function creates a windowing function using the 
  numpy kaiser function.
  
  For more information visit: 
    https://docs.scipy.org/doc/numpy-1.13.0/reference/generated/numpy.kaiser.html

  Input:
    Beta: Beta value (5 - Hamming, 6-Hanning
    Length: Length of the window [2048]

  Output:
    Window -- Vector of length 'length' 
    Win_str -- String describing the windowing function
  """
  #
  # Construct the window function
  #
  win_str ='Window function applied to data: '
  window = np.kaiser(length, beta)
  win_str += 'Kaiser Window - Beta {}'.format(str(beta))
  return window, win_str

def calcSNR(data):
  Amp = np.abs(np.real(data))
  signal = np.mean(np.amax(Amp, axis=0))
  noise = np.mean(Amp)
  SNR = 10*np.log10(np.power(signal/noise,2))
  return SNR

def rangeCompression(sci, calChirp, window, chirp='ref', fil_type='Match', diag=False):
  """
    This function performs the chirp compression on a single SHARAD record
    Inputs:
      sci - range (fast-time) record of decompressed science data
      calChirp - Cal. Chirp spectrum
      fil_type  - Match or Inverse
  """
  #
  # Get Chirp Length to allow use of the ideal chirp
  # 2048 -- reference chirp
  # 3600 -- ideal chirp
  #
  # Standard SHARAD information
  #
  dt = (3./80.)*1e-6			# 0.0375 Microseconds
  sharad_ipp = 1.0 / 700.28		# SHARAD ????
  Fc = 1 / dt				# Frequency sampling rate in Hz
  length = len(calChirp)		# Length of the calibrated chirp
  t = np.arange(0,length) * dt		# Time vector in microseconds
  if diag:
    plotTimeSeriesDiag()
  #
  # Check the length of the vectors and pad where necessary
  #
  if chirp == 'ideal' or chirp == 'UPB':
    tots_len = 3600
  elif chirp == 'ref' or chirp == 'vib':
    tots_len = 4096
  if len(sci) != tots_len:
    echoes = np.zeros(tots_len, complex)
    echoes[:len(sci)] = sci
  else:
    echoes = np.copy(sci)
  #
  # In some cases, it may be necessary to pad the chirp
  #
  if chirp == 'ref':
      #
      # Pad Chirp
      #
      dechirp = np.zeros(tots_len, complex)
      dechirp[0:2048] = np.conj(calChirp)
  else:
      dechirp = np.conj(calChirp)
  #
  # Now perform a Fourier transform to calculate the spectra
  #
  ecSpec = np.fft.fft(echoes) / len(echoes)
  ecFreq = np.fft.fftfreq(length, dt)
  #
  # Now, decompress the science data and return to time domain
  #
  if fil_type == 'Match':
#    decomp = np.fft.ifft(window*(dechirp*(ecSpec))) * len(sci)
    decomp = np.fft.ifft((dechirp*(ecSpec))) * len(ecSpec)
  elif fil_type == "Inverse":
    print('Filter type {} is currently unavailable.'.format(fil_type))
    print('Exiting...')
    sys.exit()
  if diag:
    plotRCDiagnostics()
  return decomp

##########################
#
# Ancilliary File Functions
#
##########################
def parseAncilliary(data):
  """
    This function parses the ancilliary data for an individual SHARAD EDR data record
  """
  #
  # Check that length of data is 186 bytes
  #
  if len(data) != 186:
    print('Incorrect data supplied. Please ensure data is a 186 byte string')
    return
  else:
    #
    # Set up dictionary for items
    #
    ancilliaryData = { 'SCET_BLOCK_WHOLE': bitstring.BitArray(data[0:4]).uint,
                       'SCET_BLOCK_FRAC': bitstring.BitArray(data[4:6]).uint,
                       'TLM_COUNTER': struct.unpack('>I', data[6:10])[0],
                       'FMT_LENGTH': struct.unpack('>H', data[10:12])[0],
                       'SPARE1': struct.unpack('>H', data[12:14])[0],
                       'SCET_OST_WHOLE': struct.unpack('>I', data[14:18])[0],
                       'SCET_OST_FRAC': struct.unpack('>H', data[18:20])[0],
                       'SPARE2': struct.unpack('>B', data[20:21])[0],
                       'OST_LINE_NUMBER': struct.unpack('>B', data[21:22])[0],
                       'OST_LINE': { },
                       'SPARE3': struct.unpack('>B', data[38:39])[0],
                       'DATA_BLOCK_ID': bitstring.BitArray(data[39:42]).uint,
                       'SCIENCE_DATA_SOURCE_COUNTER': struct.unpack('>H', data[42:44])[0],
                       'PACKET_SEGMENTATION_AND_FPGA_STATUS': { },
                       'SPARE4': struct.unpack('>B', data[46:47])[0],
                       'DATA_BLOCK_FIRST_PRI': bitstring.BitArray(data[47:50]).uint,
                       'TIME_DATA_BLOCK_WHOLE': struct.unpack('>I', data[50:54])[0],
                       'TIME_DATA_BLOCK_FRAC': struct.unpack('>H', data[54:56])[0],
                       'SDI_BIT_FIELD': struct.unpack('>H', data[56:58])[0],
                       'TIME_N': struct.unpack('>f', data[58:62])[0],
                       'RADIUS_N': struct.unpack('>f', data[62:66])[0],
                       'TANGENTIAL_VELOCITY_N': struct.unpack('>f', data[66:70])[0],
                       'RADIAL_VELOCITY_N': struct.unpack('>f', data[70:74])[0],
                       'TLP': struct.unpack('>f', data[74:78])[0],
                       'TIME_WPF': struct.unpack('>f', data[78:82])[0],
                       'DELTA_TIME': struct.unpack('>f', data[82:86])[0],
                       'TLP_INTERPOLATE': struct.unpack('>f', data[86:90])[0],
                       'RADIUS_INTERPOLATE': struct.unpack('>f', data[90:94])[0],
                       'TANGENTIAL_VELOCITY_INTERPOLATE': struct.unpack('>f', data[94:98])[0],
                       'RADIAL_VELOCITY_INTERPOLATE': struct.unpack('>f', data[98:102])[0],
                       'END_TLP': struct.unpack('>f', data[102:106])[0],
                       'S_COEFFS': np.array([struct.unpack('>f', data[106:110]),
                                         struct.unpack('>f', data[110:114]),
                                         struct.unpack('>f', data[114:118]),
                                         struct.unpack('>f', data[118:122]),
                                         struct.unpack('>f', data[122:126]),
                                         struct.unpack('>f', data[126:130]),
                                         struct.unpack('>f', data[130:134]),
                                         struct.unpack('>f', data[134:138])
                                        ]),
                       'C_COEFFS': np.array([struct.unpack('>f', data[138:142]),
                                         struct.unpack('>f', data[142:146]),
                                         struct.unpack('>f', data[146:150]),
                                         struct.unpack('>f', data[150:154]),
                                         struct.unpack('>f', data[154:158]),
                                         struct.unpack('>f', data[158:162]),
                                         struct.unpack('>f', data[162:166])
                                        ]),
                       'SLOPE': struct.unpack('>f', data[166:170])[0],
                       'TOPOGRAPHY': struct.unpack('>f', data[170:174])[0],
                       'PHASE_COMPENSATION_STEP': struct.unpack('>f', data[174:178])[0],
                       'RECEIVE_WINDOW_OPENING_TIME': struct.unpack('>f', data[178:182])[0],
                       'RECEIVE_WINDOW_POSITION': struct.unpack('>f', data[182:186])[0],
                     }
    #####################################################################################
    #
    # PACKET_SEGMENTATION_AND_FPGA_STATUS bit string
    #
    PSAFS = bitstring.BitArray(data[44:46])
    ancilliaryData['PACKET_SEGMENTATION_AND_FPGA_STATUS']['SCIENTIFIC_DATA_TYPE'] = PSAFS[0:1].uint
    ancilliaryData['PACKET_SEGMENTATION_AND_FPGA_STATUS']['SEGMENTATION_FLAG'] = PSAFS[1:3].uint
    ancilliaryData['PACKET_SEGMENTATION_AND_FPGA_STATUS']['SPARE1'] = PSAFS[3:8].uint
    ancilliaryData['PACKET_SEGMENTATION_AND_FPGA_STATUS']['SPARE2'] = PSAFS[8:12].uint
    ancilliaryData['PACKET_SEGMENTATION_AND_FPGA_STATUS']['DMA_ERROR'] = PSAFS[12:13].uint
    ancilliaryData['PACKET_SEGMENTATION_AND_FPGA_STATUS']['TC_OVERRUN'] = PSAFS[13:14].uint
    ancilliaryData['PACKET_SEGMENTATION_AND_FPGA_STATUS']['FIFO_FULL'] = PSAFS[14:15].uint
    ancilliaryData['PACKET_SEGMENTATION_AND_FPGA_STATUS']['TEST'] = PSAFS[15:16].uint
    #####################################################################################
    #
    # OST_LINE_NUMBER BIT STRING
    #
    OST_LINE = bitstring.BitArray(data[22:39])
    ancilliaryData['OST_LINE']['PULSE_REPETITION_INTERVAL'] = OST_LINE[0:4].uint
    ancilliaryData['OST_LINE']['PHASE_COMPENSATION_TYPE'] = OST_LINE[4:8].uint
    ancilliaryData['OST_LINE']['SPARE1'] = OST_LINE[8:10].uint
    ancilliaryData['OST_LINE']['DATA_LENGTH_TAKEN'] = OST_LINE[10:32].uint
    ancilliaryData['OST_LINE']['OPERATIVE_MODE'] = OST_LINE[32:40].uint
    ancilliaryData['OST_LINE']['MANUAL_GAIN_CONTROL'] = OST_LINE[40:48].uint
    ancilliaryData['OST_LINE']['COMPRESSION_SELECTION'] = OST_LINE[48:49].bool
    ancilliaryData['OST_LINE']['CLOSED_LOOP_TRACKING'] = OST_LINE[49:50].bool
    ancilliaryData['OST_LINE']['TRACKING_DATA_STORAGE'] = OST_LINE[50:51].bool
    ancilliaryData['OST_LINE']['TRACKING_PRE_SUMMING'] = OST_LINE[51:54].uint
    ancilliaryData['OST_LINE']['TRACKING_LOGIC_SELECTION'] = OST_LINE[54:55].uint
    ancilliaryData['OST_LINE']['THRESHOLD_LOGIC_SELECTION'] = OST_LINE[55:56].uint
    ancilliaryData['OST_LINE']['SAMPLE_NUMBER'] = OST_LINE[56:60].uint
    ancilliaryData['OST_LINE']['SPARE2'] = OST_LINE[60:61].uint
    ancilliaryData['OST_LINE']['ALPHA_BETA'] = OST_LINE[61:63].uint
    ancilliaryData['OST_LINE']['REFERENCE_BIT'] = OST_LINE[63:64].uint
    ancilliaryData['OST_LINE']['THRESHOLD'] = OST_LINE[64:72].uint
    ancilliaryData['OST_LINE']['THRESHOLD_INCREMENT'] = OST_LINE[72:80].uint
    ancilliaryData['OST_LINE']['SPARE3'] = OST_LINE[80:84].uint
    ancilliaryData['OST_LINE']['INITIAL_ECHO_VALUE'] = OST_LINE[84:87].uint
    ancilliaryData['OST_LINE']['EXPECTED_ECHO_SHIFT'] = OST_LINE[87:90].uint
    ancilliaryData['OST_LINE']['WINDOW_LEFT_SHIFT'] = OST_LINE[90:93].uint
    ancilliaryData['OST_LINE']['WINDOW_RIGHT_SHIFT'] = OST_LINE[93:96].uint
    ancilliaryData['OST_LINE']['SPARE4'] = OST_LINE[96:128].uint
  return ancilliaryData
########################
#
# Auxiliary file functions
#
#######################
def parseAuxFile(fname, df=False):
  """
  length from label filke
  This function will read in the auxilliary file and return a pandas
  dataframe with the necessary data
  Saving the below NUMPY dtype formatting in the event I can actually
  get it to work...
  dt = np.dtype([('SCET_BLOCK_WHOLE', 'u4'),
                ('SCET_BLOCK_FRAC', 'u2'),
                ('EPHEMERIS_TIME', 'f8'),
                ('GEOMETRY_EPOCH', 'U23'),
                ('SOLAR_LONGITUDE', 'f8'),
                ('ORBIT_NUMBER', 'i4'),
                ('X_MARS_SC_POSITION_VECTOR', 'f8'),
                ('Y_MARS_SC_POSITION_VECTOR', 'f8'),
                ('Z_MARS_SC_POSITION_VECTOR', 'f8'),
                ('SPACECRAFT_ALTITUDE', 'f8'),
                ('SUB_SC_EAST_LONGITUDE', 'f8'),
                ('SUB_SC_PLANETOCENTRIC_LATITUDE', 'f8'),
                ('SUB_SC_PLANETOGRAPHIC_LATITUDE', 'f8'),
                ('X_MARS_SC_VELOCITY_VECTOR', 'f8'),
                ('Y_MARS_SC_VELOCITY_VECTOR', 'f8'),
                ('Z_MARS_SC_VELOCITY_VECTOR', 'f8'),
                ('MARS_SC_RADIAL_VELOCITY', 'f8'),
                ('MARS_SC_TANGENTIAL_VELOCITY', 'f8'),
                ('LOCAL_TRUE_SOLAR_TIME', 'f8'),
                ('SOLAR_ZENITH_ANGLE', 'f8'),
                ('SC_PITCH_ANGLE', 'f8'),
                ('SC_YAW_ANGLE', 'f8'),
                ('SC_ROLL_ANGLE', 'f8'),
                ('MRO_SAMX_INNER_GIMBAL_ANGLE', 'f8'),
                ('MRO_SAMX_OUTER_GIMBAL_ANGLE', 'f8'),
                ('MRO_SAPX_INNER_GIMBAL_ANGLE', 'f8'),
                ('MRO_SAPX_OUTER_GIMBAL_ANGLE', 'f8'),
                ('MRO_HGA_INNER_GIMBAL_ANGLE', 'f8'),
                ('MRO_HGA_OUTER_GIMBAL_ANGLE', 'f8'),
                ('DES_TEMP', 'f4'),
                ('DES_5V', 'f4'),
                ('DES_12V', 'f4'),
                ('DES_2V5', 'f4'),
                ('RX_TEMP', 'f4'),
                ('TX_TEMP', 'f4'),
                ('TX_LEV', 'f4'),
                ('TX_CURR', 'f4'),
                ('CORRUPTED_DATA_FLAG', 'i2')
               ]
               )
  """
  #
  # Set up dictionary
  #
  auxData ={'SCET_BLOCK_WHOLE': [],
            'SCET_BLOCK_FRAC': [],
            'EPHEMERIS_TIME': [],
            'ELAPSED_TIME': [],
            'GEOMETRY_EPOCH': [],
            'SOLAR_LONGITUDE': [],
            'ORBIT_NUMBER': [],
            'X_MARS_SC_POSITION_VECTOR': [],
            'Y_MARS_SC_POSITION_VECTOR': [],
            'Z_MARS_SC_POSITION_VECTOR': [],
            'SPACECRAFT_ALTITUDE': [],
            'SUB_SC_EAST_LONGITUDE': [],
            'SUB_SC_PLANETOCENTRIC_LATITUDE': [],
            'SUB_SC_PLANETOGRAPHIC_LATITUDE': [],
            'X_MARS_SC_VELOCITY_VECTOR': [],
            'Y_MARS_SC_VELOCITY_VECTOR': [],
            'Z_MARS_SC_VELOCITY_VECTOR': [],
            'MARS_SC_RADIAL_VELOCITY': [],
            'MARS_SC_TANGENTIAL_VELOCITY': [],
            'LOCAL_TRUE_SOLAR_TIME': [],
            'SOLAR_ZENITH_ANGLE': [],
            'SC_PITCH_ANGLE': [],
            'SC_YAW_ANGLE': [],
            'SC_ROLL_ANGLE': [],
            'MRO_SAMX_INNER_GIMBAL_ANGLE': [],
            'MRO_SAMX_OUTER_GIMBAL_ANGLE': [],
            'MRO_SAPX_INNER_GIMBAL_ANGLE': [],
            'MRO_SAPX_OUTER_GIMBAL_ANGLE': [],
            'MRO_HGA_INNER_GIMBAL_ANGLE': [],
            'MRO_HGA_OUTER_GIMBAL_ANGLE': [],
            'DES_TEMP': [],
            'DES_5V': [],
            'DES_12V': [],
            'DES_2V5': [],
            'RX_TEMP': [],
            'TX_TEMP': [],
            'TX_LEV': [],
            'TX_CURR': [],
            'CORRUPTED_DATA_FLAG': []
          }
  #
  # Each record is composed of 267 bytes
  #
  recLen = 267
  if os.path.isfile(fname):
    _file = open(fname, 'rb')
    fsize = os.path.getsize(fname)
    for _i in range(int(fsize/recLen)): # Go through all the rows
      _file.seek(_i*recLen)
      rawData = _file.read(recLen)
      auxData['SCET_BLOCK_WHOLE'].append(struct.unpack(">I", rawData[0:4])[0])
      auxData['SCET_BLOCK_FRAC'].append(struct.unpack(">H", rawData[4:6])[0])
      auxData['EPHEMERIS_TIME'].append(struct.unpack(">d", rawData[6:14])[0])
      auxData['ELAPSED_TIME'].append(auxData['EPHEMERIS_TIME'][_i] - auxData['EPHEMERIS_TIME'][0])
      auxData['GEOMETRY_EPOCH'].append(rawData[14:37].decode("utf-8"))
      auxData['SOLAR_LONGITUDE'].append(struct.unpack(">d", rawData[37:45])[0])
      auxData['ORBIT_NUMBER'].append(struct.unpack(">i", rawData[45:49])[0])
      auxData['X_MARS_SC_POSITION_VECTOR'].append(struct.unpack(">d", rawData[49:57])[0])
      auxData['Y_MARS_SC_POSITION_VECTOR'].append(struct.unpack(">d", rawData[57:65])[0])
      auxData['Z_MARS_SC_POSITION_VECTOR'].append(struct.unpack(">d", rawData[65:73])[0])
      auxData['SPACECRAFT_ALTITUDE'].append(struct.unpack(">d", rawData[73:81])[0])
      auxData['SUB_SC_EAST_LONGITUDE'].append(struct.unpack(">d", rawData[81:89])[0])
      auxData['SUB_SC_PLANETOCENTRIC_LATITUDE'].append(struct.unpack(">d", rawData[89:97])[0])
      auxData['SUB_SC_PLANETOGRAPHIC_LATITUDE'].append(struct.unpack(">d", rawData[97:105])[0])
      auxData['X_MARS_SC_VELOCITY_VECTOR'].append(struct.unpack(">d", rawData[105:113])[0])
      auxData['Y_MARS_SC_VELOCITY_VECTOR'].append(struct.unpack(">d", rawData[113:121])[0])
      auxData['Z_MARS_SC_VELOCITY_VECTOR'].append(struct.unpack(">d", rawData[121:129])[0])
      auxData['MARS_SC_RADIAL_VELOCITY'].append(struct.unpack(">d", rawData[129:137])[0])
      auxData['MARS_SC_TANGENTIAL_VELOCITY'].append(struct.unpack(">d", rawData[137:145])[0])
      auxData['LOCAL_TRUE_SOLAR_TIME'].append(struct.unpack(">d", rawData[145:153])[0])
      auxData['SOLAR_ZENITH_ANGLE'].append(struct.unpack(">d", rawData[153:161])[0])
      auxData['SC_PITCH_ANGLE'].append(struct.unpack(">d", rawData[161:169])[0])
      auxData['SC_YAW_ANGLE'].append(struct.unpack(">d", rawData[169:177])[0])
      auxData['SC_ROLL_ANGLE'].append(struct.unpack(">d", rawData[177:185])[0])
      auxData['MRO_SAMX_INNER_GIMBAL_ANGLE'].append(struct.unpack(">d", rawData[185:193])[0])
      auxData['MRO_SAMX_OUTER_GIMBAL_ANGLE'].append(struct.unpack(">d", rawData[193:201])[0])
      auxData['MRO_SAPX_INNER_GIMBAL_ANGLE'].append(struct.unpack(">d", rawData[201:209])[0])
      auxData['MRO_SAPX_OUTER_GIMBAL_ANGLE'].append(struct.unpack(">d", rawData[209:217])[0])
      auxData['MRO_HGA_INNER_GIMBAL_ANGLE'].append(struct.unpack(">d", rawData[217:225])[0])
      auxData['MRO_HGA_OUTER_GIMBAL_ANGLE'].append(struct.unpack(">d", rawData[225:233])[0])
      auxData['DES_TEMP'].append(struct.unpack(">f", rawData[233:237])[0])
      auxData['DES_5V'].append(struct.unpack(">f", rawData[237:241])[0])
      auxData['DES_12V'].append(struct.unpack(">f", rawData[241:245])[0])
      auxData['DES_2V5'].append(struct.unpack(">f", rawData[245:249])[0])
      auxData['RX_TEMP'].append(struct.unpack(">f", rawData[249:253])[0])
      auxData['TX_TEMP'].append(struct.unpack(">f", rawData[253:257])[0])
      auxData['TX_LEV'].append(struct.unpack(">f", rawData[257:261])[0])
      auxData['TX_CURR'].append(struct.unpack(">f", rawData[261:265])[0])
      auxData['CORRUPTED_DATA_FLAG'].append(struct.unpack(">h", rawData[265:267])[0])
    #
    # Check if wanting dataframe returned
    #
    if df:
      auxData = pd.DataFrame.from_dict(auxData)
    return auxData

#########################3
#
# Chirp Functions
#
##########################
def detChirpFiles(TxTemp, RxTemp, chirp='ref'):
  """
  This function determines the appropriate calibrated chirp file to use
  for range compression and returns the decoded calibrated chirp.
  This is solely based off the temperatures of the TxTemp and RxTemp.
  """
  if chirp == 'ref' or chirp == 'vib':
    calibRoot = '../calib/'
    calibName = 'reference_chirp'
    ext = '.dat'
    TxCalNames = ['m20tx', 'm15tx', 'm10tx', 'm05tx',
                  'p00tx', 'p20tx', 'p40tx', 'p60tx']
    RxCalNames = ['m20rx', 'p00rx', 'p20rx', 'p40rx',
                  'p60rx']
    #
    # Define vectors for Tx and Rx temps
    #
    Tx = [-20, -15, -10, -5, 0, 20, 40, 60]
    Rx = [-20, 0, 20, 40, 60]
    calibChirpFiles = []
    TxDiff = []
    RxDiff = []
    #
    # Find distance
    #
    TxDiff[:] = [abs(x - TxTemp) for x in Tx]
    RxDiff[:] = [abs(x - RxTemp) for x in Rx]
    #
    # Find the indices of the closest Tx and Rx value
    #
    calibTx = TxCalNames[TxDiff.index(min(TxDiff))]
    calibRx = RxCalNames[RxDiff.index(min(RxDiff))]
    #
    # Construct File name
    #
    calChirpFile = calibRoot + calibName + '_' + \
                   TxCalNames[TxDiff.index(min(TxDiff))] + '_' + \
                   RxCalNames[RxDiff.index(min(RxDiff))] + ext
    if os.path.isfile(calChirpFile):
      calChirp = np.fromfile(calChirpFile, dtype='<f')
      if chirp == 'ref':
        real = calChirp[:2048]
        imag = calChirp[2048:]
      elif chirp == 'vib':
        #
        # Add a zero to the end to mimic missing Nyquist
        #
        real = np.zeros(4096, float)
        real[0:2048] = calChirp[:2048]
        #
        # Drop the last sample to properly stitch Nyquist and reverse
        #
        real[2049:] = np.flipud(calChirp[1:2048])
        #
        # Add a zero to the front
        #
        imag = np.zeros(4096, float)
        imag[0:2048] = calChirp[2048:]
        #
        # Drop the last sample and reverse and change sign
        #
        imag[2049:] = -1 * np.flipud(calChirp[2049:])
      calChirp = real + 1j*imag
      return calChirp
    else:
      print('Calibrated chirp file not found...exiting.')
      sys.exit()
  else:
   #
    # Use ideal chirp from UPB
    #
    delay_res = 135.00e-6 / 3600.000
    sharad_ipp = 1.0 / 700.28
    fhi = 25.00e6
    flo = 15.00e6
    plen = 85.05e-6
    nsamp = plen / delay_res
    fslope = (flo - fhi) / plen
    ctime = np.arange(0,nsamp) * delay_res
    arg = 2.0*np.pi*ctime*(fhi+fslope*ctime/2.0)
    ideal_chirp = np.zeros(3600, complex)
    ideal_chirp[:int(nsamp)] = np.sin(arg)
    ideal_chirp_FFT = np.fft.fft(ideal_chirp)
    if chirp == 'UPB':
      #
      # Load cal_filter.dat
      #
      cal_filter = np.fromfile('../calib/cal_filter.dat', '<f')
      cal_filter = cal_filter[:1800] + 1j*cal_filter[1800:]
      cal_filter = np.roll(cal_filter, 900)
      calChirp = np.zeros(3600, complex)
      calChirp[1800:] = cal_filter
      calChirp = calChirp*ideal_chirp_FFT
      return calChirp
    else:
      calChirp = ideal_chirp_FFT
      return calChirp

##########################################
#
# EDR Functions
#
##########################################
def parseFileName(_file):
  SSInstrMode = { 'SS01': {'Mode': 'SS01', 'Presum': 32, 'BitsPerSample': 8},
                  'SS02': {'Mode': 'SS02','Presum': 28, 'BitsPerSample': 6},
                  'SS03': {'Mode': 'SS03','Presum': 16, 'BitsPerSample': 4},
                  'SS04': {'Mode': 'SS04','Presum': 8, 'BitsPerSample': 8},
                  'SS05': {'Mode': 'SS05','Presum': 4, 'BitsPerSample': 6},
                  'SS06': {'Mode': 'SS06','Presum': 2, 'BitsPerSample': 4},
                  'SS07': {'Mode': 'SS07','Presum': 1, 'BitsPerSample': 8},
                  'SS08': {'Mode': 'SS08','Presum': 32, 'BitsPerSample': 6},
                  'SS09': {'Mode': 'SS09','Presum': 28, 'BitsPerSample': 4},
                  'SS10': {'Mode': 'SS10','Presum': 16, 'BitsPerSample': 8},
                  'SS11': {'Mode': 'SS11','Presum': 8, 'BitsPerSample': 6},
                  'SS12': {'Mode': 'SS12','Presum': 4, 'BitsPerSample': 4},
                  'SS13': {'Mode': 'SS13','Presum': 2, 'BitsPerSample': 8},
                  'SS14': {'Mode': 'SS14','Presum': 1, 'BitsPerSample': 6},
                  'SS15': {'Mode': 'SS15','Presum': 32, 'BitsPerSample': 4},
                  'SS16': {'Mode': 'SS16','Presum': 28, 'BitsPerSample': 8},
                  'SS17': {'Mode': 'SS17','Presum': 16, 'BitsPerSample': 6},
                  'SS18': {'Mode': 'SS18','Presum': 8, 'BitsPerSample': 4},
                  'SS19': {'Mode': 'SS19','Presum': 4, 'BitsPerSample': 8},
                  'SS20': {'Mode': 'SS20','Presum': 2, 'BitsPerSample': 6},
                  'SS21': {'Mode': 'SS21','Presum': 1, 'BitsPerSample': 4},
               }
  #
  # Receive only
  #
  ROInstrMode = { 'RO01': {'Presum': 32, 'BitsPerSample': 8},
                  'RO02': {'Presum': 28, 'BitsPerSample': 6},
                  'RO03': {'Presum': 16, 'BitsPerSample': 4},
                  'RO04': {'Presum': 8, 'BitsPerSample': 8},
                  'RO05': {'Presum': 4, 'BitsPerSample': 6},
                  'RO06': {'Presum': 2, 'BitsPerSample': 4},
                  'RO07': {'Presum': 1, 'BitsPerSample': 8},
                  'RO08': {'Presum': 32, 'BitsPerSample': 6},
                  'RO09': {'Presum': 28, 'BitsPerSample': 4},
                  'RO10': {'Presum': 16, 'BitsPerSample': 8},
                  'RO11': {'Presum': 8, 'BitsPerSample': 6},
                  'RO12': {'Presum': 4, 'BitsPerSample': 4},
                  'RO13': {'Presum': 2, 'BitsPerSample': 8},
                  'RO14': {'Presum': 1, 'BitsPerSample': 6},
                  'RO15': {'Presum': 32, 'BitsPerSample': 4},
                  'RO16': {'Presum': 28, 'BitsPerSample': 8},
                  'RO17': {'Presum': 16, 'BitsPerSample': 6},
                  'RO18': {'Presum': 8, 'BitsPerSample': 4},
                  'RO19': {'Presum': 4, 'BitsPerSample': 8},
                  'RO20': {'Presum': 2, 'BitsPerSample': 6},
                  'RO21': {'Presum': 1, 'BitsPerSample': 4}
                }
  #
  # Get the file basename
  #
  bname = os.path.basename(_file)
  #
  # Now split the file basename on _
  #
  bname = bname.split('_')
  TransID = bname[1]
  OSTLine = bname[2]
  OperMode = bname[3].upper()
  PRF = bname[4]
  Version = bname[5].split('.')[0]
  if OperMode[0:2] == 'SS':
    OperMode = SSInstrMode[OperMode]
  elif OperMode[0:2] == 'RO':
    OperMode = ROInstrMode[OperMode]
  return TransID, OSTLine, OperMode, PRF, Version

def readEDRrecord(_file, record, recLen, bps):
  """
    This function reads a single EDR record and returns the ancilliary
    and science data. Data are read and decoded based off their bit resolution
    Input:
      _file: File object identifier pointing to the data file containing the binary data
    Output:
      ancil: Data dictionary containing all the information parsed from the ancilliary data for the record
      echoes: Decoded science data
  """
  #
  # Make sure we are at the beginning of a record
  #
  _file.seek(0)
  #
  # Now fast forward to the appropriate location
  #
  _file.seek(record*recLen)
  #
  # Step 1: Read in binary data
  #
  rawData = _file.read(recLen)
  #
  # Separate ancil and echo data
  #
  # Step 2: Read and parse ancilliary data
  #
  ancil = rawData[:186]
  #
  # Step 3: Separate Actual echoes
  #
  echoes = rawData[186:]
  #
  # Okay, now decode the data
  # For 8-bit samples, there is a sample for every byte
  # For 6-bit samples, there are 2 samples for every 3 bytes
  # For 4-bit samples, there are 2 samples for every byte
  #
  # Making the vector have 4096 rows is to allow for proper decovolution of the
  # calibrated chirp
    #
  decoded_data = np.zeros(3600, int)
  #
  # Step 4: Decode the science data based on the bit resolution
  #
  # Break Byte stream into Bit Stream
  # This isn't necessary for the 8-bit resolution, but to keep the code
  # clean, split the bytes up regardless of resolution
  #
  b = bitstring.BitArray(echoes)
  for _j in range(0, len(b), bps):
    decoded_data[int(_j/bps)] = b[_j:_j+bps].int
  return decoded_data, ancil

def decompressSciData(data, compression, presum, bps, SDI):
  """
    This function decompresses the data based off page 8 of the
    SHALLOW RADAR EXPERIMENT DATA RECORD SOFTWARE INTERFACE SPECIFICATION.
    If the compression type is 'Static'
       U = C*2**S / N
         C is the Compressed Data
         S = L - R + 8
            L is base 2 log of N rounded UP to the nearest integer
            R is the bit resolution of the compressed values
         N is the number of pre-summed echoes
    If the compression type is 'dynamic'
      NOT WORKING YET
      U = C*2**S/N
        C is the compressed data
        N is the number of pre-summed echoes
        S = SDI for SDI <= 5
        S = SDI-6 for 5 < SDI <= 16
        S = SDI-16 for SDI > 16
          where SDI is the SDI_BIT_FIELD parameter from ANCILLIARY DATA
  """
  # Step 5: Decompress the data
  # Note: Only Static decompression works at this time
  #
  if compression is not True:
    compression = 'STATIC'
  else:
    compression = 'DYNAMIC'
  if compression == 'STATIC' or compression == 'DYNAMIC':
    #
    # Handle fixed scaling
    #
    if compression == 'STATIC': # Static scaling
      L = np.ceil(np.log2(int(presum)))
      R = bps
      S = L - R + 8
      N = presum
      decomp = np.power(2, S) / N
      decompressed_data = data * decomp
    elif compression == True:#dynamic scaling
      N = presum
      if SDI <= 5:
        S = SDI
      elif 5 < SDI <= 16:
        S = SDI - 6
      elif SDI > 16:
        S = SDI - 16
      decompressed_data = data * (np.power(2, S) / N)
    return decompressed_data
  else:
    print('Decompression Error: Compression Type {} not understood'.format(compression))
    return

#def calibrateData(data,
  #
  # The following factors affect the performance of SHARAD
  #     1) Transmitter and Receiver Temperatures
  #     2) Spacecraft attitude (i.e. Spacecraft roll-angle)
  #             Fabrizio says to add or subtract 180 ... and then determine which works.:q

  #     3) Spacecraft configuration (i.e. orientation of large moving parts)
  #     4) Propagation of the signal through Mars' ionosphere
  #


def detSCConf(SAMX_IG, SAPX_IG, OGA):
  print(SAMX_IG, SAPX_IG, OGA)
  return 0

########################################
#
# Label file functions (largely defunct as of 26 Sept 2018)
#
########################################
def parseLBLFile(fname):
  #
  # Define instrument modes
  #
  #
  # Subsurface sounding
  #
  SSInstrMode = { 'SS01': {'Mode': 'SS01', 'Presum': 32, 'BitsPerSample': 8},
                  'SS02': {'Mode': 'SS02','Presum': 28, 'BitsPerSample': 6},
                  'SS03': {'Mode': 'SS03','Presum': 16, 'BitsPerSample': 4},
                  'SS04': {'Mode': 'SS04','Presum': 8, 'BitsPerSample': 8},
                  'SS05': {'Mode': 'SS05','Presum': 4, 'BitsPerSample': 6},
                  'SS06': {'Mode': 'SS06','Presum': 2, 'BitsPerSample': 4},
                  'SS07': {'Mode': 'SS07','Presum': 1, 'BitsPerSample': 8},
                  'SS08': {'Mode': 'SS08','Presum': 32, 'BitsPerSample': 6},
                  'SS09': {'Mode': 'SS09','Presum': 28, 'BitsPerSample': 4},
                  'SS10': {'Mode': 'SS10','Presum': 16, 'BitsPerSample': 8},
                  'SS11': {'Mode': 'SS11','Presum': 8, 'BitsPerSample': 6},
                  'SS12': {'Mode': 'SS12','Presum': 4, 'BitsPerSample': 4},
                  'SS13': {'Mode': 'SS13','Presum': 2, 'BitsPerSample': 8},
                  'SS14': {'Mode': 'SS14','Presum': 1, 'BitsPerSample': 6},
                  'SS15': {'Mode': 'SS15','Presum': 32, 'BitsPerSample': 4},
                  'SS16': {'Mode': 'SS16','Presum': 28, 'BitsPerSample': 8},
                  'SS17': {'Mode': 'SS17','Presum': 16, 'BitsPerSample': 6},
                  'SS18': {'Mode': 'SS18','Presum': 8, 'BitsPerSample': 4},
                  'SS19': {'Mode': 'SS19','Presum': 4, 'BitsPerSample': 8},
                  'SS20': {'Mode': 'SS20','Presum': 2, 'BitsPerSample': 6},
                  'SS21': {'Mode': 'SS21','Presum': 1, 'BitsPerSample': 4},
               }
  #
  # Receive only
  #
  ROInstrMode = { 'RO01': {'Presum': 32, 'BitsPerSample': 8},
                  'RO02': {'Presum': 28, 'BitsPerSample': 6},
                  'RO03': {'Presum': 16, 'BitsPerSample': 4},
                  'RO04': {'Presum': 8, 'BitsPerSample': 8},
                  'RO05': {'Presum': 4, 'BitsPerSample': 6},
                  'RO06': {'Presum': 2, 'BitsPerSample': 4},
                  'RO07': {'Presum': 1, 'BitsPerSample': 8},
                  'RO08': {'Presum': 32, 'BitsPerSample': 6},
                  'RO09': {'Presum': 28, 'BitsPerSample': 4},
                  'RO10': {'Presum': 16, 'BitsPerSample': 8},
                  'RO11': {'Presum': 8, 'BitsPerSample': 6},
                  'RO12': {'Presum': 4, 'BitsPerSample': 4},
                  'RO13': {'Presum': 2, 'BitsPerSample': 8},
                  'RO14': {'Presum': 1, 'BitsPerSample': 6},
                  'RO15': {'Presum': 32, 'BitsPerSample': 4},
                  'RO16': {'Presum': 28, 'BitsPerSample': 8},
                  'RO17': {'Presum': 16, 'BitsPerSample': 6},
                  'RO18': {'Presum': 8, 'BitsPerSample': 4},
                  'RO19': {'Presum': 4, 'BitsPerSample': 8},
                  'RO20': {'Presum': 2, 'BitsPerSample': 6},
                  'RO21': {'Presum': 1, 'BitsPerSample': 4}
                }
  #
  # Now parse LBL File
  #
  if os.path.isfile(fname):
    #
    # Initialize dictionary
    #
    lblDic = {'INSTR_MODE_ID': [],
              'PRI': [],
              'GAIN_CONTROL': [],
              'COMPRESSION': [],
              'RECORD_BYTES': [],
              'FILE_RECORDS': []
              }
    with open(fname) as f:
      lines = f.readlines()
    #
    # Remove end of line characters from list
    #
    lines = [x.rstrip('\n') for x in lines]
    lineCount = len(lines)
    #
    # Remove all blank rows
    #
    lines = [x for x in lines if x]
    print("{} empty lines removed.".format(lineCount - len(lines)))
    lineCount = len(lines)
    #
    # Remove comments
    #
    lines = [x for x in lines if "/*" not in x]
    print("{} comment lines removed.".format(lineCount - len(lines)))
    lineCount = len(lines)
    #
    # Start parsing
    #
    print("Parsing {} lines in LBL file.".format(lineCount))
    for _i in range(lineCount):
      #
      # For this simple test all I should need out of the LBL file are:
      #  INSTRUMENT_MODE_ID
      #
      if lines[_i].split('=')[0].strip() == 'INSTRUMENT_MODE_ID':
        lblDic['INSTR_MODE_ID'] = lines[_i].split('=')[1].strip()
      if lines[_i].split('=')[0].strip() == 'MRO:PULSE_REPETITION_INTERVAL':
        lblDic['PRI'] = lines[_i].split('=')[1].strip()
      if lines[_i].split('=')[0].strip() == 'MRO:MANUAL_GAIN_CONTROL':
        lblDic['GAIN_CONTROL'] = lines[_i].split('=')[1].strip()
      if lines[_i].split('=')[0].strip() == 'MRO:COMPRESSION_SELECTION_FLAG':
        lblDic['COMPRESSION'] = lines[_i].split('=')[1].strip().strip('"')
      if lines[_i].split('=')[0].strip() == 'RECORD_BYTES':
        if lblDic['RECORD_BYTES'] == []:
          lblDic['RECORD_BYTES'] = int(lines[_i].split('=')[1].strip())
      if lines[_i].split('=')[0].strip() == 'FILE_RECORDS':
        if lblDic['FILE_RECORDS'] == []:
          lblDic['FILE_RECORDS'] = int(lines[_i].split('=')[1].strip())
    #
    # Find the instrument mode
    #
    if lblDic['INSTR_MODE_ID'][0:2] == 'SS':
      lblDic['INSTR_MODE_ID'] = SSInstrMode[lblDic['INSTR_MODE_ID']]
    elif lblDic['INSTR_MODE_ID'][0:2] == 'RO':
      lblDic['INSTR_MODE_ID'] = ROInstrMode[lblDic['INSTR_MODE_ID']]
    return lblDic
  else:
    print("{} file not found. Please check path and try again.".format(fname))
    return
###########################
#
# Plotting functions
#
###########################
def makeAuxPlots(df):
  """
    Something I can do later
  """
  f, axarr = plt.subplots(2, sharex=True)
  f.suptitle('Sharing X axis')
  X = df['ELAPSED_TIME']
  axarr[0].plot(X, df['SOLAR_LONGITUDE'], 'k.')
  axarr[1].plot(X, df['X_MARS_SC_POSITION_VECTOR'], 'k.')
  axarr[1].plot(X, df['Y_MARS_SC_POSITION_VECTOR'], 'r.')
  axarr[1].plot(X, df['Z_MARS_SC_POSITION_VECTOR'], 'b.')
  plt.show()
  return

def plotTimeSeriesDiag(dt, sci, fname='RawDiag.png', transparent=False):
  print('Producing raw data diagnostic plot...')
  t = np.arange(0, len(sci)) * dt
  spectra = np.fft.fft(sci, dt) / len(sci)
  freqs = np.fft.fftfreq(len(sci), dt)
  #
  # Plot Original Time Series
  #
  plt.subplot(3,1,1)
  plt.plot(t*1e6, np.real(sci))
  plt.title('Raw Signal')
  plt.xlabel('Time (us)')
  plt.ylabel('Amplitude')
  #
  # Plot Spectra
  #
  plt.subplot(3,1,2)
  plt.plot(freqs*1e6, np.abs(np.fft.fftshift(spectra)))
  plt.title("Modulus Amplitude Spectrum (|S|)")
  plt.xlabel('Magnitude')
  plt.ylabel('Frequency (MHz)')
  #
  # Plot Angle Spectrum
  #
  plt.subplot(3,1,3)
  plt.plot(freqs*1e6, np.unwrap(np.angle(np.fft.fftshift(spectra)))) 
  plt.title('Unwrapped Phase Spectrum')
  plt.xlabel('Angle')
  plt.ylabel('Frequency (MHz)')
  #
  # Final formatting
  #
  plt.tight_layout()
  print('Saving raw diagnostic plot as {}'.format(fname))
  plt.savefig(fname, dpi=500, transparent=transparent)
  #
  # Clear the plot
  #
  plt.clf()
  return 

def bytescl(array, mindata=None, maxdata=None, top=255):
  #
  # Byte scaling algorithm for greater contrast
  #
  if mindata is None: mindata = np.nanmin(array)
  if maxdata is None: maxdata = np.nanmax(array)
  scl = np.maximum(np.minimum(((top+0.9999)*(array-mindata)/(maxdata-mindata)).astype(np.int16), top),0)
  return scl

def rdr2san(data, fname='rdr2san', maxdb=0, top=255):
  pow_out = np.power(np.abs(data), 2)                   # Convert data to power
  db = 10 * np.log10(pow_out)                           # decibels
  maxdb = np.amax(db)
  sig = db/maxdb*255
  sig[np.where(sig < 0)] = 0.                           # Zero out values below noise floor
  sig[np.where(sig > 255)] = 255.                       # Clip values greater than MAXDB
  imName = '../runs/' + str(fname) + '.eps'
  plt.imshow(sig, cmap='gray')
  plt.savefig(imName, format='eps', dpi=1000)
  plt.show()
  return

def plotEDR(data, fname='plotEDR_new', ptype='Amp', thres=0.3, rel=False, cmap='gray', dpi=1000):
  bmpName = '../output/' + str(fname)
  pngName = '../output/' + str(fname)
  if ptype == 'Amp':
    pic = np.real(data)
    thresVal = 0.0
    bmpend = '_amp.bmp'
    pngend = '_amp.png' 
  elif ptype == 'Pow' or ptype == 'dB':
    pic = np.power(np.real(data), 2)
    if pytype == 'Pow':
      bmpend = '_pow.bmp'
      pngend = '_pow.png' 
      thresVal = 0.0
    elif ptype == 'dB':
      bmpend = '_dB.bmp'
      pngend = '_dB.png' 
      thresVal = -99.0
      if thres > 0:
        print('[WARNING] Threshold value of {} is too high for decibels...'.format(thres))
        print('[WARNING] Setting threshold to -20 dB') 
        thres = -20
  #
  # Determine scaling
  #
  if rel == True:
    mx = np.amax(pic, axis=0)
  else:
    mx = np.amax(pic)
  #
  # Scale
  # 
  if ptype == 'dB':
    pic = 10*np.log10(Pow/mx)
  pic[np.where(pic < thres)] = thresVal
  #
  # Ok...now plot
  #
  plt.imshow(pic, cmap=cmap, vmin=thres)
  plt.savefig(pngName, dpi=dpi)
  plt.imsave(bmpName, pic, cmap=cmap)
  return


def plotFirstReturn(data, type='Amp', sidelobe=False, title='First Return', fname='FirstReturn', dpi=500):
  before = 10
  after = 20
  if type == 'Amp':
    data = np.real(data)
    ylabel = 'Amplitude'
  elif type == 'Pow':
    data = np.power(np.real(data),2)
    ylabel = 'Power'
  elif type == 'dB':
    data = np.power(np.real(data), 2)
    mx = np.amax(data)
    data = 10 * np.log10(data/mx)
    ylabel = 'dB (relative to maximum)'
  idx = np.where(np.abs(data) == np.amax(np.abs(data)))
  t0 = idx[0] * 0.0375
  t = np.arange(idx[0]-before, idx[0]+after) * 0.0375
  t = t - t0
  #
  # Clear any plots
  #
  plt.clf()
  #
  # Now plot away
  #
  plt.plot(t, data[int(idx[0]-before):int(idx[0]+after)])
  if sidelobe:
    #
    # Theorectically, sidelobes should appear at around 0.24 us after the signal
    #
    plt.axvline(0.24, color='red')
  plt.xlabel('Time (us)')
  plt.ylabel(ylabel)
  plt.title(title)
  pngName = fname + '.png'
  plt.savefig(pngName, dpi=dpi)
  return

