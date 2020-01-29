import glob, sys, os, struct, argparse, bitstring
import pandas as pd
import numpy as np

###############################################
#
# Log Functions
#
###############################################
def writeLog(log, string, verb=False):
  if verb:
    print(string)
  if log is not None:
    log.write(string + '\n')
###############################################
#
# Argument Parser
#
###############################################
def parseargs(prog, vers):
  #
  # Set the default values
  #
  outDir = ['../out/']
  verb = False
  #
  # Initiate the parser
  #
  parser = argparse.ArgumentParser(description=str(prog + ' ' + vers))
  #
  # Ordered Arguments
  #
  parser.add_argument('lblFile', type=str, nargs=1,
                      help=str('Full path to PDS label file')) 
  #
  # Optional Arguments
  #
  parser.add_argument('-o', '--outDir', nargs=1, default=outDir, type=str,
                     help=str('Desired output directory'))
  #
  # Obligatory verbosity level and diagnostic options
  #
  parser.add_argument('-v', '--verbose', action="store_true", default=verb,
                     help=str('Print to screen as well as log file'))
  #
  # Print help
  #
  if len(sys.argv[1:]) == 0:
    parser.print_help()
    parser.exit()
  #
  # Parse the arguments
  #
  args = parser.parse_args()
  lblFile = args.lblFile[0]
  outDir = args.outDir[0]
  verb = args.verbose
  #
  # CHECK LBL FILE
  # 
  if not os.path.isfile(lblFile):
    print('The label file {} could not be found!'.format(lblFile))
    parser.print_help()
    parser.exit()
  elif os.path.isfile(lblFile) and not os.access(lblFile, os.R_OK):
    print('The label file {} is present, but is not readable by the current user'.format(lblFile))
    parser.print_help()
    parser.exit()
  #
  # Find science and auxiliary files
  #
  iFiles = findFiles(lblFile)
  #
  #form outDir
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
  return iFiles, oFiles, TransID, OSTLine, OperMode, verb


def findFiles(lblFile):
  #
  # Label File exists
  #
  with open(lblFile) as f:
    content = f.readlines()
    sw = 0
    length = len(content)
    for _i in range(length):
      if '^SCIENCE_TELEMETRY_TABLE' in content[_i]:
        sciFile = content[_i].split('=')[-1].strip().replace('"', '')
        sw += 1
      if '^AUXILIARY_DATA_TABLE' in content[_i]:
        auxFile = content[_i].split('=')[-1].strip().replace('"', '')
        sw += 1
      if sw == 2:
        break
  path, _ = os.path.split(lblFile)
  #
  # Find SCIENCE FILE
  #
  if os.path.isfile(path + '/' + sciFile.upper()):
    sciFile = path + '/' + sciFile.upper()
  elif os.path.isfile(path + '/' + sciFile.lower()):
    sciFile = path + '/' + sciFile.lower()
  else:
    print('[ERROR]: Science file not found...')
    exit()
  #
  # Find AUXILIARY FILE
  #
  if os.path.isfile(path + '/' + auxFile.upper()):
    auxFile = path + '/' + auxFile.upper()
  elif os.path.isfile(path + '/' + auxFile.lower()):
    auxFile = path + '/' + auxFile.lower()
  else:
    print('[ERROR]: Auxiliary file not found...')
    exit()
  iFiles = {'LABEL': lblFile, 'SCIENCE': sciFile, 'AUX': auxFile}
  return iFiles


def parseFileName(_file):
  #
  # Determine to Bits per Sample
  #
  SSInstrMode = { 'SS01': {'Mode': 'SS01', 'Presum': 32, 'BitsPerSample': 8, 'recLen': 3786},
                  'SS02': {'Mode': 'SS02','Presum': 28, 'BitsPerSample': 6, 'recLen': 2886},
                  'SS03': {'Mode': 'SS03','Presum': 16, 'BitsPerSample': 4, 'recLen': 1986},
                  'SS04': {'Mode': 'SS04','Presum': 8, 'BitsPerSample': 8, 'recLen': 3786},
                  'SS05': {'Mode': 'SS05','Presum': 4, 'BitsPerSample': 6, 'recLen': 2886},
                  'SS06': {'Mode': 'SS06','Presum': 2, 'BitsPerSample': 4, 'recLen': 1986},
                  'SS07': {'Mode': 'SS07','Presum': 1, 'BitsPerSample': 8, 'recLen': 3786},
                  'SS08': {'Mode': 'SS08','Presum': 32, 'BitsPerSample': 6, 'recLen': 2886},
                  'SS09': {'Mode': 'SS09','Presum': 28, 'BitsPerSample': 4, 'recLen': 1986},
                  'SS10': {'Mode': 'SS10','Presum': 16, 'BitsPerSample': 8, 'recLen': 3786},
                  'SS11': {'Mode': 'SS11','Presum': 8, 'BitsPerSample': 6, 'recLen': 2886},
                  'SS12': {'Mode': 'SS12','Presum': 4, 'BitsPerSample': 4, 'recLen': 1986},
                  'SS13': {'Mode': 'SS13','Presum': 2, 'BitsPerSample': 8, 'recLen': 3786},
                  'SS14': {'Mode': 'SS14','Presum': 1, 'BitsPerSample': 6, 'recLen': 2886},
                  'SS15': {'Mode': 'SS15','Presum': 32, 'BitsPerSample': 4, 'recLen': 1986},
                  'SS16': {'Mode': 'SS16','Presum': 28, 'BitsPerSample': 8, 'recLen': 3786},
                  'SS17': {'Mode': 'SS17','Presum': 16, 'BitsPerSample': 6, 'recLen': 2886},
                  'SS18': {'Mode': 'SS18','Presum': 8, 'BitsPerSample': 4, 'recLen': 1986},
                  'SS19': {'Mode': 'SS19','Presum': 4, 'BitsPerSample': 8, 'recLen': 3786},
                  'SS20': {'Mode': 'SS20','Presum': 2, 'BitsPerSample': 6, 'recLen': 2886},
                  'SS21': {'Mode': 'SS21','Presum': 1, 'BitsPerSample': 4, 'recLen': 1986},
               }
  #
  # Receive only
  #
  ROInstrMode = { 'RO01': {'Presum': 32, 'BitsPerSample': 8, 'recLen': 3786},
                  'RO02': {'Presum': 28, 'BitsPerSample': 6, 'recLen': 2886},
                  'RO03': {'Presum': 16, 'BitsPerSample': 4, 'recLen': 1986},
                  'RO04': {'Presum': 8, 'BitsPerSample': 8, 'recLen': 3786},
                  'RO05': {'Presum': 4, 'BitsPerSample': 6, 'recLen': 2886},
                  'RO06': {'Presum': 2, 'BitsPerSample': 4, 'recLen': 1986},
                  'RO07': {'Presum': 1, 'BitsPerSample': 8, 'recLen': 3786},
                  'RO08': {'Presum': 32, 'BitsPerSample': 6, 'recLen': 2886},
                  'RO09': {'Presum': 28, 'BitsPerSample': 4, 'recLen': 1986},
                  'RO10': {'Presum': 16, 'BitsPerSample': 8, 'recLen': 3786},
                  'RO11': {'Presum': 8, 'BitsPerSample': 6, 'recLen': 2886},
                  'RO12': {'Presum': 4, 'BitsPerSample': 4, 'recLen': 1986},
                  'RO13': {'Presum': 2, 'BitsPerSample': 8, 'recLen': 3786},
                  'RO14': {'Presum': 1, 'BitsPerSample': 6, 'recLen': 2886},
                  'RO15': {'Presum': 32, 'BitsPerSample': 4, 'recLen': 1986},
                  'RO16': {'Presum': 28, 'BitsPerSample': 8, 'recLen': 3786},
                  'RO17': {'Presum': 16, 'BitsPerSample': 6, 'recLen': 2886},
                  'RO18': {'Presum': 8, 'BitsPerSample': 4, 'recLen': 1986},
                  'RO19': {'Presum': 4, 'BitsPerSample': 8, 'recLen': 3786},
                  'RO20': {'Presum': 2, 'BitsPerSample': 6, 'recLen': 2886},
                  'RO21': {'Presum': 1, 'BitsPerSample': 4, 'recLen': 1986}
                }
  #
  # Get the file basename
  #
  bname = os.path.basename(_file).split('_')
  TransID = bname[1]
  OSTLine = bname[2]
  OperMode = bname[3].upper()
  PRF = detPRF(bname[4])
  if OperMode[0:2] == 'SS':
    OperMode = SSInstrMode[OperMode]
  elif OperMode[0:2] == 'RO':
    OperMode = ROInstrMode[OperMode]
  OperMode['PRF'] = PRF
  nrec = int(os.path.getsize(_file) / OperMode['recLen'])
  OperMode['nrec'] = nrec
  return TransID, OSTLine, OperMode


def parseAuxFile(fname, oFile, dic=True, df=False, csv=False, binary=False, saveNP=False):
  #
  # Set up dictionary
  #
  a ={'SCET_BLOCK_WHOLE': [],
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
    _f = open(fname, 'rb')
    fsize = os.path.getsize(fname)
    for _i in range(int(fsize/recLen)): # Go through all the rows
      _f.seek(_i*recLen)
      r = _f.read(recLen)
      a['SCET_BLOCK_WHOLE'].append(struct.unpack(">I", r[0:4])[0])
      a['SCET_BLOCK_FRAC'].append(struct.unpack(">H", r[4:6])[0])
      a['EPHEMERIS_TIME'].append(struct.unpack(">d", r[6:14])[0])
      a['ELAPSED_TIME'].append(a['EPHEMERIS_TIME'][_i] - a['EPHEMERIS_TIME'][0])
      a['GEOMETRY_EPOCH'].append(r[14:37].decode("utf-8"))
      a['SOLAR_LONGITUDE'].append(struct.unpack(">d", r[37:45])[0])
      a['ORBIT_NUMBER'].append(struct.unpack(">i", r[45:49])[0])
      a['X_MARS_SC_POSITION_VECTOR'].append(struct.unpack(">d", r[49:57])[0])
      a['Y_MARS_SC_POSITION_VECTOR'].append(struct.unpack(">d", r[57:65])[0])
      a['Z_MARS_SC_POSITION_VECTOR'].append(struct.unpack(">d", r[65:73])[0])
      a['SPACECRAFT_ALTITUDE'].append(struct.unpack(">d", r[73:81])[0])
      a['SUB_SC_EAST_LONGITUDE'].append(struct.unpack(">d", r[81:89])[0])
      a['SUB_SC_PLANETOCENTRIC_LATITUDE'].append(struct.unpack(">d", r[89:97])[0])
      a['SUB_SC_PLANETOGRAPHIC_LATITUDE'].append(struct.unpack(">d", r[97:105])[0])
      a['X_MARS_SC_VELOCITY_VECTOR'].append(struct.unpack(">d", r[105:113])[0])
      a['Y_MARS_SC_VELOCITY_VECTOR'].append(struct.unpack(">d", r[113:121])[0])
      a['Z_MARS_SC_VELOCITY_VECTOR'].append(struct.unpack(">d", r[121:129])[0])
      a['MARS_SC_RADIAL_VELOCITY'].append(struct.unpack(">d", r[129:137])[0])
      a['MARS_SC_TANGENTIAL_VELOCITY'].append(struct.unpack(">d", r[137:145])[0])
      a['LOCAL_TRUE_SOLAR_TIME'].append(struct.unpack(">d", r[145:153])[0])
      a['SOLAR_ZENITH_ANGLE'].append(struct.unpack(">d", r[153:161])[0])
      a['SC_PITCH_ANGLE'].append(struct.unpack(">d", r[161:169])[0])
      a['SC_YAW_ANGLE'].append(struct.unpack(">d", r[169:177])[0])
      a['SC_ROLL_ANGLE'].append(struct.unpack(">d", r[177:185])[0])
      a['MRO_SAMX_INNER_GIMBAL_ANGLE'].append(struct.unpack(">d", r[185:193])[0])
      a['MRO_SAMX_OUTER_GIMBAL_ANGLE'].append(struct.unpack(">d", r[193:201])[0])
      a['MRO_SAPX_INNER_GIMBAL_ANGLE'].append(struct.unpack(">d", r[201:209])[0])
      a['MRO_SAPX_OUTER_GIMBAL_ANGLE'].append(struct.unpack(">d", r[209:217])[0])
      a['MRO_HGA_INNER_GIMBAL_ANGLE'].append(struct.unpack(">d", r[217:225])[0])
      a['MRO_HGA_OUTER_GIMBAL_ANGLE'].append(struct.unpack(">d", r[225:233])[0])
      a['DES_TEMP'].append(struct.unpack(">f", r[233:237])[0])
      a['DES_5V'].append(struct.unpack(">f", r[237:241])[0])
      a['DES_12V'].append(struct.unpack(">f", r[241:245])[0])
      a['DES_2V5'].append(struct.unpack(">f", r[245:249])[0])
      a['RX_TEMP'].append(struct.unpack(">f", r[249:253])[0])
      a['TX_TEMP'].append(struct.unpack(">f", r[253:257])[0])
      a['TX_LEV'].append(struct.unpack(">f", r[257:261])[0])
      a['TX_CURR'].append(struct.unpack(">f", r[261:265])[0])
      a['CORRUPTED_DATA_FLAG'].append(struct.unpack(">h", r[265:267])[0])
    #
    # Check output
    #
    if csv == True:
      aux = pd.DataFrame.from_dict(a)
      aux.to_csv(oFile)
      return
    if dic == True:
      return a
    elif df == True:
      return pd.DataFrame.from_dict(a)

def parseAncillary(fname):
  """
    
  """
  #
  # Set up dictionary
  #
  a = {'SCET_BLOCK_WHOLE': [],
       'SCET_BLOCK_FRAC': [],
       'TLM_COUNTER': [],
       'FMT_LENGTH': [],
       'SPARE1': [],
       'SCET_OST_WHOLE': [],
       'SCET_OST_FRAC': [],
       'SPARE2': [],
       'OST_LINE_NUMBER': [], 
       'OST_PULSE_REPETITION_INTERVAL': [],
       'OST_PHASE_COMPENSATION_TYPE': [],
       'OST_SPARE1': [],
       'OST_DATA_LENGTH_TAKEN': [], 
       'OST_OPERATIVE_MODE': [],
       'OST_MANUAL_GAIN_CONTROL': [],
       'OST_COMPRESSION_SELECTION': [],
       'OST_CLOSED_LOOP_TRACKING': [], 
       'OST_TRACKING_DATA_STORAGE': [],
       'OST_TRACKING_PRE_SUMMING': [], 
       'OST_TRACKING_LOGIC_SELECTION': [],
       'OST_THRESHOLD_LOGIC_SELECTION': [],
       'OST_SAMPLE_NUMBER': [],
       'OST_SPARE2': [],
       'OST_ALPHA_BETA': [],
       'OST_REFERENCE_BIT': [],
       'OST_THRESHOLD': [],
       'OST_THRESHOLD_INCREMENT': [],
       'OST_SPARE3': [],
       'OST_INITIAL_ECHO_VALUE': [],
       'OST_EXPECTED_ECHO_SHIFT': [],
       'OST_WINDOW_LEFT_SHIFT': [],
       'OST_WINDOW_RIGHT_SHIFT': [], 
       'OST_SPARE4': [],
       'SPARE3': [],
       'DATA_BLOCK_ID': [],
       'SCIENCE_DATA_SOURCE_COUNTER': [],
       'PSFPGA_SCIENTIFIC_DATA_TYPE': [],
       'PSFPGA_SEGMENTATION_FLAG': [],
       'PSFPGA_SPARE1': [],
       'PSFPGA_SPARE2': [],
       'PSFPGA_DMA_ERROR': [],
       'PSFPGA_TC_OVERRUN': [],
       'PSFPGA_FIFO_FULL': [],
       'PSFPGA_TEST': [],
       'SPARE4': [],
       'DATA_BLOCK_FIRST_PRI': [],
       'TIME_DATA_BLOCK_WHOLE': [],
       'TIME_DATA_BLOCK_FRAC': [],
       'SDI_BIT_FIELD': [],
       'TIME_N': [],
       'RADIUS_N': [],
       'TANGENTIAL_VELOCITY_N': [],
       'RADIAL_VELOCITY_N': [],
       'TLP': [],
       'TIME_WPF': [], 
       'DELTA_TIME': [],
       'TLP_INTERPOLATE': [],
       'RADIUS_INTERPOLATE': [],
       'TANGENTIAL_VELOCITY_INTERPOLATE': [], 
       'RADIAL_VELOCITY_INTERPOLATE': [],
       'END_TLP': [],
       'S_COEFFS': [],
       'C_COEFFS': [],
       'SLOPE': [],
       'TOPOGRAPHY': [],
       'PHASE_COMPENSATION_STEP': [],
       'RECEIVE_WINDOW_OPENING_TIME': [],
       'RECEIVE_WINDOW_POSITION': [],
     }
  #
  # Ancillary rows are 186 bytes long
  #
  recLen = 186
  if os.path.isfile(fname):
    _f = open(fname, 'rb')
    fsize = os.path.getsize(fname)
    for _i in range(int(fsize/recLen)): # Go through all the rows
      _f.seek(_i*recLen)
      data = _f.read(recLen)
      #
      # Set up dictionary for items
      #
      a['SCET_BLOCK_WHOLE'].append(struct.unpack('>I', data[0:4])[0])
      a['SCET_BLOCK_FRAC'].append(struct.unpack('>H', data[4:6])[0])
      a['TLM_COUNTER'].append(struct.unpack('>I', data[6:10])[0])
      a['FMT_LENGTH'].append(struct.unpack('>H', data[10:12])[0])
      a['SPARE1'].append(struct.unpack('>H', data[12:14])[0])
      a['SCET_OST_WHOLE'].append(struct.unpack('>I', data[14:18])[0])
      a['SCET_OST_FRAC'].append(struct.unpack('>H', data[18:20])[0])
      a['SPARE2'].append(struct.unpack('>B', data[20:21])[0])
      a['OST_LINE_NUMBER'].append(struct.unpack('>B', data[21:22])[0])
      #
      # Deal with the OST_LINE Entries
      #
      OST_LINE = bitstring.BitArray(data[22:39])
      tmp = OST_LINE[0:4].uint
      if tmp == 1:
        a['OST_PULSE_REPETITION_INTERVAL'].append(1428)
      elif tmp == 2:
        a['OST_PULSE_REPETITION_INTERVAL'].append(1492)
      elif tmp == 3:
        a['OST_PULSE_REPETITION_INTERVAL'].append(1290)
      elif tmp == 4:
        a['OST_PULSE_REPETITION_INTERVAL'].append(2856)
      elif tmp == 5:
        a['OST_PULSE_REPETITION_INTERVAL'].append(2984)
      else:
        a['OST_PULSE_REPETITION_INTERVAL'].append(2580)
      a['OST_PHASE_COMPENSATION_TYPE'].append(OST_LINE[4:8].uint)
      a['OST_SPARE1'].append(OST_LINE[8:10].uint)
      a['OST_DATA_LENGTH_TAKEN'].append(OST_LINE[10:32].uint)
      a['OST_OPERATIVE_MODE'].append(OST_LINE[32:40].uint)
      a['OST_MANUAL_GAIN_CONTROL'].append(OST_LINE[40:48].uint)
      a['OST_COMPRESSION_SELECTION'].append(OST_LINE[48:49].bool)
      a['OST_CLOSED_LOOP_TRACKING'].append(OST_LINE[49:50].bool)
      a['OST_TRACKING_DATA_STORAGE'].append(OST_LINE[50:51].bool)
      tmp = OST_LINE[51:54].uint
      if tmp == 0:
        a['OST_TRACKING_PRE_SUMMING'].append(0)
      elif tmp == 1:
        a['OST_TRACKING_PRE_SUMMING'].append(2)
      elif tmp == 2:
        a['OST_TRACKING_PRE_SUMMING'].append(3)
      elif tmp == 3:
        a['OST_TRACKING_PRE_SUMMING'].append(4)
      elif tmp == 4:
        a['OST_TRACKING_PRE_SUMMING'].append(8)
      elif tmp == 5:
        a['OST_TRACKING_PRE_SUMMING'].append(16)
      elif tmp == 6:
        a['OST_TRACKING_PRE_SUMMING'].append(32)
      else:
        a['OST_TRACKING_PRE_SUMMING'].append(64)
      a['OST_TRACKING_LOGIC_SELECTION'].append(OST_LINE[54:55].uint)
      a['OST_THRESHOLD_LOGIC_SELECTION'].append(OST_LINE[55:56].uint)
      a['OST_SAMPLE_NUMBER'].append(OST_LINE[56:60].uint)
      a['OST_SPARE2'].append(OST_LINE[60:61].uint)
      a['OST_ALPHA_BETA'].append(OST_LINE[61:63].uint)
      a['OST_REFERENCE_BIT'].append(OST_LINE[63:64].uint)
      a['OST_THRESHOLD'].append(OST_LINE[64:72].uint)
      a['OST_THRESHOLD_INCREMENT'].append(OST_LINE[72:80].uint)
      a['OST_SPARE3'].append(OST_LINE[80:84].uint)
      a['OST_INITIAL_ECHO_VALUE'].append(OST_LINE[84:87].uint)
      a['OST_EXPECTED_ECHO_SHIFT'].append(OST_LINE[87:90].uint)
      a['OST_WINDOW_LEFT_SHIFT'].append(OST_LINE[90:93].uint)
      a['OST_WINDOW_RIGHT_SHIFT'].append(OST_LINE[93:96].uint)
      a['OST_SPARE4'].append(OST_LINE[96:128].uint)
      #
      # End OST LINE
      #
      a['SPARE3'].append(struct.unpack('>B', data[38:39])[0])
      a['DATA_BLOCK_ID'].append(bitstring.BitArray(data[39:42]).uint)
      a['SCIENCE_DATA_SOURCE_COUNTER'].append(struct.unpack('>H', data[42:44])[0])
      #
      # PACKET_SEGMENTATION_AND_FPGA_STATUS bit string
      #
      PSAFS = bitstring.BitArray(data[44:46])
      a['PSFPGA_SCIENTIFIC_DATA_TYPE'].append(PSAFS[0:1].uint)
      a['PSFPGA_SEGMENTATION_FLAG'].append(PSAFS[1:3].uint)
      a['PSFPGA_SPARE1'].append(PSAFS[3:8].uint)
      a['PSFPGA_SPARE2'].append(PSAFS[8:12].uint)
      a['PSFPGA_DMA_ERROR'].append(PSAFS[12:13].uint)
      a['PSFPGA_TC_OVERRUN'].append(PSAFS[13:14].uint)
      a['PSFPGA_FIFO_FULL'].append(PSAFS[14:15].uint)
      a['PSFPGA_TEST'].append(PSAFS[15:16].uint)
      a['SPARE4'].append(struct.unpack('>B', data[46:47])[0])
      a['DATA_BLOCK_FIRST_PRI'].append(bitstring.BitArray(data[47:50]).uint)
      a['TIME_DATA_BLOCK_WHOLE'].append(struct.unpack('>I', data[50:54])[0])
      a['TIME_DATA_BLOCK_FRAC'].append(struct.unpack('>H', data[54:56])[0])
      a['SDI_BIT_FIELD'].append(struct.unpack('>H', data[56:58])[0])
      a['TIME_N'].append(struct.unpack('>f', data[58:62])[0])
      a['RADIUS_N'].append(struct.unpack('>f', data[62:66])[0])
      a['TANGENTIAL_VELOCITY_N'].append(struct.unpack('>f', data[66:70])[0])
      a['RADIAL_VELOCITY_N'].append(struct.unpack('>f', data[70:74])[0])
      a['TLP'].append(struct.unpack('>f', data[74:78])[0])
      a['TIME_WPF'].append(struct.unpack('>f', data[78:82])[0])
      a['DELTA_TIME'].append(struct.unpack('>f', data[82:86])[0])
      a['TLP_INTERPOLATE'].append(struct.unpack('>f', data[86:90])[0])
      a['RADIUS_INTERPOLATE'].append(struct.unpack('>f', data[90:94])[0])
      a['TANGENTIAL_VELOCITY_INTERPOLATE'].append(struct.unpack('>f', data[94:98])[0])
      a['RADIAL_VELOCITY_INTERPOLATE'].append(struct.unpack('>f', data[98:102])[0])
      a['END_TLP'].append(struct.unpack('>f', data[102:106])[0])
      tmp = [struct.unpack('>f', data[106:110])[0],
                            struct.unpack('>f', data[110:114])[0],
                            struct.unpack('>f', data[114:118])[0],
                            struct.unpack('>f', data[118:122])[0],
                            struct.unpack('>f', data[122:126])[0],
                            struct.unpack('>f', data[126:130])[0],
                            struct.unpack('>f', data[130:134])[0],
                            struct.unpack('>f', data[134:138])[0]
            ]
      tmp = [str(i) for i in tmp]
      a['S_COEFFS'].append('"'+ ', '.join(tmp) + '"')
      tmp = [struct.unpack('>f', data[138:142])[0],
                            struct.unpack('>f', data[142:146])[0],
                            struct.unpack('>f', data[146:150])[0],
                            struct.unpack('>f', data[150:154])[0],
                            struct.unpack('>f', data[154:158])[0],
                            struct.unpack('>f', data[158:162])[0],
                            struct.unpack('>f', data[162:166])[0]
            ]
      tmp = [str(i) for i in tmp]
      a['C_COEFFS'].append('"'+ ', '.join(tmp) + '"')
      a['SLOPE'].append(struct.unpack('>f', data[166:170])[0])
      a['TOPOGRAPHY'].append(struct.unpack('>f', data[170:174])[0])
      a['PHASE_COMPENSATION_STEP'].append(struct.unpack('>f', data[174:178])[0])
      a['RECEIVE_WINDOW_OPENING_TIME'].append(struct.unpack('>f', data[178:182])[0])
      a['RECEIVE_WINDOW_POSITION'].append(struct.unpack('>f', data[182:186])[0])
    return a


def sepSAdata(iS, oS, oA, n, b, p):
  _a = open("tmp.anc", 'wb')
  _s = open("tmp.sci", 'wb')
  #
  # Now split the SCIENCE and ANCILIARY DATA
  #
  cnt = -1
  with open(iS, 'rb') as f:
    while True:
      s = f.read(n)
      cnt += 1
      if not s:
        # EOF
        break
      _a.write(s[:186])
      _s.write(s[186:])
  _a.close()
  _s.close()
  #
  # Save Ancillary data as CSV
  #
  a = parseAncillary("tmp.anc")  
  anc = pd.DataFrame.from_dict(a)
  anc.to_csv(oA)
  #
  # Deal with the science data, the goal is to save each science file are 8-bit
  # signed integers
  #
  if b == 8:
    tmp = np.fromfile('tmp.sci', dtype='int8') 
  elif b == 6:
    _f = open('tmp.sci', 'rb')
    bs = bitstring.BitArray(_f)
    tmp = np.zeros(int(len(bs)/b), int)
    for _j in range(0, len(bs), b):
      tmp[int(_j/b)] = bs[_j:_j+b].int
  elif b == 4:
    tmp = np.fromfile('tmp.sci', dtype='int4')
  #
  # TMP now holds the data, now decompress data...
  #
  decom = getDecom(a['OST_COMPRESSION_SELECTION'][0], p, b, a['SDI_BIT_FIELD'])
  data = tmp * decom
  #
  # Now save decompressed science data as int8 binary
  #
  data.astype('int8').tofile(oS)
  #
  #
  # Remove temporary files
  #
  os.remove("tmp.anc")
  os.remove("tmp.sci")
  return cnt

def detPRF(val):
  PRF = {'335': 335.12,
         '350': 350.14,
         '387': 387.60,
         '670': 670.24,
         '700': 700.28,
         '775': 775.19,
         }
  return PRF[str(val)] 


def getDecom(c, p, b, s):
  """
    This function decompresses the data based off page 8 of the
    SHALLOW RADAR EXPERIMENT DATA RECORD SOFTWARE INTERFACE SPECIFICATION.
    If the compression type is 'Static'
       U = C* ( 2**S / N )
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
  #
  # Step 5: Decompress the data
  # Note: Only Static decompression works at this time
  #
  if c is not True:
    c = 'STATIC'
  else:
    c = 'DYNAMIC'
  if c == 'STATIC' or c == 'DYNAMIC':
    #
    # Handle fixed scaling
    #
    if c == 'STATIC': # Static scaling
      L = np.ceil(np.log2(int(p)))
      R = b
      S = L - R + 8
      N = p
      return np.power(2, S) / N
    elif c == True:#dynamic scaling
      N = p
      if SDI <= 5:
        S = SDI
      elif 5 < SDI <= 16:
        S = SDI - 6
      elif SDI > 16:
        S = SDI - 16
      return np.power(2, S) / N
    return

def formOut(lblFile, outDir):
  oDirs = {'SCIENCE': outDir + '/SCIENCE/',
           'ANCILLARY': outDir + '/ANCILLARY/',
           'AUX': outDir + '/AUXILIARY/',
           'LOGS': outDir + '/LOGS/',
          }
  tmp = os.path.basename(lblFile).split('.')[0]
  oFiles = {'SCIENCE': oDirs['SCIENCE'] + tmp + '.sci',
            'ANCILLARY': oDirs['ANCILLARY'] + tmp +'.csv',
            'AUX': oDirs['AUX'] + tmp + '.csv',
            'LOG': oDirs['LOGS'] + tmp + '.log',
           }
  return oDirs, oFiles

def makeOut(outDir, oDirs):
  #
  # Make output directories
  #
  if not os.path.isdir(outDir):
    print('[WARNING] New output directory will be created at {}'.format(outDir))
    os.mkdir(outDir)
  #
  # Check output subdirectories
  #
  if not os.path.isdir(oDirs['SCIENCE']):
    os.mkdir(oDirs['SCIENCE'])
  if not os.path.isdir(oDirs['ANCILLARY']):
    os.mkdir(oDirs['ANCILLARY'])
  if not os.path.isdir(oDirs['AUX']):
    os.mkdir(oDirs['AUX'])
  if not os.path.isdir(oDirs['LOGS']):
    os.mkdir(oDirs['LOGS'])
  return
