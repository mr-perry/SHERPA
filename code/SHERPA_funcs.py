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
  procType = ['Serial']
  proc = [1]
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
  # Make outDir
  #
  oDirs = {'SCIENCE': outDir + '/SCIENCE/',
           'ANCILLARY': outDir + '/ANCILLARY/', 
           'AUX': outDir + '/AUXILIARY/',
           'LOGS': outDir + '/LOGS/',
          }
  tmp = os.path.basename(lblFile).split('.')[0]
  oFiles = {'SCIENCE': oDirs['SCIENCE'] + tmp + '.sci',
            'ANCILLARY': oDirs['ANCILLARY'] + tmp +'.anc',
            'AUX': oDirs['AUX'] + tmp + '.csv',
            'LOG': oDirs['LOGS'] + tmp + '.log',
           }
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
  PRF = bname[4]
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
    # Check output
    #
    if csv == True:
      aux = pd.DataFrame.from_dict(auxData)
      aux.to_csv(oFile)
      return
    if dic == True:
      return auxData
    elif df == True:
      return pd.DataFrame.from_dict(auxData)

def parseAncillary(data):
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
    a = { 'SCET_BLOCK_WHOLE': bitstring.BitArray(data[0:4]).uint,
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
    a['PACKET_SEGMENTATION_AND_FPGA_STATUS']['SCIENTIFIC_DATA_TYPE'] = PSAFS[0:1].uint
    a['PACKET_SEGMENTATION_AND_FPGA_STATUS']['SEGMENTATION_FLAG'] = PSAFS[1:3].uint
    a['PACKET_SEGMENTATION_AND_FPGA_STATUS']['SPARE1'] = PSAFS[3:8].uint
    a['PACKET_SEGMENTATION_AND_FPGA_STATUS']['SPARE2'] = PSAFS[8:12].uint
    a['PACKET_SEGMENTATION_AND_FPGA_STATUS']['DMA_ERROR'] = PSAFS[12:13].uint
    a['PACKET_SEGMENTATION_AND_FPGA_STATUS']['TC_OVERRUN'] = PSAFS[13:14].uint
    a['PACKET_SEGMENTATION_AND_FPGA_STATUS']['FIFO_FULL'] = PSAFS[14:15].uint
    a['PACKET_SEGMENTATION_AND_FPGA_STATUS']['TEST'] = PSAFS[15:16].uint
    #####################################################################################
    #
    # OST_LINE_NUMBER BIT STRING
    #
    OST_LINE = bitstring.BitArray(data[22:39])
    a['OST_LINE']['PULSE_REPETITION_INTERVAL'] = OST_LINE[0:4].uint
    a['OST_LINE']['PHASE_COMPENSATION_TYPE'] = OST_LINE[4:8].uint
    a['OST_LINE']['SPARE1'] = OST_LINE[8:10].uint
    a['OST_LINE']['DATA_LENGTH_TAKEN'] = OST_LINE[10:32].uint
    a['OST_LINE']['OPERATIVE_MODE'] = OST_LINE[32:40].uint
    a['OST_LINE']['MANUAL_GAIN_CONTROL'] = OST_LINE[40:48].uint
    a['OST_LINE']['COMPRESSION_SELECTION'] = OST_LINE[48:49].bool
    a['OST_LINE']['CLOSED_LOOP_TRACKING'] = OST_LINE[49:50].bool
    a['OST_LINE']['TRACKING_DATA_STORAGE'] = OST_LINE[50:51].bool
    a['OST_LINE']['TRACKING_PRE_SUMMING'] = OST_LINE[51:54].uint
    a['OST_LINE']['TRACKING_LOGIC_SELECTION'] = OST_LINE[54:55].uint
    a['OST_LINE']['THRESHOLD_LOGIC_SELECTION'] = OST_LINE[55:56].uint
    a['OST_LINE']['SAMPLE_NUMBER'] = OST_LINE[56:60].uint
    a['OST_LINE']['SPARE2'] = OST_LINE[60:61].uint
    a['OST_LINE']['ALPHA_BETA'] = OST_LINE[61:63].uint
    a['OST_LINE']['REFERENCE_BIT'] = OST_LINE[63:64].uint
    a['OST_LINE']['THRESHOLD'] = OST_LINE[64:72].uint
    a['OST_LINE']['THRESHOLD_INCREMENT'] = OST_LINE[72:80].uint
    a['OST_LINE']['SPARE3'] = OST_LINE[80:84].uint
    a['OST_LINE']['INITIAL_ECHO_VALUE'] = OST_LINE[84:87].uint
    a['OST_LINE']['EXPECTED_ECHO_SHIFT'] = OST_LINE[87:90].uint
    a['OST_LINE']['WINDOW_LEFT_SHIFT'] = OST_LINE[90:93].uint
    a['OST_LINE']['WINDOW_RIGHT_SHIFT'] = OST_LINE[93:96].uint
    a['OST_LINE']['SPARE4'] = OST_LINE[96:128].uint
    return a


def sepSAdata(iS, oS, oA, n, _log):
  _a = open(oA, 'wb')
  _s = open(oS, 'wb')
  #
  # Now split the SCIENCE and ANCILIARY DATA
  #
  cnt = -1
  with open(iS, 'rb') as f:
    while True:
      cnt += 1
      s = f.read(n)
      if not s:
        # EOF
        break
      _a.write(s[:186])
      _s.write(s[186:])
  _a.close()
  _s.close()
  return cnt
