import glob, sys, os, struct, bitstring, argparse
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

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
  bname = os.path.basename(_file).split('_')
  TransID = bname[1]
  OSTLine = bname[2]
  OperMode = bname[3].upper()
  PRF = bname[4]
  if OperMode[0:2] == 'SS':
    OperMode = SSInstrMode[OperMode]
  elif OperMode[0:2] == 'RO':
    OperMode = ROInstrMode[OperMode]
  #
  # Determine to Bits per Sample
  #
  if OperMode['BitsPerSample'] == 4:
    recLen = 1986
  elif OperMode['BitsPerSample'] == 6:
    recLen = 2886
  elif OperMode['BitsPerSample'] == 8:
    recLen = 3786
  nrec = int(os.path.getsize(_file) / recLen)
  return TransID, OSTLine, OperMode, PRF, recLen, nrec

def parseAuxFile(fname, outDir='./.', dic=True, df=False, csv=False, binary=False, saveNP=False):
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
      csvName = outDir + '/' + os.path.basename(fname).split('.')[0] + '.csv'
      aux = pd.DataFrame.from_dict(auxData)
      aux.to_csv(csvName)
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
  prog = "SHARAD EDR Processing Algorithm"
  vers = "0.1"
  path, lblFile = os.path.split(sys.argv[1])
  lblFile = path + '/' + lblFile
  outDir = sys.argv[2]
  outSci = outDir + '/SCIENCE/'
  outAnc = outDir + '/ANCILLARY/'
  outAux = outDir + '/AUXILIARY/'
  #
  # Open LBL File and get the SCIENCE_TELEMETRY_TABLE and AUXILIARY_DATA_TABLE
  # entries
  #
  if os.path.isfile(lblFile):
    #
    # File exists
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
  else:
    return 0
  if os.path.isfile(path + '/' + sciFile.upper()):
    sciFile = path + '/' + sciFile.upper()
  elif os.path.isfile(path + '/' + sciFile.lower()):
    sciFile = path + '/' + sciFile.lower()
  else:
    print('Science file not found...')
    return 0
  if os.path.isfile(path + '/' + auxFile.upper()):
    auxFile = path + '/' + auxFile.upper()
  elif os.path.isfile(path + '/' + auxFile.lower()):
    auxFile = path + '/' + auxFile.lower()
  else:
    print('Auxiliary file not found...')
    return 0
  #
  # Parse the PDS file name for important information
  #
  TransID, OSTLine, OperMode, PRF, recLen, nrec = parseFileName(sciFile)
  print('----- {} {} Processing Log -----'.format(prog, vers))
  print('PDS Label file:\t{}'.format(lblFile))
  print('Auxiliary file:\t{}'.format(auxFile))
  print('Science file:\t{}'.format(sciFile))
  print('TransID:\t{}'.format(TransID))
  print('OSTLine:\t{}'.format(OSTLine))
  print('OperMode:\t{}'.format(OperMode['Mode']))
  print('\tOn-board presumming:\t{}'.format(OperMode['Presum']))
  print('\tBits per Sample:\t{}'.format(OperMode['BitsPerSample']))
  print('PRF:\t{}'.format(PRF))
  print('Record Length:\t{}'.format(recLen))
  print('Number of Records:\t{}'.format(nrec))
  print('')
  #
  # Make output directories
  #
  if not os.path.isdir(outDir):
    print('[WARNING] New output directory will be created at {}'.format(outDir))
    os.mkdir(outDir) 
  #
  # Check output subdirectories
  #
  if not os.path.isdir(outSci):
    os.mkdir(outSci)
#  if not os.path.isdir(outAnc):
#    os.mkdir(outAnc)
  if not os.path.isdir(outAux):
    os.mkdir(outAux)
  #
  # Parse Auxiliary File
  # 
  parseAuxFile(auxFile, outDir=outAux, dic=False, csv=True)
  #
  # Now deal with the SCIENCE and ANCILIARY DATA
  #
  echoes = np.zeros([3600, nrec])
  _sciFile = open(sciFile, 'rb')
  for _i in range(0, nrec):
    sci, ancil = readEDRrecord(_sciFile, _i, recLen, OperMode['BitsPerSample'])
    #
    # Skipping ancillary for now...the important information is all captured in
    # the auxilliary file
    # 
    if _i == 0:
      ancil = parseAncillary(ancil)
    echoes[:,_i] = decompressSciData(sci,
              ancil['OST_LINE']['COMPRESSION_SELECTION'],
              OperMode['Presum'],
              OperMode['BitsPerSample'],
              ancil['SDI_BIT_FIELD'])
  csvName = outSci + '/' + os.path.basename(sciFile).split('.')[0] + '.csv'
  np.savetxt(csvName, echoes, fmt='%.4f', delimiter=',', newline='\n', header='', footer='', comments='# ', encoding=None)
  _sciFile.close() 
  return

if __name__ == '__main__':
  main()
