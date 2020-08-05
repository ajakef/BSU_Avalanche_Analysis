import re
import numpy as np
import pandas as pd
from obspy.clients.seedlink.easyseedlink import create_client
from obspy import Stream
from obspy import UTCDateTime
from obspy.signal.cross_correlation import xcorr_max, xcorr_pick_correction, xcorr
import obspy
from scipy.interpolate import interp1d

MAX_SHIFT = 40
CGRAM_DIR = '/home/ccube-admin/correlograms'
T_FMT = '%Y-%m-%dT%H:%M:%S'
UPSAMPLE = 5

def FindIP(CCUBE_ID):
    with open('/home/ccube-admin/' + str(CCUBE_ID) + '_status.txt', 'r') as file:
        ip_addr_info = file.read()
    tun0_start = ip_addr_info.find('tun0')
    inet_start = tun0_start + ip_addr_info[tun0_start:].find('inet')
    section_info = ip_addr_info[inet_start:(100+inet_start)] 
    IP = re.findall(r'(?:\d{1,3}\.)+(?:\d{1,3}\.)+(?:\d{1,3}\.)+(?:\d{1,3})', section_info)[0]
    return IP

def dig2(n):
    return '{:0>2d}'.format(n)

def fmt(n):
    return '{:0.3f}'.format(n)

def writeSplitMSEED(stream):
    stream.merge()
    for mseedTrace in stream.split():
        traceStart = mseedTrace.stats.starttime
        filename = 'CCUBE_test/' + dig2(traceStart.year) + \
            '.' + dig2(traceStart.julday) + \
            '.' + dig2(traceStart.hour) + '.' + \
            dig2(traceStart.minute) + '.' + \
            dig2(traceStart.second) + '.' + \
            mseedTrace.stats.station + '.' + \
            mseedTrace.stats.channel + '.mseed'
        mseedTrace.write(filename, format = 'MSEED')
        print("Writing " + filename)


def ParseID(CCUBE_ID):
    if len(CCUBE_ID) == 3:
        CCUBE_ID = 'CC' + CCUBE_ID
    elif len(CCUBE_ID) != 5:
        print("Bad CCUBE ID")
        exit()
    return CCUBE_ID

def CheckAnalyzeIncoming(analysisData, t0, t1):
    ## handle missing data
    analysisData.merge()
    analysisData = analysisData.split() # for some reason split doesn't save to analysisData, but merge does

    ## look for problems with the data, and skip analysis if any problems are found
    if(len(analysisData.traces) < 3):
        print('Not enough traces to analyze: ' + str(len(analysisData.traces)))
        ccfInfo = MakeEmptyCcfInfo()
    elif(len(analysisData.traces) > 3):
        print('Too many traces: ' + str(len(analysisData.traces)) + ' (could mean missing data)')
        ccfInfo = MakeEmptyCcfInfo()
    else:
        if(not(analysisData.traces[0].meta.starttime <= (t0 + 0.1) and analysisData.traces[0].meta.endtime >= (t1-0.1) and \
               analysisData.traces[1].meta.starttime <= (t0 + 0.1) and analysisData.traces[1].meta.endtime >= (t1-0.1) and \
               analysisData.traces[2].meta.starttime <= (t0 + 0.1) and analysisData.traces[2].meta.endtime >= (t1-0.1))):
            print('Not enough data to analyze')
            ccfInfo = MakeEmptyCcfInfo()
        else:
            ccfInfo = AnalyzeIncoming(analysisData, t0, t1)
    return ccfInfo

def MakeEmptyCcfInfo():
    from numpy import NaN
    xc = np.zeros(1 + 2*MAX_SHIFT) + NaN # global MAX_SHIFT
    return {'raw_consistency': NaN,
            'raw_lag': [NaN, NaN, NaN],
            'raw_r': [NaN, NaN, NaN],
            'rms': NaN,
            'xc': [xc, xc, xc]}
    
def CustomXCorr(stream, n1, n2, shift = MAX_SHIFT): # global MAX_SHIFT
    ## shift is in samples. anything>40 makes a warning message ("window is too small")
    return xcorr(stream.traces[n1], stream.traces[n2], shift, full_xcorr=True)[2]

def AnalyzeIncoming(analysisData, current_analysis_t0, new_analysis_t0):
    consistency_threshold = 3.1 # samples
    sensitivity = (2.048/64/2**23) / 46e-6
    logfile = open('log.txt','a')
    eventfile = open('events.txt', 'a')
    logfile.close()
    eventfile.close()

    #analysisData = obspy.read()

    ## pre-process the data
    analysisData.trim(current_analysis_t0-0.005, new_analysis_t0+0.005)
    analysisData.detrend('linear')
    analysisData.taper(0.05) # Hann window
    print('Analyzing')

    xc12 = CustomXCorr(analysisData, 0, 1)
    xc23 = CustomXCorr(analysisData, 1, 2)
    xc31 = CustomXCorr(analysisData, 2, 0)
    
    [lag12, r12] = xcorr_max(xc12, abs_max = False)
    [lag23, r23] = xcorr_max(xc23, abs_max = False)
    [lag31, r31] = xcorr_max(xc31, abs_max = False)
    consistency = lag12 + lag23 + lag31
    med = np.median([r12, r23, r31])
    rms = sensitivity * np.median([analysisData.traces[j].std() for j in np.arange(3)]) # can't multiply actual traces by float; do it here instead

    D = {'raw_consistency': consistency,
         'raw_lag': [lag12, lag23, lag31],
         'raw_r': [r12, r23, r31],
         'rms': rms,
         'xc': [xc12, xc23, xc31]}
    return D

def WriteCorrelogram(correlogram, stats, t, Station_ID):
    t = obspy.UTCDateTime(round(float(t))) ## round down OR up to the nearest second. note that this may cause problems if time zones are involved.
    for i in range(3):
        filename = CGRAM_DIR +'/' + t.strftime(T_FMT) + '_' + Station_ID + '_' + str(i) + '.txt' #str(i+1) + str(1 + ((i+1) % 3))
        with open(filename,'wb') as file:
            for line in np.matrix(correlogram[i]):
                np.savetxt(file, line, fmt='%.4f')
    filename = CGRAM_DIR + '/' + t.strftime(T_FMT) + '_' + Station_ID + '_stats.txt' #str(i+1) + str(1 + ((i+1) % 3))
    with open(filename,'wb') as file:
        for line in np.matrix(stats):
            np.savetxt(file, line, fmt='%.4f')

def ProcessCorrelogram(cgram):
    cgram = np.array(cgram)
    ## resample the lag axis
    old_indices = np.arange(cgram.shape[1])
    new_indices = np.arange(1 + (cgram.shape[1]-1) * UPSAMPLE)/UPSAMPLE
    interp_function = scipy.interpolate.interp1d(old_indices, cgram, 'cubic', axis = 1)
    cgram_resamp = interp_function(new_indices)
    ## apply a gaussian filter to the image
    cgram_filt = scipy.ndimage.gaussian_filter(cgram_resamp, [1/0.5, UPSAMPLE*0.1/0.01])
    ## calculate the lag and r of the cross-correlation peak at each time
    n = cgram_filt.shape[0]
    lags = np.zeros(n)
    r = np.zeros(n)
    for i in range(n):
        xc_max = obspy.signal.cross_correlation.xcorr_max(x[i,:])
        lags[i] = xc_max[0]
        r[i] = xc_max[1]
    D = {'lag': lags, 'r': r}
    return D


def ReadCorrelogram(t):
    lags = []
    stats = pd.read_csv(CGRAM_DIR + '/' +  t.strftime(T_FMT) + '_stats.txt', sep = ' ')
    for i in range(3):
        cgram = pd.read_csv(CGRAM_DIR + '/' + t.strftime(T_FMT) + '_' + str(i) + '.txt', sep = ' ')
        xc_peaks = ProcessCorrelogram(cgram)
        lags.append(xc_peaks['lag'])

    consistency = lags[0] + lags[1] + lags[2]
    








#print("imported")
#import obspy
#analysisData = obspy.read()
#current_analysis_t0 = obspy.UTCDateTime(2009,8,24,0,20,4)
#new_analysis_t0 = current_analysis_t0 + 1
#CCUBE_seedlink.AnalyzeIncoming(analysisData, current_analysis_t0, new_analysis_t0)
