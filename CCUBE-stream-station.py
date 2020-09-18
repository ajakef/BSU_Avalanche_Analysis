import numpy as np
import subprocess
from obspy.clients.seedlink.easyseedlink import create_client
from obspy import Stream
from obspy import UTCDateTime
from obspy.signal.cross_correlation import xcorr_max, xcorr_pick_correction, xcorr
import sys
sys.path.append('/home/ccube-admin/code/lib')
import CCUBE_preprocess_lib
#from networks import LCC_2020 as network
from networks import aftershocks_2020 as network

#CCUBE_ID = CCUBE_preprocess_lib.ParseID(sys.argv[1])
#DATACUBE_ID = sys.argv[2]
data_dir = '/home/ccube-admin/CCUBE_test'
Station_ID = sys.argv[1]
station_index = np.where(network.Station == Station_ID)[0]
if(len(station_index) > 1):
    print('More than one station in network table matches input station ID')
    exit()
elif(len(station_index) == 0):
    print('No stations in network table match input station ID')
    exit()
else:
    station_index = station_index[0] # make it an int instead of an array

DATACUBE_ID = network.Datacube[station_index]
CCUBE_ID = network.CCUBE[station_index]
CCUBE_IP = CCUBE_preprocess_lib.FindIP(CCUBE_ID)

print(CCUBE_ID + ' ' + DATACUBE_ID + ' ' + CCUBE_IP)
analysis_width = 2 # seconds
analysis_overlap = 0.5
archive_hour = True # if false, use minute-long files
correlogram_length = 60 # number overlapped windows, not upsampled


def handle_data(incomingData):
    global allData, mseedData, current_archive_hour, new_archive_hour, archive_dt
    global analysisData, current_analysis_t0, new_analysis_t0, analysis_width, analysis_overlap
    global t_end, t_max#, t0, t1, tc
    global logfile, eventfile
    global correlogram, cgram_count, cgram_t
    print(incomingData)

    ## add incoming data to the buffer
    allData.append(incomingData)
    allData.merge()
    t_end = min([tr.stats.endtime for tr in allData])
    t_max = max([tr.stats.endtime for tr in allData])
    ## change the station name to LCC0-2 instead of the datacube number
    for tr in allData: 
        tr.stats.station = Station_ID

    mseedData = allData.copy()
    mseedData.trim(current_archive_hour, new_archive_hour - 0.001)
    CCUBE_preprocess_lib.writeSplitMSEED(mseedData, data_dir)

    ## if we've gone over a new file start, write that file's data, then cut it from the buffer
    if(t_end > (new_archive_hour + 10)):
        allData.trim(new_archive_hour, t_max + 1) 
        current_archive_hour = current_archive_hour + archive_dt
        new_archive_hour = new_archive_hour + archive_dt
        ################

    ## check to see if we have enough data to analyze
    #print(t_end)
    while(t_end > (current_analysis_t0 + 2*analysis_width)):
        print("Trying to analyze " + str(current_analysis_t0))

        ## calculate the basic cross-correlation and signal stats
        analysisData = allData.copy()
        ccfInfo = CCUBE_preprocess_lib.CheckAnalyzeIncoming(analysisData, current_analysis_t0, current_analysis_t0 + analysis_width)

        ## save the CCFs and increment the count
        cgram_t[cgram_count] = current_analysis_t0
        #print(ccfInfo['xc'])
        stats[cgram_count,:] = np.array([ccfInfo['rms']])
        for i in range(3):
            correlogram[i][cgram_count,:] = ccfInfo['xc'][i]
        cgram_count += 1

        ## if we've hit our quota, write the correlograms to text files and reset the count
        if(cgram_count >= correlogram_length):
            cgram_count = 0
            CCUBE_preprocess_lib.WriteCorrelogram(correlogram, stats, cgram_t[0], Station_ID)
            
        ## update the analysis start time
        current_analysis_t0 = current_analysis_t0 + analysis_width * analysis_overlap
                

## create global variables        
allData = Stream()
mseedData = Stream()
analysisData = Stream()
current_archive_hour = UTCDateTime()
if(archive_hour): ## new file every hour
    current_archive_hour.minute = 0
    archive_dt = 3600
else: ## new file every minute
    archive_dt = 60
current_archive_hour.second = 0
current_archive_hour.microsecond = 0
new_archive_hour = current_archive_hour + archive_dt
current_analysis_t0 = current_archive_hour
new_analysis_t0 = current_analysis_t0 + analysis_width
t0 = current_analysis_t0 - 2.0 * analysis_width
cgram_count = 0
cgram_t = [0] * correlogram_length
stats = np.zeros((correlogram_length, 1))
baseCgram = np.zeros((correlogram_length, 1 + 2*CCUBE_preprocess_lib.MAX_SHIFT))
correlogram = [baseCgram, baseCgram, baseCgram]

## start streaming
client = create_client(CCUBE_IP, on_data=handle_data)
client.select_stream('UP', DATACUBE_ID, '???') ## specifying the datacube name is mandatory, or it will have an opaque error
client.run()
