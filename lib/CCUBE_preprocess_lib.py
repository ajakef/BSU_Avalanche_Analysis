## obspy: need to define handle_data(trace) in advance. for example...
import numpy as np
import subprocess, re
from obspy.clients.seedlink.easyseedlink import create_client
from obspy import Stream
from obspy import UTCDateTime
from obspy.signal.cross_correlation import xcorr_max, xcorr_pick_correction, xcorr
import sys
def FindConnectionIP(CCUBE_ID):
    with open(str(CCUBE_ID) + '_status.txt', 'r') as file:
        ip_addr_info = file.read()
    tun0_start = ip_addr_info.find('tun0')
    inet_start = tun0_start + ip_addr_info[tun0_start:].find('inet')
    section_info = ip_addr_info[inet_start:(100+inet_start)] 
    IP = re.findall(r'(?:\d{1,3}\.)+(?:\d{1,3}\.)+(?:\d{1,3}\.)+(?:\d{1,3})', section_info)[0]
    return IP


CCUBE_ID = sys.argv[1]
if len(CCUBE_ID) == 3:
    CCUBE_ID = 'CC' + CCUBE_ID
elif len(CCUBE_ID) != 5:
    print("Bad CCUBE ID")
    exit()

DATACUBE_ID = sys.argv[2]
CCUBE_IP = FindConnectionIP(CCUBE_ID)

print(CCUBE_IP)
print(DATACUBE_ID)
print(CCUBE_ID)
analysis_dt = 10 # seconds
archive_hour = False # if false, use minute-long files


logfile = open('log.txt','a')
eventfile = open('events.txt', 'a')
consistency_threshold = 3.1 # samples

def fmt_int(n, d=2):
    return '{:0>2d}'.format(n)

def fmt(n):
    return '{:0.3f}'.format(n)

def writeSplitMSEED(stream):
    for mseedTrace in stream.split():
        traceStart = mseedTrace.stats.starttime
        filename = 'CCUBE_test/' + fmt_int(traceStart.year) + \
            '.' + fmt_int(traceStart.julday, 3) + \
            '.' + fmt_int(traceStart.hour) + '.' + \
            fmt_int(traceStart.minute) + '.' + \
            fmt_int(traceStart.second) + '.' + \
            mseedTrace.stats.station + '.' + \
            mseedTrace.stats.channel + '.mseed'
        mseedTrace.write(filename, format = 'MSEED')
        print("Writing " + filename)

def handle_data(incomingData):
    global allData, mseedData, current_archive_hour, new_archive_hour, archive_dt
    global analysisData, current_analysis_t0, new_analysis_t0, analysis_dt
    global t0, t1, tc
    global logfile, eventfile
    print(incomingData)
    allData.append(incomingData)
    if(t0 > (new_archive_hour + 10)):
        mseedData = allData.copy()
        mseedData.trim(current_archive_hour, new_archive_hour - 0.001)
        mseedData.merge()
        if(False):
            filename = 'CCUBE_test/' + fmt_int(current_archive_hour.year) + \
                '.' + fmt_int(current_archive_hour.julday) + \
                '.' + fmt_int(current_archive_hour.hour) + '.' + \
                fmt_int(current_archive_hour.minute) + '.' + \
                fmt_int(current_archive_hour.second) + '.' + \
                mseedData.traces[0].stats.station + '.mseed'

        writeSplitMSEED(mseedData)

        allData.trim(new_archive_hour, new_archive_hour+archive_dt+86400)
        current_archive_hour = current_archive_hour + archive_dt
        new_archive_hour = new_archive_hour + archive_dt
        ################
        
    while(UTCDateTime() > (new_analysis_t0 + 6*analysis_dt)):
        print("Trying to analyze " + str(current_analysis_t0))
        t0 = current_analysis_t0 - 2.0 * analysis_dt
        t1 = new_analysis_t0 + 2.0 * analysis_dt
        tc = current_analysis_t0 + 0.5 * analysis_dt
        analysisData = allData.copy()
        analysisData.trim(t0, t1) # wide to avoid filter artifacts
        analysisData.merge()
        analysisData = analysisData.split() # for some reason split doesn't save to analysisData, but merge does
        if(len(analysisData.traces) < 3):
            print('Not enough traces to analyze: ' + str(len(analysisData.traces)))
        elif(len(analysisData.traces) > 3):
            print('Too many traces: ' + str(len(analysisData.traces)) + ' (could mean missing data)')
        else:
            if(not(analysisData.traces[0].meta.starttime <= (t0 + 0.1) and analysisData.traces[0].meta.endtime >= (t1-0.1) and \
               analysisData.traces[1].meta.starttime <= (t0 + 0.1) and analysisData.traces[1].meta.endtime >= (t1-0.1) and \
               analysisData.traces[2].meta.starttime <= (t0 + 0.1) and analysisData.traces[2].meta.endtime >= (t1-0.1))):
                print('Not enough data to analyze')
            else:
                print('Analyzing')
                analysisData.filter('highpass', freq=2)
                analysisData.trim(current_analysis_t0, new_analysis_t0)
                [lag12, r12] = xcorr_max(xcorr(analysisData.traces[0], analysisData.traces[1], 20, full_xcorr=True)[2], abs_max = False)
                [lag23, r23] = xcorr_max(xcorr(analysisData.traces[1], analysisData.traces[2], 20, full_xcorr=True)[2], abs_max = False)
                [lag31, r31] = xcorr_max(xcorr(analysisData.traces[2], analysisData.traces[0], 20, full_xcorr=True)[2], abs_max = False)
                consistency = lag12 + lag23 + lag31
                med = np.median([r12, r23, r31])
                rms = np.median([analysisData.traces[j].std() for j in np.arange(3)])
                logfile.write(str(current_analysis_t0) + ',' + fmt(consistency * 0.01) + ',' + fmt(med) + ',' + fmt(rms) + ',' + \
                      fmt(lag12*0.01) + ',' + fmt(lag23*0.01) + ',' + fmt(lag31*0.01) + '\n')
                logfile.flush()
                print(str(current_analysis_t0) + ',' + fmt(consistency * 0.01) + ',' + fmt(med) + ',' + fmt(rms) + ',' + \
                      fmt(lag12*0.01) + ',' + fmt(lag23*0.01) + ',' + fmt(lag31*0.01))
                ## write to the events list only if we detect an event
                if(np.abs(consistency) <= consistency_threshold):
                    eventfile.write(str(current_analysis_t0) + ',' + fmt(lag12*0.01) + ',' + fmt(lag23*0.01) + ',' + fmt(lag31*0.01) +'\n')
                    eventfile.flush()
                    print('Event detected! '+ str(current_analysis_t0))
                print(" ")
        current_analysis_t0 = current_analysis_t0 + analysis_dt
        new_analysis_t0 = new_analysis_t0 + analysis_dt
                
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
new_analysis_t0 = current_analysis_t0 + analysis_dt
t0 = current_analysis_t0 - 2.0 * analysis_dt
client = create_client(CCUBE_IP, on_data=handle_data)
client.select_stream('UP', DATACUBE_ID, '???') ## specifying BE4 is mandatory, or it will have an opaque error
client.run()
