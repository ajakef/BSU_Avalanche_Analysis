import numpy as np
import subprocess, re
import sys
import os
from datetime import datetime
import time

## function to scan a CCUBE status file for the LTE (wwan0) and VPN (tun0) IP addresses
def FindIP(text, connection):
    connection_start = text.find(connection)
    inet_start = connection_start + text[connection_start:].find('inet')
    section = text[inet_start:(100+inet_start)]
    IP = re.findall(r'(?:\d{1,3}\.)+(?:\d{1,3}\.)+(?:\d{1,3}\.)+(?:\d{1,3})', section)[0]
    return IP

## function to scan a CCUBE status file for parameters like temperature and battery voltage
def FindParameter(text, parameter):
    start = text.find(parameter)
    section = text[start:(30 + start)]
    param = re.findall(r'(?:\d{1,3}\.)+(?:\d{1,3})', section)[0]
    return param

## begin the code to run
CCUBE_list = ['CC125', 'CC129', 'CC147'] # list of CCUBE IDs to search

time.sleep(120) # it seems to take less than a minute for status files to post. wait 2 minutes to be safe.

voltageFile = open('/home/ccube-admin/logfiles/CCUBE_voltage.txt', 'a')
allDataFile = open('/home/ccube-admin/logfiles/CCUBE_data.txt', 'a')

t_now = datetime.now().strftime('%Y-%m-%dT%H:%M:%S')
voltageFile.write(t_now + ' ')
allDataFile.write(t_now + ' ')
for CCUBE_ID in CCUBE_list:
    filename = '/home/ccube-admin/' + CCUBE_ID + '_status.txt'

    ## Try to determine the age of the status file
    try:
        file_age = (datetime.now() - datetime.fromtimestamp(os.stat(filename).st_mtime)).seconds 
    except:
        file_age = np.nan

    ## if file cannot be opened or is old, set everything NaN. Otherwise, read the parameters from the file.
    if np.isnan(file_age) | (file_age > 600):
        IP_tun0 = 'NaN'
        IP_wwan0 = 'NaN'
        temperature = 'NaN'
        systemVoltage = 'NaN'
        backupVoltage = 'NaN'
    else:
        with open(filename, 'r') as file:
            text = file.read()
        IP_tun0 = FindIP(text, 'tun0')
        IP_wwan0 = FindIP(text, 'wwan0')
        temperature = FindParameter(text, 'Temperature')
        systemVoltage = FindParameter(text, 'System')
        backupVoltage = FindParameter(text, 'Battery')

    ## optionally print data to the screen
    if(False):    
        print(CCUBE_ID + ',' +
              systemVoltage + ',' +
              backupVoltage + ',' +
              temperature + ',' +
              IP_tun0)
        
    ## write data to the output files    
    voltageFile.write(systemVoltage + ' ')
    allDataFile.write(systemVoltage + ' ' +
                  backupVoltage + ' ' +
                  temperature + ' ' +
                  IP_tun0 + ' ' +
                  IP_wwan0 + ' ')
voltageFile.write('\n')
allDataFile.write('\n')
    

