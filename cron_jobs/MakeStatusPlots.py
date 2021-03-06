import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
from datetime import datetime
import matplotlib.dates as mdates
import sys
sys.path.append('/home/ccube-admin/code/lib')
import SendEmail

# , names = ['t', 'Batt125', 'Backup125', 'Temp125', 'tun0125', 'wwan0125','Batt129', 'Backup129', 'Temp129', 'tun0129', 'wwan0129','Batt147', 'Backup147', 'Temp147', 'tun0147', 'wwan0147']
allData = pd.read_csv('/home/ccube-admin/logfiles/CCUBE_data.txt', sep = ' ')

t = []
for i in range(len(allData.values[:,0])):
    t += [datetime.strptime(allData.values[i,0], '%Y-%m-%dT%H:%M:%S')]

t = np.array(t)

CC125 = allData.values[:,1:6]
CC129 = allData.values[:,6:11]
CC147 = allData.values[:,11:16]

plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d %H:%M'))
plt.gca().xaxis.set_major_locator(mdates.HourLocator(interval=12))
plt.plot(t, CC129[:,0])
plt.plot(t, CC147[:,0])
plt.plot(t, CC125[:,0])
plt.gcf().autofmt_xdate()
plt.legend(['LCC0: CC129, ADA', 'LCC1: CC147, BE4', 'LCC2: CC125, AD9'])
plt.ylabel('Battery Voltage (V)')
plt.ylim([10.9, 13.3])
filename = '/home/ccube-admin/status_plots/' + datetime.now().strftime('%Y-%m-%dT%H:%M:%S') + '_BatteryVoltage.png'
plt.savefig(filename)

text = """
<html>
<head>
</head>
<body>
<p><img src="cid:batt"></p>
</body>
</html>
"""
SendEmail.SendEmail('jacobanderson152@boisestate.edu,jeffreybjohnson@boisestate.edu', 'LCC CCUBE status', text, filename, 'batt')
