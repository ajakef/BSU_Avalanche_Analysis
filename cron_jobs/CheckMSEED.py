import subprocess
import sys
sys.path.append('/home/ccube-admin/code/lib')
import SendEmail

try:
    with open('/home/ccube-admin/code/cron_jobs/prev_problems.txt', 'r') as file:
        prev_problems = file.readlines()
        prev_problems = [x.strip() for x in prev_problems] # remove newlines
except:
    prev_problems = []
    
problems = []
new_problems = []
for SN in ['AD9', 'ADA', 'BE4']:
    n = float(subprocess.check_output(['/home/ccube-admin/code/cron_jobs/CountFiles.sh', SN]))
    if(n > 2):
        problems = problems # do nothing
        #print(SN + ' ok')
    else:
        problems += [SN]
        #print(SN + ' missing')
        if not (SN in prev_problems):
            new_problems += SN


with open('/home/ccube-admin/code/cron_jobs/prev_problems.txt', 'w') as file:
    file.write('\n'.join(problems))

if len(new_problems) > 0:
    text = 'apparent dropout: ' + ', '.join(problems)
    #print(text)
    #to = 'jacobanderson152@boisestate.edu,jeffreybjohnson@boisestate.edu'
    to = 'jacobanderson152@boisestate.edu'
    SendEmail.SendEmail(to, 'Dropout report', text, [], [])
