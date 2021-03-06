#!/usr/bin/python
from email.mime.image import MIMEImage
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
import time, os, datetime, smtplib

def SendEmail(to, subject, text, attachmentFiles, attachmentNames):
    ## make sure that attachment info is list, not string
    if type(attachmentFiles) == str:
        attachmentFiles = [attachmentFiles]
    if type(attachmentNames) == str:
        attachmentNames = [attachmentNames]

    username = "boise.river.monitor@gmail.com"
    password = "arthurfoote" # DoB: July 4 1963  
    msg = MIMEMultipart()
    msg['Subject'] = subject
    msg['From'] = username
    msg['To'] = to
    msg.add_header('Content-Type', 'text/html')
    msg.attach(MIMEText(text, 'html'))
    
    for (attachmentFile, attachmentName) in zip(attachmentFiles, attachmentNames):
        with open(attachmentFile, 'rb') as a:
            img = MIMEImage(a.read())
            img.add_header('Content-ID', '<'+ attachmentName + '>')
            msg.attach(img)

    server = smtplib.SMTP('smtp.gmail.com:587')
    server.ehlo()
    server.starttls()
    server.login(username,password)
    server.sendmail(msg['From'],msg['To'].split(','),msg.as_string())
    server.quit()



