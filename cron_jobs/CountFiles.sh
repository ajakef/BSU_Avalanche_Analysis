#!/bin/bash
ls -lht /home/ccube-admin/CCUBE_test/ | head -n 40 | grep $1 | wc -l
#ls -lht /home/ccube-admin/CCUBE_mseed_pre-3-24/ | head -n 20 | grep $1 | wc -l
