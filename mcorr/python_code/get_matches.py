#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 17 13:58:18 2021

@author: piccinini
"""

# legge dal file lista.matched
# ed estrae N secondi di traccia dalla stazione STAZ

from obspy.clients.fdsn import Client
from obspy import UTCDateTime
client = Client("INGV")


# ESTRAGGO IL CATALOGO DA VERNA
#client = Client('http://verna.pi.ingv.it:8080')



STA='AQU'

file1 = open('lista.template', 'r') 
Lines = file1.readlines() 
count = 0
# Strips the newline character 
for line in Lines: 
    count=count+1;
    print(count)
    t = UTCDateTime(line.strip())
    s=t-2;
    e=t+6;
    NM=line.strip()
    try:
      S=client.get_waveforms('MN',STA,'','HHZ',s, e)
      S.detrend()
      S.taper(max_percentage=0.01)
      S.filter('bandpass',freqmin=4,freqmax=20)
      S.trim(starttime=s,endtime=e)
      S[0].write(NM+'.' + STA + '.Z.sac',format='SAC')
    except:
      print('No data found')