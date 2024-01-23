#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 23 14:19:21 2022

@author: piccinini
"""

import glob
from obspy import read_events


# READ IN NLLoc and write out QML files

D=glob.glob('chianti01*..2022*hyp')

for fil in D:
  A=read_events(fil)
  TIM=str(A[0].origins[0].time)[0:-4]
  A.write('QML/' + TIM +'.xml',format='QUAKEML')


# READ IN QML frpm NLLoc and write .ms files for eac station 
# to be used with multicorr

from obspy.clients.fdsn import Client 
from obspy import UTCDateTime

client = Client('INGV')

LIST=glob.glob('*xml')

n=0;
for EVNM in LIST:
  n=n+1;
  NM=EVNM[0:-4]
  print(str(n)+ '/'+ str(len(LIST))+'    ' + NM)
  TORIG=UTCDateTime(NM)
  STIME=TORIG-10
  ETIME=TORIG+50
  EV=read_events(EVNM)
  STALIST=[]
  for pick in EV[0].picks:
    STA=pick.waveform_id.station_code
    STALIST.append(STA)
    NEWLIST=list(set(STALIST))
  for NS in NEWLIST:
    print(NS)
    try:
      INV=client.get_stations(starttime=STIME,endtime=STIME+1,station=NS,channel='?HZ',network='*',level='response') 
      INV=INV.remove(network='IV',station='AQU')
      INV=INV.remove(network='IV',station='RM??' )
      INV=INV.remove(network='*',station='*',channel='UHZ')
      INV=INV.remove(network='*',station='*',channel='BHZ')
      INV=INV.remove(network='*',station='*',channel='LHZ')
      INV=INV.remove(network='*',station='*',channel='VHZ') 
      CHN=INV[0][0].channels[0].code                       
                              
      WVF=client.get_waveforms('*',NS,'*',CHN,STIME,ETIME)
      WVF.merge(fill_value='interpolate')
      WVF.detrend()
      NMOUT=NM+'.'+NS+'.ms'
      WVF.write(NMOUT,format='MSEED')
    except:
      print('NO DATA FOR STATION:' + NS)

  
