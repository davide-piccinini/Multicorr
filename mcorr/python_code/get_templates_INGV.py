#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 17 18:22:42 2021

@author: davide
"""

# modified 6/05/2022
# Ora legge anche un archivio locale oltre a scaricare i dati da fdsn

#import os
from obspy.clients.fdsn import Client 
from obspy import UTCDateTime

START="2022-05-24T00:00:00" 
END  ="2022-05-26T00:00:00"
LATC=43.6867
LONC=11.2554 
RAGG=30 #KM dal centro per gli eventi
RSTA=70/111 # km dal centro per le stazioni

DEC=0  # decimation of data: if set to 1 decimate value is 2


RAGG=RAGG/111

client = Client('INGV')
starttime = UTCDateTime(START) 
endtime   = UTCDateTime(END)

from obspy.clients.filesystem.sds import Client as SDClient
sdsdir='/homes/piccinini/MULTICORR/FAKEARC/'
sdsclient=SDClient(sdsdir)

# END OF COMMON PART 

INV=client.get_stations(starttime=START,endtime=END,channel='?HZ',network='*',latitude=LATC,longitude=LONC,maxradius=RSTA,level='response')
 
INV=INV.remove(network='IV',station='AQU')
INV=INV.remove(network='IV',station='RM??' )
INV=INV.remove(network='*',station='*',channel='UHZ')
INV=INV.remove(network='*',station='*',channel='BHZ')
INV=INV.remove(network='*',station='*',channel='LHZ')
INV=INV.remove(network='*',station='*',channel='VHZ')

for NET in INV:
  for STA in NET:
    STAZ=STA.code
    STLA=STA.latitude
    STLO=STA.longitude
    STEL=STA.elevation
    print("%4s %8.5f %8.5f %5.3f" % (STAZ,STLA,STLO,STEL/1000))

cat = client.get_events(starttime=starttime, endtime=endtime,latitude=LATC,longitude=LONC,maxradius=RAGG) 
print('Catalog contains ', len(cat), 'events')
# NOW DOWNLOAD EVENTS WAVEFORMS FOR SELECTED STATIONS
WAV=1 # SWITCH 0/1 per scaricare le WVF o scrivere solo il catalogo
 
n=0
for ev in cat:
    u=str(ev.resource_id) 
    kk=u.find('=') 
    id=u[kk+1:]
    eq=client.get_events(eventid=id,includearrivals=True)
    if eq[0].origins[0].evaluation_mode=='manual' and eq[0].origins[0].evaluation_status=='reviewed': 
        n=n+1 
        print(n)
        TOUT=eq[0].origins[0].time 
        OUT=str(TOUT)[: -4] + '.xml'
        eqk=eq[0]
        eqk.write('TEMPLATES/'+OUT,format='QUAKEML') 
        if WAV==1:
          for NET in INV:
              for STA in NET:
                  II=0;
                  NETNAME=NET.code 
                  STANAME=STA.code 
                  CHANAME=STA.channels[0].code 
                  print(NETNAME, STANAME, CHANAME)
                  II=0
                  try:
                      #print('FDSN')
                      WVF=client.get_waveforms(NETNAME,STANAME,'*',CHANAME,starttime=TOUT-5,endtime=TOUT+35) 
                      WVF.merge(fill_value='interpolate')
                      if DEC==1:
                        WVF.decimate(2)
                      WVF.write('TEMPLATES/' + str(TOUT)[: -4] + '.' + STANAME + '.ms',format='MSEED')
                      II=1
                  except:                    
                      print('NO FDSN DATA FOR STATION TRYING LOCAL: ' + NETNAME,STANAME,CHANAME)
                      II=0
                  if II==0:
                    try:
                        #print('SDS')
                        WVF=sdsclient.get_waveforms(NETNAME,STANAME,'*',CHANAME,starttime=TOUT-5,endtime=TOUT+35) 
                        WVF.merge(fill_value='interpolate')
                        if DEC==1:
                          WVF.decimate(2)
                        WVF.write('TEMPLATES/' + str(TOUT)[: -4] + '.' + STANAME + '.ms',format='MSEED')
                        #II=1
                    except:
                        print('NO LOCAL DATA FOR STATION: ' + NETNAME,STANAME,CHANAME)
                  else:
                    print('---')
                      
                

##### DOWNLOAD DAYLONG DATA AND CREATE DAY DIRECTORIES #####
import os
from obspy.clients.fdsn import Client 
from obspy import UTCDateTime

START="2022-05-26T00:00:00"
END  ="2022-06-10T00:00:00"
starttime=UTCDateTime(START)
endtime  =UTCDateTime(END)
import numpy as np

LATC=43.6867
LONC=11.2554 
RSTA=70/111 # km dal centro per le stazioni

DEC=0
####

#####
NDAYS=(endtime-starttime)/86400
TSTA=[]
TSTO=[]
for k in np.arange(0,NDAYS,1): 
  TSTA.append(starttime+((np.round(k)+1)*86400))

#####
client = Client("INGV")

INV=client.get_stations(starttime=starttime,endtime=endtime,channel='?HZ',network='*',latitude=LATC,longitude=LONC,maxradius=RSTA,level='response')
INV=INV.remove(network='IV',station='AQU')
INV=INV.remove(network='IV',station='RM??')
INV=INV.remove(network='*',station='*',channel='BHZ')
INV=INV.remove(network='*',station='*',channel='LHZ')
INV=INV.remove(network='*',station='*',channel='VHZ')


for NET in INV:              
  for STA in NET:                  
    II=0;
    NETNAME=NET.code
    STANAME=STA.code
    CHANAME=STA.channels[0].code
    print(NETNAME, STANAME, CHANAME)
    for D in TSTA:
      DAY=str(D-86400) 
      OUTDIR=DAY[0:4]+DAY[5:7]+DAY[8:10]
      print(OUTDIR)
      try:
        os.mkdir(OUTDIR)
      except:
        print('Dir esiste')
      try:     
        WVF=client.get_waveforms(NETNAME,STANAME,'*',CHANAME,starttime=D-86400,endtime=D)
        WVF.merge(fill_value='interpolate')
        print(WVF)
        print('writing ', NETNAME,STANAME,CHANAME)
        fil=WVF[0]
        fil.write(OUTDIR+'/' + DAY[0:-7] + fil.stats.station + '.' + fil.stats.network + '.' + fil.stats.channel,format='MSEED')                            
#        if DEC==1:
#          WVF.decimate(2)
        II=1
      except:
        print('NO FDSN')
        II=0
      if II==0:
        try:              
          print('SDS REQUEST: ', NETNAME,STANAME,CHANAME)            
          WVF=sdsclient.get_waveforms(NETNAME,STANAME,'*',CHANAME,starttime=D-86400,endtime=D)
          WVF.merge(fill_value='interpolate')
          if DEC==1:
            WVF.decimate(2)
          if len(WVF) > 0:
            II=1
          else:
            print('NO SDS DATA FOR STATION:',NETNAME,STANAME,CHANAME)
            II=0
        except:
          print('------') 
      if II==1:
        print('writing ', NETNAME,STANAME,CHANAME)
        fil=WVF[0]
        fil.write(OUTDIR+'/' + DAY[0:-7] + fil.stats.station + '.' + fil.stats.network + '.' + fil.stats.channel,format='MSEED')  
                          
                          