#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Thu Nov 12 15:25:09 2020

@author: piccinini
"""
import os
from obspy.clients.fdsn import Client
from obspy import UTCDateTime
client = Client("http://verna.pi.ingv.it:8080")

CSC_INI="2016-03-31T00:00:00"  # Data Inizio Cascinelle


latmin= 42.859
lonmin= 11.710
latmax= 42.891
lonmax= 11.756

# ESTRAGGO IL CATALOGO DA VERNA
client = Client('http://verna.pi.ingv.it:8080')
ST = UTCDateTime(CSC_INI)
ED = UTCDateTime(CSC_INI)+(86400)
try:
  cat = client.get_events(starttime=ST, endtime=ED,includearrivals=True)
  print(str(len(cat)) + ' eventi scaricati')
except:
  print('No events in catalog')


n=0
for ev in cat:
#    if ev.event_type=='earthquake':
#    print(ev)
    if ev.origins[0].evaluation_mode=='manual' :

        n=n+1
        TOUT=ev.origins[0].time
        OUT=str(TOUT)[0:-4] + '.xml'
        ev.write('TEMPLATES/'+OUT,format='QUAKEML')
        WVF=client.get_waveforms('TV','*','*','?HZ',starttime=TOUT-30,endtime=TOUT+30)
        WVF.merge(fill_value='interpolate')
        for fil in WVF:
          fil.write('TEMPLATES/' +str(TOUT)[0:-4]+ '.' +fil.stats.station+ '.ms',format='MSEED')
        WVF=client.get_waveforms('IV','ARCI','*','?HZ',starttime=TOUT-30,endtime=TOUT+30)
        WVF.merge(fill_value='interpolate')
        for fil in WVF:
          fil.write('TEMPLATES/' +str(TOUT)[0:-4]+ '.' +fil.stats.station+ '.ms',format='MSEED')
        WVF=client.get_waveforms('IV','MCIV','*','?HZ',starttime=TOUT-30,endtime=TOUT+30)
        WVF.merge(fill_value='interpolate')
        for fil in WVF:
          fil.write('TEMPLATES/' +str(TOUT)[0:-4]+ '.' +fil.stats.station+ '.ms',format='MSEED')

print (str(n) + ' eventi scritti')


#CSC_INI="2016-04-02T00:00:00"  # Data Inizio Cascinelle
#client = Client('http://verna.pi.ingv.it:8080')
#ST = UTCDateTime(CSC_INI)
#ED = UTCDateTime(CSC_INI)+(86400)


OUTDIR=CSC_INI[0:4]+CSC_INI[5:7]+CSC_INI[8:10]
try:
  os.mkdir(OUTDIR)
except:
  print('Dir esiste')
WVF=client.get_waveforms('TV','*','*','?HZ',starttime=ST,endtime=ED)
WVF.merge(fill_value='interpolate')
for sta in WVF:
    #sta.merge(fill_value='interpolate')
    sta.write(OUTDIR+'/'+CSC_INI+'.'+sta.stats.station+'.'+sta.stats.network+'.'+sta.stats.channel,format='MSEED')
WVF=client.get_waveforms('IV','ARCI','*','?HZ',starttime=ST,endtime=ED)
WVF.merge(fill_value='interpolate')
for sta in WVF:
    #sta.merge(fill_value='interpolate')
    sta.write(OUTDIR+'/'+CSC_INI+'.'+sta.stats.station+'.'+sta.stats.network+'.'+sta.stats.channel,format='MSEED')
WVF=client.get_waveforms('IV','MCIV','*','?HZ',starttime=ST,endtime=ED)
WVF.merge(fill_value='interpolate')
for sta in WVF:
    #sta.merge(fill_value='interpolate')
    sta.write(OUTDIR+'/'+CSC_INI+'.'+sta.stats.station+'.'+sta.stats.network+'.'+sta.stats.channel,format='MSEED')

print("DONE")
