# -*- coding: utf-8 -*-
"""
Created on Tue Oct  5 19:13:28 2021

@author: SmirkingJay
"""
import pandas as pd
import matplotlib.pyplot as plt
import datetime as dt
import pvlib
import numpy as np

to_skip=[*range(0,10)]
knmi=pd.read_csv('Eindhoven_KNMI.csv', skiprows=to_skip)
del knmi['# STN']
knmi.columns=['date','H','wind','temp','ghi']


knmi.date=knmi.date.astype(str)

knmi.H= knmi.H.apply(lambda x:str(x).zfill(2))
knmi.H=knmi.H.replace('24','00')

knmi['datetime']=knmi.date+knmi.H+'00'
knmi['datetime']=pd.to_datetime(knmi.datetime,format='%Y%m%d%H%M')

import datetime
knmi.datetime=knmi.datetime- datetime.timedelta(hours=0.5)
knmi.index=knmi.datetime

knmi=knmi[['wind','temp','ghi']]

knmi.ghi=knmi.ghi*2.7777
knmi.temp=knmi.temp*0.1
knmi.wind=knmi.wind*0.1

knmi.index=knmi.index.tz_localize('UTC')#.tz_convert('CET')

knmi=knmi.dropna()
KNMIData=knmi
### Calculations Solar_position with data from KNMI 
Solar_Position = pvlib.solarposition.ephemeris(time = KNMIData.index, latitude = 51.451, longitude = 5.377, pressure = 101325, temperature = KNMIData.temp)
Solar_Position=Solar_Position[Solar_Position.elevation>3.6]
Solar_Zenith=Solar_Position.zenith
KNMIData=KNMIData[KNMIData.index.isin(Solar_Position.index)]

### Calculating DNI with DIRINDEX
from pvlib import clearsky

apparent_zenith = Solar_Position.apparent_zenith
airmass = pvlib.atmosphere.get_relative_airmass(apparent_zenith)
altitude=22.6
pressure = pvlib.atmosphere.alt2pres(altitude)
airmass_absolute = pvlib.atmosphere.get_absolute_airmass(airmass, pressure)
linke_turbidity = pvlib.clearsky.lookup_linke_turbidity(KNMIData.index, 51.451,5.377)
cs = clearsky.ineichen(apparent_zenith, airmass_absolute, linke_turbidity, altitude,perez_enhancement=True)

DNI = pvlib.irradiance.dirindex(KNMIData.ghi,cs.ghi,cs.dni,Solar_Zenith,KNMIData.index,101325,True,None, 0.065,87)

### Calculating DHI
### GHI = DHI + DNI x cos(Z)
### DHI = GHI - DNI x cos(Z)
DNI = pd.DataFrame(pvlib.irradiance.dirindex(KNMIData.ghi,cs.ghi,cs.dni,Solar_Zenith,KNMIData.index,101325,True,None, 0.065,87))
Solar_Zenith_rad = np.radians(Solar_Position.zenith)
DHI = KNMIData.ghi - DNI.dni * np.cos(Solar_Zenith_rad)

orientation = {'A_SE':135,'A_SW':225,'B_E':90,'B_S':180,'B_W':270,'C_N':0,'C_S':180,'D_E':90,'D_W':270,}
def POA_all(a,b,c,d,e,f,g,h,i,o,**kwargs):
    x_1=kwargs.get('x_1')
    x_2=kwargs.get('x_2')
    if x_1 and x_2:
        o = {'A_SE':135,'A_SW':225,'B_E':90,'B_S':180,'B_W':270,'C_N':0,'C_S':180,'D_E':90,'D_W':270, 'roofA':135,'roofB':180}
        tilt_1 = {'A_SE':a, 'A_SW':b, 'B_E':c, 'B_S':d, 'B_W':e, 'C_N':f, 'C_S':g, 'D_E':h, 'D_W':i, 'roofA':x_1, 'roofB':x_2}
        poa_dict={} 
        poa_dict_total={}
        poa_direct={}
        poa_diffuse={}
        for surface in ['A_SE','A_SW', 'B_E', 'B_S', 'B_W', 'C_N', 'C_S', 'D_E','D_W','roofA','roofB']:
            poa_dict[surface] = pvlib.irradiance.get_total_irradiance(surface_tilt = tilt_1[surface], surface_azimuth = o[surface], solar_zenith = Solar_Zenith, solar_azimuth = Solar_Position.azimuth, dni = DNI.dni, ghi = KNMIData.ghi, dhi = DHI)
            poa_dict_total[surface] = sum(poa_dict[surface]['poa_global'].dropna())/1000
            poa_direct[surface]= poa_dict[surface]['poa_direct'].dropna()
            poa_diffuse[surface]=poa_dict[surface]['poa_diffuse'].dropna()
    else:
        tilt_2 = {'A_SE':a, 'A_SW':b, 'B_E':c, 'B_S':d, 'B_W':e, 'C_N':f, 'C_S':g, 'D_E':h, 'D_W':i}
        poa_dict={} 
        poa_dict_total={}
        poa_direct={}
        poa_diffuse={}
        for surface in ['A_SE','A_SW', 'B_E', 'B_S', 'B_W', 'C_N', 'C_S', 'D_E','D_W']:
            poa_dict[surface] = pvlib.irradiance.get_total_irradiance(surface_tilt = tilt_2[surface], surface_azimuth = o[surface], solar_zenith = Solar_Zenith, solar_azimuth = Solar_Position.azimuth, dni = DNI.dni, ghi = KNMIData.ghi, dhi = DHI)
            poa_dict_total[surface] = sum(poa_dict[surface]['poa_global'].dropna())/1000
            poa_direct[surface]= poa_dict[surface]['poa_direct'].dropna()
            poa_diffuse[surface]=poa_dict[surface]['poa_diffuse'].dropna()
    return poa_dict, poa_dict_total, poa_direct, poa_diffuse

POA_All, POA_Total_ALL, POA_DIRECT_ALL, POA_DIFFUSE_ALL=POA_all(90,90,90,90,90,40,40,40,40,orientation)


def POA(a,b,c,d,e,f,g,h,o):
    tilt = {'10':a, '15':b, '20':c, '25':d, '30':e, '35':f, '40':g, '45':h}
    poa_dict={} 
    poa_dict_total={}
    poa_direct={}
    poa_diffuse={}
    for surb in ['10','15', '20', '25', '30', '35', '40', '45']:
        poa_dict[surb] = pvlib.irradiance.get_total_irradiance(surface_tilt = tilt[surb], surface_azimuth = o, solar_zenith = Solar_Zenith, solar_azimuth = Solar_Position.azimuth, dni = DNI.dni, ghi = KNMIData.ghi, dhi = DHI)
        poa_dict_total[surb] = sum(poa_dict[surb]['poa_global'].dropna())/1000
    return poa_dict, poa_dict_total

POA_roofB,  POA_Total_roofB=   POA(10,15,20,25,30,35,40,45,180)
POA_roofA_SE,   POA_Total_roofASE=POA(10,15,20,25,30,35,40,45,135)
POA_roofA_SW,   POA_Total_roofASW=POA(10,15,20,25,30,35,40,45,225)

POA_All, POA_Total_ALL, POA_DIRECT_ALL, POA_DIFFUSE_ALL=POA_all(90,90,90,90,90,40,40,40,40,orientation,x_1=30,x_2=30)


def PLOT(x,height,title,xl,yl,**kwargs):
    fig, ax=plt.subplots()
    Bot=kwargs.get('Bot',None)
    Top=kwargs.get('Top',None)
    x_1=kwargs.get('x_1')
    y_1=kwargs.get('y_1')
    x_2=kwargs.get('x_2')
    y_2=kwargs.get('y_2')
    leg_1=kwargs.get('leg_1')
    leg_2=kwargs.get('leg_2')
    leg_3=kwargs.get('leg_3')
    r1=ax.bar(x, height, width=0.8,align='center')
    if x_1 and y_1:
        r2=ax.bar(x_1, y_1, width=0.4,align='edge')
    if x_2 and y_2:
        r3=ax.bar(x_2, y_2, width=0.2,align='edge')
    plt.ylim(bottom=Bot,top=Top)
    plt.title(title,fontsize=13)
    plt.xlabel(xl, fontsize=10)
    plt.ylabel(yl, fontsize=10)
    if leg_1 and leg_2:
        ax.legend((leg_1,leg_2,leg_3),loc='upper left')
    fig.tight_layout()
    return plt.show()

plot_roofB=PLOT(POA_Total_roofB.keys(),POA_Total_roofB.values(),'Comparison of Tilt Angles for Building B', 'Tilt angles [Degrees]','Total Annual POA [KWh/m2]',Bot=1000,Top=1250)
plot_roofA=PLOT(POA_Total_roofASE.keys(),POA_Total_roofASE.values(),'Comparison of Tilt Angles and Orientations for Building A','Tilt angles [Degrees]','Total Annual POA [KWh/m2]',Bot=900,Top=1300,x_1=POA_Total_roofASW.keys(),y_1=POA_Total_roofASW.values(),leg_1='Roof A-SE',leg_2='Roof A-SW')
plot_all_surfaces=PLOT(POA_Total_ALL.keys(),POA_Total_ALL.values(),'Comparison of POA for all surfaces','Surfaces','Total POA [KWh/m2]',Bot=400,Top=1300)

#Question 3:

module_parameters = pd.read_excel('Module parameters.xlsx', index_col= 'Parameters')

def Temp_Cell(x,a,b,delta_T):
    temp_cell= {}
    for cell in x:
        temp_cell[cell] = pvlib.temperature.sapm_cell(POA_Total_ALL[cell],KNMIData.temp, KNMIData.wind, a = a, b = b, deltaT = delta_T)
    return temp_cell
Tempcell_surfaces=Temp_Cell(['roofA','roofB'],a=-3.5,b=-0.0672,delta_T=3)
Tempcell_surfaces.update(Temp_Cell(['A_SE','A_SW','B_E','B_S','B_W'],a=-2.81,b=-0.0455,delta_T=0))
Tempcell_surfaces.update(Temp_Cell(['C_N','C_S','D_E',"D_W"],a=-2.98,b=-0.0471,delta_T=1))

orientation = {'A_SE':135,'A_SW':225,'B_E':90,'B_S':180,'B_W':270,'C_N':0,'C_S':180,'D_E':90,'D_W':270, 'roofA':135,'roofB':180}
def AOI(a,b,c,d,e,f,g,h,i,j,k,):
    aoi = {}
    tilt = {'A_SE':a, 'A_SW':b, 'B_E':c, 'B_S':d, 'B_W':e, 'C_N':f, 'C_S':g, 'D_E':h, 'D_W':i, 'roofA':j, 'roofB':k}
    for a in ['A_SE','A_SW', 'B_E', 'B_S', 'B_W', 'C_N', 'C_S', 'D_E', 'D_W', 'roofA', 'roofB']:
        aoi[a] = pvlib.irradiance.aoi(surface_tilt = tilt[a], surface_azimuth = orientation[a], solar_zenith = Solar_Position.zenith, solar_azimuth = Solar_Position.azimuth)
    return aoi
AOI_All=AOI(90,90,90,90,90,40,40,40,40,30,30)

def EFF_IRR(x):    
    effect_irra= {}
    for irra in ['A_SE','A_SW', 'B_E', 'B_S', 'B_W', 'C_N', 'C_S', 'D_E', 'D_W', 'roofA', 'roofB']: 
        effect_irra[irra] = pvlib.pvsystem.sapm_effective_irradiance(poa_direct = POA_DIRECT_ALL[irra], poa_diffuse = POA_DIFFUSE_ALL[irra], airmass_absolute = airmass, aoi = AOI_All[irra], module = module_parameters.loc['Cells_in_Series':,x])    
    return effect_irra
Effect_Irr_HIT=EFF_IRR('HIT')
Effect_Irr_CdTe=EFF_IRR('CdTe')  
Effect_Irr_monoSi=EFF_IRR('mono-Si')   

print('eff irr ase=',sum(Effect_Irr_HIT['A_SE']))
print('eff irr asw=',sum(Effect_Irr_HIT['A_SW']))
print('direct irr ase=',sum(POA_DIRECT_ALL['A_SE'])) 
print('direct irr asw=',sum(POA_DIRECT_ALL['A_SW'])) 
print('diffuse irr asw=',sum(POA_DIFFUSE_ALL['A_SW'])) 
print('diffuse irr asw=',sum(POA_DIFFUSE_ALL['A_SW'])) 

HIT = {'C_N':95, 'C_S':95, 'D_E':95, 'D_W':95, 'A_SW':1190, 'A_SE':1428, 'B_W':214, 'B_S':214, 'B_E':357, 'roofA':1190, 'roofB':595}
CdTe = {'C_N':167, 'C_S':167, 'D_E':167, 'D_W':167, 'A_SW':2083, 'A_SE':2500, 'B_W':375, 'B_S':375, 'B_E':625, 'roofA':2083, 'roofB':1041}
mono = {'C_N':62, 'C_S':62, 'D_E':62, 'D_W':62, 'A_SW':777, 'A_SE':932, 'B_W':139, 'B_S':139, 'B_E':233, 'roofA':777, 'roofB':388}
area_surfaces = {'C_N':120.79, 'C_S':120.79, 'D_E':120.79, 'D_W':120.79, 'A_SW':1500, 'A_SE':1800, 'B_W':270, 'B_S':270, 'B_E':450, 'roofA':1500, 'roofB':750}
def Get_DC(x,y):
    dc = {}
    dc_total={}
    dc_total_yield={}
    for DC in ['A_SE','A_SW', 'B_E', 'B_S', 'B_W', 'C_N', 'C_S', 'D_E', 'D_W', 'roofA', 'roofB']:
        dc[DC] = pvlib.pvsystem.sapm(Effect_Irr_HIT[DC], Tempcell_surfaces[DC], module_parameters[x])
        dc[DC]=dc[DC]['p_mp'].dropna()
        #dc_total[DC]=sum(dc[DC]*y[DC])/(1000000)
        dc_total[DC]=sum(dc[DC])/(1000)
        dc_total_yield[DC]=sum(dc[DC]*y[DC])/(1000*area_surfaces[DC])
    return dc,dc_total,dc_total_yield
DC_HIT,  DC_Total_HIT, DC_Total_HIT_yield=Get_DC('HIT',HIT)  
DC_CdTe,    DC_Total_CdTe,DC_Total_CdTe_yield=Get_DC('CdTe',CdTe)      
DC_monoSi,  DC_Total_monoSi,DC_Total_monoSi_yield=Get_DC('mono-Si',mono)     

#plot_DC_HIT=PLOT(DC_Total_HIT.keys(),DC_Total_HIT.values(),'Comparison of all surfaces for HIT', 'Surfaces','Total DC Power [KWh/panel]',Bot=None,Top=330000)   
#plot_DC_CdTe=PLOT(DC_Total_CdTe.keys(),DC_Total_CdTe.values(),'Comparison of all surfaces for CdTe', 'Surfaces','Total DC Power [KWh/m2]',Bot=None,Top=250)    
#plot_DC_monoSi=PLOT(DC_Total_monoSi.keys(),DC_Total_monoSi.values(),'Comparison of all surfaces for monoSi', 'Surfaces','Total DC Power [KWh]',Bot=None,Top=330000)    
#plot_roofA=PLOT(DC_Total_HIT.keys(),DC_Total_HIT.values(),'Comparison of Tilt Angles and Orientations for Building A','Tilt angles [Degrees]','Total POA [KWh/m2]',x_1=DC_Total_monoSi.keys(),y_1=DC_Total_monoSi.values(),leg_1='HIT',leg_2='MonoSi', Top=330000)
plot_roofA_1=PLOT(DC_Total_HIT.keys(),DC_Total_HIT.values(),'Comparison of annual DC output for all surfaces','Surfaces','Total DC Power [MWh]',x_1=DC_Total_monoSi.keys(),y_1=DC_Total_monoSi.values(),x_2=DC_Total_CdTe.keys(),y_2=DC_Total_CdTe.values(),leg_1='HIT',leg_2='MonoSi', leg_3='CdTe', Top=330)
plot_roofA_1=PLOT(DC_Total_HIT_yield.keys(),DC_Total_HIT_yield.values(),'Comparison of annual DC yield for all surfaces','Surfaces','Total DC Power [KWh/m2]',x_1=DC_Total_monoSi_yield.keys(),y_1=DC_Total_monoSi_yield.values(),x_2=DC_Total_CdTe_yield.keys(),y_2=DC_Total_CdTe_yield.values(),leg_1='HIT',leg_2='MonoSi', leg_3='CdTe', Top=330)

dc_1={}
dc_1=dict(DC_HIT)

def AC_Eff(pdc):
    nom_eff=0.96    
    pdc0=module_parameters.loc['Wp','HIT']/nom_eff
    if pdc==0:
        pac=0
    elif pdc0>pdc>0:
        tau=pdc/pdc0
        eff=(-0.0162*tau)-(0.0059/tau)+0.9858
        pac=pdc*eff
    elif pdc>pdc0:
        pac=module_parameters.loc['Wp','HIT']
    return pac

Power_AC={}
Power_AC_Total={}
for x in ['A_SE','A_SW', 'B_E', 'B_S', 'B_W', 'C_N', 'C_S', 'D_E', 'D_W', 'roofA', 'roofB']: 
    v=pd.DataFrame(dc_1[x])
    y=[]
    for z in v.p_mp:
        p=AC_Eff(z)
        y.append(p)
    n=pd.DataFrame(y, index=Solar_Position.index)
    n.columns=['p_mp']
    Power_AC[x]=n
  

for durface in ['A_SE','A_SW', 'B_E', 'B_S', 'B_W', 'C_N', 'C_S', 'D_E', 'D_W', 'roofA', 'roofB']:
    Power_AC_Total[durface] = sum(Power_AC[durface]['p_mp']*HIT[durface])/(1000000)
    
plot_AC_total=PLOT(Power_AC_Total.keys(),Power_AC_Total.values(),'Comparison of annual AC power for all surfaces (HIT)', 'Surfaces','Total AC Power [MWh]',Bot=None,Top=330) 


Power_AC_A=Power_AC_Total['A_SE']+Power_AC_Total['A_SW']+Power_AC_Total['roofA']
Power_AC_B=Power_AC_Total['B_E']+Power_AC_Total['B_S']+Power_AC_Total['B_W']+Power_AC_Total['roofB']
Power_AC_C=Power_AC_Total['C_N']+Power_AC_Total['C_S']
Power_AC_D=Power_AC_Total['D_E']+Power_AC_Total['D_W']


AC_Buildings={}
AC_Buildings['A']=Power_AC_A
AC_Buildings['B']=Power_AC_B
AC_Buildings['C']=Power_AC_C
AC_Buildings['D']=Power_AC_D

plot_AC_Buildings=PLOT(AC_Buildings.keys(),AC_Buildings.values(),'Comparison of annual AC power for all Buildings (HIT)', 'Buildings','Total AC Power [MWh]',Bot=0,Top=800)


'''y=[]
r=pd.DataFrame(index=Solar_Position.index)
for x in ['A_SE','A_SW', 'B_E', 'B_S', 'B_W', 'C_N', 'C_S', 'D_E', 'D_W', 'roofA', 'roofB']: 
    i=pd.DataFrame(Power_AC[x])
    t=pd.DataFrame(DC_HIT[x])
    r['1']=i['p_mp']
    r['2']=t['p_mp']
    f=r['2']*0.96
    g=-r['1']+f
    for z in f:
        if z<0:
            y.append(z)
        else:
            continue'''

Irr_HIT_daily={}
for vurface in Effect_Irr_HIT:      
        Irr_HIT_daily[vurface]=Effect_Irr_HIT[vurface].resample('d').sum()   
    
#index=Irr_HIT_daily['A_SE'].index
sa_spring={}
sa_summer={}
sa_autumn={}
for turface in ['A_SE','A_SW', 'B_E', 'B_S', 'B_W', 'C_N', 'C_S', 'D_E', 'D_W', 'roofA', 'roofB']:
    sa_spring[turface]=(Power_AC[turface].loc['2019-05-15 00:00:00':'2019-05-15 23:00:00'])
    sa_summer[turface]=(Power_AC[turface].loc['2019-06-27 00:00:00':'2019-06-27 23:00:00'])
    sa_autumn[turface]=(Power_AC[turface].loc['2019-09-14 00:00:00':'2019-09-14 23:00:00'])


for turface in ['A_SE','A_SW', 'B_E', 'B_S', 'B_W', 'C_N', 'C_S', 'D_E', 'D_W', 'roofA', 'roofB']:
    sa_spring[turface].columns=[turface]
    sa_summer[turface].columns=[turface]
    sa_autumn[turface].columns=[turface]

def PLOT_ac(**kwargs):
    a=kwargs.get('a')
    b=kwargs.get('b')
    c=kwargs.get('c')
    d=kwargs.get('d')
    if a:
        ax=sa_spring[a].plot(ylabel='AC output (W)', title='Hourly Variation of AC output for a day in Spring',)
    if b:
        sa_spring[b].plot(ax=ax)
    if c:
        sa_spring[c].plot(ax=ax)
    if d:
       sa_spring[d].plot(ax=ax)   
    
    if a:
        ax=sa_summer[a].plot(ylabel='AC output (W)', title='Hourly Variation of AC output for a day in Summer')
        #plt.xlim(left='05:00',right='19:00')
    if b:    
        sa_summer[b].plot(ax=ax)
    if c:
        sa_summer[c].plot(ax=ax)
    if d:
        sa_summer[d].plot(ax=ax)
            
    if a:
        ax=sa_autumn[a].plot(ylabel='AC output (W)', title='Hourly Variation of AC output for a day in Autumn')
        #plt.xlim(left='05:00',right='19:00')
    if b:
        sa_autumn[b].plot(ax=ax)
    if c:
        sa_autumn[c].plot(ax=ax)
    if d:
        sa_autumn[d].plot(ax=ax)
    
        

A_seasons_plot=PLOT_ac(a='A_SE',b='A_SW',c='roofA')
B_seasons_plot=PLOT_ac(a='B_E',b='B_S',c='B_W',d='roofB')
C_seasons_plot=PLOT_ac(a='C_N',b='C_S')
D_seasons_plot=PLOT_ac(a='D_E',b='D_W')



    