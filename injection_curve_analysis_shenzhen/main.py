#python分为可变结构与不可变结构，不可变结构基本等于复杂结构，包括list map等。复杂结构的赋值和传参都是传递的引用
# for it in list 中，it是只读的，修改它不会改变list值
# ctrl+[,] 代码缩进
import os
import datetime
import shutil
import configparser
import sys

import numpy as np
from scipy.optimize import leastsq
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import pandas as pd

import subprocess
import time  

import threading  
import math
import xlrd

def pres_loss_friction(L,D,rou,Q,visc,e):
    if Q<1e-3 :
        return 0
    Cf = 4.343e-15 #metric
    #get reynolds number
    Cr = 0.01474 #metric
    Re = Cr*rou*Q/(D*visc)
    #get fanning friction factor
    if Re<4000:
        f = 16/Re
    else:
        t=pow(e/(3.7*D),10/9)
        t=t+6.9/Re
        t = math.log10(t)*(-3.6)
        f = 1/t/t
    return Cf*f*L*rou*Q*Q/pow(D,5)*0.1 #bar to MPa

def pres_loss_friction_array_ver(L,D,rou,Q,visc,e):
    Cf = 4.343e-15 #metric
    Cr = 0.01474 #metric
    
    rev = Q.copy()
    for i in range(len(Q)):
        if Q[i]<1e-3 :
            rev[i]=0
        else:
            #get reynolds number
            Re = Cr*rou*Q[i]/(D*visc)
            #get fanning friction factor
            if Re<4000:
                f = 16/Re
            else:
                t=pow(e/(3.7*D),10/9)
                t=t+6.9/Re
                t = math.log10(t)*(-3.6)
                f = 1/t/t
            rev[i] = Cf*f*L*rou*Q[i]*Q[i]/pow(D,5)*0.1 #bar to MPa
    
    return rev

def pres_loss_head(H,rou):
    return rou*9.8*H*1e-6


def pres_loss_singleICD(fld,q,strength):
    global fluid_wat
    dens_cal = fluid_wat.dens
    visc_cal = fluid_wat.visc
    K = strength # / dens_cal =1.0 #str:bars/(rm3/day)2
    dp = pow((dens_cal/fld.dens)*(fld.visc/visc_cal),0.25)*fld.dens/dens_cal*K*q*q
    return dp*0.1 #bar 2 MPa

class Fluid:
    def __init__(self,rou,visc):
        self.dens = rou
        self.visc = visc

class Well:
    def __init__(self,l,l2,h,d,e,wn,fluid,coff,strength=0.00021,PI=1.0): #ecl typical value
        self.length_v = l
        self.length_h = l2-l
        self.tvd = h
        self.diam = d*0.0254
        self.rough = e*0.001
        self.name = wn
        self.fluid = fluid
        self.coff = coff
        self.strength = strength
        self.PI = PI
    
    def get_pres_loss_icd(self,Q,f):
        global icd_interval
        nIcds = (int)(self.length_h/icd_interval)
        #print(Q)
        dp = pres_loss_singleICD(f,Q/nIcds,self.strength)
        return dp

    def get_bhp(self,Q,thp):
        bhp = thp+pres_loss_head( self.tvd, self.fluid.dens)- \
            self.coff*pres_loss_friction(self.length_v,self.diam,self.fluid.dens,Q,self.fluid.visc,self.rough)- \
                self.get_pres_loss_icd(Q,fluid)
                #pres_loss_local(self.fluid.dens,self.coff,Q,self.diam) #3.7
        # if bhp<resP:
        #     bhp = resP
        return bhp
    def get_effective_pres(self,Q,thp,resP):
        tmp = self.get_bhp(Q,thp)-resP
        if tmp<0:
            tmp=0
        return tmp
    
    def read_thp_curve(self,fname):
        wb = xlrd.open_workbook(fname)
        st = wb.sheet_by_name(self.name)
        self.timeT = st.col_values(0)[1:]
        self.presT = st.col_values(1)[1:]
        self.injT = st.col_values(2)[1:]  #list 不同于array，不能直接乘。乘代表重复列表多次组成新列表
        self.injT = [i*1e-3*60*24 for i in self.injT] 
        self.retT = st.col_values(3)[1:]
        self.retT = [abs(i)*1e-3*60*24 for i in self.retT] 
        #对流量做平均，避免剧烈波动。向前5步平均
        self.ave_rate(self.injT,200)
        
        self.effec_pres = [i for i in range(len(self.timeT))]

        global resp, plot_icd_pl,plot_pres_loss
        self.effec_PI = self.effec_pres.copy()
        for i in range(len(self.timeT)):
            if not plot_icd_pl and not plot_pres_loss:
                self.effec_pres[i] = self.get_effective_pres(abs(self.injT[i]),self.presT[i],resp)
                if self.effec_pres[i]>0.1:
                    self.effec_PI[i]= (self.injT[i]-self.retT[i])/self.effec_pres[i]/10
                    if self.effec_PI[i]>self.PI:
                        self.effec_PI[i]=self.PI
                    elif self.effec_PI[i]<0:
                        self.effec_PI[i]=0
                else:
                    self.effec_PI[i]= 0
                    
            elif plot_pres_loss:
                self.effec_pres[i] = self.coff*pres_loss_friction(self.length_v,self.diam,self.fluid.dens,abs(self.injT[i]),self.fluid.visc,self.rough)+ \
                    self.get_pres_loss_icd(abs(self.injT[i]),fluid)
            else:
                self.effec_pres[i] = self.get_pres_loss_icd(abs(self.injT[i]),self.fluid)
        #处理时间，标准化
        for i in range(len(self.timeT)):
            if i==0 :
                continue
            if self.timeT[i]<self.timeT[i-1]:
                self.timeT[i:] = [it+1.0 for it in self.timeT[i:]]
        self.timeT = [it-self.timeT[0] for it in self.timeT]

    def read_wat_test_dat(self,fname,ithst=0):
        wb = xlrd.open_workbook(fname)
        st = wb.sheet_by_index(ithst)
        start_line = 1
        self.wattest_tT =   st.col_values(1)[start_line:]
        self.wattest_pT =   st.col_values(2)[start_line:] #MPa
        self.wattest_inT =  st.col_values(3)[start_line:] #L/min 
        self.wattest_outT = st.col_values(4)[start_line:] #L/min 
        
        self.wattest_pT = np.array( self.wattest_pT )
        self.wattest_inT  = np.array( [i*1e-3*60*24 for i in self.wattest_inT] )#L/min to m3/day
        self.wattest_outT = np.array( [abs(i)*1e-3*60*24 for i in self.wattest_outT] )#L/min to m3/day
    
    def ave_rate(self,tarT,istep):
        #average rate value to avoid huge oscillation. for water test data only for now
        tmpl = tarT.copy()
        for i in range(len(tmpl)):
            j=i
            n=0
            total = 0.0
            while j>=0 and n<istep:
                total+=tmpl[j]
                j-=1
                n+=1
            tarT[i] = total / n
    
    def wat_test_calc_p(self):
        global resp,fluid_wat
        intoResT = self.wattest_inT - self.wattest_outT
        #intoResT[intoResT<0] = 0
        #q = pi(bhp-resp)
        bhp = intoResT/self.PI*0.1 + resp
        #thp = bhp-head+dp
        dp_icd = self.get_pres_loss_icd(self.wattest_inT,fluid_wat)
        dp_fric = self.coff*pres_loss_friction_array_ver(self.length_v,self.diam,fluid_wat.dens,self.wattest_inT,fluid_wat.visc,self.rough)
        dp_head = pres_loss_head( self.tvd, fluid_wat.dens)
        self.wattest_calcpT = bhp-dp_head+dp_icd+dp_fric
        self.wattest_DP = dp_icd+dp_fric

    def set_PI(self,pi):
        self.PI = pi

    def set_K(self,k):
        self.strength = k
    
    def set_coff(self,c):
        self.coff = c

def pres_loss(coff,L,D,rou,Q,visc,e):
    return pres_loss_friction(L,D,rou,Q,visc,e)*coff


def wat_test_thp_calc(p,well):
    #well.set_coff(p[0])
    #well.set_K(p)
    well.set_PI(p[1])
    well.wat_test_calc_p()
    return well.wattest_calcpT

def error_thp_calc(p,well,y):
    return wat_test_thp_calc(p,well)-y


#脚本使用流程：1) 使用斜井段长度，拟合A4H5实测摩阻 5.4 
#            2）使用得到的摩阻系数计算斜井段压力损失  3）认为流量沿ICD均匀分布，压差也相同（相当于定了ICD参数），算出第一个ICD压降即可
#               总注入/ICD数 即为流入ICD的流速，按这个流速算压降。
# WELL   V_MD     H_MD
# A4H5   3500     4300
# C3H4   4000     5068
# A3H3   1650     2500
# A5H1   1483     2205

#region match
# v=[8,7,6,5,4,2] #in bpm, 1 bpm = 1800 m3/day
# bpm_2_m3_day = 228.9417
# p2=[7,5.4,4.5,3.1,2,0.65] #in MPa
# v = [it*bpm_2_m3_day for it in v]

# L = 3500 #m
# D = 4.7*0.0254 #inch 2 m
# rou = 1000 # kg/m3
# visc = 0.7 #cp
# e = 0.3*0.001 #mm 2 m

# p_test = [pres_loss(5.4,L,D,rou,it,visc,e) for it in v]

# plt.plot(v,p2,'o',label='实测阻力')
# plt.plot(v,p_test,'r',label='计算阻力')
# plt.legend(loc='upper left')
# #plt.set_xlabel(xlabel='流速(m3/day)',ylabel='压力损失(MPa)')
# plt.xlabel('流速(m3/day)')
# plt.ylabel('压力损失(MPa)')
# plt.show()

# # # y = np.array(p2)
# # # x = np.array(v)
# # # pset = [L,D,rou,visc,e]
# # # para =leastsq(error, 3.7, args=(x,y)) # 进行拟合

# # # p_match = [pres_loss(para,L,D,rou,it,visc,e) for it in v]
# # # plt.plot(v,p_match,'g',label='拟合阻力')

# sys.exit()
#endregion

dat_file_nam = R"C:\Users\shixi\Desktop\workContracts\ShenZhenLh11\WorkSpace\inje_dat_for_reading.xlsx"
resp = 12.4 #Mpa
dens = 1050
visc = 2.0

fluid = Fluid(dens,visc)

fluid_wat = Fluid(1000,0.7)

#a4h5 len tvd diam_in rough_mm 
# wA4H5 = Well(3450,1240,4.7,0.3,'A4H5',fluid)
# wA4H5.read_thp_curve(dat_file_nam)

# wA5H1 = Well(3450,1240,4.7,0.3,'A5H1',fluid)
# wA5H1.read_thp_curve(dat_file_nam)

# wA4H5 = Well(3450,1240,4.7,0.3,'A4H5',fluid)
# wA4H5.read_thp_curve(dat_file_nam)

# wA4H5 = Well(3450,1240,4.7,0.3,'A4H5',fluid)
# wA4H5.read_thp_curve(dat_file_nam)

# #plot
# plt.plot(wA4H5.timeT,wA4H5.effec_pres,label='有效压力')
# plt.plot(wA4H5.timeT,wA4H5.presT,label='泵注压力')
# plt.legend(loc='upper right') 
# plt.show()

icd_interval=10 #m ?需确定 K需确定（开放注入？）
plot_icd_pl = False
plot_wat_test_match = False
plot_pres_loss = True

wlist = []
# 摩阻系数清水拟合值：5.4
#                    vlen hlen  tvd diam_in rough_mm,              coff, K  PI(40-100),  
wlist.append(Well(3500,4300, 1240,4.7,     0.3,    'A4H5',fluid, 2.37,0.021,60))
wlist.append(Well(4000,5068, 1246,4.7,     0.3,    'C3H4',fluid,   2.37, 0.021,100 )) #最大5mpa摩阻，考虑到a4h5他们估算是4mpa，c3h4较长，基本合理
wlist.append(Well(1650,2500, 1214,4.7,     0.3,    'A3H3',fluid,    2.37, 0.021,94))#1.7
wlist.append(Well(1483,2205, 1240,4.7,     0.3,    'A5H1',fluid,     2.37, 0.021,55 ))#3.4 太大了，有点扯
# wlist.append(Well(2200,3000, 1240,4.7,     0.3,    'test',fluid,     2.37, 0.021,55 ))#3.4 太大了，有点扯
# qqq = 1200*1e-3*60*24
# dp1 = -pres_loss_head( wlist[0].tvd, wlist[0].fluid.dens)
# dp2 = wlist[0].coff*pres_loss_friction(wlist[0].length_v,wlist[0].diam,wlist[0].fluid.dens,qqq,wlist[0].fluid.visc,wlist[0].rough)
# dp3 = wlist[0].get_pres_loss_icd(qqq,fluid)
# dp = dp1+dp2+dp3+resp
                
# print(dp1,'',dp2,'',dp3,'',dp)
# exit(0)
nfig = 4
fig, axs = plt.subplots(nfig)
plt.figure()
i=0

if plot_wat_test_match:
    root = R"C:\Users\shixi\Desktop\workContracts\深圳礁灰岩\资料\4连续封隔体控水\清水实验资料\\"
    root = os.path.dirname(root)
    ave_step = 1
    for w in wlist:
        fnam = os.path.join(root, w.name+'井清水测试.xlsx')
        for i in range(1): #根据是否注释掉，切换一个井多个测试数据和一个井一个测试数据的情况
            #w.read_wat_test_dat(fnam)
            #print('reading sheet 0 of c3h4/a5h1')
            w.read_wat_test_dat(fnam,2)
            print('reading sheet 1-6 of a3h3')
            #axs[i].plot(w.wattest_tT,w.wattest_inT,'go',label=w.name+'ori in')
            #axs[i].plot(w.wattest_tT,w.wattest_outT,label=w.name+'ori out')
            w.ave_rate(w.wattest_inT,ave_step)
            w.ave_rate(w.wattest_outT,ave_step)
            #axs[i].plot(w.wattest_tT,w.wattest_inT,'r',label=w.name+'ave in')
            #axs[i].plot(w.wattest_tT,w.wattest_outT,label=w.name+'ave out')
            #axs[i].legend(loc='upper left')

            axs[i].set(ylabel='MPa')
            axs[i].plot(w.wattest_tT,w.wattest_pT,label=w.name+'记录井口压力')
            w.wat_test_calc_p()
            axs[i].plot(w.wattest_tT,w.wattest_calcpT,label=w.name+'计算井口压力')
            axs[i].plot(w.wattest_tT,w.wattest_DP,label=w.name+'压力损失')
            axs[i].grid()
            #match
            # p0 = [2.37,100] # 拟合的初始参数设置
            # para =leastsq(error_thp_calc, p0, args=(w,w.wattest_pT)) # 进行拟合
            # y_fitted = wat_test_thp_calc(para[0],w) # 画出拟合后的曲线
            # axs[i].plot(w.wattest_tT,y_fitted,label=w.name+'拟合c,pi:'+str(para[0][0])+' '+str(para[0][1]))

            # w.set_K(0.015) #w.set_PI(100)
            # w.wat_test_calc_p()
            # axs[i].plot(w.wattest_tT,w.wattest_calcpT,label=w.name+'计算井口压力'+str(w.strength))

            # w.set_K(0.01) #w.set_PI(100)
            # w.wat_test_calc_p()
            # axs[i].plot(w.wattest_tT,w.wattest_calcpT,label=w.name+'计算井口压力'+str(w.strength))

            # ax2 = axs[i].twinx()
            # ax2.plot(w.wattest_tT,w.wattest_inT,'g',label=w.name+'入口流量')
            # ax2.plot(w.wattest_tT,w.wattest_inT-w.wattest_outT,'b',label=w.name+'净注入量')
            # ax2.set(ylabel='m3/day')

            axs[i].legend(loc='upper left')

            #dataframe = pd.DataFrame({'Time':w.wattest_tT,'Caculate Pres(MPa)':w.wattest_calcpT,'Record Pres(MPa)':w.wattest_pT})
            #dataframe.to_csv(os.path.join(root, w.name+'_'+str(i)+'th_wat_test.csv'),index=False,sep=',')

            i+=1
else:
    root = os.path.dirname(dat_file_nam)
    for w in wlist:
        print(w.name+'begin')
        w.read_thp_curve(dat_file_nam)
        axs[i].plot(w.timeT,w.effec_pres,'r',label=w.name+'有效压力')
        axs[i].plot(w.timeT,w.presT,'b',label=w.name+'泵注压力')
        ax2 = axs[i].twinx()
        ax2.plot(w.timeT,w.effec_PI,'g',label=w.name+'有效PI')
        axs[i].legend(loc='upper right')
        axs[i].set(ylabel='MPa')
        ax2.set(ylabel='PI')
        ax2.legend(loc='upper right')
        plt.plot(w.timeT,w.effec_pres,label=w.name+'有效压力')
        #plt.plot(w.timeT,w.effec_PI,label=w.name+'有效PI')
        plt.legend(loc='upper right')
        
        #plt.set(ylabel='MPa')
        #output csv
        #dataframe = pd.DataFrame({'Time':w.timeT,'Effective_pressure':w.effec_pres})
        #dataframe.to_csv(os.path.join(root, w.name+'_effcPres.csv'),index=False,sep=',')
        i+=1

plt.show()
