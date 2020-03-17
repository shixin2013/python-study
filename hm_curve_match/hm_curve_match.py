#python分为可变结构与不可变结构，不可变结构基本等于复杂结构，包括list map等。复杂结构的赋值和传参都是传递的引用
# for it in list 中，it是只读的，修改它不会改变list值
import os
import datetime
import shutil
import configparser

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

import subprocess
import time  

import threading  
import zipfile




OUT_FILE_EXT = ['.runover','.out','.UDR','.nsp','.stp','.sstp','.DATA_TMP','.RSF','.wrs','.wrsh','.prp','.hex','.idx','.log','.ini','.DATA_TMP']
IN_FILE_EXT = ['.DATA']

def rateCurve(hist_dat, sim_dat, curve_nam ): #要求输入的曲线只包含报告步 RPTONLY, 否则结果不准确
    # score range is 0-100 ; each data point has a score, its arithmatic average is the final score
    # point score is calculated as (1-relErr)*100, minimum is 0
    # relErr is calculated as normal, other than the history data is below a min value. This min value is 1% of the overall average in the history data 
    # 大于1%的数据, 直接算相对误差, 设置最大值 100%; 小于1%的数据, 要么不考虑, 要么用绝对误差/1%值作为误差, 同样设置最大值100%
    assert len(hist_dat)==len(sim_dat)
    thread_hold_value = np.average(hist_dat)*0.05
    final_score = 0.0
    if thread_hold_value<1e-10:
        thread_hold_value = 1.0
    for i in range(len(hist_dat)):
        if hist_dat[i]<thread_hold_value:
            if sim_dat[i]<thread_hold_value:
                rel_err = 0
            else:
                rel_err = abs(hist_dat[i]-sim_dat[i])/thread_hold_value
        else:
            rel_err = abs(hist_dat[i]-sim_dat[i])/hist_dat[i]
        rel_err = min(rel_err, 1.0)
        final_score += 1-rel_err
    final_score /= len(sim_dat)
    final_score *= 100
    return final_score

def getBestCurve(hist_dat, sim_dat_list, curve_nam):
    curve_idx = 0
    best_rate = 0
    for idx, curve in enumerate(sim_dat_list):
        rate = rateCurve(hist_dat, curve, curve_nam) 
        if(rate>best_rate):
            best_rate = rate
            curve_idx = idx

    return curve_idx, best_rate

def highlightBestCurve(time_array, hist_dat, sim_dat_list, curve_nam):
    best_curve_idx, best_score = getBestCurve(hist_dat, sim_dat_list, curve_nam)
    plt.figure()
    #fig.subplot(1,1,1)
    plt.plot(time_array,hist_dat, color='blue', marker='o', markersize=2,linewidth=0)
    plt.xlabel('Time/days')
    plt.ylabel(curve_nam)
    
    for i in range(len(sim_dat_list)):
        if i==best_curve_idx:
            color_str = 'g'
        else:
            tmpd = 150/len(sim_dat_list)*i+50
            color_str = str(tmpd)
            plt.plot(time_array,sim_dat_list[i],color=color_str, label=str(i))
    plt.plot(time_array,sim_dat_list[best_curve_idx],color='g', label=str(best_curve_idx))
    
    plt.title(curve_nam + ' Best Score = '+str(int(best_score)))
    plt.legend()
    


def findAllDATA(rootDir): #找到DATA文件
    print('\n SCANNING DATA FILES...')
    alldatalist=[]
    rundatalist=[]
    for rt, dirs, files in os.walk(rootDir):  
        for file in files:
            fnam, fext = os.path.splitext(file)
            if  fext in IN_FILE_EXT:
                alldatalist.append(os.path.join(rt,file))
                if (fnam+'.runover') not in files:
                    rundatalist.append(os.path.join(rt,file))
    print('   -Files to be tested %d' %len(alldatalist))
    return alldatalist,rundatalist


def readUDR(UDR_FileName):
    fr = open(UDR_FileName,'r')
    fr.readline()
    fr.readline()
    tablehead1 = fr.readline().split() #line3
    fr.readline()
    tablehead2 = fr.readline().split() #line5
    fr.readline()
    fr.readline()
    if len(tablehead1)!=len(tablehead2):
        exit('met error when reading %s ' %UDR_FileName)

    lines = fr.readlines()
    tableArray = [[-1 for i in range(len(lines))] for n in range(len(tablehead1))]
    for i in range(len(tablehead2)):
        if tablehead2[i]=='-':
            tablehead2[i]=''
        else:
            tablehead2[i]=':'+tablehead2[i]

    tableHead = [tablehead1[i]+tablehead2[i] for i in range(len(tablehead1))]

    j=0
    for line in lines:
        i=0
        linesplit = line.split()
        for it in linesplit:
            tableArray[i][j]= float(it)
            i+=1
        j+=1

    fr.close()

    dic_head_array = {}
#    if len(dic_head_array)!=len(tableArray):
#        exit('met error when reading %s ' %UDR_FileName)
    for i in range(len(tableArray)):
        dic_head_array[tableHead[i]] = tableArray[i]

    return dic_head_array

def curveCompare(x,f_test,f_std):
    totalInteral = 0
    diffInteral = 0
    for i in range(len(x)-1):
        dx= x[i+1]-x[i]
        fave_std = (f_std(x[i])+f_std(x[i+1]))/2
        fave_test = (f_test(x[i])+f_test(x[i+1]))/2
        totalInteral+=abs(dx*fave_std)
        diffInteral+=abs(fave_test-fave_std)*dx
    if totalInteral<1e-5:
        relErr = 0
    else:
        relErr = diffInteral/totalInteral
    return relErr


def getUDR(spcArchDir): 
    UDRList=[]
    for rt, __, files in os.walk(spcArchDir):  
        for file in files:
            __, fext = os.path.splitext(file)
            if  fext == '.UDR':
                testUDR = os.path.join(rt,file)
                UDRList.append(testUDR)

    return  UDRList

 


#======================main========================
#-------run para-------
resultRoot = R'C:\Users\shixi\Desktop\debugCases\curveRating'    #path
tarType = 'WWPR'

curve_list =[]
time_dat =[]
hist_dat =[]
tarCurveType=[]
tarHistType=[]

UdrLists = getUDR(resultRoot)
anyDic = readUDR(UdrLists[0])
for key in anyDic:
    if key.split(':')[0] == tarType:
        tarCurveType.append(key)
    elif key.split(':')[0] == tarType+'H':
        tarHistType.append(key)
assert len(tarCurveType)==len(tarHistType)

for i in range(len(tarCurveType)):
    time_dat = []
    hist_dat = []
    curve_list = []
    for udr in UdrLists:
        testUdrDic = readUDR(udr)
        if len(time_dat)==0:
            time_dat.append( testUdrDic['TIME'])
            hist_dat.append( testUdrDic[tarHistType[i]])
        curve_list.append(testUdrDic[tarCurveType[i]])

    highlightBestCurve(time_dat[0],hist_dat[0],curve_list,tarCurveType[i])    
    
plt.show()

print('\n===================================')
print('--- Regression Over')
print('===================================\n')