#对比ecl prt中trans结果与cloud trans结果,NNC要实用 'ALLNNC'
##参数

PRTPATH = r"C:\Users\shixi\Documents\newLandProj\zaoyuan\SimCases\C_project\C\C.PRT"
CLOUDTRANSPATH=r"C:\Users\shixi\Desktop\trans.test"
NX=71   
NY=140 
NZ=42

##PRTPATH = r"C:\Users\shixi\Desktop\cases\ZhanJiang\BASE_newC\BASE_newC.PRT"
##CLOUDTRANSPATH=r"C:\Users\shixi\Desktop\trans.test"
##NX=72
##NY=32
##NZ=112

#CLOUDTRANSPATH=r"C:\Users\shixi\Desktop\Cloud_VS2017_SC\CLOUD2017\transCloud.txt"
#PRTPATH = r"C:\Users\shixi\Desktop\cases\big\KZ0_ecl\KZ0_ECL.PRT"
#NX=171
#NY=77
#NZ=12
##全局变量
IJK_SEC_FLAG=r'(I,  J,  K)'
PROP_END_FLAG=r'1                             **********************************************************************'
import matplotlib.pyplot as plt
####################### One IJK Sec Read Function ########## 
def readIJKSection(ijksec,IIDX,JIDX,KIDX,VAL,offset):
    r"传入从 ...(I,  J,  K) I=  1... 开始的list，直到此段结束，将结果按I,J,K,Val存入指定位置 \
    其中ijksec是文件内容按行的list，offset是几个数组的写入起始位置"
    ijksec[0] = [ijksec[0][:15],ijksec[0][15:].split()]
    ijksec[1:] = [[ll[:12],ll[12:].split()]  for ll in ijksec[1:]]
    IRange=[0 for i in range(15)]
    idx=0
    nCol=0
    #print(ijksec[0])
    #print(ijksec[1])
    for i in ijksec[0][1]: #解析I下标,对应列
        IRange[nCol]=int(i)
        nCol+=1
    #print(nCol)
    nN = nCol*(len(ijksec)-2)
    #print(nN)
    idx=offset
    for ll in ijksec[2:]: #对每一数据行,得到JK下标
        jj = int(ll[0].split(',')[1])
        kk = int(ll[0].split(',')[2][:-1])
        nCol=0
        for i in ll[1]:
            if '--'  in i:
                #VAL[idx]=float('nan')
                VAL[idx]=0
            else:
                VAL[idx]=float(i)
            JIDX[idx]=jj
            KIDX[idx]=kk
            IIDX[idx]=IRange[nCol]
            nCol+=1
            idx+=1
    #print(idx)
    return nN
####################### Function End ############

####################### One Property multi-IJK Read Function ########## 
def readProperty(LineList,PropStartFlag,IIDX,JIDX,KIDX,VAL):
    r"读取某prt文件LineList中的属性，注意PropStartFlag要有唯一性，用于定位开始位置 "
    lineNum=0
    for ll in LineList:
        if PropStartFlag in ll:
            break
        else:
            lineNum+=1

    if lineNum==len(LineList) :
        print("Make sure ",PropStartFlag," is in the .RPT file!")
        return -1

    lineNum+=2
    ijkSl=0
    ijkEl=0
    findStartIJK=False
    findEndIJK=False
    breakflag=False
    offset=0
    while lineNum<len(LineList) and breakflag==False:
        lineNum+=1
        if lineNum==len(LineList):
            if findStartIJK==True and findEndIJK==False:
                findEndIJK==True
                ijkEl=lineNum
        else:
            ll=LineList[lineNum]
            if PROP_END_FLAG in ll:
                if findStartIJK==True and findEndIJK==False:
                    findEndIJK=True
                    ijkEl=lineNum
                breakflag=True
            elif IJK_SEC_FLAG in ll and findStartIJK==False:
                ijkSl=lineNum
                findStartIJK=True
                lineNum+=1
            elif ll[0]=='\n' and findStartIJK==True and findEndIJK==False:
                ijkEl=lineNum
                findEndIJK=True

        if findStartIJK==True and findEndIJK==True:
            findStartIJK=False
            findEndIJK=False
            tmpi = readIJKSection(LineList[ijkSl:ijkEl],IIDX,JIDX,KIDX,VAL,offset)
            offset+=tmpi
            #print(' ',tmpi,' entries read. ',offset,' entries read in total.')
            
    return 0
####################### Function End ############

########## ============NNC Reading ==========###########
def readNNC(LineList): 
    r"nnc flag 后第12行是正文，以 --- 结束 "
    NNCFlag_ = 'ALL NNCS AT      0.00  DAYS'
    NNC_END_Flag_ = r'---'
    lineNum=0
    for ll in LineList:
        if NNCFlag_ in ll:
            break
        else:
            lineNum+=1

    if lineNum==len(LineList) :
        print("Make sure ",NNCFlag_," is in the .RPT file!")
        return 0,0,0

    lineNum+=12
    ijkSl=lineNum
    for ll in LineList[ijkSl:]:
        if NNC_END_Flag_ in ll:
            break
        else:
            lineNum+=1

    ijkEl=lineNum

    I1 = [0 for i in range(ijkSl,ijkEl)]
    J1 =[0 for i in range(ijkSl,ijkEl)]
    K1=[0 for i in range(ijkSl,ijkEl)]
    I2=[0 for i in range(ijkSl,ijkEl)]
    J2=[0 for i in range(ijkSl,ijkEl)]
    K2=[0 for i in range(ijkSl,ijkEl)]
    VAL= [0.0 for i in range(ijkSl,ijkEl)]
    IJK1=[0 for i in range(ijkSl,ijkEl)]
    IJK2=[0 for i in range(ijkSl,ijkEl)]
    
    idx=0
    for ll in LineList[ijkSl:ijkEl]:
        ll=ll.split()
        I1[idx]=int(ll[0])
        J1[idx]=int(ll[1])
        K1[idx]=int(ll[2])
        I2[idx]=int(ll[3])
        J2[idx]=int(ll[4])
        K2[idx]=int(ll[5])
        VAL[idx]=float(ll[6])
        IJK1[idx]=(K1[idx]-1)*NX*NY+(J1[idx]-1)*NX+(I1[idx]-1)
        IJK2[idx]=(K2[idx]-1)*NX*NY+(J2[idx]-1)*NX+(I2[idx]-1)
        idx+=1
            
    return IJK1,IJK2,VAL

####################### Function End ############

## TRANS XYZ READ#####
def DirTransRead(RPTF,tranDir,IIDX,JIDX,KIDX,VAL,IJKIDX):
    if tranDir=='X':
        PropFlag = "TRANX    AT      0.00  DAYS"
    elif tranDir=='Y':
        PropFlag = "TRANY    AT      0.00  DAYS"
    elif tranDir=='Z':
        PropFlag = "TRANZ    AT      0.00  DAYS"
    else:
        print('Trans dirction error')
        exit(0)

    flag = readProperty(RPTF,PropFlag,IIDX,JIDX,KIDX,VAL)
    if flag<0:
        exit(0)

    for i in range(NALLCELL):
        IJKIDX[i]=(KIDX[i]-1)*NX*NY+(JIDX[i]-1)*NX+(IIDX[i]-1)
    
    print('TRANS ',tranDir,' RPT Read Complete')      
    return 1


### DIFF PLOT ######
def DiffPlot(picHandle,p2,tranDir,CloudTransDic,IJKIDX,VAL,IJKIDX2=0):
    r'计算差异,对于每个ecl的trans值，找cloud对应.xyz ecl只有IJKIDX，nnc ecl另有IJKIDX2'
    if tranDir=='X':
        di=1
    elif tranDir=='Y':
        di=NX
    elif tranDir=='Z':
        di = NX*NY 
    else:
        di = -9999

    tmpidx=0
    difflist=[0 for i in IJKIDX]
    tmpidx2=0
    difflist2=[0 for i in IJKIDX]
    for ijk in range(0,len(IJKIDX)):
        if ijk%10000==0:
            print('.',end='')
        tecl=VAL[ijk]
        if di>=0:
            keyi = (IJKIDX[ijk],IJKIDX[ijk]+di)
            keyi2 = (IJKIDX[ijk]+di,IJKIDX[ijk])
        else:
            keyi = (IJKIDX[ijk],IJKIDX2[ijk])
            keyi2 = (IJKIDX2[ijk],IJKIDX[ijk])
        if keyi in CloudTransDic:
            tcloud=CloudTransDic[keyi]
        elif keyi2 in CloudTransDic:
            tcloud=CloudTransDic[keyi2]
        else:
            tcloud=0
            
        if tecl<1e-5 :
            if tcloud<1e-5:
                continue
            else:
                difflist2[tmpidx2]=1
                tmpidx2+=1
                continue
        elif tcloud<1e-5 :
            if tecl<1e-5:
               continue
            else:
                difflist2[tmpidx2]=-1
                tmpidx2+=1
                continue
        else:
            difflist[tmpidx]=(tcloud-tecl)/tecl
        if tranDir=='X':
            maxDiff=5
        elif tranDir=='Y':
            maxDiff=5
        elif tranDir=='Z':
            maxDiff=5
        else:
            maxDiff=5
        maxDiff=0.5
        if abs(difflist[tmpidx])>maxDiff :
            difflist[tmpidx]=difflist[tmpidx]/abs(difflist[tmpidx])*maxDiff
            #print('IJK ',ijk,' I J K ',IIDX[ijk],' ',JIDX[ijk],' ',KIDX[ijk], ' TECL ',tecl,' TC ',tcloud)
        tmpidx+=1

    print('Diff Calculating End')
            
    #直方图
    picHandle.set_ylabel('NUMBER')
    if tranDir=='X':
        picHandle.set_xlabel('TRANX RELDIFF(+ means C>E)')  
    elif tranDir=='Y':
        picHandle.set_xlabel('TRANY RELDIFF(+ means C>E)') 
    elif tranDir=='Z':
        picHandle.set_xlabel('TRANZ RELDIFF(+ means C>E)') 
    else:
        picHandle.set_xlabel('TRAN NNC RELDIFF') 
    #picHandle.set_xticks([-5,-4,-3,-2,-1,0,1,2,3,4,5])
   
    picHandle.hist(
        difflist[:tmpidx],      #要统计的数据   
        #bins=[-6,-5,-4,-3,-2,-1,-0.5,-0.1,0.1,0.5,1,2,3,4,5,6],      #柱的个数 -1:0.01 1%误差以内，50%
        #bins=[-1.2,-1,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1],
        bins=[-0.5,-0.45	,-0.4,	-0.35,	-0.3,	-0.25,	-0.2,	-0.15,	-0.1,	-0.05,	-0.01,	0,	0.01,	0.05,	0.1,	0.15,	0.2	,0.25,	0.3,	0.35,	0.4,	0.45,0.5],

		#range =None,      #默认即可
        #normed=True,    #true为概率密度
        #stacked=True    #频率
		#log=True,
		histtype='bar',
		rwidth=0.8
    )

    p2.set_ylabel('NUMBER')
    if tranDir=='X':
        p2.set_xlabel('TRANX DISMATCH(+ means C>E)')  
    elif tranDir=='Y':
        p2.set_xlabel('TRANY DISMATCH(+ means C>E)') 
    elif tranDir=='Z':
        p2.set_xlabel('TRANZ DISMATCH(+ means C>E)') 
    else:
        p2.set_xlabel('TRAN NNC DISMATCH') 
    p2.set_xticks([-1,1])
   
    p2.hist(
        difflist2[:tmpidx2],      #要统计的数据   
        bins=[-1.5,0,1.5],      #柱的个数 -1:0.01 1%误差以内，50%
		#range =None,      #默认即可
        #normed=True,    #true为概率密度
        #stacked=True    #频率
		#log=True,
		histtype='bar',
		rwidth=0.8
    )
    #p2=plt.plot(range(0,tmpidx),difflist[0:tmpidx])
    #plt.show()


    # 柱状图 bar/barh
    #rects1=plt.bar(                      #(x,data) 就是所要画的二维数据
     #       e100_idxArray,                      #x 是X坐标轴数据，即每个块的x轴起始位置
      #      e100_valArray,                 #data是Y坐标轴的数据，即每个块的y轴高度
            #width=[0.1,0.2,0.3],         #每一个块的显示宽度
            #bottom=[1,2,3],              #每一个块的底部高度
    ##        color='y',                   #块的颜色
    ##        edgecolor='g',               #块的边界颜色
    ##        linewidth=2,                 #块的线条宽度
    ##        xerr=1,                      #x轴误差bar
    ##        yerr=1,                      #y轴误差bar
    ##        ecolor='r',                  #误差bar的颜色
    ##        capsize=1,                   #误差bar的线条宽度
    ##        orientation='vertical',     #块的方向  (horizontal,vertical)
    ##        align="center",              #块的位置 (center, left, right)
    ##        hold=None
         #   )
    return 1

###==========MAIN LOOP============
#读cloud文件
ftmp = open(CLOUDTRANSPATH)
flines=ftmp.readlines()
flines=[ll.split() for ll in flines]

Aidx = [int(ll[0]) for ll in flines[:]]
Bidx = [int(ll[1]) for ll in flines[:]]
Tran = [float(ll[2]) for ll in flines[:]]
A_B = [0.0 for ll in flines[:]]
for i in range(len(Aidx)):
    A_B[i]=(Aidx[i],Bidx[i])

CloudTransDic = dict(zip(A_B,Tran))
print('Cloud connList Read Complete')
ftmp.close()

#读rpt文件
ftmp = open(PRTPATH)
flines = ftmp.readlines()

#用于存储xyz trans
NALLCELL=NX*NY*NZ
IIDX=[-1 for i in range(NALLCELL)]
JIDX=[-1 for i in range(NALLCELL)]
KIDX=[-1 for i in range(NALLCELL)]
VAL=[-1 for i in range(NALLCELL)]
IJKIDX=[0 for i in range(NALLCELL)]

fig,axes = plt.subplots(2,2)
fig1,axes1 = plt.subplots(2,2)
#xyz trans read
DirTransRead(flines,'X',IIDX,JIDX,KIDX,VAL,IJKIDX)
DiffPlot(axes[0,0],axes1[0,0],'X',CloudTransDic,IJKIDX,VAL)
#将ecl trans写为cloud可识别的CONN
fconnOut=open(r"ECL_CONN",'w')
nConns=0
writeOutList=[0 for i in range(NALLCELL*4)]
for i in range(NALLCELL):
    if VAL[i]>0:
        writeOutList[nConns]=str(IJKIDX[i])+' '+str(IJKIDX[i]+1)+' '+str(VAL[i])+'\n'
        nConns=nConns+1
    
DirTransRead(flines,'Y',IIDX,JIDX,KIDX,VAL,IJKIDX)
DiffPlot(axes[0,1],axes1[0,1],'Y',CloudTransDic,IJKIDX,VAL)
#将ecl trans写为cloud可识别的CONN
for i in range(NALLCELL):
    if VAL[i]>0:
        writeOutList[nConns]=str(IJKIDX[i])+' '+str(IJKIDX[i]+NX)+' '+str(VAL[i])+'\n'
        nConns=nConns+1

DirTransRead(flines,'Z',IIDX,JIDX,KIDX,VAL,IJKIDX)
DiffPlot(axes[1,0],axes1[1,0],'Z',CloudTransDic,IJKIDX,VAL)
#将ecl trans写为cloud可识别的CONN
for i in range(NALLCELL):
    if VAL[i]>0:
        writeOutList[nConns]=str(IJKIDX[i])+' '+str(IJKIDX[i]+NX*NY)+' '+str(VAL[i])+'\n'
        nConns=nConns+1
        
#nnc trans
IJKIDX,IJKIDX2,VALNNC = readNNC(flines)
if IJKIDX==0 and IJKIDX2==0 and VALNNC==0:
    print("NO nnc\n")
else:
    DiffPlot(axes[1,1],axes1[1,1],'N',CloudTransDic,IJKIDX,VALNNC,IJKIDX2)
    #将ecl trans写为cloud可识别的CONN
    for i in range(len(IJKIDX)):
        if VALNNC[i]>0:
            writeOutList[nConns]=str(IJKIDX[i])+' '+str(IJKIDX2[i])+' '+str(VALNNC[i])+'\n'
            nConns=nConns+1
    fconnOut.write("CONN\n")
    fconnOut.write(str(nConns)+" 2\n") #2代表ecl输入的conn
    fconnOut.writelines(writeOutList[:nConns])
    fconnOut.write(r"/")
    fconnOut.close()               

fig.tight_layout()
fig1.tight_layout()
plt.show()   
print('END')
