#!/usr/bin/python3

#libraries
import xlrd, re
import glob
import os
import string
import sys,random
import numpy,scipy, statistics,math
import argparse
import time
from matplotlib import rcParams
rcParams['font.family']='monospace'
import matplotlib.pyplot as plt

class SubjectData():
    def __init__(self,filename):
        self.filename=filename
        self.name=os.path.basename(filename)
        self.directory=os.path.dirname(filename)
        print(self.name.replace(".","_").split("_")[0])
        self.SubjectName=self.name.replace(".","_").split("_")[0]
        self.readFile()
    def readFile(self):
        self.data=dict()
        self.ratio=dict()
        self.chi2=dict()
        self.siteChi2=dict()
        with open(self.filename,"r") as Fname:
            for line in Fname:
                if "mutation" in line:
                    continue
                splitline=line.split()
                self.data[splitline[0]]=splitline
                self.ratio[splitline[0]]=float(splitline[3])
                self.chi2[splitline[0]]=self.calcChi2(splitline[1],splitline[2])
                if splitline[0][0]==splitline[0][-1]:
                    self.siteChi2[splitline[0]]=self.calcChi2(1.0-float(splitline[1]),1.0-float(splitline[2]))

                
    def getMutation(self,cutoff,obsFreqcutoff):
        goodMuts=[]
        vals=[]
        for mut in self.ratio.keys():
            if abs(self.ratio[mut])>cutoff and math.inf!=self.ratio[mut] and float(self.data[mut][2])>obsFreqcutoff:
                goodMuts.append(mut)
                vals.append(abs(self.ratio[mut]))

        return goodMuts,vals
    def calcChi2(self,expected,observed):
        exp=float(expected)
        obs=float(observed)
        nom=(obs-exp)*(obs-exp)
        if obs==0.0 and exp==0.0:
            results=0.0
        elif exp==0.0 and obs!=0.0:
            results=0.0
        elif  exp>obs:
            try:
                results=-float(nom)/float(exp)
            except:
                print("error in calculating Chi2")
                results=0.0
        else:
            try:
                results=float(nom)/float(exp)
            except:
                print("error in calculating Chi2")
                results=0.0

        return results

def write_mutations_per_mouse(allData):
    selection_cutoff=1.0
    for groups in allData.keys():
        for i,subj in  enumerate(sort_keys(allData[groups].keys())):
            with open(os.path.join("MutationLists","MutationList.cutoff{}.{}.txt").format(selection_cutoff,subj),'w') as fopen:
                for mut in allData[groups][subj].data.keys():
                    if float(allData[groups][subj].data[mut][2])<=0.01:
                        continue
                    if float(allData[groups][subj].data[mut][1])<=0.0:
                        continue
                    if float(allData[groups][subj].data[mut][2])/float(allData[groups][subj].data[mut][1])>=selection_cutoff:
                        fopen.write("{} {:0.2f}\n".format(mut,float(allData[groups][subj].data[mut][2])/float(allData[groups][subj].data[mut][1])))
    
    
def plotChi2graph(allData,graphname,mutationlist,chitype="chi2"):
    subj_labels=[]
    tmp=[]
    for k in allData.keys():
        for x in allData[k].keys():
            tmp.append(x)
    fig,ax=plt.subplots(figsize=(len(tmp),len(mutationlist)/2))
    xpos=0
    for groups in allData.keys():
        for i,subj in  enumerate(sort_keys(allData[groups].keys())):
            subj_labels.append(subj)
            for n,mut in enumerate(mutationlist):
                if numpy.isnan(allData[groups][subj].__dict__[chitype][mut]):
                    ax.plot(xpos,n,'kx',ms=0)
                elif allData[groups][subj].__dict__[chitype][mut]==math.inf:
                    ax.plot(xpos,n,'rx',ms=0)
                elif allData[groups][subj].__dict__[chitype][mut]==0.0:
                    ax.plot(xpos,n,'wo',ms=0,mec='r')
                elif allData[groups][subj].__dict__[chitype][mut]<0:
                    ax.plot(xpos,n,'bo',ms=round(abs(math.log2(abs(allData[groups][subj].__dict__[chitype][mut])))))
                else:
                    #print(math.log10(allData[groups][subj].chi2[mut]))
                    ax.plot(xpos,n,'ro',ms=round(abs(math.log2(allData[groups][subj].__dict__[chitype][mut]))))
            #print(xpos)
            xpos+=1
        if xpos<len(tmp):
            ax.axvline(xpos-0.5,color='black',ls='--')
    ax.set_xticks(range(len(subj_labels)))
    ax.set_xticklabels(subj_labels)
    ax.set_yticks(range(len(mutationlist)))
    ax.set_yticklabels(mutationlist)
    ax.set_xlim([-0.5,len(subj_labels)-0.5])
    ax.set_ylim([-1,len(mutationlist)])
    plt.xlabel("Mice Labels")
    plt.ylabel("Mutations")
    plt.savefig(graphname,bbox_inches='tight')
    plt.close()
    
def plotgraph(allData,graphname,mutationlist):
    subj_labels=[]
    tmp=[]
    for k in allData.keys():
        for x in allData[k].keys():
            tmp.append(x)
    fig,ax=plt.subplots(figsize=(len(tmp),len(mutationlist)/4))
    xpos=0
    for groups in allData.keys():
        #subjs=[]
        #if "GT1.2" in groups:
        #    subjs=['V11411','10504','V11421', 'V11417', '10502', 'V11412', '10503','V11425']
        #else:
        #    subjs=['V11405', 'V11407', 'V11404']
        for i,subj in  enumerate(sort_keys(allData[groups].keys())):
        #for i,subj in enumerate(subjs):
            subj_labels.append(subj)
            for n,mut in enumerate(mutationlist):
                if numpy.isnan(allData[groups][subj].ratio[mut]):
                    ax.plot(xpos,n,'kx',ms=0)
                elif allData[groups][subj].ratio[mut]==math.inf:
                    ax.plot(xpos,n,'rx',ms=0)
                elif allData[groups][subj].ratio[mut]==0.0:
                    ax.plot(xpos,n,'wo',ms=0,mec='r')
                elif allData[groups][subj].ratio[mut]<0:
                    ax.plot(xpos,n,'bo',ms=2*round(abs(allData[groups][subj].ratio[mut])))
                else:
                    ax.plot(xpos,n,'ro',ms=2*round(abs(allData[groups][subj].ratio[mut])))
            #print(xpos)
            xpos+=1
        if xpos<len(tmp):
            ax.axvline(xpos-0.5,color='black',ls='--')
    ax.set_xticks(range(len(subj_labels)))
    ax.set_xticklabels(subj_labels)
    ax.set_yticks(range(len(mutationlist)))
    ax.set_yticklabels(mutationlist)
    ax.set_xlim([-0.5,len(subj_labels)-0.5])
    ax.set_ylim([-1,len(mutationlist)])
    plt.xlabel("Mice Labels")
    plt.ylabel("Mutations")
    plt.savefig(graphname,bbox_inches='tight')
    plt.close()

def plotgraphRatioV2(allData,graphname,mutationlist):
    subj_labels=[]
    tmp=[]
    for k in allData.keys():
        for x in allData[k].keys():
            tmp.append(x)
    fig,ax=plt.subplots(figsize=(len(tmp),10))
    xpos=0
    counts=dict()
    for groups in allData.keys():
        for i,subj in  enumerate(sort_keys(allData[groups].keys())):
            subj_labels.append(subj)
            counts[subj]=0
            for n,mut in enumerate(mutationlist):
                #print(allData[groups][subj].data[mut])
                if float(allData[groups][subj].data[mut][2])<=0.01:
                    continue
                if float(allData[groups][subj].data[mut][1])<=0.0:
                    continue
                if float(allData[groups][subj].data[mut][2])<float(allData[groups][subj].data[mut][1]):
                    ratio=(1-float(allData[groups][subj].data[mut][2]))/(1-float(allData[groups][subj].data[mut][1]))
                    ax.plot(xpos,ratio,'ro')
                    counts[subj]=counts[subj]+1

            #print(xpos)
            xpos+=1
        if xpos<len(tmp):
            ax.axvline(xpos-0.5,color='black',ls='--')

    ax.set_xticks(range(len(subj_labels)))
    ax.set_xticklabels(subj_labels)
    ax.set_xlim([-0.5,len(subj_labels)-0.5])
    ax2=ax.twinx() 
    ax2.set_ylabel('count of positive selection',color='blue')
    ax2.plot(counts.values(),'-sb',mfc='none',ms=10)
    #ax.set_yscale('log')
    plt.xlabel("Mice Labels")
    ax.set_ylabel("(1-Observed)/(1-Expected)")
    plt.savefig(graphname,bbox_inches='tight')
    plt.close()

def plotgraphSimpRatio(allData,graphname,mutationlist=[]):
    subj_labels=[]
    random.seed(20220215)#sets the seed to lock in the jiggle positions
    tmp=[]
    for k in allData.keys():
        for x in allData[k].keys():
            tmp.append(x)
    fig,ax=plt.subplots(figsize=(len(tmp),6))
    xpos=0
    counts=dict()
    mutcount=dict()
    grouplines=[0,0]
    for groups in allData.keys():
        
        for i,subj in  enumerate(sort_keys(allData[groups].keys())):
            medianValue=[]
            subj_labels.append(subj)
            counts[subj]=0
            mutcount[subj]=0
            for n,mut in enumerate(allData[groups][subj].data.keys()):
                sel=[]
                for g in allData.keys():
                    for s in allData[g].keys():
                        sel.append(allData[g][s].data[mut][2])
                if float(max(sel))<=0.01:
                    continue
                if float(allData[groups][subj].data[mut][1])<=0.0:
                    continue
                if float(allData[groups][subj].data[mut][2])>float(allData[groups][subj].data[mut][1]):
                    ratio=float(allData[groups][subj].data[mut][2])/float(allData[groups][subj].data[mut][1])
                    medianValue.append(ratio)
                    grouplines[1],=ax.plot(xpos+random.random()/4-0.15,ratio,'ro',alpha=0.75)
                    counts[subj]=counts[subj]+1
            counts[subj]=counts[subj]/float(len(allData[groups][subj].data.keys()))
            ax.plot([xpos-0.4,xpos+0.4],[numpy.median(medianValue),numpy.median(medianValue)],color='black',ls='-',zorder=1)
            #print(numpy.median(medianValue))
            ax.axvline(xpos,color='grey',ls=':',zorder=0)
            xpos+=1
        if xpos<len(tmp):
            ax.axvline(xpos-0.5,color='black',ls='--')
    ax.set_xticks(range(len(subj_labels)))
    ax.set_xticklabels(subj_labels)
    ax.set_xlim([-0.5,len(subj_labels)-0.5])
    #ax.set_ylim([0.0001,10000])
    ax.set_yscale('log')
    ax.yaxis.grid(True)
    plt.xlabel("Mice Labels")
    ax.set_ylabel("Observed/Expected")

    plt.savefig(graphname,bbox_inches='tight')
    plt.close()

def plotgraphRatio(allData,graphname,fullmutationlist=[]):
    subj_labels=[]
    tmp=[]
    for k in allData.keys():
        for x in allData[k].keys():
            tmp.append(x)
    fig,ax=plt.subplots(figsize=(len(tmp),6))
    xpos=0
    counts=dict()
    goodcounts=dict()
    mutcount=dict()
    allmutcounts=dict()
    
    grouplines=[0,0]
    contactdots=[]
    
    mutationlist=[x[:-1] for x in fullmutationlist]
    for groups in allData.keys():
        for i,subj in enumerate(sort_keys(allData[groups].keys())):
            subj_labels.append(subj)
            counts[subj]=0
            mutcount[subj]=0
            allmutcounts[subj]=0
            goodcounts[subj]=0
            for n,mut in enumerate(allData[groups][subj].data.keys()):

                if float(allData[groups][subj].data[mut][2])<=0.01:
                    continue
                if float(allData[groups][subj].data[mut][1])<=0.0:
                    continue
                if mut[:-1] in mutationlist:
                    allmutcounts[subj]=allmutcounts[subj]+1
                goodcounts[subj]=goodcounts[subj]+1
                if float(allData[groups][subj].data[mut][2])>float(allData[groups][subj].data[mut][1]):
                    ratio=float(allData[groups][subj].data[mut][2])/float(allData[groups][subj].data[mut][1])
                    counts[subj]=counts[subj]+1
                    if mut[:-1] in mutationlist:
                        contactdots.append((xpos,ratio))
                        #grouplines[0],=ax.plot(xpos,ratio,'bo')
                        mutcount[subj]=mutcount[subj]+1
                    else:
                        grouplines[1],=ax.plot(xpos,ratio,'ro',alpha=0.5)

            #counts[subj]=counts[subj]/float(len(allData[groups][subj].data.keys()))
            counts[subj]=counts[subj]/float(goodcounts[subj])
            xpos+=1
        if xpos<len(tmp):
            ax.axvline(xpos-0.5,color='black',ls='--')
    for dot in contactdots:
        grouplines[0],=ax.plot(dot[0],dot[1],'bo')
    ax.set_xticks(range(len(subj_labels)))
    ax.set_xticklabels(subj_labels)
    ax.set_xlim([-0.5,len(subj_labels)-0.5])
    ax.set_yscale('log')
    ax2=ax.twinx()
    ax2.set_ylabel('count of positive selection/Number of Mutations',color='blue')
    print("All mutations percentages")
    print(["{:0.4f}".format(float(counts[x])) for x in subj_labels])
    ax2.plot([counts[x] for x in subj_labels],'k-s',ms=10,label="All Mutations",mfc='none')
    if mutationlist:
        print("contact percentages")
        print(["{:0.4f}".format(mutcount[x]/float(allmutcounts[x])) for x in subj_labels])
        ax2.plot([mutcount[x]/float(allmutcounts[x]) for x in subj_labels],'b:s',ms=10,label="Contact Positions",mfc='none')
    ax2.legend(bbox_to_anchor=(0,1,1,0),loc="lower right", ncol=2)
    ax.legend(grouplines,["Contact Positions","All Positions"],bbox_to_anchor=(0,1,1,0),loc="lower left", ncol=2)
    plt.xlabel("Mice Labels")
    ax.set_ylabel("Observed/Expected")

    plt.savefig(graphname,bbox_inches='tight')
    plt.close()
    
def plotgraphRatioMutations(allData,graphname,mutationlist):
    subj_labels=[]
    tmp=[]
    for k in allData.keys():
        for x in allData[k].keys():
            tmp.append(x)
    fig,ax=plt.subplots(figsize=(len(tmp),6))
    xpos=0
    counts=dict()
    for groups in allData.keys():
        for i,subj in  enumerate(sort_keys(allData[groups].keys())):
            subj_labels.append(subj)
            counts[subj]=0
            for n,mut in enumerate(mutationlist):
                if float(allData[groups][subj].data[mut][2])<=0.01:
                    continue
                if float(allData[groups][subj].data[mut][1])<=0.0:
                    continue
                if float(allData[groups][subj].data[mut][2])>float(allData[groups][subj].data[mut][1]):
                    ratio=float(allData[groups][subj].data[mut][2])/float(allData[groups][subj].data[mut][1])
                    counts[subj]=counts[subj]+1
                    ax.plot(xpos,ratio,'ro')
            #print(xpos)
            xpos+=1
        if xpos<len(tmp):
            ax.axvline(xpos-0.5,color='black',ls='--')

    ax.set_xticks(range(len(subj_labels)))
    ax.set_xticklabels(subj_labels)
    ax.set_xlim([-0.5,len(subj_labels)-0.5])
    ax.set_yscale('log')
    ax2=ax.twinx()
    ax2.set_ylabel('Number of positive selection/Total contact mutations',color='blue')
    ax2.plot([counts[x]/float(len(mutationlist)) for x in subj_labels],'b:s',ms=7)
    t=["{:0.5f},{}".format(float(counts[x])/float(len(mutationlist)),x) for x in subj_labels]
    print("\n".join(t))
    plt.xlabel("Mice Labels")
    ax.set_ylabel("Observed/Expected")

    plt.savefig(graphname,bbox_inches='tight')
    plt.close()
    
def plotgraphLineRatio(allData,graphname,AreaCutoff=0.6,label=False):
    #different symbols for different set of mutations

    subj_labels=[]
    tmp=[]
    for k in allData.keys():
        for x in allData[k].keys():
            tmp.append(x)
    fig,ax=plt.subplots(figsize=(6,6))
    xpos=0
    counts=dict()
    colors=[[1,0,0],[0,0,0],[0,0,1]]
    grouplines=[0,0,0,0,0,0]
    positiveSelection=[0 for x in range(len(tmp))]
    negativeSelection=[0 for x in range(len(tmp))]
    CMpositiveSelection=[0 for x in range(len(tmp))]
    CMnegativeSelection=[0 for x in range(len(tmp))]
    goodMuts=["T30P","Y33F","W47W","W50W","N52N","N54T","S55N","G56G","G57A","T58V","N59N","Y60Y","A61A","Q62W","K63Y","T69T","R72R","T74R","E106E","W107W"]
    goodMuts2=[x[:-1] for x in goodMuts]
    areaCount=0
    alphaValue=1.0
    if label:
        alphaValue=0.5

    for g,groups in enumerate(allData.keys()):
        for i,subj in  enumerate(sort_keys(allData[groups].keys())):
            subj_labels.append(subj)
            counts[subj]=0
            for n,mut in enumerate(allData[groups][subj].data.keys()):
                if float(allData[groups][subj].data[mut][2])<=0.01:
                    continue
                if float(allData[groups][subj].data[mut][1])<=0.0:
                    continue
                if mut[0]==mut[-1] and label:
                    ax.text(float(allData[groups][subj].data[mut][1]),float(allData[groups][subj].data[mut][2]),"{}-{}".format(subj,mut),size=3)
                if mut[0]==mut[-1] and not label:
                    continue
                if float(allData[groups][subj].data[mut][1]) <= AreaCutoff and float(allData[groups][subj].data[mut][2]) >= AreaCutoff:
                    areaCount+=1
                if mut[:-1] in goodMuts2:
                    grouplines[g],=ax.plot(float(allData[groups][subj].data[mut][1]),float(allData[groups][subj].data[mut][2]),'s',color=colors[g],ms=4,fillstyle='none',alpha=alphaValue)
                    if allData[groups][subj].data[mut][1]>allData[groups][subj].data[mut][2]:
                        CMnegativeSelection[xpos]=CMnegativeSelection[xpos]+1
                        negativeSelection[xpos]=negativeSelection[xpos]+1
                    else:
                        CMpositiveSelection[xpos]=CMpositiveSelection[xpos]+1
                        positiveSelection[xpos]=positiveSelection[xpos]+1
                else:
                    #pass
                    grouplines[g+3],=ax.plot(float(allData[groups][subj].data[mut][1]),float(allData[groups][subj].data[mut][2]),'o',color=colors[g],ms=2,fillstyle='none',alpha=alphaValue/2.0)
                    if allData[groups][subj].data[mut][1]>allData[groups][subj].data[mut][2]:
                        negativeSelection[xpos]=negativeSelection[xpos]+1
                    else:
                        positiveSelection[xpos]=positiveSelection[xpos]+1
            xpos+=1
        if xpos<len(tmp):
            ax.axvline(xpos-0.5,color='black',ls='--')
    print("Count in area {}:{}".format(AreaCutoff,areaCount))
    ax.plot([-0.01,1.01],[-0.01,1.01],'-k',zorder=1)
    #ax.set_yscale('log')
    #ax.set_xscale('log')
    ax.set_ylim([-0.01,1.01])
    ax.set_xlim([-0.01,1.01])
    ax.legend(grouplines,["Cont. Pos.\nBG505","Cont. Pos.\nGT1.2 short","Cont. Pos.\nGT1.2 long","all BG505","all GT1.2 short","all GT1.2 long"],bbox_to_anchor=(1.02,1),loc='upper left')
    
    plt.ylabel("Observed Frequency")
    plt.xlabel("Expected Frequency")

    plt.savefig(graphname,bbox_inches='tight')
    plt.close()
    print("positive\t{}".format(positiveSelection))
    print("negative\t{}".format(negativeSelection))
    print("Cont. positive\t{}".format(CMpositiveSelection))
    print("Cont. negative\t{}".format(CMnegativeSelection))
    
def plotgraphViolin(allData,graphname):
    #different symbols for different set of mutations
    subj_labels=[]
    tmp=[]
    for k in allData.keys():
        for x in allData[k].keys():
            tmp.append(x)

    fig,ax=plt.subplots(figsize=(len(tmp),6))
    xpos=0
    counts=dict()

    goodMuts=["T30P","Y33F","W47W","W50W","N52N","N54T","S55N","G56G","G57A","T58V","N59N","Y60Y","A61A","Q62W","K63Y","T69T","R72R","T74R","E106E","W107W"]
    colors=[[1,0,0],[0,0,1]]
    Plotdata=[]
    for g,groups in enumerate(allData.keys()):
        for i,subj in  enumerate(sort_keys(allData[groups].keys())):
            subj_labels.append(subj)
            counts[subj]=0
            dataset=[]
            for n,mut in enumerate(allData[groups][subj].data.keys()):
                sel=[]
                for g in allData.keys():
                    for s in allData[g].keys():
                        sel.append(allData[g][s].data[mut][2])
                if float(max(sel))<=0.01:
                    continue
                
                #if float(allData[groups][subj].data[mut][2])<=0.01:
                #    continue
                if float(allData[groups][subj].data[mut][1])<=0.0:
                    continue
                if float(allData[groups][subj].data[mut][2])>float(allData[groups][subj].data[mut][1]):
                    dataset.append(float(allData[groups][subj].data[mut][2])/float(allData[groups][subj].data[mut][1]))
            #Plotdata.append([x/float(len(dataset)) for x in dataset])
            Plotdata.append(dataset)
            #print("{} {}".format(len(dataset),numpy.mean(dataset)))
            xpos+=1
        if xpos<len(tmp):
            ax.axvline(xpos-0.5+1,color='black',ls='--')
    ax.violinplot(Plotdata, points=300, widths=0.5,showmeans=False, showextrema=True, showmedians=True)

    ax.set_xticks(range(1,len(subj_labels)+1))
    ax.set_xticklabels(subj_labels)
    ax.set_xlim([-0.5+1,len(subj_labels)-0.5+1])
    ax.set_yscale('log')
    ax.yaxis.grid(True)
    
    plt.xlabel("Mice Labels")
    ax.set_ylabel("Observed/Expected")

    plt.savefig(graphname,bbox_inches='tight')
    plt.close()

def plotgraphCombineViolin(allData,graphname):
    #different symbols for different set of mutations
    random.seed(20220215)#sets the seed to lock in the jiggle positions
    subj_labels=[]
    tmp=[]
    for k in allData.keys():
        for x in allData[k].keys():
            tmp.append(x)

    fig,ax=plt.subplots(figsize=(len(tmp),6))
    xpos=0
    counts=dict()

    goodMuts=["T30P","Y33F","W47W","W50W","N52N","N54T","S55N","G56G","G57A","T58V","N59N","Y60Y","A61A","Q62W","K63Y","T69T","R72R","T74R","E106E","W107W"]
    colors=[[1,0,0],[0,0,1]]
    Plotdata=[]
    for g,groups in enumerate(allData.keys()):
        for i,subj in  enumerate(sort_keys(allData[groups].keys())):
            medianValue=[]
            subj_labels.append(subj)
            counts[subj]=0
            dataset=[]
            for n,mut in enumerate(allData[groups][subj].data.keys()):
                sel=[]
                for g in allData.keys():
                    for s in allData[g].keys():
                        sel.append(allData[g][s].data[mut][2])
                if float(max(sel))<=0.01:
                    continue
                
                #if float(allData[groups][subj].data[mut][2])<=0.01:
                #    continue
                if float(allData[groups][subj].data[mut][1])<=0.0:
                    continue
                #else:
                if float(allData[groups][subj].data[mut][2])>float(allData[groups][subj].data[mut][1]):
                    ratio=float(allData[groups][subj].data[mut][2])/float(allData[groups][subj].data[mut][1])
                    #dataset.append(float(allData[groups][subj].data[mut][2])/float(allData[groups][subj].data[mut][1]))
                    dataset.append(numpy.log10(ratio))
                    medianValue.append(ratio)
                    ax.plot(xpos+random.random()/4-0.15+1,ratio,'ro',alpha=0.5)
                    counts[subj]=counts[subj]+1

            counts[subj]=counts[subj]/float(len(allData[groups][subj].data.keys()))
            ax.plot([xpos-0.4+1,xpos+0.4+1],[numpy.median(medianValue),numpy.median(medianValue)],color='black',ls='-',zorder=1)
            #Plotdata.append([x/float(len(dataset)) for x in dataset])
            Plotdata.append(dataset)
            #Plotdata.append([random.gauss(random.random()*100,50)*10 for x in range(1000)])
            #print("{} {}".format(len(dataset),numpy.mean(dataset)))
            xpos+=1
        if xpos<len(tmp):
            ax.axvline(xpos-0.5+1,color='black',ls='--')

   

    #xpos=0
    #for g,groups in enumerate(allData.keys()):
    #    for i,subj in  enumerate(sort_keys(allData[groups].keys())):
    #        v_parts=ax.violinplot(points=300,positions=[xpos+1])#, showmeans=False, showextrema=True, showmedians=True)
    #        print(list(v_parts.keys()))
    #        for pc in v_parts['bodies']:
    #            pc.set_edgecolors('blue')
    #            pc.set_facecolors('blue')
    #        xpos+=1
    #input('e')
    ax2=ax.twinx()
    ax2.violinplot(Plotdata, points=300,widths=0.75,showmeans=False, showextrema=True, showmedians=True)
    ax2.set_yticks([])
    ax2.set_yticklabels([])
    ax.set_xticks(range(1,len(subj_labels)+1))
    ax.set_xticklabels(subj_labels)
    ax.set_xlim([-0.5+1,len(subj_labels)-0.5+1])
    ax.set_yscale('log')
    #ax.set_title(str(time.time()))
    #ax.set_ylim([0,10])
    #ax.set_xlim([10,12])
    ax.yaxis.grid(True)
    
    plt.xlabel("Mice Labels")
    ax.set_ylabel("Observed/Expected")

    plt.savefig(graphname,bbox_inches='tight')
    plt.close()
    
def writeOutData(allData,filename,mutationlist):
    mutationDict=dict()
    subj_labels=[]
    tmp=[]
    for groups in allData.keys():
        for i,subj in  enumerate(sort_keys(allData[groups].keys())):
            subj_labels.append(subj)
            for n,mut in enumerate(mutationlist):
                if mut not in mutationDict.keys():
                    mutationDict[mut]=[]
                if numpy.isnan(allData[groups][subj].ratio[mut]):
                    mutationDict[mut].append("NaN")
                elif allData[groups][subj].ratio[mut]==math.inf:
                    mutationDict[mut].append("inf")
                elif allData[groups][subj].ratio[mut]==0.0:
                    mutationDict[mut].append("0.000")
                elif allData[groups][subj].ratio[mut]<0:
                    mutationDict[mut].append("{:0.4f}".format(allData[groups][subj].ratio[mut]))
                else:
                    mutationDict[mut].append("{:0.4f}".format(allData[groups][subj].ratio[mut]))

    with open(filename,'w') as ofile:
        ofile.write("\t"+"\t".join(subj_labels)+"\n")
        for mut in mutationlist[::-1]:
            ofile.write("{}\t{}\t{}\n".format(mut,"\t".join(mutationDict[mut]),mut))
        ofile.write("\t"+"\t".join(subj_labels)+"\n")

def findException(allData,expectFreq,ratio):
    for groups in allData.keys():
        for i,subj in  enumerate(sort_keys(allData[groups].keys())):
            for mut in allData[groups][subj].data.keys():
                if float(allData[groups][subj].data[mut][1])>expectFreq and float(allData[groups][subj].data[mut][2])/float(allData[groups][subj].data[mut][1])>ratio:
                    print("{}\t{}\t{}\t{}".format(subj,mut,allData[groups][subj].data[mut][1],allData[groups][subj].data[mut][2]))

def sort_mutations(TMuts,maxN):
    sort_Muts_Tuple=sorted(TMuts)
    goodMuts=[]
    for entry in sort_Muts_Tuple[::-1]:
        #print(entry)
        if len(goodMuts)>=maxN:
            break
        else:
            if entry[1] not in goodMuts:
                goodMuts.append(entry[1])

    posMut=[(int(x[1:len(x)-1]),x) for x in goodMuts]
    sortMut=sorted(posMut)
    return [x[1] for x in sortMut[::-1]]

def sort_mutations_chi2(TMuts,maxN):
    sort_Muts_Tuple=sorted(TMuts)
    goodMuts=[]
    for entry in sort_Muts_Tuple[::-1]:
        #print(entry)
        if len(goodMuts)>=maxN:
            break
        else:
            if entry[1] not in goodMuts:
                goodMuts.append(entry[1])

    posMut=[(int(x[1:len(x)-1]),x) for x in goodMuts]
    sortMut=sorted(posMut)
    return [x[1] for x in sortMut[::-1]]

def write_out_percentages(allData,filename):
    goodMuts=["T30P","Y33F","W47W","W50W","N52N","N54T","S55N","G56G","G57A","T58V","N59N","Y60Y","A61A","Q62W","K63Y","T69T","R72R","T74R","E106E","W107W"]
    cutoff_measurement=2.0
    count_dict=dict()
    for groups in allData.keys():
        for i,subj in  enumerate(sort_keys(allData[groups].keys())):
            count_dict[subj]=[0,0]
            for n,mut in enumerate(allData[groups][subj].data.keys()):
                sel=[]
                for g in allData.keys():
                    for s in allData[g].keys():
                        sel.append(allData[g][s].data[mut][2])
                if float(max(sel))<=0.01:
                    continue
                if float(allData[groups][subj].data[mut][1])<=0.0:
                    continue
                ratio=float(allData[groups][subj].data[mut][2])/float(allData[groups][subj].data[mut][1])
                if ratio>2.0:
                    count_dict[subj][0]+=1
                    if mut in goodMuts and ratio>cutoff_measurement:
                        count_dict[subj][1]+=1

    with open(filename,"w") as f:
        f.write("\t#Muts\t#Cont.Muts\n")
        for groups in allData.keys():
            for i,subj in  enumerate(sort_keys(allData[groups].keys())):
                f.write("{}\t{}\t{}\n".format(subj,count_dict[subj][0],count_dict[subj][1]))


def parse_args():
    # Construct an argument parser
    all_args = argparse.ArgumentParser()
    
    # Add arguments to the parser
    all_args.add_argument("-d","--listdirectory",action="extend",nargs="+",required=False,help="directory with the list files")
    all_args.add_argument("-f","--listfile",action="extend",nargs="+",required=False,help="list files",default=[])
    all_args.add_argument("-o","--graphname",required=False,help="graph name to print out")
    all_args.add_argument("-v","--version",required=False,help="graph name to print out")
    all_args.add_argument("-c","--cutoff",default=2.0,required=False,help="cutoff value",type=float)
    all_args.add_argument("--observedcutoff",default=0.01,required=False,help="cutoff value",type=float)
    all_args.add_argument("--top",default=25,required=False,help="top X",type=int)
    args = vars(all_args.parse_args())
    return args

def sort_keys(keys):
    list_keys=list(keys)
    ordered_keys=["V11404","V11405","V11407","10504", "10502", "10503","V11411", "V11421", "V11417", "V11412", "V11425"]
    surviving_keys=["V11404","V11405","V11407","10504", "10502", "10503","V11411", "V11421", "V11417", "V11412", "V11425"]
    for key in ordered_keys:
        if key not in list_keys:
            surviving_keys.remove(key)

    return surviving_keys

def main(args):
    groups=dict()

    goodMut=set()
    values=set()
    subjects=dict()
    Muts_Tuple=[]
    if args["listfile"]:
        for listfile in args["listfile"]:
            subject=SubjectData(listfile)
            subjects[subject.SubjectName]=subject
        groups["individual"]=subjects

    if args["listdirectory"]:

        for subdir in args["listdirectory"]:
            subjects=dict()
            for file in glob.iglob(os.path.join(subdir,"*.txt")):
                subject=SubjectData(file)
                if subject.SubjectName not in subjects.keys():
                    subjects[subject.SubjectName]=subject
            groups[subdir]=subjects

    for group in groups.keys():
        for subject in groups[group].keys():
            muts,vals=groups[group][subject].getMutation(args["cutoff"],args["observedcutoff"])
            for i,mut in enumerate(muts):
                Muts_Tuple.append((vals[i],mut))
        sort_keys(groups[group].keys())

    goodMuts=sort_mutations(Muts_Tuple,args["top"])
    write_mutations_per_mouse(groups)

    print("Number of Mutations:{}".format(len(goodMuts)))
    #print("{} - {}".format(max(values),min(values)))
    plotgraph(groups,"log10ratio.top{}.{}.pdf".format(args["top"],args["version"]),goodMuts)
 
    writeOutData(groups,"log10ratio.top{}.{}.txt".format(args["top"],args["version"]),goodMuts)
    #findException(groups,0.01,2)

    print("Subject\tMutation\tlog10(ratio)\tchi2\texp\tobs")
    Chi2_Tuple=[]
    for group in groups.keys():
        for subject in groups[group].keys():
            for mut in groups[group][subject].chi2.keys():
                Chi2_Tuple.append((abs(groups[group][subject].chi2[mut]),mut))
                if "T58V" in mut:
                    print("{}\t{}\t{}\t{}\t{}\t{}".format(subject,mut,groups[group][subject].data[mut][3],groups[group][subject].chi2[mut],groups[group][subject].data[mut][1],groups[group][subject].data[mut][2]))

    print(len(Chi2_Tuple))
    goodMuts=sort_mutations_chi2(Chi2_Tuple,args["top"])

    print("Number of Mutations:{}".format(len(goodMuts)))
    plotChi2graph(groups,"Chi2graph.top{}.{}.pdf".format(args["top"],args["version"]),goodMuts)

    #contact mutations
    goodMuts=["T30P","Y33F","N54T","S55N","G57A","T58V","Q62W","K63Y","T74R"]
    plotChi2graph(groups,"chi2graph.ContactMutations.{}.pdf".format(args["version"]),goodMuts[::-1])

    goodMuts=["T30P","Y33F","W47W","W50W","N52N","N54T","S55N","G56G","G57A","T58V","N59N","Y60Y","A61A","Q62W","K63Y","T69T","R72R","T74R","E106E","W107W"]
    plotChi2graph(groups,"Chi2graph.ContactPositions.{}.pdf".format(args["version"]),goodMuts[::-1])
    
    plotgraphRatio(groups,"selectionGraph.{}.pdf".format(args["version"]),goodMuts)
    plotgraphSimpRatio(groups,"Simple.selectionGraph.{}.pdf".format(args["version"]))
    plotgraphLineRatio(groups,"ExpVsObs.{}.pdf".format(args["version"]),AreaCutoff=0.6)
    plotgraphLineRatio(groups,"ExpVsObs.Labeled.{}.pdf".format(args["version"]),label=True,AreaCutoff=0.6)
    plotgraphViolin(groups,"Violin.selectionRatio.{}.pdf".format(args["version"]))
    plotgraphCombineViolin(groups,"ViolinCombine.selectionRatio.{}.pdf".format(args["version"]))
    plotgraphRatioMutations(groups,"SelectionGraph.ContactMutations.{}.pdf".format(args["version"]),goodMuts)

    #redo chi2 graph using 1-(freq of no mutation at that position)
    siteChi2=[]
    for group in groups.keys():
        for subject in groups[group].keys():
            for mut in groups[group][subject].siteChi2.keys():
                siteChi2.append((abs(groups[group][subject].siteChi2[mut]),mut))
    goodMuts=sort_mutations_chi2(siteChi2,5000)
    plotChi2graph(groups,"Chi2graph.SiteSelection.AllPositions.{}.pdf".format(args["version"]),goodMuts,'siteChi2')
    plotgraphRatioV2(groups,"PositionAntiselectionGraph.{}.pdf".format(args["version"]),goodMuts)
    goodMuts=sort_mutations_chi2(siteChi2,args["top"])
    plotChi2graph(groups,"Chi2graph.SiteSelection.top{}.{}.pdf".format(args["top"],args["version"]),goodMuts,'siteChi2')

    goodMuts=["T30T","Y33Y","W47W","W50W","N52N","N54N","S55S","G56G","G57G","T58T","N59N","Y60Y","A61A","Q62Q","K63K","T69T","R72R","T74T","E106E","W107W"]
    plotChi2graph(groups,"Chi2graph.SiteSelection.contactPositions.{}.pdf".format(args["version"]),goodMuts[::-1],'siteChi2')
    
    write_out_percentages(groups,"MutationCounts.{}.txt".format(args["version"]))
    
if __name__=="__main__":
    args=parse_args()
    main(args)
