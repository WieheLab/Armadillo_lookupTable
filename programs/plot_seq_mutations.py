#!/usr/bin/python3

#libraries
import xlrd, re
import glob
import os
import string
import sys
import argparse
import numpy,scipy, statistics
import Bio
from Bio.Seq import Seq
from Bio import pairwise2
#from Bio.pairwise2 import format_alignment
from Bio import Align
from Bio import AlignIO
from matplotlib.patches import Rectangle
from matplotlib import rcParams
rcParams['font.family']='monospace'
import matplotlib.pyplot as plt

#file=sys.stderr

class SeqCompare():
    def __init__(self,name,seq,ucaseq,markup):
        self.name=name.rstrip()
        self.seq=seq
        self.uca=ucaseq
        self.markup=markup

    def checkandFix_complex(self):
        #print(self.seq)
        #print(len(self.seq))
        #print(self.uca)
        #print(len(self.uca))
        aligner=Align.PairwiseAligner()
        aligner.match_score=2
        aligner.mismatch_score=-0.5
        aligner.target_left_open_gap_score=-8
        aligner.target_right_open_gap_score=-8
        aligner.query_internal_open_gap_score=-100
        aligner.query_internal_extend_gap_score=-100
        aligner.query_left_open_gap_score=-100
        aligner.query_left_extend_gap_score=-8
        aligner.query_right_open_gap_score=-100
        aligner.query_right_extend_gap_score=-8
        
        aligner.mode="global"
        #print(aligner)
        #input('e')
        if len(self.seq)<len(self.uca):
            sequca=Seq(self.uca)
            seqseq=Seq(self.seq)
            alignments=aligner.align(seqseq,sequca)[0]
            #alignments=pairwise2.align.globalms(seqseq,sequca,2,-1,-10,-0.5)
            #print(alignments)
            #input('e')
            seqseq=alignments.__str__().split("\n")[0]
            #print(seqseq)
            #input('e')
            self.seq=seqseq
            if "-" in alignments.__str__().split("\n")[2]:
                return False
        elif len(self.seq)>len(self.uca):
            return False
        return True
    def checkandFix(self):
        if len(self.seq)<len(self.uca):
            sequca=Seq(self.uca)
            seqseq=Seq(self.seq)
            #alignments=aligner.align(seqseq,sequca)[0]
            alignments=pairwise2.align.globalms(seqseq,sequca,5,-4,-10,-0.5)[0]
            #print(alignments)
            #input('e')
            try:
                seqseq=alignments.seqA
            except:
                return False
            #print(seqseq)
            #input('e')
            self.seq=seqseq
            if "-" in alignments.seqB:
                return False
        elif len(self.seq)>len(self.uca):
            return False
        return True
    
    def translate(self,ntseq):
        dna_to_aa=DNAtoAAmap()
        tranSeq=[]
        for i in range(0, int(len(ntseq)/3)*3, 3):
            try:
                tranSeq.append(dna_to_aa[ntseq[i:i+3]])
            except:
                tranSeq.append("X")
        return tranSeq
    def getMutations(self):
        self.AAseq=self.translate(self.seq)
        self.AAuca=self.translate(self.uca)

        self.mutations=[(pos+1,mut) for pos,mut in enumerate(self.AAseq) if mut!=self.AAuca[pos]]

        self.AAdiff=[]
        for pos,AA in enumerate(self.AAseq):
            if AA==self.AAuca[pos]:
                self.AAdiff.append(".")
            else:
                self.AAdiff.append(AA)

        self.NumbMutations=len(self.mutations)
        return self.mutations
    def getDNAMutations(self):
        self.dnamutations=[(pos+1,mut) for pos,mut in enumerate(self.seq) if mut!=self.uca[pos]]

        self.NTdiff=[]
        for pos,nt in enumerate(self.seq):
            if nt==self.uca[pos]:
                self.NTdiff.append(".")
            else:
                self.NTdiff.append(nt)
        
        return len(self.dnamutations)
    def getFile(self,NumbMutations,freqTableDir):
        freqFiles=glob.glob(os.path.join(freqTableDir,'*%d.freq_table.txt'%NumbMutations))[0]
        freqFileContents=dict()
        with open(freqFiles,'r') as freqF:
            for line in freqF:
                dataLine=line.rstrip().split(',')
                try:
                    freqFileContents[dataLine[0]]=[float(x) for x in dataLine[1:]]
                except:
                    freqFileContents[dataLine[0]]=dataLine[1:]
        self.freqFileContents=freqFileContents
        self.AApos=dict()
        for i,aa in enumerate(freqFileContents['pos']):
            self.AApos[aa]=i

def plot_sequence_differences(seqs,markup,figname,label):
    c=[1,0,0]
    fig,ax=plt.subplots(2,1,figsize=(12,4))
    
    ax[0].set_xlim([-1,len(markup)])
    ax[0].set_ylim([0,1])
    ax[1].set_xlim([-1,len(markup)])
    ax[1].set_ylim([0,1])
    yrange=ax[0].get_ylim()
    #print(datanames)
    
    region_pos={"FR1":(-1,-1),"FR2":(-1,-1),"FR3":(-1,-1),"FR4":(-1,-1),"CDR1":(-1,-1),"CDR2":(-1,-1),"CDR3":(-1,-1)}
    keyhash={"FR1":"1","FR2":"2","FR3":"3","FR4":"4","CDR1":"A","CDR2":"B"}
    for n,key in enumerate(region_pos.keys()):
        if "CDR3"==key:
            continue
        r=[i for i,x in enumerate(markup) if x==keyhash[key]]
        region_pos[key]=(min(r),max(r))
        if "FR" in key:
            c=[1,1,1]
        else:
            c=[0.5,0.5,1]
        ax[0].add_patch(Rectangle((min(r), 0), max(r)-min(r)+0.9, yrange[1],facecolor=c,alpha=0.5,zorder=0))
        ax[0].text((min(r)+max(r))/2.0,1.02,key,ha='center')
        ax[1].add_patch(Rectangle((min(r), 0), max(r)-min(r)+0.9, yrange[1],facecolor=c,alpha=0.5,zorder=0))
        #ax[1].text(min(r)+5,1.02,key)
        #ax.axvline(min(r),color='black',zorder=1,linestyle=':')

    r=[i for i,x in enumerate(markup) if x not in list(keyhash.values()) and x!= "U"]
    ax[0].add_patch(Rectangle((min(r), 0), max(r)-min(r)+0.9, yrange[1],facecolor=[0.5,0.5,1],alpha=0.5,zorder=0))
    ax[0].text((min(r)+max(r))/2.0,1.02,"CDR3",ha="center")
    ax[1].add_patch(Rectangle((min(r), 0), max(r)-min(r)+0.9, yrange[1],facecolor=[0.5,0.5,1],alpha=0.5,zorder=0))
    #ax[1].text(min(r)+5,1.02,"CDR3")

    ax[0].yaxis.grid(True)
    ax[1].yaxis.grid(True)
    
    ntdiff=[{"A":0,"G":0,"C":0,"T":0,"-":0,"N":0} for x in seqs[0].NTdiff]
    for seq in seqs:
        for i,nt in enumerate(seq.NTdiff):
            if nt==".":
                continue
            else:
                ntdiff[i][nt]+=1

    ax[0].bar([i for i,x in enumerate(ntdiff)],[float(sum(pos.values()))/float(len(seqs)) for pos in ntdiff],1,facecolor=[0,0,0])

    H=[]
    nt=list( ntdiff[0].keys())
    for d in nt:
        H.append([x[d]/float(len(seqs)) for x in ntdiff])

    b=[0 for x in ntdiff]
    for i,h in enumerate(H):
        ax[1].bar([i for i,x in enumerate(ntdiff)],h,1,bottom=b,label="{}".format(nt[i].replace("-","del")))
        for n,x in enumerate(h):
            b[n]=b[n]+x

    ax[0].set_title("{} N={}".format(label,len(seqs)),y=1.1)
    ax[0].set_ylabel("Percent Mutations")
    ax[0].set_xticklabels([])

    ax[1].legend(loc='upper center', bbox_to_anchor=(0.5, 1.20), ncol=6,shadow=True)
    ax[1].set_ylabel("Percent Mutations")
    ax[1].set_xlabel("Positions")
    plt.savefig(figname,bbox_inches='tight')
    plt.close()

def plot_AAsequence_differences(seqs,markup,figname,label):
    c=[1,0,0]
    fig,ax=plt.subplots(figsize=(12,4))
    
    ax.set_xlim([-1,len(markup)/3.0])
    ax.set_ylim([0,1])
    yrange=ax.get_ylim()
    #print(datanames)
    
    region_pos={"FR1":(-1,-1),"FR2":(-1,-1),"FR3":(-1,-1),"FR4":(-1,-1),"CDR1":(-1,-1),"CDR2":(-1,-1),"CDR3":(-1,-1)}
    keyhash={"FR1":"1","FR2":"2","FR3":"3","FR4":"4","CDR1":"A","CDR2":"B"}
    for n,key in enumerate(region_pos.keys()):
        if "CDR3"==key:
            continue
        r=[i for i,x in enumerate(markup) if x==keyhash[key]]
        region_pos[key]=(min(r),max(r))
        if "FR" in key:
            c=[1,1,1]
        else:
            c=[0.5,0.5,1]
        ax.add_patch(Rectangle((min(r)/3.0, 0), max(r)/3.0-min(r)/3.0+0.9, yrange[1],facecolor=c,alpha=0.5,zorder=0))
        ax.text((min(r)+max(r))/6.0,1.02,key,ha='center')

        #ax[1].text(min(r)+5,1.02,key)
        #ax.axvline(min(r),color='black',zorder=1,linestyle=':')

    r=[i for i,x in enumerate(markup) if x not in list(keyhash.values()) and x!= "U"]
    ax.add_patch(Rectangle((min(r)/3.0, 0), max(r)/3.0-min(r)/3.0+0.9, yrange[1],facecolor=[0.5,0.5,1],alpha=0.5,zorder=0))
    ax.text((min(r)+max(r))/6.0,1.02,"CDR3",ha="center")

    ax.yaxis.grid(True)

    AAdiff=[0 for x in seqs[0].AAdiff]
    for seq in seqs:
        for i,nt in enumerate(seq.AAdiff):
            if nt==".":
                continue
            else:
                AAdiff[i]+=1

    ax.bar([i for i,x in enumerate(AAdiff)],[float(pos)/float(len(seqs)) for pos in AAdiff],1,facecolor=[0,0,0])

    ax.set_title("{} N={}".format(label,len(seqs)),y=1.1)
    ax.set_ylabel("Percent Mutations")
    ax.set_xlabel("AA Positions")
    plt.savefig(figname,bbox_inches='tight')
    plt.close()

def plot_AAnt_differences(seqs,markup,figname,label):
    c=[1,0,0]
    fig,ax=plt.subplots(figsize=(12,3))
    
    ax.set_xlim([-1,len(markup)])
    ax.set_ylim([0,1])
    yrange=ax.get_ylim()
    #print(datanames)
    
    region_pos={"FR1":(-1,-1),"FR2":(-1,-1),"FR3":(-1,-1),"FR4":(-1,-1),"CDR1":(-1,-1),"CDR2":(-1,-1),"CDR3":(-1,-1)}
    keyhash={"FR1":"1","FR2":"2","FR3":"3","FR4":"4","CDR1":"A","CDR2":"B"}
    for n,key in enumerate(region_pos.keys()):
        if "CDR3"==key:
            continue
        r=[i for i,x in enumerate(markup) if x==keyhash[key]]
        region_pos[key]=(min(r),max(r))
        if "FR" in key:
            c=[1,1,1]
        else:
            c=[0.5,0.5,1]
        ax.add_patch(Rectangle((min(r), 0), max(r)-min(r)+0.9, yrange[1],facecolor=c,alpha=0.5,zorder=0))
        ax.text((min(r)+max(r))/2.0,1.02,key,ha='center')

        #ax[1].text(min(r)+5,1.02,key)
        #ax.axvline(min(r),color='black',zorder=1,linestyle=':')

    r=[i for i,x in enumerate(markup) if x not in list(keyhash.values()) and x!= "U"]
    ax.add_patch(Rectangle((min(r), 0), max(r)-min(r)+0.9, yrange[1],facecolor=[0.5,0.5,1],alpha=0.5,zorder=0))
    ax.text((min(r)+max(r))/2.0,1.02,"CDR3",ha="center")

    ax.yaxis.grid(True)

    AAdiff=[0 for x in seqs[0].AAdiff]
    for seq in seqs:
        for i,nt in enumerate(seq.AAdiff):
            if nt==".":
                continue
            else:
                AAdiff[i]+=1

    ax.bar([3*i+1 for i,x in enumerate(AAdiff)],[float(pos)/float(len(seqs)) for pos in AAdiff],3,facecolor=[0,0,0])

    ntdiff=[{"A":0,"G":0,"C":0,"T":0,"-":0,"N":0} for x in seqs[0].NTdiff]
    for seq in seqs:
        for i,nt in enumerate(seq.NTdiff):
            if nt==".":
                continue
            else:
                ntdiff[i][nt]+=1
    H=[]
    nt=list( ntdiff[0].keys())
    for d in nt:
        H.append([x[d]/float(len(seqs)) for x in ntdiff])

        
    b=[0 for x in ntdiff]
    for i,h in enumerate(H):
        if nt[i]=="N":
            ax.bar([i for i,x in enumerate(ntdiff)],h,1,bottom=b,alpha=0.75)
        else:
            ax.bar([i for i,x in enumerate(ntdiff)],h,1,bottom=b,label="{}".format(nt[i].replace("-","del")),alpha=0.75)
        for n,x in enumerate(h):
            b[n]=b[n]+x
    ax.legend(loc='upper right', bbox_to_anchor=(1, -0.1), ncol=6)
    ax.set_title("{} N={}".format(label,len(seqs)),y=1.1)
    ax.set_ylabel("Mutation Frequency")
    ax.set_xlabel("Positions")
    plt.savefig(figname,bbox_inches='tight')
    plt.close()
    
def DNAtoAAmap():
    dna_to_aa_tranx_map=dict()
    #ALA
    dna_to_aa_tranx_map["GCT"]="A";  
    dna_to_aa_tranx_map["GCC"]="A";  
    dna_to_aa_tranx_map["GCA"]="A";  
    dna_to_aa_tranx_map["GCG"]="A";  
    #ARG
    dna_to_aa_tranx_map["CGT"]="R";  
    dna_to_aa_tranx_map["CGC"]="R";  
    dna_to_aa_tranx_map["CGA"]="R";  
    dna_to_aa_tranx_map["CGG"]="R"; 
    dna_to_aa_tranx_map["AGA"]="R"; 
    dna_to_aa_tranx_map["AGG"]="R"; 
    #ASN
    dna_to_aa_tranx_map["AAT"]="N"; 
    dna_to_aa_tranx_map["AAC"]="N"; 
    #ASP
    dna_to_aa_tranx_map["GAT"]="D"; 
    dna_to_aa_tranx_map["GAC"]="D"; 
    #CYS
    dna_to_aa_tranx_map["TGT"]="C"; 
    dna_to_aa_tranx_map["TGC"]="C"; 
    #GLN
    dna_to_aa_tranx_map["CAA"]="Q"; 
    dna_to_aa_tranx_map["CAG"]="Q"; 
    #GLU
    dna_to_aa_tranx_map["GAA"]="E"; 
    dna_to_aa_tranx_map["GAG"]="E"; 
    #GLY
    dna_to_aa_tranx_map["GGT"]="G"; 
    dna_to_aa_tranx_map["GGC"]="G"; 
    dna_to_aa_tranx_map["GGA"]="G"; 
    dna_to_aa_tranx_map["GGG"]="G"; 
    #HIS
    dna_to_aa_tranx_map["CAT"]="H"; 
    dna_to_aa_tranx_map["CAC"]="H"; 
    #ILE
    dna_to_aa_tranx_map["ATT"]="I"; 
    dna_to_aa_tranx_map["ATC"]="I"; 
    dna_to_aa_tranx_map["ATA"]="I"; 
    #LEU
    dna_to_aa_tranx_map["TTA"]="L"; 
    dna_to_aa_tranx_map["TTG"]="L"; 
    dna_to_aa_tranx_map["CTT"]="L"; 
    dna_to_aa_tranx_map["CTC"]="L"; 
    dna_to_aa_tranx_map["CTA"]="L"; 
    dna_to_aa_tranx_map["CTG"]="L"; 
    #LYS
    dna_to_aa_tranx_map["AAA"]="K"; 
    dna_to_aa_tranx_map["AAG"]="K"; 
    #MET
    dna_to_aa_tranx_map["ATG"]="M"; 
    #PHE
    dna_to_aa_tranx_map["TTT"]="F"; 
    dna_to_aa_tranx_map["TTC"]="F"; 
    #PRO
    dna_to_aa_tranx_map["CCT"]="P"; 
    dna_to_aa_tranx_map["CCC"]="P"; 
    dna_to_aa_tranx_map["CCA"]="P"; 
    dna_to_aa_tranx_map["CCG"]="P"; 
    #SER
    dna_to_aa_tranx_map["TCT"]="S"; 
    dna_to_aa_tranx_map["TCC"]="S"; 
    dna_to_aa_tranx_map["TCA"]="S"; 
    dna_to_aa_tranx_map["TCG"]="S"; 
    dna_to_aa_tranx_map["AGT"]="S"; 
    dna_to_aa_tranx_map["AGC"]="S"; 
    #THR
    dna_to_aa_tranx_map["ACT"]="T"; 
    dna_to_aa_tranx_map["ACC"]="T"; 
    dna_to_aa_tranx_map["ACA"]="T"; 
    dna_to_aa_tranx_map["ACG"]="T"; 
    #TRP
    dna_to_aa_tranx_map["TGG"]="W"; 
    #TYR
    dna_to_aa_tranx_map["TAT"]="Y"; 
    dna_to_aa_tranx_map["TAC"]="Y"; 
    #VAL
    dna_to_aa_tranx_map["GTT"]="V"; 
    dna_to_aa_tranx_map["GTC"]="V"; 
    dna_to_aa_tranx_map["GTA"]="V"; 
    dna_to_aa_tranx_map["GTG"]="V"; 
    #STOP
    dna_to_aa_tranx_map["TAA"]="*"; 
    dna_to_aa_tranx_map["TGA"]="*"; 
    dna_to_aa_tranx_map["TAG"]="*"; 
    
    return dna_to_aa_tranx_map

def parse_args():
    # Construct an argument parser
    all_args = argparse.ArgumentParser()
    
    # Add arguments to the parser
    all_args.add_argument("-u","--UCAfile",default=[],required=False,help="UCA file")
    all_args.add_argument("-s","--seqfiles",action="extend",nargs="+",required=True,help="direcotry containing tables")    #DirwithTables    
    all_args.add_argument("-m","--markupfile",default=[],required=False,help="markup file")
    all_args.add_argument("-v","--version",default="V1",required=False,help="version")
  
    args = vars(all_args.parse_args())
    return args
    

    
def main(args):
    markupfile=args["markupfile"]
    UCAfile=args["UCAfile"]
    
    with open(UCAfile,'r') as ucaFile:
        ucaname=ucaFile.readline()
        ucaSeq=ucaFile.readline().rstrip()

    with open(markupfile,'r') as mfile:
        mfilename=mfile.readline()
        markupseq=mfile.readline().strip()


    for fastafile in args["seqfiles"]:
        seqs=[]
        with open(fastafile,'r') as seqfile:
            while True:
                #read batches of lines
                name=seqfile.readline().rstrip()
                seq=seqfile.readline().rstrip()
                if not name:#break if no more lines
                    break
                #process data
                sequence=SeqCompare(name,seq,ucaSeq,markupseq)
                if sequence.checkandFix():
                    muts=sequence.getMutations()
                    dna_muts=sequence.getDNAMutations()
                    seqs.append(sequence)

                else:
                    pass
        setname=fastafile.split("/")
        setname[-1]=setname[-1].replace(".fasta","")
        
        #plot_sequence_differences(seqs,markupseq,"/".join(setname)+".mutations.{}.pdf".format(args["version"]),"/".join(setname))
        #plot_AAsequence_differences(seqs,markupseq,"/".join(setname)+".AAmutations.{}.pdf".format(args["version"]),"/".join(setname))
        plot_AAnt_differences(seqs,markupseq,"/".join(setname)+".AANTmutations.{}.pdf".format(args["version"]),"/".join(setname))


if __name__=="__main__":
    args=parse_args()
    main(args)

