#!/usr/bin/python3

#libraries
import xlrd, re, argparse
import glob
import os
import string
import sys
import numpy,scipy, statistics
import Bio
from Bio.Seq import Seq
from Bio import pairwise2
#from Bio.pairwise2 import format_alignment
from Bio import Align
from Bio import AlignIO
from matplotlib import rcParams
rcParams['font.family']='monospace'
import matplotlib.pyplot as plt

MIN_PYTHON=(3,5)
if sys.version_info < MIN_PYTHON:
    sys.exit("Python %s.%s or later is required.\n" % MIN_PYTHON)
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

        self.NumbMutations=len(self.mutations)
        return self.mutations
    def getDNAMutations(self):
        self.dnamutations=[(pos+1,mut) for pos,mut in enumerate(self.seq) if mut!=self.uca[pos]] 
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

    def printResults(self):
        counts={"20":0,"10":0,"2":0,"1":0,"0.1":0,"0.01":0,"<0.01":0}
        for mut in self.mutations:
            if mut[1]=="X" or mut[1]=="*":
                continue
            if str(mut[0]) in self.freqFileContents.keys():
                mutProbValue=self.freqFileContents["%d"%mut[0]][self.AApos[mut[1]]]
            else:
                return ""
                

            if mutProbValue>=0.20:
                counts["20"]+=1
            elif mutProbValue>=0.10:
                counts["10"]+=1
            elif mutProbValue>=0.02:
                counts["2"]+=1
            elif mutProbValue>=0.01:
                counts["1"]+=1
            elif mutProbValue>=0.001:
                counts["0.1"]+=1
            elif mutProbValue>=0.0001:
                counts["0.01"]+=1
            elif mutProbValue<0.0001:
                counts["<0.01"]+=1
        stringout=[self.name]
        for key in ["20","10","2","1","0.1","0.01","<0.01"]:
            stringout.append(str(counts[key]))
        return "\t".join(stringout)
            #add in counting here for distinguishing between types of probability mutations
    
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

def writeSMUA(outfilename,name,seq):
    with open(outfilename,'a+') as outf:
        outf.write(name+"\n")
        outf.write(seq+"\n")

def printGraph(dataset,pdfname):
    fig,ax=plt.subplots(figsize=(6,4))
    #ax.grid(visible=True,axis='x',which='both',color='k',linestyle='-',linewidth=2)
    ax.hist(dataset,bins=range(0,numpy.max(dataset)+1,1))
    ax.axvline(numpy.median(dataset),color='red',ls='-',label="median:{}".format(numpy.median(dataset)))
    ax.axvline(numpy.mean(dataset),color='blue',ls='--',label="Mean:{:0.2f}".format(numpy.mean(dataset)))
    ax.axvline(numpy.mean(dataset)+numpy.std(dataset),color='blue',ls=':')
    ax.axvline(numpy.mean(dataset)-numpy.std(dataset),color='blue',ls=':')
    try:
        ax.axvline(statistics.mode(dataset),color='black',ls='-',label="mode:{}".format(statistics.mode(dataset)))
    except:
        ax.axvline(statistics.mode(dataset),color='black',ls='-',label="mode:{}".format(statistics.mode(dataset)[0]))
    plt.xlabel("Number of Nucleotide Mutations")
    plt.ylabel("Count of Sequences")
    plt.title(pdfname.replace(".pdf",""))
    plt.legend()
    plt.savefig(pdfname,bbox_inches='tight')
    plt.close()
        
def main(UCAfile,SMUfile,DirFiles):
    with open(SMUfile+".error",'w') as SMUerror:
        pass
    with open(UCAfile,'r') as ucaFile:
        ucaname=ucaFile.readline()
        ucaSeq=ucaFile.readline().rstrip()
    markupseq=ucaSeq
    columnheaders=["name",">20","20-10","10-2","2-1","1-0.1","0.1-0.01","<0.01"]
    print("\t".join(columnheaders))
    mut_counts=[]
    with open(SMUfile,'r') as SMU:
        while True:
            #read batches of lines
            name=SMU.readline().rstrip()
            seq=SMU.readline().rstrip()


            if not name:#break if no more lines
                break
            #process data
            sequence=SeqCompare(name,seq,ucaSeq,markupseq)
            if sequence.checkandFix():
                muts=sequence.getMutations()
                dna_muts=sequence.getDNAMutations()
                mut_counts.append(dna_muts)
                sequence.getFile(dna_muts,DirFiles)
                print(sequence.printResults())
            else:
                writeSMUA(SMUfile+".error",name,seq)
    with open(SMUfile.replace(".fasta",".counts.txt"),"w") as cfiles:
        cfiles.write("{}".format(SMUfile))
        cfiles.write("\tmedian\t{}".format(statistics.median(mut_counts)))
        cfiles.write("\tmean\t{:0.2f}".format(statistics.mean(mut_counts)))
        try:
            cfiles.write("\tmode\t{}\n".format(statistics.mode(mut_counts)))
        except:
            cfiles.write(statistics.mode(mut_counts))
        
    printGraph(mut_counts,SMUfile.replace('.fasta','.counts.pdf'))
                
if __name__=="__main__":
    if len(sys.argv)<2:
        print("usage:python get_NumbMuts.py UCAfile VDJ DirWithtables")
    else:
        main(sys.argv[1],sys.argv[2],sys.argv[3])
