#!/usr/bin/python
#
#
#

import os
import argparse
import sys
import re
import string
import gzip
import math
import numpy

MIN_PYTHON=(3,5)
if sys.version_info < MIN_PYTHON:
    sys.exit("Python %s.%s or later is required.\n" %MIN_PYTHON)

class FreqTable():
    def __init__(self,name):
        self.name=name

    def readFreqTable(self,freqTablefilename):
        freqFileContents=dict()
        with open(freqTablefilename,"r") as freqF:
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
    def readuca(self,UCAfile):
        name=""
        seq=[]
        with open(UCAfile,'r') as ucafile:
            for line in ucafile:
                dataline=line.rstrip()
                if ">" in dataline and len(name)<1:
                    name=dataline[1:]
                elif ">" in dataline:
                    break
                else:
                    seq.append(dataline)
        self.ucaseq="".join(seq)
        self.ucaname=name
            
    def convertTable(self,outname):
        self.writeMatrixTable(self.freqFileContents,self.AApos,outname)

    def writeMatrixTable(self,freqFile,AApos,outname):
        if "pos" in freqFile:
            freqFile.pop("pos")
        positions=[int(x) for x in freqFile.keys()]
        AAlist=["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","X","Y","-"]
        with open(outname,"w") as outTable:
            outTable.write("\t{}\n".format("\t".join(AAlist)))
            for pos in positions:
                outTable.write("{}".format(pos))
                for AA in AAlist:
                    try:
                        p=AApos[AA]
                        outTable.write("\t{:1.6f}".format(freqFile[str(pos)][p]))
                    except:
                        outTable.write("\t0.000001")
                
                outTable.write("\n")
                
    def readLogoTable(self,LogoTable):
        logoFileContents=dict()
        with open(LogoTable,"r") as logoF:
            for line in logoF:
                dataLine=line.rstrip().split("\t")
                if not dataLine[0]:
                    dataLine[0]="pos"
                try:
                    logoFileContents[dataLine[0]]=[float(x) for x in dataLine[1:]]
                except:
                    logoFileContents[dataLine[0]]=dataLine[1:]

        logoAApos=dict()
        for i,aa in enumerate(logoFileContents['pos']):
            logoAApos[aa]=i
        self.logoFileContents=logoFileContents
        self.logoAApos=logoAApos
        return logoAApos,logoFileContents
    def compare(self,l1,l2):
        # here l1 and l2 must be lists
        #print(l1)
        #print(l2)
        #input('e')
        #if len(l1) != len(l2):
        #    return False
        l1.sort()
        l2.sort()
        if l1 == l2 or l1[-1]+1==l2[-1] or l2[-1]+1==l1[-1]:
            return True
        else:
            return False
    def compareTables(self,logoAApos,logoFileContents,outnameComparison):
        if "pos" in self.freqFileContents:
            self.freqFileContents.pop("pos")
        if "pos" in logoFileContents:
            logoFileContents.pop("pos")
        positions=[int(x) for x in self.freqFileContents.keys()]
        logoPositions=[int(x) for x in logoFileContents.keys()]
        if not self.compare(positions,logoPositions):
            print("Tables are not compatible. Make sure they are the same length.")
            print(positions)
            print(logoPositions)
            return
        AAlist=["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","X","Y","-"]
        AAdict=dict()
        for i,AA in enumerate(AAlist):
            AAdict[AA]=i
        self.LogRatio=dict()
        for pos in positions:
            ratios=[]
            for AA in AAlist:
                freqTableValue=0
                logoTableValue=0
                try:
                    freqTableValue=self.freqFileContents[str(pos)][self.AApos[AA]]
                except:
                    freqTableValue=0.0
                #if freqTableValue<self.cutoffValue:
                #    freqTableValue=0.0
                    
                try:
                    logoTableValue=logoFileContents[str(pos)][logoAApos[AA]]
                except:
                    logoTableValue=0.0
                ratio=0
                if logoTableValue==0.0 and freqTableValue==0.0:
                    ratios.append(math.nan)
                elif logoTableValue==0.0:
                    ratios.append(0.0000)
                elif freqTableValue==0.0:
                    ratios.append(math.inf)
                else:
                    print("{}{}".format(AA,pos))
                    print(logoTableValue)
                    print(freqTableValue)
                    print(logoTableValue/float(freqTableValue))
                    ratio=logoTableValue/float(freqTableValue)
                    print(math.log10(ratio))
                    ratios.append(math.log10(ratio))
                    print(ratios)
                    #input('e')
            self.LogRatio[str(pos)]=ratios

        self.AAdict=AAdict
        if outnameComparison:
            self.writeMatrixTable(self.LogRatio,self.AAdict,outnameComparison)
        return self.LogRatio,AAdict

    def rankorder(self,freqField,AApos,outname):
        #Mutlist=dict()
        Mutlist=[]
        for pos in freqField.keys():
            for AA in AApos.keys():
                Mutlist.append(("{}{}".format(AA,pos),"{:.6f}".format(freqField[pos][AApos[AA]])))

        print(Mutlist)
        input('e')
        Mutlist2=sorted(Mutlist,key=sortSecond)
        #                lambda v:float(v[1]))
        input('e')
        for M in Mutlist2:
            print(M)
        #input('e')
        
    def convert(self,ntseq):
        dna_to_aa=DNAtoAAmap()
        if "-" in ntseq:
            nt=ntseq.replace("-","")
            ntseq=nt
        AAseq=[]
        for i in range(0, int(len(ntseq)/3)*3, 3):
            try:
                AAseq.append(dna_to_aa[ntseq[i:i+3]])
            except:
                AAseq.append("X")
        return "".join(AAseq)
    
    def printList(self,freqField,AApos,outname):
        Mutlist=[]
        #G31D ARM_freq obs_Freq foldchange
        header=["mutation"]
        if "freqFileContents" in dir(self):
            header.append("ARMADiLLO Freq")
        if "logoFileContents" in dir(self):
            header.append("Observed Freq")
        header.append("Log10 foldChange")
        if "ucaseq" in dir(self):
            UCAAA_seq=self.convert(self.ucaseq)
        for pos in freqField.keys():
            for AA in AApos.keys():
                datalist=[]
                if "ucaseq" in dir(self):
                    try:
                        ucaAA=UCAAA_seq[int(pos)-1]
                    except:
                        pass
                        #print(int(pos))
                        #print(UCAAA_seq)
                        #print(len(UCAAA_seq))
                        #input('e')
                    datalist.append("{}{}{}".format(ucaAA,pos,AA))
                else:
                    datalist.append("{}{}".format(pos,AA))

                if "freqFileContents" in dir(self):
                    try:
                        datalist.append("{:.6f}".format(self.freqFileContents[pos][self.AApos[AA]]))
                    except:
                        datalist.append("{:.6f}".format(0.000001))
                if "logoFileContents" in dir(self):
                    try:
                        datalist.append("{:.6f}".format(self.logoFileContents[pos][self.logoAApos[AA]]))
                    except:
                        datalist.append("{:.6f}".format(0.000000))
                    
                datalist.append("{:.6f}".format(freqField[pos][AApos[AA]]))
                Mutlist.append(datalist)
                    
        with open(outname,"w") as ofile:
            ofile.write("{}\n".format("\t".join(header)))
            for muts in Mutlist:
                ofile.write("{}\n".format("\t".join(muts)))
        
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

def sortSecond(val):
    #print(val)
    #print(float(val[1]))
    return float(val[1])
     
def parse_args():
    all_args=argparse.ArgumentParser()
    
    all_args.add_argument("-u","--uca",required=False,help="UCA fasta file")
    all_args.add_argument("-f","--infile",required=False,help="name of input file")
    all_args.add_argument("-c","--convert",required=False,help="freq_table to convert",default=[])
    all_args.add_argument("-l","--logofile",required=False,help="logo freq_table to convert",default=[])
    all_args.add_argument("-o","--outfile",required=False,help="name of outfile for comparison",default=[])
    all_args.add_argument("--writelist",required=False,help="write ratios to list")
    #all_args.add_argument("-r","--rankorder",required=False,help="name of rank ordered file of comparison")

    args=vars(all_args.parse_args())
    return args

def main(args):

    RunBool=True
    if os.path.exists(args["infile"]):
        filename=args['infile']
    else:
        print("can not find infile: {}".format(args["infile"]))
        RunBool=False

    if not RunBool:
        sys.exit("Program failed")

    freqtable=FreqTable("test")
    freqtable.readFreqTable(args["infile"])
    if args["uca"]:
        freqtable.readuca(args["uca"])
    if args["convert"]:
        freqtable.convertTable(args["convert"])
    logoAApos,logoContents=freqtable.readLogoTable(args["logofile"])

    CompFreq,CompAA=freqtable.compareTables(logoAApos,logoContents,args["outfile"])
    if args["writelist"]:
        freqtable.printList(CompFreq,CompAA,args["writelist"])
    #if args["rankorder"]:
    #    freqtable.rankorder(CompFreq,CompAA,args["rankorder"])
    return

if __name__=="__main__":
    args=parse_args()
    main(args)
