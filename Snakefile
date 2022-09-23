import glob
import os,sys
import shutil

#Paths
PATH="/datacommons/dhvi/scripts/snakemake_pipelines/LookupTable_analysis"

#variables
UCAFILE=config["ucaFile"]
UCA_MARKUP=[]
if config["ucaMarkup"]:
    UCA_MARKUP=config["ucaMarkup"]
tmp=os.path.split(UCAFILE)[-1]
TableName=tmp.replace(".fasta","")

#getting file names
data_dir=[]
sample_dir2=[]
vdjnames=[]

if config["sample"]:
   sample=config["sample"]
   pieces=sample.split("/")
   vdjnames.append(pieces[-1].replace(".fasta",""))
   sample_dir2.append(pieces[-2])
   data_dir="/".join(pieces[:-2])
elif config["basedir"]:
   data_dir=config["basedir"]
   sample_dir2,vdjnamestmp = glob_wildcards(data_dir+"/{dir}/{file}.human_VDJ.fxnl.fasta")
   vdjnames=[x+".human_VDJ.fxnl" for x in vdjnamestmp]

with open(UCAFILE,'r') as file:
    lines=file.readlines()
UCAcount=len(lines[-1].rstrip())
UCATitle=lines[0][1:].rstrip()
paramspace = range(0,UCAcount)

#maxjobs=4

if config["ucaMarkup"]:
    rule LookupAnalysis:#the top level rule to make sure the necessary files/dirs are produced
        input:
            expand(data_dir+"/{sample_dir}/{vdjname}.aa.aligned.freq_table.txt", zip,sample_dir=sample_dir2,vdjname=vdjnames),
            expand(data_dir+"/{sample_dir}/{vdjname}.AANTmutations.V1.pdf",zip,sample_dir=sample_dir2,vdjname=vdjnames),
            expand(data_dir+"/{sample_dir}/{vdjname}.ImprobCounts.txt",zip,sample_dir=sample_dir2,vdjname=vdjnames),
            expand(data_dir+"/{sample_dir}/{vdjname}.ratios.txt",zip,sample_dir=sample_dir2,vdjname=vdjnames)
else:
    rule LookupAnalysis:#the top level rule to make sure the necessary files/dirs are produced
        input:
            expand(data_dir+"/{sample_dir}/{vdjname}.aa.aligned.freq_table.txt", zip,sample_dir=sample_dir2,vdjname=vdjnames),
            expand(data_dir+"/{sample_dir}/{vdjname}.ImprobCounts.txt",zip,sample_dir=sample_dir2,vdjname=vdjnames),
            expand(data_dir+"/{sample_dir}/{vdjname}.ratios.txt",zip,sample_dir=sample_dir2,vdjname=vdjnames)
        
rule GenerateLookupTable:
    input:
        expand(UCATitle+".N{N}.freq_table.txt",N=paramspace),
        d=PATH+"/"+TableName
    output:
        tlog=PATH+"/"+TableName+"_TableGeneration.log"
    shell:
        "mv {UCATitle}.N*.freq_table.txt {input.d}/.;"
	"rm -f {UCATitle}.ARMADiLLO.Detailed.text;"
	"rm -f {UCATitle}.ARMADiLLO.fasta;"
        "touch {output.tlog};"

rule RunArmadillo:
     input:
        uca={UCAFILE},
        d=PATH+"/"+TableName,
        armadillo_path=PATH+"/programs/armadillo",
        armadillo=PATH+"/programs/armadillo/ARMADiLLO"
     output:
        UCATitle+".N{N}.freq_table.txt"
     run:
        max_iter=1000000
        seed=202108024
        shell("{input.armadillo} -m {input.armadillo_path}/Mutability.csv -s {input.armadillo_path}/Substitution.csv -uca {input.uca} -seq {input.uca} -text -n {wildcards.N} -max_iter {max_iter} -threads 1 -random_seed {seed}")
        #shell("mv {UCATitle}.N{wildcards.N}.freq_table.txt {input.d}/.")

rule GenMutGraphs:
     input:
        program=PATH+"/programs/plot_seq_mutations.py",
        vdj=data_dir+"/{sample_dir}/{vdjname}.fasta",
        uca={UCAFILE},
        markup=config["ucaMarkup"]
     output:
        mutationpdf=data_dir+"/{sample_dir}/{vdjname}.AANTmutations.V1.pdf"
     shell:
        "python3 {input.program} -u {input.uca} -s {input.vdj} -m {input.markup} -v V1"
    
rule GenerateDir:
     input:
     output:
        directory(PATH+"/"+TableName)
     shell:
        "mkdir {output}"

rule compileTranslater:
     input:
         PATH+"/programs/translater"
     output:
        PATH+"/programs/translater/translater"
     shell:
        "cd {input}; make"
rule compileAlign:
     input:
        PATH+"/programs/align_all_to_first"
     output:
        PATH+"/programs/align_all_to_first/align_all_to_first"
     shell:
        "cd {input}; make"
rule compileMatrix:
     input:
        PATH+"/programs/positional_aa_matrix_maker"
     output:
        PATH+"/programs/positional_aa_matrix_maker/positional_aa_matrix_maker"
     shell:
        "cd {input}; make"
rule compileArmadillo:
     input:
        PATH+"/programs/armadillo"
     output:
        armadillo=PATH+"/programs/armadillo/ARMADiLLO"
     shell:
        "cd {input}; make"

rule runVDJCalc:
     input:
        uca={UCAFILE},
        vdj=data_dir+"/{sample_dir}/{vdjname}.fasta",
        d=PATH+"/"+TableName,
	glog=PATH+"/"+TableName+"_TableGeneration.log"
     output:
        counts=temp(data_dir+"/{sample_dir}/{vdjname}.counts.txt"),
	mutcounts=data_dir+"/{sample_dir}/{vdjname}.ImprobCounts.txt"
     shell:
        "python3 {PATH}/programs/get_NumbMuts_VDJ.py {input.uca} {input.vdj} {input.d} > {output.mutcounts}"

rule PredictedAAfreqTable:
     input:
        glog=PATH+"/"+TableName+"_TableGeneration.log",
        aafile=data_dir+"/{sample_dir}/{vdjname}.aa.aligned.freq_table.txt",
        countfile=data_dir+"/{sample_dir}/{vdjname}.counts.txt",
        uca={UCAFILE},
        d=PATH+"/"+TableName
     output:
        writelist=data_dir+"/{sample_dir}/{vdjname}.ratios.txt",
	junkfile=temp(data_dir+"/{sample_dir}/{vdjname}.junk.txt"),
        olog=temp(data_dir+"/{sample_dir}/{vdjname}.freq.log")
     run:
        counts=[]
        with open(input.countfile,"r") as f:
             counts=f.readline().rstrip().split()[2]
        Muttable=glob.glob(os.path.join(input.d,UCATitle+".N{}.freq_table.txt".format(int(round(float(counts))))))
        shell("python3 {PATH}/programs/process_Freq_Table.py -u {input.uca} -l {input.aafile} -f {Muttable} --writelist {output.writelist} -o {output.junkfile} > {output.olog}")

rule translate:
     input:
        file=data_dir+"/{sample_dir}/{vdjname}.fasta",
        translater=PATH+"/programs/translater/translater"	
     output:
        temp(data_dir+"/{sample_dir}/{vdjname}.aa.tmp")
     shell:
        "{input.translater} -i {input.file} > {output}"

rule alignAAFasta:
     input:
        program=PATH+"/programs/align_all_to_first/align_all_to_first",
        aafile=data_dir+"/{sample_dir}/{vdjname}.aa.tmp"
     output:
        alignedFasta=temp(data_dir+"/{sample_dir}/{vdjname}.aa.aligned"),
        junkfile=temp(data_dir+"/{sample_dir}/{vdjname}.junk.fasta")
     shell:
        "{input.program} -i {input.aafile} -matrix {PATH}/programs/align_all_to_first/BLOSUM62.txt -retain_ins_to_file {output.junkfile} -go 11 -ge 1 > {output.alignedFasta}"

rule ObservedAAfreqTable:
     input:
        program=PATH+"/programs/positional_aa_matrix_maker/positional_aa_matrix_maker",
	alignedFasta=data_dir+"/{sample_dir}/{vdjname}.aa.aligned"
     output:
        aatable=data_dir+"/{sample_dir}/{vdjname}.aa.aligned.freq_table.txt"
     shell:
        "{input.program} -i {input.alignedFasta} -o {output.aatable}"

rule cleanTranslater:
     input:
        d=PATH+"/programs/translater",
        p=PATH+"/programs/translater/translater"
     output:
        temp(PATH+"/translater.clean")
     shell:
        "cd {input.d}; make clean;"
	"touch {output}"
rule cleanAlign:
     input:
        d=PATH+"/programs/align_all_to_first",
        p=PATH+"/programs/align_all_to_first/align_all_to_first",

     output:
        temp(PATH+"/align.clean")
     shell:
        "cd {input.d}; make clean;"
	"touch {output}"	
rule cleanMatrix:
     input:
        d=PATH+"/programs/positional_aa_matrix_maker",
        p=PATH+"/programs/positional_aa_matrix_maker/positional_aa_matrix_maker"
     output:
        temp(PATH+"/matrix.clean")
     shell:
        "cd {input.d}; make clean;"
	"touch {output}"
rule cleanArmadillo:
     input:
        d=PATH+"/programs/armadillo",
        armadillo=PATH+"/programs/armadillo/ARMADiLLO"
     output:
        temp(PATH+"/armadillo.clean")
     shell:
        "cd {input.d}; make clean;"
	"touch {output}"

rule clean: #resets the whole proccess by removing the outputs
    message:
        "cleaning the generated files"
    run:
        for s_dir in set(sample_dir2):
             for file in glob.iglob(os.path.join(data_dir,s_dir,"*.AANTmutations.V1.pdf")):
                  os.remove(file)
             for file in glob.iglob(os.path.join(data_dir,s_dir,"*.ImprobCounts.txt")):
                  os.remove(file)
             for file in glob.iglob(os.path.join(data_dir,s_dir,"*.aa.aligned.freq_table.txt")):
                  os.remove(file)
             for file in glob.iglob(os.path.join(data_dir,s_dir,"*.ratios.txt")):
                  os.remove(file)

rule partialclean: #resets the whole proccess by removing the outputs
    message:
        "partial cleaning of the generated files"
    run:
        for s_dir in set(sample_dir2):
             for file in glob.iglob(os.path.join(data_dir,s_dir,"*.AANTmutations.V1.pdf")):
                  os.remove(file)
             for file in glob.iglob(os.path.join(data_dir,s_dir,"*.ImprobCounts.txt")):
                  os.remove(file)
             for file in glob.iglob(os.path.join(data_dir,s_dir,"*.ratios.txt")):
                  os.remove(file)
	
rule fullclean: #resets the whole proccess by removing the outputs
    message:
        "complete cleaning of compiled programs, lookup table and generated files"
    input:
        PATH+"/armadillo.clean",
        PATH+"/matrix.clean",
        PATH+"/align.clean",
        PATH+"/translater.clean"
    run:
        shell("rm -rf {PATH}/{TableName} *.log")
        for s_dir in set(sample_dir2):
             for file in glob.iglob(os.path.join(data_dir,s_dir,"*.AANTmutations.V1.pdf")):
                  os.remove(file)
             for file in glob.iglob(os.path.join(data_dir,s_dir,"*.ImprobCounts.txt")):
                  os.remove(file)
             for file in glob.iglob(os.path.join(data_dir,s_dir,"*.aa.aligned.freq_table.txt")):
                  os.remove(file)
             for file in glob.iglob(os.path.join(data_dir,s_dir,"*.ratios.txt")):
                  os.remove(file)
             for file in glob.iglob(os.path.join(data_dir,s_dir,"*.counts.pdf")):
                  os.remove(file)
