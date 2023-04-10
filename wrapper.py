#! /usr/bin/env python3

import subprocess
import argparse
import os
import re
import glob
import zipfile

parser = argparse.ArgumentParser()

parser.add_argument("-b",type=str,required=True,help="Base name used in all the filenames, probably the name of the assembled FASTA.")
parser.add_argument("-f",type=str,required=True,help="Fullpath to fasta.")
parser.add_argument("-q",type=str,required=True,help="Fullpath to raw reads.")
parser.add_argument("-g",type=str,required=True,help="Fullpath to reference gtf.")
parser.add_argument("-o",type=str,required=True,help="Fullpath to desired output directory.")
parser.add_argument("-h",type=str,required=True,help="Fullpath to directory containing Hisat2 databases.")
parser.add_argument("-s",type=str,required=True,help="Base name of Hisat2 databases (i.e. Dmel).")
parser.add_argument("-r",type=str,required=True,help="Fullpath to the directory that contains ref blast databases. You may need to adjust the versions in the script if they have been updated.")
parser.add_argument("-d",type=str,required=True,help="Fullpath to the directory that contains additional desired blast databases.")
parser.add_argument("-t",type=str,default=1,required=False,help="TPM cutoff.")
parser.add_argument("-p",type=int,default=80,required=False,help="Percent alignment...") 
parser.add_argument("-v",type=int,default=100,required=False,help="... over how many base pairs.")
parser.add_argument("-k",type=int,default=5000,required=False,help="Desired length (bp) around target region for synteny check.")

args = parser.parse_args()

def makebamswhisat(reads, hisatdb, hbase): # Make BAM files for each dataset with Hisat 2. 
    makebamdircmd = "mdkir bamdir"
    subprocess.run(makebamdircmd,shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    location = reads
    for file in os.listdir(reads):
        if file.endswith(".fastq.gz"):
            if "Undetermined" not in file:
                base = file.split('.')[0]
                print(base)
                base = '_'.join(base.split('_')[:3])
                hisatcmd = "hisat2 -p 20 --dta -q -x {3}{4} -1 {0}{1} -2 {0}{1} 2> bamdir/{2}.align.out | samtools view -b - 2> bamdir/{2}.view.out | samtools sort -o bamdir/{2}.bam -O bam 2> bamdir/{2}.sort.out".format(reads, file, base, hisatdbs, hbase)
                subprocess.run(hisatcmd,shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    
def customrunblast(fasta, base, output, refdbs): # Begun Lab / fly specific blast loops done here. 
    
    blastdbs=["dmel-chromosome-r6.41", "dsim-rna-V3"] 
    
    for species in ("mel","sim","ana","yak"): # Make lists of all of the reference dbs - some of these will be empty.
        if species == "mel":
            suffix = "r6.41" 
        elif species == "sim": 
            suffix = "r2.02"
        elif species == "ana": 
            suffix = "r1.06"
        else:
            suffix = "r1.05"

        for reference in ("-3prime-","-5prime-","-intergenic-","-intron-","-miRNA-","-miscRNA-","-ncRNA-","-pseudogene-","-transposon-","-tRNA-", "-CDS-"):
            ref = "d"+species+reference+suffix
            blastdbs.append(ref)
        
    for each in blastdbs:
        blastcmd = "blastn -query {0} -db {4}{1} -out {3}{2}.{1}.outfmt6.txt -outfmt 6".format(fasta, each, base, output, refdbs)
        subprocess.run(blastcmd,shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
def runblast(fasta, base, db, output): #General blast command - in this case, used for the YO databases.

    yoblastdbs= []
    
    # make list of all of the YO reference dbs using os scandir
    yotoblast = os.scandir(db)
    for entry in yotoblast:
        if entry.is_file():
            if ".nto" in str(entry) and "fasta" not in str(entry):
                if "melanogaster" not in str(entry):
                    yostuff = entry.name[0:-4]
                    yoblastdbs.append(yostuff)
    onlyyobds = list(set(yoblastdbs))

    for yo in onlyyobds:
        yoblastcmd = "blastn -query {0} -db {1}{2} -out {3}{4}_{2}_outfmt6.txt -outfmt 6".format(fasta, db, yo, output, base)
        subprocess.run(yoblastcmd,shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

def normalizeblastouts(output): #Adjust filenames for future steps.
    
    os.chdir(output)
    list = os.listdir(output)
    for each in list: 
        if "outfmt6.txt" not in each:
            list.remove(each)
    for every in list:
        for species in ["Drosophila_melanogaster","Drosophila_ananassae","Drosophila_mojavensis","Drosophila_persimilis","Drosophila_pseudoobscura","Drosophila_virilis","Drosophila_willistoni","Drosophila_yakuba"]:
            if species in every:
                content = every.split("_outfmt6.txt")[0].split("dmel_ALL_")[1]
                content = content.replace("_","-")
                content = content.replace(".","-")
                newstring = "dmel_ALL."+content+".outfmt6.txt"
                print(newstring)
                os.rename(every,newstring)

        else:
            for ref in ["dmel-","dana-","dyak-","dsim-"]:
                if ref in every:
                    content = every.split(".outfmt6.txt")[0].split("dmel_ALL.")[1]
                    content = content.replace("_","-")
                    content = content.replace(".","-")
                    newstring2 = "dmel_ALL."+content+".outfmt6.txt"
                    print(newstring2)
                    os.rename(every,newstring2)  
    
def updatesetparameter(output, base, perc, over, fasta): #Update set parameter to move through the tissues and species. Default is 80% alignment over 100bp. Produces .screened, .positions, and .transcripts outputs. 
    uspcmd = "perl /home/kdlombardo/code/kaescripts/femalerepo/updated_parse_trinity_transcripts.pl {0} {2} {3} {1} {0} {4}".format(output, base, perc, over, fasta) 
    subprocess.run(uspcmd,shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

def sortnscreen(base, perc, over, output): #Sort by chromosome and check that they're 500 bp from OTHER things- "base" refers to the name of the dataset that carries across files, from assembly onward. In this case, "dmel_ALL" is the "base" of the name for the pooled assembly. 
    screencmd = "perl /home/kdlombardo/hayleys/assemblies_and_blastout/scripts/screen_candidates_by_gene_distance.pl /data/julie/ref/Dmel_gene_500extend {3}{0}.{1}_{2}.positions {3}{0}.{1}_{2}.positions.500screened".format(base, perc, over, output)
    sortcmd = "perl /home/kdlombardo/hayleys/assemblies_and_blastout/scripts/sort_RALL_positions.pl {3}{0}.{1}_{2}.positions.500screened {3}{0}.{1}_{2}.positions.500screened.sorted".format(base, perc, over, output)
    subprocess.run(screencmd,shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    subprocess.run(sortcmd,shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

def maketable(base, perc, over, ouput): # Make a table with the exon and strand information. Length screening done here too.
    tablecmd = "perl /home/juliecridland/scripts/make_RALL_transcript_sets.pl {3}{0}.{1}_{2}.positions.500screened {3}{0}.{1}_{2}.positions.500screened.record".format(base, perc, over, ouput)
    subprocess.run(tablecmd,shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
def updategtfs(base, perc, over, output): # Make an updated set of gtf inputs to add to the previous gtf file to be used for a new Stringtie estimation.
    ugtfcmd = "perl /home/kdlombardo/hayleys/assemblies_and_blastout/scripts/make_updated_gtf_from_RALL.pl {3}{0}.{1}_{2}.positions.500screened.record T {3}{0}.{1}_{2}.gtf".format(base, perc, over, output)
    subprocess.run(ugtfcmd,shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

def catgtfs(base, rgtf, perc, over, output): # Concatenate to previously made gtf. Be sure to use gtf's fullpath.
        catgtfcmd = "cat {4}{0}.{2}_{3}.gtf {1} > {4}{0}.{2}_{3}.joined.gtf".format(base, rgtf, perc, over, output)
        subprocess.run(catgtfcmd,shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

def runstringtie(base, perc, over, output): # Use the new joined GTF and the bam files created above to estimate TPMs of the newly identified candidates.
    for tissue in ["PV","SR","ST"]:
        for sample in ["dmelf1", "justdmel"]: 
            bamdir = sample+"_"+tissue+"_bam/"
            thisbam =(glob.glob("bamdir/"))
            for each in thisbam:
                if ".bam" in each:
                    strain=each.split("/")[5]
                    strain=strain.split("_")[0]
                    thisabund =("abund/"+strain+"_"+bamdir[:-1]+ '.abund.tab')
                    stringtiecmd = "stringtie -p 20 -e {0} -G {6}{3}.{4}_{5}.joined.gtf -A {1} -o {6}{2}_{7}_{3}.final.gtf".format(each, thisabund, bamdir[:-5], base, perc, over, output, strain)
                    print(stringtiecmd)
                    abcmd = "mkdir abund"
                    subprocess.run(abcmd,shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                    subprocess.run(stringtiecmd,shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    
def maketabletwo(base, tpm, output): # This table uses TPM cutoffs - defaut 1.
    tabledir = "{1}tables_{0}".format(base, output)
    mttcmd = "perl /home/juliecridland/scripts/make_TPM_table_denovo.pl abund/ {2} {1}/{0}_TPM.{3}.table".format(base, tabledir, tpm) #provide path to directory that contains all abundance files
    if os.path.isdir(tabledir) == False:
        print("Making directory...")
        mktabledircmd = "mkdir {0}".format(tabledir)
        subprocess.run(mktabledircmd,shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print("Making tables...")
        subprocess.run(mttcmd,shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    else:
        print("{0} already exists, writing files to this directory will overwrite any relevant files there already. Proceed? (y/n)".format(tabledir))
        yes = {'yes', 'y'}
        no = {'no', 'n'}  
        done = False
        while not done:
            choice = input().lower()
            if choice in yes:
                print("Making tables...")
                subprocess.run(mttcmd,shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                print(mttcmd)
                return True
            elif choice in no:
                print("Cancelling...")
                return False
            else:
                print("Please respond with yes or no.") 
    
def makerecordtable(base, perc, over, tpm, output):  # Take the screened record and subset by things only in the TPM > 1 table. 
            rtcmd = "perl /home/juliecridland/scripts/subset_records_by_TPM_matches.pl {4}tables_{0}/{0}_TPM.{3}.table {4}{0}.{1}_{2}.positions.500screened.record {4}tables_{0}/{0}.{1}_{2}_candidates.record.table".format(base, perc, over, tpm, output)
            subprocess.run(rtcmd,shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # TABLE COLUMNS FROM LEFT TO RIGHT: AG, Transcript ID, Gene Length, Chromosome, Start, Stop, Exon Start, Exon Stop
    
def checkanddrop(base, perc, over, output): # Ensure that candidates are 500bp away from ONE ANOTHER and shorter than 300bp part - drop if not.
    tabledir = "tables_{0}".format(base)
    cdcmd = "perl /home/juliecridland/scripts/sort_distances_and_count_chrom_Female_reproductive.pl {4}{0}/{1}.{2}_{3}_candidates.record.table test_inter {4}{0}/sorted_{1}.{2}_{3}_candidates".format(tabledir, base, perc, over, output)
    subprocess.run(cdcmd,shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    
def makenewfasta(base, perc, over, output): #Create a new fasta.
    tabledir = "tables_{0}".format(base)
    mnfcmd = "perl /home/juliecridland/scripts/make_AG_denovo_fasta.pl {4}{1}/sorted_{0}.{2}_{3}_candidates {4}{0}.{2}_{3}.transcripts {4}{0}.{2}_{3}.Intergenic.TPM.1.fasta".format(base,tabledir, perc, over, output)
    subprocess.run(mnfcmd,shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

def finalblast(base, perc, over, output, refdbs): #Blast new fasta back to the mel genome reference and double check that there are no unsual BLAST hits. 
    dcdir = "final_{0}_blasts".format(base)
    if os.path.isdir(dcdir) == False:
        print("Making directory for blast results of final fasta...")
        mkdcdircmd = "mkdir {1}{0}".format(dcdir, output)
        subprocess.run(mkdcdircmd,shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
    blastdbs=["dmel-chromosome-r6.41", "dsim-rna-V3"] 
    species = "mel"
    suffix = "r6.41" 
    for reference in ("-3prime-","-5prime-","-intergenic-","-intron-","-miRNA-","-miscRNA-","-ncRNA-","-pseudogene-","-transposon-","-tRNA-", "-CDS-"):
        ref = "d"+species+reference+suffix
        blastdbs.append(ref)

        for each in blastdbs:
            fbcmd = "blastn -query {5}{1}.{3}_{4}.Intergenic.TPM.1.fasta -db {6}{0} -out {5}{2}/final_{1}.{0}.outfmt6.txt -outfmt 6".format(each, base, dcdir, perc, over, output, refdbs)
            subprocess.run(fbcmd,shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    
def syntenycheck(base, perc, over, output): # Uses the 2021 Flybase table of orthologs (the tsv). Intergenic ONLY. Goes through the sorted file, finds the nearest upstead and downstream genes and their orthologs, then collects those orthologs' FBGNs. (Align to other species' genomes.)
    syntcmd = "perl /home/juliecridland/scripts/get_synteny_positions_dmel.pl {3}tables_{0}/sorted_{0}.{1}_{2}_candidates /data/FlyRef/dmel_orthologs_in_drosophila_species_fb_2021_02.tsv {3}{0}_synteny_positions_and_orthologs".format(base, perc, over, output)
    subprocess.run(syntcmd,shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # Look up positions of orthologs in yak and sim.
    simcheckcmd = "perl /home/juliecridland/scripts/get_ortholog_positions_for_synteny.pl sim {1}{0}_synteny_positions_and_orthologs /data/FlyRef/Dsim_gene_ranges {1}{0}_Dsim_synteny_positions {1}{0}_Dsim_synteny_manual_check".format(base, output)
    yakcheckcmd = "perl /home/juliecridland/scripts/get_ortholog_positions_for_synteny.pl yak {1}{0}_synteny_positions_and_orthologs /data/FlyRef/Dyak_gene_ranges {1}{0}_Dyak_synteny_positions {1}{0}_Dyak_synteny_manual_check".format(base, output)
    subprocess.run(simcheckcmd,shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    subprocess.run(yakcheckcmd,shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    
def pullregionaroundcandidates(base, perc, over, tlen, output): # Makes a fasta by pulling out the regions around the candidates plus 5kb on each side. The fasta is from the mel choromosome file.
    pullcmd = "perl /home/juliecridland/scripts/get_region_X_around_candidates.pl {4}tables_{0}/sorted_{0}.{1}_{2}_candidates {3} /data/julie/AG_Denovo/dmel-all-chromosome-r6.19.concat.fasta {4}{0}_plus{3}bp.fasta {4}{0}_all_positions".format(base, perc, over, tlen, output)
    subprocess.run(pullcmd,shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    
def blastdoublecheck(base, tlen, output): # Checks that the candidates fall within the regions specified.
    sim2cmd = "blastn -db /data/FlyRef/blast/dsim-all-chromosome-r2.02 -query {2}{0}_plus{1}bp.fasta -out {2}{0}_plus{1}bp.to.Dsim.align -outfmt 6 -evalue 1e-6".format(base, tlen, output)
    yak2cmd = "blastn -db /data/FlyRef/blast/dyak-chromosome-r1.05 -query {2}{0}_plus{1}bp.fasta -out {2}{0}_plus{1}bp.to.Dyak.align -outfmt 6 -evalue 1e-6".format(base, tlen, output)
    subprocess.run(sim2cmd,shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    subprocess.run(yak2cmd,shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    
def getoutgrouppositions(base, output): # Get positions in outgroups of the orthologous matches.
    getposcmd = "perl /home/juliecridland/scripts/get_synteny_regions.pl {1}{0}_all_positions {1} {1}{0}_Synteny_ortholog_positions > {1}{0}_Synteny_ortholog_no_match".format(base, output)
    subprocess.run(getposcmd,shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    
def checknesting(base, output): # Check if the syntenic positions in the orthologs are where they're supposed to be.
    simnestcmd = "perl /home/juliecridland/scripts/compare_ortholog_positions_to_candidate_ortho_Female_reproductive.pl {1}{0}_Dsim_synteny_positions {1}{0}_Synteny_ortholog_positions {1}{0}_Dsim_synteny_confirmed.txt".format(base, output)
    yaknestcmd = "perl /home/juliecridland/scripts/compare_ortholog_positions_to_candidate_ortho_Female_reproductive.pl {1}{0}_Dyak_synteny_positions {1}{0}_Synteny_ortholog_positions {1}{0}_Dyak_synteny_confirmed.txt".format(base, output)
    subprocess.run(yaknestcmd,shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    subprocess.run(simnestcmd,shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

makebamswithhisat(args.r, args.h, args.s)                                
customrunblast(args.f, args.b, args.o, args.r)
runblast(args.f, args.b, args.r, args.d)                               
normalizeblastouts(args.o)
updatesetparameter(args.o, args.b, args.p, args.v, args.f)
sortnscreen(args.b, args.p, args.v, args.o)
maketable(args.b, args.p, args.v, args.o)
updategtfs(args.b, args.p, args.v, args.o)
catgtfs(args.b, args.g, args.p, args.v, args.o)
runstringtie(args.b, args.p, args.v, args.o)
maketabletwo(args.b, args.t, args.o)
makerecordtable(args.b, args.p, args.v, args.t, args.o)
checkanddrop(args.b, args.p, args.v, args.o)
makenewfasta(args.b, args.p, args.v, args.o)
finalblast(args.b, args.p, args.v, args.o, args.r)
syntenycheck(args.b, args.p, args.v, args.o)
pullregionaroundcandidates(args.b, args.p, args.v, args.k, args.o)
blastdoublecheck(args.b, args.k, args.o)
getoutgrouppositions(args.b, args.o)
checknesting(args.b, args.o)
