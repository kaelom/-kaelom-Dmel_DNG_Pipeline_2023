# De Novo Gene Identification Pipeline

This repository contains the perl scripts, written by Julie M. Cridland, used in 'Identifying candidate de novo genes expressed in the somatic female reproductive tract of *Drosophila melanogaster*.' The pipeline summarized in the wrapper is a modified version of the one used in ['Population biology of accessory gland-expressed de novo genes in Drosophila melanogaster'](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8733444/pdf/iyab207.pdf).
Scripts are numbered in the order they are used in the python wrapper, __wrapper.py__. This wrapper includes some intermediate steps that are not included in the scripts themselves.

## ![Pipeline Flowchart](https://raw.githubusercontent.com/kaelom/Dmel_DNG_Pipeline_2023/main/flowchartwithlegend.png)

1. updated_parse_trinity_transcripts.pl
2. screen_candidates_by_gene_distance.pl
3. sort_RALL_positions.pl
4. make_RALL_transcript_sets.pl
5. make_updated_gtf_from_RALL.pl
6. make_TPM_table_denovo.pl
7. subset_records_by_TPM_matches.pl
8. sort_distances_and_count_chrom_Female_reproductive.pl
9. make_AG_denovo_fasta.pl
10. get_synteny_positions_dmel.pl
11. get_ortholog_positions_for_synteny.pl
12. get_region_X_around_candidates.pl
13. get_synteny_regions.pl
14. compare_ortholog_positions_to_candidate_ortho_Female_reproductive.pl


While this pipeline is focused on processing *Drosophila melanogaster* transcriptome data specifically, minor alterations, such as adjusting necessary BLAST databases, can be made to use this pipeline on other species. 

## Wrapper Manual

The wrapper was utilized by calling it from the command line like so:

``` 
./wrapper.py -b basename -f path_to_query.fasta -q path_to_raw_read_dir/ -g path_to_reference.gtf -o output_directory/ -h path_to_hisat_db_dir/ -s basename_of_hisat2_dbs -r path_to_ref_blast_db_dir/ -d path_to_Yang_et_al_db_dir/

``` 

In addition to the perl scripts above, the following programs were also used:
bash, BLAST, HISAT2, Stringtie

Other input files:
Assembled transcript.
BAM files for the transcripts of interest.
All required databases used with BLAST/HISAT2.


### Required Flags:

* -b: Base name used in all the filenames, probably the name of the assembled FASTA. String.
* -f: Fullpath to fasta. String.
* -q: Fullpath to the directory that contains the raw reads. String.
* -g: Fullpath to reference gtf. String.
* -o: Fullpath to desired output directory. String.
* -r: Fullpath to the directory that contains mel ref sequences. You may need to adjust the versions in the script if they have been updated. String.
* -d: Fullpath to the directory that contains additional desired blast databases. String.
* -a: Fullpath to abundance file directory. String.

### Optional Flags:

* -t: TPM cutoff. String. 
* -p: Percent alignment used. Default = 80. Integer.
* -v: Number of base pairs over which the the percent alignment is used. Default = 100. Integer. 
* -k: Desired length (bp) around target region for synteny check. Default = 5kb. Integer
