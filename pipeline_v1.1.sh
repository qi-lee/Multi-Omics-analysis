#!/bin/bash
###########################################################
#    File Name: pipeline.sh
#    Author: liqi
#    Mail: liqi@ihb.ac.cn
#    Copyright: Lab of Algal Genomics
#    Created Time: Fri 10 May 2020 01:55:54 PM CST
###########################################################

# default params
threads=20; # number of threads (default=4)
samplename="test"
workdir="Test"
out="Out"

if [ ! -d $out ]; then mkdir $out; fi
if [ ! -d $out/1_preprocess ]; then mkdir $out/1_preprocess; fi
if [ ! -d $out/2_assembly ]; then mkdir $out/2_assembly; fi
if [ ! -d $out/2_assembly/tmp ]; then mkdir $out/2_assembly/tmp; fi
if [ ! -d $out/3_kraken ]; then mkdir $out/3_kraken; fi
if [ ! -d $out/4_initial_bin ]; then mkdir $out/4_initial_bin; fi
if [ ! -d $out/4_bin_refinement ]; then mkdir $out/4_bin_refinement; fi
if [ ! -d $out/5_taxonomy ]; then mkdir $out/5_taxonomy; fi
if [ ! -d $out/6_visualize ]; then mkdir $out/6_visualize; fi
if [ ! -d $out/7_annotation ]; then mkdir $out/7_annotation; fi


for F in $name 
do
# Preprocess
    RAWREAD1=$workdir/RawReads/${F}_1.fq.gz
    RAWREAD2=$workdir/RawReads/${F}_2.fq.gz
    if [[ ! -s $out/1_preprocess/$F/${F}_final_R1.fastq ]]; then
        metawrap read_qc -1 $RAWREAD1 -2 $RAWREAD2 -t 20 --skip-bmtagger -o $out/1_preprocess/$F > $out/1_preprocess/${F}_readQC.log 2>&1

        cd $workdir/1_preprocess/$F/
        pigz -p 20 -k ${F}_final_R1.fastq
        pigz -p 20 -k ${F}_final_R2.fastq
        rm *_unpaired_trimmed.fq
        cd $workdir/
    else
        echo "Looks like the assembly was already existed. Skipping..."
    fi
    READ1=$workdir/1_preprocess/$F/${F}_final_R1.fastq
    READ2=$workdir/1_preprocess/$F/${F}_final_R2.fastq

# Assembly
    if [[ ! -s $out/2_assembly/$F/$F.contigs.fa ]]; then
        rm -r $out/2_assembly/$F
        megahit -1 $READ1 -2 $READ2 -o $out/2_assembly/$F -t $threads -m 0.5 --out-prefix $F --min-contig-len 500 --tmp-dir $out/2_assembly/tmp > $out/2_assembly/${F}_assembly.log 2>&1
        rm -r $out/2_assembly/tmp/*
        rm -r $out/2_assembly/$F/intermediate_contigs
    else
        echo "Looks like the assembly was already existed. Skipping..."
    fi

    Script/rm_short_contigs.pl -i $out/2_assembly/$F/$F.contigs.fa -l 1500 -o $out/2_assembly/$F/${F}_scaffolds_lg1500.fasta
    
    if [[ ! -s $out/2_assembly/$F/${F}_scaffolds_lg1500.fasta ]]; then 
        echo "#######################################"
        echo "Something went wrong with the assembly."
        echo "#######################################"
        continue
    fi

    Ctgs=$workdir/$out/2_assembly/$F/${F}_scaffolds_lg1500.fasta

# kraken2
    if [ ! -d $out/3_kraken/$F ]; then mkdir $out/3_kraken/$F; fi
        ln -s $READ1 $out/3_kraken/$F/${F}_1.fastq
        ln -s $READ2 $out/3_kraken/$F/${F}_2.fastq
    kraken2 --db minikraken2 --threads $threads --output $out/3_kraken/$F/${F}_reads.kraken2 --confidence 0.05 --report $out/3_kraken/$F/${F}_reads.kreport2 --paired $out/3_kraken/$F/${F}_1.fastq $out/3_kraken/$F/${F}_2.fastq
    ktImportTaxonomy -q 2 -t 3 -o $out/3_kraken/$F/${F}_kraken2.html $out/3_kraken/$F/${F}_reads.kraken2
    rm -rf $out/3_kraken/$F/*.html.files
    est_abundance.py -i $out/3_kraken/$F/${F}_reads.kreport2 -k database150mers.kmer_distrib -o $out/3_kraken/$F/${F}_reads.O.txt -l O
    est_abundance.py -i $out/3_kraken/$F/${F}_reads.kreport2 -k database150mers.kmer_distrib -o $out/3_kraken/$F/${F}_reads.P.txt -l P
    rm $out/3_kraken/$F/*.kraken2
    rm $out/3_kraken/$F/*.fastq

# Binning
    if [[ ! -s $out/4_initial_bin/$F/metabat2_bins/bin.unbinned.fa ]]; then
        mkdir $out/4_initial_bin/$F
        ln -s $READ1 $out/4_initial_bin/$F/${F}_1.fastq
        ln -s $READ2 $out/4_initial_bin/$F/${F}_2.fastq
        metawrap binning -o $out/4_initial_bin/$F -t $threads -a $Ctgs -m 350 -l 1500 --metabat2 --maxbin2 $out/4_initial_bin/$F/${F}_1.fastq $out/4_initial_bin/$F/${F}_2.fastq > $out/4_initial_bin/${F}_ibinning.log 2>&1
        rm $out/4_initial_bin/$F/${F}_1.fastq
        rm $out/4_initial_bin/$F/${F}_2.fastq
        rm $out/4_initial_bin/$F/work_files/${F}.bam*
        rm $out/4_initial_bin/$F/work_files/assembly*

        metawrap bin_refinement -o $out/4_bin_refinement/$F -A $out/4_initial_bin/$F/metabat2_bins/ -B $out/4_initial_bin/$F/maxbin2_bins/ -c 50 -x 10 -t 24 -m 350 > $out/4_bin_refinement/${F}_bin_refinement.log 2>&1
        rm -r $out/4_bin_refinement/$F/work_files
        rm -r $out/4_bin_refinement/$F/maxbin2_bins
        rm -r $out/4_bin_refinement/$F/metabat2_bins
    else
        echo "Looks like the bins were already existed. Skipping..."
    fi

    if [[ ! -s $out/4_bin_refinement/$F/metawrap_50_10_bins.stats ]]; then 
        echo "######################################"
        echo "Something went wrong with the binning."
        echo "######################################"
        continue
    fi

# Taxonomy
    metawrap classify_bins -b $out/4_bin_refinement/$F/metawrap_50_10_bins -o $out/5_taxonomy/$F -t $threads > $out/5_taxonomy/${F}_Classify_bins.log 2>&1
    mkdir $out/5_taxonomy/$F
    rm -r $out/5_taxonomy/$F/bins
    cp -rf $out/4_bin_refinement/$F/metawrap_50_10_bins $out/5_taxonomy/$F/bins
    for B in `ls 5_taxonomy/$F/bins/*.fa`
    do
        sample=${B##*/}
        mv $B $out/5_taxonomy/$F/bins/${F}_$sample
    done
    awk '{print "'$F'_"$0}' $out/5_taxonomy/$F/bin_taxonomy.tab > $out/5_taxonomy/$F/${F}_bin_taxonomy.txt
    rm $out/5_taxonomy/$F/all_contigs.fa
    rm $out/5_taxonomy/$F/mapping.tax

# Blobology
    metawrap blobology -a $Ctgs -t $threads -o $out/6_visualize/$F --bins $out/4_bin_refinement/$F/metawrap_50_10_bins $out/4_initial_bin/$F/${F}_1.fastq $out/4_initial_bin/$F/${F}_2.fastq > $out/6_visualize/${F}_blobology.log 2>&1
    metawrap quant_bins -a $Ctgs -t $threads -o $out/6_visualize/$F -b $out/4_bin_refinement/$F/metawrap_50_10_bins $out/4_initial_bin/$F/${F}_1.fastq $out/4_initial_bin/$F/${F}_2.fastq > $out/6_visualize/${F}_quant.log 2>&1
    rm $out/6_visualize/$F/*.bt2
    rm $out/6_visualize/$F/*.fasta
    rm -r $out/6_visualize/$F/assembly_index
    rm -r $out/6_visualize/$F/alignment_files
    awk '{print "'$F'_"$0}' $out/6_visualize/$F/bin_abundance_table.tab > $out/6_visualize/$F/${F}_bin_abundance.txt

# Annotation
    if [[ ! -d $out/7_annotation/$F/bin_translated_genes ]]; then
        metaWRAP annotate_bins -o $out/7_annotation/$F -t $threads -b $out/5_taxonomy/$F/bins -s $F > $out/7_annotation/${F}_annotation.log 2>&1
        Script/faa_files_rename.pl -i $workdir/$out/7_annotation/$F/bin_translated_genes -o $workdir/$out/7_annotation/$F/bin_faa_rename
        mv $workdir/$out/7_annotation/$F/bin_faa_rename/gene_stats.txt $workdir/$out/7_annotation/$F/${F}_bin_genes.txt
    fi

    #rm $READ1
    #rm $READ2
done

