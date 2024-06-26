# A flexible pipeline for multi-omics data analysis

## This package is built around a collection of publicly available tools and personal scripts tied together for analyzing multi-omics datasets.

### The pipeline was developed by QiLi (liqi at ihb.ac.cn, Lab of Algal Genomics).<br> External software used by this pipeline are copyright respective authors.
<br>

### The pipeline can be broadly separated into seven main sections：
* 1. Preprocess
* 2. Metagenomics
* 3. Metatranscriptomics
* 4. Metaproteomics
* 5. Metabonomics
* 6. Phylogenetic analysis
* 7. Functional analysis
<br>

### * 1. Preprocess
  The presence of poor quality or technical sequences such as adapters in the sequencing data can easily result in suboptimal downstream analyses. There are many useful read preprocessing tools to perform the quality control (FastQC, Trimmomatic). Here we choose Trimmomatic to clean our sequencing datasets, for example:
  
    java -jar xx/software/trimmomatic-0.xx.jar PE -threads 30 -phred64 ${SAMPLE}_R1.fq ${SAMPLE}_R2.fq ${SAMPLE}_R1_paired_trimmed.fq ${SAMPLE}_R1_unpaired_trimmed.fq ${SAMPLE}_R2_paired_trimmed.fq ${SAMPLE}_R2_unpaired_trimmed.fq ILLUMINACLIP:xx/Trimmomatic-0.xx/adapters/TruSeqxxx.fa:2:30:10 LEADING:10 TRAILING:10 SLIDINGWINDOW:4:20 MINLEN:20

<br>

### * 2.1 Metagenomics/Assembly
  According to our experience, we use SPAdes Genomes Assembler to assembly our data. Executinng SPAdes as the following command:
  
    spades.py -o output -k 21,33,55 --pe1-1 ${SAMPLE}_R1_paired_trimmed.fq --pe1-2 ${SAMPLE}_R2_paired_trimmed.fq --pe1-s ${SAMPLE}_R1R2_unpaired_trimmed.fq --mp1-1 ${SAMPLE}_MP_R1_paired_trimmed.fq --mp1-2 ${SAMPLE}_MP_R2_paired_trimmed.fq --mp1-s ${SAMPLE}_MP_R1R2_unpaired_trimmed.fq --pacbio ${SAMPLE}_filtered_subreads.fastq --careful --cov-cutoff 3 --disable-gzip-output -t 20 -m 500

<br>

### * 2.2 Metagenomics/Binning
  In order to bin metagenomics sequence, we have to collect some fetures about the assembled fragments. The output of some assembler, such as SPAdes, Velvet and AByss, provide the coverage value of each cotigs in the header. If there isn’t coverage information in the assembly data, we can map reads to contigs using BWA or Bowtie to get coverage information instead. For each contigs, GC content, length and coverage are extracted from assembly file using my perl script. The script outputs the file contigs.info including 4 column, that is contig name, contig length, GC content and coverage. Additionly, we utilize the essential gene sets to draw back some assembly fragements from the unassigned sequences. The process of obtaining the information of essential gene contained in the assembly was doned by the same perl script.

    calc_contiginfo_essentialgene_v3.pl -r ${SAMPLE}_scaffolds.fasta -1 ${SAMPLE}_R1_paired_trimmed.fq -2 ${SAMPLE}_R2_paired_trimmed.fq -m 500 -I 0 -X 600 -p phred64 -t $THREADS
Binning the metagenomic assembly using selected features.

    Meta-binning.R
The process of selecting the contigs in each area may be time-consuming, it totally depends on the complexity of the metagenomics project. Sometimes, you need to combine more than one level value to get the optimal inital group.

<br>

### * 2.3 Metagenomics/Bins Evaluation
  Since genomes reconstructed from metagenomic data usually vary substantially in their qualities, we proposed a set of quality criteria with quantitative thresholds to evaluate the quality of these genomes for subsequent analyses.

![image](https://github.com/qi-lee/Meta-Microbiome/blob/master/5_Bins_Evaluation/quality.JPG)
Other related software: CheckM.    
<br>

### * 2.4 Taxonomy
  We implemented taxonomic assignments of the genome bins using TAXAassign with some modiﬁed codes for efﬁciency and accuracy. We used deduced amino acid sequence information through DIAMOND BLASTP searches, instead of nucleotide sequences through BLASTN searches, to produce a protein sequence alignment against the NCBI non-redundant (nr) protein database. 
    
    TAXAassign_prot.sh -c 30 -r 100 -t 60 -m 60 -q 50 -a "60,65,70,80,95,95" -f ${SAMPLE}_All_bins.faa

<br>


### * 3.1 Metatranscriptomics/RemovingrRNA
  Removing rRNA Sequences with SortMeRNA.
  <br>
  "SortMeRNA is a program tool for filtering, mapping and OTU-picking NGS reads in metatranscriptomic and metagenomic data. The core algorithm is based on approximate seeds and allows for fast and sensitive analyses of nucleotide sequences. The main application of SortMeRNA is filtering ribosomal RNA from metatranscriptomic data."

    conda install -c bioconda sortmerna --yes
    
    sortmerna -h
    sortmerna --ref rfam-5.8s-database-id98.fasta,rfam-5.8s-database-id98.idx:silva-bac-16s-id90.fasta,silva-bac-16s-id90.idx:silva-bac-23s-id98.fasta,silva-bac-23s-id98.idx --reads ${SAMPLE}_trimmed.fq --aligned ${SAMPLE}_aligned --other ${SAMPLE}_other --fastx --blast 1 --log --num_alignments 1 --paired_in -v -a 40 -e 1e-5 --num_seeds 3

<br>

### * 3.2 Metatranscriptomics/Aligning
  Aligning to Genome with STAR-aligner.
  <br>
  "The STAR aligner is a very fast and efficent spliced aligner tools for aligning RNAseq data to genomes."
    
    conda install -c bioconda star --yes
    
    STAR -h
    STAR --genomeDir star_index --readFilesIn ${SAMPLE}_filtered.fq  --runThreadN 20 --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts

<br>

### * 3.3 Metatranscriptomics/GeneCounts
  Summarizing Gene Counts with featureCounts.
  <br>
  "featureCounts is a highly efficient general-purpose read summarization program that counts mapped reads for genomic features such as genes, exons, promoter, gene bodies, genomic bins and chromosomal locations. It can be used to count both RNA-seq and genomic DNA-seq reads. featureCounts takes as input SAM/BAM files and an annotation file including chromosomal coordinates of features."
    
    conda install -c bioconda subread --yes
    
    featureCounts -h
    featureCounts -T 20 -p -t gene -a genomic.gtf -o ${SAMPLE}_final_counts.txt ${SAMPLE}.bam

<br>

### * 3.4 Metatranscriptomics/DEG
  Importing Gene Counts into RStudio.
  <br>
  "you can now use the gene count table as an input into DESeq2 for statistical analysis using the R-programming language."

    library(DESeq2)
    library(ggplot2)
    library(pheatmap)
    library(dplyr)

    countdata <- read.table("${SAMPLE}_final_counts.txt", header = TRUE, skip = 1, row.names = 1)
    metadata <- read.delim("metadata.txt", row.names = 1)
    ddsMat <- DESeqDataSetFromMatrix(countData = countdata,colData = metadata,design = ~Group)
    ddsMat <- DESeq(ddsMat)
    results <- results(ddsMat, pAdjustMethod = "fdr", alpha = 0.05)
    summary(results)
    write.table(x = as.data.frame(counts(ddsMat), normalized = T), file = 'normalized_counts.txt', sep = '\t', quote = F,col.names = NA)
    write.table(x = counts(ddsMat[row.names(results_sig)], normalized = T), file = 'normalized_counts_significant.txt', sep = '\t', quote = F, col.names = NA)

<br>

### * 4.1 Metaproteomics/MetaProteomeAnalyzer
  MetaProteomeAnalyzer (MPA)
  <br>
  "MetaProteomeAnalyzer (MPA) software for analyzing and visualizing MS-based metaproteomics data."

    conda install mpa-portable -c bioconda
    java -cp mpa-portable-X.Y.Z.jar de.mpa.cli.CmdLineInterface [parameters]
    java -cp mpa-portable-X.Y.Z.jar de.mpa.cli.CmdLineInterface -spectrum_files spectrum_file.mgf -database uniprot_sprot.fasta -missed_cleav 1 -prec_tol 10ppm -frag_tol 0.5Da -output_folder output

<br>

### * 4.2 Metaproteomics/ProteoStorm
  ProteoStorm
  <br>
  "ProteoStorm: An Ultrafast Metaproteomics Database Search Framework."
    
    wget https://repo.continuum.io/archive/Anaconda2-5.1.0-Linux-x86_64.sh
    bash Anaconda2-5.1.0-Linux-x86_64.sh

    python -u ProteoStorm.py --Database ./fasta --Spectra ./mgf --RemoveSpectra ./HS_matched_spectra.txt --SpectralDataset "demo" --output ./ProteoStorm_Out --PrecursorMassTolerance 10 --FragmentMassTolerance 0.015 --InstrumentID 3 --FragmentMethodID 3

<br>

### * 5.1 Metabonomics/NOREVA
  NOREVA
  <br>
  "R Package for Systematic Optimization of Metabolomic Data Processing."

    install.packages("devtools")
    devtools::install_github("idrblab/NOREVA")
    library(NOREVA)
    PrepareInuputFiles(dataformat, rawdata, label)
    normulticlassqcall(fileName, SAalpha="Y", SAbeta="Y", SAgamma="Y")
    normulticlassnoall(fileName, SAalpha="Y", SAbeta="Y", SAgamma="Y")
    normulticlassisall(fileName, IS)

<br>

### * 5.2 Metabonomics/MetaboAnalyst
  MetaboAnalyst
  <br>
  "MetaboAnalyst is a web-based platform dedicated for comprehensive metabolomics data analysis, interpretation and integration with other omics data."
  
    www.metaboanalyst.ca/
    new.metaboanalyst.ca/docs/Tutorials.xhtml

<br>


### * 6. Phylogenetic analysis
To infer phylogenetic relationships among bacteria, a whole-genome based and alignment-free Composition Vector Tree (CVTree) method was applied to the comparison and clustering of the 68 genomes we extracted from our assemblies.

    cvtree -i species.list -p data -o CVTree_k6.txt -k 6
    neighbor
Other related software: PHYLIP, MEGAN.
<br>

### * 7. Functional analysis

Functional characterization and annotation of protein encoding genes were performed by MOCAT2. To further compare the functional potential of each group, the predicted ORFs were analyzed using the GhostKOALA service on the KEGG website. When the results were returned through your email, open the links and then compare your interesting pathway.
    
    http://www.kegg.jp/ghostkoala/
Other related software: MEGAN.

<br>


### Dependencies

    FastQC
    Trimmomatic
    SPAdes
    Meta-IDBA
    CheckM
    TAXAassign
    SortMeRNA
    STAR
    featureCounts
    DESeq2
    MetaProteomeAnalyzer
    ProteoStorm
    NOREVA
    MetaboAnalyst
    CVTree
    phylip
    GhostKOALA
    KEGG
    MOCAT2
    MEGAN
    

### Reference:
    Andrews S. FastQC: a quality control tool for high throughput sequence data[J]. 2010.
    Bankevich A, Nurk S, Antipov D, et al. SPAdes: a new genome assembly algorithm and its applications to single-cell sequencing[J]. Journal of computational biology, 2012, 19(5): 455-477.
    Beyter and Lin et al. (2018). ProteoStorm: An Ultrafast Metaproteomics Database Search Framework. Cell Systems. 7, 463–467
    Bolger A M, Lohse M, Usadel B. Trimmomatic: a flexible trimmer for Illumina sequence data[J]. Bioinformatics, 2014, 30(15): 2114-2120.
    Buchfink B, Xie C, Huson D H. Fast and sensitive protein alignment using DIAMOND[J]. Nature methods, 2015, 12(1): 59-60.
    Fu J, Zhang Y, Wang Y, et al. Optimization of metabolomic data processing using NOREVA[J]. Nature protocols, 2022, 17(1): 129-151.
    Huson D H, Beier S, Flade I, et al. MEGAN community edition-interactive exploration and analysis of large-scale microbiome sequencing data[J]. PLoS computational biology, 2016, 12(6): e1004957.
    Ijaz U, Quince C. TAXAassign v0. 4[J]. 2013.
    Kanehisa M, Goto S. KEGG: kyoto encyclopedia of genes and genomes[J]. Nucleic acids research, 2000, 28(1): 27-30.
    Kanehisa M, Sato Y, Morishima K. BlastKOALA and GhostKOALA: KEGG tools for functional characterization of genome and metagenome sequences[J]. Journal of molecular biology, 2016, 428(4): 726-731.
    Kultima J R, Coelho L P, Forslund K, et al. MOCAT2: a metagenomic assembly, annotation and profiling framework[J]. Bioinformatics, 2016, 32(16): 2520-2523.
    Muth T, Behne A, Heyer R, et al. The MetaProteomeAnalyzer: a powerful open-source software suite for metaproteomics data analysis and interpretation[J]. Journal of proteome research, 2015, 14(3): 1557-1565.
    Muth T, Kohrs F, Heyer R, et al. MPA portable: a stand-alone software package for analyzing metaproteome samples on the go[J]. Analytical chemistry, 2018, 90(1): 685-689.
    Pang Z, Lu Y, Zhou G, et al. MetaboAnalyst 6.0: towards a unified platform for metabolomics data processing, analysis and interpretation[J]. Nucleic Acids Research, 2024: gkae253.
    Parks D H, Imelfort M, Skennerton C T, et al. CheckM: assessing the quality of microbial genomes recovered from isolates, single cells, and metagenomes[J]. Genome research, 2015, 25(7): 1043-1055.
    Peng Y, Leung H C M, Yiu S M, et al. Meta-IDBA: a de Novo assembler for metagenomic data[J]. Bioinformatics, 2011, 27(13): i94-i101.
    Plotree D, Plotgram D. PHYLIP-phylogeny inference package (version 3.2)[J]. cladistics, 1989, 5(163): 6.
    Qi J, Luo H, Hao B. CVTree: a phylogenetic tree reconstruction tool based on whole genomes[J]. Nucleic acids research, 2004, 32(suppl_2): W45-W47.
    Quinlan A R, Hall I M. BEDTools: a flexible suite of utilities for comparing genomic features[J]. Bioinformatics, 2010, 26(6): 841-842.
    Xie M, Ren M, Yang C, et al. Metagenomic analysis reveals symbiotic relationship among bacteria in microcystis-dominated community[J]. Frontiers in microbiology, 2016, 7.
    Xu Z, Hao B. CVTree update: a newly designed phylogenetic study platform using composition vectors and whole genomes[J]. Nucleic acids research, 2009, 37(suppl_2): W174-W178.
<br>

#### Questions or Comments, please contact: liqi at ihb.ac.cn.
