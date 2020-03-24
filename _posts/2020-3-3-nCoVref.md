---
layout: post
title: Using SARS-CoV-2 Reference Genome to Align Outbreak Sequencing Data
---

In this post we use the reference genome established by the CDC: [Severe acute respiratory syndrome coronavirus 2 isolate Wuhan-Hu-1, complete genome.](https://www.ncbi.nlm.nih.gov/nuccore/1798174254) To begin we need to grab just the *FASTA* data and save it as a file under something like **nCoVref.fasta** 

Next we will use samtools to index this reference file.

    git clone https://github.com/samtools/samtools.git
    
    cd samtools
    ./configure
    make
    make install
    
    samtools faidx nCoVref.fasta

This creates our nCoVref.fasta.fai index file, which will be useful later. To perform our actual alignment we are going to use the Burrow-Wheeler Aligner for short-reads, BWA. 

    git clone https://github.com/lh3/bwa.git
    cd bwa; make
    ./bwa index nCoVref.fasta
    
BWA index created some additional files, like .bwt, .pac, .ann, and others. Now our reference genome is ready to align reads from sequencing data of patients from outbreaks. We are going to use isolate raw read files from the Washington SARS-CoV-2 outbreak uploaded to the [Sequence Read Archive.](https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR11278092)To retrieve this data we need to use the SRA Toolkit, whose compiled binaries/install scripts of February 26, 2020 can be found [Here](https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software). Once installed we grab our reads:

    ./fastq-dump --outdir /path SRR11278092



![_config.yml]({{ site.baseurl }}/images/config.png)

