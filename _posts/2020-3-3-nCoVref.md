---
layout: post
title: Using SARS-CoV-2 Reference Genome to Align Outbreak Sequencing Data
---

In this post we take the reference genome established by the CDC: [Severe acute respiratory syndrome coronavirus 2 isolate Wuhan-Hu-1, complete genome.](https://www.ncbi.nlm.nih.gov/nuccore/1798174254) To begin we need to grab just the *FASTA* data and save it as a file under something like **nCoVref.fasta** 

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
    
BWA index created some additional files, like .bwt, .pac, .ann, and others. Now our reference genome is ready to align reads from sequencing data of patients from outbreaks.

![_config.yml]({{ site.baseurl }}/images/config.png)

The easiest way to make your first post is to edit this one. Go into /_posts/ and update the Hello World markdown file. For more instructions head over to the [Jekyll Now repository](https://github.com/barryclark/jekyll-now) on GitHub.
