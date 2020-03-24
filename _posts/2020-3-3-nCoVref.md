---
layout: post
title: Using SARS-CoV-2 Reference Genome to Align Outbreak Sequencing Data
---

In this post we take the reference genome established by the CDC at <https://www.ncbi.nlm.nih.gov/nuccore/1798174254> 
To begin we need to grab just the *FASTA* data and save it as a file under something like **nCoVref.fasta**
Next we will use samtools to index this reference file.

git clone https://github.com/samtools/samtools.git
samtools faidx nCoVref.fasta

This creates our nCoVref.fasta.fai index.

![_config.yml]({{ site.baseurl }}/images/config.png)

The easiest way to make your first post is to edit this one. Go into /_posts/ and update the Hello World markdown file. For more instructions head over to the [Jekyll Now repository](https://github.com/barryclark/jekyll-now) on GitHub.
