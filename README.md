# **Introduction to ChIP-seq analysis**

Welcome in this tutorial, you will find here a detailed ChIP-seq workflow, starting from sequencing read to the final coverage tracks and differentially accessible genomic regions.

This tutorial is using the Galaxy platform to perform the data download, quality control, mapping and peak calling. We will then explore the result via IGV and RSAT.

A small guide for this course :
  * ⚡️ : time to shine, this is your hands-on objective
  * 👀 : unscroll a help note if you're stuck
  * ❓ : quiz time!
  * 🪐 : this is a how-to referring to Galaxy
  * 🪩 : feeling on fire? try this out

Usefull links :
  * the [Myers *et al.*](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1003565) article
  * The [NCBI website](https://www.ncbi.nlm.nih.gov/)
  * The [GATK help page](https://gatk.broadinstitute.org/hc/en-us/articles/360035890791-SAM-or-BAM-or-CRAM-Mapped-sequence-data-formats) on alignment formats

## 📍 **1. A brief note on ChIP-seq**

ChIP-seq stands for **Ch**romatin **I**mmuno **P**recipitation followed by **seq**uencing.

## 📍 **2. A brief note on Galaxy**

Today, we will work on the Galaxy platform. It's simple, free and open-source.

## 📍 **3. Let's start the analysis : loading the raw data**

### 🔸 **3.a Find the identifier**

We will work on the study from [Myers *et al.*](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1003565).
[ADD STUDY DETAILS]

Upon publication of their work, authors should deposit their raw data on a publicly available repositories. You can access and download these archives via two main platforms, the Sequence Read Archive (SRA) from NCBI (US) and the European Nucleotide Archive (ENA) from EBI (EU). Both platforms regularly cross-update each other.

To find the correct accession ID in a study, you should look for the following :
* A **BioProject accession**, starting with `PRJ` (*e.g.* `PRJNA176146`), that will link to the complete project archive
* A **GEO identifier**, starting with `GSE` (*e.g.* `GSE41186`), that will link to a specific experiments, in our case ChIP-seq.

⚡️ Your turn : Find the accession identifiers for raw sequencing data from X et al.

<details>
  <summary>Tips 👀</summary>

  > You are looking for a code starting with `GSE`. You usually find it in the *Data accessibility* section of an article, else you can try to `Ctrl+F` for `GSE` in the paper.

  </details>
<br>

### 🔸 **3.b Load the raw data to Galaxy**

You can see in this project that multiple experiments were performed (ChIP-on-ChIP, ChIP-seq, RNA-seq). For the sake of time and simplicity, we will focus our tutorial on the following two ChIP-seq samples :

* ChIP-seq of the FNR protein in anaerobic condition, sample A
* Input DNA in anaerobic condition


❓ Can you guess why we are selecting this pair of dataset?

❓ Are the sequencing data single-end or paired-end?

⚡️ Find the SRA identifier (starts with `SRR`) of these two samples and upload them to Galaxy. Assign them a clear name (*e.g.* **FNR** and **Input**).

* 🪐 **Get Data** : *Download and Extract Reads in FASTA/Q*
* 🪐 Assign a new name to a sample :  edit the **Name** attribute via the ✏️`Edit attributes` link and save. <br><div style="text-align:center"><img src="image/chap3/edit_name.png" width="220"/></div>

<details>
  <summary>Tips 👀</summary>

  > The two sample's identifier are `SRR576933` (FNR ChIP) & `SRR576938` (Input).<br>
  ><div style="text-align:center"><img src="image/chap3/find_srr1.png" width="400"/></div>
  ><div style="text-align:center"><img src="image/chap3/find_srr2.png" width="400"/></div>
  >🪐 Paste the SSR identifier in Galaxy's tool and click `Execute`. The job will start running and turn green once finished.<br>
  > <div style="text-align:center"><img src="image/chap3/load_data.png" width="400"/></div>
  > Once finished, edit the name for both and group them as a collection (see below).

  </details>

<br>

### 🔸 **3.c Group into a Dataset Collection**

We will group both sample into a **Dataset Collection**. Working on a collection allows to perform the same type of command on both samples in parrallel in a single go.

⚡️ Make a Dataset Collection with both samples

* 🪐 Click the `Operations on multiple datasets` button in the History and select the two renamed input files
<div style="text-align:center"><img src="image/chap3/select_dc.png" width="200"/></div>

* 🪐 Select the `Build Dataset List` command and name the collection
<div style="text-align:center"><img src="image/chap3/build_dc.png" width="200"/></div>

* 🪐 Now you can run all tools on the two files in batch using the `Dataset collection` input 🎉
<div style="text-align:center"><img src="image/chap3/run_dc.png" width="300"/></div>



## 📍 **4. Sweep up the dust : cleaning raw data**

We need to assess the quality of the sequencing reads. We need to know the confidence in the sequence call and the potential presence of contaminants (index or other organism).

### 🔸 **4.a Quality Check of raw data**

FastQC in the most common tool used to get an overview of a fastq file's quality. We will run it with default settings.

⚡️ Run FASTQC on the Dataset Collection with both samples

* 🪐 **FASTQC** : *Reads Quality Report*
* 🪐 Select the Dataset Collection as input
* 🪐 Explore the result with the webpage output (`view data`)

❓ Do you see any problem with the quality of the dataset?

<details>
  <summary>Tips 👀</summary>

  > FastQC output comprises several parts, you can refer to the [tool's website](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) for detailed information. Take a close look to the **Overrepresented sequences** in FNR sample.

  </details>
<br>

### 🔸 **4.b Trimming contaminants**

We have seen in the previous chapter that **29%** of our FNR raw data correspond an Illumina Adapter sequence.

❓ What's an adapter sequence again ?

From here we have two choices :

* **Leave the contaminant sequence** in the dataset and assuming correct read mapping will filter them out
* **Remove the contaminant sequence** from out dataset set prior to mapping

We will hop for the **2nd** choice, as it will provide more acurrate mapping statitics that could help us detect other potential problem in the dataset that prevent a high fraction of mapped reads.

 ⚡️ Trim the contaminant Index sequence on both samples with Trimmomatic

* 🪐 Copy the Index sequence from the FASTQ output
* 🪐 **Trimmomatic** : *flexible read trimming tool for Illumina NGS data*
* 🪐 Use the correct data type (single/paired end)
* 🪐 Set the accuracy of the match between adapter to 8.
* 🪐 Save the log output <div style="text-align:center"><img src="image/chap4/out_trim.png" width="250"/></div>

High accuracy thresholds for short reads will remove adapter dimers but adapter contamination at the 3'end of the reads will remain undetected, this is why we lower the accuracy to 8. A threshold of 8 corresponds to 12 perfect matches between read and adapter, so adapter contamination of 12bp and longer will be detected.

<details>
  <summary>Tips 👀</summary>

  > The overrepresented Adapter Sequence is an Illumina Index 5 :
  ```
>TruSeq Adapter
GATCGGAAGAGCACACGTCTGAACTCCAGTCACA
```
> We can specifically select it as a custom fasta filter :
> <div style="text-align:center"><img src="image/chap4/trim_run.png" width="450"/></div>
> <div style="text-align:center"><img src="image/chap4/trim_param.png" width="450"/></div>

  </details>
<br>

### 🔸 **4.c Checking trimming**

You can access the trimming statistics in the output log file. Take a look.

❓ How many reads have been discared for each dataset ? How many reads do we have left?

An important question to ask ourself is whehter we have sequenced our sample deep enough. This will depend on two modalities :

* Have we saturated the sequencing library?
* Do we expect to cover the majority of the genome?

❓ How can we solve these two above questions?

<details>
  <summary>Tips 👀</summary>

  > Library is saturated when you have sequenced (almost) all the spot in the lane and you mostly increase duplication level.<br>

  > We must compare the number of good quality reads with the genome size of our organism (*e.g.* for the 3 Gb human genome, 10 million reads are considered sufficient). To do so, take a look at *Escherichia coli*'s genome information on NCBI.

  </details>
<br>

Everything is almost ready, make sure your data now looks all right with a new FASTQC run on trimmed data.

⚡️ Run FASTQC on the Dataset Collection **after** trimming

❓ Is the Adapter sequenced properly filtered out?

🪩 You can merge multiple FASTQC output in a fancy way with **MultiQC**, try it out!

## 📍 **5. Gotta map them all : read alignment**

Mapping is a step where you can explore a very large parameter space. We need to be careful about our selected settings.

### 🔸 **5.a Loading a reference genome**

First thing first, what map do we use? We can check that in the article's text.

❓ What is the model organism used?

❓ What is the reference genome of the study?

🪩 Try to fetch its information back from the [NCBI website](https://www.ncbi.nlm.nih.gov/)

<details>
  <summary>Tips 👀</summary>

  > You can see the sentence *Resulting reads were aligned to the published E. coli K-12 MG1655 genome (U00096.2)* in the Materials and Methods section

  </details>
<br>

This reference comes from NCBI Nucleotide database, which is not ideal (large amount of errors in sequences and annotations). It is now highly recommended to work with **NCBI Reference Sequence** (RefSeq) data, the curated subset of NCBI’s Nucleotide database. Similar to sample identifier, reference genome also have a unique ID.

⚡️ Find E.Coli RefSeq genome ID and upload it to Galaxy
* 🪐 **NCBI Accession Download** : *Download sequences from GenBank/RefSeq by accession through the NCBI ENTREZ API*
* Paste the RefSeq ID as a direct entry, keep output as fasta

<details>
  <summary>Tips 👀</summary>

  > The RefSeq ID for the E.Coli strain K-12 MG1655 is `NC_000913.3`. You can find it by searching this strain on NCBI's nucleotide database (1), filtering for RefSeq match only (2) and looking at the 1st entry (3).
  ><div style="text-align:center"><img src="image/chap5/get_refseq.png" width="500"/></div>
  > You can then paste it to Galaxy genome loader
  ><div style="text-align:center"><img src="image/chap5/load_refseq.png" width="500"/></div>

  </details>
<br>


### 🔸 **5.b Selecting mapping parameters**

There's plenty of mappers available, you probably have already heard of STAR, Bowtie and BWA for example. We need to select the appropriate one for our data.

We have **short** (36bp) and **single-end** reads, this is best suited for **Bowtie**.

Now the tricly part : understanding and playing with the parameters.

⚡️ Take a look at the Bowtie **-v** and **-n** parameters in Galaxy parameter list.
* 🪐 **Map with Bowtie for Illumina**

❓ What is the `-v 2` parameter doing?

❓ What are the `-n 2 -l 35` parameters doing?

❓ Make a guess : which of the two options would you pick here?

🪩 Do you see any other parameter that would have been useful to us?

<details>
  <summary>Tips 👀</summary>

  > The Bowtie **v-mode** performs end-to-end alignment of the read and filter out the ones with more than *n* mismatches. <br>
  > The Bowtie **n-mode** aligns a seed of *l* basepairs from the high quality 5'end and remove the ones with more than *n* mismatches. <br>
  > Remember your FASTQC reports, does the **Per base sequence quality** graph hint a something you would want to exclude? <br>

  </details>
<br>

These two parameters modify the **alignment** strategy, we also need to take care of the **reads reporting** strategy.

⚡️ Take a look at the Bowtie **-m** and **--un** parameters in Galaxy parameter list.

❓ What are the `-m 1` or `-m 3` parameters doing?

❓ What is the `--un` parameter doing?

❓ Make a guess : which parameter set would you pick here?

<details>
  <summary>Tips 👀</summary>

  > **-m** discards reads that map to more than *n* reported genome locations <br>
  > **--un** write the reads that failed to align in a separate file, think about why this would be useful. <br>

  </details>
<br>

### 🔸 **5.c Mapping the reads with Bowtie**

Now that we have explored the Bowtie parameter space, let's map our reads to the E.Coli reference genome.

⚡️ Map both trimmed samples with Bowtie on E.Coli RefSeq genome
* 🪐 **Map with Bowtie for Illumina**
* 🪐 Genome index is done before mapping, keep default settings
* 🪐 Expand Bowtie full parameter list
* 🪐 Perform a **seed-based** alignment with a seed of **35bp** and a maximum of **2 mismatches**.
* 🪐 **Discard** all the reads that map to multiple positions
* 🪐  Save the bowtie mapping statistics to the history

<details>
  <summary>Tips 👀</summary>

  > Select the tool **Map with Bowtie for Illumina** and fill the parameter form. Select the RefSeq genome as a single file and not a collection.
  ><div style="text-align:center"><img src="image/chap5/map_reads.png" width="500"/></div>
  > Expand Bowtie parameter form
  ><div style="text-align:center"><img src="image/chap5/full_list.png" width="500"/></div>
  > Select n-mode mapping with a 35bp seed and max 2 mismatches
  ><div style="text-align:center"><img src="image/chap5/n2_l36.png" width="500"/></div>
  > Remove multi-mapping reads
  ><div style="text-align:center"><img src="image/chap5/m1.png" width="500"/></div>
  > Save the mapping output to Galaxy history
  ><div style="text-align:center"><img src="image/chap5/save_hist.png" width="300"/></div>

  </details>
<br>


⚡️ After reads have finished mapping, check the mapping statistics output.

❓ How many alignment are reported per sample?

❓ Does the mapping ratio seem good enough?

🪩 Check the mapping statistics with **Samtools flagstat** : *tabulate descriptive stats for BAM datset*

🪩 Run another mapping job with other Bowtie parameters or the untrimmed samples, how does is affect the mapping performances?

<details>
  <summary>Tips 👀</summary>

  > Take a look at the *mapping stats* output of Bowtie, scroll down the short read warning down to the global stats.
  > You should obtain the following result :
  > * Input sample has **6217172** alignments with **94.67%** mapped reads.
  > * FNR sample has **2173747** alignments with **89.98%** mapped reads<br>

  </details>
<br>

Everything good on the mapping statistics side, let's take a closer look a the *mapped reads* output now.

❓ What is the alignment file format?

❓ What is the **SO** (sorting order)?

<details>
  <summary>Tips 👀</summary>

  > The output mapped reads are stored in a **SAM** file format. You can see this in the **Attributes** section of the data. For more detail on this format, check the [GATK help page](https://gatk.broadinstitute.org/hc/en-us/articles/360035890791-SAM-or-BAM-or-CRAM-Mapped-sequence-data-formats) and the image below (credit : D. Caetano-Anolles).
  ><div style="text-align:center"><img src="image/chap5/sam_bam.png" width="500"/></div>
  > The SAM files are **Unsorted** (*i.e.* reads are unordered), you can see this is the SAM header.
  ><div style="text-align:center"><img src="image/chap5/sam_header.png" width="500"/></div>

  </details>
<br>

⚡️ Sort the mapping output by genomic coordinates and save it as a compressed BAM output

* 🪐 **Samtools sort** : *order of storing aligned sequences*
* Use the *coordinate* sorting key
* Make sure the output is a coordinate-sorted **BAM** file

## 📍 **6. Reaching the summits : Peak calling**

### 🔸 **6.a A brief note on Peak Calling algorithm**

Here we are!

We have our data ready, we can now get to the detection of biological signal by detecting genomic position with enriched ChIP signal.

The most widely used tools for this is the good old **MACS** (**M**odel-based **A**nalysis for **C**hIP-**S**eq) tool. It is used for ChIP-seq as well as for ATAC-seq. It can robustly detect peaks of different profiles (*e.g.* broad or narrow). Current version is MACS2, beta MACS3 is on its way. Check its [GitHub repo](https://github.com/macs3-project/MACS) for more details.

❓ What's behind a peak calling algorithm?

The goal of MACS2 is to **detect region significantly enriched for chromatin immunoprecipitation signal** (our FNR IP) **compared to a given background**. The input sample does exactly that : it gives a measure of the background signal obtain without a targetted IP.


### 🔸 **6.b Extracting FNR and Input sample from the Collection**

Now that you know what MACS2 need to do, select the proper input files. We need to specifically assign both samples as **treatment** and **control** files, so we need to extract them separately from the collection.

⚡️ Extract FNR and Input sorted BAM files from the Dataset collection
* 🪐 **Extract Dataset** *from a list*
* Select the element identifier (*i.e* the name you assigned)
* Run once per sample and check the output file attributes

<div style="text-align:center"><img src="image/chap6/extract_collection.png" width="500"/></div>

### 🔸 **6.c Running MACS2 callpeak**

We have our BAM files ready, now let's take a look at the required parameters of MACS2 on the Galaxy platform.

* 🪐 **MACS2 callpeak** : *Call peaks from alignment results*

❓ What is the effective genome size and why is it needed?

❓ Should we build a shifting model?
