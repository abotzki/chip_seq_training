# **Introduction to ChIP-seq analysis**

Welcome in this tutorial, you will find here a detailed ChIP-seq workflow, starting from sequencing read to the final coverage tracks and differentially accessible genomic regions.

This tutorial is using the Galaxy platform to perform the data download, quality control, mapping and peak calling. We will then explore the result via IGV and RSAT.

A small guide for this course :
  * ⚡️ : time to shine, this is your hands-on objective
  * 👀 : unscroll a help note if you're stuck
  * ❓ : quiz time!
  * 🪩 : feeling on fire? try this out
  * 🪐 : this is referring to Galaxy
  * 🧬 : this is referring to IGV
  * 🔮 : this is referring to RSAT

Usefull links :
  * the [Myers *et al.*](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1003565) article
  * The [NCBI website](https://www.ncbi.nlm.nih.gov/)
  * The [GATK help page](https://gatk.broadinstitute.org/hc/en-us/articles/360035890791-SAM-or-BAM-or-CRAM-Mapped-sequence-data-formats) on alignment formats
  * The [UCSC Genome Browser User Guide](https://genome.ucsc.edu/goldenPath/help/hgTracksHelp.html)

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

⚡️ Your turn : Find the accession identifiers for raw sequencing data from Myers *et al.*.

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

❓ What are duplicate tags and how MACS2 treat them? Remember you also have duplicate informatiom from the FastQC output.

❓ How a change in the p-value cutoff impact the result?

<details>
  <summary>Tips 👀</summary>

  > The **effective genome size**  is the portion of the genome that is mappable. There is no predefined value for E.Coli, we will use the complete genome size of the reference.<br>
  > The **shifting model** is doable when we have enough reads and enough detected peaks, or ideally paired-end data. This is not the case in our sample. We will not use a model.<br>
  > **Duplicates** are read the uniquely align to the exact same position that might come from PCR duplicates. We can ignore them with the default parameter `--keep-dup 1` in the Advanced Options.<br>
  > If you want to remove your duplicates rather than ignoring them, use the Picard tool 🪐 **MarkDuplicates** : *examine aligned records in BAM datasets to locate duplicate molecules*. 🪩 You can give it a try on the sorted BAM files.<br>
  > The **p-value (or p-adjusted) threshold** set the limit for calling a peak, the lower the threshold, the higher the number of reported peaks. This is a critical parameter to take into account.
  </details>
<br>

All onboard, we're ready to launch MACS2 run!

⚡️ Run MACS2 callpeak on FNR ChIP-seq with input as control.
* 🪐 **MACS2 callpeak** : *Call peaks from alignment results*
* Set E.Coli effective genome size to 4641652bp
* Set the MACS2 peak calling to run without model, with a 0 shift and a 200bp extsize
* Ignore duplicates
* Save a BedGraph output

❓ What does MACS2 when given `--no-model --shift 0 --extsize 200`?

❓ How many region have been reported, what is the output file format?

<details>
  <summary>Tips 👀</summary>

  > Select FNR IP as **treatment** and input as **control**
  ><div style="text-align:center"><img src="image/chap6/macs2_input.png" width="500"/></div>
  > Do not build a model
  ><div style="text-align:center"><img src="image/chap6/macs2_model.png" width="500"/></div>
  > Select an additional BedGraph output
  ><div style="text-align:center"><img src="image/chap6/macs2_output.png" width="500"/></div><br>

  > The output file is a **BED** file, we report **213** peaks enriched in FNR IP.
  ><div style="text-align:center"><img src="image/chap6/peak_info.png" width="200"/></div>
  </details>
<br>

🪩 Explore the different output files on your own, can you guess what they refer to? Expand the **Tool Standard Output** section in Galaxy output info, do you understand the log message?

🪩 Run another MACS callpeak with a different p-value or with the modelling on, what happens?

## 📍 **7. Enjoying the landscape : Data visualization**

### 🔸 **7.a Converting tracks to BigWig**

We have a **BAM** output for read alignment signal and a **BED** output for peaks. For signal browsing and large dataset, it is more efficient to work with **BigWig** and **BigBed** file format. **BigWig** fill will summarise the individual read information from **BAM** into a sum of read mapping on each position of the genome. **BigBed** is the indexed binary equivalent of **BED**. For a detail explanation of the BigWig formats, see the **UCSC Genome Browser** help on [bigWig](https://genome.ucsc.edu/goldenpath/help/bigWig.html) and [bigBed](https://genome.ucsc.edu/goldenpath/help/bigBed.html). As our BED file is small, we will not convert it to BigBed here.

⚡️ Convert both FNR IP and input BAM file to BigWig with **DeepTools**.
* 🪐  **bamCoverage** : *generates a coverage bigWig file from a given BAM or CRAM file*
* Use the **sorted BAM files** as input.
* Keep the default **50bp** granularity
* Use the **1x** coverage normalisation
* Set E.Coli effective genome size to **4641652bp**
* We are aiming to generate the track visualisation that most closely match our MACS2 analysis, so we need to also **ignore the Duplicates** and extend our single-end reads to **200bp** (see *Advanced Options*).

❓ What is the 1x normalisation? Why do we need it?

<details>
  <summary>Tips 👀</summary>

  > Select FNR IP and Input **sorted BAM** as input
  ><div style="text-align:center"><img src="image/chap7/bamcov_input.png" width="500"/></div>
  > Extend single-end reads to 200bp and ignore duplicates
  ><div style="text-align:center"><img src="image/chap7/bamcov_param.png" width="500"/></div>

  </details>
<br>

### 🔸 **7.b Assessing the quality of the peaks and FNR IP**

We now have the two required inputs for DeepTools graphs :
* A set of regions (ours MACS2 peaks as BED file)
* A set of signal tracks (our FNR and input BigWig files)

DeepTools is a very efficient tool for signal exploration, and it comes with an easy-to-use [online manual](https://deeptools.readthedocs.io/en/develop/index.html). It is usually done in two steps :
1. **pre-compute** the signal for the region of interest
2. **plot** the signal in your requested format

In this tutorial, we will do a **Profile** plot, but you can explore other possibilities there (*e.g.* Heatmaps)!

⚡️ Compute the signal matrix using **DeepTools**.
* 🪐  **computeMatrix** : *prepares data for plotting a heatmap or a profile of given regions*
* Select our MACS2 peaks as regions to plot
* Select the **BamCoverge** collection output as score files. Do not use the BAM file directl, as they are not 1x normalised it will make them harder to contrast.
* Use the **reference-point** mode centered as region's center and **1kb** window up/downstream.

❓ Why don't we use the **scale-region** mode? See [DeepTools manual](https://deeptools.readthedocs.io/en/develop/content/tools/computeMatrix.html#details) for help.

<details>
  <summary>Tips 👀</summary>

  > Select regions and signal track inputs
  ><div style="text-align:center"><img src="image/chap7/compmat_input.png" width="500"/></div>
  > Use the reference-point mode
  ><div style="text-align:center"><img src="image/chap7/compmat_refpoint.png" width="500"/></div>

  </details>
<br>

Now the visualization part!

⚡️ Generate the Profile plot using **DeepTools**.
* 🪐  **plotProfile** : *creates a profile plot for score distributions across genomic regions*
* Use our freshly-computed matrix as input
* Our tracks are 1x normalised, we can group them in the same plot

<details>
  <summary>Tips 👀</summary>

  > Select matrix as input
  ><div style="text-align:center"><img src="image/chap7/plotprof_input.png" width="500"/></div>
  > Group the two signal (FNR and input) in the same plot
  ><div style="text-align:center"><img src="image/chap7/plotprof_param.png" width=200"/></div>

  </details>
<br>

❓ Is the output graph consistent with our biological hypothesis?

Ok, the signal in our detected MACS2 peaks seems good! Yet, this plot is only focusing on our **selected regions**, how can we make sure our ChIP-seq worked? To answer that, we should look at the general distribution of ChIP-seq signal across the genome :
* If the FNR ChIP **worked**, we expect *a small number of peaks with high signal* and the rest with low signal
* If the FNR ChIP **failed**, we expect it to ressemble the non-specific Input, meaning *a uniform low signal in every region*

To assess where we are standing, we need to compute the **Fingerprint** of our data (*aka* the **Lorenz curve**), it's a quality-check plot that is also supported in DeepTools (told you this tool is great). You also have a nice explanatory section about it on the [DeepTools manual](https://deeptools.readthedocs.io/en/develop/content/tools/plotFingerprint.html?highlight=plotFingerprint#what-the-plots-tell-you).

⚡️ Generate the Lorenz curve for our FNR IP and Input data
* 🪐  **plotFingerprint** : *plots profiles of BAM files; useful for assessing ChIP signal strength*
* Select our Collection of **sorted BAM** files as input
* For time sake (and using E.Coli), reduce the **sampling rate** to 1000
* Same as the other plot, extend reads to **200bp** and **ignore duplicates**

<details>
  <summary>Tips 👀</summary>

  > Select matrix as input
  ><div style="text-align:center"><img src="image/chap7/lorenz_input.png" width="500"/></div>
  > Group the two signal (FNR and input) in the same plot
  ><div style="text-align:center"><img src="image/chap7/lorenz_param.png" width=500"/></div>
  ><div style="text-align:center"><img src="image/chap7/lorenz_out.png" width=300"/></div>
  </details>
<br>

❓ On a scale of Challenger Deep to Mt Everest, how *peaky* is our IP?

### 🔸 **7.c Browsing the signal in the Genome**

We have our data, we have visualised it's overall quality at full genome scale and at MASC2 peak levels. Now, we want to go look a specific examples of regions. To do so, we will use a **Genome Browser** that will serves as a search tool.

We will use the 🧬 **I**ntegrative **G**enome **V**iewer (**IGV**) that you should have installed locally on your machine. Another widely used browser is the [**UCSC Genome Browser**](https://genome.ucsc.edu/), which can be used online.

To get IGV up and running, we need to download locally some of our Galaxy output files :
* The E.Coli genome sequence as a **FASTA** files and its gene annotation as a **GFF3** file.
* The MACS peaks of FNR signal as **BED** file
* The normalised FNR and input signal as **BigWig** files

⚡️ Save the FASTA, BED and BigWig files locally from Galaxy history
* 🪐 Scroll down your history and simply use the *download* button, unzip the bigwig collection

<div style="text-align:center"><img src="image/chap7/save_bigwig.png" width=200"/></div>
<div style="text-align:center"><img src="image/chap7/save_bed.png" width=200"/></div>
<div style="text-align:center"><img src="image/chap7/save_fasta.png" width=200"/></div>

We're just missing the GFF3 file of gene annotations. We haven't loaded to Galaxy. To retrieve it, we must go back to where we found the E.Coli fasta sequence.

⚡️ Find the correct gene annotation file for our refrence genome
* Remember that our *E.Coli genome* is from the RefSeq database (chapter **5.a**)

<details>
  <summary>Tips 👀</summary>

  > The RefSeq ID for the E.Coli strain K-12 MG1655 is `NC_000913.3`. You can find it by searching this strain on NCBI's nucleotide database (1), filtering for **RefSeq** match only (2) and looking at the **1st** entry (3).
  ><div style="text-align:center"><img src="image/chap5/get_refseq.png" width="500"/></div>

  > On the RefSeq page, you can ask for a local data drop via the **Send to** tab.<br>
  ><div style="text-align:center"><img src="image/chap7/get_gff3.png" width="500"/></div>

  </details>
<br>

Ready to go? Let's turn IGV on and load all our files!
* 🧬 **Genomes** : *Load Genome from file...*
* Select the RefSeq fasta sequence
* 🧬 **File** : *Load from file...*
* Select the BED, GFF3 and BigWig files
* For better visualisation, **comment** (add `##` at the beginning of the line) the third line of the GFF3 in a raw text editor (*i.e.* not word) before loading.
<div style="text-align:center"><img src="image/chap7/comment_gff3.png" width=500"/></div><br>

To be able to compare our FNR IP with input, we need to look at them with the **same y-axis** and comparable **normalisation strategies**! Both tracks are 1x-normalised, so we should just make sure the y-axis range is always the same for the two tracks.
* 🧬 **Select BigWig tracks** : Right-clik and *Group Autoscale*
* 🧬 **Select BigWig tracks** : Right-clik and *Show Data Range*
* Start browsing, zoom in, zoom out, customise your tracks, ...
* Go to the following genes : *pepT*, *ycfP*

❓ What would have happen if we did not took care of group auto-scaling?

❓ Anything interesting to see ?

🪩 IGV also support BAM file as input, load the FNR sorted BAM from Galaxy to IGV

🪩 Wanna explore other genes? Sort the MACS2 BED file on Galaxy by p-adjusted values or log2 Fold Change to get the top and bottom quality peaks of our set.




IGV (GFF3, fasta, scaling, ...)

 Going further : QC = Lorenz ; Annot= ChIPseeker
