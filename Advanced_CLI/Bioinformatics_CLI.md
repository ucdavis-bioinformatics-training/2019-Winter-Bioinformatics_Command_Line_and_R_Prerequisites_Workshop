# Project/Experiment Structure and a simple bioinformatics task

## Experiment structure

Goals:
* Entire contained within 1 folder.
* Folder structure is intuitive and self explanatory.
* Scripts are apparent
* Another user (collaborator) can intuit your analysis just from the folders
* A year from now (hopefully sooner) when you publish, you can quickly review your analysis path.
* Your analysis is Reproducible


Align E. coli reads (PacBio) to E. coli Reference
---------------------------------------------------

Let's align some PacBio reads from the E. coli K12, MG1655 strain to an E. coli K12, MG1655 reference genome.
If you google for "Koren MG1655 wgs-assembler" you should find [this](http://wgs-assembler.sourceforge.net/wiki/index.php/Escherichia_coli_K12_MG1655,_using_uncorrected_PacBio_reads,_with_CA8.1) page, where three xzipped fastq files can be downloaded. Take only the first ~4000 reads from the read set, and align them to the reference, to produce a SAM file.


Let's set up a fake sequencing project directory, and talk a bit about project philosophy..

**1\.** First, create a directory for your user and the example project in the workshop directory:

    cd
    mkdir -p /share/workshop/$USER/bioinfo_example

---

**2a\.** Next, go into that directory, create a raw data directory (we are going to call this 00-RawData) and cd into that directory. Lets then create symbolic links to the sample directories that contains the raw data (borrowed from a previous RNAseq workshop).

    cd /share/workshop/$USER/bioinfo_example
    mkdir 00-RawData
    cd 00-RawData/
    curl -L -o escherichia_coli.fastq.xz http://sourceforge.net/projects/wgs-assembler/files/wgs-assembler/wgs-8.0/datasets/escherichia_coli_k12_mg1655.m130404_014004_sidney_c100506902550000001823076808221337_s1_p0.1.fastq.xz/download
    unxz escherichia_coli.fastq.xz

This directory now contains the fastq file for one sample.

**2b\.** Let's create a sample sheet for the project and store sample names in a file called samples.txt

    ls > ../samples.txt
    cat ../samples.txt

---
**3a\.** Now, take a look at the raw data directory.

    ls /share/workshop/$USER/bioinfo_example/00-RawData


**3b\.** Lets get a better look at the files in the directory.

    ls -lah *

---

### fastq files
fastq files combine the sequence and quality scores into 1 file. Each sequence here has 4 lines (should be enforced strictly), header, sequence, historical '+', and quality.

<img src="filetypes_figures/filetypes_figure3.png" alt="figure3" width="800px"/>

CASAVA 1.8 Read IDs

@EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG
* EAS139 the unique instrument name
* 136 the run id
* FC706VJ the flowcell id
* 2 flowcell lane
* 2104 tile number within the flowcell lane
* 15343 ’x’-coordinate of the cluster within the tile
* 197393 ’y’-coordinate of the cluster within the tile
* 1 the member of a pair, 1 or 2 (paired-end or mate-pair reads only)
* Y Y if the read fails filter (read is bad), N otherwise
* 18 0 when none of the control bits are on, otherwise it is an even number
* ATCACG index sequence

#### Quality scores
Quality scores are paired 1 to 1 with sequence characters.

Each quality character has a numerical value associated with it (ASCII value). In Illumina 1.8+ you subtract 33 from the ascii value associated with the quality character to get the quality score.

<img src="filetypes_figures/filetypes_figure5.png" alt="figure5" width="800px"/>

<img src="filetypes_figures/filetypes_figure4b.png" alt="figure4b" width="400px"/>

<img src="filetypes_figures/filetypes_figure4a.png" alt="figure4a" width="400px"/>

----

**4\.** View the contents of the fastq file using the 'less' command, when gzipped used 'zless' (which is just the 'less' command for gzipped files):

    less escherichia_coli.fastq

Make sure you can identify which lines correspond to a single read and which lines are the header, sequence, and quality values. Press 'q' to exit this screen. Then, let's figure out the number of reads in this file. A simple way to do that is to count the number of lines and divide by 4 (because the record of each read uses 4 lines). In order to do this use cat to output the uncompressed file and pipe that to "wc" to count the number of lines:

    cat escherichia_coli.fastq | wc -l

Divide this number by 4 and you have the number of reads in this file. One more thing to try is to figure out the length of the reads without counting each nucleotide. First get the first 4 lines of the file (i.e. the first record):

    cat escherichia_coli.fastq  | head -2 | tail -1

Note the header lines (1st and 3rd line) and sequence and quality lines (2nd and 4th) in each 4-line fastq block. Then, copy and paste the DNA sequence line into the following command (replace [sequence] with the line):

---

**5\.** Now go back to your 'bioinfo_example' directory and create two directories called 'slurmout' and '01-HTS_Preproc':

    cd /share/workshop/$USER/bioinfo_example
    mkdir scriptout
    mkdir 01-Preproc

The results of all our scripts will output to .out and .err files into the scriptout folder. The results of our preprocessing steps will be put into the 01-Preproc directory. The next step after that will go into a "02-..." directory, etc. You can collect scripts that perform each step, and notes and metadata relevant for each step, in the directory for that step. This way anyone looking to replicate your analysis has limited places to search for the commands you used. In addition, you may want to change the permissions on your original 00-RawData directory to "read only", so that you can never corrupt your raw data. (We won't worry about this here, because we've linked in sample folders).

---

### fasta
The fasta format uses the '>' to indicate a new sequence followed by the name of the sequence on the same line. The following line(s) are the DNA sequence and may be split on multiple lines (wrapped), until the next '>' is reached. Genome and transcriptome files are most often in fasta format.

<img src="filetypes_figures/filetypes_figure2.png" alt="figure2" width="800px"/>


Qual files are so rarely used these days and so are not discussed.

---

**6\.** Lets go get the reference sequence from Ensembl, extract, and view.

    mkdir References
    cd References
    wget ftp://ftp.ensemblgenomes.org/pub/release-45/bacteria//fasta/bacteria_88_collection/escherichia_coli_str_k_12_substr_mg1655_gca_000801205/dna/Escherichia_coli_str_k_12_substr_mg1655_gca_000801205.ASM80120v1.dna.toplevel.fa.gz
    gunzip Escherichia_coli_str_k_12_substr_mg1655_gca_000801205.ASM80120v1.dna.toplevel.fa.gz
    less Escherichia_coli_str_k_12_substr_mg1655_gca_000801205.ASM80120v1.dna.toplevel.fa

---

**7\.** Lets 'preprocess' our data by just taking the first 4000 reads (16,000) lines from the fasta file placing the result it into the 01-Preproc folder. Often preprocessing will include evaluating the raw data, removing poor quality bases, reads, removing primers/adapters, PCR duplicates, etc. CREATE a script (preproc.sh) that does the following. And then run the script.


<pre class="prettyprint"><code class="language-bash" style="background-color:333333">
#!/bin/bash

cd /share/workshop/$USER/bioinfo_example
head -n 16000 00-RawData/escherichia_coli.fastq > 01-Preproc/escherichia_coli.4k.fastq
</code></pre>

    bash preproc.sh
    cd /share/workshop/$USER/bioinfo_example
    ll 01-Preproc
    less 01-Preproc/escherichia_coli.4k.fastq
    wc -l 00-RawData/escherichia_coli.fastq
    wc -l 01-Preproc/escherichia_coli.4k.fastq


**8\.** Lets align the reads to the genome using [minimap2](https://github.com/lh3/minimap2), then sort to bam format with samtools. Create a script (mapping.sh) to align the reads. View the help documentation.

<pre class="prettyprint"><code class="language-bash" style="background-color:333333">
#!/bin/bash

module load minimap2/2.7
module load samtools/1.9
minimap2 References/Escherichia_coli_str_k_12_substr_mg1655_gca_000801205.ASM80120v1.dna.toplevel.fa 01-Preproc/escherichia_coli.4k.fastq -a | samtools sort -o 02-Mapping/escherichia_coli.4k.bam
samtools index 02-Mapping/escherichia_coli.4k.bam
</code></pre>

    mkdir 02-Mapping
    bash mapping.sh
    ll 02-Mapping
    samtools view 02-Mapping/escherichia_coli.4k.bam

For those who are interested, here's how you'd create a BAM file that could be viewed in [IGV](https://software.broadinstitute.org/software/igv/download)

**9\.** Lets get some stats on the mapping results, using samtools.

    module load samtools/1.9
    samtools flagstat 02-Mapping/escherichia_coli.4k.bam > 02-Mapping/escherichia_coli.4k.flagstat
    samtools depth 02-Mapping/escherichia_coli.4k.bam > 02-Mapping/escherichia_coli.4k.depth
    samtools stats 02-Mapping/escherichia_coli.4k.bam > 02-Mapping/escherichia_coli.4k.stats

View the output files to see the outcome of mapping results.

**10\.** Transfer the depth file to your computer and analyze in **R**, perform some summary statistics, plot the distribution of depth over the genome, plot the coverage depth over the whole genome
