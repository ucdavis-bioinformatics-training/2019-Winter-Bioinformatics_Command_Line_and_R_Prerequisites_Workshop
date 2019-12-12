Wednesday Warmup Exercise
===========================

Align E. coli reads (PacBio) to E. coli Reference
---------------------------------------------------

Let's align some PacBio reads from the E. coli K12, MG1655 strain to an E. coli K12, MG1655 reference genome. If you google for "Koren MG1655 wgs-assembler" you should find [this](http://wgs-assembler.sourceforge.net/wiki/index.php/Escherichia_coli_K12_MG1655,_using_uncorrected_PacBio_reads,_with_CA8.1) page, where three xzipped fastq files can be downloaded. Grab just one, and uncompress it. Then, grab the appropriate reference genome (as we did yesterday). Finally, take only the first ~4000 reads from the read set, and align them to the reference, to produce a SAM file, as we did yesterday. Does everything look as expected?

SOLUTION
----------

From the provided link, download the first fastq.xz file ... the direct links aren't very direct, but the first of the three provided curl commands works nicely:

    curl -L -o escherichia_coli_k12_mg1655.m130404_014004_sidney_c100506902550000001823076808221337_s1_p0.1.fastq.xz http://sourceforge.net/projects/wgs-assembler/files/wgs-assembler/wgs-8.0/datasets/escherichia_coli_k12_mg1655.m130404_014004_sidney_c100506902550000001823076808221337_s1_p0.1.fastq.xz/download

Then, a little googling around will tell you the command to uncompress a .xz file:

    unxz escherichia_coli_k12_mg1655.m130404_014004_sidney_c100506902550000001823076808221337_s1_p0.1.fastq.xz

Some more googling should get you to the Illumina iGenomes [site](https://support.illumina.com/sequencing/sequencing_software/igenome.html):

    https://support.illumina.com/sequencing/sequencing_software/igenome.html

And we'll download the only MG1655 archive:

    wget ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/Escherichia_coli_K_12_MG1655/NCBI/2001-10-15/Escherichia_coli_K_12_MG1655_NCBI_2001-10-15.tar.gz

And un-archive it using the 'tar' command:

    tar -xzvf Escherichia_coli_K_12_MG1655_NCBI_2001-10-15.tar.gz

Now we want to align the first 4,000 reads to the genome. Since each sequence takes up four lines, we want the first 16k lines:

    head -n 16000 escherichia_coli_k12_mg1655.m130404_014004_sidney_c100506902550000001823076808221337_s1_p0.1.fastq > ecoli.4k.fastq

Then, we *could* use the pre-indexed files included in the archive, but for simplicity's sake, let's do it ourselves. We'll link in the genome, load the BWA module, then index:

    ln -s Escherichia_coli_K_12_MG1655/NCBI/2001-10-15/Sequence/WholeGenomeFasta/genome.fa .
    module load bwa/0.7.16a
    bwa index genome.fa

Finally, we'll align using BWA's MEM algorithm:

    bwa mem genome.fa ecoli.4k.fastq > ecoli.PB.MEM.sam

For those who are interested, here's how you'd create a BAM file that could be viewed in [IGV](https://software.broadinstitute.org/software/igv/download):

    module load samtools/1.8
    samtools view -uhS ecoli.PB.MEM.sam | samtools sort -o ecoli.PB.MEM.bam -
    samtools index ecoli.PB.MEM.bam
    samtools flagstat ecoli.PB.MEM.bam  # for some summary stats

Note that the trailing '-' in the second line of the preceding code block is a substitute for a file name, when the command is receiving data from standard in (STDIN), but the command format requires a file name.





