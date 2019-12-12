Wednesday Warmup Exercise
===========================

Align E. coli reads (PacBio) to E. coli Reference
---------------------------------------------------

Let's align some PacBio reads from the E. coli K12, MG1655 strain to an E. coli K12, MG1655 reference genome. If you google for "Koren MG1655 wgs-assembler" you should find [this](http://wgs-assembler.sourceforge.net/wiki/index.php/Escherichia_coli_K12_MG1655,_using_uncorrected_PacBio_reads,_with_CA8.1) page, where three xzipped fastq files can be downloaded. Grab just one, and uncompress it. Then, grab the appropriate reference genome (as we did yesterday). Finally, take only the first ~4000 reads from the read set, and align them to the reference, to produce a SAM file, as we did yesterday. Does everything look as expected?


