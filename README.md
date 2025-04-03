# transcriptome2star

Convert a FASTA file of a transcriptome into a FASTA and GTF suitable for a STAR database. Particularly useful for single-cell analysis when no genome assembly is available.  

A very simple, very short script that takes a FASTA file as input and outputs a new FASTA file with a ".ref" suffix added to seach sequence name and GTF with a single exon transcript covering the entire length of each sequence, retaining the original transcript name. This output is suitable for creating a STAR database as required by many single cell and RNAseq analysis pipelines.  

Originally written for Benham-Pyle 2021, PMID: 34475533  

Usage: `python transcriptome2star.py <fasta_in> <gtf_out> <fasta_out>`  

Example STAR indexing command:  
`STAR --runThreadN 8 --runMode genomeGenerate --genomeDir star_index --genomeFastaFiles fasta.ref.fa --sjdbGTFfile fasta.ref.gtf --sjdbOverhang 99 --genomeSAindexNbases 11`  



## FAQS  

Question: Why does this 100 line script need it's own repository?   

Answer: I have been asked for this script or output from it dozens of times over the last decade and shared it ad hoc. When it was suggested that I make a repository, I wondered why I hadn't made it years ago. While this solution is simple, it is often non-obvious.  
