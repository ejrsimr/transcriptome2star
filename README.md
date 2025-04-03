# transcriptome2star

Convert a FASTA file of a transcriptome into a FASTA and GTF suitable for a STAR database. Particularly useful for single-cell analysis.

A very simple, very short script that takes a FASTA file as input and outputs a new FASTA file with renamed sequences and GTF suitable for creating a STAR database as required by many single cell pipelines. Most of the single cell pipelines we worked with require a STAR database indexed with a GTF, we have found this simple solution works.  

Originally written for Benham-Pyle 2021, PMID: 34475533

Usage: `python transcriptome2star.py <fasta_in> <gtf_out> <fasta_out>`

Example STAR indexing command:  
`STAR --runThreadN 8 --runMode genomeGenerate --genomeDir star_index --genomeFastaFiles fasta.ref.fa --sjdbGTFfile fasta.ref.gtf --sjdbOverhang 99 --genomeSAindexNbases 11`



## FAQS

Question: Why does this 100 line script need it's own repository?  
Answer: I have been asked for this script or output from it dozens of times over the last decade and shared it ad hoc. When it was suggested that I make a repository, I wondered why I hadn't made it years ago. While this solution is simple, it is often non-obvious. 
