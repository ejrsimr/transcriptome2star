# transcriptome2star

Convert a FASTA file of a transcriptome into a FASTA and GTF suitable for a STAR database. Particularly useful for single-cell analysis.

A very simple, very short script that takes a FASTA file as input and outputs a new FASTA file with renamed sequences and GTF suitable for creating a STAR database as required by many single cell pipelines.

Originally written for Benham-Pyle 2021, PMID: 34475533

Usage: `python transcriptome2star.py <fasta_in> <gtf_out> <fasta_out>`

Example STAR indexing command:  
`STAR --runThreadN 8 --runMode genomeGenerate --genomeDir star_index --genomeFastaFiles fasta.ref.fa --sjdbGTFfile fasta.ref.gtf --sjdbOverhang 99 --genomeSAindexNbases 11`
