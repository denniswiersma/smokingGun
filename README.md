# Smoking Gun
The goal of this assignment is to predict whether a Single Nucleotide Polymorphism (SNP) is deleterious for a given gene sequence.

# Usage
`galectin3/homologeneFasta.fasta -p 28 G -s galectin3/galectin3_DNA.fasta`

# Assignment execution
## Finding a protein
Since I'm currently doing a project related to heart disease I decided to look for a protein linked to this condition.
After a quick google (or rather duckduckgo) search I found just what I was looking for on [this](https://www.nih.gov/news-events/news-releases/protein-linked-increased-risk-heart-failure-death-older-adults) page.
It describes how the galectin-3 protein can identify people at higher risk of heart failure.
So, naturally I went ahead and did a protein search on NCBI to get to [galectin-3's page](https://www.ncbi.nlm.nih.gov/protein/NP_001183972.1). 

Next stop: the protein's [HomoloGene page](https://www.ncbi.nlm.nih.gov/homologene?LinkName=protein_homologene&from_uid=308081799) where we can find a [Multiple Sequence Alignment](https://www.ncbi.nlm.nih.gov/homologene?cmd=Retrieve&dopt=MultipleAlignment&list_uids=37608) (MSA).
Unfortunately we run into a small issue at this step. When we try to download this MSA we simply get a non-aligned fasta file (which can be found in `galectin3/HomologeneFasta.fasta`). 
To get an alignment out of this fasta file I fed it into [Clustal Omega](http://www.clustal.org/omega/) using default settings, resulting in the `galectin3/msa.fasta`.
