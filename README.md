# funcmeta_libseq
### Characterization of functional metagenomics libraries with long-read sequencing
Sequencing of libraries used in functional metagenomic screenings with long-read technologies, such as nanopore sequencing, can provide information of library quality and diversity as well as alleviate the need for time consuming primer walking determination of hit sequences. The scripts that will be collected here can be used to analyse such nanopore sequencing data of metagenomic libraries as, for example, used in Steffi's PhD thesis.


## Notes on library generation

Functional metagenomic libraries were fragmented to 2-5 kb inserts and cloned into expressable plasmids for functional screening. The generated libaries were sequenced on ONT's MinION at ~10-fold expected coverage, simply prepared from the plasmid stock with the RAD-003 transposon based library sequencing kit.


## Pre-requisites

Clone the repository and run the python scripts. The scripts require python 3 with Biopython, scikit-bio, numpy and matplotlib.


## Insert extraction

The plasmid library is randomly broken and sequenced. We have to remove the plasmid backbone from the reads to extract the metagenomic insert.
*all input files as .fasta*
Prepare the following files:
- library_reads.fasta: The nanopore reads of the library. Run nanofilt to beforehand for QC of appropriate approximate length and read quality.
- backbone.fasta: The plasmid backbone sequence
- Optional: potential_contaminants.fasta: A few sequences that you might want to check for as contamination in your library.

```python removeBB.py -r library_reads.fasta -bb backbone.fasta --contaminants potential_contaminants.fasta```

TODO: HELP PRINT HERE. 

TODO: ADD EXAMPLE FIGURES.

TODO: Size estimation! Diversity metrics!

The extracted inserts can be further analysed for diversity and will be used to identify and polish hits in the second step.

## Find and polish hit
Once hits are found via functional screening, the hits only need to be Sanger sequenced once at one end. We can look for matches to this hit within the nanopore insert and polish it via consensus assembly, which should be very accurate if a few reads can be collected (*see Zurek et al. [2020](https://www.nature.com/articles/s41467-020-19687-9)*) .
Input files:
- query.ab1 / query.fasta: Parial sequence of hit after screening obtained by Sanger in ab1 or fasta format
- extracted_inserts.fasta: Insert reads extracted in the previous step

```python findfamiliar.py -q query.ab1 -db extracted_inserts.fasta ```

TODO: PRINT HELP.

TODO: ADD EXAMPLE FIGURES

The alignment threshold to collect matching reads for consensus assembly should be set stringently with the --aln_thresh parameter.

The resulting insert hits are collected and output as hits.fasta. For polishing of those, racon and medaka could be used. A shell script is provided to automate this.

```medaka_consensus.sh hits.fasta```

Finally, an accurate sequence of the screened metagenomic hit is generated.

