# funcmeta_libseq
### Characterization of functional metagenomics libraries with long-read sequencing
Sequencing of libraries used in functional metagenomic screenings with long-read technologies, such as nanopore sequencing, can provide information of library quality and diversity as well as alleviate the need for time consuming primer walking determination of hit sequences. The scripts that will be collected here can be used to analyse such nanopore sequencing data of metagenomic libraries as, for example, used in Steffi's PhD thesis.


## Notes on library generation

Functional metagenomic libraries were fragmented to 2-5 kb inserts and cloned into expressable plasmids for functional screening. The generated libaries were sequenced on ONT's MinION at ~10-fold expected coverage, simply prepared from the plasmid stock with the transposase-based rapid library preparation kit SQK-RAD003.


## Pre-requisites

Clone the repository and run the python scripts. The scripts require python 3 with Biopython, scikit-bio, numpy and matplotlib.


## Insert extraction

The plasmid library is randomly broken and sequenced. We have to remove the plasmid backbone from the reads to extract the metagenomic insert.

Prepare the following files:
- library_reads.fasta: The nanopore reads of the library. Run nanofilt to beforehand for QC of appropriate approximate length and read quality.
- backbone.fasta: The plasmid backbone sequence
- Optional: potential_contaminants.fasta: A few sequences that you might want to check for as contamination in your library.
Furthermore, a size selection of the final inserts can be performet with the --lengthfilter argument, e.g. keeping only inserts between 400 and 10000 bp in size. If intermediate files are to be saved to disk, the --keep_intermediate flag can be set.

The user will be asked to provide a backbone alignment threshold, which can be approximated from the backbone alignment histogram (see figure below). The value chosen should lie just below the majority of alignment scores, in this example 2500 works well. 

Example usage:
```python removeBB.py -r library_reads.fasta -bb backbone.fasta --lengthfilter 400 10000 --contaminants potential_contaminants.fasta```


```
usage: removeBB.py [-h] -r RECORDS -bb BACKBONE
                   [--lengthfilter LENGTHFILTER LENGTHFILTER]
                   [--contaminants CONTAMINANTS] [--keep_intermediates]

Script to remove plasmid backbone from sequence collections. Author:
pauljannis (paul_jannis@web.de). Version 0.2.

optional arguments:
  -h, --help            show this help message and exit
  -r RECORDS, --records RECORDS
                        Fasta file containing the target sequence collection.
  -bb BACKBONE, --backbone BACKBONE
                        Single record fasta file containing the plasmid
                        backbone sequence that will be removed from the
                        sequence collection.
  --lengthfilter LENGTHFILTER LENGTHFILTER
                        Optional: Filter final sequences for insert length.
  --contaminants CONTAMINANTS
                        Optional: Fasta file containing potential contaminants
                        to be detected.
  --keep_intermediates  Set flag to write out intermediate files.
```

Example output:
```
Figure 1: Histogram of backbone alignment scores
Figure 2: Histogram of insert lengths
```
1: <img src="/example/alignmenthistogram_example.png" width="384" height="288">
2: <img src="/example/lengthhistogram_example.png" width="384" height="288">


The extracted inserts can be further analysed for diversity and will be used to identify and polish hits.


## Diversity analysis

The resulting insert sequences can be analysed for library diversity.

Example usage:
```analyse_lib.py -i extracted_inserts.fasta```

The script will approximate the over-sequencing rate and provide estimates on the library diversity. 


## Find and polish hit

Once hits are found via functional screening, the hits only need to be Sanger sequenced once at one end. We can look for matches to this hit within the nanopore insert and polish it via consensus assembly, which should be very accurate if a few reads can be collected (*see Zurek et al. [2020](https://www.nature.com/articles/s41467-020-19687-9)*) .

Input files:
- query.ab1 / query.fasta: Parial sequence of hit after screening obtained by Sanger in ab1 or fasta format
- extracted_inserts.fasta: Insert reads extracted in the previous step

```python findfamiliar.py -q query.ab1 -db extracted_inserts.fasta ```


```
usage: findfamiliar.py [-h] [-v] [-T THREADS] [--aln_thresh ALN_THRESH] -q
                       QUERY -db DATABASE

Script for nanopore metagenomics. Author: Paul Zurek (pjz26@cam.ac.uk).
Version 0.1

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
  -T THREADS, --threads THREADS
                        Number of threads to execute in parallel. Defaults to
                        CPU count.
  --aln_thresh ALN_THRESH
                        Alingment score threshold for extraction.
  -q QUERY, --query QUERY
                        Query sequence (.ab1 or .fasta)
  -db DATABASE, --database DATABASE
                        Database (.fasta)
```

The alignment threshold to collect matching reads for consensus assembly should be set stringently with the --aln_thresh parameter and a suitable value can be identified with the alignment score histogram. Almost all of the sequences in the database are expected to not match the query and have a very low alignment score, while a few hits should have a high alignment score. The threshold should be stringently set in between those two cases (for example 1000 for the data below).

Example alignment score histogram:
<img src="/example/queryalignment_histogram_example.png" width="341" height="259">


The resulting insert hits are collected and output as hits.fasta. For polishing of those, racon and medaka could be used. A shell script is provided to automate this.

```medaka_consensus.sh hits.fasta```

Finally, an accurate sequence of the screened metagenomic hit is generated.
