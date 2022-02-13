"""
Looks for query sequence in database (bidirectional) and outputs alignment scores as well as hits

author: pauljannis
email: paul_jannis@web.de
date: 16/03/20202

Code partially taken from UMIC-seq (https://github.com/fhlab/UMIC-seq)
"""
__version__ = 0.1

import sys
from Bio import SeqIO
import multiprocessing
from skbio.alignment import StripedSmithWaterman
import argparse
import numpy as np
from tqdm import tqdm
from matplotlib import pyplot as plt


def parse_args(sysargs):
    parser = argparse.ArgumentParser(description=f"""Script for nanopore metagenomics.
                                    Author: Paul Zurek (pjz26@cam.ac.uk).
                                    Version {__version__}""")
    parser.add_argument('-v', '--version', action='version', version=__version__)
    parser.add_argument('-T', '--threads', type=int, default=0, help='Number of threads to execute in parallel. Defaults to CPU count.')
    parser.add_argument('--aln_thresh', type=int, default=200, help='Alingment score threshold for extraction.')

    parser.add_argument('-q', '--query', help='Query sequence (.ab1 or .fasta)', required=True)
    parser.add_argument('-db', '--database', help='Database (.fasta)', required=True)
    
    args = parser.parse_args(sysargs)
    return args


def main(args):
    """Find query hits in the database and save those to a new file"""
    threads = args.threads
    if threads == 0:
        threads = multiprocessing.cpu_count()

    #Loading query: Fasta or sanger sequence
    query_file = args.query
    if '.ab1' in query_file:
        query_rec = SeqIO.read(query_file, "abi-trim")
        query_file = query_file[:-4]
    elif ".fasta" in query_file:
        query_rec = SeqIO.read(query_file, "fasta")
        query_file = query_file[:-6]
    else:
        raise Exception('File format error. Please provide .ab1 or .fasta')

    #Generate SSW queries: Need to align the reverse complement too!
    ssw_queries = [StripedSmithWaterman(str(query_rec.seq), score_only=True), StripedSmithWaterman(str(query_rec.reverse_complement().seq), score_only=True)]

    #Aligning function: align all queries and return best alignment score
    def align(rec_str):
        aln_scores = []
        for query in ssw_queries:
            aln = query(rec_str)
            aln_scores.append(aln.optimal_alignment_score)
        return max(aln_scores)

    #Count records (only needed for tqdm)
    database_file = args.database
    N_seq = len([1 for line in open(database_file) if line.startswith(">")])
    print(f"Working with {N_seq} database sequences.")

    """
    #Generating scores
    print(f"Aligning ({threads} threads)...")
    pool = multiprocessing.Pool(processes=threads)
    database_gen = (str(rec.seq) for rec in SeqIO.parse(database_file, 'fasta'))  #Generator to save memory when aligning
    score_lst = list(tqdm(pool.imap(align, database_gen, chunksize=1000), total=N_seq))
    #pool imap takes a generator and returns a generator (thus turning to list here)! saves memory and enables results to be used as they come (not needed here)!
    pool.close()
    pool.join()
    """
    print("Aligning ...")
    database_gen = (str(rec.seq) for rec in SeqIO.parse(database_file, 'fasta'))
    score_lst = list(tqdm(map(align, database_gen), total=N_seq))


    #Alignment score histogram
    print("Plotting alignment scores")
    plt.figure()
    binwidth = 50
    plt.hist(score_lst, bins=np.arange(min(score_lst), max(score_lst)+binwidth, binwidth))
    plt.xlabel("SSW alignment score")
    plt.ylabel("Count")
    plt.yscale("log")
    name = query_file + "_alnscores-hist.png"
    plt.savefig(name, dpi=600, bbox_inches='tight')


    #Selecting good alignments: Print and write to file
    print("Selecting good alignments")
    aln_thresh = args.aln_thresh
    records = SeqIO.parse(database_file, "fasta")
    hits = []
    for score, rec in zip(score_lst, records):
        if score > aln_thresh:
            print("ID: %s" % rec.id)
            print("Length: %d" % len(rec.seq))
            print(rec.seq[:50] + " ... " + rec.seq[-50:])
            print("Score: %s" % str(score))
            hits.append(rec)

    print(f"Idenfitied {len(hits)} hit sequences")
    if len(hits) > 0:
        name = query_file + "_hits.fasta"
        SeqIO.write(hits, name, "fasta")


if __name__ == "__main__":
    args = parse_args(sys.argv[1:])
    main(args)
