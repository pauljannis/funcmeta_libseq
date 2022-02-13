#v0.2
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#FUNCTIONS TO REMOVE THE PLASMID BACKBONE FORM NANOPORE READS
#USED FOR METAGENOMIC LIBRARY SEQUENCING WITH STEFFI AND MIRIAM
#PZ 08/08/2019

## FLIP NANOPORE READS TO BE IN BACKBONE ORIENTATION 
# fix_orientation(query_rec, targets, low_tresh, high_thresh)
#query_rec: Takes the plasmid backbone sequence, linearized at the integration cut site, as a biopython record
#targets: Takes the reads as a python _list_ of biopython records
#low_thresh, high_thresh: One orientation must have a lower alnignment score than low_thresh and the other must have a higher score than high_thresh
#                         to flip the sequence accordingly. Otherwise discarded (no good alignment at all or both orientations are aligned well -> ambiguous)
#RETURNS: 2 lists: list of records with correct orientation, list of records discarded due to ambiguous alignment
#Suggested: Test alignment scores first with alnscore_histogram

## REMOVE BACKBONE FROM READS
# remove_backbone(query_seq, targets, copy, thresh, verbose=False)
#query_seq: string sequence that is to be removed (e.g. plasmid sequence, after orienting the reads)
#targets: list of biopython records
#copy: Flag to fix broken inserts: Should be true on first run, false on consecutive runs.
#      -> True: Copies everything of the target sequence before the alignment start to the end of the target sequence,
#               thus fixing the read if it is split in the insert
#      -> False: On the second run this should be false and will then remove the "second half" of the query if the query is split
# Thresh: Alignment score threshold for removal
#RETURNS: List of targets with removed query_seq


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
from matplotlib import pyplot as plt
import numpy as np
from Bio.Seq import Seq
from skbio.alignment import StripedSmithWaterman
import multiprocessing

## HELPER FUNCTIONS

#Alignment sub-function, returning full alignment info to the multiprocessing function (SSW_alignment object can't be pickled)
def calc_score(rec):
    global SSWquery
    aln = SSWquery(rec)
    return aln.optimal_alignment_score, aln.query_begin, aln.query_end, aln.target_begin, aln.target_end_optimal

#Calculates alignment values on multiple cores
#One sequence (string, query_seq) vs a list of biopython records
#Returns: [score, query_begin, query_end, target_begin, target_end] list of lists (or score only list of int if passed)
def calc_aln_lst(query_seq, records, score_only=False):  #Uses calc_score
    rec_lst = [str(rec.seq) for rec in records]
    global SSWquery
    SSWquery = StripedSmithWaterman(query_seq)
    threads = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(processes=threads)
    aln_values = pool.map(calc_score, rec_lst)
    pool.close()
    if score_only:
        aln_values = [val[0] for val in aln_values]
    return aln_values


## MAIN FUNCTIONS

#Fixes read orientation based on query sequence alignment
#Targets needs to be a list not a generator ~~ list(SeqIO.parse())
def fix_orientation(query_rec, targets, low_tresh, high_thresh):  #Uses calc_aln_lst
    fwd_seq = str(query_rec.seq)
    rev_seq = str(query_rec.reverse_complement().seq)

    print("Fixing read orientation (%d cores)" % multiprocessing.cpu_count())
    #Calculate all alingment scores (scores only)
    alnF_lst = calc_aln_lst(fwd_seq, targets, True)
    alnR_lst = calc_aln_lst(rev_seq, targets, True)

    flipped_tar, disc_tar = [], []
    for i in range(len(targets)):
        scoreF = alnF_lst[i]
        scoreR = alnR_lst[i]
        #print("Fwd: %d\t Rev: %d" % (scoreF, scoreR))
        if scoreF > high_thresh and scoreR < low_tresh:
            flipped_tar.append(targets[i])
        elif scoreF < low_tresh and scoreR > high_thresh:
            flip_seq = targets[i]
            flip_seq.seq = targets[i].seq.reverse_complement()
            flipped_tar.append(flip_seq)
        else:
            #print("Fwd: %d\t Rev: %d" % (scoreF, scoreR))
            disc_tar.append(targets[i])

    print("Flipping completed")
    print("Discarded %d sequences because of ambiguous alignment" % len(disc_tar))
    return flipped_tar, disc_tar
        

#Removes query_seq from all targets if alignment score is > thresh. Copies beginning of target to its end if copy = true.
def remove_backbone(query_seq, targets, copy, thresh, verbose=False):  #Uses calc_aln_lst
    print("Removing backbone (%d cores)" % multiprocessing.cpu_count())
    aln_lst = calc_aln_lst(query_seq, targets, False)

    removed = []
    for i in range(len(targets)):
        tar_seq = str(targets[i].seq)
        score, q_begin, q_end, t_begin, t_end = aln_lst[i]
        if verbose:
            print("Score: %d Readlength: %d" % (score, len(tar_seq)))
            print("QUERY: Start: %d End: %d" % (q_begin, q_end))
            print("TARGET: Start: %d End: %d" % (t_begin, t_end))
            print("")
        if score > thresh:
            if copy:
                new_tar_seq = tar_seq[t_end:] + tar_seq[:t_begin]
            else:
                new_tar_seq = tar_seq[:t_begin] + tar_seq[t_end:]
            new_tar = targets[i]
            new_tar.seq = Seq(new_tar_seq)
            removed.append(new_tar)
        else:
            removed.append(targets[i])
    print("Backbone removed")
    return removed

