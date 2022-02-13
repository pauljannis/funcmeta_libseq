#v0.2
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#FUNCTIONS TO ANALYZE AND FILTER READS
#PZ 08/08/2019
#Added sequence reversing (bidirectionality) to alignments in similarity and alignment score histogram

## SEQUENCE LENGTH HISTOGRAM
# length_histogram(records, name)
#recods: list of biopython records
#name: exact filename of figure (example.png)

## SEQUENCE LENGTH FILTER
# length_filter(records, low, high)
#recods: list of biopython records
#low, high: low and high threshold for sequence length in bp
#           -> Sequences with length L will be kept if low < L < high
#Returns: Two lists: list of records passing length criterion, list of records failing length criterion

## SIMILARITY HISTOGRAM
#Useful for clustering
#Samples (size) sequences from (records) and aligns them against all (records)
# similarity_histogram(records, size, name, verbose=False)
#Records: list of biopython records
#Size: Number of sequences to be sampled from records and aligned vs records
#name: name of plots (without ".png")

## ALIGNMENT SCORE HISTOGRAM OF ONE SEQUENCE VS ALL READS
#Useful to search for a sequence in a collection
# alnscore_histogram(query_rec, targets, name)
#query_rec: Record to be aligned to all targets
#targets: list of biopython records
#Name: exact name of figure file (example.png)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
from matplotlib import pyplot as plt
import numpy as np
from skbio.alignment import StripedSmithWaterman
import multiprocessing as mlt
import itertools

threads = mlt.cpu_count()
#threads = 6


## HELPER FUNCTIONS

#Sub-function for multiprocessing alignments. Calculates alignment and returns different alignment metrics:
# [SSW score, % identity in alignment, % identity over shorter sequence, % coverage, % relative SSW score over shorter sequence]
def all_scores(query_nr, record, bidirect):   
    global query_lst
    global rev_query_lst
    aln = query_lst[query_nr](record)
    if bidirect:  #Do reverse alignment and decide based on alignment score
        alnR = rev_query_lst[query_nr](record)
        if alnR.optimal_alignment_score > aln.optimal_alignment_score:
            aln = alnR
    #Alignment score
    score = aln.optimal_alignment_score
    #Calculate query coverage
    query_length = len(aln.query_sequence)
    aligned_query_length = aln.query_end - aln.query_begin + 1
    coverage = aligned_query_length / query_length * 100
    #Calculate %identity in alignment
    aln_query = aln.aligned_query_sequence
    aln_target = aln.aligned_target_sequence
    aln_length = len(aln_query)
    same_aa = sum(e1 == e2 for e1, e2 in zip(aln_query, aln_target))
    ALNident = same_aa / aln_length * 100
    #Calculate global %identity (%matches relative to full length of shorter sequence (like CD-HIT))
    target_length = len(aln.target_sequence)
    if query_length < target_length:
        glo_len = query_length
    else:
        glo_len = target_length
    GLOident = same_aa / glo_len * 100
    #Relative alignment score
    max_score = glo_len * 2   #match_score is set to default 2
    rel_score = score / max_score * 100
    return [score, ALNident, GLOident, coverage, rel_score]


#Multiprocessing few vs all alignments
#Calculates (full) alignment scores of a list of records vs a list of records
#RETURNS a list of [score, %ALNident, %GLOident, %covery, %score] for each record for each records
def multi_aln_scores(records, queries, bidirect=True, remove_high=False, print_thresh=0):  #Uses all_scores
    rec_lst = [str(rec.seq) for rec in records]
    N_queries = len(queries)
    global query_lst
    global rev_query_lst
    if bidirect:  #Make the reverse complement query list
        rev_query_lst = [StripedSmithWaterman(str(q.reverse_complement().seq)) for q in queries]
    query_lst = [StripedSmithWaterman(str(q.seq)) for q in queries]
    pool = mlt.Pool(processes=threads)
    print('Started processing with %d threads' % threads)
    count = 0
    score_lst_lst = []
    for i in range(N_queries):
        score_lst = pool.starmap(all_scores, zip(itertools.repeat(i), rec_lst, itertools.repeat(bidirect)))
        score_score_lst = [v[0] for v in score_lst]
        if print_thresh > 0:  #Print record if alignment score is over threshold
            for j in range(len(score_score_lst)):
                if score_score_lst[j] > print_thresh:
                    print("ID: %s" % records[j].id)
                    print("Length: %d" % len(records[j].seq))
                    print(records[j].seq[:100] + " ... " + records[j].seq[-100:])
                    print("Scores: %s" % str(score_lst[j]))
        if remove_high:  #Remove the sequence with the highest score (e.g. if the query is in the records itself)
            maxind = score_score_lst.index(max(score_score_lst))
            del score_lst[maxind]  #Remove the self-alignment score
        score_lst_lst.append(score_lst)
        count += 1    
        if N_queries > 1: print("iteration %d of %d" % (count, N_queries), end='\r')
    pool.close()
    print('\nfinished generating the alignment scores')
    if N_queries == 1:
        return score_lst_lst[0]
    else:
        return score_lst_lst





## MAIN FUNCTIONS

#Calculates sequence length histogram
def length_histogram(records, name):
    lengths = [len(rec.seq) for rec in records]
    plt.figure()
    binwidth = 50
    plt.hist(lengths, bins=np.arange(min(lengths), max(lengths)+binwidth, binwidth))
    plt.xlabel("Sequence length")
    plt.ylabel("Count")
    plt.savefig(name, dpi=600)


#Filters by sequences by length
def length_filter(records, low, high):
    records_filt, records_disc = [], []
    for rec in records:
        if low < len(rec) < high:
            records_filt.append(rec)
        else:
            records_disc.append(rec)
    #records_filt = [rec for rec in records if low < len(rec.seq) < high]
    print("%d sequences of %d retained after length filtering (%.2f%%)" % (len(records_filt), len(records), (len(records_filt)/len(records)*100)))
    return records_filt, records_disc


#Calculate similarity histogram (one histogram for each metric)
#Sample x sequences from fasta for calculation and sums the result into one histogram
def sum_similarity_histogram(records, size, name, bidirect):  #Uses multi_aln_scores
    N_rec = len(records)
    #sample some sequneces
    sample = np.random.choice(records, size=size, replace=False)
    #Calculate aln scores
    scores = multi_aln_scores(records, sample, bidirect=bidirect, remove_high=True)
    scores_f = [x for y in scores for x in y]  #Flatten scores list
    print("Generating histogram of %d alignment scores from %d sequences vs %d data" % (len(scores_f), size, N_rec))
    #Plotting
    metric = ["SSW alignment score", "Alignment identity [%]", "Global Identity [%]", "Query coverage [%]", "Relative SSW alignment score [%]"]
    binwidths = [50, 1, 1, 1, 1]
    for i in range(len(metric)):
        values = [sc[i] for sc in scores_f]
        plt.figure()
        binwidth = binwidths[i]
        plt.hist(values, bins=np.arange(min(values), max(values)+binwidth, binwidth))
        plt.xlabel(metric[i])
        plt.ylabel("Count")
        plt.yscale("log")
        plt.title("Summed similarity histogram")
        plt.savefig(name+str(i)+".png", dpi=600)


#Calculate similarity histogram (one histogram for each metric)
#Generates 25 one sequence vs all plot
def single_similiarity_histogram(records, name, bidirect):  #Uses multi_aln_scores
    N_rec = len(records)
    #sample 25 sequneces
    sample = np.random.choice(records, 25, replace=False)
    #Calculate aln scores
    scores = multi_aln_scores(records, sample, bidirect=bidirect, remove_high=True)
    #Plotting
    metric = ["SSW alignment score", "Alignment identity [%]", "Global Identity [%]", "Query coverage [%]", "Relative SSW alignment score [%]"]
    binwidths = [50, 1, 1, 1, 1]
    for met_i in range(len(metric)):
        f, axarr = plt.subplots(5, 5, sharex=True, sharey=True, figsize=(12,8))
        f.subplots_adjust(hspace = 0.2, wspace = 0.1)
        f.suptitle("Individual similarity histograms", fontsize=12)
        f.text(0.5, 0.04, metric[met_i], ha='center')
        f.text(0.04, 0.5, 'Count', va='center', rotation='vertical')
        i = 0
        for a in range(5):
            for b in range(5):
                values = [lst[met_i] for lst in scores[i]]   #Get the right metric for each dataset
                axarr[a,b].hist(values, bins=np.arange(0, max(values)+binwidths[met_i], binwidths[met_i]))
                #axarr[a,b].set_yticklabels([])
                #axarr[a,b].set_xticklabels([])
                axarr[a,b].set_yscale('log')
                i += 1
        plt.savefig(name+str(met_i)+".png", dpi=600)


#Calculates all alignments of one sequence to a list of records and plots the scores as a histogram
def alnscore_histogram(query_rec, records, name, bidirect=True, print_thresh=0):  #uses calc_aln_lst
    print("Calculating alignment scores for histogram (%d cores)" % mlt.cpu_count())
    ALL_scores = multi_aln_scores(records, [query_rec], bidirect=bidirect, remove_high=False, print_thresh=print_thresh)
    aln_scores = [sc[0] for sc in ALL_scores]
    print("Histogram of alignment scores for %d sequences generated" % len(records))
    plt.figure()
    binwidth = 50
    plt.hist(aln_scores, bins=np.arange(min(aln_scores), max(aln_scores)+binwidth, binwidth))
    plt.xlabel("SSW alignment score")
    plt.ylabel("Count")
    plt.savefig(name, dpi=600)
    return aln_scores

