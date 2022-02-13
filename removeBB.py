"""
Remove plasmid backbone from sequence collections.
- Remove plasmid backbone in two stages.
- check for contamination

author: pauljannis
email: paul_jannis@web.de
date: 08/08/2019

version history
0.2 - added argparse
0.1 - initial commit
"""
import sys
import argparse
from Bio import SeqIO

from metalibseq.filteranalysis import alnscore_histogram, length_histogram, length_filter
from metalibseq.removal import fix_orientation, remove_backbone

__version__ = 0.2

def parse_args(sysargs):
    """Argument parser"""
    parser = argparse.ArgumentParser(description=f"""Script to remove plasmid backbone from sequence collections.
                                 Author: pauljannis (paul_jannis@web.de).
                                 Version {__version__}.""")
    parser.add_argument("-r", "--records", required=True,
                        help="Fasta file containing the target sequence collection.")
    parser.add_argument("-bb", "--backbone", required=True,
                        help="""Single record fasta file containing the plasmid backbone sequence
                        that will be removed from the sequence collection.""")
    parser.add_argument("--lengthfilter", nargs=2, type=int,
                        help="""Optional: Filter final sequences for insert length.""")
    parser.add_argument("--contaminants", help="Optional: Fasta file containing potential contaminants to be detected.")
    parser.add_argument("--keep_intermediates", action="store_true",
                        help="""Set flag to write out intermediate files.""")

    args = parser.parse_args(sysargs)
    return args


def main(args):
    """Pipeline removing plasmid backbone and checking for contamination"""
    # Read in records.
    records = list(SeqIO.parse(args.records, "fasta"))
    query_rec = SeqIO.read(args.backbone, "fasta")
    query_seq = str(query_rec.seq)

    # Get alingment score histogramm to see if backbone can be found (backbone alignment score vs all sequences histogram)
    alnscore_histogram(query_rec, records, "histo.alignment.plamid_backbone.png", bidirect=True)

    # Check for contaminants
    if args.contaminants is not None:
        contaminants = list(SeqIO.parse(args.contaminants, "fasta"))
        for cont in contaminants:
            alnscore_histogram(cont, records, f"hiso.alignment.contaminant_{cont.id}.png", bidirect=True)

    # lengthhistoram before removal
    length_histogram(records, "histo.length.pre-removal.png")

    # Flip reads to all have the same orientation as the backbone (needed for backbone removal to script to work)
    print("Fixing read orientation: please provide threshold")
    print("Note: The threshold to identify backbone alignment should be identify from the backbone alignment histogram.")
    print("      A suitable value lies just below the the main alignment scores")
    thresh = int(input("backbone threshold: "))
    flipped, discarded = fix_orientation(query_rec, records, int(thresh/10), thresh)  # used 500 and 2500 before
    if args.keep_intermediates:
        SeqIO.write(flipped, "records.corrected_orientation.fasta", "fasta")
        SeqIO.write(discarded, "records.discarded_orientation.fasta", "fasta")

    # Remove backbone and copy beginning to fix insert cuts (round 1)
    remove_thresh = 200   # remove thresh low value that just makes sure no small unspecific regions are removed
    print(f"Removing backbone above {remove_thresh} alignment threshold")
    removed = remove_backbone(query_seq, flipped, copy=True, thresh=remove_thresh)
    if args.keep_intermediates:
        SeqIO.write(removed, "records.removedBB_first-pass.fasta", "fasta")

    # Check how much BB is left after first removal (could be only half of a split backbone)
    alnscore_histogram(query_rec, removed, "histo.alignment.removedBB_first-pass.png", bidirect=False)   # bidirect not necessary, all records are oriented

    #Remove backbone round two (backbone could be split)
    removed = remove_backbone(query_seq, removed, copy=False, thresh=remove_thresh)
    if args.keep_intermediates or args.lengthfilter is None:
        SeqIO.write(removed, "records.removedBB_second-pass.fasta", "fasta")

    #Check that no BB is left finally
    alnscore_histogram(query_rec, removed, "histo.alignment.removedBB_second-pass.png", bidirect=False)

    # Length histogram of all inserts (BB should be gone now)
    length_histogram(removed, "histo.length.post-removal.png")

    # Filter by length: Good values min 400 bp, max 10000 bp
    if args.lengthfilter is not None:
        filtered, discarded = length_filter(removed, args.lengthfilter[0], args.lengthfilter[1])
        length_histogram(filtered, "histo.length.post-filter.png")
        SeqIO.write(filtered, "records.length_filter.fasta", "fasta")
        if args.keep_intermediates:
            SeqIO.write(discarded, "records.length_filter_discarded.fasta", "fasta")


if __name__ == "__main__":
    args = parse_args(sys.argv[1:])
    main(args)

