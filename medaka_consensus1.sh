#!/bin/bash
echo 'Runs minimap, racon and medaka on the reads provided as the first argument!'
#First argument needs to be the sequences to assemble
hits_file=$1

#Pomoxis miniassemble generates a draft assembly with miniasm and racon
#mini_assemble -i hits_file -o draft_assm -p assm
##--> ALWAYS EMPTY! -> SELF:
racon_file='racon-draftassm_'$hits_file

minimap2 -x ava-ont $hits_file $hits_file > overlaps.paf
racon $hits_file overlaps.paf $hits_file > $racon_file
rm overlaps.paf

medaka_consensus -i $hits_file -d $racon_file -m r941_min_high_g344 -o .
rm calls_to_draft.bam
rm calls_to_draft.bam.bai
rm consensus_probs.hdf
rm $racon_file'.fai'
rm $racon_file'.mmi'

consensus_file='medaka-consensus_'$hits_file
mv consensus.fasta $consensus_file
echo ""
echo "Assembly and polishing complete. Moved consensus to "$consensus_file
