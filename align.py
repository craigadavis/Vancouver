"""
This is an example script that takes a test file containing some
sequences, aligns everything against the first sequence, and 
clips out insertions relative to this reference.  It generates
a FASTA-formatted output file that will contain these aligned
sequences.
"""

# load modules
import HyPhy
import os
import hyphyAlign
from seqUtils import convert_fasta

# start an instance of HyPhy to communicate with
hyphy = HyPhy._THyPhy (os.getcwd(), 1)


# set HyPhy for nucleotide alignment
hyphyAlign.change_settings(hyphy,
                           alphabet=hyphyAlign.nucAlphabet,
                           scoreMatrix=hyphyAlign.nucScoreMatrix,
                           gapOpen=20, gapOpen2=20,
                           gapExtend=10, gapExtend2=10,
                           noTerminalPenalty=1)

# open the FASTA file and convert into a Python object
# handle = open('test full1302.fasta', 'rU')
handle = open('Vancouver_Bref_1302.fasta', 'rU')
fasta = convert_fasta(handle)
handle.close()


# use B sequence as reference
nameref, refseq = fasta[0]

# prepare file to write results
# outfile = open('align-out.fa', 'w')
outfile = open('Vancouver_Bref_1302_aligned-out.fa', 'w')

for header, sequence in fasta[1:]:
    # remove any gap characters from the original sequences
    new_seq = sequence.replace('-', '')

    # align this sequence against the reference
    aquery, aref, ascore = hyphyAlign.pair_align(hyphy, refseq, new_seq)

    # see where in the reference the sequence aligned
    left, right = hyphyAlign.get_boundaries(aref)
    
    # ignore insertions relative to reference
    new_seq2 = ''  # initialize a new empty string
    
    # iterate over every base in the aligned sequence
    for i, bp in enumerate(aquery):
        
        # examine the corresponding base in the aligned reference
        if aref[i] == '-':
            # if the reference contains a gap then this is an insertion
            continue
        
        # otherwise concatenate this base to a new sequence
        new_seq2 += aquery[i]

    # write the result to our file
    # outfile.write('>%s\n%s\n' % (header, new_))
    outfile.write('>%s\n%s\n' % (header, new_seq2))

outfile.close()
