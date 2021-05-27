#!/usr/bin/env python

import sys
from Bio import SeqIO


def get_ids(fname):
    """ """
    return( [ L.strip() for L in open(fname).readlines()] )

def get_recs(seqs, id_list, out_handle = sys.stdout):
    to_write = (r for r in recs if r.id in id_list)
    n = SeqIO.write(to_write, out_handle, "fasta")
    return(n)

if __name__ == "__main__":
    try:
        seq_file, id_file =  sys.argv[1:]
    except ValueError:
        print("Usage get_fasta_subset.py [seq.fa] [ids.list]")
        sys.exit(1)
    recs = SeqIO.parse(seq_file, "fasta")
    ids = get_ids(id_file)
    n = get_recs(recs, ids)
    sys.stderr.write("Wrote {} records\n".format(n))
    sys.exit(0)
