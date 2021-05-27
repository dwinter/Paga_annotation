import sys
import re
from Bio import SeqIO

rxlr_motif = re.compile("R\SLR")

tricky_motif = re.compile("R\SLR\S{1,40}[ED][ED][KR]")
crinkler_motif = re.compile("^\S{30,70}LFLA[RK]")

def find_motif(motif, seq, start=29, end=60):
    m = motif.search(seq[start:end])
    if m:
        return( [True, m.start() + start +1 ] )
    return( [False, None] )

def find_rxlr(motif, seq, cleavage, sp_cutoff= 40, end=100):
   if cleavage > sp_cutoff:
       return( [False, None])
   return(find_motif(motif, seq, start=cleavage, end=cleavage+end))


def get_cleaveage_site(row):
    res = row.split()
    res[1] = int(res[1].split("-")[0])
    return(res)

def write_list(gene_id_list, out_file):
    with open(out_file, "w") as out:
        for g in gene_id_list:
            out.write(g+"\n")
    


#start by doing the simples seaches, but store the sequences as dict for later
#more compelex seearch


recs = SeqIO.parse("data/secreted_proteins.faa", "fasta")
seq_dict = SeqIO.to_dict(recs)



#set up emptpy lists for genes that meet criteria
rxlrs = []
crinkler = []
nlp = []

for r in seq_dict.values():
    present, start = find_motif(rxlr_motif, str(r.seq))
    if present:
        rxlrs.append( [r, start] )
    crinkler_match = crinkler_motif.search(str(r.seq))
    if crinkler_match:
        crinkler.append(r.id)
    if r.seq.find("GHRHDWE") >= 0:
        nlp.append(r.id)


#write everything out. rxlrs are a bit mroe complex because I was orginally
# storing the sequences for futher analysis. Could be streamlined
write_list(crinkler, "results/crinkler.list")
sys.stdout.write("wrote {} putative crinklers\n".format(len(crinkler)))
write_list(nlp, "results/NLP.list\n")
sys.stdout.write("wrote {} putative NLPs\n".format(len(nlp)))

with open("results/simple_RE.tsv", "w") as out:
    for rec, s in rxlrs:
        out.write("{}\t{}\n".format(rec.id, s))

sys.stdout.write("wrote {} rxlrs for the simple RE approach\n".format(len(rxlrs)))


# THen the more complex one, needing to take cleavage sites into account

to_keep = []
sp = open("data/cleavage_sites.tsv")
for line in sp:                     
    gene, site  = get_cleaveage_site(line)
    rxlr, site = find_rxlr(tricky_motif, str(seq_dict[gene].seq),site)
    if rxlr:
        to_keep.append([gene, site])

with open("results/complex_RE.tsv", "w") as out:
    for g, s in to_keep:
        out.write("{}\t{}\n".format(g, s))
sys.stdout.write("wrote {} rxlr for the complex RE approach\n".format(len(to_keep)))


