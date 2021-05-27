# annotation of rxlr effecotrs from a set of proteins

# before this analysis

We need to run augustus (or whatever) to get our protein coding genes and
singalP to identify. The files should be in the data directory:

```sh
data/proteins.faa
data/singal_p_res.tsv
```

The data directory also includes the published hmm profile.

# rxlr hunting

We are going to use three approaches to ID putative rxlrs, two based on simle
sequence matching and another based on an hmm. 

## hmm

Start by getting a file with just the putative secreted proteins.

```
awk '{ if ($2 != "OTHER") { print $1 } } ' data/signal_p_res.tsv  > data/secreted.list
scripts/get_fasta_subset.py ../proteins.faa data/secreted.list  > data/secreted_proteins.faa
```
Then you just need to run hmmscan to search for the motif. The basis of this hmm
is from doi:10.1105/tpc.107.051037.

```
hmmscan --tblout results/secrted_rxlr_hmm.tsv --cpu 10 hmm/rxlr.hmm data/proteins.faa > hmm.out
```

And the number can of hits can be counted in the results

```
grep -v "#" results/secrted_rxlr_hmm.tsv | wc
```

NOTE a number of the hits discovred by this approach seem to not actually have a
RxLR domain (but score well elsewhere in the alignment. Should these be manually
checked?

## regex

To do the regex steps we first need to record the cleavage sites of the
predicted secrted proteins (as the RxLRs must fall within 100bp of this site for
one approach). We can extract that we awk easily enough 

```
awk '{ if ($2 != "OTHER") { print $1,$7 } } ' data/signal_p_res.tsv  > data/cleavage_sites.tsv
```

The actually logic of the calling the regexs is handed off to a python script:

## Compare them

Finally, we can use a little bit of R to see how much the different call-sets
overlap

```r
library(UpSetR)
simple_RE <- read.table("results/simple_RE.tsv", sep="\t")
complex_RE <- read.table("results/complex_RE.tsv", sep="\t")
hmm <- read.table("results/hmm.tsv")
US <- fromList( list(simple_RE=simple_RE$V1, complex_RE=complex_RE$V1, hmm=hmm$V3))
upset(US)
```



