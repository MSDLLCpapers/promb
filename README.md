# promb: protein humanness evaluation toolkit

```
              *           █
▄▄▄▄    ▄▄▄  ▄█▄   ▄▄▄▄   █▄▄▄  
█   █  █    █▓███  █ █ █  █   █ 
█▄▄▄▀  █    ▀███▀  █   █  █▄▄▄▀ 
█                           
▀       protein mutation burden 

# 1          2          3          4          5          6        
# CVQLQQSGA  VQLQQSGAE  QLQQSGAEL  LQQSGAELA  QQSGAELAR  QSGAELARP queries
# QVQLQQSGP  VQLMQSGAE  QLMQSGAEV  LHQSGSELA  QQSGSELAL  PSAAELARP hits
# ^       ^     ^         ^     +   ^   +         +   ^  ^ ^      
```

<p>
<a href="https://pypi.org/project/promb/">
    <img src="https://img.shields.io/pypi/dm/promb"
        alt="Pip Install"></a>
<a href="https://github.com/MSDLLCPapers/promb/releases">
    <img src="https://img.shields.io/pypi/v/promb"
        alt="Latest release"></a>
<a href="https://huggingface.co/spaces/prihodad/promb-humanness">
    <img src="https://img.shields.io/badge/🤗%20Spaces-prihodad/promb--humanness-blue"
        alt="Hugging Face Spaces"></a>
</p>

Python library for computing protein humanness (or mouse-ness, dog-ness, cow-ness, ...) based on average number of mutations 
to nearest peptide in a reference proteome (Human SwissProt, human antibody repertoires in OAS, or a custom reference DB).

## How does this work? What is protein humanness?

A database of human peptides (or other reference proteome) is loaded in-memory into a big python set for fast lookups.
Near matches are found using an optimized lookup strategy implemented with `cdist` function in `scipy`.

Human SwissProt and Human OAS are stored as part of this package. 
A custom reference DB can be provided as a gzipped fasta or text file.

> [!NOTE]
> More human doesn't always mean "good". Peptides like "GGGGGGGGG" or "EKEKEKEKE" are human but you don't necessarily want your protein to contain those. Make sure to check for sequence entropy, AlphaFold2 confidence, or other quality metrics.


## Usage

Try out promb on HuggingFace Spaces: https://huggingface.co/spaces/prihodad/promb-humanness


Install using pip:

```bash
pip install promb
```

### Command-line interface

Check which peptides in a sequence are present in the DB (exact match, very fast):

```bash
# FASTA file(s)
promb exact -d human-reference -l 9 -o output_matches.csv query.fasta
# PDB file(s)
promb exact -d human-reference -l 9 -o output_matches.csv --chain A ./pdb_directory/
# Output:
# File,ID,Position,Peptide,Found
# test_antibody.fa,Pembrolizumab,1,EIVLTQSPA,True
# test_antibody.fa,Pembrolizumab,2,IVLTQSPAT,True
# test_antibody.fa,Pembrolizumab,3,VLTQSPATL,True

```

Find nearest human peptides (a lot slower than exact match search). Peptides are ranked by identity (number of mutations), ties are resolved by BLOSUM similarity. The search is exhaustive - we don't use any pre-filtering step, all peptides from reference database are compared to the query.

```bash
# FASTA file(s)
promb nearest -d human-reference -l 9 -o output_matches.csv query.fasta
# PDB file(s)
promb nearest -d human-reference -l 9 -o output_matches.csv --chain A ./pdb_directory/
# Output:
# File,ID,Position,Peptide,Nearest,Mutations
# test_antibody.fa,Pembrolizumab,1,EIVLTQSPA,EIVLTQSPA,0
# test_antibody.fa,Pembrolizumab,2,IVLTQSPAT,IVLTQSPAT,0
# test_antibody.fa,Pembrolizumab,3,VLTQSPATL,VLTQSPATL,0
```

You can use promb for very fast computation of the OASis antibody humanness score from BioPhi:

```bash
# FASTA file(s)
promb oasis -o output_scores.csv ./input/sequences/
# Output:
# File,ID,Identity
# test_antibody.fa,Pembrolizumab,0.7572815533980582
# test_antibody.fa,Pembrolizumab,0.7142857142857143
```

### Python library

First, initialize the in-memory peptide database using `init_db`:

```python
from promb import init_db
# Use pre-defined Human reference proteome (UP000005640)
db = init_db('human-reference', 9)
# Use pre-defined Human SwissProt DB
db = init_db('human-swissprot', 9)
# Use pre-defined Human OAS DB (antibody peptides found in >=10% of subjects)
db = init_db('human-oas')
# Use custom path with gzipped fasta file
db = init_db('path/to/db.fasta.gz', 9)
# Use custom path with gzipped text file, one peptide per line
db = init_db('path/to/db.txt.gz')
```

Check which peptides in a sequence are present in the DB (exact match, very fast):

```python
from promb import init_db

db = init_db('human-oas')

db.contains('HWGRRKAWC')
# True
db.contains('MGPLHQFLL')
# False

# Use 'compute_peptide_content' to compute OASis-like score (% of peptides that are human, exact match)
db.compute_peptide_content('QVQLVQSGVEVKKPGASVKVSCKASGYTFTNYYMYWVRQAPGQGLEWMGGINPSNGGTNFNEKFKNRVTLTTDSSTTTAYMELKSLQFDDTAVYYCARRDYRFDMGFDYWGQGTTVTVSS')
# 0.7142857142857143
```

Find nearest human peptides (a lot slower than exact match search). Peptides are ranked by identity (number of mutations), ties are resolved by BLOSUM similarity. The search is exhaustive - we don't use any pre-filtering step, all peptides from reference database are compared to the query.

```python
from promb import init_db, print_nearest

db = init_db('human-reference', 9)

peptides = db.chop_seq_peptides('CVQLQQSGAELARPPASVKMSCKAS')
# ['CVQLQQSGA', 'VQLQQSGAE', 'QLQQSGAEL', ...]
nearest = db.find_nearest_peptides(peptides, 1)
# [['QVQLQQSGP'], ['VQLMQSGAE'], ['QLMQSGAEV'], ...]

print_nearest(peptides, nearest)
# ^ = mutation, + = similar residue
#
# 1          2          3          4          5          6          7          8          9          10       
# CVQLQQSGA  VQLQQSGAE  QLQQSGAEL  LQQSGAELA  QQSGAELAR  QSGAELARP  SGAELARPP  GAELARPPA  AELARPPAS  ELARPPASV queries
# QVQLQQSGP  VQLMQSGAE  QLMQSGAEV  LHQSGSELA  QQSGSELAL  PSAAELARP  SGAELRQPP  GAELLEPPA  LELQRPPAS  ELQRPPAST hits
# ^       ^     ^         ^     +   ^   +         +   ^  ^ ^             ^+        ^^     ^  ^         ^     ^

# 11         12         13         14         15         16         17       
# LARPPASVK  ARPPASVKM  RPPASVKMS  PPASVKMSC  PASVKMSCK  ASVKMSCKA  SVKMSCKAS queries
# LAMPPASVK  AMPPASVKV  RPPASVGMS  PPASYKSSC  VASTKMSCK  ASVKVSCKA  SVKVSCKAS hits
#   ^         ^      +        ^        ^ ^    ^  ^           +         +     

```

Compute protein mutation burden (PMB) = average number of non-human mutations per peptide:

```python
pmb = db.compute_pmb('CVQLQQSGAELARPPASVKMSCKAS')
pmb
# 1.176 mutations per peptide
```

PMB is slow for peptides with more than 2 mutations, you can cap the number of mutations per peptide to speed up computation:

```python
pmb = db.compute_pmb('ELVISISALIVE', max=3)
pmb
# 2.75 mutations per peptide (capped at 3 - peptides with more mutations will be counted as 3 mutations)
```

## Acknowledgements

Promb is building upon our previous work - OASis antibody humanness metric:

> David Prihoda, Jad Maamary, Andrew Waight, Veronica Juan, Laurence Fayadat-Dilman, Daniel Svozil & Danny A. Bitton (2022) 
> BioPhi: A platform for antibody design, humanization, and humanness evaluation based on natural antibody repertoires and deep learning, mAbs, 14:1, DOI: https://doi.org/10.1080/19420862.2021.2020203

OASis is based on antibody repertoires from the Observed Antibody Space:

> Kovaltsuk, A., Leem, J., Kelm, S., Snowden, J., Deane, C. M., & Krawczyk, K. (2018). Observed Antibody Space: A Resource for Data Mining Next-Generation Sequencing of Antibody Repertoires. The Journal of Immunology, 201(8), 2502–2509. https://doi.org/10.4049/jimmunol.1800708
