# *protein-feature-vectors*

Code to produce fixed-length feature vectors from protein sequences based on iFeatureOmega.

## Introduction

This package is derived from [iFeatureOmega-CLI](https://github.com/Superzchen/iFeatureOmega-CLI), with the following differences:

- It only makes feature vectors using protein sequence (not based on DNA, RNA, or compounds)
- It only makes fixed-length feature vectors
- It has fewer dependencies than iFeatureOmega-CLI
- It just creates feature vectors and has no other plotting or analytical capabilities
- It has additional variants of specific algorithms (e.g. TPC_type_3) created for faster Deep Learning

## Installation

```sh
  git clone git@github.com:bosborne/protein-feature-vectors.git
  cd protein-feature-vectors
  pip3 install .
```

## Usage

Sequences from a Fasta file:

```python
from protein_feature_vectors import Calculator
calc = Calculator()
calc.display_feature_types()
calc.get_feature_vectors("CTriad", file='data_examples/multi.fa')
print(calc.encodings)
calc.to_csv("CTriad.csv", "index=False", header=False)
```

Sequences from a dict:

```python
from protein_feature_vectors import Calculator
seqs = {"A1": "MLVTIKIQLKEAL...LDTGADVLEDM", "B2": "MHLPGKWMIGGIG....GFIKVRQYDEICGH"}
calc = Calculator()
calc.get_feature_vectors("CTriad", pdict=seqs)
print(calc.encodings)
```

## Feature Vectors

| Name | Description | Length |
|------|-------------|--------|
| AAC | Amino acid composition | 20 |
| AC | Auto covariance | 24 |
| ACC | Auto-cross covariance | 192 |
| APAAC | Amphiphilic PAAC | 26 |
| ASDC | Adaptive skip dipeptide composition | 400 |
| CC | Cross covariance | 168 |
| CKSAAGP_type_1 | Composition of k-spaced amino acid group pairs - normalized | 100 |
| CKSAAGP_type_2 | Composition of k-spaced amino acid group pairs - raw count | 100 |
| CKSAAP_type_1 | k-spaced amino acid pairs, k = 0,1,2,3 - normalized | 1600 |
| CKSAAP_type_2 | k-spaced amino acid pairs, k = 0,1,2,3 - raw count | 1600 |
| CKSAAP_type_3 | k-spaced amino acid pairs, k = 0,1,2,3 - normalized, rounded count | 1600 |
| CTDC | Composition | 39 |
| CTDD | Distribution | 195 |
| CTDT | Transition | 39 |
| CTriad | Triads of 7 groups based on dipoles and volumes | 343 |
| DDE | Composition and distribution of dipeptides | 400 |
| DistancePair | PseAAC of distance-pairs and reduced alphabet | 20 |
| DPC_type_1 | Dipeptide composition - normalized | 400 |
| DPC_type_2 | Dipeptide composition - raw count | 400 |
| GAAC | Grouped amino acid composition | 5 |
| GDPC_type_1 | Grouped dipeptide composition - normalized | 25 |
| GDPC_type_2 | Grouped dipeptide composition - raw count | 25 |
| Geary | Spatial autocorrelation of physicochemical properties | 24 |
| GTPC_type_1 | Grouped tripeptide composition - normalized | 125 |
| GTPC_type_2 | Grouped tripeptide composition - raw count | 125 |
| KSCTriad | k-spaced CTriad, k = 0,1,2,3 | 1372 |
| Moran | Spatial autocorrelation of physicochemical properties | 24 |
| NMBroto | Normalized Moreau-Broto | 24 |
| PAAC | Pseudo-amino acid composition | 23 |
| PseKRAAC_type_1 | Pseudo K-tuple reduced amino acids composition | 4 |
| PseKRAAC_type_2 | Pseudo K-tuple reduced amino acids composition | 4 |
| PseKRAAC_type_3A | Pseudo K-tuple reduced amino acids composition | 4 |
| PseKRAAC_type_3B | Pseudo K-tuple reduced amino acids composition | 4 |
| PseKRAAC_type_4 | Pseudo K-tuple reduced amino acids composition | 25 |
| PseKRAAC_type_5 | Pseudo K-tuple reduced amino acids composition | 9 |
| PseKRAAC_type_6A | Pseudo K-tuple reduced amino acids composition | 16 |
| PseKRAAC_type_6B | Pseudo K-tuple reduced amino acids composition | 25 |
| PseKRAAC_type_6C | Pseudo K-tuple reduced amino acids composition | 25 |
| PseKRAAC_type_7 | Pseudo K-tuple reduced amino acids composition | 4 |
| PseKRAAC_type_8 | Pseudo K-tuple reduced amino acids composition | 4 |
| PseKRAAC_type_9 | Pseudo K-tuple reduced amino acids composition | 4 |
| PseKRAAC_type_10 | Pseudo K-tuple reduced amino acids composition | 4 |
| PseKRAAC_type_11 | Pseudo K-tuple reduced amino acids composition | 4 |
| PseKRAAC_type_12 | Pseudo K-tuple reduced amino acids composition | 4 |
| PseKRAAC_type_13 | Pseudo K-tuple reduced amino acids composition | 16 |
| PseKRAAC_type_14 | Pseudo K-tuple reduced amino acids composition | 4 |
| PseKRAAC_type_15 | Pseudo K-tuple reduced amino acids composition | 4 |
| PseKRAAC_type_16 | Pseudo K-tuple reduced amino acids composition | 4 |
| QSOrder | Quasi-sequence-order descriptors (nlag=13) | 46 |
| SOCNumber | Sequence-order-coupling number | 6 |
| TPC_type_1 | Tripeptide composition - normalized | 8000  |
| TPC_type_2 | Tripeptide composition - raw count | 8000 |
| TPC_type_3 | Tripeptide composition - normalized, rounded count | 8000 |
