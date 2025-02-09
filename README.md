# *protein-feature-vectors*

Code to produce fixed-length feature vectors from protein sequences based on iFeatureOmega.

## Introduction

This package is derived from [iFeatureOmega-CLI](https://github.com/Superzchen/iFeatureOmega-CLI), with the following differences:

- It only makes feature vectors using protein sequence (not based on DNA, RNA, or chemical structures)
- It only makes fixed-length feature vectors
- It has fewer dependencies than iFeatureOmega-CLI
- It just creates feature vectors and has no other plotting or analytical capabilities

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
proteins = Calculator()
proteins.display_feature_types()
proteins.get_feature_vectors("CTriad", file='data_examples/multi.fa')
print(proteins.encodings)
proteins.to_csv("CTriad.csv", "index=False", header=False)
```

Sequences from a dict:

```python
from protein_feature_vectors import Calculator
seqs = {"A1": "MLVTIKIGGQLKEAL...LDTGADDTVLEDM", "B2": "MHLPGKWKPKMIGGIG....GFIKVRQYDQILVEICGH"}
proteins = Calculator()
proteins.get_feature_vectors("CTriad", pdict=seqs)
print(proteins.encodings)
```

## Feature Vectors

| Name | Description | Length |
|------|-------------|--------|
| AAC | Amino acid composition | 20 |
| CKSAAP_type_1 | Composition of k-spaced amino acid pairs type 1 - normalized | 1600 |
| CKSAAP_type_2 | Composition of k-spaced amino acid pairs type 2 - raw count | 1600 |
| DPC_type_1 | Dipeptide composition type 1 - normalized | 400 |
| DPC_type_2 | Dipeptide composition type 2 - raw count | 400 |
| TPC_type_1 | Tripeptide composition type 1 - normalized | 8000  |
| TPC_type_2 | Tripeptide composition type 1 - raw count | 8000  |
| CTDC | Composition | 39 |
| CTDT | Transition | 39 |
| CTDD | Distribution | 195 |
| CTriad | Conjoint triad | 343 |
| KSCTriad | Conjoint k-spaced triad | 1372 |
| ASDC | Adaptive skip dipeptide composition | 400 |
| DistancePair | PseAAC of distance-pairs and reduced alphabe | 20 |
| GAAC | Grouped amino acid composition | 5 |
| CKSAAGP_type_1 | Composition of k-spaced amino acid group pairs type 1- normalized | 100 |
| CKSAAGP_type_2 | Composition of k-spaced amino acid group pairs type 2- raw count | 100 |
| GDPC_type_1 | Grouped dipeptide composition type 1 - normalized | 25 |
| GDPC_type_2 | Grouped dipeptide composition type 2 - raw count | 25 |
| GTPC_type_1 | Grouped tripeptide composition type 1 - normalized | 125 |
| GTPC_type_2 | Grouped tripeptide composition type 1 - raw count | 125 |
| Moran | Spatial autocorrelation of physicochemical properties | 24 |
| Geary | Spatial autocorrelation of physicochemical properties | 24 |
| NMBroto | Normalized Moreau-Broto | 24 |
| AC | Auto covariance | 24 |
| CC | Cross covariance | 168 |
| ACC | Auto-cross covariance | 192 |
| SOCNumber | Sequence-order-coupling number | 6 |
| QSOrder | Quasi-sequence-order descriptors | 46 |
| PAAC | Pseudo-amino acid composition | 23 |
| APAAC | Amphiphilic PAAC | 26 |
| PseKRAAC_type_1 | Pseudo K-tuple reduced amino acids composition type 1 | 4 |
| PseKRAAC_type_2 | Pseudo K-tuple reduced amino acids composition type 2 | 4 |
| PseKRAAC_type_3A | Pseudo K-tuple reduced amino acids composition type 3A | 4 |
| PseKRAAC_type_3B | Pseudo K-tuple reduced amino acids composition type 3B | 4 |
| PseKRAAC_type_4 | Pseudo K-tuple reduced amino acids composition type 4 | 25 |
| PseKRAAC_type_5 | Pseudo K-tuple reduced amino acids composition type 5 | 9 |
| PseKRAAC_type_6A | Pseudo K-tuple reduced amino acids composition type 6A | 16 |
| PseKRAAC_type_6B | Pseudo K-tuple reduced amino acids composition type 6B | 25 |
| PseKRAAC_type_6C | Pseudo K-tuple reduced amino acids composition type 6C | 25 |
| PseKRAAC_type_7 | Pseudo K-tuple reduced amino acids composition type 7 | 4 |
| PseKRAAC_type_8 | Pseudo K-tuple reduced amino acids composition type 8 | 4 |
| PseKRAAC_type_9 | Pseudo K-tuple reduced amino acids composition type 9 | 4 |
| PseKRAAC_type_10 | Pseudo K-tuple reduced amino acids composition type 10 | 4 |
| PseKRAAC_type_11 | Pseudo K-tuple reduced amino acids composition type 11 | 4 |
| PseKRAAC_type_12 | Pseudo K-tuple reduced amino acids composition type 12 | 4 |
| PseKRAAC_type_13 | Pseudo K-tuple reduced amino acids composition type 13 | 16 |
| PseKRAAC_type_14 | Pseudo K-tuple reduced amino acids composition type 14 | 4 |
| PseKRAAC_type_15 | Pseudo K-tuple reduced amino acids composition type 15 | 4 |
| PseKRAAC_type_16 | Pseudo K-tuple reduced amino acids composition type 16 | 4 |
