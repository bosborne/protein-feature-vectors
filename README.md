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

```python
  $ python3
  >>> from ProteinFeatureVectors import Protein
  >>> proteins = Protein('data_examples/multi.fa')
  >>> proteins.display_feature_types()
  >>> proteins.get_feature_vectors("CTriad")
  >>> print(proteins.encodings)
  >>> proteins.to_csv("CTriad.csv", "index=False", header=False)
```

## Feature Vectors

| Name | Description | Length |
|------|-------------|--------|
| AAC | Amino acid composition |   |
| CKSAAP_type_1 | Composition of k-spaced amino acid pairs type 1 - normalized |   |
| CKSAAP_type_2 | Composition of k-spaced amino acid pairs type 2 - raw count |   |
| DPC_type_1 | Dipeptide composition type 1 - normalized |   |
| DPC_type_2 | Dipeptide composition type 2 - raw count |   |
| TPC_type_1 | Tripeptide composition type 1 - normalized |   |
| TPC_type_2 | Tripeptide composition type 1 - raw count |   |
| CTDC | Composition |   |
| CTDT | Transition |   |
| CTDD | Distribution |   |
| CTriad | Conjoint triad |   |
| KSCTriad | Conjoint k-spaced triad |   |
| ASDC | Adaptive skip dipeptide composition |   |
| DistancePair | PseAAC of distance-pairs and reduced alphabe |   |
| GAAC | Grouped amino acid composition |   |
| CKSAAGP_type_1 | Composition of k-spaced amino acid group pairs type 1- normalized |   |
| CKSAAGP_type_2 | Composition of k-spaced amino acid group pairs type 2- raw count |   |
| GDPC_type_1 | Grouped dipeptide composition type 1 - normalized |   |
| GDPC_type_2 | Grouped dipeptide composition type 2 - raw count |   |
| GTPC_type_1 | Grouped tripeptide composition type 1 - normalized |   |
| GTPC_type_2 | Grouped tripeptide composition type 1 - raw count |   |
| Moran | Spatial autocorrelation of physicochemical properties |   |
| Geary | Spatial autocorrelation of physicochemical properties |   |
| NMBroto | Normalized Moreau-Broto |   |
| AC | Auto covariance |   |
| CC | Cross covariance |   |
| ACC | Auto-cross covariance |   |
| SOCNumber | Sequence-order-coupling number |   |
| QSOrder | Quasi-sequence-order descriptors |   |
| PAAC | Pseudo-amino acid composition |   |
| APAAC | Amphiphilic PAAC |   |
| PseKRAAC_type_1 | Pseudo K-tuple reduced amino acids composition type 1 |   |
| PseKRAAC_type_2 | Pseudo K-tuple reduced amino acids composition type 2 |   |
| PseKRAAC_type_3A | Pseudo K-tuple reduced amino acids composition type 3A |   |
| PseKRAAC_type_3B | Pseudo K-tuple reduced amino acids composition type 3B |   |
| PseKRAAC_type_4 | Pseudo K-tuple reduced amino acids composition type 4 |   |
| PseKRAAC_type_5 | Pseudo K-tuple reduced amino acids composition type 5 |   |
| PseKRAAC_type_6A | Pseudo K-tuple reduced amino acids composition type 6A |   |
| PseKRAAC_type_6B | Pseudo K-tuple reduced amino acids composition type 6B |   |
| PseKRAAC_type_6C | Pseudo K-tuple reduced amino acids composition type 6C |   |
| PseKRAAC_type_7 | Pseudo K-tuple reduced amino acids composition type 7 |   |
| PseKRAAC_type_8 | Pseudo K-tuple reduced amino acids composition type 8 |   |
| PseKRAAC_type_9 | Pseudo K-tuple reduced amino acids composition type 9 |   |
| PseKRAAC_type_10 | Pseudo K-tuple reduced amino acids composition type 10 |   |
| PseKRAAC_type_11 | Pseudo K-tuple reduced amino acids composition type 11 |   |
| PseKRAAC_type_12 | Pseudo K-tuple reduced amino acids composition type 12 |   |
| PseKRAAC_type_13 | Pseudo K-tuple reduced amino acids composition type 13 |   |
| PseKRAAC_type_14 | Pseudo K-tuple reduced amino acids composition type 14 |   |
| PseKRAAC_type_15 | Pseudo K-tuple reduced amino acids composition type 15 |   |
| PseKRAAC_type_16 | Pseudo K-tuple reduced amino acids composition type 16 |   |
