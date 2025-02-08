# *protein-feature-vectors*

Code to produce fixed-length feature vectors from protein sequences.

## Introduction

This package is derived from [iFeatureOmega-CLI](https://github.com/Superzchen/iFeatureOmega-CLI), with the following differences:

- It only makes feature vectors using protein sequence (not based on DNA, RNA, or chemical structures)
- The feature vectors are all fixed-length
- It has fewer dependencies than iFeatureOmega-CLI
- It just creates feature vectors and has no other plotting or analytical capabilities

## Installation

```sh  
  git clone git@github.com:bosborne/protein-feature-vectors.git
  cd protein-feature-vectors
  pip3 install .
```

### Usage

```python
  $ python3
  >>> from ProteinFeatureVectors import Protein
  >>> proteins.display_feature_types()
  >>> proteins = Protein('data_examples/multi.fa')
  >>> proteins.get_feature_vectors("CTriad")
  >>> print(proteins.encodings)
```
