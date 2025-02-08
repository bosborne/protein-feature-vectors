import sys
import os

srcdir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
sys.path.insert(0, srcdir)

from ProteinFeatureVectors import Protein  # noqa: E402

proteins = Protein(f"{srcdir}/data_examples/multi.fa")
proteins.display_feature_types()
proteins.get_feature_vectors("CTriad")
print(proteins.encodings)
proteins.get_feature_vectors("TPC_type_2")
print(proteins.encodings)
