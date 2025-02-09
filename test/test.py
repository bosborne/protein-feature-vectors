import sys
import os

srcdir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
sys.path.insert(0, srcdir)

from ProteinFeatureVectors import Protein  # noqa: E402

proteins = Protein()
proteins.display_feature_types()
proteins.get_feature_vectors("CTriad", file=f"{srcdir}/data_examples/multi.fa")
print(proteins.encodings)

seq = {
    "A12345": "LVTIKIGGQLKEALLDTGADDTVLEDMHLPGKWKPKMIGGIGGFIKVRQYDQILVEICGH"
    + "KAIGTVLVGPTPVNIIGRNLLTQIGCTLNFEMEKEGKISKIGPENPYNTPIFAIKKKDST"
    + "KWRKLVDFRELNKRTQDFWEVQLGIPHPAGLKKKKSVTVLDVGDAYFSVPLDEDFRKYTA"
    + "FTIPSTNNETPGIRYQYNVLPQGWKGSPAIFQSSMTKILEPFRKQNPDIVIYQYMDDLYV"
    + "GSDLEIGQHRIKVEELRQHLLRWGLTTPDKKHQKEPPFLWMG"
}

proteins.get_feature_vectors("TPC_type_2", pdict=seq)
print(proteins.encodings)
