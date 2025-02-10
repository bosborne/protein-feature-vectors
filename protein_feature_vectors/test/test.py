import sys
import os
import pandas as pd

srcdir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
sys.path.insert(0, srcdir)

from protein_feature_vectors import Calculator  # noqa: E402

proteins = Calculator()
proteins.display_feature_types()
proteins.get_feature_vectors(
    "QSOrder", file=f"{srcdir}/data_examples/multi.fa"
)
assert isinstance(proteins.encodings, pd.DataFrame)

seqs = {
    "A12345": "LVTIKIGGQLKEALLDTGADDTVLEDMHLPGKWKPKMIGGIGGFIKVRQYDQILVEICGH"
    + "KAIGTVLVGPTPVNIIGRNLLTQIGCTLNFEMEKEGKISKIGPENPYNTPIFAIKKKDST"
    + "KWRKLVDFRELNKRTQDFWEVQLGIPHPAGLKKKKSVTVLDVGDAYFSVPLDEDFRKYTA"
    + "FTIPSTNNETPGIRYQYNVLPQGWKGSPAIFQSSMTKILEPFRKQNPDIVIYQYMDDLYV"
    + "GSDLEIGQHRIKVEELRQHLLRWGLTTPDKKHQKEPPFLWMG",
    "B3QQG6": "MQETSATITAERKHSHVDVCLNRPVCFDGQDTGLDAWRFEHNAAPEIDFAEIDLTA"
    + "EFLGHAIGMPLMISSMTGGYGDALALNRTLAEAAERFRIPLGVGSMRQALEGNSHRESFSIVRSSA"
    + "PSVPIFANIGAPEVAAGLSREQLSTLVELIEADGLIVHLNPAQELFQPEGSTNFRGFLDRLHDITA"
    + "TINVPVIAKEVGCGISAPLASKLADAGVKAIDVAGAGGISWQKVEECRYLDRFGNEERFSPSALDEF"
    + "LNWGIPTAECLTGIAALKEKSPEYGSLAVISSGGIRNGLDVAKSIALGADIAASAQHLLKALRAGTLEE"
    + "TIRTWANDLRAAMFLTGSATTAQLKHAPIYRKP",
    # "short": "M",
}
for alg in [
    "AAC",
    "AC",
    "ACC",
    "APAAC",
    "ASDC",
    "CC",
    "CKSAAGP_type_1",
    "CKSAAGP_type_2",
    "CKSAAP_type_1",
    "CKSAAP_type_2",
    "CTDC",
    "CTDD",
    "CTDT",
    "CTriad",
    "DistancePair",
    "DDE",
    "DPC_type_1",
    "DPC_type_2",
    "GAAC",
    "GDPC_type_1",
    "GDPC_type_2",
    "Geary",
    "GTPC_type_1",
    "GTPC_type_2",
    "KSCTriad",
    "Moran",
    "NMBroto",
    "PAAC",
    "PseKRAAC_type_1",
    "PseKRAAC_type_2",
    "PseKRAAC_type_3A",
    "PseKRAAC_type_3B",
    "PseKRAAC_type_4",
    "PseKRAAC_type_5",
    "PseKRAAC_type_6A",
    "PseKRAAC_type_6B",
    "PseKRAAC_type_6C",
    "PseKRAAC_type_7",
    "PseKRAAC_type_8",
    "PseKRAAC_type_9",
    "PseKRAAC_type_10",
    "PseKRAAC_type_11",
    "PseKRAAC_type_12",
    "PseKRAAC_type_13",
    "PseKRAAC_type_14",
    "PseKRAAC_type_15",
    "PseKRAAC_type_16",
    "QSOrder",
    "SOCNumber",
    "TPC_type_1",
    "TPC_type_2",
]:
    proteins.get_feature_vectors(alg, pdict=seqs)
    print(proteins.encodings)
    assert isinstance(proteins.encodings, pd.DataFrame)
    assert proteins.encodings.shape[0] == 2
