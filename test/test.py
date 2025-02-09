import sys
import os

srcdir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
sys.path.insert(0, srcdir)

from ProteinFeatureVectors import Calculator  # noqa: E402

proteins = Calculator()
proteins.display_feature_types()
proteins.get_feature_vectors("CTriad", file=f"{srcdir}/data_examples/multi.fa")
print(proteins.encodings)

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
}

proteins.get_feature_vectors("TPC_type_2", pdict=seqs)
print(proteins.encodings)
