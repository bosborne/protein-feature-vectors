import os
import sys
import json
import math
import pickle
import numpy as np
import pandas as pd
import fastapy
from collections import Counter


class Calculator:
    """
    from protein_feature_vectors import Calculator

    # Create a instance
    calc = Calculator()

    # Display available methods
    calc.display_feature_types()

    # Calculate feature vectors
    calc.get_feature_vectors("AAC", file="./data_examples/multi.fa")

    # Display the feature vectors
    print(calc.encodings)

    # Get ids and values from the *encodings* DataFrame
    protein_ids = [x[0] for x in calc.encodings.iterrows()]
    values = [x[1].tolist() for x in calc.encodings.iterrows()]

    # Save feature vectors
    calc.to_csv("AAC.csv", index=False, header=False)
    """

    def __init__(self, verbose=False):
        self.verbose = verbose
        self.AA = "ACDEFGHIKLMNPQRSTVWY"
        self.non_AA = "BJOUXZ-_"
        self.seq_list = list()
        self.vector_length = dict()
        self.datadir = os.path.join(
            os.path.dirname(os.path.realpath(__file__)), "data"
        )
        self.import_parameters("Protein_parameters_setting.json")
        self.__default_para = {
            "sliding_window": 5,
            "kspace": 3,
            "nlag": 3,
            "weight": 0.05,
            "lambdaValue": 3,
            "PseKRAAC_model": "g-gap",
            "g-gap": 2,
            "k-tuple": 2,
            "RAAC_clust": 1,
        }
        self.__cmd_dict = {
            "AAC": "self._AAC()",
            "AC": "self._AC()",
            "ACC": "self._ACC()",
            "APAAC": "self._APAAC()",
            "ASDC": "self._ASDC()",
            "CC": "self._CC()",
            "CKSAAGP_type_1": "self._CKSAAGP(normalized=True)",
            "CKSAAGP_type_2": "self._CKSAAGP(normalized=False)",
            "CKSAAP_type_1": "self._CKSAAP(type=1)",  # normalized=True
            "CKSAAP_type_2": "self._CKSAAP(type=2)",  # normalized=False
            "CKSAAP_type_3": "self._CKSAAP(type=3)",
            "CTDC": "self._CTDC()",
            "CTDD": "self._CTDD()",
            "CTDT": "self._CTDT()",
            "CTriad": "self._CTriad()",
            "DDE": "self._DDE()",
            "DistancePair": "self._DistancePair()",
            "DPC_type_1": "self._DPC(normalized=True)",
            "DPC_type_2": "self._DPC(normalized=False)",
            "GAAC": "self._GAAC()",
            "GDPC_type_1": "self._GDPC(normalized=True)",
            "GDPC_type_2": "self._GDPC(normalized=False)",
            "Geary": "self._Geary()",
            "GTPC_type_1": "self._GTPC(normalized=True)",
            "GTPC_type_2": "self._GTPC(normalized=False)",
            "KSCTriad": "self._KSCTriad()",
            "Moran": "self._Moran()",
            "NMBroto": "self._NMBroto()",
            "PAAC": "self._PAAC()",
            "PseKRAAC_type_1": "self._PseKRAAC_type_1()",
            "PseKRAAC_type_2": "self._PseKRAAC_type_2()",
            "PseKRAAC_type_3A": "self._PseKRAAC_type_3A()",
            "PseKRAAC_type_3B": "self._PseKRAAC_type_3B()",
            "PseKRAAC_type_4": "self._PseKRAAC_type_4()",
            "PseKRAAC_type_5": "self._PseKRAAC_type_5()",
            "PseKRAAC_type_6A": "self._PseKRAAC_type_6A()",
            "PseKRAAC_type_6B": "self._PseKRAAC_type_6B()",
            "PseKRAAC_type_6C": "self._PseKRAAC_type_6C()",
            "PseKRAAC_type_7": "self._PseKRAAC_type_7()",
            "PseKRAAC_type_8": "self._PseKRAAC_type_8()",
            "PseKRAAC_type_9": "self._PseKRAAC_type_9()",
            "PseKRAAC_type_10": "self._PseKRAAC_type_10()",
            "PseKRAAC_type_11": "self._PseKRAAC_type_11()",
            "PseKRAAC_type_12": "self._PseKRAAC_type_12()",
            "PseKRAAC_type_13": "self._PseKRAAC_type_13()",
            "PseKRAAC_type_14": "self._PseKRAAC_type_14()",
            "PseKRAAC_type_15": "self._PseKRAAC_type_15()",
            "PseKRAAC_type_16": "self._PseKRAAC_type_16()",
            "QSOrder": "self._QSOrder()",
            "SOCNumber": "self._SOCNumber()",
            "TPC_type_1": "self._TPC(type=1)",
            "TPC_type_2": "self._TPC(type=2)",
            "TPC_type_3": "self._TPC(type=3)",
        }

    def read_fasta(self):
        """read_fasta
        >>>calc.fasta_list
        [['H6SH45', 'MFFRENLAFQQREARKFSS....DLIAEIQKQGQGQWTYQIYQE'],
        ['9HIV1', 'MDFWEVQLGIPHP...YQEIVTLTEEAELELAENREI'],
        ['H6SH45_9HIV1', 'MFFRENLAFQQREARKF...WTYQIYQEIELELAENREI']]
        """
        fasta_sequences = []
        if not os.path.exists(self.fasta_file):
            sys.exit(f"Error: cannot find Fasta file '{self.fasta_file}'")
        for record in fastapy.parse(self.fasta_file):
            fasta_sequences.append([record.id, str(record.seq).upper()])

        self.seq_list = fasta_sequences

    def import_parameters(self, file):
        """import_parameters

        Parameters
        ----------
        file : str
            Input JSON settings file
        """
        settings = os.path.join(
            os.path.dirname(os.path.realpath(__file__)), file
        )
        if os.path.exists(settings):
            with open(settings) as f:
                records = f.read().strip()
            self.__default_para_dict = json.loads(records)
        else:
            sys.exit(f"Parameter file '{settings}' not found")

    def get_feature_vectors(self, descriptor, file=None, pdict=None):
        """get_feature_vectors

        Parameters
        ----------
        descriptor : string
            Name of feature vector
        file : sequence file name, optional
            Fasta file name
        pdict : dict, optional
            Key is an id, value is a sequence string
        """
        if descriptor not in self.__cmd_dict:
            sys.exit(f"No such algorithm: {descriptor}")

        if file is not None:
            self.fasta_file = file
            self.read_fasta()
        elif pdict is not None and isinstance(pdict, dict):
            self.seq_list = [[id, sequence] for id, sequence in pdict.items()]
        elif len(self.seq_list) == 0:
            sys.exit("No sequence supplied")
        # Remove sequences with invalid chars
        self.validate()

        self.sequence_number = len(self.seq_list)
        self.encodings = None
        # Copy parameters
        if descriptor in self.__default_para_dict:
            for key in self.__default_para_dict[descriptor]:
                self.__default_para[key] = self.__default_para_dict[
                    descriptor
                ][key]
        # Run the command
        eval(self.__cmd_dict[descriptor])
        if len(self.encodings > 0):
            self.vector_length[descriptor] = self.encodings.shape[1]
            if self.verbose:
                print(
                    f"{descriptor} vector length: {self.vector_length[descriptor]}"
                )

    def validate(self):
        """validate
        Remove sequences with invalid characters and convert to uppercase
        """
        validated = list()
        for fasta in self.seq_list:
            seqstr = fasta[1].upper()
            invalid = self.check_for_nonstandard(seqstr)
            if len(invalid) == 0:
                validated.append([fasta[0], seqstr])
            else:
                if self.verbose:
                    print(f"Skipping {fasta[0]} non-standard aa: {invalid}")
        self.seq_list = validated

    def check_for_nonstandard(self, seqstr):
        """check_for_nonstandard
        '&' is intersection

        Parameters
        ----------
        seqstr : str
            AA sequence

        Returns
        -------
        set
            Set of invalid characters
        """
        return set(seqstr) & set(self.non_AA)

    def display_feature_types(self):
        info = """        
        AAC                                                Amino acid composition
        AC                                                 Auto covariance
        ACC                                                Auto-cross covariance
        APAAC                                              Amphiphilic PAAC
        ASDC                                               Adaptive skip dipeptide composition
        CC                                                 Cross covariance
        CKSAAGP_type_1                                     Composition of k-spaced amino acid group pairs type 1 - normalized
        CKSAAGP_type_2                                     Composition of k-spaced amino acid group pairs type 2 - raw count
        CKSAAP_type_1                                      Composition of k-spaced amino acid pairs - normalized
        CKSAAP_type_2                                      Composition of k-spaced amino acid pairs - raw count
        CKSAAP_type_3                                      Composition of k-spaced amino acid pairs - normalized, rounded count
        CTDC                                               Composition
        CTDD                                               Distribution
        CTDT                                               Transition
        CTriad                                             Conjoint triad
        DistancePair                                       PseAAC of distance-pairs and reduced alphabet
        DDE                                                Composition and distribution of dipeptides
        DPC_type_1                                         Dipeptide composition type 1 - normalized
        DPC_type_2                                         Dipeptide composition type 2 - raw count
        GAAC                                               Grouped amino acid composition
        GDPC_type_1                                        Grouped dipeptide composition type 1 - normalized
        GDPC_type_2                                        Grouped dipeptide composition type 2 - raw count
        Geary                                              Geary
        GTPC_type_1                                        Grouped tripeptide composition type 1 - normalized
        GTPC_type_2                                        Grouped tripeptide composition type 1 - raw count
        KSCTriad                                           Conjoint k-spaced triad
        Moran                                              Correlation between neighboring residues based on physicochemical properties
        NMBroto                                            Correlation between neighboring residues based on physicochemical properties
        PAAC                                               Pseudo-amino acid composition
        PseKRAAC_type_1                                    Pseudo K-tuple reduced amino acids composition type 1
        PseKRAAC_type_2                                    Pseudo K-tuple reduced amino acids composition type 2
        PseKRAAC_type_3A                                   Pseudo K-tuple reduced amino acids composition type 3A
        PseKRAAC_type_3B                                   Pseudo K-tuple reduced amino acids composition type 3B
        PseKRAAC_type_4                                    Pseudo K-tuple reduced amino acids composition type 4
        PseKRAAC_type_5                                    Pseudo K-tuple reduced amino acids composition type 5
        PseKRAAC_type_6A                                   Pseudo K-tuple reduced amino acids composition type 6A
        PseKRAAC_type_6B                                   Pseudo K-tuple reduced amino acids composition type 6B
        PseKRAAC_type_6C                                   Pseudo K-tuple reduced amino acids composition type 6C
        PseKRAAC_type_7                                    Pseudo K-tuple reduced amino acids composition type 7
        PseKRAAC_type_8                                    Pseudo K-tuple reduced amino acids composition type 8
        PseKRAAC_type_9                                    Pseudo K-tuple reduced amino acids composition type 9
        PseKRAAC_type_10                                   Pseudo K-tuple reduced amino acids composition type 10
        PseKRAAC_type_11                                   Pseudo K-tuple reduced amino acids composition type 11
        PseKRAAC_type_12                                   Pseudo K-tuple reduced amino acids composition type 12
        PseKRAAC_type_13                                   Pseudo K-tuple reduced amino acids composition type 13
        PseKRAAC_type_14                                   Pseudo K-tuple reduced amino acids composition type 14
        PseKRAAC_type_15                                   Pseudo K-tuple reduced amino acids composition type 15
        PseKRAAC_type_16                                   Pseudo K-tuple reduced amino acids composition type 16        
        QSOrder                                            Quasi-sequence-order descriptors
        SOCNumber                                          Sequence-order-coupling number
        TPC_type_1                                         Tripeptide composition - normalized
        TPC_type_2                                         Tripeptide composition - raw count
        TPC_type_3                                         Tripeptide composition - normalized, rounded count
        """
        print(info)

    def _AAC(self):
        try:
            header = ["SampleName"]
            encodings = []
            for i in self.AA:
                header.append("AAC_{0}".format(i))
            encodings.append(header)
            for i in self.seq_list:
                name, sequence = i[0], i[1]
                count = Counter(sequence)
                for key in count:
                    count[key] = count[key] / len(sequence)
                code = [name]
                for aa in self.AA:
                    code.append(count[aa])
                encodings.append(code)
            encodings = np.array(encodings)
            self.encodings = pd.DataFrame(
                encodings[1:, 1:].astype(float),
                columns=encodings[0, 1:],
                index=encodings[1:, 0],
            )
            return True
        except Exception as e:
            sys.exit(f"Error: {e}")

    def _CKSAAP(self, type=None):
        try:
            encodings = []
            aaPairs = []
            for aa1 in self.AA:
                for aa2 in self.AA:
                    aaPairs.append(aa1 + aa2)
            header = ["SampleName"]
            gap = self.__default_para["kspace"]
            for g in range(gap + 1):
                for aa in aaPairs:
                    header.append(f"CKSAAP{type}_" + aa + ".gap" + str(g))
            encodings.append(header)
            for i in self.seq_list:
                name, sequence = i[0], i[1]
                code = [name]
                for g in range(gap + 1):
                    myDict = {}
                    for pair in aaPairs:
                        myDict[pair] = 0
                    sum = 0
                    for index1 in range(len(sequence)):
                        index2 = index1 + g + 1
                        if (
                            index1 < len(sequence)
                            and index2 < len(sequence)
                            and sequence[index1] in self.AA
                            and sequence[index2] in self.AA
                        ):
                            myDict[sequence[index1] + sequence[index2]] = (
                                myDict[sequence[index1] + sequence[index2]] + 1
                            )
                            sum = sum + 1
                    for pair in aaPairs:
                        if type == 1:
                            code.append(myDict[pair] / sum)
                        elif type == 2:
                            code.append(myDict[pair])
                        # If the multiplier is 1000 the F1 score drops dramaticaly.
                        # Example values now: 2.8, 1.7, 0.4
                        elif type == 3:
                            code.append(
                                round(
                                    (myDict[pair] * (100 / len(sequence))),
                                    1,
                                )
                            )
                encodings.append(code)
            encodings = np.array(encodings)
            self.encodings = pd.DataFrame(
                encodings[1:, 1:].astype(float),
                columns=encodings[0, 1:],
                index=encodings[1:, 0],
            )
            return True
        except Exception as e:
            sys.exit(f"Error: {e}")

    def _DPC(self, normalized=True):
        try:
            encodings = []
            diPeptides = [
                "DPC_" + aa1 + aa2 for aa1 in self.AA for aa2 in self.AA
            ]
            header = ["SampleName"] + diPeptides
            encodings.append(header)
            AADict = {}
            for i in range(len(self.AA)):
                AADict[self.AA[i]] = i
            for i in self.seq_list:
                name, sequence = i[0], i[1]
                code = [name]
                tmpCode = [0] * 400
                for j in range(len(sequence) - 2 + 1):
                    tmpCode[
                        AADict[sequence[j]] * 20 + AADict[sequence[j + 1]]
                    ] = (
                        tmpCode[
                            AADict[sequence[j]] * 20 + AADict[sequence[j + 1]]
                        ]
                        + 1
                    )
                if sum(tmpCode) != 0:
                    if normalized:
                        tmpCode = [i / sum(tmpCode) for i in tmpCode]
                code = code + tmpCode
                encodings.append(code)
            encodings = np.array(encodings)
            self.encodings = pd.DataFrame(
                encodings[1:, 1:].astype(float),
                columns=encodings[0, 1:],
                index=encodings[1:, 0],
            )
            return True
        except Exception as e:
            sys.exit(f"Error: {e}")

    def _DDE(self):
        try:
            myCodons = {
                "A": 4,
                "C": 2,
                "D": 2,
                "E": 2,
                "F": 2,
                "G": 4,
                "H": 2,
                "I": 3,
                "K": 2,
                "L": 6,
                "M": 1,
                "N": 2,
                "P": 4,
                "Q": 2,
                "R": 6,
                "S": 6,
                "T": 4,
                "V": 4,
                "W": 1,
                "Y": 2,
            }
            encodings = []
            diPeptides = [aa1 + aa2 for aa1 in self.AA for aa2 in self.AA]
            header = ["SampleName"] + ["DDE_{0}".format(i) for i in diPeptides]
            encodings.append(header)
            myTM = []
            for pair in diPeptides:
                myTM.append(
                    (myCodons[pair[0]] / 61) * (myCodons[pair[1]] / 61)
                )

            AADict = {}
            for i in range(len(self.AA)):
                AADict[self.AA[i]] = i
            for i in self.seq_list:
                name, sequence = i[0], i[1]
                code = [name]
                tmpCode = [0] * 400
                for j in range(len(sequence) - 2 + 1):
                    tmpCode[
                        AADict[sequence[j]] * 20 + AADict[sequence[j + 1]]
                    ] = (
                        tmpCode[
                            AADict[sequence[j]] * 20 + AADict[sequence[j + 1]]
                        ]
                        + 1
                    )
                if sum(tmpCode) != 0:
                    tmpCode = [i / sum(tmpCode) for i in tmpCode]

                myTV = []
                for j in range(len(myTM)):
                    myTV.append(myTM[j] * (1 - myTM[j]) / (len(sequence) - 1))
                for j in range(len(tmpCode)):
                    tmpCode[j] = (tmpCode[j] - myTM[j]) / math.sqrt(myTV[j])
                code = code + tmpCode
                encodings.append(code)
            encodings = np.array(encodings)
            self.encodings = pd.DataFrame(
                encodings[1:, 1:].astype(float),
                columns=encodings[0, 1:],
                index=encodings[1:, 0],
            )
            return True
        except Exception as e:
            sys.exit(f"Error: {e}")

    def _TPC(self, type=None):
        encodings = []
        triPeptides = [
            f"TPC{type}_" + aa1 + aa2 + aa3
            for aa1 in self.AA
            for aa2 in self.AA
            for aa3 in self.AA
        ]
        header = ["SampleName"] + triPeptides
        encodings.append(header)
        AADict = {}
        for i in range(len(self.AA)):
            AADict[self.AA[i]] = i
        for i in self.seq_list:
            try:
                name, sequence = i[0], i[1]
                if len(sequence) < 3:
                    print(f"Skipping {name} TPC sequence is less than 3")
                    continue
                code = [name]
                tmpCode = [0] * 8000
                for j in range(len(sequence) - 3 + 1):
                    tmpCode[
                        AADict[sequence[j]] * 400
                        + AADict[sequence[j + 1]] * 20
                        + AADict[sequence[j + 2]]
                    ] = (
                        tmpCode[
                            AADict[sequence[j]] * 400
                            + AADict[sequence[j + 1]] * 20
                            + AADict[sequence[j + 2]]
                        ]
                        + 1
                    )
                if sum(tmpCode) != 0:
                    if type == 1:
                        tmpCode = [i / sum(tmpCode) for i in tmpCode]
                    elif type == 3:
                        """
                        Type 3 is count of trimer per 1000 aa
                        """
                        tmpCode = [
                            round((i * (1000 / len(sequence))), 1)
                            for i in tmpCode
                        ]
                code = code + tmpCode
                encodings.append(code)
            except Exception as e:
                print(f"Error: {e} {name}")
        encodings = np.array(encodings)
        self.encodings = pd.DataFrame(
            encodings[1:, 1:].astype(float),
            columns=encodings[0, 1:],
            index=encodings[1:, 0],
        )
        return True

    def _GAAC(self):
        try:
            group = {
                "alphatic": "GAVLMI",
                "aromatic": "FYW",
                "postivecharge": "KRH",
                "negativecharge": "DE",
                "uncharge": "STCPNQ",
            }
            groupKey = group.keys()
            encodings = []
            header = ["SampleName"]
            for key in groupKey:
                header.append("GAAC_{0}".format(key))
            encodings.append(header)
            for i in self.seq_list:
                name, sequence = i[0], i[1]
                code = [name]
                count = Counter(sequence)
                myDict = {}
                for key in groupKey:
                    for aa in group[key]:
                        myDict[key] = myDict.get(key, 0) + count[aa]
                for key in groupKey:
                    code.append(myDict[key] / len(sequence))
                encodings.append(code)
            encodings = np.array(encodings)
            self.encodings = pd.DataFrame(
                encodings[1:, 1:].astype(float),
                columns=encodings[0, 1:],
                index=encodings[1:, 0],
            )
            return True
        except Exception as e:
            sys.exit(f"Error: {e}")

    def generateGroupPairs(self, groupKey):
        gPair = {}
        for key1 in groupKey:
            for key2 in groupKey:
                gPair[key1 + "." + key2] = 0
        return gPair

    def _CKSAAGP(self, normalized=True):
        try:
            gap = self.__default_para["kspace"]
            group = {
                "alphaticr": "GAVLMI",
                "aromatic": "FYW",
                "postivecharger": "KRH",
                "negativecharger": "DE",
                "uncharger": "STCPNQ",
            }
            groupKey = group.keys()
            index = {}
            for key in groupKey:
                for aa in group[key]:
                    index[aa] = key
            gPairIndex = []
            for key1 in groupKey:
                for key2 in groupKey:
                    gPairIndex.append(key1 + "." + key2)
            encodings = []
            header = ["SampleName"]
            for g in range(gap + 1):
                for p in gPairIndex:
                    header.append("CKSAAGP_" + p + ".gap" + str(g))
            encodings.append(header)
            for i in self.seq_list:
                name, sequence = i[0], i[1]
                code = [name]
                for g in range(gap + 1):
                    gPair = self.generateGroupPairs(groupKey)
                    sum = 0
                    for p1 in range(len(sequence)):
                        p2 = p1 + g + 1
                        if (
                            p2 < len(sequence)
                            and sequence[p1] in self.AA
                            and sequence[p2] in self.AA
                        ):
                            gPair[
                                index[sequence[p1]] + "." + index[sequence[p2]]
                            ] = (
                                gPair[
                                    index[sequence[p1]]
                                    + "."
                                    + index[sequence[p2]]
                                ]
                                + 1
                            )
                            sum = sum + 1
                    if sum == 0:
                        for gp in gPairIndex:
                            code.append(0)
                    else:
                        for gp in gPairIndex:
                            if normalized:
                                code.append(gPair[gp] / sum)
                            else:
                                code.append(gPair[gp])
                encodings.append(code)
            encodings = np.array(encodings)
            self.encodings = pd.DataFrame(
                encodings[1:, 1:].astype(float),
                columns=encodings[0, 1:],
                index=encodings[1:, 0],
            )
            return True
        except Exception as e:
            sys.exit(f"Error: {e}")

    def _GDPC(self, normalized=True):
        try:
            group = {
                "alphaticr": "GAVLMI",
                "aromatic": "FYW",
                "postivecharger": "KRH",
                "negativecharger": "DE",
                "uncharger": "STCPNQ",
            }
            groupKey = group.keys()
            dipeptide = [g1 + "." + g2 for g1 in groupKey for g2 in groupKey]
            index = {}
            for key in groupKey:
                for aa in group[key]:
                    index[aa] = key
            encodings = []
            header = ["SampleName"] + ["GDPC_{0}".format(i) for i in dipeptide]
            encodings.append(header)
            for i in self.seq_list:
                name, sequence = i[0], i[1]
                code = [name]
                myDict = {}
                for t in dipeptide:
                    myDict[t] = 0
                sum = 0
                for j in range(len(sequence) - 2 + 1):
                    myDict[
                        index[sequence[j]] + "." + index[sequence[j + 1]]
                    ] = (
                        myDict[
                            index[sequence[j]] + "." + index[sequence[j + 1]]
                        ]
                        + 1
                    )
                    sum = sum + 1
                if sum == 0:
                    for t in dipeptide:
                        code.append(0)
                else:
                    for t in dipeptide:
                        if normalized:
                            code.append(myDict[t] / sum)
                        else:
                            code.append(myDict[t])
                encodings.append(code)
            encodings = np.array(encodings)
            self.encodings = pd.DataFrame(
                encodings[1:, 1:].astype(float),
                columns=encodings[0, 1:],
                index=encodings[1:, 0],
            )
            return True
        except Exception as e:
            sys.exit(f"Error: {e}")

    def _GTPC(self, normalized=True):
        try:
            group = {
                "alphaticr": "GAVLMI",
                "aromatic": "FYW",
                "postivecharger": "KRH",
                "negativecharger": "DE",
                "uncharger": "STCPNQ",
            }
            groupKey = group.keys()
            triple = [
                g1 + "." + g2 + "." + g3
                for g1 in groupKey
                for g2 in groupKey
                for g3 in groupKey
            ]
            index = {}
            for key in groupKey:
                for aa in group[key]:
                    index[aa] = key
            encodings = []
            header = ["SampleName"] + ["GTPC_{0}".format(i) for i in triple]
            encodings.append(header)
            for i in self.seq_list:
                name, sequence = i[0], i[1]
                code = [name]
                myDict = {}
                for t in triple:
                    myDict[t] = 0
                sum = 0
                for j in range(len(sequence) - 3 + 1):
                    myDict[
                        index[sequence[j]]
                        + "."
                        + index[sequence[j + 1]]
                        + "."
                        + index[sequence[j + 2]]
                    ] = (
                        myDict[
                            index[sequence[j]]
                            + "."
                            + index[sequence[j + 1]]
                            + "."
                            + index[sequence[j + 2]]
                        ]
                        + 1
                    )
                    sum = sum + 1
                if sum == 0:
                    for t in triple:
                        code.append(0)
                else:
                    for t in triple:
                        if normalized:
                            code.append(myDict[t] / sum)
                        else:
                            code.append(myDict[t])
                encodings.append(code)
            encodings = np.array(encodings)
            self.encodings = pd.DataFrame(
                encodings[1:, 1:].astype(float),
                columns=encodings[0, 1:],
                index=encodings[1:, 0],
            )
            return True
        except Exception as e:
            sys.exit(f"Error: {e}")

    def _NMBroto(self):
        try:
            props = self.__default_para["aaindex"].split(";")
            nlag = self.__default_para["nlag"]
            fileAAidx = os.path.join(self.datadir, "AAidx.txt")
            with open(fileAAidx) as f:
                records = f.readlines()[1:]
            myDict = {}
            for i in records:
                array = i.rstrip().split("\t")
                myDict[array[0]] = array[1:]
            AAidx = []
            AAidxName = []
            for i in props:
                if i in myDict:
                    AAidx.append(myDict[i])
                    AAidxName.append(i)
                else:
                    print('"' + i + '" properties not exist.')
                    return None
            AAidx1 = np.array([float(j) for i in AAidx for j in i])
            AAidx = AAidx1.reshape((len(AAidx), 20))
            pstd = np.std(AAidx, axis=1)
            pmean = np.average(AAidx, axis=1)
            for i in range(len(AAidx)):
                for j in range(len(AAidx[i])):
                    AAidx[i][j] = (AAidx[i][j] - pmean[i]) / pstd[i]
            index = {}
            for i in range(len(self.AA)):
                index[self.AA[i]] = i
            encodings = []
            header = ["SampleName"]
            for p in props:
                for n in range(1, nlag + 1):
                    header.append("NMBroto_" + p + ".lag" + str(n))
            encodings.append(header)
            for i in self.seq_list:
                name, sequence = i[0], i[1]
                if len(sequence) <= nlag + 1:
                    if self.verbose:
                        print(
                            f"Skipping {name} NMBroto requires sequence length > nlag+1: {str(nlag + 1)}"
                        )
                        continue
                code = [name]
                N = len(sequence)
                for prop in range(len(props)):
                    for n in range(1, nlag + 1):
                        if len(sequence) > nlag:
                            # if key is '-', then the value is 0
                            rn = sum(
                                [
                                    AAidx[prop][index.get(sequence[j], 0)]
                                    * AAidx[prop][
                                        index.get(sequence[j + n], 0)
                                    ]
                                    for j in range(len(sequence) - n)
                                ]
                            ) / (N - n)
                        else:
                            rn = "NA"
                        code.append(rn)
                encodings.append(code)
            encodings = np.array(encodings)
            self.encodings = pd.DataFrame(
                encodings[1:, 1:].astype(float),
                columns=encodings[0, 1:],
                index=encodings[1:, 0],
            )
            return True
        except Exception as e:
            sys.exit(f"Error: {e}")

    def _Moran(self):
        try:
            props = self.__default_para["aaindex"].split(";")
            nlag = self.__default_para["nlag"]
            fileAAidx = os.path.join(self.datadir, "AAidx.txt")
            with open(fileAAidx) as f:
                records = f.readlines()[1:]
            myDict = {}
            for i in records:
                array = i.rstrip().split("\t")
                myDict[array[0]] = array[1:]
            AAidx = []
            AAidxName = []
            for i in props:
                if i in myDict:
                    AAidx.append(myDict[i])
                    AAidxName.append(i)
                else:
                    print('"' + i + '" properties not exist.')
                    return None
            AAidx1 = np.array([float(j) for i in AAidx for j in i])
            AAidx = AAidx1.reshape((len(AAidx), 20))
            propMean = np.mean(AAidx, axis=1)
            propStd = np.std(AAidx, axis=1)
            for i in range(len(AAidx)):
                for j in range(len(AAidx[i])):
                    AAidx[i][j] = (AAidx[i][j] - propMean[i]) / propStd[i]
            index = {}
            for i in range(len(self.AA)):
                index[self.AA[i]] = i
            encodings = []
            header = ["SampleName"]
            for p in props:
                for n in range(1, nlag + 1):
                    header.append("Moran_" + p + ".lag" + str(n))
            encodings.append(header)
            for i in self.seq_list:
                name, sequence = i[0], i[1]
                if len(sequence) <= nlag + 1:
                    if self.verbose:
                        print(
                            f"Skipping {name} Moran requires sequence length > nlag+1: {str(nlag + 1)}"
                        )
                        continue
                code = [name]
                N = len(sequence)
                for prop in range(len(props)):
                    xmean = (
                        sum([AAidx[prop][index[aa]] for aa in sequence]) / N
                    )
                    for n in range(1, nlag + 1):
                        if len(sequence) > nlag:
                            # if key is '-', then the value is 0
                            fenzi = sum(
                                [
                                    (
                                        AAidx[prop][index.get(sequence[j], 0)]
                                        - xmean
                                    )
                                    * (
                                        AAidx[prop][
                                            index.get(sequence[j + n], 0)
                                        ]
                                        - xmean
                                    )
                                    for j in range(len(sequence) - n)
                                ]
                            ) / (N - n)
                            fenmu = (
                                sum(
                                    [
                                        (
                                            AAidx[prop][
                                                index.get(sequence[j], 0)
                                            ]
                                            - xmean
                                        )
                                        ** 2
                                        for j in range(len(sequence))
                                    ]
                                )
                                / N
                            )
                            rn = fenzi / fenmu
                        else:
                            rn = "NA"
                        code.append(rn)
                encodings.append(code)
            encodings = np.array(encodings)
            self.encodings = pd.DataFrame(
                encodings[1:, 1:].astype(float),
                columns=encodings[0, 1:],
                index=encodings[1:, 0],
            )
            return True
        except Exception as e:
            sys.exit(f"Error: {e}")

    def _Geary(self):
        try:
            props = self.__default_para["aaindex"].split(";")
            nlag = self.__default_para["nlag"]
            fileAAidx = os.path.join(
                self.datadir,
                "AAidx.txt",
            )
            with open(fileAAidx) as f:
                records = f.readlines()[1:]
            myDict = {}
            for i in records:
                array = i.rstrip().split("\t")
                myDict[array[0]] = array[1:]
            AAidx = []
            AAidxName = []
            for i in props:
                if i in myDict:
                    AAidx.append(myDict[i])
                    AAidxName.append(i)
                else:
                    print('"' + i + '" properties not exist.')
                    return None
            AAidx1 = np.array([float(j) for i in AAidx for j in i])
            AAidx = AAidx1.reshape((len(AAidx), 20))
            propMean = np.mean(AAidx, axis=1)
            propStd = np.std(AAidx, axis=1)
            for i in range(len(AAidx)):
                for j in range(len(AAidx[i])):
                    AAidx[i][j] = (AAidx[i][j] - propMean[i]) / propStd[i]
            index = {}
            for i in range(len(self.AA)):
                index[self.AA[i]] = i
            encodings = []
            header = ["SampleName"]
            for p in props:
                for n in range(1, nlag + 1):
                    header.append("Geary_" + p + ".lag" + str(n))
            encodings.append(header)
            for i in self.seq_list:
                name, sequence = i[0], i[1]
                code = [name]
                N = len(sequence)
                for prop in range(len(props)):
                    xmean = (
                        sum([AAidx[prop][index[aa]] for aa in sequence]) / N
                    )
                    for n in range(1, nlag + 1):
                        if len(sequence) > nlag:
                            # if key is '-', then the value is 0
                            rn = (
                                (N - 1)
                                / (2 * (N - n))
                                * (
                                    (
                                        sum(
                                            [
                                                (
                                                    AAidx[prop][
                                                        index.get(
                                                            sequence[j], 0
                                                        )
                                                    ]
                                                    - AAidx[prop][
                                                        index.get(
                                                            sequence[j + n], 0
                                                        )
                                                    ]
                                                )
                                                ** 2
                                                for j in range(
                                                    len(sequence) - n
                                                )
                                            ]
                                        )
                                    )
                                    / (
                                        sum(
                                            [
                                                (
                                                    AAidx[prop][
                                                        index.get(
                                                            sequence[j], 0
                                                        )
                                                    ]
                                                    - xmean
                                                )
                                                ** 2
                                                for j in range(len(sequence))
                                            ]
                                        )
                                    )
                                )
                            )
                        else:
                            rn = "NA"
                        code.append(rn)
                encodings.append(code)
            encodings = np.array(encodings)
            self.encodings = pd.DataFrame(
                encodings[1:, 1:].astype(float),
                columns=encodings[0, 1:],
                index=encodings[1:, 0],
            )
            return True
        except Exception as e:
            sys.exit(f"Error: {e}")

    def generatePropertyPairs(self, myPropertyName):
        pairs = []
        for i in range(len(myPropertyName)):
            for j in range(i + 1, len(myPropertyName)):
                pairs.append([myPropertyName[i], myPropertyName[j]])
                pairs.append([myPropertyName[j], myPropertyName[i]])
        return pairs

    def _AC(self):
        try:
            property_name = self.__default_para["aaindex"].split(";")
            nlag = self.__default_para["nlag"]
            try:
                data_file = os.path.join(self.datadir, "AAindex.data")
                with open(data_file, "rb") as handle:
                    property_dict = pickle.load(handle)
            except Exception as e:
                sys.exit(
                    f"Could not find the physicochemical properties file: {e}"
                )
                return False
            for p_name in property_name:
                tmp = np.array(property_dict[p_name], dtype=float)
                pmean = np.average(tmp)
                pstd = np.std(tmp)
                property_dict[p_name] = [(elem - pmean) / pstd for elem in tmp]
            AA_order_dict = {}
            for i in range(len(self.AA)):
                AA_order_dict[self.AA[i]] = i
            encodings = []
            header = ["SampleName"]
            for p_name in property_name:
                for i in range(nlag):
                    header.append("AC_%s.lag%s" % (p_name, i + 1))
            encodings.append(header)
            for i in self.seq_list:
                name, sequence = i[0], i[1]
                if len(sequence) < nlag + 1:
                    if self.verbose:
                        f"Skipping {name} AC requires sequence length > nlag + 1: {str(nlag + 1)}"
                    continue
                code = [name]
                L = len(sequence)
                for p_name in property_name:
                    xmean = (
                        sum(
                            [
                                property_dict[p_name][AA_order_dict[aa]]
                                for aa in sequence
                            ]
                        )
                        / L
                    )
                    for lag in range(1, nlag + 1):
                        ac = 0
                        try:
                            ac = sum(
                                [
                                    (
                                        property_dict[p_name][
                                            AA_order_dict[sequence[j]]
                                        ]
                                        - xmean
                                    )
                                    * (
                                        property_dict[p_name][
                                            AA_order_dict[sequence[j + lag]]
                                        ]
                                        - xmean
                                    )
                                    for j in range(L - lag)
                                ]
                            ) / (L - lag)
                        except Exception as e:
                            ac = 0
                            print(f"Error: {e}")
                        code.append(ac)
                encodings.append(code)
            encodings = np.array(encodings)
            self.encodings = pd.DataFrame(
                encodings[1:, 1:].astype(float),
                columns=encodings[0, 1:],
                index=encodings[1:, 0],
            )
            return True
        except Exception as e:
            sys.exit(f"Error: {e}")

    def _CC(self):
        try:
            property_name = self.__default_para["aaindex"].split(";")
            if len(property_name) < 2:
                self.error_msg = "More than two property should be selected for this descriptor."
                return False
            nlag = self.__default_para["nlag"]
            try:
                data_file = os.path.join(self.datadir, "AAindex.data")
                with open(data_file, "rb") as handle:
                    property_dict = pickle.load(handle)
            except Exception as e:
                sys.exit(
                    f"Could not find the physicochemical properties file: {e}"
                )
                return False
            for p_name in property_name:
                tmp = np.array(property_dict[p_name], dtype=float)
                pmean = np.average(tmp)
                pstd = np.std(tmp)
                property_dict[p_name] = [(elem - pmean) / pstd for elem in tmp]
            AA_order_dict = {}
            for i in range(len(self.AA)):
                AA_order_dict[self.AA[i]] = i
            property_pairs = self.generatePropertyPairs(property_name)
            encodings = []
            header = ["SampleName"]
            header += [
                "CC_" + p[0] + "_" + p[1] + "_lag." + str(lag)
                for p in property_pairs
                for lag in range(1, nlag + 1)
            ]
            encodings.append(header)
            for i in self.seq_list:
                name, sequence = i[0], i[1]
                if len(sequence) < nlag + 1:
                    if self.verbose:
                        f"Skipping {name} CC requires sequence length > nlag + 1: {str(nlag + 1)}"
                    continue
                code = [name]
                L = len(sequence)
                for pair in property_pairs:
                    mean_p1 = (
                        sum(
                            [
                                property_dict[pair[0]][AA_order_dict[aa]]
                                for aa in sequence
                            ]
                        )
                        / L
                    )
                    mean_p2 = (
                        sum(
                            [
                                property_dict[pair[1]][AA_order_dict[aa]]
                                for aa in sequence
                            ]
                        )
                        / L
                    )
                    for lag in range(1, nlag + 1):
                        cc = 0
                        try:
                            cc = sum(
                                [
                                    (
                                        property_dict[pair[0]][
                                            AA_order_dict[sequence[j]]
                                        ]
                                        - mean_p1
                                    )
                                    * (
                                        property_dict[pair[1]][
                                            AA_order_dict[sequence[j + lag]]
                                        ]
                                        - mean_p2
                                    )
                                    for j in range(L - lag)
                                ]
                            ) / (L - lag)
                        except Exception as e:
                            cc = 0
                            print(f"Error: {e}")
                        code.append(cc)
                encodings.append(code)
            encodings = np.array(encodings)
            self.encodings = pd.DataFrame(
                encodings[1:, 1:].astype(float),
                columns=encodings[0, 1:],
                index=encodings[1:, 0],
            )
            return True
        except Exception as e:
            sys.exit(f"Error: {e}")

    def _ACC(self):
        try:
            property_name = self.__default_para["aaindex"].split(";")
            if len(property_name) < 2:
                self.error_msg = "More than two property should be selected for this descriptor."
                return False
            nlag = self.__default_para["nlag"]
            try:
                data_file = os.path.join(self.datadir, "AAindex.data")
                with open(data_file, "rb") as handle:
                    property_dict = pickle.load(handle)
            except Exception as e:
                sys.exit(
                    f"Could not find the physicochemical properties file: {e}"
                )
                return False
            for p_name in property_name:
                tmp = np.array(property_dict[p_name], dtype=float)
                pmean = np.average(tmp)
                pstd = np.std(tmp)
                property_dict[p_name] = [(elem - pmean) / pstd for elem in tmp]
            AA_order_dict = {}
            for i in range(len(self.AA)):
                AA_order_dict[self.AA[i]] = i
            property_pairs = self.generatePropertyPairs(property_name)
            encodings = []
            header = ["SampleName"]
            for p_name in property_name:
                for i in range(nlag):
                    header.append("ACC_%s.lag%s" % (p_name, i + 1))
            header += [
                "ACC_" + p[0] + "_" + p[1] + "_lag." + str(lag)
                for p in property_pairs
                for lag in range(1, nlag + 1)
            ]
            encodings.append(header)
            for i in self.seq_list:
                name, sequence = i[0], i[1]
                if len(sequence) < nlag + 1:
                    if self.verbose:
                        print(
                            f"Skipping {name} ACC requires sequence length > nlag + 1: {str(nlag + 1)}"
                        )
                        continue
                code = [name]
                L = len(sequence)
                for p_name in property_name:
                    xmean = (
                        sum(
                            [
                                property_dict[p_name][AA_order_dict[aa]]
                                for aa in sequence
                            ]
                        )
                        / L
                    )
                    for lag in range(1, nlag + 1):
                        ac = 0
                        try:
                            ac = sum(
                                [
                                    (
                                        property_dict[p_name][
                                            AA_order_dict[sequence[j]]
                                        ]
                                        - xmean
                                    )
                                    * (
                                        property_dict[p_name][
                                            AA_order_dict[sequence[j + lag]]
                                        ]
                                        - xmean
                                    )
                                    for j in range(L - lag)
                                ]
                            ) / (L - lag)
                        except Exception as e:
                            print(f"Error: {e}")
                            ac = 0
                        code.append(ac)
                for pair in property_pairs:
                    mean_p1 = (
                        sum(
                            [
                                property_dict[pair[0]][AA_order_dict[aa]]
                                for aa in sequence
                            ]
                        )
                        / L
                    )
                    mean_p2 = (
                        sum(
                            [
                                property_dict[pair[1]][AA_order_dict[aa]]
                                for aa in sequence
                            ]
                        )
                        / L
                    )
                    for lag in range(1, nlag + 1):
                        cc = 0
                        try:
                            cc = sum(
                                [
                                    (
                                        property_dict[pair[0]][
                                            AA_order_dict[sequence[j]]
                                        ]
                                        - mean_p1
                                    )
                                    * (
                                        property_dict[pair[1]][
                                            AA_order_dict[sequence[j + lag]]
                                        ]
                                        - mean_p2
                                    )
                                    for j in range(L - lag)
                                ]
                            ) / (L - lag)
                        except Exception as e:
                            print(f"Error: {e}")
                            cc = 0
                        code.append(cc)
                encodings.append(code)
            encodings = np.array(encodings)
            self.encodings = pd.DataFrame(
                encodings[1:, 1:].astype(float),
                columns=encodings[0, 1:],
                index=encodings[1:, 0],
            )
            return True
        except Exception as e:
            sys.exit(f"Error: {e}")

    def Count(self, seq1, seq2):
        sum = 0
        for aa in seq1:
            sum = sum + seq2.count(aa)
        return sum

    def _CTDC(self):
        try:
            group1 = {
                "hydrophobicity_PRAM900101": "RKEDQN",
                "hydrophobicity_ARGP820101": "QSTNGDE",
                "hydrophobicity_ZIMJ680101": "QNGSWTDERA",
                "hydrophobicity_PONP930101": "KPDESNQT",
                "hydrophobicity_CASG920101": "KDEQPSRNTG",
                "hydrophobicity_ENGD860101": "RDKENQHYP",
                "hydrophobicity_FASG890101": "KERSQD",
                "normwaalsvolume": "GASTPDC",
                "polarity": "LIFWCMVY",
                "polarizability": "GASDT",
                "charge": "KR",
                "secondarystruct": "EALMQKRH",
                "solventaccess": "ALFCGIVW",
            }
            group2 = {
                "hydrophobicity_PRAM900101": "GASTPHY",
                "hydrophobicity_ARGP820101": "RAHCKMV",
                "hydrophobicity_ZIMJ680101": "HMCKV",
                "hydrophobicity_PONP930101": "GRHA",
                "hydrophobicity_CASG920101": "AHYMLV",
                "hydrophobicity_ENGD860101": "SGTAW",
                "hydrophobicity_FASG890101": "NTPG",
                "normwaalsvolume": "NVEQIL",
                "polarity": "PATGS",
                "polarizability": "CPNVEQIL",
                "charge": "ANCQGHILMFPSTWYV",
                "secondarystruct": "VIYCWFT",
                "solventaccess": "RKQEND",
            }
            group3 = {
                "hydrophobicity_PRAM900101": "CLVIMFW",
                "hydrophobicity_ARGP820101": "LYPFIW",
                "hydrophobicity_ZIMJ680101": "LPFYI",
                "hydrophobicity_PONP930101": "YMFWLCVI",
                "hydrophobicity_CASG920101": "FIWC",
                "hydrophobicity_ENGD860101": "CVLIMF",
                "hydrophobicity_FASG890101": "AYHWVMFLIC",
                "normwaalsvolume": "MHKFRYW",
                "polarity": "HQRKNED",
                "polarizability": "KMHFRYW",
                "charge": "DE",
                "secondarystruct": "GNPSD",
                "solventaccess": "MSPTHY",
            }
            groups = [group1, group2, group3]
            property = (
                "hydrophobicity_PRAM900101",
                "hydrophobicity_ARGP820101",
                "hydrophobicity_ZIMJ680101",
                "hydrophobicity_PONP930101",
                "hydrophobicity_CASG920101",
                "hydrophobicity_ENGD860101",
                "hydrophobicity_FASG890101",
                "normwaalsvolume",
                "polarity",
                "polarizability",
                "charge",
                "secondarystruct",
                "solventaccess",
            )
            encodings = []
            header = ["SampleName"]
            for p in property:
                for g in range(1, len(groups) + 1):
                    header.append("CTDC_" + p + ".G" + str(g))
            encodings.append(header)
            for i in self.seq_list:
                name, sequence = i[0], i[1]
                code = [name]
                for p in property:
                    c1 = self.Count(group1[p], sequence) / len(sequence)
                    c2 = self.Count(group2[p], sequence) / len(sequence)
                    c3 = 1 - c1 - c2
                    code = code + [c1, c2, c3]
                encodings.append(code)
            encodings = np.array(encodings)
            self.encodings = pd.DataFrame(
                encodings[1:, 1:].astype(float),
                columns=encodings[0, 1:],
                index=encodings[1:, 0],
            )
            return True
        except Exception as e:
            sys.exit(f"Error: {e}")

    def _CTDT(self):
        try:
            group1 = {
                "hydrophobicity_PRAM900101": "RKEDQN",
                "hydrophobicity_ARGP820101": "QSTNGDE",
                "hydrophobicity_ZIMJ680101": "QNGSWTDERA",
                "hydrophobicity_PONP930101": "KPDESNQT",
                "hydrophobicity_CASG920101": "KDEQPSRNTG",
                "hydrophobicity_ENGD860101": "RDKENQHYP",
                "hydrophobicity_FASG890101": "KERSQD",
                "normwaalsvolume": "GASTPDC",
                "polarity": "LIFWCMVY",
                "polarizability": "GASDT",
                "charge": "KR",
                "secondarystruct": "EALMQKRH",
                "solventaccess": "ALFCGIVW",
            }
            group2 = {
                "hydrophobicity_PRAM900101": "GASTPHY",
                "hydrophobicity_ARGP820101": "RAHCKMV",
                "hydrophobicity_ZIMJ680101": "HMCKV",
                "hydrophobicity_PONP930101": "GRHA",
                "hydrophobicity_CASG920101": "AHYMLV",
                "hydrophobicity_ENGD860101": "SGTAW",
                "hydrophobicity_FASG890101": "NTPG",
                "normwaalsvolume": "NVEQIL",
                "polarity": "PATGS",
                "polarizability": "CPNVEQIL",
                "charge": "ANCQGHILMFPSTWYV",
                "secondarystruct": "VIYCWFT",
                "solventaccess": "RKQEND",
            }
            group3 = {
                "hydrophobicity_PRAM900101": "CLVIMFW",
                "hydrophobicity_ARGP820101": "LYPFIW",
                "hydrophobicity_ZIMJ680101": "LPFYI",
                "hydrophobicity_PONP930101": "YMFWLCVI",
                "hydrophobicity_CASG920101": "FIWC",
                "hydrophobicity_ENGD860101": "CVLIMF",
                "hydrophobicity_FASG890101": "AYHWVMFLIC",
                "normwaalsvolume": "MHKFRYW",
                "polarity": "HQRKNED",
                "polarizability": "KMHFRYW",
                "charge": "DE",
                "secondarystruct": "GNPSD",
                "solventaccess": "MSPTHY",
            }
            # groups = [group1, group2, group3]
            property = (
                "hydrophobicity_PRAM900101",
                "hydrophobicity_ARGP820101",
                "hydrophobicity_ZIMJ680101",
                "hydrophobicity_PONP930101",
                "hydrophobicity_CASG920101",
                "hydrophobicity_ENGD860101",
                "hydrophobicity_FASG890101",
                "normwaalsvolume",
                "polarity",
                "polarizability",
                "charge",
                "secondarystruct",
                "solventaccess",
            )
            encodings = []
            header = ["SampleName"]
            for p in property:
                for tr in ("Tr1221", "Tr1331", "Tr2332"):
                    header.append("CTDT_" + p + "." + tr)
            encodings.append(header)
            for i in self.seq_list:
                name, sequence = i[0], i[1]
                code = [name]
                aaPair = [
                    sequence[j : j + 2] for j in range(len(sequence) - 1)
                ]
                for p in property:
                    c1221, c1331, c2332 = 0, 0, 0
                    for pair in aaPair:
                        if (pair[0] in group1[p] and pair[1] in group2[p]) or (
                            pair[0] in group2[p] and pair[1] in group1[p]
                        ):
                            c1221 = c1221 + 1
                            continue
                        if (pair[0] in group1[p] and pair[1] in group3[p]) or (
                            pair[0] in group3[p] and pair[1] in group1[p]
                        ):
                            c1331 = c1331 + 1
                            continue
                        if (pair[0] in group2[p] and pair[1] in group3[p]) or (
                            pair[0] in group3[p] and pair[1] in group2[p]
                        ):
                            c2332 = c2332 + 1
                    code = code + [
                        c1221 / len(aaPair),
                        c1331 / len(aaPair),
                        c2332 / len(aaPair),
                    ]
                encodings.append(code)
            encodings = np.array(encodings)
            self.encodings = pd.DataFrame(
                encodings[1:, 1:].astype(float),
                columns=encodings[0, 1:],
                index=encodings[1:, 0],
            )
            return True
        except Exception as e:
            sys.exit(f"Error: {e}")

    def Count1(self, aaSet, sequence):
        number = 0
        for aa in sequence:
            if aa in aaSet:
                number = number + 1
        cutoffNums = [
            1,
            math.floor(0.25 * number),
            math.floor(0.50 * number),
            math.floor(0.75 * number),
            number,
        ]
        cutoffNums = [i if i >= 1 else 1 for i in cutoffNums]

        code = []
        for cutoff in cutoffNums:
            myCount = 0
            for i in range(len(sequence)):
                if sequence[i] in aaSet:
                    myCount += 1
                    if myCount == cutoff:
                        code.append((i + 1) / len(sequence) * 100)
                        break
            if myCount == 0:
                code.append(0)
        return code

    def _CTDD(self):
        try:
            group1 = {
                "hydrophobicity_PRAM900101": "RKEDQN",
                "hydrophobicity_ARGP820101": "QSTNGDE",
                "hydrophobicity_ZIMJ680101": "QNGSWTDERA",
                "hydrophobicity_PONP930101": "KPDESNQT",
                "hydrophobicity_CASG920101": "KDEQPSRNTG",
                "hydrophobicity_ENGD860101": "RDKENQHYP",
                "hydrophobicity_FASG890101": "KERSQD",
                "normwaalsvolume": "GASTPDC",
                "polarity": "LIFWCMVY",
                "polarizability": "GASDT",
                "charge": "KR",
                "secondarystruct": "EALMQKRH",
                "solventaccess": "ALFCGIVW",
            }
            group2 = {
                "hydrophobicity_PRAM900101": "GASTPHY",
                "hydrophobicity_ARGP820101": "RAHCKMV",
                "hydrophobicity_ZIMJ680101": "HMCKV",
                "hydrophobicity_PONP930101": "GRHA",
                "hydrophobicity_CASG920101": "AHYMLV",
                "hydrophobicity_ENGD860101": "SGTAW",
                "hydrophobicity_FASG890101": "NTPG",
                "normwaalsvolume": "NVEQIL",
                "polarity": "PATGS",
                "polarizability": "CPNVEQIL",
                "charge": "ANCQGHILMFPSTWYV",
                "secondarystruct": "VIYCWFT",
                "solventaccess": "RKQEND",
            }
            group3 = {
                "hydrophobicity_PRAM900101": "CLVIMFW",
                "hydrophobicity_ARGP820101": "LYPFIW",
                "hydrophobicity_ZIMJ680101": "LPFYI",
                "hydrophobicity_PONP930101": "YMFWLCVI",
                "hydrophobicity_CASG920101": "FIWC",
                "hydrophobicity_ENGD860101": "CVLIMF",
                "hydrophobicity_FASG890101": "AYHWVMFLIC",
                "normwaalsvolume": "MHKFRYW",
                "polarity": "HQRKNED",
                "polarizability": "KMHFRYW",
                "charge": "DE",
                "secondarystruct": "GNPSD",
                "solventaccess": "MSPTHY",
            }
            # groups = [group1, group2, group3]
            property = (
                "hydrophobicity_PRAM900101",
                "hydrophobicity_ARGP820101",
                "hydrophobicity_ZIMJ680101",
                "hydrophobicity_PONP930101",
                "hydrophobicity_CASG920101",
                "hydrophobicity_ENGD860101",
                "hydrophobicity_FASG890101",
                "normwaalsvolume",
                "polarity",
                "polarizability",
                "charge",
                "secondarystruct",
                "solventaccess",
            )
            encodings = []
            header = ["SampleName"]
            for p in property:
                for g in ("1", "2", "3"):
                    for d in ["0", "25", "50", "75", "100"]:
                        header.append("CTDD_" + p + "." + g + ".residue" + d)
            encodings.append(header)
            for i in self.seq_list:
                name, sequence = i[0], i[1]
                code = [name]
                for p in property:
                    code = (
                        code
                        + self.Count1(group1[p], sequence)
                        + self.Count1(group2[p], sequence)
                        + self.Count1(group3[p], sequence)
                    )
                encodings.append(code)
            encodings = np.array(encodings)
            self.encodings = pd.DataFrame(
                encodings[1:, 1:].astype(float),
                columns=encodings[0, 1:],
                index=encodings[1:, 0],
            )
            return True
        except Exception as e:
            sys.exit(f"Error: {e}")

    def CalculateKSCTriad(self, sequence, gap, features, AADict):
        res = []
        for g in range(gap + 1):
            myDict = {}
            for f in features:
                myDict[f] = 0

            for i in range(len(sequence)):
                if i + g + 1 < len(sequence) and i + 2 * g + 2 < len(sequence):
                    fea = (
                        AADict[sequence[i]]
                        + "."
                        + AADict[sequence[i + g + 1]]
                        + "."
                        + AADict[sequence[i + 2 * g + 2]]
                    )
                    myDict[fea] = myDict[fea] + 1

            maxValue, minValue = np.max(list(myDict.values())), np.min(
                list(myDict.values())
            )
            for f in features:
                res.append((myDict[f] - minValue) / maxValue)
        return res

    def _CTriad(self):
        AAGroup = {
            "g1": "AGV",
            "g2": "ILFP",
            "g3": "YMTS",
            "g4": "HNQW",
            "g5": "RK",
            "g6": "DE",
            "g7": "C",
        }
        myGroups = sorted(AAGroup.keys())
        AADict = {}
        for g in myGroups:
            for aa in AAGroup[g]:
                AADict[aa] = g
        features = [
            f1 + "." + f2 + "." + f3
            for f1 in myGroups
            for f2 in myGroups
            for f3 in myGroups
        ]
        encodings = []
        header = ["SampleName"]
        for f in features:
            header.append("CTriad_{0}".format(f))
        encodings.append(header)
        for i in self.seq_list:
            name, sequence = i[0], i[1]
            if len(sequence) < 30:
                if self.verbose:
                    print(
                        f"Skipping {name} CTriad needs sequence with length > 30 {name}"
                    )
                continue
            code = [name]
            code = code + self.CalculateKSCTriad(sequence, 0, features, AADict)
            encodings.append(code)
        encodings = np.array(encodings)
        self.encodings = pd.DataFrame(
            encodings[1:, 1:].astype(float),
            columns=encodings[0, 1:],
            index=encodings[1:, 0],
        )
        return True

    def _KSCTriad(self):
        try:
            gap = self.__default_para["kspace"]
            AAGroup = {
                "g1": "AGV",
                "g2": "ILFP",
                "g3": "YMTS",
                "g4": "HNQW",
                "g5": "RK",
                "g6": "DE",
                "g7": "C",
            }
            myGroups = sorted(AAGroup.keys())
            AADict = {}
            for g in myGroups:
                for aa in AAGroup[g]:
                    AADict[aa] = g
            features = [
                f1 + "." + f2 + "." + f3
                for f1 in myGroups
                for f2 in myGroups
                for f3 in myGroups
            ]
            encodings = []
            header = ["SampleName"]
            for g in range(gap + 1):
                for f in features:
                    header.append("KSCTriad_" + f + ".gap" + str(g))
            encodings.append(header)
            for i in self.seq_list:
                name, sequence = i[0], i[1]
                code = [name]
                if len(sequence) < 2 * gap + 3:
                    if self.verbose:
                        print(
                            f"Skipping {name} KSCTriand sequence should be > (2*gap+3)"
                        )
                    continue
                code = code + self.CalculateKSCTriad(
                    sequence, gap, features, AADict
                )
                encodings.append(code)
            encodings = np.array(encodings)
            self.encodings = pd.DataFrame(
                encodings[1:, 1:].astype(float),
                columns=encodings[0, 1:],
                index=encodings[1:, 0],
            )
            return True
        except Exception as e:
            sys.exit(f"Error: {e}")

    def _SOCNumber(self):
        try:
            nlag = self.__default_para["nlag"]
            dataFile = os.path.join(self.datadir, "Schneider-Wrede.txt")
            dataFile1 = os.path.join(self.datadir, "Grantham.txt")
            DictAA = {}
            for i in range(len(self.AA)):
                DictAA[self.AA[i]] = i
            DictAA1 = {}
            for i in range(len(self.AA)):
                DictAA1[self.AA[i]] = i
            with open(dataFile) as f:
                records = f.readlines()[1:]
            AADistance = []
            for i in records:
                array = i.rstrip().split()[1:] if i.rstrip() != "" else None
                AADistance.append(array)
            AADistance = np.array(
                [
                    float(AADistance[i][j])
                    for i in range(len(AADistance))
                    for j in range(len(AADistance[i]))
                ]
            ).reshape((20, 20))
            with open(dataFile1) as f:
                records = f.readlines()[1:]
            AADistance1 = []
            for i in records:
                array = i.rstrip().split()[1:] if i.rstrip() != "" else None
                AADistance1.append(array)
            AADistance1 = np.array(
                [
                    float(AADistance1[i][j])
                    for i in range(len(AADistance1))
                    for j in range(len(AADistance1[i]))
                ]
            ).reshape((20, 20))
            encodings = []
            header = ["SampleName"]
            for n in range(1, nlag + 1):
                header.append("SOCNumber_Schneider.lag" + str(n))
            for n in range(1, nlag + 1):
                header.append("SOCNumber_gGrantham.lag" + str(n))
            encodings.append(header)
            for i in self.seq_list:
                name, sequence = i[0], i[1]
                code = [name]
                for n in range(1, nlag + 1):
                    code.append(
                        sum(
                            [
                                AADistance[DictAA[sequence[j]]][
                                    DictAA[sequence[j + n]]
                                ]
                                ** 2
                                for j in range(len(sequence) - n)
                            ]
                        )
                        / (len(sequence) - n)
                    )
                for n in range(1, nlag + 1):
                    code.append(
                        sum(
                            [
                                AADistance1[DictAA1[sequence[j]]][
                                    DictAA1[sequence[j + n]]
                                ]
                                ** 2
                                for j in range(len(sequence) - n)
                            ]
                        )
                        / (len(sequence) - n)
                    )
                encodings.append(code)
            encodings = np.array(encodings)
            self.encodings = pd.DataFrame(
                encodings[1:, 1:].astype(float),
                columns=encodings[0, 1:],
                index=encodings[1:, 0],
            )
            return True
        except Exception as e:
            sys.exit(f"Error: {e}")

    def _QSOrder(self):
        nlag = self.__default_para["nlag"]
        w = self.__default_para["weight"]
        dataFile0 = os.path.join(self.datadir, "Schneider-Wrede.txt")
        dataFile1 = os.path.join(self.datadir, "Grantham.txt")
        DictAA = {}
        for i in range(len(self.AA)):
            DictAA[self.AA[i]] = i
        DictAA1 = {}
        for i in range(len(self.AA)):
            DictAA1[self.AA[i]] = i
        with open(dataFile0) as f:
            records = f.readlines()[1:]
        AADistance = []
        for i in records:
            array = i.rstrip().split()[1:] if i.rstrip() != "" else None
            AADistance.append(array)
        AADistance = np.array(
            [
                float(AADistance[i][j])
                for i in range(len(AADistance))
                for j in range(len(AADistance[i]))
            ]
        ).reshape((20, 20))
        with open(dataFile1) as f:
            records = f.readlines()[1:]
        AADistance1 = []
        for i in records:
            array = i.rstrip().split()[1:] if i.rstrip() != "" else None
            AADistance1.append(array)
        AADistance1 = np.array(
            [
                float(AADistance1[i][j])
                for i in range(len(AADistance1))
                for j in range(len(AADistance1[i]))
            ]
        ).reshape((20, 20))
        encodings = []
        header = ["SampleName"]
        for aa in self.AA:
            header.append("QSOrder_Schneider.Xr." + aa)
        for aa in self.AA:
            header.append("QSOrder_Grantham.Xr." + aa)
        for n in range(1, nlag + 1):
            header.append("QSOrder_Schneider.Xd." + str(n))
        for n in range(1, nlag + 1):
            header.append("QSOrder_Grantham.Xd." + str(n))
        encodings.append(header)
        try:
            for i in self.seq_list:
                name, sequence = i[0], i[1]
                if nlag > len(sequence) - 1:
                    print(
                        f"Skipping {name} QSOrder requires nlag > sequence length - 1: {str(nlag)}"
                    )
                    continue
                code = [name]
                arraySW = []
                arrayGM = []
                for n in range(1, nlag + 1):
                    arraySW.append(
                        sum(
                            [
                                AADistance[DictAA[sequence[j]]][
                                    DictAA[sequence[j + n]]
                                ]
                                ** 2
                                for j in range(len(sequence) - n)
                            ]
                        )
                    )
                    arrayGM.append(
                        sum(
                            [
                                AADistance1[DictAA1[sequence[j]]][
                                    DictAA1[sequence[j + n]]
                                ]
                                ** 2
                                for j in range(len(sequence) - n)
                            ]
                        )
                    )
                myDict = {}
                for aa in self.AA:
                    myDict[aa] = sequence.count(aa)
                for aa in self.AA:
                    code.append(myDict[aa] / (1 + w * sum(arraySW)))
                for aa in self.AA:
                    code.append(myDict[aa] / (1 + w * sum(arrayGM)))
                for num in arraySW:
                    code.append((w * num) / (1 + w * sum(arraySW)))
                for num in arrayGM:
                    code.append((w * num) / (1 + w * sum(arrayGM)))
                encodings.append(code)
        except Exception as e:
            print(f"Error: {e} {name}")

        encodings = np.array(encodings)
        self.encodings = pd.DataFrame(
            encodings[1:, 1:].astype(float),
            columns=encodings[0, 1:],
            index=encodings[1:, 0],
        )
        return True

    def Rvalue(self, aa1, aa2, AADict, Matrix):
        return sum(
            [
                (Matrix[i][AADict[aa1]] - Matrix[i][AADict[aa2]]) ** 2
                for i in range(len(Matrix))
            ]
        ) / len(Matrix)

    def _PAAC(self):
        try:
            lambdaValue = self.__default_para["lambdaValue"]
            w = self.__default_para["weight"]
            dataFile = os.path.join(self.datadir, "PAAC.txt")
            with open(dataFile) as f:
                records = f.readlines()
            AA = "".join(records[0].rstrip().split()[1:])
            AADict = {}
            for i in range(len(AA)):
                AADict[AA[i]] = i
            AAProperty = []
            AAPropertyNames = []
            for i in range(1, len(records)):
                array = (
                    records[i].rstrip().split()
                    if records[i].rstrip() != ""
                    else None
                )
                AAProperty.append([float(j) for j in array[1:]])
                AAPropertyNames.append(array[0])
            AAProperty1 = []
            for i in AAProperty:
                meanI = sum(i) / 20
                fenmu = math.sqrt(sum([(j - meanI) ** 2 for j in i]) / 20)
                AAProperty1.append([(j - meanI) / fenmu for j in i])
            encodings = []
            header = ["SampleName"]
            for aa in AA:
                header.append("PAAC_Xc1." + aa)
            for n in range(1, lambdaValue + 1):
                header.append("PAAC_Xc2.lambda" + str(n))
            encodings.append(header)
            for i in self.seq_list:
                name, sequence = i[0], i[1]
                if len(sequence) < lambdaValue + 1:
                    if self.verbose:
                        print(
                            f"Skipping {name} PAAC requires sequence length > lambdaValue+1: {str(lambdaValue + 1)}"
                        )
                        continue
                code = [name]
                theta = []
                for n in range(1, lambdaValue + 1):
                    theta.append(
                        sum(
                            [
                                self.Rvalue(
                                    sequence[j],
                                    sequence[j + n],
                                    AADict,
                                    AAProperty1,
                                )
                                for j in range(len(sequence) - n)
                            ]
                        )
                        / (len(sequence) - n)
                    )
                myDict = {}
                for aa in AA:
                    myDict[aa] = sequence.count(aa)
                code = code + [myDict[aa] / (1 + w * sum(theta)) for aa in AA]
                code = code + [(w * j) / (1 + w * sum(theta)) for j in theta]
                encodings.append(code)
            encodings = np.array(encodings)
            self.encodings = pd.DataFrame(
                encodings[1:, 1:].astype(float),
                columns=encodings[0, 1:],
                index=encodings[1:, 0],
            )
            return True
        except Exception as e:
            sys.exit(f"Error: {e}")

    def _APAAC(self):
        try:
            lambdaValue = self.__default_para["lambdaValue"]
            w = self.__default_para["weight"]
            dataFile = os.path.join(self.datadir, "PAAC.txt")
            with open(dataFile) as f:
                records = f.readlines()
            AA = "".join(records[0].rstrip().split()[1:])
            AADict = {}
            for i in range(len(AA)):
                AADict[AA[i]] = i
            AAProperty = []
            AAPropertyNames = []
            for i in range(1, len(records) - 1):
                array = (
                    records[i].rstrip().split()
                    if records[i].rstrip() != ""
                    else None
                )
                AAProperty.append([float(j) for j in array[1:]])
                AAPropertyNames.append(array[0])
            AAProperty1 = []
            for i in AAProperty:
                meanI = sum(i) / 20
                fenmu = math.sqrt(sum([(j - meanI) ** 2 for j in i]) / 20)
                AAProperty1.append([(j - meanI) / fenmu for j in i])
            encodings = []
            header = ["SampleName"]
            for i in AA:
                header.append("APAAC_Pc1." + i)
            for j in range(1, lambdaValue + 1):
                for i in AAPropertyNames:
                    header.append("APAAC_Pc2." + i + "." + str(j))
            encodings.append(header)
            for i in self.seq_list:
                name, sequence = i[0], i[1]
                if lambdaValue > len(sequence) - 1:
                    if self.verbose:
                        print(
                            f"Skipping {name} APAAC requires sequence length - 1 > lambdaValue: {str(lambdaValue)}"
                        )
                    continue
                code = [name]
                theta = []
                for n in range(1, lambdaValue + 1):
                    for j in range(len(AAProperty1)):
                        theta.append(
                            sum(
                                [
                                    AAProperty1[j][AADict[sequence[k]]]
                                    * AAProperty1[j][AADict[sequence[k + n]]]
                                    for k in range(len(sequence) - n)
                                ]
                            )
                            / (len(sequence) - n)
                        )
                myDict = {}
                for aa in AA:
                    myDict[aa] = sequence.count(aa)

                code = code + [myDict[aa] / (1 + w * sum(theta)) for aa in AA]
                code = code + [
                    w * value / (1 + w * sum(theta)) for value in theta
                ]
                encodings.append(code)
            encodings = np.array(encodings)
            self.encodings = pd.DataFrame(
                encodings[1:, 1:].astype(float),
                columns=encodings[0, 1:],
                index=encodings[1:, 0],
            )
            return True
        except Exception as e:
            sys.exit(f"Error: {e}")

    def _ASDC(self):
        try:
            encodings = []
            aaPairs = []
            for aa1 in self.AA:
                for aa2 in self.AA:
                    aaPairs.append(aa1 + aa2)
            header = ["SampleName"]
            header += [
                "ASDC_" + aa1 + aa2 for aa1 in self.AA for aa2 in self.AA
            ]
            encodings.append(header)
            for i in self.seq_list:
                name, sequence = i[0], i[1]
                code = [name]
                sum = 0
                pair_dict = {}
                for pair in aaPairs:
                    pair_dict[pair] = 0
                for j in range(len(sequence)):
                    for k in range(j + 1, len(sequence)):
                        if sequence[j] in self.AA and sequence[k] in self.AA:
                            pair_dict[sequence[j] + sequence[k]] += 1
                            sum += 1
                for pair in aaPairs:
                    code.append(pair_dict[pair] / sum)
                encodings.append(code)
            encodings = np.array(encodings)
            self.encodings = pd.DataFrame(
                encodings[1:, 1:].astype(float),
                columns=encodings[0, 1:],
                index=encodings[1:, 0],
            )
            return True
        except Exception as e:
            sys.exit(f"Error: {e}")

    def _DistancePair(self):
        try:
            cp20_dict = {
                "A": "A",
                "C": "C",
                "D": "D",
                "E": "E",
                "F": "F",
                "G": "G",
                "H": "H",
                "I": "I",
                "K": "K",
                "L": "L",
                "M": "M",
                "N": "N",
                "P": "P",
                "Q": "Q",
                "R": "R",
                "S": "S",
                "T": "T",
                "V": "V",
                "W": "W",
                "Y": "Y",
            }
            cp19_dict = {
                "A": "A",
                "C": "C",
                "D": "D",
                "E": "E",
                "F": "F",
                "G": "G",
                "H": "H",
                "I": "I",
                "K": "K",
                "L": "L",
                "M": "M",
                "N": "N",
                "P": "P",
                "Q": "Q",
                "R": "R",
                "S": "S",
                "T": "T",
                "V": "V",
                "W": "W",
                "Y": "F",  # YF
            }
            cp14_dict = {
                "A": "A",
                "C": "C",
                "D": "D",
                "E": "E",
                "F": "F",
                "G": "G",
                "H": "H",
                "I": "I",
                "K": "H",  # HRKQ
                "L": "L",
                "M": "I",  # IMV
                "N": "N",
                "P": "P",
                "Q": "H",  # HRKQ
                "R": "H",  # HRKQ
                "S": "S",
                "T": "T",
                "V": "I",  # IMV
                "W": "W",
                "Y": "W",  # WY
            }
            cp13_dict = {
                "A": "A",
                "C": "C",
                "D": "D",
                "E": "E",
                "F": "F",
                "G": "G",
                "H": "H",
                "I": "I",
                "K": "K",
                "L": "I",  # IL
                "M": "F",  # FM
                "N": "N",
                "P": "H",  # HPQWY
                "Q": "H",  # HPQWY
                "R": "K",  # KR
                "S": "S",
                "T": "T",
                "V": "V",
                "W": "H",  # HPQWY
                "Y": "H",  # HPQWY
            }
            cp20_AA = "ACDEFGHIKLMNPQRSTVWY"
            cp19_AA = "ACDEFGHIKLMNPQRSTVW"
            cp14_AA = "ACDEFGHILNPSTW"
            cp13_AA = "ACDEFGHIKNSTV"
            distance = self.__default_para["distance"]
            cp = self.__default_para["cp"]
            AA = cp20_AA
            AA_dict = cp20_dict
            if cp == "cp(19)":
                AA = cp19_AA
                AA_dict = cp19_dict
            if cp == "cp(14)":
                AA = cp14_AA
                AA_dict = cp14_dict
            if cp == "cp(13)":
                AA = cp13_AA
                AA_dict = cp13_dict
            encodings = []
            pair_dict = {}
            single_dict = {}
            for aa1 in AA:
                single_dict[aa1] = 0
                for aa2 in AA:
                    pair_dict[aa1 + aa2] = 0
            header = ["SampleName"]
            for d in range(distance + 1):
                if d == 0:
                    for key in sorted(single_dict):
                        header.append("DP_{0}".format(key))
                else:
                    for key in sorted(pair_dict):
                        header.append("DP_%s.distance%s" % (key, d))
            encodings.append(header)
            for elem in self.seq_list:
                name, sequence = elem[0], elem[1]
                if len(sequence) < distance + 1:
                    if self.verbose:
                        print(
                            f"Skipping {name} DistancePair requires sequence length > distance+1: {str(distance + 1)}"
                        )
                    continue
                code = [name]
                for d in range(distance + 1):
                    if d == 0:
                        tmp_dict = single_dict.copy()
                        for i in range(len(sequence)):
                            tmp_dict[AA_dict[sequence[i]]] += 1
                        for key in sorted(tmp_dict):
                            code.append(tmp_dict[key] / len(sequence))
                    else:
                        tmp_dict = pair_dict.copy()
                        for i in range(len(sequence) - d):
                            tmp_dict[
                                AA_dict[sequence[i]] + AA_dict[sequence[i + d]]
                            ] += 1
                        for key in sorted(tmp_dict):
                            code.append(tmp_dict[key] / (len(sequence) - d))
                encodings.append(code)
            encodings = np.array(encodings)
            self.encodings = pd.DataFrame(
                encodings[1:, 1:].astype(float),
                columns=encodings[0, 1:],
                index=encodings[1:, 0],
            )
            return True
        except Exception as e:
            sys.exit(f"Error: {e}")

    def gapModel(self, fastas, myDict, gDict, gNames, ktuple, glValue, ttype):
        encodings = []
        header = ["SampleName"]
        if ktuple == 1:
            header = header + [
                ttype + "_" + g + "_gap" + str(glValue) for g in gNames
            ]
            encodings.append(header)
            for i in fastas:
                name, sequence = i[0], i[1]
                code = [name]
                numDict = {}
                for j in range(0, len(sequence), glValue + 1):
                    numDict[gDict[myDict[sequence[j]]]] = (
                        numDict.get(gDict[myDict[sequence[j]]], 0) + 1
                    )

                for g in gNames:
                    code.append(numDict.get(g, 0))
                encodings.append(code)
        if ktuple == 2:
            header = header + [
                ttype + "_" + g1 + "_" + g2 + "_gap" + str(glValue)
                for g1 in gNames
                for g2 in gNames
            ]
            encodings.append(header)
            for i in fastas:
                name, sequence = i[0], i[1]
                code = [name]
                numDict = {}
                for j in range(0, len(sequence), glValue + 1):
                    if j + 1 < len(sequence):
                        numDict[
                            gDict[myDict[sequence[j]]]
                            + "_"
                            + gDict[myDict[sequence[j + 1]]]
                        ] = (
                            numDict.get(
                                gDict[myDict[sequence[j]]]
                                + "_"
                                + gDict[myDict[sequence[j + 1]]],
                                0,
                            )
                            + 1
                        )

                for g in [g1 + "_" + g2 for g1 in gNames for g2 in gNames]:
                    code.append(numDict.get(g, 0))
                encodings.append(code)
        if ktuple == 3:
            header = header + [
                ttype + "_" + g1 + "_" + g2 + "_" + g3 + "_gap" + str(glValue)
                for g1 in gNames
                for g2 in gNames
                for g3 in gNames
            ]
            encodings.append(header)
            for i in fastas:
                name, sequence = i[0], i[1]
                code = [name]
                numDict = {}
                for j in range(0, len(sequence), glValue + 1):
                    if j + 1 < len(sequence) and j + 2 < len(sequence):
                        numDict[
                            gDict[myDict[sequence[j]]]
                            + "_"
                            + gDict[myDict[sequence[j + 1]]]
                            + "_"
                            + gDict[myDict[sequence[j + 2]]]
                        ] = (
                            numDict.get(
                                gDict[myDict[sequence[j]]]
                                + "_"
                                + gDict[myDict[sequence[j + 1]]]
                                + "_"
                                + gDict[myDict[sequence[j + 2]]],
                                0,
                            )
                            + 1
                        )
                for g in [
                    g1 + "_" + g2 + "_" + g3
                    for g1 in gNames
                    for g2 in gNames
                    for g3 in gNames
                ]:
                    code.append(numDict.get(g, 0))
                encodings.append(code)
        return encodings

    def lambdaModel(
        self, fastas, myDict, gDict, gNames, ktuple, glValue, ttype
    ):
        encodings = []
        header = ["SampleName"]
        if ktuple == 1:
            header = header + [
                ttype + "_" + g + "_LC" + str(glValue) for g in gNames
            ]
            encodings.append(header)
            for i in fastas:
                name, sequence = i[0], i[1]
                code = [name]
                numDict = {}
                for j in range(0, len(sequence)):
                    numDict[gDict[myDict[sequence[j]]]] = (
                        numDict.get(gDict[myDict[sequence[j]]], 0) + 1
                    )

                for g in gNames:
                    code.append(numDict.get(g, 0))
                encodings.append(code)
        if ktuple == 2:
            header = header + [
                ttype + "_" + g1 + "_" + g2 + "_LC" + str(glValue)
                for g1 in gNames
                for g2 in gNames
            ]
            encodings.append(header)
            for i in fastas:
                name, sequence = i[0], i[1]
                code = [name]
                numDict = {}
                for j in range(0, len(sequence)):
                    if j + glValue < len(sequence):
                        numDict[
                            gDict[myDict[sequence[j]]]
                            + "_"
                            + gDict[myDict[sequence[j + glValue]]]
                        ] = (
                            numDict.get(
                                gDict[myDict[sequence[j]]]
                                + "_"
                                + gDict[myDict[sequence[j + glValue]]],
                                0,
                            )
                            + 1
                        )

                for g in [g1 + "_" + g2 for g1 in gNames for g2 in gNames]:
                    code.append(numDict.get(g, 0))
                encodings.append(code)
        if ktuple == 3:
            header = header + [
                ttype + "_" + g1 + "_" + g2 + "_" + g3 + "_LC" + str(glValue)
                for g1 in gNames
                for g2 in gNames
                for g3 in gNames
            ]
            encodings.append(header)
            for i in fastas:
                name, sequence = i[0], i[1]
                code = [name]
                numDict = {}
                for j in range(0, len(sequence)):
                    if j + glValue < len(sequence) and j + 2 * glValue < len(
                        sequence
                    ):
                        numDict[
                            gDict[myDict[sequence[j]]]
                            + "_"
                            + gDict[myDict[sequence[j + glValue]]]
                            + "_"
                            + gDict[myDict[sequence[j + 2 * glValue]]]
                        ] = (
                            numDict.get(
                                gDict[myDict[sequence[j]]]
                                + "_"
                                + gDict[myDict[sequence[j + glValue]]]
                                + "_"
                                + gDict[myDict[sequence[j + 2 * glValue]]],
                                0,
                            )
                            + 1
                        )
                for g in [
                    g1 + "_" + g2 + "_" + g3
                    for g1 in gNames
                    for g2 in gNames
                    for g3 in gNames
                ]:
                    code.append(numDict.get(g, 0))
                encodings.append(code)
        return encodings

    def _PseKRAAC_type_1(self):
        try:
            AAGroup = {
                2: ["CMFILVWY", "AGTSNQDEHRKP"],
                3: ["CMFILVWY", "AGTSP", "NQDEHRK"],
                4: ["CMFWY", "ILV", "AGTS", "NQDEHRKP"],
                5: ["WFYH", "MILV", "CATSP", "G", "NQDERK"],
                6: ["WFYH", "MILV", "CATS", "P", "G", "NQDERK"],
                7: ["WFYH", "MILV", "CATS", "P", "G", "NQDE", "RK"],
                8: ["WFYH", "MILV", "CA", "NTS", "P", "G", "DE", "QRK"],
                9: ["WFYH", "MI", "LV", "CA", "NTS", "P", "G", "DE", "QRK"],
                10: [
                    "WFY",
                    "ML",
                    "IV",
                    "CA",
                    "TS",
                    "NH",
                    "P",
                    "G",
                    "DE",
                    "QRK",
                ],
                11: [
                    "WFY",
                    "ML",
                    "IV",
                    "CA",
                    "TS",
                    "NH",
                    "P",
                    "G",
                    "D",
                    "QE",
                    "RK",
                ],
                12: [
                    "WFY",
                    "ML",
                    "IV",
                    "C",
                    "A",
                    "TS",
                    "NH",
                    "P",
                    "G",
                    "D",
                    "QE",
                    "RK",
                ],
                13: [
                    "WFY",
                    "ML",
                    "IV",
                    "C",
                    "A",
                    "T",
                    "S",
                    "NH",
                    "P",
                    "G",
                    "D",
                    "QE",
                    "RK",
                ],
                14: [
                    "WFY",
                    "ML",
                    "IV",
                    "C",
                    "A",
                    "T",
                    "S",
                    "NH",
                    "P",
                    "G",
                    "D",
                    "QE",
                    "R",
                    "K",
                ],
                15: [
                    "WFY",
                    "ML",
                    "IV",
                    "C",
                    "A",
                    "T",
                    "S",
                    "N",
                    "H",
                    "P",
                    "G",
                    "D",
                    "QE",
                    "R",
                    "K",
                ],
                16: [
                    "W",
                    "FY",
                    "ML",
                    "IV",
                    "C",
                    "A",
                    "T",
                    "S",
                    "N",
                    "H",
                    "P",
                    "G",
                    "D",
                    "QE",
                    "R",
                    "K",
                ],
                17: [
                    "W",
                    "FY",
                    "ML",
                    "IV",
                    "C",
                    "A",
                    "T",
                    "S",
                    "N",
                    "H",
                    "P",
                    "G",
                    "D",
                    "Q",
                    "E",
                    "R",
                    "K",
                ],
                18: [
                    "W",
                    "FY",
                    "M",
                    "L",
                    "IV",
                    "C",
                    "A",
                    "T",
                    "S",
                    "N",
                    "H",
                    "P",
                    "G",
                    "D",
                    "Q",
                    "E",
                    "R",
                    "K",
                ],
                19: [
                    "W",
                    "F",
                    "Y",
                    "M",
                    "L",
                    "IV",
                    "C",
                    "A",
                    "T",
                    "S",
                    "N",
                    "H",
                    "P",
                    "G",
                    "D",
                    "Q",
                    "E",
                    "R",
                    "K",
                ],
                20: [
                    "W",
                    "F",
                    "Y",
                    "M",
                    "L",
                    "I",
                    "V",
                    "C",
                    "A",
                    "T",
                    "S",
                    "N",
                    "H",
                    "P",
                    "G",
                    "D",
                    "Q",
                    "E",
                    "R",
                    "K",
                ],
            }
            fastas = self.seq_list
            subtype = self.__default_para["PseKRAAC_model"]
            raactype = self.__default_para["RAAC_clust"]
            ktuple = self.__default_para["k-tuple"]
            glValue = (
                self.__default_para["g-gap"]
                if subtype == "g-gap"
                else self.__default_para["lambdaValue"]
            )
            if raactype not in AAGroup:
                print(raactype)
                self.error_msg = "PseKRAAC descriptor value error."
                return False
            # index each amino acids to their group
            myDict = {}
            for i in range(len(AAGroup[raactype])):
                for aa in AAGroup[raactype][i]:
                    myDict[aa] = i
            gDict = {}
            gNames = []
            for i in range(len(AAGroup[raactype])):
                gDict[i] = "T1.G." + str(i + 1)
                gNames.append("T1.G." + str(i + 1))
            encodings = []
            if subtype == "g-gap":
                encodings = self.gapModel(
                    fastas, myDict, gDict, gNames, ktuple, glValue, "type1"
                )
            else:
                encodings = self.lambdaModel(
                    fastas, myDict, gDict, gNames, ktuple, glValue, "type1"
                )
            encodings = np.array(encodings)
            self.encodings = pd.DataFrame(
                encodings[1:, 1:].astype(float),
                columns=encodings[0, 1:],
                index=encodings[1:, 0],
            )
            return True
        except Exception as e:
            sys.exit(f"Error: {e}")

    def _PseKRAAC_type_2(self):
        try:
            AAGroup = {
                2: ["LVIMCAGSTPFYW", "EDNQKRH"],
                3: ["LVIMCAGSTP", "FYW", "EDNQKRH"],
                4: ["LVIMC", "AGSTP", "FYW", "EDNQKRH"],
                5: ["LVIMC", "AGSTP", "FYW", "EDNQ", "KRH"],
                6: ["LVIM", "AGST", "PHC", "FYW", "EDNQ", "KR"],
                8: ["LVIMC", "AG", "ST", "P", "FYW", "EDNQ", "KR", "H"],
                15: [
                    "LVIM",
                    "C",
                    "A",
                    "G",
                    "S",
                    "T",
                    "P",
                    "FY",
                    "W",
                    "E",
                    "D",
                    "N",
                    "Q",
                    "KR",
                    "H",
                ],
                20: [
                    "L",
                    "V",
                    "I",
                    "M",
                    "C",
                    "A",
                    "G",
                    "S",
                    "T",
                    "P",
                    "F",
                    "Y",
                    "W",
                    "E",
                    "D",
                    "N",
                    "Q",
                    "K",
                    "R",
                    "H",
                ],
            }
            fastas = self.seq_list
            subtype = self.__default_para["PseKRAAC_model"]
            raactype = self.__default_para["RAAC_clust"]
            ktuple = self.__default_para["k-tuple"]
            glValue = (
                self.__default_para["g-gap"]
                if subtype == "g-gap"
                else self.__default_para["lambdaValue"]
            )
            if raactype not in AAGroup:
                self.error_msg = "PseKRAAC descriptor value error."
                return False
            # index each amino acids to their group
            myDict = {}
            for i in range(len(AAGroup[raactype])):
                for aa in AAGroup[raactype][i]:
                    myDict[aa] = i
            gDict = {}
            gNames = []
            for i in range(len(AAGroup[raactype])):
                gDict[i] = "T1.G." + str(i + 1)
                gNames.append("T1.G." + str(i + 1))
            encodings = []
            if subtype == "g-gap":
                encodings = self.gapModel(
                    fastas, myDict, gDict, gNames, ktuple, glValue, "type2"
                )
            else:
                encodings = self.lambdaModel(
                    fastas, myDict, gDict, gNames, ktuple, glValue, "type2"
                )
            encodings = np.array(encodings)
            self.encodings = pd.DataFrame(
                encodings[1:, 1:].astype(float),
                columns=encodings[0, 1:],
                index=encodings[1:, 0],
            )
            return True
        except Exception as e:
            sys.exit(f"Error: {e}")

    def _PseKRAAC_type_3A(self):
        try:
            AAGroup = {
                2: ["AGSPDEQNHTKRMILFYVC", "W"],
                3: ["AGSPDEQNHTKRMILFYV", "W", "C"],
                4: ["AGSPDEQNHTKRMIV", "W", "YFL", "C"],
                5: ["AGSPDEQNHTKR", "W", "YF", "MIVL", "C"],
                6: ["AGSP", "DEQNHTKR", "W", "YF", "MIL", "VC"],
                7: ["AGP", "DEQNH", "TKRMIV", "W", "YF", "L", "CS"],
                8: ["AG", "DEQN", "TKRMIV", "HY", "W", "L", "FP", "CS"],
                9: ["AG", "P", "DEQN", "TKRMI", "HY", "W", "F", "L", "VCS"],
                10: [
                    "AG",
                    "P",
                    "DEQN",
                    "TKRM",
                    "HY",
                    "W",
                    "F",
                    "I",
                    "L",
                    "VCS",
                ],
                11: [
                    "AG",
                    "P",
                    "DEQN",
                    "TK",
                    "RI",
                    "H",
                    "Y",
                    "W",
                    "F",
                    "ML",
                    "VCS",
                ],
                12: [
                    "FAS",
                    "P",
                    "G",
                    "DEQ",
                    "NL",
                    "TK",
                    "R",
                    "H",
                    "W",
                    "Y",
                    "IM",
                    "VC",
                ],
                13: [
                    "FAS",
                    "P",
                    "G",
                    "DEQ",
                    "NL",
                    "T",
                    "K",
                    "R",
                    "H",
                    "W",
                    "Y",
                    "IM",
                    "VC",
                ],
                14: [
                    "FA",
                    "P",
                    "G",
                    "T",
                    "DE",
                    "QM",
                    "NL",
                    "K",
                    "R",
                    "H",
                    "W",
                    "Y",
                    "IV",
                    "CS",
                ],
                15: [
                    "FAS",
                    "P",
                    "G",
                    "T",
                    "DE",
                    "Q",
                    "NL",
                    "K",
                    "R",
                    "H",
                    "W",
                    "Y",
                    "M",
                    "I",
                    "VC",
                ],
                16: [
                    "FA",
                    "P",
                    "G",
                    "ST",
                    "DE",
                    "Q",
                    "N",
                    "K",
                    "R",
                    "H",
                    "W",
                    "Y",
                    "M",
                    "L",
                    "I",
                    "VC",
                ],
                17: [
                    "FA",
                    "P",
                    "G",
                    "S",
                    "T",
                    "DE",
                    "Q",
                    "N",
                    "K",
                    "R",
                    "H",
                    "W",
                    "Y",
                    "M",
                    "L",
                    "I",
                    "VC",
                ],
                18: [
                    "FA",
                    "P",
                    "G",
                    "S",
                    "T",
                    "DE",
                    "Q",
                    "N",
                    "K",
                    "R",
                    "H",
                    "W",
                    "Y",
                    "M",
                    "L",
                    "I",
                    "V",
                    "C",
                ],
                19: [
                    "FA",
                    "P",
                    "G",
                    "S",
                    "T",
                    "D",
                    "E",
                    "Q",
                    "N",
                    "K",
                    "R",
                    "H",
                    "W",
                    "Y",
                    "M",
                    "L",
                    "I",
                    "V",
                    "C",
                ],
                20: [
                    "F",
                    "A",
                    "P",
                    "G",
                    "S",
                    "T",
                    "D",
                    "E",
                    "Q",
                    "N",
                    "K",
                    "R",
                    "H",
                    "W",
                    "Y",
                    "M",
                    "L",
                    "I",
                    "V",
                    "C",
                ],
            }
            fastas = self.seq_list
            subtype = self.__default_para["PseKRAAC_model"]
            raactype = self.__default_para["RAAC_clust"]
            ktuple = self.__default_para["k-tuple"]
            glValue = (
                self.__default_para["g-gap"]
                if subtype == "g-gap"
                else self.__default_para["lambdaValue"]
            )
            if raactype not in AAGroup:
                self.error_msg = "PseKRAAC descriptor value error."
                return False
            # index each amino acids to their group
            myDict = {}
            for i in range(len(AAGroup[raactype])):
                for aa in AAGroup[raactype][i]:
                    myDict[aa] = i
            gDict = {}
            gNames = []
            for i in range(len(AAGroup[raactype])):
                gDict[i] = "T1.G." + str(i + 1)
                gNames.append("T1.G." + str(i + 1))
            encodings = []
            if subtype == "g-gap":
                encodings = self.gapModel(
                    fastas, myDict, gDict, gNames, ktuple, glValue, "type3A"
                )
            else:
                encodings = self.lambdaModel(
                    fastas, myDict, gDict, gNames, ktuple, glValue, "type3A"
                )
            encodings = np.array(encodings)
            self.encodings = pd.DataFrame(
                encodings[1:, 1:].astype(float),
                columns=encodings[0, 1:],
                index=encodings[1:, 0],
            )
            return True
        except Exception as e:
            sys.exit(f"Error: {e}")

    def _PseKRAAC_type_3B(self):
        try:
            AAGroup = {
                2: ["HRKQNEDSTGPACVIM", "LFYW"],
                3: ["HRKQNEDSTGPACVIM", "LFY", "W"],
                4: ["HRKQNEDSTGPA", "CIV", "MLFY", "W"],
                5: ["HRKQNEDSTGPA", "CV", "IML", "FY", "W"],
                6: ["HRKQNEDSTPA", "G", "CV", "IML", "FY", "W"],
                7: ["HRKQNEDSTA", "G", "P", "CV", "IML", "FY", "W"],
                8: ["HRKQSTA", "NED", "G", "P", "CV", "IML", "FY", "W"],
                9: ["HRKQ", "NED", "ASTG", "P", "C", "IV", "MLF", "Y", "W"],
                10: [
                    "RKHSA",
                    "Q",
                    "NED",
                    "G",
                    "P",
                    "C",
                    "TIV",
                    "MLF",
                    "Y",
                    "W",
                ],
                11: [
                    "RKQ",
                    "NG",
                    "ED",
                    "AST",
                    "P",
                    "C",
                    "IV",
                    "HML",
                    "F",
                    "Y",
                    "W",
                ],
                12: [
                    "RKQ",
                    "ED",
                    "NAST",
                    "G",
                    "P",
                    "C",
                    "IV",
                    "H",
                    "ML",
                    "F",
                    "Y",
                    "W",
                ],
                13: [
                    "RK",
                    "QE",
                    "D",
                    "NG",
                    "HA",
                    "ST",
                    "P",
                    "C",
                    "IV",
                    "ML",
                    "F",
                    "Y",
                    "W",
                ],
                14: [
                    "R",
                    "K",
                    "QE",
                    "D",
                    "NG",
                    "HA",
                    "ST",
                    "P",
                    "C",
                    "IV",
                    "ML",
                    "F",
                    "Y",
                    "W",
                ],
                15: [
                    "R",
                    "K",
                    "QE",
                    "D",
                    "NG",
                    "HA",
                    "ST",
                    "P",
                    "C",
                    "IV",
                    "M",
                    "L",
                    "F",
                    "Y",
                    "W",
                ],
                16: [
                    "R",
                    "K",
                    "Q",
                    "E",
                    "D",
                    "NG",
                    "HA",
                    "ST",
                    "P",
                    "C",
                    "IV",
                    "M",
                    "L",
                    "F",
                    "Y",
                    "W",
                ],
                17: [
                    "R",
                    "K",
                    "Q",
                    "E",
                    "D",
                    "NG",
                    "HA",
                    "S",
                    "T",
                    "P",
                    "C",
                    "IV",
                    "M",
                    "L",
                    "F",
                    "Y",
                    "W",
                ],
                18: [
                    "R",
                    "K",
                    "Q",
                    "E",
                    "D",
                    "NG",
                    "HA",
                    "S",
                    "T",
                    "P",
                    "C",
                    "I",
                    "V",
                    "M",
                    "L",
                    "F",
                    "Y",
                    "W",
                ],
                19: [
                    "R",
                    "K",
                    "Q",
                    "E",
                    "D",
                    "NG",
                    "H",
                    "A",
                    "S",
                    "T",
                    "P",
                    "C",
                    "I",
                    "V",
                    "M",
                    "L",
                    "F",
                    "Y",
                    "W",
                ],
                20: [
                    "R",
                    "K",
                    "Q",
                    "E",
                    "D",
                    "N",
                    "G",
                    "H",
                    "A",
                    "S",
                    "T",
                    "P",
                    "C",
                    "I",
                    "V",
                    "M",
                    "L",
                    "F",
                    "Y",
                    "W",
                ],
            }
            fastas = self.seq_list
            subtype = self.__default_para["PseKRAAC_model"]
            raactype = self.__default_para["RAAC_clust"]
            ktuple = self.__default_para["k-tuple"]
            glValue = (
                self.__default_para["g-gap"]
                if subtype == "g-gap"
                else self.__default_para["lambdaValue"]
            )
            if raactype not in AAGroup:
                self.error_msg = "PseKRAAC descriptor value error."
                return False
            # index each amino acids to their group
            myDict = {}
            for i in range(len(AAGroup[raactype])):
                for aa in AAGroup[raactype][i]:
                    myDict[aa] = i
            gDict = {}
            gNames = []
            for i in range(len(AAGroup[raactype])):
                gDict[i] = "T1.G." + str(i + 1)
                gNames.append("T1.G." + str(i + 1))
            encodings = []
            if subtype == "g-gap":
                encodings = self.gapModel(
                    fastas, myDict, gDict, gNames, ktuple, glValue, "type3B"
                )
            else:
                encodings = self.lambdaModel(
                    fastas, myDict, gDict, gNames, ktuple, glValue, "type3B"
                )
            encodings = np.array(encodings)
            self.encodings = pd.DataFrame(
                encodings[1:, 1:].astype(float),
                columns=encodings[0, 1:],
                index=encodings[1:, 0],
            )
            return True
        except Exception as e:
            sys.exit(f"Error: {e}")

    def _PseKRAAC_type_4(self):
        try:
            AAGroup = {
                5: ["G", "IVFYW", "ALMEQRK", "P", "NDHSTC"],
                8: ["G", "IV", "FYW", "ALM", "EQRK", "P", "ND", "HSTC"],
                9: ["G", "IV", "FYW", "ALM", "EQRK", "P", "ND", "HS", "TC"],
                11: [
                    "G",
                    "IV",
                    "FYW",
                    "A",
                    "LM",
                    "EQRK",
                    "P",
                    "ND",
                    "HS",
                    "T",
                    "C",
                ],
                13: [
                    "G",
                    "IV",
                    "FYW",
                    "A",
                    "L",
                    "M",
                    "E",
                    "QRK",
                    "P",
                    "ND",
                    "HS",
                    "T",
                    "C",
                ],
                20: [
                    "G",
                    "I",
                    "V",
                    "F",
                    "Y",
                    "W",
                    "A",
                    "L",
                    "M",
                    "E",
                    "Q",
                    "R",
                    "K",
                    "P",
                    "N",
                    "D",
                    "H",
                    "S",
                    "T",
                    "C",
                ],
            }
            fastas = self.seq_list
            subtype = self.__default_para["PseKRAAC_model"]
            raactype = self.__default_para["RAAC_clust"]
            ktuple = self.__default_para["k-tuple"]
            glValue = (
                self.__default_para["g-gap"]
                if subtype == "g-gap"
                else self.__default_para["lambdaValue"]
            )
            if raactype not in AAGroup:
                self.error_msg = "PseKRAAC descriptor value error."
                return False
            # index each amino acids to their group
            myDict = {}
            for i in range(len(AAGroup[raactype])):
                for aa in AAGroup[raactype][i]:
                    myDict[aa] = i
            gDict = {}
            gNames = []
            for i in range(len(AAGroup[raactype])):
                gDict[i] = "T1.G." + str(i + 1)
                gNames.append("T1.G." + str(i + 1))
            encodings = []
            if subtype == "g-gap":
                encodings = self.gapModel(
                    fastas, myDict, gDict, gNames, ktuple, glValue, "type4"
                )
            else:
                encodings = self.lambdaModel(
                    fastas, myDict, gDict, gNames, ktuple, glValue, "type4"
                )
            encodings = np.array(encodings)
            self.encodings = pd.DataFrame(
                encodings[1:, 1:].astype(float),
                columns=encodings[0, 1:],
                index=encodings[1:, 0],
            )
            return True
        except Exception as e:
            sys.exit(f"Error: {e}")

    def _PseKRAAC_type_5(self):
        try:
            AAGroup = {
                3: ["FWYCILMVAGSTPHNQ", "DE", "KR"],
                4: ["FWY", "CILMV", "AGSTP", "EQNDHKR"],
                8: ["FWY", "CILMV", "GA", "ST", "P", "EQND", "H", "KR"],
                10: [
                    "G",
                    "FYW",
                    "A",
                    "ILMV",
                    "RK",
                    "P",
                    "EQND",
                    "H",
                    "ST",
                    "C",
                ],
                15: [
                    "G",
                    "FY",
                    "W",
                    "A",
                    "ILMV",
                    "E",
                    "Q",
                    "RK",
                    "P",
                    "N",
                    "D",
                    "H",
                    "S",
                    "T",
                    "C",
                ],
                20: [
                    "G",
                    "I",
                    "V",
                    "F",
                    "Y",
                    "W",
                    "A",
                    "L",
                    "M",
                    "E",
                    "Q",
                    "R",
                    "K",
                    "P",
                    "N",
                    "D",
                    "H",
                    "S",
                    "T",
                    "C",
                ],
            }
            fastas = self.seq_list
            subtype = self.__default_para["PseKRAAC_model"]
            raactype = self.__default_para["RAAC_clust"]
            ktuple = self.__default_para["k-tuple"]
            glValue = (
                self.__default_para["g-gap"]
                if subtype == "g-gap"
                else self.__default_para["lambdaValue"]
            )
            if raactype not in AAGroup:
                self.error_msg = "PseKRAAC descriptor value error."
                return False
            # index each amino acids to their group
            myDict = {}
            for i in range(len(AAGroup[raactype])):
                for aa in AAGroup[raactype][i]:
                    myDict[aa] = i
            gDict = {}
            gNames = []
            for i in range(len(AAGroup[raactype])):
                gDict[i] = "T1.G." + str(i + 1)
                gNames.append("T1.G." + str(i + 1))
            encodings = []
            if subtype == "g-gap":
                encodings = self.gapModel(
                    fastas, myDict, gDict, gNames, ktuple, glValue, "type5"
                )
            else:
                encodings = self.lambdaModel(
                    fastas, myDict, gDict, gNames, ktuple, glValue, "type5"
                )
            encodings = np.array(encodings)
            self.encodings = pd.DataFrame(
                encodings[1:, 1:].astype(float),
                columns=encodings[0, 1:],
                index=encodings[1:, 0],
            )
            return True
        except Exception as e:
            sys.exit(f"Error: {e}")

    def _PseKRAAC_type_6A(self):
        try:
            AAGroup = {
                4: ["AGPST", "CILMV", "DEHKNQR", "FYW"],
                5: ["AHT", "CFILMVWY", "DE", "GP", "KNQRS"],
                20: [
                    "A",
                    "C",
                    "D",
                    "E",
                    "F",
                    "G",
                    "H",
                    "I",
                    "K",
                    "L",
                    "M",
                    "N",
                    "P",
                    "Q",
                    "R",
                    "S",
                    "T",
                    "V",
                    "W",
                    "Y",
                ],
            }
            fastas = self.seq_list
            subtype = self.__default_para["PseKRAAC_model"]
            raactype = self.__default_para["RAAC_clust"]
            ktuple = self.__default_para["k-tuple"]
            glValue = (
                self.__default_para["g-gap"]
                if subtype == "g-gap"
                else self.__default_para["lambdaValue"]
            )

            if raactype not in AAGroup:
                self.error_msg = "PseKRAAC descriptor value error."
                return False

            # index each amino acids to their group
            myDict = {}
            for i in range(len(AAGroup[raactype])):
                for aa in AAGroup[raactype][i]:
                    myDict[aa] = i

            gDict = {}
            gNames = []
            for i in range(len(AAGroup[raactype])):
                gDict[i] = "T1.G." + str(i + 1)
                gNames.append("T1.G." + str(i + 1))

            encodings = []
            if subtype == "g-gap":
                encodings = self.gapModel(
                    fastas, myDict, gDict, gNames, ktuple, glValue, "type6A"
                )
            else:
                encodings = self.lambdaModel(
                    fastas, myDict, gDict, gNames, ktuple, glValue, "type6A"
                )
            encodings = np.array(encodings)
            self.encodings = pd.DataFrame(
                encodings[1:, 1:].astype(float),
                columns=encodings[0, 1:],
                index=encodings[1:, 0],
            )
            return True
        except Exception as e:
            sys.exit(f"Error: {e}")

    def _PseKRAAC_type_6B(self):
        try:
            AAGroup = {
                5: ["AEHKQRST", "CFILMVWY", "DN", "G", "P"],
            }
            fastas = self.seq_list
            subtype = self.__default_para["PseKRAAC_model"]
            raactype = self.__default_para["RAAC_clust"]
            ktuple = self.__default_para["k-tuple"]
            glValue = (
                self.__default_para["g-gap"]
                if subtype == "g-gap"
                else self.__default_para["lambdaValue"]
            )

            if raactype not in AAGroup:
                self.error_msg = "PseKRAAC descriptor value error."
                return False

            # index each amino acids to their group
            myDict = {}
            for i in range(len(AAGroup[raactype])):
                for aa in AAGroup[raactype][i]:
                    myDict[aa] = i

            gDict = {}
            gNames = []
            for i in range(len(AAGroup[raactype])):
                gDict[i] = "T1.G." + str(i + 1)
                gNames.append("T1.G." + str(i + 1))

            encodings = []
            if subtype == "g-gap":
                encodings = self.gapModel(
                    fastas, myDict, gDict, gNames, ktuple, glValue, "type6B"
                )
            else:
                encodings = self.lambdaModel(
                    fastas, myDict, gDict, gNames, ktuple, glValue, "type6B"
                )
            encodings = np.array(encodings)
            self.encodings = pd.DataFrame(
                encodings[1:, 1:].astype(float),
                columns=encodings[0, 1:],
                index=encodings[1:, 0],
            )
            return True
        except Exception as e:
            sys.exit(f"Error: {e}")

    def _PseKRAAC_type_6C(self):
        try:
            AAGroup = {
                5: ["AG", "C", "DEKNPQRST", "FILMVWY", "H"],
            }
            fastas = self.seq_list
            subtype = self.__default_para["PseKRAAC_model"]
            raactype = self.__default_para["RAAC_clust"]
            ktuple = self.__default_para["k-tuple"]
            glValue = (
                self.__default_para["g-gap"]
                if subtype == "g-gap"
                else self.__default_para["lambdaValue"]
            )
            if raactype not in AAGroup:
                self.error_msg = "PseKRAAC descriptor value error."
                return False
            # index each amino acids to their group
            myDict = {}
            for i in range(len(AAGroup[raactype])):
                for aa in AAGroup[raactype][i]:
                    myDict[aa] = i
            gDict = {}
            gNames = []
            for i in range(len(AAGroup[raactype])):
                gDict[i] = "T1.G." + str(i + 1)
                gNames.append("T1.G." + str(i + 1))
            encodings = []
            if subtype == "g-gap":
                encodings = self.gapModel(
                    fastas, myDict, gDict, gNames, ktuple, glValue, "type6C"
                )
            else:
                encodings = self.lambdaModel(
                    fastas, myDict, gDict, gNames, ktuple, glValue, "type6C"
                )
            encodings = np.array(encodings)
            self.encodings = pd.DataFrame(
                encodings[1:, 1:].astype(float),
                columns=encodings[0, 1:],
                index=encodings[1:, 0],
            )
            return True
        except Exception as e:
            sys.exit(f"Error: {e}")

    def _PseKRAAC_type_7(self):
        try:
            AAGroup = {
                2: ["C", "MFILVWYAGTSNQDEHRKP"],
                3: ["C", "MFILVWYAKR", "GTSNQDEHP"],
                4: ["C", "KR", "MFILVWYA", "GTSNQDEHP"],
                5: ["C", "KR", "MFILVWYA", "DE", "GTSNQHP"],
                6: ["C", "KR", "WYA", "MFILV", "DE", "GTSNQHP"],
                7: ["C", "KR", "WYA", "MFILV", "DE", "QH", "GTSNP"],
                8: ["C", "KR", "WYA", "MFILV", "D", "E", "QH", "GTSNP"],
                9: ["C", "KR", "WYA", "MFILV", "D", "E", "QH", "TP", "GSN"],
                10: [
                    "C",
                    "KR",
                    "WY",
                    "A",
                    "MFILV",
                    "D",
                    "E",
                    "QH",
                    "TP",
                    "GSN",
                ],
                11: [
                    "C",
                    "K",
                    "R",
                    "WY",
                    "A",
                    "MFILV",
                    "D",
                    "E",
                    "QH",
                    "TP",
                    "GSN",
                ],
                12: [
                    "C",
                    "K",
                    "R",
                    "WY",
                    "A",
                    "MFILV",
                    "D",
                    "E",
                    "QH",
                    "TP",
                    "GS",
                    "N",
                ],
                13: [
                    "C",
                    "K",
                    "R",
                    "W",
                    "Y",
                    "A",
                    "MFILV",
                    "D",
                    "E",
                    "QH",
                    "TP",
                    "GS",
                    "N",
                ],
                14: [
                    "C",
                    "K",
                    "R",
                    "W",
                    "Y",
                    "A",
                    "FILV",
                    "M",
                    "D",
                    "E",
                    "QH",
                    "TP",
                    "GS",
                    "N",
                ],
                15: [
                    "C",
                    "K",
                    "R",
                    "W",
                    "Y",
                    "A",
                    "FILV",
                    "M",
                    "D",
                    "E",
                    "Q",
                    "H",
                    "TP",
                    "GS",
                    "N",
                ],
                16: [
                    "C",
                    "K",
                    "R",
                    "W",
                    "Y",
                    "A",
                    "FILV",
                    "M",
                    "D",
                    "E",
                    "Q",
                    "H",
                    "TP",
                    "G",
                    "S",
                    "N",
                ],
                17: [
                    "C",
                    "K",
                    "R",
                    "W",
                    "Y",
                    "A",
                    "FI",
                    "LV",
                    "M",
                    "D",
                    "E",
                    "Q",
                    "H",
                    "TP",
                    "G",
                    "S",
                    "N",
                ],
                18: [
                    "C",
                    "K",
                    "R",
                    "W",
                    "Y",
                    "A",
                    "FI",
                    "LV",
                    "M",
                    "D",
                    "E",
                    "Q",
                    "H",
                    "T",
                    "P",
                    "G",
                    "S",
                    "N",
                ],
                19: [
                    "C",
                    "K",
                    "R",
                    "W",
                    "Y",
                    "A",
                    "F",
                    "I",
                    "LV",
                    "M",
                    "D",
                    "E",
                    "Q",
                    "H",
                    "T",
                    "P",
                    "G",
                    "S",
                    "N",
                ],
                20: [
                    "C",
                    "K",
                    "R",
                    "W",
                    "Y",
                    "A",
                    "F",
                    "I",
                    "L",
                    "V",
                    "M",
                    "D",
                    "E",
                    "Q",
                    "H",
                    "T",
                    "P",
                    "G",
                    "S",
                    "N",
                ],
            }
            fastas = self.seq_list
            subtype = self.__default_para["PseKRAAC_model"]
            raactype = self.__default_para["RAAC_clust"]
            ktuple = self.__default_para["k-tuple"]
            glValue = (
                self.__default_para["g-gap"]
                if subtype == "g-gap"
                else self.__default_para["lambdaValue"]
            )

            if raactype not in AAGroup:
                self.error_msg = "PseKRAAC descriptor value error."
                return False

            # index each amino acids to their group
            myDict = {}
            for i in range(len(AAGroup[raactype])):
                for aa in AAGroup[raactype][i]:
                    myDict[aa] = i

            gDict = {}
            gNames = []
            for i in range(len(AAGroup[raactype])):
                gDict[i] = "T1.G." + str(i + 1)
                gNames.append("T1.G." + str(i + 1))

            encodings = []
            if subtype == "g-gap":
                encodings = self.gapModel(
                    fastas, myDict, gDict, gNames, ktuple, glValue, "type7"
                )
            else:
                encodings = self.lambdaModel(
                    fastas, myDict, gDict, gNames, ktuple, glValue, "type7"
                )
            encodings = np.array(encodings)
            self.encodings = pd.DataFrame(
                encodings[1:, 1:].astype(float),
                columns=encodings[0, 1:],
                index=encodings[1:, 0],
            )
            return True
        except Exception as e:
            sys.exit(f"Error: {e}")

    def _PseKRAAC_type_8(self):
        try:
            AAGroup = {
                2: ["ADEGKNPQRST", "CFHILMVWY"],
                3: ["ADEGNPST", "CHKQRW", "FILMVY"],
                4: ["AGNPST", "CHWY", "DEKQR", "FILMV"],
                5: ["AGPST", "CFWY", "DEN", "HKQR", "ILMV"],
                6: ["APST", "CW", "DEGN", "FHY", "ILMV", "KQR"],
                7: ["AGST", "CW", "DEN", "FY", "HP", "ILMV", "KQR"],
                8: ["AST", "CG", "DEN", "FY", "HP", "ILV", "KQR", "MW"],
                9: ["AST", "CW", "DE", "FY", "GN", "HQ", "ILV", "KR", "MP"],
                10: [
                    "AST",
                    "CW",
                    "DE",
                    "FY",
                    "GN",
                    "HQ",
                    "IV",
                    "KR",
                    "LM",
                    "P",
                ],
                11: [
                    "AST",
                    "C",
                    "DE",
                    "FY",
                    "GN",
                    "HQ",
                    "IV",
                    "KR",
                    "LM",
                    "P",
                    "W",
                ],
                12: [
                    "AST",
                    "C",
                    "DE",
                    "FY",
                    "G",
                    "HQ",
                    "IV",
                    "KR",
                    "LM",
                    "N",
                    "P",
                    "W",
                ],
                13: [
                    "AST",
                    "C",
                    "DE",
                    "FY",
                    "G",
                    "H",
                    "IV",
                    "KR",
                    "LM",
                    "N",
                    "P",
                    "Q",
                    "W",
                ],
                14: [
                    "AST",
                    "C",
                    "DE",
                    "FL",
                    "G",
                    "H",
                    "IV",
                    "KR",
                    "M",
                    "N",
                    "P",
                    "Q",
                    "W",
                    "Y",
                ],
                15: [
                    "AST",
                    "C",
                    "DE",
                    "F",
                    "G",
                    "H",
                    "IV",
                    "KR",
                    "L",
                    "M",
                    "N",
                    "P",
                    "Q",
                    "W",
                    "Y",
                ],
                16: [
                    "AT",
                    "C",
                    "DE",
                    "F",
                    "G",
                    "H",
                    "IV",
                    "KR",
                    "L",
                    "M",
                    "N",
                    "P",
                    "Q",
                    "S",
                    "W",
                    "Y",
                ],
                17: [
                    "AT",
                    "C",
                    "DE",
                    "F",
                    "G",
                    "H",
                    "IV",
                    "K",
                    "L",
                    "M",
                    "N",
                    "P",
                    "Q",
                    "R",
                    "S",
                    "W",
                    "Y",
                ],
                18: [
                    "A",
                    "C",
                    "DE",
                    "F",
                    "G",
                    "H",
                    "IV",
                    "K",
                    "L",
                    "M",
                    "N",
                    "P",
                    "Q",
                    "R",
                    "S",
                    "T",
                    "W",
                    "Y",
                ],
                19: [
                    "A",
                    "C",
                    "D",
                    "E",
                    "F",
                    "G",
                    "H",
                    "IV",
                    "K",
                    "L",
                    "M",
                    "N",
                    "P",
                    "Q",
                    "R",
                    "S",
                    "T",
                    "W",
                    "Y",
                ],
                20: [
                    "A",
                    "C",
                    "D",
                    "E",
                    "F",
                    "G",
                    "H",
                    "I",
                    "V",
                    "K",
                    "L",
                    "M",
                    "N",
                    "P",
                    "Q",
                    "R",
                    "S",
                    "T",
                    "W",
                    "Y",
                ],
            }
            fastas = self.seq_list
            subtype = self.__default_para["PseKRAAC_model"]
            raactype = self.__default_para["RAAC_clust"]
            ktuple = self.__default_para["k-tuple"]
            glValue = (
                self.__default_para["g-gap"]
                if subtype == "g-gap"
                else self.__default_para["lambdaValue"]
            )
            if raactype not in AAGroup:
                self.error_msg = "PseKRAAC descriptor value error."
                return False
            # index each amino acids to their group
            myDict = {}
            for i in range(len(AAGroup[raactype])):
                for aa in AAGroup[raactype][i]:
                    myDict[aa] = i

            gDict = {}
            gNames = []
            for i in range(len(AAGroup[raactype])):
                gDict[i] = "T1.G." + str(i + 1)
                gNames.append("T1.G." + str(i + 1))

            encodings = []
            if subtype == "g-gap":
                encodings = self.gapModel(
                    fastas, myDict, gDict, gNames, ktuple, glValue, "type8"
                )
            else:
                encodings = self.lambdaModel(
                    fastas, myDict, gDict, gNames, ktuple, glValue, "type8"
                )
            encodings = np.array(encodings)
            self.encodings = pd.DataFrame(
                encodings[1:, 1:].astype(float),
                columns=encodings[0, 1:],
                index=encodings[1:, 0],
            )
            return True
        except Exception as e:
            sys.exit(f"Error: {e}")

    def _PseKRAAC_type_9(self):
        try:
            AAGroup = {
                2: ["ACDEFGHILMNPQRSTVWY", "K"],
                3: ["ACDFGMPQRSTW", "EHILNVY", "K"],
                4: ["AGPT", "CDFMQRSW", "EHILNVY", "K"],
                5: ["AGPT", "CDQ", "EHILNVY", "FMRSW", "K"],
                6: ["AG", "CDQ", "EHILNVY", "FMRSW", "K", "PT"],
                7: ["AG", "CDQ", "EHNY", "FMRSW", "ILV", "K", "PT"],
                8: ["AG", "C", "DQ", "EHNY", "FMRSW", "ILV", "K", "PT"],
                9: ["AG", "C", "DQ", "EHNY", "FMW", "ILV", "K", "PT", "RS"],
                10: [
                    "A",
                    "C",
                    "DQ",
                    "EHNY",
                    "FMW",
                    "G",
                    "ILV",
                    "K",
                    "PT",
                    "RS",
                ],
                11: [
                    "A",
                    "C",
                    "DQ",
                    "EHNY",
                    "FM",
                    "G",
                    "ILV",
                    "K",
                    "PT",
                    "RS",
                    "W",
                ],
                12: [
                    "A",
                    "C",
                    "DQ",
                    "EHNY",
                    "FM",
                    "G",
                    "IL",
                    "K",
                    "PT",
                    "RS",
                    "V",
                    "W",
                ],
                13: [
                    "A",
                    "C",
                    "DQ",
                    "E",
                    "FM",
                    "G",
                    "HNY",
                    "IL",
                    "K",
                    "PT",
                    "RS",
                    "V",
                    "W",
                ],
                14: [
                    "A",
                    "C",
                    "D",
                    "E",
                    "FM",
                    "G",
                    "HNY",
                    "IL",
                    "K",
                    "PT",
                    "Q",
                    "RS",
                    "V",
                    "W",
                ],
                15: [
                    "A",
                    "C",
                    "D",
                    "E",
                    "FM",
                    "G",
                    "HNY",
                    "IL",
                    "K",
                    "PT",
                    "Q",
                    "R",
                    "S",
                    "V",
                    "W",
                ],
                16: [
                    "A",
                    "C",
                    "D",
                    "E",
                    "F",
                    "G",
                    "HNY",
                    "IL",
                    "K",
                    "M",
                    "PT",
                    "Q",
                    "R",
                    "S",
                    "V",
                    "W",
                ],
                17: [
                    "A",
                    "C",
                    "D",
                    "E",
                    "F",
                    "G",
                    "HNY",
                    "IL",
                    "K",
                    "M",
                    "P",
                    "Q",
                    "R",
                    "S",
                    "T",
                    "V",
                    "W",
                ],
                18: [
                    "A",
                    "C",
                    "D",
                    "E",
                    "F",
                    "G",
                    "HNY",
                    "I",
                    "K",
                    "L",
                    "M",
                    "P",
                    "Q",
                    "R",
                    "S",
                    "T",
                    "V",
                    "W",
                ],
                19: [
                    "A",
                    "C",
                    "D",
                    "E",
                    "F",
                    "G",
                    "HN",
                    "I",
                    "K",
                    "L",
                    "M",
                    "P",
                    "Q",
                    "R",
                    "S",
                    "T",
                    "V",
                    "W",
                    "Y",
                ],
                20: [
                    "A",
                    "C",
                    "D",
                    "E",
                    "F",
                    "G",
                    "H",
                    "N",
                    "I",
                    "K",
                    "L",
                    "M",
                    "P",
                    "Q",
                    "R",
                    "S",
                    "T",
                    "V",
                    "W",
                    "Y",
                ],
            }
            fastas = self.seq_list
            subtype = self.__default_para["PseKRAAC_model"]
            raactype = self.__default_para["RAAC_clust"]
            ktuple = self.__default_para["k-tuple"]
            glValue = (
                self.__default_para["g-gap"]
                if subtype == "g-gap"
                else self.__default_para["lambdaValue"]
            )
            if raactype not in AAGroup:
                self.error_msg = "PseKRAAC descriptor value error."
                return False
            # index each amino acids to their group
            myDict = {}
            for i in range(len(AAGroup[raactype])):
                for aa in AAGroup[raactype][i]:
                    myDict[aa] = i
            gDict = {}
            gNames = []
            for i in range(len(AAGroup[raactype])):
                gDict[i] = "T1.G." + str(i + 1)
                gNames.append("T1.G." + str(i + 1))
            encodings = []
            if subtype == "g-gap":
                encodings = self.gapModel(
                    fastas, myDict, gDict, gNames, ktuple, glValue, "type9"
                )
            else:
                encodings = self.lambdaModel(
                    fastas, myDict, gDict, gNames, ktuple, glValue, "type9"
                )
            encodings = np.array(encodings)
            self.encodings = pd.DataFrame(
                encodings[1:, 1:].astype(float),
                columns=encodings[0, 1:],
                index=encodings[1:, 0],
            )
            return True
        except Exception as e:
            sys.exit(f"Error: {e}")

    def _PseKRAAC_type_10(self):
        try:
            AAGroup = {
                2: ["CMFILVWY", "AGTSNQDEHRKP"],
                3: ["CMFILVWY", "AGTSP", "NQDEHRK"],
                4: ["CMFWY", "ILV", "AGTS", "NQDEHRKP"],
                5: ["FWYH", "MILV", "CATSP", "G", "NQDERK"],
                6: ["FWYH", "MILV", "CATS", "P", "G", "NQDERK"],
                7: ["FWYH", "MILV", "CATS", "P", "G", "NQDE", "RK"],
                8: ["FWYH", "MILV", "CA", "NTS", "P", "G", "DE", "QRK"],
                9: ["FWYH", "ML", "IV", "CA", "NTS", "P", "G", "DE", "QRK"],
                10: [
                    "FWY",
                    "ML",
                    "IV",
                    "CA",
                    "TS",
                    "NH",
                    "P",
                    "G",
                    "DE",
                    "QRK",
                ],
                11: [
                    "FWY",
                    "ML",
                    "IV",
                    "CA",
                    "TS",
                    "NH",
                    "P",
                    "G",
                    "D",
                    "QE",
                    "RK",
                ],
                12: [
                    "FWY",
                    "ML",
                    "IV",
                    "C",
                    "A",
                    "TS",
                    "NH",
                    "P",
                    "G",
                    "D",
                    "QE",
                    "RK",
                ],
                13: [
                    "FWY",
                    "ML",
                    "IV",
                    "C",
                    "A",
                    "T",
                    "S",
                    "NH",
                    "P",
                    "G",
                    "D",
                    "QE",
                    "RK",
                ],
                14: [
                    "FWY",
                    "ML",
                    "IV",
                    "C",
                    "A",
                    "T",
                    "S",
                    "NH",
                    "P",
                    "G",
                    "D",
                    "QE",
                    "R",
                    "K",
                ],
                15: [
                    "FWY",
                    "ML",
                    "IV",
                    "C",
                    "A",
                    "T",
                    "S",
                    "N",
                    "H",
                    "P",
                    "G",
                    "D",
                    "QE",
                    "R",
                    "K",
                ],
                16: [
                    "W",
                    "FY",
                    "ML",
                    "IV",
                    "C",
                    "A",
                    "T",
                    "S",
                    "N",
                    "H",
                    "P",
                    "G",
                    "D",
                    "QE",
                    "R",
                    "K",
                ],
                17: [
                    "W",
                    "FY",
                    "ML",
                    "IV",
                    "C",
                    "A",
                    "T",
                    "S",
                    "N",
                    "H",
                    "P",
                    "G",
                    "D",
                    "Q",
                    "E",
                    "R",
                    "K",
                ],
                18: [
                    "W",
                    "FY",
                    "M",
                    "L",
                    "IV",
                    "C",
                    "A",
                    "T",
                    "S",
                    "N",
                    "H",
                    "P",
                    "G",
                    "D",
                    "Q",
                    "E",
                    "R",
                    "K",
                ],
                19: [
                    "W",
                    "F",
                    "Y",
                    "M",
                    "L",
                    "IV",
                    "C",
                    "A",
                    "T",
                    "S",
                    "N",
                    "H",
                    "P",
                    "G",
                    "D",
                    "Q",
                    "E",
                    "R",
                    "K",
                ],
                20: [
                    "W",
                    "F",
                    "Y",
                    "M",
                    "L",
                    "I",
                    "V",
                    "C",
                    "A",
                    "T",
                    "S",
                    "N",
                    "H",
                    "P",
                    "G",
                    "D",
                    "Q",
                    "E",
                    "R",
                    "K",
                ],
            }
            fastas = self.seq_list
            subtype = self.__default_para["PseKRAAC_model"]
            raactype = self.__default_para["RAAC_clust"]
            ktuple = self.__default_para["k-tuple"]
            glValue = (
                self.__default_para["g-gap"]
                if subtype == "g-gap"
                else self.__default_para["lambdaValue"]
            )
            if raactype not in AAGroup:
                self.error_msg = "PseKRAAC descriptor value error."
                return False
            # index each amino acids to their group
            myDict = {}
            for i in range(len(AAGroup[raactype])):
                for aa in AAGroup[raactype][i]:
                    myDict[aa] = i
            gDict = {}
            gNames = []
            for i in range(len(AAGroup[raactype])):
                gDict[i] = "T1.G." + str(i + 1)
                gNames.append("T1.G." + str(i + 1))
            encodings = []
            if subtype == "g-gap":
                encodings = self.gapModel(
                    fastas, myDict, gDict, gNames, ktuple, glValue, "tpye10"
                )
            else:
                encodings = self.lambdaModel(
                    fastas, myDict, gDict, gNames, ktuple, glValue, "type10"
                )
            encodings = np.array(encodings)
            self.encodings = pd.DataFrame(
                encodings[1:, 1:].astype(float),
                columns=encodings[0, 1:],
                index=encodings[1:, 0],
            )
            return True
        except Exception as e:
            sys.exit(f"Error: {e}")

    def _PseKRAAC_type_11(self):
        try:
            AAGroup = {
                2: ["CFYWMLIV", "GPATSNHQEDRK"],
                3: ["CFYWMLIV", "GPATS", "NHQEDRK"],
                4: ["CFYW", "MLIV", "GPATS", "NHQEDRK"],
                5: ["CFYW", "MLIV", "G", "PATS", "NHQEDRK"],
                6: ["CFYW", "MLIV", "G", "P", "ATS", "NHQEDRK"],
                7: ["CFYW", "MLIV", "G", "P", "ATS", "NHQED", "RK"],
                8: ["CFYW", "MLIV", "G", "P", "ATS", "NH", "QED", "RK"],
                9: ["CFYW", "ML", "IV", "G", "P", "ATS", "NH", "QED", "RK"],
                10: [
                    "C",
                    "FYW",
                    "ML",
                    "IV",
                    "G",
                    "P",
                    "ATS",
                    "NH",
                    "QED",
                    "RK",
                ],
                11: [
                    "C",
                    "FYW",
                    "ML",
                    "IV",
                    "G",
                    "P",
                    "A",
                    "TS",
                    "NH",
                    "QED",
                    "RK",
                ],
                12: [
                    "C",
                    "FYW",
                    "ML",
                    "IV",
                    "G",
                    "P",
                    "A",
                    "TS",
                    "NH",
                    "QE",
                    "D",
                    "RK",
                ],
                13: [
                    "C",
                    "FYW",
                    "ML",
                    "IV",
                    "G",
                    "P",
                    "A",
                    "T",
                    "S",
                    "NH",
                    "QE",
                    "D",
                    "RK",
                ],
                14: [
                    "C",
                    "FYW",
                    "ML",
                    "IV",
                    "G",
                    "P",
                    "A",
                    "T",
                    "S",
                    "N",
                    "H",
                    "QE",
                    "D",
                    "RK",
                ],
                15: [
                    "C",
                    "FYW",
                    "ML",
                    "IV",
                    "G",
                    "P",
                    "A",
                    "T",
                    "S",
                    "N",
                    "H",
                    "QE",
                    "D",
                    "R",
                    "K",
                ],
                16: [
                    "C",
                    "FY",
                    "W",
                    "ML",
                    "IV",
                    "G",
                    "P",
                    "A",
                    "T",
                    "S",
                    "N",
                    "H",
                    "QE",
                    "D",
                    "R",
                    "K",
                ],
                17: [
                    "C",
                    "FY",
                    "W",
                    "ML",
                    "IV",
                    "G",
                    "P",
                    "A",
                    "T",
                    "S",
                    "N",
                    "H",
                    "Q",
                    "E",
                    "D",
                    "R",
                    "K",
                ],
                18: [
                    "C",
                    "FY",
                    "W",
                    "M",
                    "L",
                    "IV",
                    "G",
                    "P",
                    "A",
                    "T",
                    "S",
                    "N",
                    "H",
                    "Q",
                    "E",
                    "D",
                    "R",
                    "K",
                ],
                19: [
                    "C",
                    "F",
                    "Y",
                    "W",
                    "M",
                    "L",
                    "IV",
                    "G",
                    "P",
                    "A",
                    "T",
                    "S",
                    "N",
                    "H",
                    "Q",
                    "E",
                    "D",
                    "R",
                    "K",
                ],
                20: [
                    "C",
                    "F",
                    "Y",
                    "W",
                    "M",
                    "L",
                    "I",
                    "V",
                    "G",
                    "P",
                    "A",
                    "T",
                    "S",
                    "N",
                    "H",
                    "Q",
                    "E",
                    "D",
                    "R",
                    "K",
                ],
            }
            fastas = self.seq_list
            subtype = self.__default_para["PseKRAAC_model"]
            raactype = self.__default_para["RAAC_clust"]
            ktuple = self.__default_para["k-tuple"]
            glValue = (
                self.__default_para["g-gap"]
                if subtype == "g-gap"
                else self.__default_para["lambdaValue"]
            )
            if raactype not in AAGroup:
                self.error_msg = "PseKRAAC descriptor value error."
                return False
            # index each amino acids to their group
            myDict = {}
            for i in range(len(AAGroup[raactype])):
                for aa in AAGroup[raactype][i]:
                    myDict[aa] = i
            gDict = {}
            gNames = []
            for i in range(len(AAGroup[raactype])):
                gDict[i] = "T1.G." + str(i + 1)
                gNames.append("T1.G." + str(i + 1))
            encodings = []
            if subtype == "g-gap":
                encodings = self.gapModel(
                    fastas, myDict, gDict, gNames, ktuple, glValue, "type11"
                )
            else:
                encodings = self.lambdaModel(
                    fastas, myDict, gDict, gNames, ktuple, glValue, "type11"
                )
            encodings = np.array(encodings)
            self.encodings = pd.DataFrame(
                encodings[1:, 1:].astype(float),
                columns=encodings[0, 1:],
                index=encodings[1:, 0],
            )
            return True
        except Exception as e:
            sys.exit(f"Error: {e}")

    def _PseKRAAC_type_12(self):
        try:
            AAGroup = {
                2: ["IVMLFWYC", "ARNDQEGHKPST"],
                3: ["IVLMFWC", "YA", "RNDQEGHKPST"],
                4: ["IVLMFW", "C", "YA", "RNDQEGHKPST"],
                5: ["IVLMFW", "C", "YA", "G", "RNDQEHKPST"],
                6: ["IVLMF", "WY", "C", "AH", "G", "RNDQEKPST"],
                7: ["IVLMF", "WY", "C", "AH", "GP", "R", "NDQEKST"],
                8: ["IVLMF", "WY", "C", "A", "G", "R", "Q", "NDEHKPST"],
                9: ["IVLMF", "WY", "C", "A", "G", "P", "H", "K", "RNDQEST"],
                10: [
                    "IVLM",
                    "F",
                    "W",
                    "Y",
                    "C",
                    "A",
                    "H",
                    "G",
                    "RN",
                    "DQEKPST",
                ],
                11: [
                    "IVLMF",
                    "W",
                    "Y",
                    "C",
                    "A",
                    "H",
                    "G",
                    "R",
                    "N",
                    "Q",
                    "DEKPST",
                ],
                12: [
                    "IVLM",
                    "F",
                    "W",
                    "Y",
                    "C",
                    "A",
                    "H",
                    "G",
                    "N",
                    "Q",
                    "T",
                    "RDEKPS",
                ],
                13: [
                    "IVLM",
                    "F",
                    "W",
                    "Y",
                    "C",
                    "A",
                    "H",
                    "G",
                    "N",
                    "Q",
                    "P",
                    "R",
                    "DEKST",
                ],
                14: [
                    "IVLM",
                    "F",
                    "W",
                    "Y",
                    "C",
                    "A",
                    "H",
                    "G",
                    "N",
                    "Q",
                    "P",
                    "R",
                    "K",
                    "DEST",
                ],
                15: [
                    "IVLM",
                    "F",
                    "W",
                    "Y",
                    "C",
                    "A",
                    "H",
                    "G",
                    "N",
                    "Q",
                    "P",
                    "R",
                    "K",
                    "D",
                    "EST",
                ],
                16: [
                    "IVLM",
                    "F",
                    "W",
                    "Y",
                    "C",
                    "A",
                    "H",
                    "G",
                    "N",
                    "Q",
                    "P",
                    "R",
                    "K",
                    "S",
                    "T",
                    "DE",
                ],
                17: [
                    "IVL",
                    "M",
                    "F",
                    "W",
                    "Y",
                    "C",
                    "A",
                    "H",
                    "G",
                    "N",
                    "Q",
                    "P",
                    "R",
                    "K",
                    "S",
                    "T",
                    "DE",
                ],
                18: [
                    "IVL",
                    "M",
                    "F",
                    "W",
                    "Y",
                    "C",
                    "A",
                    "H",
                    "G",
                    "N",
                    "Q",
                    "P",
                    "R",
                    "K",
                    "S",
                    "T",
                    "D",
                    "E",
                ],
                20: [
                    "I",
                    "V",
                    "L",
                    "M",
                    "F",
                    "W",
                    "Y",
                    "C",
                    "A",
                    "H",
                    "G",
                    "N",
                    "Q",
                    "P",
                    "R",
                    "K",
                    "S",
                    "T",
                    "D",
                    "E",
                ],
            }
            fastas = self.seq_list
            subtype = self.__default_para["PseKRAAC_model"]
            raactype = self.__default_para["RAAC_clust"]
            ktuple = self.__default_para["k-tuple"]
            glValue = (
                self.__default_para["g-gap"]
                if subtype == "g-gap"
                else self.__default_para["lambdaValue"]
            )
            if raactype not in AAGroup:
                self.error_msg = "PseKRAAC descriptor value error."
                return False
            # index each amino acids to their group
            myDict = {}
            for i in range(len(AAGroup[raactype])):
                for aa in AAGroup[raactype][i]:
                    myDict[aa] = i
            gDict = {}
            gNames = []
            for i in range(len(AAGroup[raactype])):
                gDict[i] = "T1.G." + str(i + 1)
                gNames.append("T1.G." + str(i + 1))
            encodings = []
            if subtype == "g-gap":
                encodings = self.gapModel(
                    fastas, myDict, gDict, gNames, ktuple, glValue, "type12"
                )
            else:
                encodings = self.lambdaModel(
                    fastas, myDict, gDict, gNames, ktuple, glValue, "type12"
                )
            encodings = np.array(encodings)
            self.encodings = pd.DataFrame(
                encodings[1:, 1:].astype(float),
                columns=encodings[0, 1:],
                index=encodings[1:, 0],
            )
            return True
        except Exception as e:
            sys.exit(f"Error: {e}")

    def _PseKRAAC_type_13(self):
        try:
            AAGroup = {
                4: ["ADKERNTSQ", "YFLIVMCWH", "G", "P"],
                12: [
                    "A",
                    "D",
                    "KER",
                    "N",
                    "TSQ",
                    "YF",
                    "LIVM",
                    "C",
                    "W",
                    "H",
                    "G",
                    "P",
                ],
                17: [
                    "A",
                    "D",
                    "KE",
                    "R",
                    "N",
                    "T",
                    "S",
                    "Q",
                    "Y",
                    "F",
                    "LIV",
                    "M",
                    "C",
                    "W",
                    "H",
                    "G",
                    "P",
                ],
                20: [
                    "A",
                    "D",
                    "K",
                    "E",
                    "R",
                    "N",
                    "T",
                    "S",
                    "Q",
                    "Y",
                    "F",
                    "L",
                    "I",
                    "V",
                    "M",
                    "C",
                    "W",
                    "H",
                    "G",
                    "P",
                ],
            }
            fastas = self.seq_list
            subtype = self.__default_para["PseKRAAC_model"]
            raactype = self.__default_para["RAAC_clust"]
            ktuple = self.__default_para["k-tuple"]
            glValue = (
                self.__default_para["g-gap"]
                if subtype == "g-gap"
                else self.__default_para["lambdaValue"]
            )

            if raactype not in AAGroup:
                self.error_msg = "PseKRAAC descriptor value error."
                return False

            # index each amino acids to their group
            myDict = {}
            for i in range(len(AAGroup[raactype])):
                for aa in AAGroup[raactype][i]:
                    myDict[aa] = i

            gDict = {}
            gNames = []
            for i in range(len(AAGroup[raactype])):
                gDict[i] = "T1.G." + str(i + 1)
                gNames.append("T1.G." + str(i + 1))

            encodings = []
            if subtype == "g-gap":
                encodings = self.gapModel(
                    fastas, myDict, gDict, gNames, ktuple, glValue, "type13"
                )
            else:
                encodings = self.lambdaModel(
                    fastas, myDict, gDict, gNames, ktuple, glValue, "type13"
                )
            encodings = np.array(encodings)
            self.encodings = pd.DataFrame(
                encodings[1:, 1:].astype(float),
                columns=encodings[0, 1:],
                index=encodings[1:, 0],
            )
            return True
        except Exception as e:
            sys.exit(f"Error: {e}")

    def _PseKRAAC_type_14(self):
        try:
            AAGroup = {
                2: ["ARNDCQEGHKPST", "ILMFWYV"],
                3: ["ARNDQEGHKPST", "C", "ILMFWYV"],
                4: ["ARNDQEGHKPST", "C", "ILMFYV", "W"],
                5: ["AGPST", "RNDQEHK", "C", "ILMFYV", "W"],
                6: ["AGPST", "RNDQEK", "C", "H", "ILMFYV", "W"],
                7: ["ANDGST", "RQEK", "C", "H", "ILMFYV", "P", "W"],
                8: ["ANDGST", "RQEK", "C", "H", "ILMV", "FY", "P", "W"],
                9: ["AGST", "RQEK", "ND", "C", "H", "ILMV", "FY", "P", "W"],
                10: [
                    "AGST",
                    "RK",
                    "ND",
                    "C",
                    "QE",
                    "H",
                    "ILMV",
                    "FY",
                    "P",
                    "W",
                ],
                11: [
                    "AST",
                    "RK",
                    "ND",
                    "C",
                    "QE",
                    "G",
                    "H",
                    "ILMV",
                    "FY",
                    "P",
                    "W",
                ],
                12: [
                    "AST",
                    "RK",
                    "ND",
                    "C",
                    "QE",
                    "G",
                    "H",
                    "IV",
                    "LM",
                    "FY",
                    "P",
                    "W",
                ],
                13: [
                    "AST",
                    "RK",
                    "N",
                    "D",
                    "C",
                    "QE",
                    "G",
                    "H",
                    "IV",
                    "LM",
                    "FY",
                    "P",
                    "W",
                ],
                14: [
                    "AST",
                    "RK",
                    "N",
                    "D",
                    "C",
                    "Q",
                    "E",
                    "G",
                    "H",
                    "IV",
                    "LM",
                    "FY",
                    "P",
                    "W",
                ],
                15: [
                    "A",
                    "RK",
                    "N",
                    "D",
                    "C",
                    "Q",
                    "E",
                    "G",
                    "H",
                    "IV",
                    "LM",
                    "FY",
                    "P",
                    "ST",
                    "W",
                ],
                16: [
                    "A",
                    "RK",
                    "N",
                    "D",
                    "C",
                    "Q",
                    "E",
                    "G",
                    "H",
                    "IV",
                    "LM",
                    "F",
                    "P",
                    "ST",
                    "W",
                    "Y",
                ],
                17: [
                    "A",
                    "R",
                    "N",
                    "D",
                    "C",
                    "Q",
                    "E",
                    "G",
                    "H",
                    "IV",
                    "LM",
                    "K",
                    "F",
                    "P",
                    "ST",
                    "W",
                    "Y",
                ],
                18: [
                    "A",
                    "R",
                    "N",
                    "D",
                    "C",
                    "Q",
                    "E",
                    "G",
                    "H",
                    "IV",
                    "LM",
                    "K",
                    "F",
                    "P",
                    "S",
                    "T",
                    "W",
                    "Y",
                ],
                19: [
                    "A",
                    "R",
                    "N",
                    "D",
                    "C",
                    "Q",
                    "E",
                    "G",
                    "H",
                    "IV",
                    "L",
                    "K",
                    "M",
                    "F",
                    "P",
                    "S",
                    "T",
                    "W",
                    "Y",
                ],
                20: [
                    "A",
                    "R",
                    "N",
                    "D",
                    "C",
                    "Q",
                    "E",
                    "G",
                    "H",
                    "I",
                    "V",
                    "L",
                    "K",
                    "M",
                    "F",
                    "P",
                    "S",
                    "T",
                    "W",
                    "Y",
                ],
            }
            fastas = self.seq_list
            subtype = self.__default_para["PseKRAAC_model"]
            raactype = self.__default_para["RAAC_clust"]
            ktuple = self.__default_para["k-tuple"]
            glValue = (
                self.__default_para["g-gap"]
                if subtype == "g-gap"
                else self.__default_para["lambdaValue"]
            )

            if raactype not in AAGroup:
                self.error_msg = "PseKRAAC descriptor value error."
                return False

            # index each amino acids to their group
            myDict = {}
            for i in range(len(AAGroup[raactype])):
                for aa in AAGroup[raactype][i]:
                    myDict[aa] = i

            gDict = {}
            gNames = []
            for i in range(len(AAGroup[raactype])):
                gDict[i] = "T1.G." + str(i + 1)
                gNames.append("T1.G." + str(i + 1))

            encodings = []
            if subtype == "g-gap":
                encodings = self.gapModel(
                    fastas, myDict, gDict, gNames, ktuple, glValue, "type14"
                )
            else:
                encodings = self.lambdaModel(
                    fastas, myDict, gDict, gNames, ktuple, glValue, "type14"
                )
            encodings = np.array(encodings)
            self.encodings = pd.DataFrame(
                encodings[1:, 1:].astype(float),
                columns=encodings[0, 1:],
                index=encodings[1:, 0],
            )
            return True
        except Exception as e:
            sys.exit(f"Error: {e}")

    def _PseKRAAC_type_15(self):
        try:
            AAGroup = {
                2: ["MFILVAW", "CYQHPGTSNRKDE"],
                3: ["MFILVAW", "CYQHPGTSNRK", "DE"],
                4: ["MFILV", "ACW", "YQHPGTSNRK", "DE"],
                5: ["MFILV", "ACW", "YQHPGTSN", "RK", "DE"],
                6: ["MFILV", "A", "C", "WYQHPGTSN", "RK", "DE"],
                7: ["MFILV", "A", "C", "WYQHP", "GTSN", "RK", "DE"],
                8: ["MFILV", "A", "C", "WYQHP", "G", "TSN", "RK", "DE"],
                9: ["MF", "ILV", "A", "C", "WYQHP", "G", "TSN", "RK", "DE"],
                10: [
                    "MF",
                    "ILV",
                    "A",
                    "C",
                    "WYQHP",
                    "G",
                    "TSN",
                    "RK",
                    "D",
                    "E",
                ],
                11: [
                    "MF",
                    "IL",
                    "V",
                    "A",
                    "C",
                    "WYQHP",
                    "G",
                    "TSN",
                    "RK",
                    "D",
                    "E",
                ],
                12: [
                    "MF",
                    "IL",
                    "V",
                    "A",
                    "C",
                    "WYQHP",
                    "G",
                    "TS",
                    "N",
                    "RK",
                    "D",
                    "E",
                ],
                13: [
                    "MF",
                    "IL",
                    "V",
                    "A",
                    "C",
                    "WYQHP",
                    "G",
                    "T",
                    "S",
                    "N",
                    "RK",
                    "D",
                    "E",
                ],
                14: [
                    "MF",
                    "I",
                    "L",
                    "V",
                    "A",
                    "C",
                    "WYQHP",
                    "G",
                    "T",
                    "S",
                    "N",
                    "RK",
                    "D",
                    "E",
                ],
                15: [
                    "MF",
                    "IL",
                    "V",
                    "A",
                    "C",
                    "WYQ",
                    "H",
                    "P",
                    "G",
                    "T",
                    "S",
                    "N",
                    "RK",
                    "D",
                    "E",
                ],
                16: [
                    "MF",
                    "I",
                    "L",
                    "V",
                    "A",
                    "C",
                    "WYQ",
                    "H",
                    "P",
                    "G",
                    "T",
                    "S",
                    "N",
                    "RK",
                    "D",
                    "E",
                ],
                20: [
                    "M",
                    "F",
                    "I",
                    "L",
                    "V",
                    "A",
                    "C",
                    "W",
                    "Y",
                    "Q",
                    "H",
                    "P",
                    "G",
                    "T",
                    "S",
                    "N",
                    "R",
                    "K",
                    "D",
                    "E",
                ],
            }
            fastas = self.seq_list
            subtype = self.__default_para["PseKRAAC_model"]
            raactype = self.__default_para["RAAC_clust"]
            ktuple = self.__default_para["k-tuple"]
            glValue = (
                self.__default_para["g-gap"]
                if subtype == "g-gap"
                else self.__default_para["lambdaValue"]
            )

            if raactype not in AAGroup:
                self.error_msg = "PseKRAAC descriptor value error."
                return False

            # index each amino acids to their group
            myDict = {}
            for i in range(len(AAGroup[raactype])):
                for aa in AAGroup[raactype][i]:
                    myDict[aa] = i

            gDict = {}
            gNames = []
            for i in range(len(AAGroup[raactype])):
                gDict[i] = "T1.G." + str(i + 1)
                gNames.append("T1.G." + str(i + 1))

            encodings = []
            if subtype == "g-gap":
                encodings = self.gapModel(
                    fastas, myDict, gDict, gNames, ktuple, glValue, "type15"
                )
            else:
                encodings = self.lambdaModel(
                    fastas, myDict, gDict, gNames, ktuple, glValue, "type15"
                )
            encodings = np.array(encodings)
            self.encodings = pd.DataFrame(
                encodings[1:, 1:].astype(float),
                columns=encodings[0, 1:],
                index=encodings[1:, 0],
            )
            return True
        except Exception as e:
            sys.exit(f"Error: {e}")

    def _PseKRAAC_type_16(self):
        try:
            AAGroup = {
                2: ["IMVLFWY", "GPCASTNHQEDRK"],
                3: ["IMVLFWY", "GPCAST", "NHQEDRK"],
                4: ["IMVLFWY", "G", "PCAST", "NHQEDRK"],
                5: ["IMVL", "FWY", "G", "PCAST", "NHQEDRK"],
                6: ["IMVL", "FWY", "G", "P", "CAST", "NHQEDRK"],
                7: ["IMVL", "FWY", "G", "P", "CAST", "NHQED", "RK"],
                8: ["IMV", "L", "FWY", "G", "P", "CAST", "NHQED", "RK"],
                9: ["IMV", "L", "FWY", "G", "P", "C", "AST", "NHQED", "RK"],
                10: [
                    "IMV",
                    "L",
                    "FWY",
                    "G",
                    "P",
                    "C",
                    "A",
                    "STNH",
                    "RKQE",
                    "D",
                ],
                11: [
                    "IMV",
                    "L",
                    "FWY",
                    "G",
                    "P",
                    "C",
                    "A",
                    "STNH",
                    "RKQ",
                    "E",
                    "D",
                ],
                12: [
                    "IMV",
                    "L",
                    "FWY",
                    "G",
                    "P",
                    "C",
                    "A",
                    "ST",
                    "N",
                    "HRKQ",
                    "E",
                    "D",
                ],
                13: [
                    "IMV",
                    "L",
                    "F",
                    "WY",
                    "G",
                    "P",
                    "C",
                    "A",
                    "ST",
                    "N",
                    "HRKQ",
                    "E",
                    "D",
                ],
                14: [
                    "IMV",
                    "L",
                    "F",
                    "WY",
                    "G",
                    "P",
                    "C",
                    "A",
                    "S",
                    "T",
                    "N",
                    "HRKQ",
                    "E",
                    "D",
                ],
                15: [
                    "IMV",
                    "L",
                    "F",
                    "WY",
                    "G",
                    "P",
                    "C",
                    "A",
                    "S",
                    "T",
                    "N",
                    "H",
                    "RKQ",
                    "E",
                    "D",
                ],
                16: [
                    "IMV",
                    "L",
                    "F",
                    "W",
                    "Y",
                    "G",
                    "P",
                    "C",
                    "A",
                    "S",
                    "T",
                    "N",
                    "H",
                    "RKQ",
                    "E",
                    "D",
                ],
                20: [
                    "I",
                    "M",
                    "V",
                    "L",
                    "F",
                    "W",
                    "Y",
                    "G",
                    "P",
                    "C",
                    "A",
                    "S",
                    "T",
                    "N",
                    "H",
                    "R",
                    "K",
                    "Q",
                    "E",
                    "D",
                ],
            }
            fastas = self.seq_list
            subtype = self.__default_para["PseKRAAC_model"]
            raactype = self.__default_para["RAAC_clust"]
            ktuple = self.__default_para["k-tuple"]
            glValue = (
                self.__default_para["g-gap"]
                if subtype == "g-gap"
                else self.__default_para["lambdaValue"]
            )

            if raactype not in AAGroup:
                self.error_msg = "PseKRAAC descriptor value error."
                return False

            # index each amino acids to their group
            myDict = {}
            for i in range(len(AAGroup[raactype])):
                for aa in AAGroup[raactype][i]:
                    myDict[aa] = i

            gDict = {}
            gNames = []
            for i in range(len(AAGroup[raactype])):
                gDict[i] = "T1.G." + str(i + 1)
                gNames.append("T1.G." + str(i + 1))

            encodings = []
            if subtype == "g-gap":
                encodings = self.gapModel(
                    fastas, myDict, gDict, gNames, ktuple, glValue, "type16"
                )
            else:
                encodings = self.lambdaModel(
                    fastas, myDict, gDict, gNames, ktuple, glValue, "type16"
                )
            encodings = np.array(encodings)
            self.encodings = pd.DataFrame(
                encodings[1:, 1:].astype(float),
                columns=encodings[0, 1:],
                index=encodings[1:, 0],
            )
            return True
        except Exception as e:
            sys.exit(f"Error: {e}")

    def to_csv(self, file="encodings.csv", index=False, header=False):
        try:
            self.encodings.to_csv(file, index=index, header=header)
        except Exception as e:
            sys.exit(f"Error: {e}")
        return True

    def to_tsv(self, file="encodings.tsv", index=False, header=False):
        try:
            self.encodings.to_csv(file, sep="\t", index=index, header=header)
        except Exception as e:
            sys.exit(f"Error: {e}")
        return True

    def to_svm(self, file="encodings.svm"):
        try:
            with open(file, "w") as f:
                for line in self.encodings.values:
                    f.write("1")
                    for i in range(len(line)):
                        f.write("  %d:%s" % (i + 1, line[i]))
                    f.write("\n")
            return True
        except Exception as e:
            sys.exit(f"Error: {e}")

    def to_arff(self, file="encodings.arff"):
        with open(file, "w") as f:
            f.write("@relation descriptor\n\n")
            for i in range(1, len(self.encodings.values[0]) + 1):
                f.write("@attribute f.%d numeric\n" % i)
            f.write("@attribute play {yes, no}\n\n")
            f.write("@data\n")
            for line in self.encodings.values:
                line = line
                for fea in line:
                    f.write("%s," % fea)
                f.write("yes\n")
