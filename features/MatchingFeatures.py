from abc import ABC
from itertools import product
from typing import Dict

from pandas import Series

from consts.features import AU, GC, GU, HOT_ENCODING_LEN, MM
from duplex.utils import mix_inter_bulge_seq
from features.Features import Features



class MatchingFeatures(Features):
    def extract_features(self):
        self._features_dict["miRNA_match_position"] = self.miRNA_match_position()
        self._features_dict["miRNA_pairing_count"] = self.miRNA_pairing_count()


    def miRNA_match_position(self):  # 20
        mmp_dic = {}
        # print(list(self._duplex.mir_iterator()))
        # pair_list = [mir+mix_inter_bulge_seq(self._duplex._mrna_inter[i], self._duplex._mrna_bulge[i])
        #              for i, mir in self._duplex.mir_iterator()] # bug ?
        basepair_list = [mir + self._duplex._mrna_inter[i] for i, mir in self._duplex.mir_iterator()]
        misspair_list = [mir + self._duplex._mrna_bulge[i] for i, mir in self._duplex.mir_iterator()]

        # print(basepair_list)
        # print(misspair_list)

        for i in range (22):
            key = 'miRNAMatchPosition_' + str(i+1)
            try:
                basepair = basepair_list[i] #bug: was i-1
                misspair = misspair_list[i]
                # print(f"{i}  {pair}")

                if basepair in AU:
                    mmp_dic[key] = "AU"
                elif basepair in GC:
                    mmp_dic[key] = "GC"
                elif basepair in GU:
                    mmp_dic[key] = "GU"
                elif ' ' in basepair:
                    if ' ' in misspair:
                        mmp_dic[key] = "BB"
                    else:
                        mmp_dic[key] = "MM"
                else:
                    raise Exception(f"It shouldn't be here. i={i}, basepair={basepair} misspair={misspair} "
                                    f"{self._duplex}")
            except IndexError:
                mmp_dic[key] = ""

        return mmp_dic

    def miRNA_pairing_count(self):  # 6*3=18
        mpc_dic = {'miRNAPairingCount_Seed_GC': 0,
                   'miRNAPairingCount_Seed_AU': 0,
                   'miRNAPairingCount_Seed_GU': 0,
                   'miRNAPairingCount_Seed_mismatch': 0,
                   'miRNAPairingCount_Seed_bulge': 0,
                   'miRNAPairingCount_Seed_bulge_nt': 0,
                   'miRNAPairingCount_Total_GC': 0,
                   'miRNAPairingCount_Total_AU': 0,
                   'miRNAPairingCount_Total_GU': 0,
                   'miRNAPairingCount_Total_mismatch': 0,
                   'miRNAPairingCount_Total_bulge': 0,
                   'miRNAPairingCount_Total_bulge_nt': 0,
                   'miRNAPairingCount_X3p_GC': 0,
                   'miRNAPairingCount_X3p_AU': 0,
                   'miRNAPairingCount_X3p_GU': 0,
                   'miRNAPairingCount_X3p_mismatch': 0,
                   'miRNAPairingCount_X3p_bulge': 0,
                   'miRNAPairingCount_X3p_bulge_nt': 0}

        i = 0 # seed index
        for k in range(len(self._duplex.mir_bulge)):
            mir = mix_inter_bulge_seq(self._duplex.mir_bulge[k], self._duplex.mir_inter[k])
            mrna = mix_inter_bulge_seq(self._duplex.mrna_bulge[k], self._duplex.mrna_inter[k])
            pair = mir + mrna


            if mir != "":
                i += 1

            if i == 0:
                continue

            if pair in AU:
                mpc_dic['miRNAPairingCount_Total_AU'] += 1
                if 0 < i < 9:
                    mpc_dic['miRNAPairingCount_Seed_AU'] += 1
                if i >=9 :
                    mpc_dic['miRNAPairingCount_X3p_AU'] += 1
            elif pair in GC:
                mpc_dic['miRNAPairingCount_Total_GC'] += 1
                if 0 < i < 9:
                    mpc_dic['miRNAPairingCount_Seed_GC'] += 1
                if i >= 9:
                    mpc_dic['miRNAPairingCount_X3p_GC'] += 1
            elif pair in GU:
                mpc_dic['miRNAPairingCount_Total_GU'] += 1
                if 0 < i < 9:
                    mpc_dic['miRNAPairingCount_Seed_GU'] += 1
                if i >= 9:
                    mpc_dic['miRNAPairingCount_X3p_GU'] += 1
            elif pair in MM:
                mpc_dic['miRNAPairingCount_Total_mismatch'] += 1
                if 0 < i < 9:
                    mpc_dic['miRNAPairingCount_Seed_mismatch'] += 1
                if i >= 9:
                    mpc_dic['miRNAPairingCount_X3p_mismatch'] += 1
            elif len(pair.strip()) == 1:
                mpc_dic['miRNAPairingCount_Total_bulge_nt'] += 1
                if 0 < i < 9:
                    mpc_dic['miRNAPairingCount_Seed_bulge_nt'] += 1
                if i >= 9:
                    mpc_dic['miRNAPairingCount_X3p_bulge_nt'] += 1

        mpc_dic['miRNAPairingCount_Total_bulge'] = self._duplex.mir_bulge_count
        mpc_dic['miRNAPairingCount_Seed_bulge'] = self._duplex.seed.mir_bulge_count
        mpc_dic['miRNAPairingCount_X3p_bulge'] = mpc_dic['miRNAPairingCount_Total_bulge'] - mpc_dic['miRNAPairingCount_Seed_bulge']

        return mpc_dic

