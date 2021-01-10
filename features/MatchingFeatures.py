from abc import ABC
from itertools import product
from typing import Dict

from pandas import Series

from consts.features import HOT_ENCODING_LEN
from features.Features import Features



class MatchingFeatures(Features):
    def extract_features(self):
        self._features_dict["miRNA_match_position"] = self.miRNA_match_position()
        self._features_dict["miRNA_pairing_count"] = self.miRNA_pairing_count()


    def miRNA_match_position(self):  # 20
        mmp_dic = {}
        pair_list = list(self._duplex.mir_pairing_iterator())
        for i in range (1, 21):
            key = 'miRNAMatchPosition_' + str(i)
            try:
                pair = pair_list[i-1]
            except IndexError:
                pair = '  '

            if pair in AU:
                mmp_dic[key] = 2
            elif pair in GC:
                mmp_dic[key] = 1
            elif pair in GU:
                mmp_dic[key] = 3
            elif ' ' in pair:
                mmp_dic[key] = 5
            else:
                mmp_dic[key] = 4

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

        i = 0
        for pair in self._duplex.mir_pairing_iterator():
            i += 1
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

