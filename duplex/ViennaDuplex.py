from typing import Tuple

import RNA

from duplex.Duplex import Duplex


def find_pairing(s, ch):
    return [i for i, ltr in enumerate(s) if ltr == ch]




class ViennaDuplex(Duplex):
    @classmethod
    def createDuplex(cls, mirna: str, target: str) -> Tuple[str, str, str, str]:
        duplex = RNA.duplexfold(mirna, target)

        (mir_pairing, mrna_pairing) = duplex.structure.split('&')
        mir_coor = (duplex.i - len(mir_pairing), duplex.i)
        mrna_coor = (duplex.j - 1, duplex.j + len(mrna_pairing) - 1)
        active_mir = mirna[mir_coor[0]:mir_coor[1]]
        active_mrna = target[mrna_coor[0]:mrna_coor[1]]

        mir_idx = find_pairing(mir_pairing, '(')
        mrna_idx = find_pairing(mrna_pairing, ')')
        mrna_idx = mrna_idx[::-1]
        mrna_len = len(active_mrna)

        mrna_bulge = ""
        mrna_inter = ""
        mir_inter = ""
        mir_bulge = ""

        mir_i = 0
        mrna_i = mrna_len - 1
        if (mir_coor[0] > 0):
            mir_bulge += mirna[:mir_coor[0]]
            mir_inter += " " * mir_coor[0]
            mrna_inter += " " * mir_coor[0]
            mrna_addition_len = mrna_coor[1] + mir_coor[0] - mrna_coor[1]
            mrna_bulge_additon = target[mrna_coor[1]:mrna_coor[1] + mir_coor[0]]
            mrna_bulge_additon = mrna_bulge_additon + "#" * (mrna_addition_len - len(mrna_bulge_additon))
            mrna_bulge += mrna_bulge_additon[::-1]

        for i in range(len(mir_idx)):
            # deal with the bulge
            mir_bulge_idx = range(mir_i, mir_idx[i])
            mir_bulge += active_mir[mir_i:mir_idx[i]]
            mrna_bulge_idx = range(mrna_i, mrna_idx[i], -1)
            mrna_bulge += active_mrna[mrna_i:mrna_idx[i]: -1]
            c_pos = max(len(mrna_bulge_idx), len(mir_bulge_idx))
            mrna_inter += " " * c_pos
            mir_inter += " " * c_pos
            mrna_bulge += " " * (c_pos - len(mrna_bulge_idx))
            mir_bulge += " " * (c_pos - len(mir_bulge_idx))

            # deal with the interaction
            mir_bulge += " "
            mir_inter += active_mir[mir_idx[i]]
            mrna_bulge += " "
            mrna_inter += active_mrna[mrna_idx[i]]
            # update the idx
            mir_i = mir_idx[i] + 1
            mrna_i = mrna_idx[i] - 1


        mir_i += mir_coor[0]
        mir_bulge_additon=mirna[mir_i:]
        mir_bulge+=mir_bulge_additon

        mrna_addition = target[max(0, mrna_coor[0]-len(mir_bulge_additon)+mrna_i+1):mrna_coor[0]+mrna_i+1]
        mrna_bulge += mrna_addition[::-1]



        # mir_i += mir_coor[0] # new. bug fix
        # print(mir_i)
        # if (mir_i <= len(mirna)):
        #     mir_bulge += mirna[mir_i:]
        #     addition = mirna[mir_i:]
        #     # s = max(0, mrna_coor[0] + 1 - len(addition))
        #     # e = mrna_coor[0] + 1
        #     s = mrna_i - len(addition)
        #     e = mrna_i
        #     mrna_bulge += target[s:e][::-1]

        return mrna_bulge, mrna_inter, mir_inter, mir_bulge

if __name__ == '__main__':
    mirna = 'UGAGGUAGUAGGUUGUAUAGUU'
    target = 'CAACAAGAACCAACACUACUGCCCAACCGUCAACGUCG'
    dp = ViennaDuplex.fromChimera(mirna, target)
    print(dp)


