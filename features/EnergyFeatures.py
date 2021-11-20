import tempfile
from pathlib import Path
from subprocess import call
import RNA
import re

from duplex.Duplex import Duplex
from features.Features import Features
from utils.utils import get_subsequence_by_coordinates


class EnergyFeatures(Features):
    def extract_features(self):
        self._features_dict["minimum_free_energy"] = self.minimum_free_energy()

    def minimum_free_energy(self):
        mrna_full: str = self._region_sequence
        mrna_surrounding100: str = get_subsequence_by_coordinates(mrna_full, self._start, self._end, extra_chars=50)
        seed: Duplex = self._duplex.seed
        mrna_seed = seed.site
        mrna_site = self._duplex.site[::-1]
        mrna_site_3p = mrna_site[len(mrna_seed):]



        MFE_Seed = RNA.fold(mrna_seed)
        MFE_surrounding100 = RNA.fold(mrna_surrounding100)
        MFE_surrounding100_norm = MFE_surrounding100[1] / len(mrna_surrounding100)
        duplex = RNA.duplexfold(self._miRNA_sequence, self._site)
        MEF_duplex = duplex.energy

        if mrna_site_3p=="":
            MFE_3p = (0, 0)
        else:
            MFE_3p = RNA.fold(mrna_site_3p)

        constraint_low = "."*min(self._start - 1, 50)
        constraint_site = "x"*min(len(mrna_site), len(self._region_sequence) - self._start + 1)
        constraint_high = "."*min(len(self._region_sequence) - self._end, 50)
        constraint = constraint_low + constraint_site + constraint_high
        # print(len(self._region_sequence))
        # print(self._start)
        # print(len(self._region_sequence) - self._start + 1)
        # print(len(mrna_site))
        # print(len(constraint_site))
        assert len(constraint) == len(mrna_surrounding100), \
            f"""constraint and mrna_surrounding100 are not in the same length
mrna_site: {mrna_site}
constraint:          {constraint}
mrna_surrounding100: {mrna_surrounding100}"""

        with tempfile.TemporaryDirectory() as tmpdirname:
            tmp_dir = Path(tmpdirname)

            mrna100file_in = tmp_dir / 'mrna100_with_constraints.fa'
            mrna100file_out = tmp_dir /'mrna100_with_constraints.result'
            with mrna100file_in.open('w') as f:
                f.write(mrna_surrounding100 + "\n" + constraint + "\n")

            cmd = "RNAfold -C {infile} > {outfile}".format(infile=mrna100file_in, outfile=mrna100file_out)
            status = call(cmd, cwd=tmp_dir.resolve(), shell=True)

            with mrna100file_out.open() as f:
                twolines = f.readlines()
            cmfe =float(re.findall("\d+\.\d+", twolines[1])[0])*(-1)

        MFE = {'Energy_MEF_Seed': round(MFE_Seed[1], 4),
               'Energy_MEF_3p': round(MFE_3p[1], 4),
               'Energy_MEF_local_target': round(MFE_surrounding100[1], 4),
               'Energy_MEF_local_target_normalized': round(MFE_surrounding100_norm, 4),
               'Energy_MEF_Duplex': round(MEF_duplex, 4),
               'Energy_MEF_cons_local_target': round(cmfe, 4),
               'Energy_MEF_cons_local_target_normalized': round(cmfe/len(mrna_surrounding100), 4)
               }
        return MFE


