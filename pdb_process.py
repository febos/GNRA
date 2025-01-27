
import os, glob

from urslib2 import SplitmmCIF


if __name__ == "__main__":

    source   = "PDB"
    firsts   = "PDB1"
    crystals = "PDB1cc"

    os.makedirs(firsts,   exist_ok = True)
    os.makedirs(crystals, exist_ok = True)

    #SplitmmCIF.All(source,firsts, m1only = True)

    cnt = 0

    for first in sorted(glob.glob(os.path.join(firsts,"*.cif1"))):

        crystal = os.path.join(crystals,os.path.basename(first)[:-1])
        cnt += 1

        xray = False
        with open(first) as inp:
            for line in inp:
                if line.startswith("_exptl.method") and "X-RAY DIFFRACTION" in line.upper():
                    xray = True
                    break

        # Don't derive crystal contacts for files larger 10MB (~ribosomes)
        if os.path.getsize(first) > 10**7:
            xray = False
                
        if os.path.exists(crystal) and os.path.getsize(crystal) > 0:
            print(cnt, first, '->', end=' ')
            print(crystal, "(skipped)")
        else:
            print(cnt, first, '->', end=' ')
            os.system('echo "open {} format mmcif; {}'.format(first,
                                                              "crystalcontacts #1 distance 6;"*xray)+\
                      ' save {} format mmcif" | chimerax --nogui'.format(crystal))
            print(crystal)

        # COPY THE MISSING ENTRIES, e.g. 180D 181D
        if not os.path.exists(crystal):
            os.system("cp {} {}".format(first,crystal))
