from ase import Atoms
from ase.io import read, write
from ase.calculators.vasp import Vasp
from ase.optimize.bfgs import BFGS
import numpy as np


def get_vasp_setting():
    prec   = "normal"
    encut  = 400
    xc     = "pbe"
    ivdw   = 0
    nsw    = 0  # will be overwritten by steps
    nelm   = 30
    nelmin = 3
    ibrion = -1
    potim  = 0.2
    algo   = "Fast"  # sometimes VeryFast fails
    ismear = 0
    sigma  = 0.1
    ediff  = 1.0e-4
    ediffg = -0.1
    kpts   = [1, 1, 1]
    ispin  = 1
    pp     = "potpaw_PBE.54"
    npar   = 6
    nsim   = npar
    isym   = 0
    lreal  = True
    lorbit = 10  # to avoid error
    lwave  = False
    lcharg = False

    ldipol = True
    idipol = 3

    calc = Vasp(prec=prec, encut=encut, xc=xc, ivdw=ivdw, algo=algo, ediff=ediff, ediffg=ediffg,
                ibrion=ibrion, potim=potim, nsw=nsw, nelm=nelm, nelmin=nelmin, kpts=kpts, ismear=ismear,
                kgamma=True, ispin=ispin, pp=pp, npar=npar, nsim=nsim, isym=isym, lreal=lreal, lwave=lwave,
                lcharg=lcharg, sigma=sigma, lorbit=lorbit, ldipol=ldipol, idipol=idipol)

    return calc


surf_and_ads = read("surf_plus_ads.db")

ads  = list(filter(lambda atom: atom.tag == -1, surf_and_ads))
surf = list(filter(lambda atom: atom.tag != -1, surf_and_ads))
ads  = Atoms(ads)
surf = Atoms(surf)

E = np.zeros(3)

steps = 10
for i, mol in enumerate([ads, surf, surf_and_ads]):
    calc = get_vasp_setting()

    # optimization
    calc.int_params["ibrion"] = 2
    calc.int_params["nsw"] = steps

    mol.set_calculator(calc)

    #E[i] = mol.get_potential_energy()

Eads = E[2] - (E[0] + E[1])
print("Adsorption energy (eV) = {}".format(Eads))

