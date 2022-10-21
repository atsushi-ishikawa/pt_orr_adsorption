from ase import Atoms
from ase.io import read, write
from ase.calculators.vasp import Vasp
from ase.optimize.bfgs import BFGS
import numpy as np
import os


def get_vasp_setting(directory=None):
    prec   = "normal"
    encut  = 400
    xc     = "pbe"
    ivdw   = 12
    nsw    = 0  # will be overwritten by steps
    nelm   = 30
    nelmin = 3
    ibrion = -1
    potim  = 0.2
    algo   = "VeryFast"  # sometimes VeryFast fails
    ismear = 0
    sigma  = 0.1
    ediff  = 1.0e-4
    ediffg = -5.0e-2
    kpts   = [1, 1, 1]
    ispin  = 1
    pp     = "potpaw_PBE.54"
    npar   = 10
    nsim   = npar
    isym   = 0
    lreal  = True
    lorbit = 10  # to avoid error
    lwave  = False
    lcharg = False

    ldipol = True
    idipol = 3

    calc = Vasp(directory=directory, prec=prec, encut=encut, xc=xc, ivdw=ivdw, algo=algo, ediff=ediff, ediffg=ediffg,
                ibrion=ibrion, potim=potim, nsw=nsw, nelm=nelm, nelmin=nelmin, kpts=kpts, ismear=ismear,
                kgamma=True, ispin=ispin, pp=pp, npar=npar, nsim=nsim, isym=isym, lreal=lreal, lwave=lwave,
                lcharg=lcharg, sigma=sigma, lorbit=lorbit, ldipol=ldipol, idipol=idipol)

    return calc

workdir = "work"
os.makedirs(os.path.join(os.getcwd(), workdir), exist_ok=True)

surf_and_ads = read("surf_plus_ads.db")

ads  = list(filter(lambda atom: atom.tag == -1, surf_and_ads))
surf = list(filter(lambda atom: atom.tag != -1, surf_and_ads))
ads  = Atoms(ads)
surf = Atoms(surf)

vacuum = 10.0
ads.set_cell([(vacuum, 0, 0), (0, vacuum, 0), (0, 0, vacuum)])
ads.center()
ads.set_pbc(True)

surf.cell = surf_and_ads.cell
surf.set_pbc(True)

constraints = surf_and_ads.constraints
surf.set_constraint(constraints)

E = np.zeros(3)

steps = 100
for i, mol in enumerate([ads, surf, surf_and_ads]):
    directory = os.path.join(workdir, mol.get_chemical_formula())
    calc = get_vasp_setting(directory=directory)

    # optimization
    calc.int_params["ibrion"] = 2
    calc.int_params["nsw"] = steps

    mol.set_calculator(calc)

    E[i] = mol.get_potential_energy()

Eads = E[2] - (E[0] + E[1])
print("Adsorption energy (eV) = {}".format(Eads))

