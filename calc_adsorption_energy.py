from ase import Atoms
from ase.io import read, write
from ase.calculators.vasp import Vasp
from ase.optimize.bfgs import BFGS
from ase.db import connect
import numpy as np
import os
import argparse
import json


def get_vasp_setting(directory=None, gas_phase=False):
    prec   = "normal"
    encut  = 400
    xc     = "pbe"
    ivdw   = 12
    nsw    = 0  # will be overwritten by steps
    nelm   = 40
    nelmin = 5
    ibrion = -1
    potim  = 0.2
    algo   = "VeryFast"  # sometimes VeryFast fails
    ismear = 0
    sigma  = 0.1
    ediff  = 1.0e-5
    ediffg = -3.0e-2
    kpts   = [1, 1, 1]
    ispin  = 1
    lasph  = True
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

    if gas_phase:
        calc = Vasp(directory=directory, prec=prec, encut=encut, xc=xc, ivdw=ivdw, algo=algo, ediff=ediff, ediffg=ediffg,
                    ibrion=ibrion, potim=potim, nsw=nsw, nelm=nelm, nelmin=nelmin, kpts=[1, 1, 1], ismear=0,
                    ispin=ispin, pp=pp, npar=npar, nsim=nsim, isym=isym, lreal=True, lwave=lwave,
                    lcharg=lcharg, sigma=sigma, lorbit=lorbit, lasph=lasph)
    else:
        calc = Vasp(directory=directory, prec=prec, encut=encut, xc=xc, ivdw=ivdw, algo=algo, ediff=ediff, ediffg=ediffg,
                    ibrion=ibrion, potim=potim, nsw=nsw, nelm=nelm, nelmin=nelmin, kpts=kpts, ismear=ismear,
                    ispin=ispin, pp=pp, npar=npar, nsim=nsim, isym=isym, lreal=lreal, lwave=lwave,
                    lcharg=lcharg, sigma=sigma, lorbit=lorbit, lasph=lasph, ldipol=ldipol, idipol=idipol)

    return calc

parser = argparse.ArgumentParser()
parser.add_argument("--basedir", default=None, type=str)
parser.add_argument("--workdir", default="work", type=str)
parser.add_argument("--worksubdir", default="", type=str)
parser.add_argument("--jsonfile", default="result.json", type=str)
parser.add_argument("--steps", default=200, type=int)
args = parser.parse_args()

# directory
if args.basedir is None:
    basedir = os.getcwd()
else:
    basedir = args.basedir

workdir = os.path.join(basedir, args.workdir, args.worksubdir)
if not os.path.isdir(workdir):
    os.makedirs(workdir)
os.chdir(workdir)

# json file to write output
jsonfile = os.path.join(basedir, args.jsonfile)
if not os.path.isfile(jsonfile):
    with open(jsonfile, "w") as f:
        f.write("")

db = connect("surf_and_ads.json")
surf_and_ads = db.get_atoms(id=1)
row = db.get(id=1)
smiles = row.data.smiles
anchor = row.data.anchor_atom

ads  = list(filter(lambda atom: atom.tag == -1, surf_and_ads))
surf = list(filter(lambda atom: atom.tag != -1, surf_and_ads))
ads  = Atoms(ads)
surf = Atoms(surf)

vacuum = 15.0
ads.set_cell([(vacuum, 0, 0), (0, vacuum, 0), (0, 0, vacuum)])
ads.center()
ads.set_pbc(True)

surf.cell = surf_and_ads.cell
surf.set_pbc(True)

constraints = surf_and_ads.constraints
surf.set_constraint(constraints)

E = np.zeros(3)

for i, mol in enumerate([ads, surf, surf_and_ads]):
    directory = os.path.join(workdir, mol.get_chemical_formula())

    if i == 0:
        gas_phase = True
    else:
        gas_phase = False

    calc = get_vasp_setting(directory=directory, gas_phase=gas_phase)

    # optimization
    calc.int_params["ibrion"] = 2
    calc.int_params["nsw"] = args.steps

    mol.set_calculator(calc)

    print("calculating {:s}".format(mol.get_chemical_formula()), flush=True)

    E[i] = mol.get_potential_energy()

Eads = E[2] - (E[0] + E[1])
print("Adsorption energy (eV) = {}".format(Eads))

data = {"adsorbate_formula": ads.get_chemical_formula(), "adsorbate_smiles": smiles,
        "anchor_atom": anchor, "adsorption_energy": Eads}
        
with open(jsonfile, "a") as f:
    json.dump(data, f, indent=4)

