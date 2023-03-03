from ase.build import fcc111, surface, sort, niggli_reduce, add_adsorbate
from ase.calculators.emt import EMT
from ase.db import connect
from ase.io import read, write
from ase.visualize import view
from ase import Atoms
import os
import numpy as np
import collections
import aseplus.tools as tools
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--cif", default=None, type=str)
parser.add_argument("--adsorbate_smiles", default=None, type=str, help="smiles string for adsorbate")
parser.add_argument("--rotate", default=None, type=str, help="x|y|z_degree")
parser.add_argument("--rotate2", default=None, type=str, help="x|y|z_degree")
parser.add_argument("--height", default=None, type=float)
parser.add_argument("--nlayer", default=3, type=int)
parser.add_argument("--vacuum", default=10.0, type=float)
parser.add_argument("--basedir", default="", type=str)
parser.add_argument("--workdir", default="work", type=str)
parser.add_argument("--worksubdir", default="", type=str)

args = parser.parse_args()

cif = args.cif
adsorbate_smiles = args.adsorbate_smiles

workdir = os.path.join(args.basedir, args.workdir, args.worksubdir)
if not os.path.isdir(workdir):
    os.makedirs(workdir)
os.chdir(workdir)

if args.rotate is None:
    rotate_dir_and_angle = ["x", 0]
else:
    rotate_dir_and_angle = [args.rotate[0], args.rotate[2:]]

if args.rotate2 is None:
    rotate_dir_and_angle2 = ["x", 0]
else:
    rotate_dir_and_angle2 = [args.rotate2[0], args.rotate2[2:]]

if args.height is None:
    height = 4.0
else:
    height = args.height

nlayer = args.nlayer
vacuum = args.vacuum

lattice = "fcc"
facet   = "111"
#lattice = "hcp"
#lattice = "sp15"
#facet = "010"

indices = []
for c in facet:
    indices.append(int(c))
bulk = read(cif)
surf = surface(lattice=bulk, indices=indices, layers=nlayer, vacuum=vacuum, periodic=True)

surf = surf*[3, 3, 1]
#surf = surf*[3, 2, 1]
surf = sort(surf)
surf = tools.sort_atoms_by(surf, xyz="z")

formula = surf.get_chemical_formula()

offset_fac = (2.1, 1.5)
offset_fac = np.array(offset_fac)

surf.translate([0, 0, -vacuum+0.5])
surf = tools.fix_lower_surface(surf)
#
# prepare adsorbate
#
os.system('obabel -:"{0:s}" -oxyz -h --gen3D -O tmp.xyz'.format(adsorbate_smiles))
adsorbate = read("tmp.xyz")
adsorbate.set_tags([-1]*len(adsorbate))
adsorbate.center()
adsorbate.rotate(v=rotate_dir_and_angle[0], a=int(rotate_dir_and_angle[1]))
adsorbate.rotate(v=rotate_dir_and_angle2[0], a=int(rotate_dir_and_angle2[1]))

# shift
min_ind = min(adsorbate.positions[2, :])
adsorbate.translate([0, 0, -min_ind])

#
# adsorb on surface
#
if adsorbate is not None:
    offset = (0.16, 0.33)  # x1y1
    offset = np.array(offset)
    add_adsorbate(surf, adsorbate, height=height, position=(0, 0), offset=offset*offset_fac)

write("POSCAR", surf)

db = connect("surf_and_ads.json")
db.write(surf, data={"smiles": adsorbate_smiles})

