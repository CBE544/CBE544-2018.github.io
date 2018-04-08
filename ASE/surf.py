#!/usr/bin/env python
import numpy as np

from ase import io
from espresso import espresso
from ase.build import surface
from ase.geometry import get_layers
import sys
#from ase.optimize import BFGS
#sys.stdout = open('output.dat', 'w')
atoms =  io.read('SrTiO3-conv.cif')
a=3.936
cell=atoms.set_cell([(a,0,0),(0,a,0), (0,0,a)], scale_atoms=True)
traj=io.Trajectory('sto001.traj','w')
#traj.write(atoms)
s1=surface(atoms, (0,0,1), 4)
s1.center(vacuum=10, axis=2)

traj.write(s1)
#(lay,dist)=get_layers(s1, (0,0,1))
#sym=atoms.get_chemical_symbols()
#print("atom")
#print(sym)
#print ("distance between layers")
#np.set_printoptions(precision=3)
#print(dist)
#print ("atom layer")
#print(lay)



