#!/usr/bin/env python
import numpy as np

from ase import io
from ase import Atoms
from espresso import espresso
from ase.optimize import QuasiNewton
from ase.constraints import FixAtoms
import sys
import pickle

slab =  io.read('opt.traj')
slab.set_pbc([True,True,True])


calc = espresso(pw=500,             #plane-wave cutoff
                dw=5000,                    #density cutoff
                xc='PBE',          #exchange-correlation functional
                kpts=(5,5,1),       #k-point sampling;
                nbands=-20,             #20 extra bands besides the bands needed to hold
                sigma=0.1,
                #mode = 'vc-relax',
                #cell_dynamics = 'bfgs',
                #opt_algorithm = 'bfgs',
                #fmax = 0.05,
		nosym=True,
                convergence= {'energy':1e-5,
                    'mixing':0.1,
                    'nmix':10,
                    'mix':4,
                    'maxsteps':500,
                    'diag':'david'
                    },  #convergence parameters
                 dipole={'status':True}, #dipole correction to account for periodicity in z
                 output = {'avoidio':False,
                    'removewf':True,
                    'wf_collect':False},
                 spinpol=False,
                 parflags='-npool 2',
                 outdir='calcdir')   #output directory for Quantum Espresso files

#calc.calculation_required = lambda x, y: True
slab.set_calculator(calc)
slab.get_potential_energy()
traj=io.Trajectory('opt.traj','w')
traj.write(slab)
#qn=QuasiNewton(slab,trajectory='opt.traj',logfile='opt.log')
#qn.run(fmax=0.03)
calc.save_flev_chg('efchg.tgz')
potential = calc.extract_total_potential()

potential_file=open('potential.pickle','w')
pickle.dump(potential,potential_file)
potential_file.close()

posin= io.read('opt.traj')
p=posin.copy()
p.calc=calc
p.calc.load_flev_chg('efchg.tgz')
dos = calc.calc_pdos(nscf=True,
		     kpts=(10,10,1),
		     Emin=-15.0,
		     Emax=15.0,
	             DeltaE=0.01,
                     sigma = 0.1, 
                     ngauss =0, 
                     tetrahedra=False,
                     slab=True)
f = open('dos.pickle', 'w')
pickle.dump(dos, f)
f.close()
