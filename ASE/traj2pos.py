#!/usr/bin/env python

from ase.io import read,write
from ase.build import sort
import sys

rin=sys.argv[1]

p=read(rin)
p=sort(p)

write('POSCAR',p)

