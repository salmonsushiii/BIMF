#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 14:33:16 2024

@author: selinayao & bilinzhuang
"""

# Calculate volume of BIMF
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors

# Setting a random starting point for the simulation
import random
import time
random.seed(round(time.time()*1000))  # randomize starting point of random numbers


# for automating the simulation through python inputs
import sys

if len(sys.argv)>1:
    jobname = sys.argv[1]
else:
    jobname = '1_bimf_0.50_water_RunTest'
    
variables = jobname.split('_')
n_bimf = int(variables[0])
volumeFrac_water = float(variables[2])
joblabel = variables[4]



# Approximate Box Size
L = 50
#n_bimf = 1
#volumeFrac_water = 0.0


# Estimate van der Waals volume of BIMF
mol = Chem.MolFromPDBFile("bimf.pdb")
v_bimf = AllChem.ComputeMolVolume(mol) # return molecular van der Waal's volume 


# Calculate molecular volume from density
def molecularVolume(density,molarMass):
    # return molecular volume in Angstrom**3
    return molarMass/density*10/6.022
    
# Calculate water and ethanol molecular volume
mol = Chem.MolFromPDBFile("water.pdb")
m_water = Descriptors.ExactMolWt(mol)
rho_water = 0.9970  # g/mL
v_water = molecularVolume(rho_water,m_water)


mol = Chem.MolFromPDBFile("ethanol.pdb")
m_ethanol = Descriptors.ExactMolWt(mol)
rho_ethanol = 0.7852  # g/mL
v_ethanol = molecularVolume(rho_ethanol,m_ethanol)



# Calculate number of molecules
V = L**3
V_solvent = V - n_bimf*v_bimf
V_water = V_solvent*volumeFrac_water
n_water = int(V_water/v_water)
n_ethanol = int((V_solvent-V_water)/v_ethanol)

V_actual = n_bimf*v_bimf + n_water*v_water + n_ethanol*v_ethanol
L_actual = V_actual**(1/3)
print(n_water)
print(n_ethanol)

# Pdb files
f_bimf = "bimf.pdb"
f_water = "water.pdb"
f_ethanol = "ethanol.pdb"


# Write packmol input file
PackmolInputFile = 'mixture_%ibimf_%.2fwater_%.2fEtOH_%s' %(n_bimf, volumeFrac_water, 1-volumeFrac_water,joblabel)
f = open(PackmolInputFile+'.inp','w')
f.write("""#
# A mixture of bimf, water, and ethanol
#
""")
f.write("""tolerance 2.5
filetype pdb
output """+PackmolInputFile+"_packmol.pdb\n\n")
f.write("""structure %s 
  number %i 
  center
  fixed  %.3f %.3f %.3f 0. 0. 0. 
end structure\n\n""" %(f_bimf,n_bimf,L_actual/2.0,L_actual/2.0,L_actual/2.0))

if n_water != 0:
    f.write("""structure %s 
      number %i 
      inside box 0. 0. 0. %.3f %.3f %.3f  
    end structure\n\n""" %(f_water,n_water,L_actual,L_actual,L_actual))

if n_ethanol != 0:
    f.write("""structure %s 
      number %i 
      inside box 0. 0. 0. %.3f %.3f %.3f  
    end structure\n\n""" %(f_ethanol,n_ethanol,L_actual,L_actual,L_actual))
    
f.close()

# Run packmol
packmolfile = PackmolInputFile+".inp"
#!./Packmol/packmol-20.14.4/packmol < $packmolfile


import os
os.system("./Packmol/packmol-20.14.4/packmol < "+packmolfile)


######
# Fixing bond information in the pdb file
######

def grep_conect_lines(pdb_file):
    conect_lines = []
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith("CONECT"):
                conect_lines.append(line)
    return conect_lines

def count_atoms(pdb_file):
    atom_count = 0
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                atom_count += 1
    return atom_count

def write_connect_lines(f_pdb,numMol,start_register):
    pdb_file = f_pdb
    num_atoms_in_mol = count_atoms(f_pdb)
    conect_lines = grep_conect_lines(pdb_file)
    
    for i in range(numMol):
        for line in conect_lines:
            atom_numbers = [start_register + i*num_atoms_in_mol + int(k) for k in line.strip().split()[1:]]
            connect_line = ('CONECT'+'%5i'*len(atom_numbers)+'\n') %tuple(atom_numbers)
            f_output.write(connect_line)
    

    
    

input_file = PackmolInputFile + "_packmol.pdb"
output_file = PackmolInputFile + ".pdb"



with open(input_file, 'r') as f_input, open(output_file, 'w') as f_output:
    for line in f_input:
        if line.startswith("TITLE"):
            f_output.write(line)
            f_output.write('CRYST1%9.3f%9.3f%9.3f  90.00  90.00  90.00 P 1           1 \n'%(L_actual,L_actual,L_actual))
        elif not line.startswith("CONECT"):
            f_output.write(line)

    start_register = 0
    
    write_connect_lines(f_bimf,n_bimf,start_register)
    start_register += n_bimf*count_atoms(f_bimf)
    
    write_connect_lines(f_water,n_water,start_register)
    start_register += n_water*count_atoms(f_water)
    
    write_connect_lines(f_ethanol,n_ethanol,start_register)
    start_register += n_ethanol*count_atoms(f_ethanol)
    




######
# Perform simulation
######

from openmm.app import *
from openmm import *
from openmm.unit import *

from sys import stdout



# Create an OpenFF Molecule object for benzene from SMILES
from openff.toolkit import Molecule
molecule = Molecule.from_file('BIMF.mol')

# Create the GAFF template generator
from openmmforcefields.generators import GAFFTemplateGenerator
gaff = GAFFTemplateGenerator(molecules=molecule)

# Create an OpenMM ForceField object with AMBER ff14SB and TIP3P with compatible ions
forcefield = ForceField('tip3p.xml','oplsaa-ethanol.xml')

# Register the GAFF template generator
forcefield.registerTemplateGenerator(gaff.generator)

# You can now parameterize an OpenMM Topology object that contains the specified molecule.
# forcefield will load the appropriate GAFF parameters when needed, and antechamber
# will be used to generate small molecule parameters on the fly.

#from openmm.app import PDBFile
#pdbfile = PDBFile('ligand_wH.pdb')
#system = forcefield.createSystem(pdbfile.topology)


# simulation parameters
equilibration_time = 10*nanosecond
simulation_time    = 50*nanosecond
t_step   = 0.001*picosecond
t_report = 0.01*nanosecond
DCD_output_file = "BIMF_trial.dcd"
final_pdb_file  = "BIMF_trial.pdb"


# import the BIMF molecule
pdbfile = PackmolInputFile+".pdb"
pdb = PDBFile(pdbfile)

#pdb = PDBFile("bimf.pdb")



# fixing OPLS force rules
def OPLS_LJ(system):
    forces = {system.getForce(index).__class__.__name__: system.getForce(
        index) for index in range(system.getNumForces())}
    nonbonded_force = forces['NonbondedForce']
    lorentz = CustomNonbondedForce(
        '4*epsilon*((sigma/r)^12-(sigma/r)^6); sigma=sqrt(sigma1*sigma2); epsilon=sqrt(epsilon1*epsilon2)')
    lorentz.setNonbondedMethod(nonbonded_force.getNonbondedMethod())
    lorentz.addPerParticleParameter('sigma')
    lorentz.addPerParticleParameter('epsilon')
    lorentz.setCutoffDistance(nonbonded_force.getCutoffDistance())
    system.addForce(lorentz)
    LJset = {}
    for index in range(nonbonded_force.getNumParticles()):
        charge, sigma, epsilon = nonbonded_force.getParticleParameters(index)
        LJset[index] = (sigma, epsilon)
        lorentz.addParticle([sigma, epsilon])
        nonbonded_force.setParticleParameters(
            index, charge, sigma, epsilon * 0)
    for i in range(nonbonded_force.getNumExceptions()):
        (p1, p2, q, sig, eps) = nonbonded_force.getExceptionParameters(i)
        # ALL THE 12,13 and 14 interactions are EXCLUDED FROM CUSTOM NONBONDED
        # FORCE
        lorentz.addExclusion(p1, p2)
        if eps._value != 0.0:
            #print p1,p2,sig,eps
            sig14 = sqrt(LJset[p1][0] * LJset[p2][0])
            eps14 = sqrt(LJset[p1][1] * LJset[p2][1])
            nonbonded_force.setExceptionParameters(i, p1, p2, q, sig14, eps)
    return system





print ("Building Model...")
modeller = Modeller(pdb.topology, pdb.positions) 

print ("Creating System...")
system = forcefield.createSystem(modeller.topology, nonbondedMethod=PME, nonbondedCutoff=1.0*nanometer, constraints=HBonds)
#system = OPLS_LJ(system)
integrator = LangevinIntegrator(300*kelvin, 1/picosecond, t_step)  # step size: 0.001*picoseconds
 
print ("Set up the simulation...")
simulation = Simulation(modeller.topology, system, integrator)
simulation.context.setPositions(modeller.positions)

# Minimizing the energy 
print("Minimizing Energy...")
simulation.minimizeEnergy(maxIterations=1000) 

# Equilibrate the System
print("Equilibrating...")
# simulation.step(10000000) # equilibrate for 10 ns (use on the cluster)  
simulation.step(int(equilibration_time/t_step))  #

# Adding Reporters 
print("Collecting Data...")
Nsteps_report = int(t_report/t_step)
simulation.reporters.append(DCDReporter(DCD_output_file, Nsteps_report)) #Every 100ps record one traj data
simulation.reporters.append(StateDataReporter(stdout, Nsteps_report, step=True, potentialEnergy=True, temperature=True, separator=' ')) 
simulation.reporters.append(PDBReporter(final_pdb_file, Nsteps_report)) 

# Running simulation 
# simulation.step(100000000) 
simulation.step(int(simulation_time/t_step))
