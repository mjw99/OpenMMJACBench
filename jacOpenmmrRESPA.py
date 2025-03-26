#!/usr/bin/env python
# DHFR rRESPA Benchmark using OpenMM
#
# 12 fs outer timestep
#  6 fs non-bonded far
#           long-range electrostatic interaction
#  4 fs non-bonded near
#           van der Waals and short-range electrostatic interactions
#  2 fs bonded
#
# refs:
# http://docs.openmm.org/7.0.0/api-python/generated/simtk.openmm.mtsintegrator.MTSIntegrator.html
# http://dx.doi.org/10.1063/1.460004
# http://dx.doi.org/10.1063/1.463137


from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit
from sys import stdout
import time


## Platform
################

#platform = mm.Platform.getPlatformByName("OpenCL")
platform = mm.Platform.getPlatformByName("CUDA")
#platform = mm.Platform.getPlatformByName("Reference")


platformProperties = {}
## Precision
################

# OpenCL
#platformProperties['OpenCLPrecision'] = 'mixed'
# CUDA 
platformProperties['CudaPrecision'] = 'mixed'

## Parallel GPUs
################

# Run on multiple cards

# OpenCL parallel
#platformProperties['OpenCLDeviceIndex'] = '0,1,2'
#platformProperties['OpenCLDeviceIndex'] = '1'
#platformProperties['OpenCLDeviceIndex'] = '0'

# CUDA parallel
#platformProperties['CudaDeviceIndex'] = '0,1,2'
#platformProperties['CudaDeviceIndex'] = '1'
platformProperties['CudaDeviceIndex'] = '0'

prmtop = app.AmberPrmtopFile('prmtop')
inpcrd = app.AmberInpcrdFile('inpcrd',  loadVelocities=True, loadBoxVectors=True)

system = prmtop.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=0.8*unit.nanometer, constraints=app.HBonds)

# Set the COM Removal to something sensible
for i in range(system.getNumForces()):
   if (type(system.getForce(i)) == mm.CMMotionRemover):
      system.getForce(i).setFrequency(1000)

# Match the PME settings of the original benchmark
   if (type(system.getForce(i)) == mm.NonbondedForce):
      # NFFT1 =   64       NFFT2 =   64       NFFT3 =   64
      # Ewald Coefficient =  0.39467 (A^-1)
      system.getForce(i).setPMEParameters(3.9467, 64, 64, 64)

# rRESPA
# By default, all forces are in force group 0, hence no need to mark
# the bonded forces
 
# Non-bonded near: van der Waals and short-range electrostatic interactions
near = [f for f in system.getForces() if isinstance(f, mm.NonbondedForce)][0]
near.setForceGroup(1)

# Non-bonded far: long-range electrostatic interactions
far = [f for f in system.getForces() if isinstance(f, mm.NonbondedForce)][0]
far.setReciprocalSpaceForceGroup(2)

#                                             far    near   bonded
integrator = mm.MTSIntegrator(12*unit.femtoseconds, [(2,2), (1,3), (0,6)])


simulation = app.Simulation(prmtop.topology, system, integrator, platform, platformProperties)

print("OpenMM version: %s" % simulation.context.getPlatform().getOpenMMVersion())
print("Platform: %s" % simulation.context.getPlatform().getName())
for item in simulation.context.getPlatform().getPropertyNames():
  print("%s: %s" % (item, simulation.context.getPlatform().getPropertyValue(simulation.context,item)))

print()
print("Number of atoms %i"      % len(inpcrd.positions))
print("Number of velocities %i" % len(inpcrd.velocities))

simulation.context.setPositions(inpcrd.positions)
simulation.context.setVelocities(inpcrd.velocities)

simulation.reporters.append(app.PDBReporter('output.pdb', 1000))
simulation.reporters.append(app.StateDataReporter(stdout, 1000, step=True, time=True,  totalEnergy=True, kineticEnergy=True, potentialEnergy=True, temperature=True, progress=True, totalSteps=10000, speed=True))

start_time = time.time()
simulation.step(10000) # i.e. 120,000 fs == 120 ps == 0.12 ns

# If it takes 100 seconds to run 0.12 ns,
# then it will take 1/0.12  * 100s to run one ns.
totalDynamicsRunTimeInSeconds = time.time() - start_time

timeNeedToRunOneNsinSeconds = (1/0.12) * totalDynamicsRunTimeInSeconds

NsPerDay = 86400 / timeNeedToRunOneNsinSeconds


print("%f seconds" % totalDynamicsRunTimeInSeconds)
print("%f is the time needed to run 1 ns" % timeNeedToRunOneNsinSeconds)
print("%f ns/day" % NsPerDay)

