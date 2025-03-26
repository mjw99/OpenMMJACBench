#!/usr/bin/env python
# DHFR Benchmark using OpenMM

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

# Remember, this is being run NVE
integrator = mm.VerletIntegrator(2*unit.femtoseconds)

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
simulation.step(10000) # i.e. 20,000 fs == 20 ps == 0.02 ns

# If it takes 100 seconds to run 0.02 ns,
# then it will take 1/0.02 (==50)  * 100s to run one ns.
totalDynamicsRunTimeInSeconds = time.time() - start_time

timeNeedToRunOneNsinSeconds = (1/0.02) * totalDynamicsRunTimeInSeconds

NsPerDay = 86400 / timeNeedToRunOneNsinSeconds


print("%f seconds" % totalDynamicsRunTimeInSeconds)
print("%f is the time needed to run 1 ns" % timeNeedToRunOneNsinSeconds)
print("%f ns/day" % NsPerDay)

# Refs
# Python API docs
# https://simtk.org/api_docs/openmm/api6_1/python/



# http://wiki.simtk.org/openmm/BenchmarkOpenMMDHRF
#  1xC2070 = 30.9 ns/day
#
# http://ambermd.org/gpus/benchmarks.htm
#  1xM2090 = 43.74 ns/day
#  



###########
# Results #
###########

# OpenMM 6.2.0/CUDA 6.5, K40c, ecc on		66.58 ns/day	(25.95 run time)
# OpenMM 6.2.0/OpenCL 6.5, K40c, ecc on		48.22 ns/day	(35.83 run time)

# AMBER 14.0.1 / CUDA 6.5, K40c, ecc on         108.69 ns/day   (15.93 run time)
