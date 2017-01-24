#!/usr/bin/python2.7
# DHFR rRESPA Benchmark using OpenMM
#
# 12 fs outer timestep
#  6 fs non-bonded far
#           long-range electrostatic interaction
#  4 fs non-bonded near
#           van der Waals and short-range electrostatic interactions
#  2 bonded
#
# refs:
# http://docs.openmm.org/7.0.0/api-python/generated/simtk.openmm.mtsintegrator.MTSIntegrator.html


from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import time


## Platform
################

#platform = openmm.Platform_getPlatformByName("OpenCL")
platform = openmm.Platform_getPlatformByName("CUDA")
#platform = openmm.Platform_getPlatformByName("Reference")


platformProperties = {}
## Precision
################

# OpenCL
#platformProperties['OpenCLPrecision'] = 'mixed'
# CUDA 
platformProperties['CudaPrecision'] = 'mixed'

## Parallel GPUs
################

# Run on multiple cards; current setup on vertex

#OpenCL parallel
#platformProperties['OpenCLDeviceIndex'] = '0,1,2'
#platformProperties['OpenCLDeviceIndex'] = '1'
#platformProperties['OpenCLDeviceIndex'] = '0'

# CUDA parallel
#platformProperties['CudaDeviceIndex'] = '0,1,2'
#platformProperties['CudaDeviceIndex'] = '1'
platformProperties['CudaDeviceIndex'] = '0'

prmtop = AmberPrmtopFile('prmtop')
inpcrd = AmberInpcrdFile('inpcrd',  loadVelocities=True, loadBoxVectors=True)

system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=0.8*nanometer, constraints=HBonds)

# Set the COM Removal to something sensible
for i in range(system.getNumForces()):
   if (type(system.getForce(i)) == openmm.CMMotionRemover):
      system.getForce(i).setFrequency(1000)

   if (type(system.getForce(i)) == openmm.NonbondedForce):
      # NFFT1 =   64       NFFT2 =   64       NFFT3 =   64
      # Ewald Coefficient =  0.39467 (A^-1)
      system.getForce(i).setPMEParameters(3.9467, 64, 64, 64)

# rRESPA
# By default, all forces are in force group 0, hence no need to mark
# the bonded forces

# Non-bonded near: van der Waals and short-range electrostatic interactions
near = [f for f in system.getForces() if isinstance(f, NonbondedForce)][0]
near.setForceGroup(1)

# Non-bonded far: long-range electrostatic interactions
far = [f for f in system.getForces() if isinstance(f, NonbondedForce)][0]
far.setReciprocalSpaceForceGroup(2)

#                                             far    near   bonded
integrator = MTSIntegrator(12*femtoseconds, [(2,2), (1,3), (0,6)])



simulation = Simulation(prmtop.topology, system, integrator, platform, platformProperties)

print "OpenMM version: %s" % (simulation.context.getPlatform().getOpenMMVersion())
print "Platform: %s" % (simulation.context.getPlatform().getName())
for item in simulation.context.getPlatform().getPropertyNames():
  print "%s: %s" % (item, simulation.context.getPlatform().getPropertyValue(simulation.context,item))

print ""
print "Number of atoms %i"      % len(inpcrd.positions)
print "Number of velocities %i" % len(inpcrd.velocities)

simulation.context.setPositions(inpcrd.positions)
simulation.context.setVelocities(inpcrd.velocities)

simulation.reporters.append(PDBReporter('output.pdb', 1000))
simulation.reporters.append(StateDataReporter(stdout, 1000, step=True, time=True,  totalEnergy=True, kineticEnergy=True, potentialEnergy=True, temperature=True, progress=True, totalSteps=10000, speed=True))

start_time = time.time()
simulation.step(10000) # i.e. 120,000 fs == 120 ps == 0.12 ns

# If it takes 100 seconds to run 0.12 ns,
# then it will take 1/0.12  * 100s to run one ns.
totalDynamicsRunTimeInSeconds = time.time() - start_time

timeNeedToRunOneNsinSeconds = (1/0.12) * totalDynamicsRunTimeInSeconds

NsPerDay = 86400 / timeNeedToRunOneNsinSeconds


print str(totalDynamicsRunTimeInSeconds) + " seconds"
print str(timeNeedToRunOneNsinSeconds) + " is the time needed to run 1 ns"
print str(NsPerDay)  + " ns/day"
