#!/usr/bin/python2.7
# DHFR HMR (5fs step) Benchmark using OpenMM


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
platformProperties['CudaDeviceIndex'] = '1'

prmtop = AmberPrmtopFile('prmtop')
inpcrd = AmberInpcrdFile('inpcrd',  loadVelocities=True, loadBoxVectors=True)

hydrogenMass = 4*amu

system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=0.8*nanometer, constraints=AllBonds, hydrogenMass=hydrogenMass)

# Set the COM Removal to something sensible
for i in range(system.getNumForces()):
   if (type(system.getForce(i)) == openmm.CMMotionRemover):
      system.getForce(i).setFrequency(1000)

   if (type(system.getForce(i)) == openmm.NonbondedForce):
      # NFFT1 =   64       NFFT2 =   64       NFFT3 =   64
      # Ewald Coefficient =  0.39467 (A^-1)
      system.getForce(i).setPMEParameters(3.9467, 64, 64, 64)

# Remember, this is being run NVE
integrator = VerletIntegrator(5*femtoseconds)

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
simulation.step(10000) # i.e. 50,000 fs == 50 ps == 0.05 ns

# If it takes 100 seconds to run 0.05 ns,
# then it will take 1/0.05 (==20)  * 100s to run one ns.
totalDynamicsRunTimeInSeconds = time.time() - start_time

timeNeedToRunOneNsinSeconds = (1/0.05) * totalDynamicsRunTimeInSeconds

NsPerDay = 86400 / timeNeedToRunOneNsinSeconds


print str(totalDynamicsRunTimeInSeconds) + " seconds"
print str(timeNeedToRunOneNsinSeconds) + " is the time needed to run 1 ns"
print str(NsPerDay)  + " ns/day"