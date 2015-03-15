#!/usr/bin/python2.7
# DHFR Benchmark using OpenMM


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

# Remember, this is being run NVE
integrator = VerletIntegrator(2*femtoseconds)

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
simulation.reporters.append(StateDataReporter(stdout, 1000, step=True, time=True,  totalEnergy=True, kineticEnergy=True, potentialEnergy=True, temperature=True))

start_time = time.time()
simulation.step(10000) # i.e. 20,000 fs == 20 ps == 0.02 ns

# If it takes 100 seconds to run 0.02 ns,
# then it will take 1/0.02 (==50)  * 100s to run one ns.
totalDynamicsRunTimeInSeconds = time.time() - start_time

timeNeedToRunOneNsinSeconds = (1/0.02) * totalDynamicsRunTimeInSeconds

NsPerDay = 86400 / timeNeedToRunOneNsinSeconds


print str(totalDynamicsRunTimeInSeconds) + " seconds"
print str(timeNeedToRunOneNsinSeconds) + " is the time needed to run 1 ns"
print str(NsPerDay)  + " nS/day"

# Refs
# Python API docs
# https://simtk.org/api_docs/openmm/api5_1/python/



# http://wiki.simtk.org/openmm/BenchmarkOpenMMDHRF
#  1xC2070 = 30.9 ns/day
#
# http://ambermd.org/gpus/benchmarks.htm
#  1xM2090 = 43.74 ns/day
#  



###########
# Results #
###########

# AMBER11, M2090,  ecc on 			38.12 ns/day	(45.28s run time)
# AMBER11, M2090,  ecc off 			42.48 ns/day	(40.67s run time)

# AMBER12 (Bugfix 15), M2090,  ecc on		41.92 ns/day	(41.40s run time)

# OpenMM 5.0/Reference				x ns/day    (x run time)
# OpenMM 5.0/OpenCL, M2090,  ecc on 		17.27 ns/day  (	100.04 run time)
# OpenMM 5.0/CUDA, M2090,  ecc on 		19.56 ns/day  	(88.30 run time)
# OpenMM 5.0/OpenCL, C2075,  ecc on           	12.18 ns/day  	(141.21 run time)

# OpenMM 5.1/OpenCL (CUDA 5.5), M2090, ecc on	17.47 ns/day	(98.89 run time)
# OpenMM 5.1/CUDA 5.5, M2090, ecc on		21.53 ns/day	(85.44 run time)
# OpenMM 5.1/OpenCL (CUDA 5.5), C2075, ecc on	14.76 ns/day	(117.02 run time)
# OpenMM 5.1/CUDA 5.5, C2075, ecc on		16.90 ns/day	(102.22 run time)


# OpenMM 6.0/OpenCL 5.5, M2090, ecc on          17.72 ns/day    (97.49 run time)
# OpenMM 6.0/CUDA 5.5, M2090, ecc on            20.30 ns/day    (85.10 run time)
# OpenMM 6.0/OpenCL 5.5, C2075, ecc on		14.64 ns/day	(118.02 run time)
# OpenMM 6.0/CUDA 5.5, C2075, ecc on		16.91 ns/day	(102.16 run time)
# OpenMM 6.0/CUDA 5.5, K40c, ecc on		51.37 ns/day	(33.63 run time)
# OpenMM 6.0/OpenCL (CUDA 5.5), K40c, ecc on	36.61 ns/day	(47.19 run time)

# OpenMM 6.1/CUDA 6.0, K40c, ecc on		49.73 ns/day	(34.75 run time)
# OpenMM 6.1/OpenCL 6.0, K40c, ecc on		40.17 ns/day	(43.02 run time)

# OpenMM 6.2.0/CUDA 6.5, K40c, ecc on		57.30 ns/day	(30.15 run time)
# OpenMM 6.2.0/OpenCL 6.5, K40c, ecc on		47.23 ns/day	(36.59 run time)
