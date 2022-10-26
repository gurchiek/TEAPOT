# TEAPOT

TEAPOT is an acronym for Track Everything Available, Predict Other Things
This describes a general approach for sparse-signal tracking simulations.
A simple example can be implemented from this code using OpenSim Moco.
You will need to have already [setup](https://simtk-confluence.stanford.edu:8443/display/OpenSim/Scripting+with+Matlab) use of the OpenSim API within MATLAB.

## Step 1
Run s1_baselineTrackingSimulation.
This generates a physically consistent 2D walking gait by tracking reference coordinate and ground reaction force data. A single step is simulated and constraints are imposed to force bilateral symmetry in every state and the lumbar actuator control signal. An OUTPUT_s1_* folder will be created in which simulation results are stored.

## Step 2
Run s2_predictiveSimulation.
This is a predictive simulation of 2D walking gait. No reference data are tracked. Constraints are imposed to force bilateral symmetry as in the step 1 tracking simulation and to match the gait speed of the reference data. The only objective function term represented control effort: the sum of the squared control signals of all system actuators (including muscles). Once the simulation has finished, the function synchronizeWithBaseline is called to time-synchronize the predicitve simulation trajectory with the reference tracking simulation from step 1 (using the GRF signal). An OUTPUT_s2_* folder will be created in which simulation results are stored (all files ending with _unsynchronized* are not synchronized with the reference tracking simulation; else they are).

## Step 3
Run s3_trackOnlyRightShankIMU.
This is a tracking simulation where only (simulated) signals from a right shank-worn gyroscope (angular velocity) and accelerometer (specific force) are tracked. The IMU signals are generated using the reference tracking simulation. The simulated signals are noise- and bias-free. A full stride was simulated instead of a single step as in the simulations of steps 1 and 2 because no bilateral symmetry constraints were enforced in this simulation. A symmetry constraint was placed on each state with itself which models periodicity of the stride (the initial state is exactly the same as the last). No constraint was used to match the gait speed of the reference simulation of step 1. 

## Step 4
Run s4_compareSimulations.
The compares the difference between simulations 2 and 3 with simulation 1. The results demonstrate the improved performance facilitated by tracking the signals from even a single IMU under ideal conditions (noise-/bias-free signals). Characterization of mechanical variables are generally better for the right leg (the side on which the sensor signals were tracked).
