%% SCRIPT3: track only right shank IMU

% information used:
% 1. time bounds: assumes we can accurately get foot contact/off events
%       from shank IMU
% 2. accelerometer and gyroscope signals from shank
% 3. pelvis displacement: assumes we can accurately get pelvis displacement
%       from shank IMU: we assume stride symmetry, thus joint angles same
%       at beginning and end, thus cartesian displacement of every segment
%       is same as pelvis displacement, which could be calculated from IMU
%       signals several ways (e.g., gait speed * stride time, shank
%       tracking algos from pedestrian tracking, etc.)

% initial guess is solution from predictive simulation

% track shank accel/gyro signals generated from track all solution for full
% stride, same simple bounds placed on this solution as for predictive sim,
% symmetry was enforced on every state variable and control since
% simulating a full stride

% tracking mocap-tracking solution

clear
close all
clc

import org.opensim.modeling.*

% logistics
[scriptDir,scriptName] = fileparts(mfilename('fullpath'));
dateStr = datestr(datetime,'yyyymmddTHHMMSS');
sessionName = [scriptName,'_',dateStr];
trackDir = dir('OUTPUT_s1*');
trackDir = trackDir.name;
predictDir = dir('OUTPUT_s2*');
predictDir = predictDir.name;

% track IMUs from these frames
trackTheseFrames = {'tibia_r'};

% cost term weights
accelTrackingWeight = 10;
gyroTrackingWeight = 10;
controlEffortWeight = 1;

% get initial/final time and pelvis_tx value
trackSolutionFile = dir(fullfile(trackDir,'*solution_stride.sto'));
trackSolutionFile = fullfile(trackDir,trackSolutionFile.name);
trackSolutionTrajectory = MocoTrajectory(trackSolutionFile);
trackSolutionStates = trackSolutionTrajectory.exportToStatesTable;
trackSolutionControls = trackSolutionTrajectory.exportToControlsTable;
initialTime = trackSolutionStates.getIndependentColumn.get(0);
finalTime = trackSolutionStates.getIndependentColumn.get(trackSolutionStates.getNumRows-1);
initialPelvisX = trackSolutionStates.getDependentColumn('/jointset/groundPelvis/pelvis_tx/value').get(0);
finalPelvisX = trackSolutionStates.getDependentColumn('/jointset/groundPelvis/pelvis_tx/value').get(trackSolutionStates.getNumRows-1);

% get predictive solution file (used as initial guess here)
predictSolutionFile = dir(fullfile(predictDir,'*solution_stride.sto'));
predictSolutionFile = fullfile(predictDir,predictSolutionFile.name);

%% MODEL

% load 2D_gait model from MocoExample
model = Model(fullfile('NecessaryFilesFrom_Moco_example2DWalking','2D_gait.osim'));

% add IMU frames
imuFramePaths = StdVectorString;
for k = 1:length(trackTheseFrames)
    body = model.getBodySet.get(trackTheseFrames{k});
    imuname = [trackTheseFrames{k} '_imu'];
    imuframe = PhysicalOffsetFrame(imuname, body, Transform);
    imuframe.set_translation(Vec3(0,0,0));
    imuframe.set_orientation(Vec3(0,0,0));
    body.addComponent(imuframe);
    imuFramePaths.add(['/bodyset/' trackTheseFrames{k} '/' imuname]);
end
model.finalizeConnections;

% add model IMUs
OpenSenseUtilities.addModelIMUs(model, imuFramePaths);
model.initSystem;

%% STUDY

% init study
study = MocoStudy;
study.setName(sessionName);

%% PROBLEM

% get problem
problem = study.updProblem;

% set model
problem.setModel(model);

%% SIMULATE IMU SIGNALS FROM REFERENCE TRAJECTORY

% reference data

% accel
accelSignalPath = StdVectorString;
accelSignalPath.add('.*accelerometer_signal');
accelSignals = opensimSimulation.analyzeVec3(model,trackSolutionStates,trackSolutionControls,accelSignalPath);

% average variance of acceleration magnitude for weight normalization
accel_var = 0;
for k = 0:accelSignals.getNumColumns-1
    accel_var = accel_var + var(vecnorm(accelSignals.getDependentColumnAtIndex(k).getAsMat')) / accelSignals.getNumColumns;
end

% gyro
gyroSignalPath = StdVectorString;
gyroSignalPath.add('.*gyroscope_signal');
gyroSignals = opensimSimulation.analyzeVec3(model,trackSolutionStates,trackSolutionControls,gyroSignalPath);

% average var of gyroscope magnitude for weight normalization
gyro_var = 0;
for k = 0:gyroSignals.getNumColumns-1
    gyro_var = gyro_var + std(vecnorm(gyroSignals.getDependentColumnAtIndex(k).getAsMat')) / gyroSignals.getNumColumns;
end

% update column labels to match frame paths
accelSignals.setColumnLabels(imuFramePaths);
gyroSignals.setColumnLabels(imuFramePaths);

%% COST TERM: TRACK ACCELEROMETER SIGNALS

accelTracking = MocoAccelerationTrackingGoal('accelerometer_tracking',accelTrackingWeight / accel_var);
accelTracking.setFramePaths(imuFramePaths);
accelTracking.setAccelerationReference(accelSignals);
accelTracking.setGravityOffset(true);
accelTracking.setExpressAccelerationsInTrackingFrames(true);
problem.addGoal(accelTracking);

%% COST TERM: TRACK GYROSCOPE SIGNALS

gyroTracking = MocoAngularVelocityTrackingGoal('gyroscope_tracking',gyroTrackingWeight / gyro_var);
gyroTracking.setFramePaths(imuFramePaths);
gyroTracking.setAngularVelocityReference(gyroSignals);
problem.addGoal(gyroTracking);

%%  COST TERM: CONTROL EFFORT

problem.addGoal(MocoControlGoal('control_effort',controlEffortWeight));

%% CONSTRAINTS: BILATERAL SYMMETRY

% init symmetry goal enabling simulation of only single step from which a
% full stride can be reconstructed
symmetryGoal = MocoPeriodicityGoal('symmetryGoal');

% enforce symmetry on every state variable except pelvis_tx/value
% if contains activation then leading part of path is just /(musclename),
% use this to enforce symmetry on muscle controls
model.initSystem;
for k = 1:model.getNumStateVariables
    
    % get this state name + contralateral (used only if bilateral)
    thisStateName = char(model.getStateVariableNames().getitem(k-1));
    if ~endsWith(thisStateName,'pelvis_tx/value')
        symmetryGoal.addStatePair(MocoPeriodicityGoalPair(thisStateName));
        if endsWith(thisStateName,'/activation')
            symmetryGoal.addControlPair(MocoPeriodicityGoalPair(replace(thisStateName,'/activation','')));
        end
    end
end

% enforce symmetry in lumbar actuation
symmetryGoal.addControlPair(MocoPeriodicityGoalPair('/lumbarAct'));
    
% add goal
problem.addGoal(symmetryGoal);

%% CONSTRAINTS: SIMPLE BOUNDS

% set time bounds
problem.setTimeBounds(initialTime, finalTime);

% set bound
problem.setStateInfo('/jointset/groundPelvis/pelvis_tx/value', MocoBounds(initialPelvisX-0.1,initialPelvisX+2.0), MocoInitialBounds(initialPelvisX));

% simple bounds from example2DWalking
problem.setStateInfo('/jointset/groundPelvis/pelvis_tilt/value', [-20*pi/180, -10*pi/180]);
problem.setStateInfo('/jointset/groundPelvis/pelvis_ty/value', [0.75, 1.25]);
problem.setStateInfo('/jointset/hip_l/hip_flexion_l/value', [-10*pi/180, 60*pi/180]);
problem.setStateInfo('/jointset/hip_r/hip_flexion_r/value', [-10*pi/180, 60*pi/180]);
problem.setStateInfo('/jointset/knee_l/knee_angle_l/value', [-50*pi/180, 0]);
problem.setStateInfo('/jointset/knee_r/knee_angle_r/value', [-50*pi/180, 0]);
problem.setStateInfo('/jointset/ankle_l/ankle_angle_l/value', [-15*pi/180, 25*pi/180]);
problem.setStateInfo('/jointset/ankle_r/ankle_angle_r/value', [-15*pi/180, 25*pi/180]);
problem.setStateInfo('/jointset/lumbar/lumbar/value', [0, 20*pi/180]);
    
%% SOLVER

% solver settings
solver = MocoCasADiSolver.safeDownCast(study.updSolver);
solver.resetProblem(problem);
solver.set_optim_solver('ipopt');
solver.set_verbosity(2); % 0 for none, 1 for moco only, 2 for casadi output
solver.set_optim_ipopt_print_level(5);
solver.set_num_mesh_intervals(100);
solver.set_transcription_scheme('hermite-simpson'); % trapezoidal or hermite-simpson
solver.set_optim_convergence_tolerance(1e-3);
solver.set_optim_constraint_tolerance(1e-3);
solver.set_optim_max_iterations(10000);
solver.set_optim_hessian_approximation('limited-memory');
solver.set_optim_finite_difference_scheme('forward'); % forward, central, or backward
solver.set_parallel(1);

% initialize at with random
solver.setGuessFile(predictSolutionFile);

%% SOLVE

solution = study.solve;

% write solution
resultsDir = fullfile(scriptDir,['OUTPUT_',sessionName]);
mkdir(resultsDir)
solution.write(fullfile(resultsDir,[sessionName,'_solution_stride.sto']));

% write study
study.print(fullfile(resultsDir,[sessionName,'_study.xml']));

% save script
copyfile(fullfile(scriptDir,[scriptName,'.m']),fullfile(resultsDir,[sessionName,'_script.m']))

% write grf
forceNamesRightFoot = StdVectorString();
forceNamesRightFoot.add('contactHeel_r');
forceNamesRightFoot.add('contactFront_r');
forceNamesLeftFoot = StdVectorString();
forceNamesLeftFoot.add('contactHeel_l');
forceNamesLeftFoot.add('contactFront_l');
externalForcesTableFlatStep = opensimMoco.createExternalLoadsTableForGait(model,solution,forceNamesRightFoot,forceNamesLeftFoot);
STOFileAdapter.write(externalForcesTableFlatStep,fullfile(resultsDir,[sessionName,'_grf_stride.sto']));

% write muscle outputs
outputs = StdVectorString;
outputs.add('.*excitation');
outputs.add('.*activation');
outputs.add('.*pennation_angle');
outputs.add('.*tendon_length');
outputs.add('.*normalized_fiber_length');
outputs.add('.*fiber_length_along_tendon');
outputs.add('.*tendon_strain');
outputs.add('.*passive_force_multiplier');
outputs.add('.*active_force_length_multiplier');
outputs.add('.*normalized_fiber_velocity');
outputs.add('.*fiber_velocity_along_tendon');
outputs.add('.*force_velocity_multiplier');
outputs.add('.*fiber_force');
outputs.add('.*fiber_force_along_tendon');
outputs.add('.*active_fiber_force_along_tendon');
outputs.add('.*passive_fiber_force_along_tendon');
outputs.add('.*tendon_force');
outputs.add('.*tendon_power');
outputs.add('.*muscle_power');
STOFileAdapter.write(study.analyze(solution, outputs),fullfile(resultsDir,[sessionName,'_outputs_stride.sto']));

