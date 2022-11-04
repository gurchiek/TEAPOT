%% SCRIPT2: predictive simulation

% objective is minimize excitation squared

% only constraints are that initial/final time and pelvis_tx/value match
% track all solution (effectively ensuring equal velocity), bilateral
% symmetry, and the simple bounds placed on coordinate values (except
% pelvis_tx) from the example2DWalking predictive simulation script that
% comes with the OpenSim download: https://github.com/opensim-org/opensim-core/blob/36f88ea00bf8fd76d18c6fb4f1e87f4e1801da60/Bindings/Java/Matlab/examples/Moco/example2DWalking/example2DWalking.m#L310

% this solution initialized with a solution to the same problem computed
% with a 2021 iMac Apple M1 initialized with setGuess('bounds'). I have not
% been able to solve this problem (init with setGuess('bounds')) on my 2016
% MacBook pro in <10000 iterations.

clear
close all
clc

import org.opensim.modeling.*

% logistics
trackDir = dir('OUTPUT_s1*');
trackDir = trackDir.name;
[scriptDir,scriptName] = fileparts(mfilename('fullpath'));
dateStr = datestr(datetime,'yyyymmddTHHMMSS');
sessionName = [scriptName,'_',dateStr];

% cost term weights
controlEffortWeight = 1;

% get initial/final time and pelvis_tx value
trackSolutionFile = dir(fullfile(trackDir,'*solution_step.sto'));
trackSolutionFile = fullfile(trackDir,trackSolutionFile.name);
trackSolution = TimeSeriesTable(trackSolutionFile);
initialTime = trackSolution.getIndependentColumn.get(0);
finalTime = trackSolution.getIndependentColumn.get(trackSolution.getNumRows-1);
initialPelvisX = trackSolution.getDependentColumn('/jointset/groundPelvis/pelvis_tx/value').get(0);
finalPelvisX = trackSolution.getDependentColumn('/jointset/groundPelvis/pelvis_tx/value').get(trackSolution.getNumRows-1);

%% MODEL

% load 2D_gait model from MocoExample
model = Model(fullfile('NecessaryFilesFrom_Moco_example2DWalking','2D_gait.osim'));

%% STUDY

% init study
study = MocoStudy;
study.setName(sessionName);

%% PROBLEM

% get problem
problem = study.updProblem;

% set model
problem.setModel(model);

%%  COST: CONTROL EFFORT

problem.addGoal(MocoControlGoal('control_effort',controlEffortWeight));

%% CONSTRAINTS: BILATERAL SYMMETRY

% init symmetry goal enabling simulation of only single step from which a
% full stride can be reconstructed
symmetryGoal = MocoPeriodicityGoal('symmetryGoal');

% enforce symmetry s.t. initial value right = final value left (and vv) for
% all states (since model only 2D) except pelvis_tx/value
model.initSystem;
for k = 1:model.getNumStateVariables
    
    % get this state name + contralateral (used only if bilateral)
    thisStateName = char(model.getStateVariableNames().getitem(k-1));
    contralateralStateName = replace(thisStateName,{'_r/','_l/'},{'_l/','_r/'});
    isBilateral = ~strcmp(thisStateName,contralateralStateName);
    if isBilateral
        symmetryGoal.addStatePair(MocoPeriodicityGoalPair(thisStateName, contralateralStateName));
    else
        if ~endsWith(thisStateName,'pelvis_tx/value')
            symmetryGoal.addStatePair(MocoPeriodicityGoalPair(thisStateName));
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
problem.setStateInfo('/jointset/groundPelvis/pelvis_tx/value', MocoBounds(), MocoInitialBounds(initialPelvisX), MocoFinalBounds(finalPelvisX));

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

% initialize at with bounds
solver.setGuess('bounds');

%% SOLVE

step_unsync_traj = study.solve;

%% OUTPUTS

% write unsynhronized solution for step and stride
resultsDir = fullfile(scriptDir,['OUTPUT_',sessionName]);
mkdir(resultsDir)
step_unsync_traj.write(fullfile(resultsDir,[sessionName,'_solution_step_unsynchronized.sto']));
stride_unsync_traj = opensimMoco.createPeriodicTrajectory(step_unsync_traj);
stride_unsync_traj.write(fullfile(resultsDir,[sessionName,'_solution_stride_unsynchronized.sto']));

% get unsynchronized grf
forceNamesRightFoot = StdVectorString();
forceNamesRightFoot.add('contactHeel_r');
forceNamesRightFoot.add('contactFront_r');
forceNamesLeftFoot = StdVectorString();
forceNamesLeftFoot.add('contactHeel_l');
forceNamesLeftFoot.add('contactFront_l');
grf_step_unsync = opensimMoco.createExternalLoadsTableForGait(model,step_unsync_traj,forceNamesRightFoot,forceNamesLeftFoot);
grf_stride_unsync = opensimMoco.createExternalLoadsTableForGait(model,stride_unsync_traj,forceNamesRightFoot,forceNamesLeftFoot);
STOFileAdapter.write(grf_step_unsync,fullfile(resultsDir,[sessionName,'_grf_step_unsynchronized.sto']));
STOFileAdapter.write(grf_stride_unsync,fullfile(resultsDir,[sessionName,'_grf_stride_unsynchronized.sto']));

% synchronize stride + write
stride_unsync_table = TimeSeriesTable(fullfile(resultsDir,[sessionName,'_solution_stride_unsynchronized.sto']));
stride_sync_table = synchronizeWithBaseline(stride_unsync_table,grf_stride_unsync,trackDir);
STOFileAdapter.write(stride_sync_table,fullfile(resultsDir,[sessionName,'_solution_stride.sto']));
stride_sync_traj = MocoTrajectory(fullfile(resultsDir,[sessionName,'_solution_stride.sto']));

% write grf
grf_stride_sync = opensimMoco.createExternalLoadsTableForGait(model,stride_sync_traj,forceNamesRightFoot,forceNamesLeftFoot);
STOFileAdapter.write(grf_stride_sync,fullfile(resultsDir,[sessionName,'_grf_stride.sto']));

% write study
study.print(fullfile(resultsDir,[sessionName,'_study.xml']));

% save script
copyfile(fullfile(scriptDir,[scriptName,'.m']),fullfile(resultsDir,[sessionName,'_script.m']))

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
STOFileAdapter.write(study.analyze(stride_sync_traj, outputs),fullfile(resultsDir,[sessionName,'_outputs_stride.sto']));
