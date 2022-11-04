%% SCRIPT1: baseline tracking simulation, track all data

% uses data and model from the MocoExamples/example2DWalking

% tracks all coordinates except pelvis_ty, this allowed to be adjusted in
% order to better match grf. Enforcing bilateral symmetry as constraint.
% Enforcing initial conditions from reference data as constraint. Initial
% pelvis_tx set to 0 and final set to net displacement from the reference
% data.

clear
close all
clc

import org.opensim.modeling.*

% logistics
mocoExamplesDir = 'NecessaryFilesFrom_Moco_example2DWalking';
[scriptDir,scriptName] = fileparts(mfilename('fullpath'));
dateStr = datestr(datetime,'yyyymmddTHHMMSS');
sessionName = [scriptName,'_',dateStr];

% track only coordinate values but not pelvis_ty
trackStatesWithTheseExpressions = {'value'}; % state name must contain these
doNotTrackStatesWithTheseExpressions = {'pelvis_ty'}; % AND state name must not contain these

% cost term weights
stateTrackingWeight = 10;
controlEffortWeight = 1;
grfTrackingWeight = 1;

% intial and final time for first step
initialTime = 0.0;
finalTime = 0.47008941;

%% MODEL

% load 2D_gait model from MocoExample
model = Model(fullfile(mocoExamplesDir,'2D_gait.osim'));

%% STUDY

% init study
study = MocoStudy;
study.setName(sessionName);

%% PROBLEM

% get problem
problem = study.updProblem;

% set model
problem.setModel(model);

%% COST: STATE TRACKING

% init state tracking cost
stateTrackingGoal = MocoStateTrackingGoal;
stateTrackingGoal.setName('state_tracking');
stateTrackingGoal.setAllowUnusedReferences(true);

% get reference
referenceCoordinatesProcessor = TableProcessor(fullfile(mocoExamplesDir,'referenceCoordinates.sto'));
referenceCoordinatesProcessor.append(TabOpLowPassFilter(6));
stateTrackingGoal.setReference(referenceCoordinatesProcessor);

% set weight for each state
model.initSystem;
stateNames = model.getStateVariableNames;
for k = 1:model.getNumStateVariables
    if contains(char(stateNames.get(k-1)),trackStatesWithTheseExpressions) && ~contains(char(stateNames.get(k-1)),doNotTrackStatesWithTheseExpressions)
        stateTrackingGoal.setWeightForState(stateNames.get(k-1),stateTrackingWeight);
    else
        stateTrackingGoal.setWeightForState(stateNames.get(k-1),0);
    end
end

% scale weights with range
stateTrackingGoal.setScaleWeightsWithRange(true);

% add to problem
problem.addGoal(stateTrackingGoal);

%%  COST: CONTROL EFFORT

problem.addGoal(MocoControlGoal('control_effort',controlEffortWeight));

%% COST: GRF TRACKING

% track the right and left vertical and fore-aft ground reaction forces
grfTrackingGoal = MocoContactTrackingGoal('grf_tracking', grfTrackingWeight);
grfTrackingGoal.setExternalLoadsFile(fullfile(mocoExamplesDir,'referenceGRF.xml'));
forceNamesRightFoot = StdVectorString();
forceNamesRightFoot.add('contactHeel_r');
forceNamesRightFoot.add('contactFront_r');
grfTrackingGoal.addContactGroup(forceNamesRightFoot, 'Right_GRF');
forceNamesLeftFoot = StdVectorString();
forceNamesLeftFoot.add('contactHeel_l');
forceNamesLeftFoot.add('contactFront_l');
grfTrackingGoal.addContactGroup(forceNamesLeftFoot, 'Left_GRF');
grfTrackingGoal.setProjection('plane');
grfTrackingGoal.setProjectionVector(Vec3(0, 0, 1));
problem.addGoal(grfTrackingGoal);

%% CONSTRAINTS: BILATERAL SYMMETRY

% init symmetry goal enabling simulation of only single step from which a
% full stride can be reconstructed
symmetryGoal = MocoPeriodicityGoal('symmetryGoal');

% enforce symmetry s.t. initial value right = final value left (and vv) for
% all states (since model only 2D) except pelvis_tx/value
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

% process referenceCoordinates to use as reference
referenceCoordinates = referenceCoordinatesProcessor.process(model);

% filtering adds data, trim
referenceCoordinates.trim(initialTime-1e-6,finalTime+1e-6);

% get pelvis displacement for pelvis_tx final bound
pelvisDisplacement = referenceCoordinates.getDependentColumn('/jointset/groundPelvis/pelvis_tx/value').get(referenceCoordinates.getNumRows-1) - ...
                     referenceCoordinates.getDependentColumn('/jointset/groundPelvis/pelvis_tx/value').get(0);

% set bound
problem.setStateInfo('/jointset/groundPelvis/pelvis_tx/value', MocoBounds(), MocoInitialBounds(0), MocoFinalBounds(pelvisDisplacement));

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

% initialize guess with bound midpoint + use coords/speeds from ref data
initialGuess = solver.createGuess('bounds');

% add speeds to coordinates table
n = referenceCoordinates.getNumRows;
t = zeros(1,n);
for k = 1:n; t(k) = referenceCoordinates.getIndependentColumn.get(k-1); end
referenceColumnLabels = referenceCoordinates.getColumnLabels;
for k = 0:referenceCoordinates.getNumColumns-1
    thisColumnName = char(referenceColumnLabels.get(k));
    if endsWith(thisColumnName,'/value')
        thisColumn = referenceCoordinates.getDependentColumn(thisColumnName);
        q = zeros(1,n);
        for j = 1:n; q(j) = thisColumn.get(j-1); end
        qdot = fdiff5(q,t);
        qdotVec = Vector(n,0.0);
        for j = 1:n; qdotVec.set(j-1,qdot(j)); end
        speedColumnName = replace(thisColumnName,'/value','/speed');
        referenceCoordinates.appendColumn(speedColumnName,qdotVec);
    end
end

% insert referenceCoordinates for initial guess
initialGuess.insertStatesTrajectory(referenceCoordinates,true);

% set guess
solver.setGuess(initialGuess);

%% SOLVE

solution = study.solve;

% get full stride
fullStride = opensimMoco.createPeriodicTrajectory(solution);

% write solution
resultsDir = fullfile(scriptDir,['OUTPUT_',sessionName]);
mkdir(resultsDir)
fullStride.write(fullfile(resultsDir,[sessionName,'_solution_stride.sto']));
solution.write(fullfile(resultsDir,[sessionName,'_solution_step.sto']));

% write study
study.print(fullfile(resultsDir,[sessionName,'_study.xml']));

% save script
copyfile(fullfile(scriptDir,[scriptName,'.m']),fullfile(resultsDir,[sessionName,'_script.m']))

% write grf
externalForcesTableFlatStride = opensimMoco.createExternalLoadsTableForGait(model,fullStride,forceNamesRightFoot,forceNamesLeftFoot);
STOFileAdapter.write(externalForcesTableFlatStride,fullfile(resultsDir,[sessionName,'_grf_stride.sto']));
externalForcesTableFlatStep = opensimMoco.createExternalLoadsTableForGait(model,solution,forceNamesRightFoot,forceNamesLeftFoot);
STOFileAdapter.write(externalForcesTableFlatStep,fullfile(resultsDir,[sessionName,'_grf_step.sto']));

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
STOFileAdapter.write(study.analyze(fullStride, outputs),fullfile(resultsDir,[sessionName,'_outputs_stride.sto']));
STOFileAdapter.write(study.analyze(solution, outputs),fullfile(resultsDir,[sessionName,'_outputs_step.sto']));
