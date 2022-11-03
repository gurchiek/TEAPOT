function [solution_sync] = synchronizeWithBaseline(solution,grf,baselinedir)

% adjust time to better match baseline simulation

import org.opensim.modeling.*

% clone unsynchronized solution and grf
solution_sync = solution.clone;

% get GRF at first instant of baseline solution, use this to sync
baselinename = replace(baselinedir,'OUTPUT_','');
baselinegrf = TimeSeriesTable(fullfile(baselinedir,[baselinename '_grf_step.sto']));
grfthresh = baselinegrf.getDependentColumn('ground_force_r_vy').get(0);

% get instant closest to grfthresh
[~,icontact] = min(abs(grf.getDependentColumn('ground_force_r_vy').getAsMat - grfthresh));
if icontact == grf.getNumRows; return; end

% adjust
iunsync = icontact;
sync_pelvis_tx = 0;
for isync = 0:solution.getNumRows-2
    solution_sync.setRowAtIndex(isync,solution.getRowAtIndex(iunsync));
    solution_sync.getDependentColumn('/jointset/groundPelvis/pelvis_tx/value').set(isync,sync_pelvis_tx);
    change_pelvis_tx = solution.getDependentColumn('/jointset/groundPelvis/pelvis_tx/value').get(iunsync+1) - solution.getDependentColumn('/jointset/groundPelvis/pelvis_tx/value').get(iunsync);
    sync_pelvis_tx = sync_pelvis_tx + change_pelvis_tx;
    iunsync = iunsync+1;
    if iunsync > solution.getNumRows-2; iunsync = 0; end
end

% copy first row to last
solution_sync.setRowAtIndex(solution_sync.getNumRows-1,solution_sync.getRowAtIndex(0));

% adjust pelvis at end 
final_dt = solution_sync.getIndependentColumn.get(solution_sync.getNumRows-1) - solution_sync.getIndependentColumn.get(solution_sync.getNumRows-2);
final_change_pelvis_tx = final_dt * solution_sync.getDependentColumn('/jointset/groundPelvis/pelvis_tx/speed').get(solution_sync.getNumRows-2);
final_pelvis_tx = solution_sync.getDependentColumn('/jointset/groundPelvis/pelvis_tx/value').get(solution_sync.getNumRows-2) + final_change_pelvis_tx;
solution_sync.getDependentColumn('/jointset/groundPelvis/pelvis_tx/value').set(solution_sync.getNumRows-1,final_pelvis_tx);

end

