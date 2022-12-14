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
icontact = crossing0(grf.getDependentColumn('ground_force_r_vy').getAsMat - grfthresh,{'n2p','z2p','n2z'});
icontact = icontact(1);
if icontact == grf.getNumRows; return; end
icontact = icontact - 1;

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

function [ i, type ] = crossing0( x, type0 )
%Reed Gurchiek,
%   crossing0 finds the zero crossings in the 1D array x.
%
%---------------------------INPUTS-----------------------------------------
%
%   x:
%       1D array
%       
%   type0 (optional):
%       cell array, type of zero crossing.  
%           1) 'p2z': positive to zerO
%           2) 'p2n': positive to negative
%           3) 'z2n': zero to negative
%           4) 'n2z': negative to zero
%           5) 'n2p': negative to positive
%           6) 'z2p': zero to positive
%           7) 'all': 1) through 6)
%               -combos allowed (e.g. type = {'p2n' 'n2p'});
%
%--------------------------OUTPUTS-----------------------------------------
%
%   i:
%       1xp array of zero crossing indices
%
%   type: 
%       1xp cell array of type of crossing associated with the indices in i
%       (see INPUT type for description of names
%
%--------------------------------------------------------------------------
%% crossing0

%type
if nargin == 2
    if ~iscell(type0)
        if strcmpi(type0,'all')
            type0 = {'p2z' 'p2n' 'z2n' 'n2z' 'n2p' 'z2p'};
        else
            type0{1} = type0;
        end
    elseif any(strcmpi('all',type0))
        type0 = {'p2z' 'p2n' 'z2n' 'n2z' 'n2p' 'z2p'};
    end
else
    type0 = {'p2z' 'p2n' 'z2n' 'n2z' 'n2p' 'z2p'};
end

%for each sample
ct = 0;
sgn0 = sign(x(1));
i = [];
type = {''};
for k = 2:length(x)
    
    %current sign
    sgn = sign(x(k));
    
    %if crossed
    if sgn ~= sgn0
        
        %positive to zero
        if sgn0 == 1 && sgn == 0
            if any(strcmpi('p2z',type0))
                ct = ct + 1;
                i(ct) = k;
                type{ct} = 'p2z';
            end
        %positive to negative
        elseif sgn0 == 1 && sgn == -1
            if any(strcmpi('p2n',type0))
                ct = ct + 1;
                i(ct) = k;
                type{ct} = 'p2n';
            end
        %zero to negative
        elseif sgn0 == 0 && sgn == -1
            if any(strcmpi('z2n',type0))
                ct = ct + 1;
                i(ct) = k;
                type{ct} = 'z2n';
            end
        %negative to zero
        elseif sgn0 == -1 && sgn == 0
            if any(strcmpi('n2z',type0))
                ct = ct + 1;
                i(ct) = k;
                type{ct} = 'n2z';
            end
        %negative to positive
        elseif sgn0 == -1 && sgn == 1
            if any(strcmpi('n2p',type0))
                ct = ct + 1;
                i(ct) = k;
                type{ct} = 'n2p';
            end
        %zero to positive
        elseif sgn0 == 0 && sgn == 1
            if any(strcmpi('z2p',type0))
                ct = ct + 1;
                i(ct) = k;
                type{ct} = 'z2p';
            end
        end
        
    end
    
    %update previous sign
    sgn0 = sgn;
    
end
    

end