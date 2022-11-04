%% compare muscle outputs

clear
close all
clc

import org.opensim.modeling.*

% analyze stride
suffix = '_stride';

% load 2d gait model
model = Model(fullfile('NecessaryFilesFrom_Moco_example2DWalking','2D_gait.osim'));

% output folders
f1 = dir('OUTPUT_s1*'); f1 = f1(end).name;
f2 = dir('OUTPUT_s2*'); f2 = f2(end).name;
f3 = dir('OUTPUT_s3*'); f3 = f3(end).name;
folder = {f1, f2, f3};
labels = {'baseline','predict','IMU-tracking'};
n = 3;

% color
color = [0 0 0; 140 21 21; 45 136 246]/255;
alpha = [0.9; 1.0; 1.0];
color_spec = [color alpha];

% line style
linestyle = {'-','-','-','-'};

% line width
linewidth = [1.5,1.5,1.5,1.5];

% get left/right vgrf
lfc = zeros(n,1);
lfo = zeros(n,1);
rfo = zeros(n,1);
x = 0:0.25:100;
vgrf.r = zeros(n,length(x));
vgrf.l = zeros(n,length(x));
for k = 1:n
    
    % get raw right/left vgrf + time
    grf = TimeSeriesTable(fullfile(folder{k},[replace(folder{k},'OUTPUT_','') '_grf' suffix '.sto']));
    t = zeros(1,grf.getNumRows);
    for j = 1:grf.getNumRows; t(j) = grf.getIndependentColumn.get(j-1); end
    r = grf.getDependentColumn('ground_force_r_vy').getAsMat;
    l = grf.getDependentColumn('ground_force_l_vy').getAsMat;
    
    % get left foot contact
    [~,ind_min] = min(l);
    ind_contact = find(l > 30);
    lfc(k) = ind_contact(1);
    
    % filter grf
    dt = mean(diff(t));
    
    % resample
    vgrf.r(k,:) = interp1(linspace(0,100,length(r)),r,x,'pchip');
    vgrf.l(k,:) = interp1(linspace(0,100,length(l)),l,x,'pchip');
    
end

% get solutions
for k = 1:n
    
    % import sto
    solution(k) = TimeSeriesTable(fullfile(folder{k},[replace(folder{k},'OUTPUT_','') '_solution' suffix '.sto']));
    
end

% get muscle outputs
muscles = {'bifemsh','gastroc','glut_max','hamstrings','iliopsoas','rect_fem','soleus','tib_ant','vasti'};
outputs = {'tendon_force','muscle_power','activation'};
side = {'r','l'};
for k = 1:n
    
    % import sto
    out = TimeSeriesTable(fullfile(folder{k},[replace(folder{k},'OUTPUT_','') '_outputs' suffix '.sto']));
    t = zeros(1,grf.getNumRows);
    for j = 1:out.getNumRows; t(j) = out.getIndependentColumn.get(j-1); end
    
    % for each muscle
    for m = 1:length(muscles)
        
        % for each output
        for o = 1:length(outputs)
            
            % for each side
            for s = 1:2
                
                % get signal
                v = out.getDependentColumn(['/' muscles{m} '_' side{s} '|' outputs{o}]).getAsMat;
                
                % store interpolation
                muscle.([muscles{m} '_' side{s}]).(outputs{o})(k,:) = interp1(linspace(0,100,length(v)),v,x,'pchip');
                
                % calculate work if power requested
                if strcmp(outputs{o},'muscle_power')
                    
                    con_power = v;
                    ecc_power = v;
                    
                    con_power(con_power < 0) = 0;
                    ecc_power(ecc_power > 0) = 0;
                    
                    muscle.([muscles{m} '_' side{s}]).net_eccentric_work_stance(k) = -trapz(t(1:length(v)),ecc_power);
                    muscle.([muscles{m} '_' side{s}]).net_concentric_work_stance(k) = trapz(t(1:length(v)),con_power);
                    
                end
                
            end
            
        end
        
    end
    
end

% get moment arms
state = model.initSystem;
for k = 1:n
    
    % load solution file
    traj = MocoTrajectory(fullfile(folder{k},[replace(folder{k},'OUTPUT_','') '_solution' suffix '.sto']));
    statevars = traj.exportToStatesTable;
    columnLabels = statevars.getColumnLabels;
    
    numrows = statevars.getNumRows;
    if k == 1; numrows = numrows - 1; end
    
    % for each state
    for j = 1:numrows
        
        % set the state
        for i = 0:statevars.getNumColumns-1
            model.setStateVariableValue(state,columnLabels.get(i),statevars.getDependentColumn(columnLabels.get(i)).get(j-1));
        end
        model.realizePosition(state);
        
        % for each muscle
        for m = 1:length(muscles)
            
            % for each side
            for s = 1:2
                
                % get muscle and ankle/knee coords
                musc = DeGrooteFregly2016Muscle.safeDownCast(model.getComponent([muscles{m} '_' side{s}]));
                knee = model.getCoordinateSet.get(['knee_angle_' side{s}]);
                ankle = model.getCoordinateSet.get(['ankle_angle_' side{s}]);
                
                % get moment arm
                moment_arm.([muscles{m} '_' side{s}]).(['knee_angle_' side{s}])(j) = musc.computeMomentArm(state,knee);
                moment_arm.([muscles{m} '_' side{s}]).(['ankle_angle_' side{s}])(j) = musc.computeMomentArm(state,ankle);
                
            end
            
        end
        
    end
    
    % for each muscle
    for m = 1:length(muscles)

        % for each side
        for s = 1:2
            
            % get knee moment arm
            v = moment_arm.([muscles{m} '_' side{s}]).(['knee_angle_' side{s}]);

            % store interpolation
            muscle.([muscles{m} '_' side{s}]).(['knee_angle_' side{s} '_momentArm'])(k,:) = interp1(linspace(0,100,length(v)),v,x,'pchip');
            
            % get ankle moment arm
            v = moment_arm.([muscles{m} '_' side{s}]).(['ankle_angle_' side{s}]);

            % store interpolation
            muscle.([muscles{m} '_' side{s}]).(['ankle_angle_' side{s} '_momentArm'])(k,:) = interp1(linspace(0,100,length(v)),v,x,'pchip');
            
        end
        
    end
    
end

% get knee/ankle joint moment
moment.rknee = zeros(n,length(x));
moment.lknee = zeros(n,length(x));
moment.rankle = zeros(n,length(x));
moment.lankle = zeros(n,length(x));
for k = 1:n
    
    % for each muscle
    for m = 1:length(muscles)
        
        % for each side
        for s = 1:2
            
            % knee moment
            f = muscle.([muscles{m} '_' side{s}]).tendon_force(k,:);
            r = muscle.([muscles{m} '_' side{s}]).(['knee_angle_' side{s} '_momentArm'])(k,:);
            moment.([side{s} 'knee'])(k,:) = moment.([side{s} 'knee'])(k,:) + r .* f;
            
            % ankle moment
            f = muscle.([muscles{m} '_' side{s}]).tendon_force(k,:);
            r = muscle.([muscles{m} '_' side{s}]).(['ankle_angle_' side{s} '_momentArm'])(k,:);
            moment.([side{s} 'ankle'])(k,:) = moment.([side{s} 'ankle'])(k,:) + r .* f;
            
        end
        
    end
    
end
   
%% TIME SERIES PLOTS: GRF
        
fig = figure;
fig.Position = [1137 546 426 180];

% for each side
for s = 1:2

    sp = subplot(1,2,s);

    % for each simulation
    for k = 1:n

        plot(x,vgrf.(side{s})(k,:),'Color',color_spec(k,:),'LineStyle',linestyle{k},'LineWidth',linewidth(k));
        hold on
        
        if k > 1; fprintf('%s %s vGRF RMSE: %f N\n',labels{k},side{s},rms(vgrf.(side{s})(k,:)-vgrf.(side{s})(1,:))); end

    end

    sp.Box = 'off';
    if s == 1; ylabel('Force (N)'); end
    xlabel('% Stride')
    title([side{s} ' vGRF'])

end
   
%% TIME SERIES PLOTS: MUSCLE OUTPUTS

% muscles to plot
muscplot = {'vasti','rect_fem','gastroc','soleus'};

% variables to plot
varplot = {'tendon_force','muscle_power'};

% titles
vtitle = {'Ft','Pm'};
mtitle = {'Vasti','RF','Gastroc','Soleus'};

% axis labels
ylab = {'Force (N)','Power (W)'};

% units
units = {'N','W'};
        
% for each variable
for v = 1:length(varplot)

    % for each muscle
    for m = 1:length(muscplot)
        
        fig = figure;
        fig.Position = [1137 546 426 180];

        % for each side
        for s = 1:2

            sp = subplot(1,2,s);

            % for each simulation
            for k = 1:n

                musc = [muscplot{m} '_' side{s}];
                var = varplot{v};

                plot(x,muscle.(musc).(var)(k,:),'Color',color_spec(k,:),'LineStyle',linestyle{k},'LineWidth',linewidth(k));
                hold on

                if k > 1; fprintf('%s, %s %s, %s RMSE: %f %s\n',varplot{v},side{s},muscplot{m},labels{k},rms(muscle.(musc).(var)(k,:)-muscle.(musc).(var)(1,:)),units{v}); end

            end

            sp.Box = 'off';
            title([vtitle{v} ': ' side{s} ' ' mtitle{m}])
            if s == 1; ylabel(ylab{v}); end
            xlabel('% Stride')

        end

    end

end
   
%% TIME SERIES PLOTS: JOINT MOMENT

% joints to plot
jointplot = {'knee','ankle'};

% titles
jtitle = {'Ft','Pm'};
        
% for each joint
for j = 1:length(jointplot)
        
    fig = figure;
    fig.Position = [1137 546 426 180];

    % for each side
    for s = 1:2

        sp = subplot(1,2,s);

        % for each simulation
        for k = 1:n

            thisMoment = [side{s} jointplot{j}];

            plot(x,moment.(thisMoment)(k,:),'Color',color_spec(k,:),'LineStyle',linestyle{k},'LineWidth',linewidth(k));
            hold on
            
            if k > 1; fprintf('Joint moment, %s %s, %s RMSE: %f Nm\n',side{s},jointplot{j},labels{k},rms(moment.(thisMoment)(k,:)-moment.(thisMoment)(1,:))); end

        end
        
        sp.Box = 'off';
        if s == 1; ylabel('Moment (Nm)'); end
        xlabel('% Stride')
        title([side{s} ' ' jointplot{j}])

    end

end
   
%% DISCRETE VAR HISTOGRAMS: MUSCLE WORK

% muscles to plot
muscplot = {'vasti','rect_fem','gastroc','soleus'};

% variables to plot
varplot = {'net_eccentric_work_stance','net_concentric_work_stance'};

% titles
vtitle = {'Ecc','Con'};
mtitle = {'Vasti','RF','Gastroc','Soleus'};

% axis labels
ylab = {'Work (J)','Work (J)'};
        
% for each variable
for v = 1:length(varplot)
    
    var = varplot{v};

    % for each muscle
    for m = 1:length(muscplot)
        
        fig = figure;
        fig.Position = [954 227 238 127];

        % for each side
        for s = 1:2
            
            sp = subplot(1,2,s);
            
            musc = [muscplot{m} '_' side{s}];
            
            % for each simulation
            discvar = zeros(1,n);
            for k = 1:n
                discvar(k) = muscle.(musc).(var)(k); 
                if k > 1
                    err = discvar(k) - discvar(1);
                    fprintf('%s, %s %s, %s Error = %f J (%f%%)\n',varplot{v},side{s},muscplot{m},labels{k},err,err/discvar(1)*100)
                end
            end
            
            b = bar(1:n,discvar);
            b.FaceColorMode = 'manual';
            b.FaceColor = 'flat';
            for k = 1:n; b.CData(k,:) = color_spec(k,1:3); end
            b.FaceAlpha = 0.7;
            b.EdgeColor = [1 1 1];
            sp.Box = 'off';
            sp.XTickLabel = labels;
            sp.XTickLabelRotation = 45;
            sp.YAxis.LineWidth = 1.5;
            sp.XAxis.LineWidth = 1.5;
            sp.XAxis.TickLength = [0 0];
            sp.YAxis.LineWidth = 1.5;
            sp.XAxis.LineWidth = 1.5;
            title([vtitle{v} ': ' side{s} ' ' mtitle{m}])
            ylabel(ylab{v})
            
        end
        
    end

end
