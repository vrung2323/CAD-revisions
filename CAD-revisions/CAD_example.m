%%%%%%%%%%%%%%%%%%%%%%%%%% TUTORIAL ON CELL ASSEMBLY DETECTION %%%%%%%%%%%%%%%%%%%%%%%%%
close all; 
%load('testData/test_spike_data.mat');  % Contains spM

addpath('CADfunc/');
addpath('CADfunc/CreateTestData');
addpath('CADfunc/ExtraPlotFn');
%%

nneu = size(spM,1);
BinSizes=[0.015 0.025 0.04 0.06 0.085 0.15 0.25 0.4 0.6 0.85 1.5];


MaxLags=10*ones(size(BinSizes));
% Restrict MaxLags so that total assembly length is less than ~=2 s
%MaxLags = min(MaxLags,ceil(2./BinSizes));


% This is a threshold on E[AB]. Below this threshold, pairs will NOT be
%   tested for correlations.
Exp_thres = 1;



%% %%%%%%%%%%%%%%%%%%%% ASSEMBLY DETECTION %%%%%%%%%%%%%%%%%%%%%%%

% 
% Use "plotOnly" if you have already run the algorithm, and just 
%  want to experiment with the visualizations/pruning options.
%
plotOnly = 0;
    if (plotOnly)
        load('test_CAD_example.mat','assembly');
    else
        [assembly]=Main_assemblies_detection(spM,MaxLags,BinSizes, Exp_thres,...
            0,0.05,100,i);

        %% If you want to save your results for future reference
        save('test_CAD_example.mat','assembly');
    end

    %% %%%%%%%%%%%%%%%%%%%%%%%% VISUALIZATION %%%%%%%%%%%%%%%%%%%%%%%%

    % ASSEMBLY REORDERING
    figure
    [As_across_bins,As_across_bins_index]=assemblies_across_bins(assembly,BinSizes);
    display='raw';
    %display='clustered';
    % VISUALIZATION
    [Amatrix,Binvector,Unit_order,As_order]=...
    assembly_assignment_matrix(As_across_bins, nneu, BinSizes, display,1);


%% %%%%%%%%%%%%%%%%%%%%%%%% PRUNING %%%%%%%%%%%%%%%%%%%%%%%%

prune = input('Type B to prune by size, P to prune by pvalue, or O to prune by occurrences: ','s');
while prune ~= 'B' && prune ~= 'P' && prune ~= 'O'
    prune = input('ERROR: Type B to prune by size, P to prune by pvalue, or O to prune by occurrences: ','s');
end
figure(2)

if prune == 'B'
    % PRUNING: criteria = 'biggest';
    criteria = 'biggest';
    [As_across_bins_pr,As_across_bins_index_pr]=...
    pruning_across_bins(As_across_bins,As_across_bins_index,nneu,criteria);
    display='raw';
    [Amatrix,Binvector,Unit_order,As_order]=...
    assembly_assignment_matrix(As_across_bins_pr, nneu, BinSizes, display,1);
    title(sprintf('Pruning: %s',criteria));
elseif prune == 'P'
    % PRUNING: criteria = 'distance', style = 'pvalue';
    criteria = 'distance'; style = 'pvalue'; th = 0.3;
    [As_across_bins_pr,As_across_bins_index_pr]=...
    pruning_across_bins(As_across_bins,As_across_bins_index,nneu,criteria,th,style);
    display='raw';
    [Amatrix,Binvector,Unit_order,As_order]=...
    assembly_assignment_matrix(As_across_bins_pr, nneu, BinSizes, display,1);
    title(sprintf('Pruning: %s, %s',criteria,style));
elseif prune == 'O'
    % PRUNING: criteria = 'distance', style = 'occ';
    criteria = 'distance'; style = 'occ'; th = 0.3;
    [As_across_bins_pr,As_across_bins_index_pr]=...
    pruning_across_bins(As_across_bins,As_across_bins_index,nneu,criteria,th,style);
    display='raw';
    [Amatrix,Binvector,Unit_order,As_order]=...
    assembly_assignment_matrix(As_across_bins_pr, nneu, BinSizes, display,1);
    title(sprintf('Pruning: %s, %s',criteria,style));
end

clc
disp('This is the PRUNED Assembly Assignment Matrix.');
disp('It contains a subset of assemblies; in this case ');
disp('  it discards any assemblies which are a proper ');
disp('  subset of another (B), or in the case of');
disp('  pruning by distance, groups assemblies within');
disp('  a cosine distance, either choosing the most');
disp('  significant assembly (P) or the one with the');
disp('  most occurrences (O). ');
disp('If the same cell configuration is identified at multiple');
disp('  time scales, the most significant is maintained');
disp(' ');
disp('Press ENTER to continue');
pause;




%% %%%%%%%%%%%%%%%%%%%%%%%% ASSEMBLY ACTIVATION %%%%%%%%%%%%%%%%%%%%%%%%
figure(3)

lagChoice = 'beginning';
% lagChoice='duration';

act_count = 'full';
[assembly_activity]=Assembly_activity_function(As_across_bins_pr,assembly,...
    spM,BinSizes,lagChoice,act_count);

%% Plot assembly activity over time
for i=1:min(length(assembly_activity),5)
    subplot(5,1,i)
    plot(assembly_activity{i}(:,1),assembly_activity{i}(:,2));
    hold on
end

clc
disp('This is (pruned) assembly activity as a function of time');
disp('  "Activity" means how many times the assembly was activated');
disp('  in each time bin');
disp(' ');
disp('Press ENTER to continue');
pause;


%% %%%% Another style of raster plot: view each assembly individually %%

figure(4)
nr = numel(As_across_bins_pr);
for jj = 1:numel(As_across_bins_pr)
    subplot(nr,1,jj);
    %if numel(As_across_bins_pr{jj}.elements)==5
     view_assembly_activity_fn(As_across_bins_pr(jj),spM);
   % end
end

clc
disp('Here we view each (pruned) assembly as separate raster plot');
disp('  This can be useful when assemblies overlap in cell membership.');
disp('  ');
disp('Press ENTER to continue');
pause;

%% %%%%%%%%%%%%%%% RASTER PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    First use RD17 function to extract assembly-linked spikes
[AAspM,AAasspikes]=assembly_rasterplot(As_across_bins_pr,assembly_activity,spM,0);
    
%% Set up parameters for raster plots
rasterpm.UseUnitOrder = 1;
rasterpm.Unit_order = Unit_order;
rasterpm.PlotParticipatingUnitsOnly = 0;
rasterpm.color_assign_type = 1;
rasterpm.PauseBetweenAssemb = 0;
% 
% % This will rotate through colors: each assembly in a different color
% hr1  = raster_all_assemblies_fn( spM,AAasspikes,As_across_bins_pr,...
%     As_order,BinSizes,rasterpm);
% set(gca,'FontSize',14); title('Without adjustment');
% axis tight; ylim([0.5,nneu+0.5]);

%% %%%%%%%%%%%%%%% RASTER PLOT: ADJUSTED ACTIVITY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Subtract off "expected" assembly activity due to 
% background firing rate
[assembly_activity_adj]=Adjust_AsAct(assembly_activity,As_across_bins_pr,assembly,...
    spM,BinSizes);

% Find spikes associated w/ the assembly
[AAspMAdj,AAasspikesAdj]=assembly_rasterplot(As_across_bins_pr,...
    assembly_activity_adj,spM,0);

hr2  = raster_all_assemblies_fn( spM,AAasspikesAdj,As_across_bins_pr,...
    As_order,BinSizes,rasterpm);
set(gca,'FontSize',14); title('With adjustment');
axis tight; ylim([0.5,nneu+0.5]);

clc
disp('This is a raster plot showing the activity of each (pruned)'); 
disp('   assembly with a different color');
disp('  ');

disp('Press ENTER to continue');

pause;

xlim([50,100]);
clc
disp('Zoom in using xlim to get a closer look at the assembly structure');