%%%%%%%%%%%%%%%%%%%%%%%%%% TUTORIAL ON CELL ASSEMBLY DETECTION %%%%%%%%%%%%%%%%%%%%%%%%%
close all; 
%load('testData/test_spike_data.mat');  % Contains spM

addpath('CADfunc/');
addpath('CADfunc/CreateTestData');
addpath('CADfunc/ExtraPlotFn');
%%

nneu = size(spM,1);
BinSizes=[0.015 0.025 0.04 0.06 0.085 0.15 0.25 0.4 0.6 0.85 1.5];
BinPhases = linspace(-0.4,0.5,10);
%BinPhases = [-0.1 0 0.1];
As_across_bins_saved = cell(length(BinPhases),1);


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
for i = 1:length(BinPhases)
    if (plotOnly)
        load('test_CAD_example.mat','assembly');
    else
        [assembly]=Main_assemblies_detection(spM,MaxLags,BinSizes, Exp_thres,...
            0,0.05,100,BinPhases(i));

        %% If you want to save your results for future reference
        save('test_CAD_example.mat','assembly');
    end

    %% %%%%%%%%%%%%%%%%%%%%%%%% VISUALIZATION %%%%%%%%%%%%%%%%%%%%%%%%

    % ASSEMBLY REORDERING
    figure;
    [As_across_bins,As_across_bins_index]=assemblies_across_bins(assembly,BinSizes);
    display='raw';
    %display='clustered';
    % VISUALIZATION
    [Amatrix,Binvector,Unit_order,As_order]=...
    assembly_assignment_matrix(As_across_bins, nneu, BinSizes, display,1);
    title(sprintf('Bin Phase: %.1f',BinPhases(i)));
    As_across_bins_saved{i} = As_across_bins;
end

output_data = lag_search(As_across_bins_saved);