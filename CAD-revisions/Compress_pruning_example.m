close all; clear; clc;

addpath('CADfunc/');
addpath('CADfunc/CreateTestData');
addpath('CADfunc/ExtraPlotFn');
%%
%% generate data with compressed pattern
n = [100 100];
T = [1 0.2];
temp = generating_ISI_fn(10,500,0.5);
data = embed_alt_data_fn(temp,n,T,[1 2 3 4 5],'twoblock');
spM = data.spM;


%% parameters
nneu=size(spM,1);  % nneu is number of recorded units
BinSizes=[0.015 0.025 0.04 0.06 0.085 0.15 0.25...
    0.4 0.85 1 1.5];
reflag = 0;
ExpThres=1;
MaxLags= 10* ones(size(BinSizes));


%% compute assembly before pruning
% ASSEMBY DETECTION
[assembly] = Main_assemblies_detection(spM,MaxLags,BinSizes,...
  ExpThres,reflag);

fprintf('Complete Detecting Assemblies \n');

% ASSEMBLY REORDERING
[As_across_bins,As_across_bins_index] = ...
  assemblies_across_bins(assembly,BinSizes);

fprintf('Complete Reordering \n');

%% Pruning
% original pruning
criteria = 'biggest';
% criteria = 'distance';
[As_across_bins_pr,As_across_bins_index_pr]= ...
  pruning_across_bins(As_across_bins,As_across_bins_index,nneu,...
  criteria);

fprintf('Complete original pruning \n');

%compress pruning
PrComprAs = compress_pruning_fn(As_across_bins,spM,BinSizes);

fprintf('Complete compress_pruning \n');

%% %%%%%%%% Produce the results %%%%%%%%%%%%%%%%%%%
if 1
% plot the true activity, USE Scratch.m file
% acttime = Data.act_time;
% actlabel = (acttime(5,:) - acttime(1,:))<1;
% rasterplot_fn(spM,{'spk',acttime});

% plot Amatrix
% VISUALIZATION
figure(1)
[Amatrix,Binvector,Unit_order,As_order]= ...
  assembly_assignment_matrix(As_across_bins, nneu, BinSizes,'raw',1);
% set(gcf, 'Position', get(0, 'Screensize'));


% plot original algorithm result
figure(2)
dis='raw';
[Amatrix,Binvector,Unit_order,As_order]= ...
  assembly_assignment_matrix(As_across_bins_pr, nneu, BinSizes, dis,1);
% set(gcf, 'Position', get(0, 'Screensize'));
% export_fig(sprintf('%sop%d.pdf',resfld,id),'-native');
% close all;

figure(3)
for jj = 1:numel(As_across_bins_pr)
    if numel(As_across_bins_pr{jj}.elements)==5
        view_assembly_activity_fn(As_across_bins_pr(jj),spM)
    end
end
title('Activities from original method')
% % set(gcf, 'Position', get(0, 'Screensize'));
% export_fig(sprintf('%sorgact%d.pdf',resfld,id),'-native');
% close all;

% plot compress pruning result
figure(4)
Amat = zeros(size(spM,1),numel(PrComprAs));
for jj = 1:numel(PrComprAs)
    Amat(PrComprAs{jj}.elements,jj)=1;
    if numel(PrComprAs{jj}.elements)==...
        max(cellfun(@(c) numel(c.elements),PrComprAs))
         %export activity of the LARGEST compressed assembly
        rasterplot_fn(spM(PrComprAs{jj}.elements,:),...
            {'spk',PrComprAs{jj}.Act});
            %     set(gcf, 'Position', get(0, 'Screensize'));
            %     export_fig(sprintf('%scompressact%d.pdf',resfld,id),'-native');
            %     close all;
        %     pause;
    end
end
title('Activities from compress pruning');

end
