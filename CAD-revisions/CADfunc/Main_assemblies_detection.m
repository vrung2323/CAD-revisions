function [assembly]=Main_assemblies_detection(spM,MaxLags,BinSizes,...
  Exp_Thres, ref_lag,alph,Dc,p_shift,No_th,O_th,bytelimit)
% this function returns cell assemblies detected in spM spike matrix binned
% at a temporal resolution specified in 'BinSizes' vector and testing for all
% lags between '-MaxLags(i)' and 'MaxLags(i)'
%
% USAGE: [assembly]=Main_assemblies_detection(spM, MaxLags, BinSizes,
% Exp_Thres,
% ref_lag, alph, Dc, No_th, O_th, bytelimit)
%
% ARGUMENTS:
% spM       := matrix with population spike trains; each row is the spike
%       train (time stamps, not binned) relative to a unit.
% MaxLags   := vector of maximal lags to be tested. For a binning dimension of
%       BinSizes(i) the program will test all pairs configurations with a time
%       shift between -MaxLags(i) and MaxLags(i);
% BinSizes  := vector of bin sizes to be tested;
% Exp_Thres := threshold on E[AB]
% (optional) ref_lag    := reference lag. Default value 0
% (optional) alph       := alpha level. Default value 0.05
% (optional) Dc         := # chunks for coincidence detection. 
%                               Default value 100
% (optional) p_shift    := phase shift for bins to check. 
%                               Default value 0
% (optional) No_th      := minimal number of occurrences required for an
%       assembly (all assemblies, even if significant, with fewer occurrences
%       than No_th are discarded). Default value 0.
% (optional) O_th       := maximal assembly order (the algorithm will return
%       assemblies of composed by maximum O_th elements).
% (optional) bytelimit  := maximal size (in bytes) allocated for all
%       assembly structures detected with a bin dimension. When the size limit is
%       reached the algorithm stops adding new units.
%
% RETURNS:
% assembly - structure containing assembly information:
%     assembly.parameters       - parameters used to run Main_assemblies_detection
%     assembly.bin{i} contains information about assemblies detected with
%                     'BinSizes(i)' bin size tested for all lags between
%                     '-MaxLags(i)' and 'MaxLags(i)'
%
%        assembly.bin{i}.bin_edges - bin edges (common to all assemblies in
%        assembly.bin{i})
%        assembly.bin{i}.n{j} information about the j-th assembly detected
%        with BinSizes(i) bin size
%                   elements: vector of units taking part to the assembly
%                   (unit order correspond to the agglomeration order)
%                        lag: vector of time lags. '.lag(z)' is the
%                        activation delay between .elements(1) and
%                        .elements(z+1)
%                         pr: vector of pvalues. '.pr(z)' is the pvalue of
%                         the statistical test between performed adding
%                         .elements(z+1) to the structure .elements(1:z)
%                       Time: assembly activation time. It reports how many
%                       times the complete assembly activates in that bin.
%                       .Time always refers to the activation of the first
%                       listed assembly element (.elements(1)), that
%                       doesn't necessarily corresponds to the first unit
%                       firing.
%               Noccurrences: number of assembly occurrence.
%               '.Noccurrences(z)' is the occurrence number of the
%               structure composed by the units .elements(1:z+1)
%
%
%
%  ?? 2016 Russo, Durstewitz.
%  for information please contact eleonora.russo@zi-mannheim.de; daniel.durstewitz@zi-mannheim.de.
%
%  last update 11/01/2016
%  last update 13/02/2016 TestPair.m update
%
%  2020 Modifications: Truong, Barreiro, Gande/Habig
%

% NOTE: Very important this be 0, not 2
if nargin<5 || isempty(ref_lag), ref_lag=0; end

if nargin<6 || isempty(alph), alph=0.05; end
if nargin<7 || isempty(Dc), Dc=100; end
if nargin<8 || isempty(p_shift), p_shift=0; end
if nargin<9 || isempty(No_th), No_th=0; end      % no limitation on the number of assembly occurrences
if nargin<10 || isempty(O_th), O_th=Inf; end     % no limitation on the assembly order (=number of elements in the assembly)
if nargin<11 || isempty(bytelimit), bytelimit=Inf; end     % no limitation on assembly dimension

nneu=size(spM,1); % number of units
assemblybin=cell(1,length(BinSizes));

% We now send this in as an optional argument
%Dc=100; %length (in # bins) of the segments in which the 
%       spike train is divided to compute #abba variance (parameter k).

% Make a folder to save temporary files
foldname = num2str(randi([100 9999]));
while exist(foldname,'dir')
    foldname = num2str(randi([100 9999]));
end
mkdir(foldname);

%start calculating
fprintf('Bin Phase = %.1f \n', p_shift);
for gg=1:length(BinSizes)
%parfor gg=1:length(BinSizes)
    int=BinSizes(gg);
    maxlag=MaxLags(gg);
    fprintf('%d - testing: bin size=%f sec; max tested lag=%d \n',...
        gg, int, maxlag);

    %%%%%%%%%%%%%%%%%%%%% Binning  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % the spike train is binned with temporal resolution 'int'
    tb=min(spM(:)):int:max(spM(:)); % keeps bin edges consistent according to histcounts
    tb = tb + p_shift*BinSizes(gg);
    binM=zeros(nneu,length(tb)-1,'uint8');
  
    for n=1:nneu
        [ binM(n,:),~] = histcounts(spM(n,:),tb);
    end
  
    if size(binM,2)-MaxLags(gg)<100
        fprintf('Warning: testing bin size=%f. The time series is too short, consider taking a longer portion of spike train or diminish the bin size to be tested \n', int);
    else
        %%%%%%%%%%%%%%%%%%%%%%% Analysis  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        [assemblybin{gg}]=FindAssemblies_recursive(binM,maxlag,Exp_Thres,alph,...
            gg,Dc,No_th,O_th,bytelimit,ref_lag,foldname);  
            % it returns assemblies detected at specific temporal resolution 'int'
    
        if ~isempty(assemblybin{gg})
            assemblybin{gg}.bin_edges=tb;
        end
        fprintf('%d - testing done\n',gg);
        fname = sprintf('%s/assembly%d.mat',foldname, gg);
        parsave(fname,assemblybin{gg})
    end
end
assembly.bin=assemblybin;
assembly.parameters.alph=alph;
assembly.parameters.Dc=Dc;
assembly.parameters.No_th=No_th;
assembly.parameters.O_th=O_th;
assembly.parameters.bytelimit=bytelimit;
assembly.parameters.Exp_Thres=Exp_Thres;

fprintf('\n');
%remove temporary folders
rmdir(foldname,'s');
end

function parsave(fname,aus)
save(fname,'aus')
end
