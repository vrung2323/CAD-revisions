function [ h ] = plot_assembly_assignment_matrix( Amatrix, Binvector, ...
    Unit_order, As_order, BinSizes, AAMpm)
%  
% - assembly assigment matrix: plot ONLY
% - REMOVE any cells which are not participating in any assembly
%
%   INPUTS:
%   Amatrix, Binvector, Unit_order, As_order: as produced by assembly_assignment_matrix
%   BinSizes:   Needed to set range for pcolor
%
%   AAMpm:  Parameters (for the Assembly Assignment Matrix i.e. AAM) to
%   govern plotting


%% Default values
% Put assembly matrix and time bins on separate plot?
SeparatePlotsFlag = 0;

% Also plot Idx  
%   A vector of cell classificiations
%   (For example, in Xu Lab project: is the cell R, DR, or C?)
PlotIdxFlag = 0;

% Only plot cells that appear in at least one assembly?
PlotParticipatingUnitsOnly = 0;

if (isfield(AAMpm,'SeparatePlotsFlag')); SeparatePlotsFlag = AAMpm.SeparatePlotsFlag; end;
if (isfield(AAMpm,'PlotIdxFlag'))
    PlotIdxFlag = AAMpm.PlotIdxFlag; 
    if (PlotIdxFlag)
        try
            AllIdx = AAMpm.AllIdx;
        catch 
            error('plot_assembly_assignment_matrix: PlotIdx=1 but AllIdx is not present');
        end
    end
end
if (isfield(AAMpm,'PlotParticipatingUnitsOnly')); PlotParticipatingUnitsOnly = AAMpm.PlotParticipatingUnitsOnly; end;
% if (nargin > 5)
%     narg = nargin-5;
%     while (narg > 0)
%         keyw   = varargin{1};
%         kvalue = varargin{2};
%         if (strcmp(keyw,'SeparatePlots'))
%             SeparatePlotsFlag = kvalue;
%         elseif (strcmp(keyw,'AllIdx'))
%             PlotIdxFlag=1;
%             AllIdx = kvalue;
%         else
%             warning('plot_assembly_assignment_matrix: keyword not recognized'); 
%         end
%         narg = narg-2;
%         varargin=varargin(3:end);
%     end
% end

nneu = size(Amatrix,1);

if (PlotParticipatingUnitsOnly)
    % Find all units that participate in at least 1 assembly
    keep_vector     = sum(~isnan(Amatrix),2);
    keep_index      = find(keep_vector);
else
    keep_index      = 1:nneu;  % Keep all units
end

% Don't plot the others!
nneu_temp  =  length(keep_index);

Amatrix_temp    = Amatrix(keep_index,:);
Unit_order_temp = Unit_order(keep_index);

% IF include ID, do it now
if (PlotIdxFlag)
    AllIdx_temp=AllIdx(Unit_order_temp);
    % Amatrix_temp = [Amatrix_temp AllIdx(Unit_order_temp)];
    figure;
    AllIdx_temp(end+1,end+1)=0;
    pcolor(AllIdx_temp);
    caxis([0, max(AllIdx_temp(:))])
    colormap jet
    ylabel('Cell class (R, DR, C, etc.)');
    colorbar
end
        
binmat=log10(Binvector);

% This is needed to make "pcolor" work, for some reason
Amatrix_temp(end+1,end+1)=0; 
binmat(end+1,end+1)=0;


if (SeparatePlotsFlag)   
    figure;
else
    subplot(2,1,1);
end



h=pcolor(Amatrix_temp);
set(h, 'EdgeColor', [0.8 0.8 0.8]);
caxis([0, max(Amatrix_temp(:))])
colormap jet
ylabel('Unit #')
set(gca, 'XTick', []);
hcb=colorbar;
hcb.Label.String = 'Time lag \it l \rm (# bins)';
yy = 1:size(Amatrix_temp,1);
set(gca,'YTick',yy(1)+0.5:yy(end),'Yticklabel',Unit_order_temp) 
set(gcf, 'Color', [1,1,1]);

        
        
        

if (SeparatePlotsFlag)
    figure;colormap(jet);
else
    subplot(2,1,2)
end

pcolor(binmat)
yy = 1:size(Amatrix_temp,2);
set(gca,'XTick',yy(1)+0.5:2:yy(end)+1+0.05,'Xticklabel',As_order(1:2:end)) 
xlabel('Assembly #')
set(gca, 'YTick', []);
hC = colorbar;
if length(BinSizes)>1
    caxis([log10(min(BinSizes)), log10(max(BinSizes))]);
else
    caxis([log10(BinSizes-0.001), log10(BinSizes+0.001)]);
end
L=[0.001,0.002,0.003,0.004,0.005, 0.006,0.007,0.008,0.009,0.01,0.02,0.03,...
    0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,...
    1,2,3,4,5,6,7,8,9,10];
set(hC,'Ytick',log10(L),'YTicklabel',L);
set(hC,'Location','southoutside')
hC.Label.String = 'Temporal precision \Delta (sec)';

 

end

