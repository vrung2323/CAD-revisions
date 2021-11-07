function [ h ] = raster_all_assemblies_fn( spM,AAasspikes,As_across_bins_pr,...
    As_order,BinSizes,rasterpm)
% raster_all_assemblies_fn: Raster plot with each assembly in a different
%               color
%   INPUTS:
%   spM:                All spikes
%   AAasspikes:         Spikes associated with each assembly
%   As_across_bins_pr: 
%   As_order:
%   BinSizes:   
%
%   rasterpm:   Parameters
%
%

%
% Default values
UseUnitOrder        = 0;           % Use provided unit order?
PlotParticipatingUnitsOnly = 0;    % Plot only units that participate in at least one assembly?
PauseBetweenAssemb  = 0;           % Pause between asemblies?
% color_assign_type
%    1 =        iterate through assemblies
%    2 =        by Bin Size
color_assign_type = 1;

if (isfield(rasterpm,'UseUnitOrder')); UseUnitOrder = rasterpm.UseUnitOrder; end;
if (isfield(rasterpm,'PlotParticipatingUnitsOnly')); PlotParticipatingUnitsOnly = rasterpm.PlotParticipatingUnitsOnly; end;
if (isfield(rasterpm,'PauseBetweenAssemb')); PauseBetweenAssemb = rasterpm.PauseBetweenAssemb; end;
if (isfield(rasterpm,'color_assign_type'));color_assign_type = rasterpm.color_assign_type;end;
if (isfield(rasterpm,'axisHandle'))
    figAx = rasterpm.axisHandle;  axes(figAx);
else
    % Make new figure
    figure;
end

if (UseUnitOrder || PlotParticipatingUnitsOnly)
    try
        Unit_order = rasterpm.Unit_order;
    catch
        error('raster_all_assemblies_fn: Unit_order is required for UseUnitOrder and PlotParticipatingUnitsOnly');
    end
    % FOR NOW:  PlotParticipatingUnitsOnly => UseUnitOrder
    UseUnitOrder = 1;
end

% How many assemblies?
nAssemb = size(As_across_bins_pr,2);

% Assign colors to each assembly
color_list = zeros(nAssemb,1);
cmap = colormap(lines);
if (color_assign_type == 1)
    if (nAssemb < 64)
        % Interpolate
        nColors    = nAssemb;
        color_list = [1:nAssemb]';
       
        intval = floor(64/(nColors-1));
        indvec = 1:intval:64; 
        if (length(indvec)<nColors); indvec = [indvec 64]; end;
    else
        % Cycle through 64 colors
        color_list = [1:nAssemb]';
        indvec     = mod(color_list,64);
        indvec(find(indvec==0))=64;
    end
elseif (color_assign_type == 2)
    nColors     = length(BinSizes);
    
    % Pick colors to be the same as assembly matrix plot
    Binlog10=log10(BinSizes);
    minBinlog10=log10(min(BinSizes)); maxBinlog10=log10(max(BinSizes)); 
    
    indvec = round(((Binlog10-minBinlog10)/(maxBinlog10-minBinlog10))*(64-1))+1;
   
    for j1=1:nAssemb
        color_list(j1) = find(BinSizes==As_across_bins_pr{As_order(j1)}.bin);
    end
end
color_vecs = cmap(indvec,:);

% We begin to plot
h=gca;

if (PlotParticipatingUnitsOnly)
    nCell = find(diff(Unit_order)<0,1,'last');
else
    nCell = size(spM,1);
end

for j1=1:nCell
    j1ind = j1;
    if (UseUnitOrder); j1ind = Unit_order(j1); end;
    aus=spM(j1ind,:);
    
    aus(isnan(aus))=[];
    plot(aus,j1*ones(size(aus)),'.','color',[0.5 0.5 0.5],'MarkerSize',6);
    hold on;    
end

for k1=1:nAssemb
    j1 = As_order(k1);
    if (UseUnitOrder)
        spikes_out_all = raster_single_assembly_fn(As_across_bins_pr{j1},...
        AAasspikes{j1},1,...
        'Color',color_vecs(color_list(k1),:),...
        'UseUnitOrder',Unit_order);
    else
        spikes_out_all = raster_single_assembly_fn(As_across_bins_pr{j1},...
        AAasspikes{j1},1,...
        'Color',color_vecs(color_list(k1),:));
    end 
    
    % Should print some info Re: assmebly as well...
    if (PauseBetweenAssemb); pause; end;
end

end
