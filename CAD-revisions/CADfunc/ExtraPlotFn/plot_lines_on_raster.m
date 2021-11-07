function [ flag ] = plot_lines_on_raster( vlines,vtype,spM,rasterpm )
% plot_lines_on_raster
%
%   INPUTS:
%       vlines  =   x-location of vertical lines
%       vtype   =   integer index classifying each line
%       spM     =   spike matrix (used to decide on height of plot)
%       rasterpm    = optional parameters
%                  


%
% Default values
UseUnitOrder        = 0;           % Use provided unit order?
PlotParticipatingUnitsOnly = 0;    % Plot only units that participate in at least one assembly?
LineColorVec        = [0.5 0.5 0.5];  % Matrix giving colors of vertical lines

if (isfield(rasterpm,'UseUnitOrder')); UseUnitOrder = rasterpm.UseUnitOrder; end
if (isfield(rasterpm,'PlotParticipatingUnitsOnly')); PlotParticipatingUnitsOnly = rasterpm.PlotParticipatingUnitsOnly; end;
if (isfield(rasterpm,'LineColorVec')); LineColorVec = rasterpm.LineColorVec; end
%if (isfield(rasterpm,'color_assign_type'));color_assign_type = rasterpm.color_assign_type;end;


if (UseUnitOrder || PlotParticipatingUnitsOnly)
    try
        Unit_order = rasterpm.Unit_order;
    catch
        warning('raster_all_assemblies_fn: Unit_order is required for UseUnitOrder and PlotParticipatingUnitsOnly')
    end
    % FOR NOW:  PlotParticipatingUnitsOnly => UseUnitOrder
    UseUnitOrder = 1;
end

if (PlotParticipatingUnitsOnly)
    nCell = find(diff(Unit_order)<0,1,'last');
else
    nCell = size(spM,1);
end

for j1=1:length(vlines) 
    xv      = vlines(j1);
    plot([xv xv],[0.5 nCell+0.5],'-','Color',LineColorVec(vtype(j1),:));
end

flag = 0;
end

