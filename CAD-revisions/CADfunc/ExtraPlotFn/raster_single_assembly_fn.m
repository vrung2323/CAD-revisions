function [spike_out_all ] = raster_single_assembly_fn( As_across_bins_single, ...
    AAasspikes_single,show,varargin )
% raster_single_assembly_fn: plot spike times of a single assembly
%
%   INPUTS:  
%   As_across_bins_single:    one member of an assembly structure cell, 
%                               as output by  "pruning_across_bins"
%                             (i.e. this argument should be 
%                             As_across_bins_pr{j}, for some j)
%   AAasspikes_single:         one member of a spike train cell, as produced
%                               by "assembly_raster_plots"
%                                (i.e. this argument should be 
%                                AAaspikes{j}, for the same j)
%                                       
%   show:                      [0/1] whether to plot
%   
%   OPTIONAL:  
%     Keyword/Argument pairs
%   'Color':                    ['r'] Color to pass to the plot function
% 

% Process optional arguments

% Color for raster plot
spcolor = 'r';

% Use supplied unit order?
UseUnitOrder = 0;

if (nargin > 3)
    narg = nargin-3;
    while (narg > 0)
        keyw   = varargin{1};
        kvalue = varargin{2};
        if (strcmp(keyw,'Color'))
            spcolor = kvalue;
        elseif (strcmp(keyw,'UseUnitOrder'))
            UseUnitOrder = 1;
            Unit_order = kvalue;
        else
            warning('raster_single_assembly_fn: keywork not recognized'); 
        end
        narg = narg-2;
        varargin=varargin(3:end);
    end
end

CellUnits   = As_across_bins_single.elements;
nCell       = length(CellUnits);

% Save spikes in the format: [unit# time]
spike_out_all = [];
for j1=1:nCell
    aus = AAasspikes_single(j1,:);
    aus(isnan(aus))=[];
    
    if (UseUnitOrder)
        % Get CellID from 
        CellID = find(Unit_order==CellUnits(j1),1);
    else
        CellID = CellUnits(j1);
    end
    spike_out_all = [spike_out_all; ...
        [CellID*ones(length(aus),1) aus'] ];
end

if (show)
    plot(spike_out_all(:,2),spike_out_all(:,1),'.',...
        'Color',spcolor,'MarkerSize',10);
end


end

