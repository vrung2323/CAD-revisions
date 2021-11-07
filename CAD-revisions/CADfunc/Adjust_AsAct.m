function [assembly_activity_adj] = Adjust_AsAct(assembly_activity, As_across_bins, ...
    assembly, spM, BinSizes)
% Adjust assembly activity to compensate for baseline firing rate
% 
%ARGUMENTS:
% assembly_activity := output of Assembly_activity_function
% As_across_bins := structure containing assembly information output of ASSEMBLIES_ACROSS_BINS;
%       assembly :=  structure containing assembly information output of Main_assemblies_detection.m
%            spM := matrix with population spike trains; each row is the spike train (time stamps, not binned) relative to a unit. 
%       BinSizes := vector of bin sizes to be tested;
% 


as=As_across_bins;
assembly_activity_adj=cell(length(as),1);
 
for n=1:length(as)
   % Time bin for assembly
   tbin=as{n}.bin;
   
   % Find EXPECTED # of coincidences, based on
   %  avg spikes per bin
   nCell  = length(as{n}.elements);
   avgSp  = zeros(nCell,1);
   for j1=1:nCell
       ee = as{n}.elements(j1);
       avgSp(j1) = sum(~isnan(spM(ee,:)))/length(as{n}.Time);
   end
   %If p1,p2,...pn are Poisson with mean u1,u2,...un,
   %   what is E(min(p1,p2,...,pn))
   %
   %AS A PROXY: USE min(u1,u2,...un)
   %avgSp'
   mAdj = min(avgSp);
   temp = assembly_activity{n};
   temp(:,2)=max(temp(:,2)-mAdj,0);
   assembly_activity_adj{n}=temp;
   
    
   % Time bin edges
   %tb=assembly.bin{find(BinSizes==as{n}.bin)}.bin_edges; 
   
end

end

