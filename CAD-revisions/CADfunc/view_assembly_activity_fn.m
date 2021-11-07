function view_assembly_activity_fn(Assembly,spM)
% plot activity of each assemblies in cell-array Assembly
%spM is the spike stamp of all neurons

% f = figure('Position',[0 0 1000 1000]);
for i = 1:numel(Assembly)
    A = Assembly{i};
    tb=min(spM(:)):A.bin:max(spM(:));
    rasterplot_fn(spM(A.elements,:),{'bin',tb,A.Time,A.lag});
    m = numel(A.elements);
    str = sprintf('bin = %f - Elements = ',A.bin);
    for j = 1:m
        str = sprintf('%s %d',str,A.elements(j));
    end
    str = sprintf('%d/%d %s',i,numel(Assembly),str);
    title(str); 
    %   name = sprintf('Results/assembly%d.png',i);
    %   export_fig(name);
    %   pause;clf;
end