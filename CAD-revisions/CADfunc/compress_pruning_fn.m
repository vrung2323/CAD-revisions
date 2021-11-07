function PrComprAs = compress_pruning_fn(As_across_bins,spM,BinSizes)

nneu = size(spM,1);

[Amatrix,binvec] = ...
  assembly_assignment_matrix(As_across_bins, nneu, BinSizes, 'raw',0);

%using Amatrix to find compressed assembly
% now we sort the assembly by size
[~,nAs] = size(Amatrix);

[~,As_order] = sort(sum(~isnan(Amatrix),1),'ascend'); %rearrange the assembly based on size

sAmat = ~isnan(Amatrix(:,As_order)); % sorted Amatrix (lag matrix)
Assorted = As_across_bins(:,As_order);
binvec = binvec(:,As_order);

As_el = cell(1,size(Amatrix,2));

%collect assembly elements
for i = 1:nAs
    As_el{i} = find(sAmat(:,i));
end

%% now we loop through As_el to find any intersection
%%    between pairs of assemblies
ComprAs = [];
elcell = [];
counter = 0;
for i = 1:nAs-1
    for j = i+1:nAs
        a = intersect(As_el{i},As_el{j});
        if numel(a)>=2
            c = [];
            counter = counter + 1;
            c.elements = a;
            c.AsId = [i j];
            elcell{1,counter} = a;
            ComprAs{1,counter} = c;
        end
    end
end

%% now prune collected compressed Assembly
% Auidx is labeled vector for assemblies in ComprAs
[Au,Auidx] = uniquecell(elcell); %get unique groups based on elements
PrComprAs = cell(1, numel(Au));% pruned ComprAs vector

for j = 1:numel(Au)
    c = [];
    c.elements = Au{j};
    c.Act = cell(1,numel(c.elements));
    c.binlag = [];
    for i = 1:numel(Auidx)
        if Auidx(i) ==j %collect lag and bin to reconstruct the activity
      
            A1 = Assorted{1,ComprAs{i}.AsId(1)};
            [~,eleid] = ismember(c.elements,A1.elements);
            templag = A1.lag(eleid);
            templag = templag - min(templag);
            tempbinlag = [A1.bin,templag];
            c.binlag = [c.binlag {tempbinlag} ];
      
            A2 = Assorted{1,ComprAs{i}.AsId(2)};
            [~,eleid] = ismember(c.elements,A2.elements);
            templag = A2.lag(eleid);
            templag = templag - min(templag);
            tempbinlag = [A2.bin,templag];
            c.binlag = [c.binlag {tempbinlag}];
        end
    end
    %% collect activity
    if 1
        c.binlag = uniquecell(c.binlag);
        %now we collect activity across binlag
        for ii = 1:numel(c.binlag)
            binw = c.binlag{ii}(1);
            lag = c.binlag{ii}(2:end);
            S = spM(c.elements,:);
            minspM = min(spM(:));
            maxspM = max(spM(:));
            tempact = collectsubSpM(S,lag,binw,[minspM maxspM]);
            for jj = 1:numel(c.Act)
                c.Act{jj} = sort(unique([c.Act{jj} tempact{jj}]),'ascend');
            end
        end
        %convert c.Act to matrix
        tempmat = nan(numel(c.Act),max(cellfun(@(c) numel(c),c.Act)));
        for jj = 1:numel(c.Act)
            tempmat(jj,1:numel(c.Act{jj})) = c.Act{jj};
        end
        c.Act = tempmat;
    end
    %%
    PrComprAs{j} = c;
  
end