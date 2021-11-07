function [assembly_output]=FindAssemblies_recursive(binM,maxlag,Exp_Thres,...
    alph,gg,Dc,No_th,O_th,bytelimit,reference_lag,foldname)
% The function agglomerate pairs of units (or a unit and a preexisting
% assembly), tests their significance and stop when the detected assemblies
% reach their maximal dimention.
%
%
%
%  © 2016 Russo, Durstewitz.
%  for information please contact eleonora.russo@zi-mannheim.de; daniel.durstewitz@zi-mannheim.de.
%
%  last update 11/01/2016
%% zero order
nu=size(binM,1);
for w1=1:nu
    assembly_in.n{w1}.elements=w1;
    assembly_in.n{w1}.lag=[];
    assembly_in.n{w1}.pr=[];
    assembly_in.n{w1}.Time=binM(w1,:)';
    assembly_in.n{w1}.Noccurrences=sum(binM(w1,:));
end

ANin=ones(nu,nu);
% fprintf('Order %d\n',0);

%% first order: test over pairs
clear assembly_out ANfo

alpha=alph*2/(nu*(nu-1)*(2*maxlag+1));      %% bonferroni correction
n_as=size(ANin,1);
nu=size(ANin,2);

prANout=ones(nu,nu);           % max possible dimension, then I cut


ANfo=zeros(nu,nu);
assembly_out.n=cell(1,n_as*nu);

nns=1;
for w1=1:(nu-1)
    for w2=w1+1:nu
        spikeTrain2=binM(w2,:)';
        n2=w2;
        isSignificant=0;
    
        [assemD]=TestPair(assembly_in.n{w1},spikeTrain2,n2,...
            maxlag,Dc,reference_lag,Exp_Thres);
    
        if (assemD.pr(end)<alpha &&  assemD.Noccurrences(end)>No_th)
            assembly_out.n{nns}=assemD;
            prANout(w1,w2)=assemD.pr;
            ANfo(w1,w2)=1;
            isSignificant=1;
        end
    
        if isSignificant
            nns=nns+1;
        end
        clear spikeTrain2
    end
  
end
assembly_out.n(nns:end)=[];

clear assemS1 assembly_in assemD

% ANfo(w1,w2) = 1 if w2>w1 and (w1,w2) are correlated
aus=ANfo+ANfo';                 % I make ANfo symmetric
aus(aus==2)=1;                  % This should never happen    
ANfo=aus;


assembly=assembly_out;
clear  assembly_out
if isempty(assembly.n)
    assembly_output=[];
end


% Assembly_O1_... contains significant pairs
fname = sprintf('%s/Assembly_O%d_b%d.mat',foldname,1,gg);
% save(fname,'-v6','assembly');
save(fname,'-v7.3','assembly');




%%  second order and more: increase the assembly size by adding a new unit

O=1;
Oincrement=1;
%     fname = sprintf('Assembly_O%d.mat', O);
%     assembly=load(fname);
while Oincrement  && O<(O_th-1)   
  % this loop stops when no more units can be added to the existing
  % assemblies or if we force a maximal assembly size (O_th:= max O_th+1
  % elements in the assembly )
  
  %     fname = sprintf('Assembly_O%d.mat', O);
  %     assembly=load(fname);
  
    Oincrement=0;
    n_as=length(assembly.n);           % number of groups previously found
    assembly_out.n=cell(1,n_as*nu);    % max possible dimension, then I cut
  
    nns=1;
    for w1=1:n_as                         % it runs over the existing assemblies
        w1_elements=assembly.n{w1}.elements;
        [~, w2_to_test]=find(ANfo(w1_elements,:)==1);   
        % I try to add only neurons that have significant first order cooccurrences with members of the assembly
        w2_to_test(ismember(w2_to_test,w1_elements))=[];  
        % I erase the one that are already in the assembly
        w2_to_test=unique(w2_to_test);
    
        alpha=alph/(length(w2_to_test)*n_as*(2*maxlag+1));  %% bonferroni correction only for the test that I actually perform
    
        for ww2=1:length(w2_to_test)
      
            w2=w2_to_test(ww2);
            spikeTrain2=binM(w2,:)';
            %             n2=w2;
            true=0;
      
            [assemD]=TestPair(assembly.n{w1},spikeTrain2,w2,maxlag,Dc,...
                reference_lag,Exp_Thres);
      
            if assemD.pr(end)<alpha && assemD.Noccurrences(end)>No_th
                assembly_out.n{nns}=assemD;
                % (AKB) This cannot actually be what is intended: 
                % Earlier line makes clear ANfo is nu x nu i.e. only cell
                %  pair significances are tested. 
                %
                % (ORIGINAL LINE)
                % ANfo(w1,w2)=1;
                %
                % (AKB) Instead, the following lines now flag a "pairwise" association
                % between w2 and each element in w1. At worst, this will
                % increase the number of tests that need to be performed in
                % the future. 
                ANfo(w1_elements,w2)=1;
                ANfo(w2,w1_elements)=1;
               
                true=1;
                Oincrement=1;
            end
      
            if true
                nns=nns+1;
            end
      
            clear spikeTrain2 assemD
        end
    end
    assembly_out.n(nns:end)=[];
    
    % If nns>1, we have assemblies to write out.
    %       Oincrement = 1 and so the outer loop will continue
    % If nns = 1, NO assembly was found in this round
    %      Also, Oincrement = 0: the outer loop will terminate
    if nns>1
        O=O+1;
        assembly=assembly_out;
        clear assembly_out
    
        %%% pruning step 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % between two assemblies with the same unit set arranged into different
        % configurations I choose the most significant one
    
        na=length(assembly.n);                      % number assemblies
        nelement=O+1;    % number elements for assembly
        selection=nan(na,nelement+1+1);
        assembly_final.n=cell(1,na);     %max possible dimension
        nns=1;
        for i=1:na
            elem=sort(assembly.n{i}.elements);
            [ism,indx]=ismember(elem,selection(:,1:nelement),'rows');
            if ~ism
                assembly_final.n{nns}=assembly.n{i};
                selection(nns,1:nelement)=elem;
                selection(nns,nelement+1)=assembly.n{i}.pr(end);
                selection(nns,nelement+2)=i;
                nns=nns+1;
            else
                % There is already an assembly with the same elements.
                % Keep only themost significant. Discard the other
                if selection(indx,nelement+1)>assembly.n{i}.pr(end)
                    assembly_final.n{indx}=assembly.n{i};
                    selection(indx,nelement+1)=assembly.n{i}.pr(end);
                    selection(indx,nelement+2)=i;
                end
            end
        end
        assembly_final.n(nns:end)=[];
        assembly=assembly_final;
        clear assembly_final
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
    end
    clear  assemS2  assemS1
  
  
    % Q: If nns = 1, then we did not update O->O+1
    %    Won't this overwrite the previous _O%d file, with an 
    %    empty matrix?
    % A: No! If I didn't execute the (nns>1) loop,
    %    assembly is unchanged. So at worst I overwrote the file with the
    %    SAME info.
    fname = sprintf('%s/Assembly_O%d_b%d.mat',foldname, O,gg);
    save(fname,'-v6','assembly');
    % fname = sprintf('Assembly_2O%d.mat', O);
    % save(fname,'-v6','assembly');
  
    bytesize=whos('assembly');
    if bytesize.bytes>bytelimit
        fprintf('The algorithm has been interrupted because assembly structures reached a global size of %d bytes, this limit can be changed in size or removed with the ''bytelimit'' option\n',bytelimit);
        O=O_th;
    end
end

maxOrder=O;

%% % pruning step 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% I remove assemblies whom elements are already ALL included in a bigger
%%% assembly

nns=1;
for o=1:maxOrder
  
    fname = sprintf('%s/Assembly_O%d_b%d.mat',foldname,o,gg);
    load(fname);
    minor=assembly;
    clear assembly
    %     minor=assembly.order{o};
    no=length(minor.n);                      % number assemblies
    selection_o=ones(no,1);
  
    for O=maxOrder:-1:o+1
        fname = sprintf('%s/Assembly_O%d_b%d.mat',foldname, O,gg);
        load(fname);
        major=assembly;
        clear assembly
        %          major=assembly.order{O};
        nO=length(major.n);                      % number of (larger)assemblies
    
        index_elemo=find(selection_o==1);        % Which smaller assemblies to check?
        for i=1:sum(selection_o)
            elemo=minor.n{index_elemo(i)}.elements;
      
            for j=1:nO
                elemO=major.n{j}.elements;
                if ismember(elemo, elemO)
                    % Is the minor assembly a proper subset of the major
                    % assembly? If so,"zero" out the smaller one. 
                    % We will no longer compare it, because we are pruning
                    % it.
                    selection_o(index_elemo(i),end)=0;
                    
                    % We can set j=n0, because there is no longer any need
                    % to continue this loop. We only need to run this until
                    % we find a single major assembly which encloses the
                    % smaller one.
                    j=nO;
                end
            end
        end
    
        % This returns "0" iff ALL entries in selection_o=0
        % This would mean ALL minor assemblies of this size have been
        % removed.
        % We can terminate this loop (for O=maxOrder:-1:o+1), 
        % and move on to the next iteration of the outer loop (for
        % o=1:maxOrder)
        if ~selection_o
            O=0;
        end
    
    end
  
    %%% give as output just the maximal groups at order "o"
    % I am now finished with this order, do not need to consider further
    % At the end of the routine assembly_output will contain all
    %   (non-pruned) assemblies, in order of increasing size
    %    
    index_elemo=find(selection_o==1);
    for i=1:sum(selection_o)
        assembly_output.n{nns}=minor.n{index_elemo(i)};
        nns=nns+1;
    end
   
    % Ensures permanent deletion
    recycle('off')
    fname = sprintf('%s/Assembly_O%d_b%d.mat',foldname,o,gg);
    delete(fname)
    recycle('on')
    
end

end










