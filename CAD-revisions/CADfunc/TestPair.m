function [assemD] = TestPair(ensemble,spikeTrain2,n2,maxlag,Dc,...
  reference_lag,Exp_Thres)

% this function tests if the two spike trains have repetitive patterns occurring more
% frequently than chance.

% ensemble      := structure with the previously formed assembly and its spike train
% spikeTrain2   := spike train of the new unit to be tested for significance (candidate to be a new assembly member)
% n2            := new unit tested
% maxlag        := maximum lag to be tested
% Dc            := size (in bins) of chunks in which I divide the spike trains to compute the variance (to reduce non stationarity effects on variance estimation)
% reference_lag := lag of reference; if zero or negative reference lag=-l
% Exp_Thres     := Threshold on E[AB]
%
%
%  ?? 2016 Russo, Durstewitz.
%  for information please contact eleonora.russo@zi-mannheim.de; daniel.durstewitz@zi-mannheim.de.
%
%  last update 11/01/2016
%  last update 13/02/2016 added correction for continuity approximation (ExpAB limit and Yates correction)
%
%  2020 Truong, Barreiro 
%  

couple=[ensemble.Time-min(ensemble.Time(:)) spikeTrain2-min(spikeTrain2(:))]; % spike train pair I am going to test

nu=2;
ntp=size(couple,1);  %%%% trial length

% I divide in parallel trials with 0/1 elements
maxrate=min(max(couple)); % it should be min(max()) (p.9)
%maxrate is M in the paper

Zaa=cell(1,maxrate); %use this to store data from subprocesses

if maxrate == 0 
    ExpABi = 0 ; %EpxABi is the expAB for each i
else
    ExpABi = zeros(maxrate,1);
    for i=1:maxrate %alpha
        Zaa{i}=zeros(size(couple), 'uint8');
        Zaa{i}(couple>=i)=1;
        ExpABi(i)=prod(sum(Zaa{i},1))/size(couple,1); 
        % ExpABi does not depend on optimal lag
    end
end


% % decide which is the lag with most coincidences (l_:=best lag)
ctAB=nan(1,maxlag+1);
ctAB_=nan(1,maxlag+1);

for l=0:maxlag
    %couple for positive lag (trAB) and negative lag (trBA)
    trAB=[couple(1:end-maxlag,1), couple(l+1:end-maxlag+l,2)];
    trBA=[couple(l+1:end-maxlag+l,1), couple(1:end-maxlag,2)];
 
    %compute joint spike count for positive lag (ctAB) and negative lag
    %(ctAB_)
    ctAB(l+1)=sum(min(trAB,[],2));
    ctAB_(l+1)=sum(min(trBA,[],2));
end

% % now choosing the optimal lag
if reference_lag<=0 % Use "-l" as a comparison with l
    aus=[ctAB; ctAB_];
    [a,b]=max(aus(:)); %a is maximum value and b is location in 1-D
    [I,J] = ind2sub(size(aus),b); %convert b to 2-D location
    % I = row index (whether optimal lag is + or -)
    % J = magnitude of lag. J=1=> lag = 0, J=2=> lag=+1 or -1, etc.
    
    l_=(I==1)*(J-1)-(I==2)*(J-1);  %% identify optimal l
    % If I==1, then lag = J-1 
    % If I==2, we should pick lag = -(J-1)
else %for reference_lag>0
    Hab_l=[ctAB_(end:-1:2), ctAB];
    [a,b]=max(Hab_l(:));
    lags=-maxlag:maxlag;
    l_=lags(b);
    % Find lag with the max count
    %  Now make sure you keep the value at the desired
    %  reference lag, which is given by an offset 
    Hab=a;
    if l_<0
        l_ref=l_+reference_lag;
        Hab_ref=Hab_l(lags==l_ref);
    else
        l_ref=l_-reference_lag;
        Hab_ref=Hab_l(lags==l_ref);
    end
end



%% HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

%% 
ExpAB=sum(ExpABi);

%% 

% if  a==0  || ExpAB<=5  || ExpAB>=(min(sum(couple))-5)     
%case of no coincidences or limit for the F asymptotical distribution (too few coincidences)
if a == 0 || ExpAB<Exp_Thres
    assemD.elements=[ensemble.elements n2];
    assemD.lag=[ensemble.lag, 99];
    assemD.pr=[ensemble.pr 1];  % setting pr=1 the tested pair will be discarded as assembly
    assemD.Time=[];
    assemD.Noccurrences=[ensemble.Noccurrences 0];
  
else
  
    % % % construct the activation time series for the couple
    len=size(couple,1);        %%%% trial length
    Time=uint8(zeros(len,1));  %%%% activation vector
    if reference_lag<=0
        if l_==0 %zero lag
            for i=1:maxrate
                sA=Zaa{i}(:,1);
                sB=Zaa{i}(:,2);
                Time(1:len)=Time(1:len)+sA(1:end).*sB(1:end);
                %so "Time" variable computes number of activations in each timebin
            end
      
            TPrMTot=[0 ctAB(1); ctAB_(3) 0]; % matrix with #AB and #BA
            %this says that if the l_=0, then use the reference_lag = -2
        elseif l_>0 %positive lag
            for i=1:maxrate
                sA=Zaa{i}(:,1);
                sB=Zaa{i}(:,2);
                Time(1:len-l_)=Time(1:len-l_)+sA(1:end-l_).*sB(l_+1:end);
            end
            TPrMTot=[0 ctAB(J); ctAB_(J) 0]; % matrix with #AB and #BA
        else % negative lag
            for i=1:maxrate
                sA=Zaa{i}(:,1);
                sB=Zaa{i}(:,2);
                Time(-l_+1:end)=Time(-l_+1:end)+sA(-l_+1:end).*sB(1:end+l_);
            end
            TPrMTot=[0 ctAB(J); ctAB_(J) 0]; % matrix with #AB and #BA
        end
    
    else
        %%%%%%%%%%%% only for reference_lag > 0 %%%%%%%%%%%%%%%%%%%
        if l_==0
            for i=1:maxrate
                sA=Zaa{i}(:,1);
                sB=Zaa{i}(:,2);
                Time(1:len)=Time(1:len)+sA(1:end).*sB(1:end);
            end
        elseif l_>0
            for i=1:maxrate
                sA=Zaa{i}(:,1);
                sB=Zaa{i}(:,2);
                Time(1:len-l_)=Time(1:len-l_)+sA(1:end-l_).*sB(l_+1:end);
            end
      
        else
            for i=1:maxrate
                sA=Zaa{i}(:,1);
                sB=Zaa{i}(:,2);
                Time(-l_+1:end)=Time(-l_+1:end)+sA(-l_+1:end).*sB(1:end+l_);
            end
        end
        TPrMTot=[0 Hab; Hab_ref 0]; % matrix with #AB and #BA
    end
  
    %% --------------------------------------------------------------------%
    % % % HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
    % I cut the spike train in stationary segments
    %%%%%%%%%%%%%%%%%%%%% chunking  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nch=ceil((size(couple,1)-maxlag)/Dc);
    Dc=floor((size(couple,1)-maxlag)/nch); %% new chunk size, this is to have all chunks of rougly the same size
    chunked=cell(nch,1);
  
    % in order to take into account the same time series part that I used
    % for MargPr I cut the time series according to l_:
    couple_cut=nan((size(couple,1)-maxlag),2);
  
    if l_==0
        couple_cut=couple(1:end-maxlag,:);
    elseif l_>0
        couple_cut(:,1)=couple(1:end-maxlag,1);
        couple_cut(:,2)=couple(l_+1:end-maxlag+l_,2);
    else
        couple_cut(:,1)=couple(1-l_:end-maxlag-l_,1);
        couple_cut(:,2)=couple(1:end-maxlag,2);
    end
  
  
    for iii=1:nch-1
        chunked{iii}=couple_cut((1+Dc*(iii-1)):Dc*iii,:);
    end
    chunked{nch}=couple_cut((1+Dc*(nch-1)):end,:); % last chunk can be of slightly different size
  
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
    MargPr_t=cell(nch,maxrate);   %% number of spikes in each parallel process in each unit
    maxrate_t=nan(nch,1);
    ch_nn=nan(nch,1);
    for iii=1:nch
        couple_t=chunked{iii};
        maxrate_t(iii)=max(max(couple_t));
        ch_nn(iii)=size(chunked{iii},1);
        Zaa_t=cell(1,maxrate_t(iii));
        for i=1:maxrate_t(iii)
            Zaa_t{i}=zeros(size(couple_t), 'uint8');
            Zaa_t{i}(couple_t>=i)=1;
        end
    
        for i=1:maxrate_t(iii)
            sA=Zaa_t{i}(:,1);
            sB=Zaa_t{i}(:,2);
            MargPr_t{iii,i}=[sum(sA); sum(sB)];
        end
    end
  
    %%%%%%%%%%%%%%%%%%%%% chunks  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
    n=ntp-maxlag;
  
    Mx0=cell(nch,maxrate);
    covABAB=cell(1,nch);
    covABBA=cell(1,nch);
    varT=cell(1,nch);
    covX=cell(1,nch);
    varX=cell(1,nch);
    varXtot=zeros(2);
    for iii=1:nch
        maxrate_t=max(chunked{iii}(:));
        ch_n=ch_nn(iii);
        % % % HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
        % % % evaluation of  #AB
        for i=1:maxrate_t
            if  ~isempty(MargPr_t{iii,i})
                Mx0{iii,i}= MargPr_t{iii,i} *ones(1,2);
            end
        end
    
        % % variance & covariance
        varT{iii}=zeros(nu);
        covABAB{iii}=cell(maxrate_t,maxrate_t);
        for i=1:maxrate_t
            Mx0{iii,i}= MargPr_t{iii,i} *ones(1,2);
            covABAB{iii}{i,i}=(Mx0{iii,i}.*Mx0{iii,i}'./ch_n).*(ch_n-Mx0{iii,i}).*(ch_n-Mx0{iii,i}')./(ch_n*(ch_n-1));
            varT{iii}=varT{iii}+covABAB{iii}{i,i};
            for j=(i+1):maxrate_t
                covABAB{iii}{i,j}=2*(Mx0{iii,j}.*Mx0{iii,j}'./ch_n).*(ch_n-Mx0{iii,i}).*(ch_n-Mx0{iii,i}')./(ch_n*(ch_n-1));
                varT{iii}=varT{iii}+covABAB{iii}{i,j};
            end
        end
    
        %HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
        % % evaluation of  X=#AB-#BA
        covX{iii}=zeros(nu);
        covABBA{iii}=cell(maxrate_t,maxrate_t);
    
        for i=1:maxrate_t
            covABBA{iii}{i,i}=(Mx0{iii,i}.*Mx0{iii,i}'./ch_n).*(ch_n-Mx0{iii,i}).*(ch_n-Mx0{iii,i}')./(ch_n*(ch_n-1)^2);
            covX{iii}=covX{iii}+covABBA{iii}{i,i};
            for j=(i+1):maxrate_t
                covABBA{iii}{i,j}=2*(Mx0{iii,j}.*Mx0{iii,j}'./ch_n).*(ch_n-Mx0{iii,i}).*(ch_n-Mx0{iii,i}')./(ch_n*(ch_n-1)^2);
                covX{iii}=covX{iii}+covABBA{iii}{i,j};
            end
        end
    
        varX{iii}=varT{iii}+varT{iii}'-covX{iii}-covX{iii}';
        varXtot=varXtot+varX{iii};
    end
    %HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
    %HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
    %HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
    X=TPrMTot-TPrMTot';
    if abs(X(1,2))>0
        X=abs(TPrMTot-TPrMTot')-0.5; %Yates correction
    end
  
    if varXtot(1,2)==0  % if variance is zero
        prF=1;
    else
        F=X.^2./varXtot;
        d2 = n;
        % Adjust d2 IF ExpAB <= 5.5
        % See Truong et al, 2020 (PhD thesis)
        if (ExpAB <=5.5)
            d2 = d2/(2*(6-ExpAB));
        end
        prF=fcdf(F(1,2),1,d2,'upper');
        %         prF=fcdf(F(1,2),1,2*double(maxrate)*n-1,'upper');
    end
  
  
  
     %%
    %HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
    %HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
    %All information about the assembly and test are returned
  
    assemD.elements=[ensemble.elements n2];
    assemD.lag=[ensemble.lag, l_];
    assemD.pr=[ensemble.pr prF];
    assemD.Time=Time;
    assemD.Noccurrences=[ensemble.Noccurrences sum(Time)];
    
end
%HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
%HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
%HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
end






