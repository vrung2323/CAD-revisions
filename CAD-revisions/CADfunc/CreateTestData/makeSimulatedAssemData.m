% Create Poisson process as described in RD17

nC = 10;            % Total # of cells
T  = 600;          % Total time (s)
lam = 3;          % Avg firing rate for background process
fracVar = 0;      % Fraction which is variable, vs. constant rate Poisson

% Smaller frequency allows full assembly to be detected relative to the
% background firing rate. Similarly, a higher background firing rate
% relative to the frequency makes it difficult to detect the full assembly.

%% Assembly parameters
assemb_type = 2; % vector of types of assembly (1-5)
nAsC = 5; % vector corresponding to the number of cells in corresp. assembly
fC = 1; % vector corresponding to the first cells in each assembly
%%

if (fracVar > 0)
    % If we need to generate a slowly varying rate
    S = makeAR(nC,T,lam);
    S = fracVar*S + lam*(1-fracVar);
end

% Estimate size of spike matrix
maxD   = round(lam*T*1.5);
bgSpM  = nan(nC,maxD);

if (fracVar >0)
    % Now create spikes, using a call to 
    % genNHPP:   creates a spike train for a nonhomogeneous Poisson process
    for j1=1:nC
        rvec = S(j1,:); tvec = 0:T;
        ETtemp  = genNHPP(@(x)eval_pw_const(x,rvec,tvec),T,1);

        nEv     = length(ETtemp);
        bgSpM(j1,1:nEv)=ETtemp;   
    end
else
    % Background process is homogeneous. Use an exponential distribution
    % for ISIs   
    for j1=1:nC
        isitemp = (-1/lam)*log(rand(1,maxD));
        sptemp  = cumsum(isitemp);
        
        % Cut off at T
        eind = find(sptemp>T,1);
        bgSpM(j1,1:(eind-1)) = sptemp(1:eind-1);
    end
end
    
% Delete unneccessary columns
bgSpM(:, all(isnan(bgSpM)))=[];

% Add on refractory period as suggested by RD17
tref  = 0.002;
maxSp = size(bgSpM,2);
bgSpM = bgSpM+ones(nC,1)*[0:maxSp-1]*tref;

% Depending on how strongly you feel about T, you might cut off the 
% spike times again here.

if (0)
% Plot
figure;
subplot(2,1,1);hold on;
for j1=1:nC
    aus = bgSpM(j1,:); aus(isnan(aus))=[];
    plot(aus,j1*ones(size(aus)),'.','color',[0.5 0.5 0.5]);
end
subplot(2,1,2);plot(0:T,S');
end

% Make assembly
for i = 1:length(assemb_type)
    switch assemb_type(i)
        case 1
            % TYPE I
            freq = 20; % <= 50 to detect full assembly @ lam = 0.1
            spAs   = genASp_Type1_fn(nAsC(i),freq,T);
        case 2
            % TYPE II
            freq = 3; % <= 15 to detect full assembly @ lam = 0.1
            spAs   = genASp_Type2_fn(nAsC(i),freq,T);
        case 3
            % TYPE III
            % The larger "stretch", the more likely cells will appear in 
            %   the indicated order
            %
            % Possible issue w/ this structure: patterns are up to 200 ms long
            %      high frequency can lead to "run-in" and violations of 
            %      tref. Need to decide what to do...
            freq    = 0.5;
            stretch = 2;
            spAs   = genASp_Type3_fn(nAsC(i),freq,T,stretch);
        case 4
            % TYPE IV
            % The larger "stretch", the more likely cells will appear in 
            %   the indicated order
            freq    = 0.2;
            stretch = 2;
            spAs   = genASp_Type4_fn(nAsC(i),freq,T,stretch); 
        case 5
            % TYPE V
            freq    = 0.2;
            spAs   = genASp_Type5_fn(nAsC(i),freq,T);
        case 10
            % Theta modulated place field spikes
            freq = 0.04;
            pm=[]; 
            pm.s    = 0.99;      % theta modulation parameter
            pm.tref = 0.002;  % refractory period
            pm.sig  = 8;
            spAs   = genASp_LinTrackTheta_fn(nAsC(i),freq,T,pm); 
    end
    
    % Am I doing refrac right?
    min(diff(spAs,[],2),[],2) 


    % Embed into background noise
    %cellIDs     = randperm(nC,nAsC(i)); % Random
    cellIDs     = (1:nAsC(i)) + fC(i) - 1; % Sequential
    spM       = embedAssembly_fn(spAs,bgSpM,cellIDs);
    bgSpM       = spM;
end


%% Save data
% Some possible naming conventions:
%
% If looking at theta procession in place fields
%   fname = sprintf('test_data_s=%4.2f_lam=%4.2f.mat',pm.s,lam);
%
% If doing a "type X" assembly
%   fname = sprintf('test_data_type=%d.mat,assemb_type);
fname = 'my_spike_trains.mat';
save(fname,'spM','spAs','cellIDs');

figure;
%Plot all spikes
subplot(211);hold on;
for j1=1:nC
    aus = spM(j1,:); aus(isnan(aus))=[];
    plot(aus,j1*ones(size(aus)),'.','color',[0.5 0.5 0.5]);
end

% Distinguish assemblies
%subplot(212);hold on;
%for j1=1:nC
    %aus = newSp(j1,:); aus(isnan(aus))=[];
    %plot(aus,j1*ones(size(aus)),'.','color',[0.5 0.5 0.5]);
%end
%for j1=1:length(cellIDs)
    %aus = spAs(j1,:); aus(isnan(aus))=[];
    %plot(aus,cellIDs(j1)*ones(size(aus)),'.','MarkerSize',10);
%end

%% Am I doing refractory period right?
min(diff(spM,[],2),[],2)  % Should be <= tref