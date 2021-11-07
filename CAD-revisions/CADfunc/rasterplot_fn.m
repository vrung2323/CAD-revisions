function rasterplot_fn(SpM,varargin)
%%
% The varargin variable is used for plotting pattern activity. It should be
% in form of a pair of string cell and an array cell
%

nneu = size(SpM,1); %number of neuron

if nargin==1
    for j = 1:nneu
        x = SpM(j,:);
        y = j*ones(1,size(SpM(j,:),2));
        scatter(x,y,20,'MarkerFaceColor','b',...
            'MarkerEdgeColor','none','MarkerFaceAlpha',1);hold on;
    end
else
    % need 1 extra variables to plot assembly activity
    var1 = varargin{1};
    %SpkTime
    if strcmp('spk',var1{1})
        SpT = var1{2};
        if size(SpT,1)~=nneu
            error('The dimension does not match');
        end
    
    % BinTime
    elseif strcmp('bin',var1{1})
        tb = var1{2};
        Time = var1{3};
        lag = var1{4};
        temp = cell(nneu,1);
        %compute SpT
        for j = 1:length(Time)
            if Time(j) ~=0
                for jj = 1:nneu
                    curM = SpM(jj,tb(j+lag(jj))<SpM(jj,:) & SpM(jj,:)<tb(j+lag(jj)+1));
                    temp{jj} = [temp{jj} curM];
                end
            end
        end
        SpT = nan(nneu,max(cellfun(@(c) numel(c),temp)));
        for jj = 1:nneu
            SpT(jj,1:numel(temp{jj})) = temp{jj};
        end
    else
        error('Option is either SpikeTime or BinTime');
    end
  
    %plot
    for j = 1:nneu
        x = SpM(j,:);
        t = SpT(j,:);
        y = j*ones(1,size(x,2));
        yt = j*ones(1,size(t,2));
        scatter(x,y,20,'MarkerFaceColor',[0.5 0.5 0.5],...
            'MarkerEdgeColor','none','MarkerFaceAlpha',0.2);hold on;
        scatter(t,yt,10,'MarkerFaceColor','r',...
            'MarkerEdgeColor','none','MarkerFaceAlpha',1);hold on;
        ylim([0 nneu+1]);
        xlabel('time');ylabel('units')
        set(gca,'FontSize',15,'FontName','Time New Roman');
    end
end