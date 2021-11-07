function AspM = collectsubSpM(spM,lag,binw,mspM)
% min/maxSpM is of the whole population,
% SpM is only of assembly elements (reduce data transfer)
% one element of lag should be zero

minspM = mspM(1);
maxspM = mspM(2);
nneu = size(spM,1);
tb = minspM:binw:maxspM;
binM = zeros(nneu,numel(tb)-1);
maxlag = max(lag);

for n=1:nneu
    [ binM(n,:),~] = histcounts(spM(n,:),tb);
end

binM = binM>0; %convert to binary
%shift with lag
I = 1;
for n = 1: nneu
    I = I.*binM(n,1+lag(n):end-maxlag+lag(n));
end

%collect activity
temp = cell(nneu,1);
for i = 1:numel(I)
    if I(i)==1
        for jj = 1:nneu
            curM = spM(jj,tb(i+lag(jj))<spM(jj,:) & spM(jj,:)<tb(i+lag(jj)+1));
            temp{jj} = [temp{jj} curM];
        end
    end
end
AspM = temp;



