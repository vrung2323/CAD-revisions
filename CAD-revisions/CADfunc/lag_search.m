function output_data = lag_search(saved_data)

%Search Script for finding what the differance in lags are for the CAD bin
%shifting algorithim 
output_data = {}; % initialize cell array of structures

t = linspace(-0.4,0.5,10);

for i = 1:length(saved_data)-1
   for j = i:length(saved_data{i})
       for m = i+1:length(saved_data)
            for n = m:length(saved_data{m})
                if isequal(sort(saved_data{i}{j}.elements),sort(saved_data{m}{n}.elements)) && ...
                        saved_data{i}{j}.bin == saved_data{m}{n}.bin % checks if the binwidth and elements of two assemblies are equal
                        
                    [temp_elements_ij,idx_ij] = sort(saved_data{i}{j}.elements); % sorts elements with indices
                    [~,idx_mn] = sort(saved_data{m}{n}.elements); % sorts elements with indices
                    temp_diff = saved_data{i}{j}.lag(idx_ij) - saved_data{m}{n}.lag(idx_mn); % calculates lag difference
                    
                    if sum(abs(temp_diff)) % only adds to cell array if lag difference is nonzero (relevant)
                        data_struct.elements = temp_elements_ij;
                        data_struct.bin = saved_data{i}{j}.bin;
                        data_struct.shift1 = t(i); data_struct.shift2 = t(m);
                        data_struct.lag_diff = temp_diff;
                        output_data = [output_data; {data_struct}];    
                    end       
                end
            end
       end
   end
end