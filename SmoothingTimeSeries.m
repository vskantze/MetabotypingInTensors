function smoothed_data = SmoothingTimeSeries(data,permute_or_not)
%SmoothingData - Smoothing metabolomics data with moving average.
%
%
% Syntax:  smoothed_data = SmoothingData(data)
%
% Inputs:
%    data - Metabolomics tensor data.
%    permute_or_not - permuting the tensor.
%
% Outputs:
%    smoothed_data - The smoothed metabolomics tensor data.                 
%
% Example: 
%    smoothed_data = SmoothingData(data)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
% CSV-files required: none
%
% See also: 
% Author: Viktor Skantze
% email: viktor.skantze@fcc.chalmers.se
% August 2020
%------------- BEGIN CODE --------------


smoothed_data = zeros(size(data));                                         % Preallocate.
for m = 1:size(data,4)                                                     % Go through the entire data by metabolite and diet.
    for d = 1:3
        met_matrix = squeeze(data(d,:,:,m));                               % Collect all data from each metabolite and diet into a matrix.
        smoothed_data(d,:,:,m) = smoothdata(met_matrix,2);                 % Smooth and save data.
    end
end

% for m = 1:size(data,4)                                                     % Go through the entire data by metabolite and diet.
%     for d = 1:3
%         for i = 1:17
%         met_matrix = squeeze(data(d,i,:,m));                               % Collect all data from each metabolite and diet into a matrix.
%         f = fit([1:size(data(d,i,:,m),3)]',met_matrix,'smoothingspline');
%         smoothed_data(d,i,:,m) = f([1:size(data(d,i,:,m),3)]');                 % Smooth and save data.
%         end
%     end
% end

if permute_or_not
    smoothed_data = permute(smoothed_data,[2 3 4 1]);                                            % d,i,t,m
end
end
%------------- END CODE --------------