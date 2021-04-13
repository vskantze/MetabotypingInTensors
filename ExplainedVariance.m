function exp_var = ExplainedVariance(model,data)
%ExplainedVariance - check what each component explains.
%
% Syntax:  exp_var = ExplainedVariance(model,data)
%
% Inputs:
%    model - tensor decomp model ,data- metabolomics data
%
% Outputs:
%    exp_var - explained variation per component               
%
% Example: 
%    exp_var = ExplainedVariance(model,data)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required:
% CSV-files required: 
%
% See also: 
% Author: Viktor Skantze
% email: viktor.skantze@fcc.chalmers.se
% August 2020
%------------- BEGIN CODE --------------
number_components = size(model{1},2);                                      % preallocate
exp_var = zeros(1,number_components);                                       
for i = 1:number_components                                                % looping through components
    reconstructed_data = nmodel([{model{1}(:,i)} {model{2}(:,i)} {model{3}(:,i)} {model{4}(:,i)}]); % reconstruct data from model
    exp_var(i)=100 * (1 - (sum((data(:)-reconstructed_data(:)).^2,'omitnan') / sum(data(:).^2,'omitnan'))); % save explained variation
end
end
%------------- END CODE --------------