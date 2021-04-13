function classification_list = ClassificationOfMetabolites(data)
%ClassificationOfMetabolites - Checking if the metabolite mostly has its peak 
% before half of the time span or not. Late metabolites are classified if
% the peak happens after the half of the time span and the inverse case for
% the fast metabolites.
%
% Syntax:  classification_list = ClassificationOfMetabolites(data)
%
% Inputs:
%    data - metabolomics tensor data
%
% Outputs:
%    classification_list - list of indeces indicating if the metabolite is 
% of fast.
%
% Example: 
%    classification_list = ClassificationOfMetabolites(data)
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


data = permute(data, [4 3 2 1]);                                           % d m t i
classification_list = zeros(size(data,2),1);
for i = 1:size(data,2)                                                     % looping through metabolites
    mean_diet = squeeze(mean(squeeze(data(:,i,:,:)),1))';                  % average curves over diets
    if  mean(mean_diet(:,1:4),'all') > mean(mean_diet(:,5:end),'all')      % classifying if slow of fast by averaging over individuals
        classification_list(i) =  1;                                       % fast
    else
        classification_list(i) =  2;                                       % slow
    end
end

end