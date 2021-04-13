function [mean_similarity,min_similarity, similarity_save] = ModelSimilarity(model_1,model_2)
%ModelSimilarity - Measure model similirty with cosine similarity.
%
% Syntax:  [mean_correlation,min_correlation, correlation_save] = ModelSimilarity(model_1,model_2)
%
% Inputs:
%    model_1,model_2 - two CP tensor decomposition models
%
% Outputs:
%    mean_correlation,min_correlation, correlation_save               
%
% Example: 
%    [mean_similarity,min_similarity, similarity_save] = ModelSimilarity(model_1,model_2)
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
nr_components = size(model_1{1},2);                                        % nr components in the models
nr_modes = size(model_1,2);                                                % nr modes in the models
similarity_save = zeros(nr_modes,nr_components);                           % preallocation

for m = 1:nr_modes                                                         % looping of modes
    mode_current_1 = model_1{m};
    mode_current_2 = model_2{m};                                           % extracting one mode at a time
    for i = 1:nr_components                                                % checking one component at a time
        factor_array_current = mode_current_1(:,i);                        % extracting factor array for comparison
        tmp = zeros(1,nr_components);
        for j = 1:nr_components                                            % checking all other components of the other model becuase models can be permuted
            
            factor_array_compare = mode_current_2(:,j);
            tmp(j) = getCosineSimilarity(factor_array_current, factor_array_compare); % check cosine similarity
        end
        similarity_save(m,i) = max(tmp);                                   % save results
    end
end
min_similarity = min(similarity_save,[],'all');
mean_similarity = mean(similarity_save,'all');
end
%------------- END CODE --------------