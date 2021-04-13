function t_anova = AssociationByANOVA(clinical_data, var_names, index_cluster_ind, plt)
%AssociationByANOVA - Evaluating the clusters found in metabolomics data
%                       with the clinical parameters. ANOVA tests and
%                       regression tests.
%
% Syntax:  [t_reg, t_anova] = AssociationByANOVA(data, var_names, index_cluster_ind, individual_comp, plt)
%
% Inputs:
%       clinical_data - The clinical data to validate potential metabotypes
%       var_names - Variable names for the clinical data.
%       index_cluster_ind - The cluster assignment in the form of an array
%                           [17,1].
%       plt - [0 1] To plot the ANOVA and the regression plots.
%
%
% Outputs:
%       t_reg - Table of regression p-values.
%       t_anova - Table of ANOVA p-values.
%
% Example: 
%    [t_reg, t_anova] = AssociationByANOVA(clinical_data, var_names, index_cluster_ind, individual_comp, 0)
%
% Other m-files required: 
% Subfunctions: none
% MAT-files required: none
% CSV-files required: none
%
% See also: 
% Author: Viktor Skantze
% email: viktor.skantze@fcc.chalmers.se
% August 2020
%------------- BEGIN CODE --------------

%% ANOVA clinical parameters
p_anova_clinic = zeros(1,size(clinical_data,2));                           % Preallocate
for i = 1:size(clinical_data,2)                                            % Looping through the clinical variables for ANOVA test
    if plt
        figure()
    end
    p_anova_clinic(i) = ValidationAnova(index_cluster_ind,clinical_data(:,i),plt,var_names{i}); % Using a one-way ANOVA test to get a p-value on the cluster assignment
end
%% Table

t_anova = array2table(p_anova_clinic,'VariableNames',var_names);           % Making table of the ANOVA p-values


end

function p_anova = ValidationAnova(cluster_assignments,clinic_param,plt,str)
%ValidationAnova - One-way ANOVA test on clinical parameters with the
%                  cluster assignment from the metabotyping algorithm.
%
% Syntax:  p_anova = ValidationAnova(cluster_assignments,clinic_param,plt,str)
%
% Inputs:
%       cluster_assignments - The cluster index every individual is
%                             assigned. Array of [17,1].
%       clinic_param - Clinical parameter to perform one-way ANOVA test on.
%       plt - To plot the result or not [0 1].
%       str - The clinical parameter name.
%
% Outputs:
%       p_anova - The p-value of the ANOVA test.
%
% Example: 
%    p_anova = ValidationAnova(cluster_assignments,clinic_param,plt,str)
%
% Other m-files required: distinguishable_colors.m
% Subfunctions: none
% MAT-files required: none
% CSV-files required: none
%
% See also: 
% Author: Viktor Skantze
% email: viktor.skantze@fcc.chalmers.se
% August 2020
%------------- BEGIN CODE --------------
if plt                                                                     % If plotting is wanted.
    %% Extract index
    colors = distinguishable_colors(4);                                    % Using different colors for different clusters     
    for i = 1:max(cluster_assignments)                                                            % For three clusters
        scatter(zeros(size(find(cluster_assignments==i))), clinic_param(cluster_assignments==i),'MarkerEdgeColor','black', 'MarkerFaceColor',colors(i,:));hold on;
    end
    for k = 1:17
        text(0.1,clinic_param(k),num2str(k),'Color','k')
    end
    title(append('Cluster assignment on ',str))         
end
%% ANOVA
if plt 
    DISPLAYOPT = 'on';                                                     % Plot the ANOVA result or not
else 
    DISPLAYOPT = 'off';
end
p_anova = anova1(clinic_param,cluster_assignments,DISPLAYOPT);             % The one-way ANOVA test on the clinical parameter
end
%------------- END CODE --------------
