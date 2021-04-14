%preprocessing_ranking_bootstrap - Printing Table 2 in the supplementary material.  Note that 
% the figures and tables will not reproduce the exact results of the paper
% due to the way the simulated data (metabolomics_data.mat, clinical_data.mat) 
% was produced. The original data can not be included publically in this code due to 
% privacy reasons. However, the attached simulated data can be used as a
% means to know how the workflow is constructed and should be used.
%
% Syntax:  preprocessing_ranking_bootstrap()
%
% Inputs:
%    
%
% Outputs:
%    table and plot of matches between clinical data and clusters found in the meta-
% -bolomics data with different pre-processing methods. The resulting table 
% will be very different from the one displayed in the paper due to the procedure of
% making the simulated data used to run this script.
%
% Example: 
%    preprocessing_ranking_bootstrap()
%
% Other m-files required: 
% Subfunctions: 
% MAT-files required:
% CSV-files required: 
%
% See also: 
% Author: Viktor Skantze
% email: viktor.skantze@fcc.chalmers.se
% August 2020
%------------- BEGIN CODE --------------
clc
clear all
%% Load data
load metabolomics_data.mat
load metabolite_list.mat
load variable_names.mat
load clinical_data.mat
%% Permute
tensor_n= permute(tensor_n, [4 3 2 1]);                                    % rearranging order
%% Subtracting baseline
data_subtracted = tensor_n - tensor_n(:,:,1,:);                                            
%% Removing first index
data_subtracted(:,:,1,:) = [];
data_subtracted = permute(data_subtracted, [1 4 3  2]);                    % doing the same as in main.m
data_subtracted = SmoothingTimeSeries(data_subtracted,0);
data_subtracted = permute(data_subtracted,[1 4 3 2]);
%% Included data
data_included = tensor_n;
%% Data type
data_type = 2;                                                              % Included = 1, subtracted = 2
%% Tensor settings
Options(2) = 0;
Options(3) = 0;
Options(1) = 1e-6;
Options(6) = 1e6;
Options(5) = NaN;
const = [0 0 0 0];
components_for_decomposition = 2;
%% Combinations of components
list_combinations = {};                                            % Preallocate
iter = 0;
for i = 1:components_for_decomposition                     % Making the list of all combinations from 1 to number_of_components_2_analyse
    tmp = nchoosek(1:components_for_decomposition,i);
    for ii = 1:size(tmp,1)
        iter = iter+1;
        list_combinations{iter} = tmp(ii,:);
    end
end
%% Preallocation
table_preprocess = array2table(zeros(7,9),'VariableNames',var_names);      % clinical parameter names
table_preprocess.Properties.RowNames = {'A','B','C','Da','Db','Ea','Eb'};  % pre-processing names
%% Number of iterations
var_num = 100;
variable_list = cell(var_num,1);                                           % Preallocate
%% Cut-off limit for signifigance - 95%
cut_off_limit_p = 0.05;                                                    % Cut-off limit for the anova p-values for the clinical parameters
silhouette_cut_off = 0.65;                                                 % Cut-off for sihlouette values of clusters
max_clust_nr = 4;                                                          % Maximal number of clusters to try
%% Main loop for normalizations
for i_preprocess = 1:7
    %% Predefine
    res_anova = cell(var_num,2);                                           % Preallocate
    clusters_save = cell(var_num,2);                                       % Preallocate
    stack_all_near_solutions = {};
    data_near_stacked_names = {};    
    %% Main loop
    for var = 1:var_num
        if data_type == 1                                                  % Load time data with all baseline included
            data_2B_processed = data_included;
        elseif data_type == 2                                              % Load time data with all baseline subtracted
            data_2B_processed = data_subtracted;
        end
        %% Variable sampling
        variables = randperm(size(data_2B_processed,2));                   % random sampling of metabolites
        variables = variables(1:randi([10 size(data_2B_processed,2)],1));  % random number of samples between 10 and 79
        metabolite_list_sampled = metabolite_list(variables);              % new list of metabolites
        data_2B_processed = data_2B_processed(:,variables,:,:);            % bootstrapped data without repeats
        %% Normalizations
        if i_preprocess == 1
            data = ( data_2B_processed - mean(data_2B_processed,4)) ./ std(data_2B_processed,0,4); % Normalization A
        elseif i_preprocess == 2
            data = (data_2B_processed - mean(data_2B_processed,[1 3 4])) ./ std(data_2B_processed,0,[1 3 4]); % Standardize the tensor according to normalization B
        elseif i_preprocess == 3
            data = (data_2B_processed - mean(data_2B_processed,[1 4])) ./ std(data_2B_processed,0,[1 4]); % Normalization C
        elseif i_preprocess == 4
            data = (data_2B_processed - mean(data_2B_processed,[1 3])) ./ std(data_2B_processed,0,4); % Normalization Da
        elseif i_preprocess == 5
            data = (data_2B_processed - mean(data_2B_processed,[1 3])) ./ std(data_2B_processed,0,[ 1 3 4]); % Normalization Db
        elseif i_preprocess == 6
            norm_E_tmp = zeros(size(data_2B_processed));                                  % Preallocate
            for m = 1:size(data_2B_processed,2)                                                      % Subtracting the baseline from each curve
                for i = 1:17
                    norm_E_tmp(:,m,:,i) = data_2B_processed(:,m,:,i) - mean( data_2B_processed(:,m,1,i) );
                end
            end
            data = norm_E_tmp ./std(data_2B_processed,0,4);                    % Normalization Ea
        elseif i_preprocess == 7
            norm_E_tmp = zeros(size(data_2B_processed));                                  % Preallocate
            for m = 1:size(data_2B_processed,2)                                                      % Subtracting the baseline from each curve
                for i = 1:17
                    norm_E_tmp(:,m,:,i) = data_2B_processed(:,m,:,i) - mean( data_2B_processed(:,m,1,i) );
                end
            end
            data = norm_E_tmp ./ std(data_2B_processed,0,[1 3 4]);             % Normalization Eb
        end
        %% Tensor decomposition
        [model,~,~,~] = parafac(data,components_for_decomposition, Options, const); % first model
        if i_preprocess == 2
            variable_list{var,1} = variables;
            for i_p = 1:5                                                        % repeating the decomposition for 5 different random initializations
                [factor_matrices,~,error,corcondia,factor_degeneracy] = ParafacTwoFactorDegeneracy(data,components_for_decomposition, Options, const); % tensor decomp
                corecondia_save(i_p, var) = corcondia;                         % saving core consistency
                error_save(i_p,var) = error;                                   % saving error
                [~,similarity_save(i_p,var), ~] = ModelSimilarity(model,factor_matrices); % saving similarity
                factor_degen(i_p,var) = factor_degeneracy;                     % saving factor degenaracy
            end
        end
        individual_component = model{4};                         % The individual mode to be clustered
        %% Cluster
        mean_silho = zeros(size(list_combinations,2), max_clust_nr-1);     % Preallocate
        index_clusters_cell = cell(size(list_combinations,2), max_clust_nr-1);% Preallocate
        loading_mets = cell(size(list_combinations,2), max_clust_nr-1);    % Preallocate
        singleton_cluster = zeros(size(list_combinations,2), max_clust_nr-1);% Preallocate
        clusters = cell(size(list_combinations,2), max_clust_nr-1);        % Preallocate
        
        for i = 1:size(list_combinations,2)                                % Going through all combinations of components and clustering each one of them
            comp_interesting_clust = list_combinations{i};                 % Current components to cluster
            for c = 2:max_clust_nr                                         % Looping over the number of clusters
                cluster_nr = c;
                [clusters{i,c-1},mean_silho(i,c-1),index_cluster_ind, centers] = ClusterComponents(individual_component(:,comp_interesting_clust),cluster_nr,comp_interesting_clust,0);% Cluster the current components with c number of clusters
                index_clusters_cell{i,c-1} =  index_cluster_ind;           % Saving the indices of the clustering
                [occurences,~] = groupcounts(index_cluster_ind);           % If the solution contains a cluster with only one individual then we mark that with a 1;
                if any(occurences == 1)
                    singleton_cluster(i,c-1) = 1;
                else
                    singleton_cluster(i,c-1) = 0;
                end
            end
        end
        %% Evaluation
        table_regression = [];                                             % Preallocate
        table_anova = [];                                                  % Preallocate
        for i = 1:size(list_combinations,2)
            comp_interesting_clust = list_combinations{i};
            for j = 1:max_clust_nr-1
                t_anova = AssociationByANOVA(clinical_data_n, var_names,index_clusters_cell{i,j}, 0);% Evaluate the clusters against clinical parameters
                t_anova.Properties.RowNames = {['Component(s): ',num2str(comp_interesting_clust),' Clusters: ',num2str(j+1)]}; % Write information to table
                t_silho = table(mean_silho(i,j),'VariableName',{'Mean_silhouette'}); % Write information to table
                if singleton_cluster(i,j) == 1                             % If a cluster contains only one individual then the solution is not valid and therefore the p-values are corrected to 1.
                    t_anova{1,:} = 1;
                end
                table_anova = [table_anova; t_anova, t_silho];             % Update and merge tables
            end
        end
        %% Sorting out insignificant solutions
        table_anova{table_anova.Mean_silhouette<silhouette_cut_off,1:9} = 1;% Discarding solutions that don't reach the cut-off value for silhouette
        tmp_matrix = table_anova.Variables;                                % Cutting the solutions that lie above the cut-off limit
        tmp_matrix(tmp_matrix>cut_off_limit_p) = 1;
        table_anova.Variables = tmp_matrix;
        %% Save results
        res_anova{var,data_type} = table_anova;                          % Saving the anova tables for each data type
        clusters_save{var,data_type} = reshape(clusters',numel(clusters),1); % Saving the clusters for each data type
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%% Extract recurring clusters %%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Looping through each clinical parameter for recurring cluster solutions
    for i_clinical = 1:9                                                   % Going through clinical parameters
        stacked_solutions_current = {};                                    % Preallocate
        row_stacked = [];                                                  % Preallocate
        for var = 1:var_num
            table_current = res_anova{var,data_type}.Variables;
            clusters_current = clusters_save{var,data_type};
            table_current = table_current(:,1:9);                          % Only use the clinical parameters, skip the mean silhouette values
            [row_hits, col_hits] = find(table_current(:,i_clinical)<1);    % Finding the valid solutions below the cut-off p-value
            
            for i = 1:length(row_hits)
                stacked_solutions_current{end+1} = clusters_current{row_hits(i)};% Stacking the cluster solutions that are valid
                row_stacked(end+1,:) = [row_hits(i), data_type, var];   % Stacking the indices and the data type
            end
        end
        %% Analyse cluster consensus
        if length(stacked_solutions_current)>1                             % If there are several valid solutions, they are compared for similarity.
            cell_groups_of_equal_solutions = EqualSolutions(stacked_solutions_current);
            %% Display near equal solutions
            for i = 1:length(cell_groups_of_equal_solutions)               % Looping through the stacked equal solutions
                current_near_solutions = cell_groups_of_equal_solutions{i};% Current solution that
                for n = 1:length(current_near_solutions)                   % Looping through the data types where the clusters were found
                    stack_all_near_solutions{end+1} = stacked_solutions_current{current_near_solutions(n)}; % Stack all valid solutions for comparison later
                    data_near_stacked_names{end+1} = var_names{i_clinical};% Stack the clinical names for display later
                end
                
            end
        end
        
    end
    %% Check similarity between clinical params
    cell_groups_of_equal_solutions = EqualSolutions(stack_all_near_solutions);
    for c_i = 1:length(cell_groups_of_equal_solutions)                     % Saving all the repeated cluster solutions and to what clinical parameter they mathced to.
        
        current_same_sol = cell_groups_of_equal_solutions{c_i};
        current_clinical_param = data_near_stacked_names(current_same_sol);
        
        
        index_clinical_parameter = find(contains(table_preprocess.Properties.VariableNames,current_clinical_param));
        for i_index = index_clinical_parameter
            if table_preprocess{i_preprocess,i_index}<length(current_same_sol)
                table_preprocess{i_preprocess,i_index} = length(current_same_sol); % Saving for each pre-processing
            end
        end
        
    end
    
end
%% Displaying results
cell_table = table2cell(table_preprocess);
Xt = cell2table(cell_table','RowNames',table_preprocess.Properties.VariableNames,'VariableNames',table_preprocess.Properties.RowNames) % Displaying table
%% Make negative core consistency to zero
corecondia_save( corecondia_save<0 ) = 0;                                                      
%% Plot the results
factor_degen = sort(factor_degen,2);
for j = 1:var_num
    for i= 1:5
        if factor_degen(i,j) == 0
            subplot(1,2,1)
            scatter(cellfun(@length, variable_list(j)),similarity_save(i,j), 'ko','filled');hold on
            ylim([0 1])
            ylabel('\textbf{Stability}','interpreter','latex')
            xlabel('\textbf{Number of metabolites}','interpreter','latex')
            subplot(1,2,2)
            scatter(cellfun(@length, variable_list(j)),corecondia_save(i,j), 'ko','filled');hold on
            ylabel('\textbf{Core consistency}','interpreter','latex')
            xlabel('\textbf{Number of metabolites}','interpreter','latex')
            ylim([0 100])
        else
            subplot(1,2,1)
            plot(cellfun(@length, variable_list(j,data_analysed)),similarity_save(i,j), 'rx','linewidth',1.5);hold on
            ylim([0 1])
            ylabel('\textbf{Stability}','interpreter','latex')
            xlabel('\textbf{Number of metabolites}','interpreter','latex')
            subplot(1,2,2)
            plot(cellfun(@length, variable_list(j,data_analysed)),corecondia_save(i,j), 'rx','linewidth',1.5);hold on
            ylabel('\textbf{Core consistency}','interpreter','latex')
            xlabel('\textbf{Number of metabolites}','interpreter','latex')
            ylim([0 100])
        end
    end
end

%------------- END CODE --------------