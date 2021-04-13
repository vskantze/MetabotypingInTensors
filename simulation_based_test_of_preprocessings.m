%simulation_based_test_of_preprocessings - Simulation of metabolomics data to test
%                                 different normalization methods for metabotyping.
% A test is made with augmented groups in the simulated metabolomics data
% and a metabotyping algorithm is supposed to find the simulated groups.
% The test is made n (rec. 150) times and a score is given to each normalization
% method.
%
% Syntax:  simulation_based_test_of_preprocessings.m
%
% Inputs:
%  
%
% Outputs:
%    Table of the scores for the different normalization methods                
%
% Example: 
%    simulation_based_test_of_preprocessings.m
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: 'metabolomics_data.mat'
%
% See also: 
% Author: Viktor Skantze
% email: viktor.skantze@fcc.chalmers.se
% April 2021
%------------- BEGIN CODE --------------
%% 
clc
clear all
%% Load data
load metabolomics_data.mat
%% Permute
tensor_n= permute(tensor_n, [4 3 2 1]);                                    % Permuting order to  diets, metabolites, time, inidividuals
%% Subtracting baseline
data_subtracted = tensor_n - tensor_n(:,:,1,:);                            
%% Removing first index
data_subtracted(:,:,1,:) = [];
data_2B_processed = permute(data_subtracted, [1 4 3  2]);                  % Permuting order to match later operations
%% Noise reduction  
data_2B_processed = SmoothingTimeSeries(data_2B_processed,0);              % Reducing noise with moving average filter
data_2B_processed = permute(data_2B_processed,[1 4 3 2]);                  % Permuting order to match later operations
%% Combinations of components
list_combinations = {};                                                    % Preallocate
iter = 0;
nr_components = 5;                                                     % Number of components to use for clustering
for i = 1:nr_components                                                % Making the list of all combinations from 1 to number_of_components_2_analyse
    tmp = nchoosek(1:nr_components,i);
    for ii = 1:size(tmp,1)
        iter = iter+1;
        list_combinations{iter} = tmp(ii,:);
    end
end
%% Select metabolites
nr_iterations = 150;                                                               
nr_datasets = 5;
score_res = zeros(nr_iterations,nr_datasets,7);                            % Preallocate
nr_metabolites = 30;                                                       % 30 metabolites per dataset                                              
diet_cluster = 3;                                                          % The diet which is augmeneted with clusters
percent_noise = 0.3;                                                       % noise for augmentation
metabolites_cluster = [1 2 12 13 22 23];                                   % Indeces of augmeneted metabolites
silhouette = score_res;                                                    % Preallocate
group_indices = [1 3 7 15];                                                % The deviating group indices
curve_1 = squeeze(data_2B_processed(1,7,:,1));                             % Individual curves to make simulated metabolites from
curve_1 = (curve_1-min(curve_1))/(max(curve_1) -min(curve_1));
curve_2 = squeeze(data_2B_processed(1,2,:,1));                             % Different curves, fast, slower, slow
curve_2 = (curve_2-min(curve_2))/(max(curve_2) -min(curve_2));
curve_3 = squeeze(data_2B_processed(1,4,:,1));
curve_3 = (curve_3-min(curve_3))/(max(curve_3) -min(curve_3));
%% Tensor settings
Options(2) = 0;
Options(3) = 0;
Options(1) = 1e-6;
Options(6) = 1e6;
Options(5) = NaN;
const = [0 0 0 0];
%% Big loop
for i_repeat = 1:nr_iterations                                             % Repeating the test for nr_iterations times
    simulated_data = zeros(nr_metabolites,17,7,3);                         % The simulated data tensor containing the groups
    for i_data = 1:nr_datasets                                             % Running the nr_datasets datasets                    
    scales = randi(nr_iterations,nr_metabolites,1);                        % Scaling the simulated metabolites
        for d = 1:3                                                        % Looping over diets
            for m = 1:nr_metabolites                                       % Looping over metabolites
                noise = rand(17,7);                                        % Individual noise over all time points
                if m <11                                                   % Copying the different curves to make similar metabolites
                    met_own = repmat(curve_1',17,1);
                    met_own_noise = met_own + mean(curve_1)*noise;
                elseif (10<m && m<21)
                    met_own = repmat(curve_2',17,1);
                    met_own_noise = met_own + mean(curve_2)*noise;
                elseif m>20
                    met_own = repmat(curve_3',17,1);
                    met_own_noise = met_own + mean(curve_3)*noise;
                end
                met_own_noise = scales(m).*met_own_noise;                  % Scale metabolites
                
                if any(m== metabolites_cluster) && d == diet_cluster       % Adding the deviating groups in the augmented metabolites and in diet 3
                    if i_data == 1                                         % Add amplitude difference to all time points, alterating with negative deviation
                        if any(m == metabolites_cluster(1:3))              % Adding the deviation to the group in the data
                            met_own_noise(group_indices,:) = met_own_noise(group_indices,:) + std(met_own_noise(group_indices,:),[],'all').*randn(1,7) + mean(met_own_noise(group_indices,:),'all')*percent_noise;
                        else
                            met_own_noise(group_indices,:) = met_own_noise(group_indices,:) - (std(met_own_noise(group_indices,:),[],'all').*randn(1,7) + mean(met_own_noise(group_indices,:),'all')*percent_noise);
                        end
                    elseif i_data == 2                                     % Add amplitude difference to the first 3 time points, alterating with negative deviation
                        if any(m == metabolites_cluster(1:3))              % Adding the deviation to the group in the data
                            met_own_noise(group_indices,1:3) = met_own_noise(group_indices,1:3) + std(met_own_noise(group_indices,1:3),[],'all').*randn(1,3) + mean(met_own_noise(group_indices,1:3),'all')*percent_noise;
                        else
                            met_own_noise(group_indices,1:3) = met_own_noise(group_indices,1:3) - (std(met_own_noise(group_indices,1:3),[],'all').*randn(1,3) + mean(met_own_noise(group_indices,1:3),'all')*percent_noise);
                        end
                    elseif i_data == 3                                     % Add amplitude difference to the last 4 time points, alternating with positive and negative
                        if any(m == metabolites_cluster(1:3))              % Adding the deviation to the group in the data
                            met_own_noise(group_indices,4:end) = met_own_noise(group_indices,4:end) + std(met_own_noise(group_indices,4:end),[],'all').*randn(1,4) + mean(met_own_noise(group_indices,4:end),'all')*percent_noise;
                        else
                            met_own_noise(group_indices,4:end) = met_own_noise(group_indices,4:end) - (std(met_own_noise(group_indices,4:end),[],'all').*randn(1,4) + mean(met_own_noise(group_indices,4:end),'all')*percent_noise);
                        end
                    elseif i_data == 4                                     % Add amplitude difference as bi-modal peaks
                        if any(m == metabolites_cluster(1:3))              % Adding the deviation to the group in the data
                            met_own_noise(group_indices,1:3) = met_own_noise(group_indices,1:3) + std(met_own_noise(group_indices,1:3),[],'all').*randn(1,3) + mean(met_own_noise(group_indices,1:3),'all')*percent_noise;
                            met_own_noise(group_indices,5:end) = met_own_noise(group_indices,5:end) + std(met_own_noise(group_indices,5:end),[],'all').*randn(1,3) + mean(met_own_noise(group_indices,5:end),'all')*percent_noise;
                        else
                            met_own_noise(group_indices,1:3) = met_own_noise(group_indices,1:3) - (std(met_own_noise(group_indices,1:3),[],'all').*randn(1,3) + mean(met_own_noise(group_indices,1:3),'all')*percent_noise);
                            met_own_noise(group_indices,5:end) = met_own_noise(group_indices,5:end) - (std(met_own_noise(group_indices,5:end),[],'all').*randn(1,3) + mean(met_own_noise(group_indices,5:end),'all')*percent_noise);
                        end
                    elseif i_data == 5                                     % Unimodal peak deviations, positive and negative
                        if any(m == metabolites_cluster(1:3))              % Adding the deviation to the group in the data
                            met_own_noise(group_indices,1:3) = met_own_noise(group_indices,1:3) + std(met_own_noise(group_indices,1:3),[],'all').*randn(1,3) + mean(met_own_noise(group_indices,1:3),'all')*percent_noise;
                            met_own_noise(group_indices,5:end) = met_own_noise(group_indices,5:end) - (std(met_own_noise(group_indices,5:end),[],'all').*randn(1,3) + mean(met_own_noise(group_indices,5:end),'all')*percent_noise);
                        else
                            met_own_noise(group_indices,1:3) = met_own_noise(group_indices,1:3) - (std(met_own_noise(group_indices,1:3),[],'all').*randn(1,3) + mean(met_own_noise(group_indices,1:3),'all')*percent_noise);
                            met_own_noise(group_indices,5:end) = met_own_noise(group_indices,5:end) + std(met_own_noise(group_indices,5:end),[],'all').*randn(1,3) + mean(met_own_noise(group_indices,5:end),'all')*percent_noise;
                        end
                    end
                end    
                simulated_data(m,:,:,d) = met_own_noise;                   % Saving the simulated data
            end
        end  
        %% Normalizations / pre-processings
        for i_decomp = 1:7                                                 % Decomposing the different normalizations of the data
            if i_decomp == 1
                norm_A = ( simulated_data - mean(simulated_data,2)) ./ std(simulated_data,0,2); % Normalization A
                tensor_2_decomp = norm_A;                                   
                pre_pro = 'A';
            elseif i_decomp == 2
                norm_B = (simulated_data - mean(simulated_data,[ 2 3 4])) ./ std(simulated_data,0,[2 3 4]); % Normalization B
                tensor_2_decomp = norm_B;
                pre_pro = 'B';
            elseif i_decomp == 3
                norm_C = (simulated_data - mean(simulated_data,[2 4])) ./ std(simulated_data,0,[2 4]); % Normalization C
                tensor_2_decomp = norm_C;
                pre_pro = 'C';
            elseif i_decomp == 4
                norm_Da = (simulated_data - mean(simulated_data,[3 4])) ./ std(simulated_data,0,2); % Normalization Da
                tensor_2_decomp = norm_Da;
                pre_pro = 'Da';
            elseif i_decomp == 5
                norm_Db = (simulated_data - mean(simulated_data,[3 4])) ./ std(simulated_data,0,[ 2 3 4]); % Normalization Db
                tensor_2_decomp = norm_Db;
                pre_pro = 'Db';
            elseif i_decomp == 6
                norm_E_try = zeros(size(norm_A));                          % Preallocate
                for m = 1:nr_metabolites                                   % Subtracting the baseline from each curve
                    for i = 1:17
                        norm_E_try(m,i,:,:) = simulated_data(m,i,:,:) - mean( simulated_data(m,i,1,:) );
                    end
                end
                norm_Ea = norm_E_try ./std(simulated_data,0,2);            % Normalization Ea
                tensor_2_decomp = norm_Ea;
                pre_pro = 'Ea';
            elseif i_decomp == 7
                norm_E_try = zeros(size(norm_A));                          % Preallocate
                for m = 1:nr_metabolites                                   % Subtracting the baseline from each curve
                    for i = 1:17
                        norm_E_try(m,i,:,:) = simulated_data(m,i,:,:) - mean( simulated_data(m,i,1,:) );
                    end
                end
                
                norm_Eb = norm_E_try ./ std(simulated_data,0,[2 3 4]);     % Normalization Eb
                tensor_2_decomp = norm_Eb;
                pre_pro = 'Eb';
            end

            %% Tensor decomposition
            [model,~] = parafac(tensor_2_decomp,nr_components, Options,const); % Decomposing the tensor
            individual_mode = model{2};                           % The individual mode for clustering
            %% Cluster to find the groups
            best_result = 0;
                for i_comb = 1:length(list_combinations)
                    current_components = list_combinations{i_comb};
                    [clusters,mean_silho,index_cluster, centers] = ClusterComponents(individual_mode(:,current_components),2,current_components,0); % Cluster 1 and 2 components because the groups are visible there
                    true_cluster_perm1 = [1 2 1 2 2 2 1 2 2 2 2 2 2 2 1 2 2]';     % True cluster assignment
                    true_cluster_perm2 = [2 1 2 1 1 1 2 1 1 1 1 1 1 1 2 1 1]';     % The inverse cluster assignment because k-means can produce both
                    
                    result1 = sum( index_cluster == true_cluster_perm2 );          % How many correct assignments is the score
                    result2 = sum( index_cluster == true_cluster_perm1 );          % How many correct assignments is the score
                    tmp_result = max([result1,result2]);
                    if tmp_result>best_result
                        best_result = tmp_result;
                        best_silho = mean_silho;
                    end
                end
            score_res(i_repeat, i_data,i_decomp) = best_result;            % The max of the two possible scores
            silhouette(i_repeat, i_data,i_decomp) = best_silho;            % The mean silhouette score for the clustering                
        end
    end
end
%% Result table
var_names = {'A','B', 'C', 'Da', 'Db', 'Ea', 'Eb'};                        % The variable names for the normalizations
runs = {'Data 1', 'Data 2', 'Data 3', 'Data 4','Data 5', 'Mean'};          % The different data sets
score_res_01 = score_res./17;                                              % Normalizing between 0 and 1
score_mean = squeeze(mean(score_res_01,1));                                % Making mean matrix                                       
score_mean(end+1,:) = mean(score_mean,1);                                  % Adding mean to matrix
sill_mean = squeeze(mean(silhouette,1));                                   % Making mean matrix

format shortg
result_table = array2table(score_mean,'VariableNames',var_names,'RowNames',runs);% Making table
result_table.Variables = round(result_table.Variables,3);

disp 'Score per normalization method'
disp(result_table)

sill_table = array2table(sill_mean,'VariableNames',var_names,'RowNames',runs(1:end-1)); % Making table
disp 'Mean silhoette per normalization method'
disp(sill_table)
%------------- END CODE --------------