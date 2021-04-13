%main - Printing most of the manuscript figures and conclusions. Note that 
% the figures and tables will not reproduce the exact results of the paper
% due to the way the simulated data (metabolomics_data.mat, clinical_data.mat) 
% were produced. The original data can not be included publically in this code due to 
% privacy reasons. However, the attached simulated data can be used as a
% means to know how the workflow is constructed and should be used.
%
% Syntax:  main()
%
% Inputs:
%    
%
% Outputs:
%    plots and prints of the conclusions of the paper.               
%
% Example: 
%    main()
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
clear all
clc
%% Load data
load metabolomics_data.mat
load metabolite_list.mat
load variable_names.mat
load clinical_data.mat
%% Subtracting baseline
data_subtracted = tensor_n - tensor_n(:,1,:,:);
%% Removing first index
data_subtracted(:,1,:,:) = [];
%% Outlier rejection of whole curve
data_subtracted = OutlierRejectionCurve(data_subtracted);                  % Augmenting outliers that span the whole time series
%% Outlier rejection spikes
data_subtracted = OutlierRejectionSpike(data_subtracted);                  % Augmenting outlier spikes in time series
%% Filter data
data_subtracted = SmoothingTimeSeries(data_subtracted,1);                  % Moving average smoothing
%% Normalize per metabolite
data_subtracted_normalized = (data_subtracted - mean(data_subtracted,[1 2 4])) ./ std(data_subtracted,0,[1 2 4]); % Standardize the tensor according to normalization B
%% Extract metabolite labels
classification_list = ClassificationOfMetabolites(data_subtracted_normalized); % Extract the labels of fast and slow metabolites
%% Tensor settings
Options(2) = 0;
Options(3) = 0;
Options(1) = 1e-8;
Options(6) = 1e6;
Options(5) = 10;
const = [0 0 0 0];
%% Plot PARAFAC
compi = 2;                                                                 % Components
reference_indexing = [1 2 1 2 1 1 1 2 1 1 1 1 1 2 2 1 2];                  % Given by kmeans on the CP components
[model,~,err,corcondia] = parafac(data_subtracted_normalized,compi,Options,const); % Tensor decomp
amino_acids = [11,12,13,16,17,20];                                         % Amino acid indeces that are predominant in the metabotype clusters
[~,indeces] = sort(model{3}(:,2),'descend');                               % Sort the metabolites according to loading
ComponentPlotter(compi, model{1},reference_indexing, model{3}, classification_list, model{2}, model{4},indeces(amino_acids)) % plot components
exp_var = ExplainedVariance(model,data_subtracted_normalized);             % Show the explained variance
clc
disp('Explained variation in % component 1 and 2 resp. ')
disp(exp_var)
%% Amino acids
figure()
col = distinguishable_colors(4);                                           % colors
diet_inspect = 2;                                                          % beef
subtracted_data = permute(data_subtracted,[4 3 2 1]);                      % rearranging order
iter = 0;
names = {'\textbf{Taurine}', '\textbf{Phenylalanine}', '\textbf{L-allothreonine}', '\textbf{Threonine}','\textbf{Proline}', '\textbf{Tyrosine}'}; % amino acid names
for i = amino_acids                                                        % plot the amino acids, labeled by metabotype clusters
    iter = iter +1;
    subplot(3,2,iter)
    for c = 1:2
        if c == 1
            plot(1:7,squeeze(subtracted_data(diet_inspect,indeces(i),:,c==reference_indexing)),'--','color',col(c,:),'LineWidth',1.5 );hold on
        else
            plot(1:7,squeeze(subtracted_data(diet_inspect,indeces(i),:,c==reference_indexing)),'color',col(c,:),'LineWidth',1.5 );hold on
        end
    end
    title(names(iter), 'Interpreter', 'latex')
    
    set(gca,'xtick',[])
    if (i == amino_acids(end)) || (i == amino_acids(end-1))
        xlabel('\textbf{Time (h)}', 'Interpreter', 'latex')
        xticks([1:7])
    end
    
    if rem(iter,2) == 1
        ylabel('\textbf{Amplitude}', 'Interpreter', 'latex')
    end
    xlim([1 7])
end
%% Concatenate diet - unfold to 3D tensor
tmp_cat = cat(1,squeeze(subtracted_data(1,:,:,:)), squeeze(subtracted_data(2,:,:,:)), squeeze(subtracted_data(3,:,:,:)));
class_list_cat = [classification_list', classification_list', classification_list'];
%% Concatenate metabolite - unfold to 2D matrix
concat_tensor = squeeze(tmp_cat(1,:,:));
class_list_cat_time = class_list_cat(1)*ones(1,7);
for i = 2:size(tmp_cat,1)
    concat_tensor = cat(1,concat_tensor, squeeze(tmp_cat(i,:,:)));
    class_list_cat_time = [class_list_cat_time, class_list_cat(i)*ones(1,7)];
end
concat_tensor = concat_tensor';
normalized_matrix = normalize(concat_tensor,1);                            % normalize the matrix
%% Plot concat PCA solution
[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED] = pca(normalized_matrix);      % Do PCA on the unfolded normalized tensor
compi = 2;                                                                 % Show two components
ComponentPlotter(compi, SCORE,reference_indexing, COEFF, class_list_cat_time, [], [])   % Plot results
%% Validation
table_anova = AssociationByANOVA(clinical_data_n, var_names, reference_indexing, 0); % Print association with ANOVA to the clusters
fprintf('ANOVA association to clinical data from clusters found in metabolomics data using CP \n')
disp(table_anova)
%------------- END CODE --------------