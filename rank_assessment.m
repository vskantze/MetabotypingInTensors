%rank_assessment - Printing the rank assessment figure from the paper. 
% Note that the figures and tables will not reproduce the exact results of the paper
% due to the way the simulated data (metabolomics_data.mat, clinical_data.mat) 
% was produced. The original data can not be included publically in this code due to 
% privacy reasons. However, the attached simulated data can be used as a
% means to know how the workflow is constructed and should be used.
%
% Syntax:  rank_assessment()
%
% Inputs:
%    
%
% Outputs:
%    plot of the rank assessment methods             
%
% Example: 
%    rank_assessment()
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
%% Loading metabolomics data
load metabolomics_data.mat
%% Subtracting baseline
data_subtracted = tensor_n - tensor_n(:,1,:,:);
%% Removing first index
data_subtracted(:,1,:,:) = [];
%% Outlier rejection of whole curve
data_subtracted = OutlierRejectionCurve(data_subtracted);                  % Augmenting outliers that span the whole time series
%% Outlier rejection spikes
data_subtracted = OutlierRejectionSpike(data_subtracted);                  % Augmenting outlier spikes in time series
%% Filter data
data_subtracted = SmoothingTimeSeries(data_subtracted,1);
%% Normalize per metabolite
data = (data_subtracted - mean(data_subtracted,[1 2 4])) ./ std(data_subtracted,0,[1 2 4]); % Standardize the tensor according to normalization B
%% Tensor settings
Options(2) = 0;
Options(3) = 0;
Options(1) = 1e-8;
Options(6) = 1e6;
Options(5) = 10;
const = [0 0 0 0];
%% Permute
data = permute(data,[4 3  2 1]);
%% Plot decomposition
compi = 2;                                                                 % Choose two components
model= parafac(data,compi,Options,const);                                  % Compute tensor decomposition
ComponentPlotter(compi, model{4},[], model{2}, [], model{3}, model{1},[])  % Plot components
%% Loop for similarity, core consistency and error
components = 10;                                                           % number of components to run
number_repitions = 10;                                                     % reptitions per iteration
similarity_mean = zeros(components,number_repitions);                      % preallocate
similarity_min = zeros(components,number_repitions);                       % preallocate
error = zeros(components,number_repitions);                                % preallocate
corri = zeros(components,number_repitions);                                % preallocate
sparsity = zeros(components,number_repitions);                             % preallocate
factor_degen = zeros(components,number_repitions);                         % preallocate

for i = 1:components                                                       % loop for each number of components
    [model_current,~,~,~] = parafac(data,i,Options,const);                 % tensor decomp
    for j = 1:number_repitions                                             % repitions
        [model_compare,~,err,corcondia,factor_degeneracy] = ParafacTwoFactorDegeneracy(data,i,Options,const);  % tensor decomp
        [mean_correlation,min_correlation] = ModelSimilarity(model_current,model_compare); % similarity
        disp(min_correlation)
        similarity_mean(i,j) = mean_correlation+rand(1)*0.001;             % saving similarity with noise to make visible in plot
        similarity_min(i,j) = min_correlation+rand(1)*0.001;               % saving similarity with noise to make visible in plot
        error(i,j) = err;                                                  % saving error with noise to make visible in plot
        corri(i,j) = corcondia;                                            % saving core consistency with noise to make visible in plot
        factor_degen(i,j) = factor_degeneracy;                             % saving factor degeneracy with noise to make visible in plot
    end 
end
%% Make negative core consistency to zero
corri( corri<0 ) = 0;                                                      
%% Plot the results
factor_degen_sorted = sort(factor_degen,2);
for i = 1:components
    for j = 1:number_repitions
        if factor_degen_sorted(i,j)==0
            subplot(221)
            plot(i,corri(i,j),'ko','markerfacecolor','k'); hold on
            set(gca,'xtick',[])
            ylabel('\textbf{Core consistency}','interpreter','latex')
            ylim([0 100])
            xlim([1 components])
            xticks([1 5 10 15])
            
            subplot(222)
            plot(i,error(i,j),'ko','markerfacecolor','k'); hold on
            set(gca,'xtick',[])
            ylabel('\textbf{Error}','interpreter','latex')
            ylim([0 inf])
            xlim([1 components])
            xticks([1 5 10 15])
            
            
            subplot(2,2,[3 4])
            plot(i,similarity_min(i,j),'ko','markerfacecolor','k');hold on
            ylim([-1 inf])
            xlim([1 components])
            set(gca,'xtick',[])
            xlabel('\textbf{Components R}','interpreter','latex')
            ylabel('\textbf{Stability}','interpreter','latex')
            xticks([1 5 10 15])
        else
            subplot(221)
            plot(i,corri(i,j),'rx','markeredgecolor','r','linewidth',1.5); hold on
            set(gca,'xtick',[])
            ylabel('\textbf{Core consistency}','interpreter','latex')
            ylim([0 100])
            xlim([1 components])
            xticks([1 5 10 15])
            
            subplot(222)
            plot(i,error(i,j),'rx','markeredgecolor','r','linewidth',1.5); hold on
            set(gca,'xtick',[])
            ylabel('\textbf{Error}','interpreter','latex')
            ylim([0 inf])
            xlim([1 components])
            xticks([1 5 10 15])
            
            
            subplot(2,2,[3 4])
            plot(i,similarity_min(i,j),'rx','markeredgecolor','r','linewidth',1.5);hold on
            ylim([-1 inf])
            xlim([1 components])
            set(gca,'xtick',[])
            xlabel('\textbf{Components R}','interpreter','latex')
            ylabel('\textbf{Stability}','interpreter','latex')
            xticks([1 5 10 15])
        end
    end
end
%------------- END CODE --------------