function [clusters ,mean_silhouette, indices_cluster, centers_clusters] = ClusterComponents(individual_mode_components, number_cluster,components_label,plt)
%ClusterComponents - Cluster individual mode components.
%
% Syntax:  [clusters,mean_silhouette,indices_cluster, centers_clusters] = ClusterComponents(individual_mode_components, number_cluster,components_label,plt)
%
% Inputs:
%    individual_mode_components - Array containing the components from the
%                                 individual mode, used for clustering.
%    number_cluster - The number of clusters that will be used for the 
%                     kmeans clustering.
%    components_label - The labels of the components used, for plotting 
%                       porpuses.
%    plt - If the components clustering should be clustered or not [0 1].
%
% Outputs:
%    clusters - Cell with the cluster from the kmeans algorithm.
%    mean_silhouette - The mean silhouette value from the cluster solution.
%    indices_cluster - The indices of the clusters.
%    centers_clusters - The centers of the clusters in components space.
%
% Example: 
%    [clusters ,mean_silhouette, indices_cluster, centers_clusters] = ClusterComponents(individual_mode_components, 3,[1 2],0)
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

%% Clustering k-means
data_kmeans = individual_mode_components;                                  % Data to be clustered. 
[indices_cluster, centers_clusters] = kmeans(data_kmeans, number_cluster,'Distance','sqeuclidean','Display','off','Replicates',10); % Clustering with kmeans
silho_value = silhouette(data_kmeans,indices_cluster,'sqeuclidean');       % Silhouette values for all individual points
mean_silhouette = mean(silho_value);                                       % The mean silhouette value for the clustering solution. 


clusters = cell(1,number_cluster);                                         % Preallocate
for i = 1:size(data_kmeans,1)                                              % Creating cells of clusters with indices.
    clusters{indices_cluster(i)} = [clusters{indices_cluster(i)},i];
end

if plt                                                                     % Plot the component clustering if wanted
    if size(data_kmeans,2) == 2
        figure();
        col = {'r','b','k','g','y','m','c'};
        for i = 1:number_cluster
            scatter(data_kmeans(indices_cluster==i,1), data_kmeans(indices_cluster==i,2),'MarkerEdgeColor','black', 'MarkerFaceColor',col{i}); hold on;
        end
        
        delta = 0.01;
        for k=1:size(data_kmeans,1)
            text(data_kmeans(k,1),data_kmeans(k,2) + delta,num2str(k),'Color','k');
        end
        
        max_value = max(data_kmeans(:));
        min_value = min(data_kmeans(:));
        
        line([(min_value+max_value)/2,(min_value+max_value)/2],[-1,1],'Color', 'black', 'LineStyle', '--'); hold on;
        line([min_value,max_value],[0,0],'Color', 'black', 'LineStyle', '--'); hold off;
        
        list = 1:number_cluster;
        legend(strcat('Cluster ',num2str(list')))
        xlabel(['Component ',num2str(components_label(1))],'interpreter','latex')
        ylabel(['Component ',num2str(components_label(2))],'interpreter','latex')
        title ('Cluster Assignments','interpreter','latex')
        hold off
        
    elseif size(data_kmeans,2) == 1
        figure;
        col = {'r','b','k','g','y','m','c'};
        for i = 1:number_cluster
            cluster_length = length(data_kmeans(indices_cluster==i));
            scatter(data_kmeans(indices_cluster==i),zeros(1,cluster_length),'MarkerEdgeColor','black', 'MarkerFaceColor',col{i}); hold on;
        end
        indi_1 = data_kmeans(:);
        delta = 0.1;
        for k=1:length(indi_1)
            text(indi_1(k),0 + delta,num2str(k),'Color','k');
        end
        
        max_value = max(data_kmeans);
        min_value = -max_value;
        
        line([(min_value+max_value)/2,(min_value+max_value)/2],[-1,1],'Color', 'black', 'LineStyle', '--'); hold on;
        line([min_value,max_value],[0,0],'Color', 'black', 'LineStyle', '--'); hold off;
        
        list = 1:number_cluster;
        legend(strcat('Cluster ',num2str(list')))
        
         title('Cluster Assignments','interpreter','latex')
        hold off
        
    elseif size(data_kmeans,2) == 3
        figure();
        col = {'r','b','k','g','y','m','c'};
        for i = 1:number_cluster
            scatter3(data_kmeans(indices_cluster==i,1), data_kmeans(indices_cluster==i,2), data_kmeans(indices_cluster==i,3), 'MarkerEdgeColor','black', 'MarkerFaceColor',col{i}); hold on;
        end
        

        
        
        list = 1:number_cluster;
        legend(strcat('Cluster ',num2str(list')))
        xlabel(['Component ',num2str(components_label(1))],'interpreter','latex')
        xlabel(['Component ',num2str(components_label(2))],'interpreter','latex')
        xlabel(['Component ',num2str(components_label(3))],'interpreter','latex')
         title ('Cluster Assignments','interpreter','latex')
        hold off
    end
end

end

%------------- END CODE --------------