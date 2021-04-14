function ComponentPlotter(components, individual, cluster_ind, metabolite, clusters_mets, time, diet, var_selected)
%ComponentPlotter - Plots the tensor decomposition model by components and 
% and modes.
%
% Syntax:  ComponentPlotter(components, individual, cluster_ind, metabolite, clusters_mets, time, diet, var_selected)
%
% Inputs:
%    components - nr components
%    individual - individual mode of model
%    cluster_ind - cluster indeces from kmeans of the individual mode
%    metabolite - metabolite mode of model
%    clusters_mets - classified metabolites like fast and slow
%    time - time mode of model
%    diet - diet mode of model
%    var_selected - selected metabolites for visualization as stars in the
%    plot
% Outputs:
%    component plot
%
% Example: 
%    ComponentPlotter(components, individual, cluster_ind, metabolite, clusters_mets, time, diet, var_selected)
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


%% Normalization to fit in the plot window by y = 0.5 -0,5
individual = individual/max(abs(individual),[],'all')/2;
metabolite = metabolite/max(abs(metabolite),[],'all')/2;
time = time/max(abs(time),[],'all')/2;
diet = diet/max(abs(diet),[],'all')/2;


figure()
col = distinguishable_colors(9);                                           % colors
if ~isempty(diet)
    iter = 0;
    for c = 1:components                                                   % looping for each component
        %% Individual mode
        iter = iter + 1;        
        subplot(components,4,iter)
        if ~isempty(cluster_ind)
            for clus = 1:max(cluster_ind)                                  % plotting the individual mode labeled by the cluster indeces
                scatter(find(cluster_ind==clus), individual(cluster_ind==clus,c),'^','filled','MarkerFaceColor',col(clus,:),'MarkerEdgeColor','k');hold on;
            end
        else
            scatter(1:17,individual(:,c),'^','filled','MarkerEdgeColor','k');hold on;
        end
        if c == 1
            title('\bf Individual mode','interpreter','latex')
        end
        ylabel(['\bf Component ', num2str(c)],'interpreter','latex')
        
        ylim([-0.5 0.52])
        xlim([0.5 size(individual,1)])
        xticks([1:4:17])
        if ~(c == components)
            set(gca,'xtick',[])
        else
            xlabel(['\bf Individual index'],'interpreter','latex')
        end
            
        box off
        iter = iter + 1;
        %% Metabolite mode
        subplot(components,4,iter)
        if isempty(clusters_mets)
            scatter(1:size(metabolite,1),metabolite(:,c),'ko','filled');hold on;
        else
            plot(1:size(metabolite,1),metabolite(:,c),'o');hold on;        % plotting the metabolite mode labeled by the metabolite classification 
            for i = 1:max(clusters_mets)
                scatter(find(clusters_mets==i),metabolite(clusters_mets==i,c),'o','filled','MarkerFaceColor',col(i+3,:),'MarkerEdgeColor','k');hold on;
            end
        end
        if ~isempty(var_selected)
            scatter(var_selected,metabolite(var_selected,c),50,'o','MarkerFaceColor',col(3,:),'MarkerEdgeColor','k');hold on;
        end
        
        ylim([-0.5 0.5])
        xlim([1 size(metabolite,1)])
%        xticks([1 20 40 60 size(metabolite,1)])
        if ~(c == components)
            set(gca,'xtick',[])
        else
            xlabel(['\bf Metabolite index'],'interpreter','latex')
        end
        box off
        if c == 1
            title('\bf Metabolite mode','interpreter','latex')
        end
        %% Time mode
        iter = iter + 1;
        subplot(components,4,iter)
        plot(1:size(time,1),time(:,c),'k','LineWidth',1.5);hold on         % plotting the time mode as line to illustrate dynamics
        ylim([-0.5 0.5])
        xlim([0.5 size(time,1)])
        xticks([1:7])
        box off
        if ~(c == components)
            set(gca,'xtick',[])
        else
            xlabel(['\bf Time index'],'interpreter','latex')
        end
        if c == 1
            title('\bf Time mode','interpreter','latex')
        end
        %% Diet mode
        iter = iter + 1;
        set(groot,'defaultAxesTickLabelInterpreter','latex');  
%         if c == 1
%             diet = diet /(2*max(abs(diet),[],'all'));
%         end
        subplot(components,4,iter)                                         % plotting the diet mode and naming the diets
        scatter(1:size(diet,1),diet(:,c),'square','filled','MarkerFaceColor',col(6,:),'MarkerEdgeColor', 'k');hold on
        
        diet_list = {'\textbf{Baked Herring}','\textbf{Beef}','\textbf{Pickled Herring}'};
        displacement_list = [0.1 0;
                             0.1 0.05;
                             -2.7 -0.05];
        for i = 1:size(diet,1)
            message_to_be_printed = diet_list{i};
            %text(i + displacement_list(i,1),diet(i,c) + displacement_list(i,2), message_to_be_printed,'FontSize',12,'Color',[1.00,0.65,0.00]);hold on;
        end
        
        ylim([-0.5 0.5])
        xlim([0.5 3])

        box off
        if ~(c == components)
            set(gca,'xtick',[])
        else
            %xlabel(['\bf Diet index'],'interpreter','latex')
            set(gca,'xtick',1:3,'xticklabel',diet_list)
            xtickangle(45)
        end
        if c == 1
            title('\bf Diet mode','interpreter','latex')
        end
    end
elseif (isempty(diet) && isempty(time))                                    % doing the same as above but if tensor is not 4th order but a matrix instead as for the unfolded tensor PCA case
        iter = 0;
    for c = 1:components
        %% Individual mode
        iter = iter + 1;
        subplot(components,2,iter)
        if ~isempty(cluster_ind)
            for clus = 1:max(cluster_ind)
                scatter(find(cluster_ind==clus), individual(cluster_ind==clus,c),'^','filled','MarkerFaceColor',col(clus,:),'MarkerEdgeColor','k');hold on;
            end
        else
            plot(1:size(individual,1),individual(:,c),'o');hold on;
        end
%         for i = 1:size(individual,1)
%             text(i,individual(i,c), num2str(i),'FontSize',13)
%         end
        if c == 1
            title('\textbf{Individual mode}','interpreter','latex')
        end
        ylabel(['\textbf{Component ', num2str(c),'}'],'interpreter','latex')
        
        if c == components
            xticks([1:4:17])
            xlabel(['\bf Individual index'],'interpreter','latex')
        else
            set(gca,'xtick',[])
        end
        xlim([1 17])
        box off
        %% The 'variable' mode
        iter = iter + 1;
        subplot(components,2,iter)
        if isempty(clusters_mets)
            plot(1:size(metabolite,1),metabolite(:,c),'o');hold on;
            for i = 1:size(metabolite,1)
                text(i,metabolite(i,c), num2str(i),'FontSize',13)
            end
        else
            %plot(1:size(metabolite,1),metabolite(:,c),'o');hold on;
            for i = 1:max(clusters_mets)
                scatter(find(clusters_mets==i),metabolite(clusters_mets==i,c),'o','filled','MarkerFaceColor',col(i+3,:), 'MarkerEdgeColor','k');hold on;
            end
%             for i = 1:size(metabolite,1)
%                 text(i,metabolite(i,c), num2str(i),'FontSize',13)
%             end
        end
        if c == 1
            title('\textbf{Variable mode}','interpreter','latex')
            
        end
        if (c == components)
            xticks([1 round(size(metabolite,1)/2) size(metabolite,1)])
            xlabel(['\bf Individual index'],'interpreter','latex')
        else
            set(gca,'xtick',[])
        end
        box off
        xlim([1 size(metabolite,1)])
    end
    
elseif (isempty(diet))
    iter = 0;
    for c = 1:components
        iter = iter + 1;
        subplot(components,3,iter)
        if ~isempty(cluster_ind)
            for clus = 1:2
                scatter(find(cluster_ind==clus), individual(cluster_ind==clus,c),'o','filled');hold on;
            end
        else
            plot(1:size(individual,1),individual(:,c),'o');hold on;
        end
        for i = 1:size(individual,1)
            text(i,individual(i,c), num2str(i),'FontSize',13)
        end
        title('Individual mode')
        ylabel(['Component ', num2str(c)])
        iter = iter + 1;
        subplot(components,3,iter)
        if isempty(clusters_mets)
            scatter(1:size(metabolite,1),metabolite(:,c),'o','filled');hold on;
            for i = 1:size(metabolite,1)
                text(i,metabolite(i,c), num2str(i),'FontSize',13)
            end
        else
            plot(1:size(metabolite,1),metabolite(:,c),'o');hold on;
            for i = 1:max(clusters_mets)
                plot(find(clusters_mets==i),metabolite(clusters_mets==i,c),'o','filled');hold on;
            end
            for i = 1:size(metabolite,1)
                text(i,metabolite(i,c), num2str(i),'FontSize',13)
            end
        end
        iter = iter + 1;
        title('Metabolite mode')
        subplot(components,3,iter)
        plot(1:size(time,1),time(:,c),'LineWidth',2);hold on
        title('Time mode')
    end
end

end
