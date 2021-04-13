function data_imputed = OutlierRejectionSpike(data_original)
%OutlierRejectionSpike - Spike outlier identification and attenuation.
%
% Syntax:  data_imputed = OutlierRejectionSpike(data_original)
%
% Inputs:
%    data_original - The metabolomics data that contains outliers and needs
%    to be cleaned up
%
% Outputs:
%    data_imputed - The cleaned data.               
%
% Example: 
%    data_imputed = OutlierRejectionSpike(data_original)
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
%% Imputed data 
data_imputed = data_original;                                              % Preallocate
%% Outlier rejection
imputed_index = [];
for m = 1:size(data_original,4)                                            % Looping through the metabolites
    for d = 1:3                                                            % Looping through the diets
        current_data = squeeze(data_original(d,:,:,m));                    % The current metabolite and diet
        for i = 1:size(data_original,2)                                    % One individual curve at a time
            tmp_data = current_data;                                       
            tmp_data(i,:) = [];
            mean_curve_current = mean(tmp_data,1);                         % Mean of the entire population except the current individual
            factor3 = 3*std(tmp_data,[],1);                                % 3 Standard deviations of the population
            index_outlier_over = (mean_curve_current + factor3) < current_data(i,:); % Is the current individual showing a peak at only one time point?
            index_outlier_under = (mean_curve_current - factor3) > current_data(i,:);
            if sum(index_outlier_over) == 1  && find(index_outlier_over)>1 && find(index_outlier_over)<8 % If the peak is over the population, impute it down to the mean point of the population
                    index_list = [find(index_outlier_over)-1 find(index_outlier_over)+1]; % The two points around the peak
                    
                    if any(index_list == size(data_original,3)+1)          % Fixing out of bounds index
                        index_list(index_list == size(data_original,3)+1) = size(data_original,3);
                    elseif any(index_list == 0)
                        index_list(index_list == 0) = 1;
                    end
                    
                    data_imputed(d,i,index_outlier_over,m) = mean(current_data(i,index_list)); % Mean of the two points will be the imputation
                    imputed_index = [imputed_index; d m i];
            elseif sum(index_outlier_under) == 1 && find(index_outlier_under)>1 && find(index_outlier_under)<8 % If the peak is under the population, impute it down to the mean point of the population
                index_list = [find(index_outlier_under)-1 find(index_outlier_under)+1];% Mean of the two points will be the imputation
                
                if any(index_list == size(data_original,3)+1)          % Fixing out of bounds index
                    index_list(index_list == size(data_original,3)+1) = size(data_original,3);
                elseif any(index_list == 0)
                    index_list(index_list == 0) = 1;
                end
                
                data_imputed(d,i,index_outlier_under,m) = mean(current_data(i,index_list));% Mean of the two points will be the imputation
                imputed_index = [imputed_index; d m i];
            end  
        end
    end
end

%------------- END CODE --------------
% %% Plot to check
% iter = 0;
% figure()
% x_time = 1:size(data_original,3);
% disp(['Nr of imputed cases ', num2str(length(imputed_index))])
% for m = imputed_index(:,2)'
%     iter = iter + 1;
%     disp(['Nr ', num2str(iter)])
%     current_original = squeeze(data_original(imputed_index(iter,1),:,:,m));
%     current_imputed = squeeze(data_imputed(imputed_index(iter,1),:,:,m));
%     subplot(1,2,1)
%     plot(x_time,current_original(imputed_index(iter,3),:),'--','LineWidth',4);hold on
%     plot(x_time,current_original);hold off;
%     title('Original','Interpreter','latex')
%     xlabel('Tid (h)','Interpreter','latex')
%     ylabel('Intensitet','Interpreter','latex')
%     legend('Altered data')
%     subplot(1,2,2)
%     plot(x_time,current_imputed(imputed_index(iter,3),:),'--','LineWidth',4);hold on;
%     plot(x_time,current_imputed);hold off;
%     title('Imputed','Interpreter','latex')
%     xlabel('Tid (h)','Interpreter','latex')
%     ylabel('Intensitet','Interpreter','latex')
%     legend('Altered data')
%     pause
% end

end