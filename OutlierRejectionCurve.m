function data_imputed = OutlierRejectionCurve(data_original)
%OutlierRejectionCurve - Entire curve outlier identification and attenuation.
%
% Syntax:  data_imputed = OutlierRejectionCurve(data_original)
%
% Inputs:
%    data_original - The metabolomics data that contains outliers and needs
%    to be cleaned up
%
% Outputs:
%    data_imputed - The cleaned data.               
%
% Example: 
%    data_imputed = OutlierRejectionCurve(data_original)
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

%% Outlier rejection
data_original = permute(data_original,[4 1 2 3]);                                            % d,i,t,m
[data_imputed, imputed_index_1] = ScalingFunc(data_original);              % Scaling down the outlier curves
[data_imputed, imputed_index_2] = ScalingFunc(data_imputed);               % Doing it again for the remainers
%% Plot to check
% iter = 0;
% x_time = 1:size(data_original,3);
% figure()
% disp(['Nr of imputed cases ', num2str(length(imputed_index_2))])
% for m = imputed_index_2(:,2)'
%     iter = iter + 1;
%     disp(['Nr ', num2str(iter),' ', num2str(m)])
%     current_original = squeeze(data_original(imputed_index_2(iter,1),:,:,m));
%     current_imputed = squeeze(data_imputed(imputed_index_2(iter,1),:,:,m));
%     subplot(1,2,1)
%     plot(x_time,current_original(imputed_index_2(iter,3),:),'--','LineWidth',4);hold on
%     plot(x_time,current_original);hold off;
%     title('Original','Interpreter','latex')
%     xlabel('Tid (h)','Interpreter','latex')
%     ylabel('Intensitet','Interpreter','latex')
%     legend('Altered data')
%     
%     subplot(1,2,2)
%     plot(x_time,current_imputed(imputed_index_2(iter,3),:),'--','LineWidth',4);hold on;
%     plot(x_time,current_imputed);hold off;
%     title('Imputed','Interpreter','latex')
%     xlabel('Tid (h)','Interpreter','latex')
%     ylabel('Intensitet','Interpreter','latex')
%     legend('Altered data')
%     pause
% end


end

function [data_imputed, imputed_index] = ScalingFunc(data_original)
%% Outlier rejection
imputed_index = [];                                                        % Preallocate
data_imputed = data_original;                                              % Preallocate
x_time = 1:size(data_original,3);
for m = 1:size(data_original,4)                                            % Looping through the metabolites
    for d = 1:3                                                            % Looping through diets
        current_data = squeeze(data_original(d,:,:,m));                    % Current metabolite and diet
        for i = 1:size(data_original,2)                                    % Going through every individual curve
            tmp_data = current_data;
            tmp_data(i,:) = [];
            mean_curve_current = mean(tmp_data,1);                         % The mean of the population at every time point except the current individual
            factor3 = 3*std(tmp_data,[],1);                                % 3 Standard deviations of the population
            
            index_outlier_over = (mean_curve_current + factor3) < current_data(i,:); % Checking if over or under 3 standard deviations
            index_outlier_under = (mean_curve_current - factor3) > current_data(i,:);
            
            if sum(index_outlier_over) > 2                                 % If more hits than two it over is considered a outlier curve
                scaling_factor = mean(mean_curve_current)/mean(current_data(i,:)); % The scaling factor to attenuate the curve to the population level
                data_imputed(d,i,:,m) = current_data(i,:) * scaling_factor;% Scaling
                imputed_index = [imputed_index; d m i];
            elseif sum(index_outlier_under) > 3                            % If more hits than three under, is considered a outlier curve
                scaling_factor = mean(mean_curve_current)/mean(current_data(i,:));
                data_imputed(d,i,:,m) = current_data(i,:) * scaling_factor;% Scale and add
                imputed_index = [imputed_index; d m i];
            end
        end
    end
end
end
%------------- END CODE --------------