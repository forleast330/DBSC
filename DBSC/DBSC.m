function labels = DBSC(data, k, p, gama)
    if nargin < 3
        p = min(1000, size(data, 1));
    end
    
    if nargin < 4 || gama==0
        gama = mean(min_dists);
    end

    p = min(size(data, 1), p);
    if p < size(data, 1)
        % [~, sdata] = kmeans(data, p, 'MaxIter', 20);
         CurtIndex = datasample(1:size(data, 1), p, 'Replace', false);
         sdata = data(CurtIndex,:);

    [~, W] = knnsearch(sdata, sdata, 'K', 2);
    min_dists = W(:, 2); 
    

    % Perform spectral clustering on sdata
    clusters = spectralcluster(sdata, k, 'KernelScale', gama);

    % dist_temp = pdist(sdata);
    % dist_temp = squareform(dist_temp);
    % min_dists_temp = min_dists;
    % dist_temp = dist_temp./(min_dists_temp*min_dists_temp');
    % dist_temp = dist_temp./gama ;
    % S = exp(-dist_temp.^2);
    % issymmetric(S);
    % clusters = spectralcluster(S,k,'Distance','precomputed');

    % gscatter(sdata(:,1), sdata(:,2),clusters);
    % Calculate distances to find the nearest points from the original data to sdata
    labels = zeros(size(data, 1), 1);

    % Pre-compute distances from all data points to sdata
    % [all_dists, I] = pdist2(sdata, data,'minkowski',2,'Smallest',2);
    % all_dists = all_dists';
    % I = I';

    % all_dists = pdist2(data, sdata);
    [all_dists, I] = pdist2(sdata, data,'minkowski',2,'Smallest',1);
    all_dists = all_dists';
    I = I';
    for i = 1:size(data, 1)
        % Get distances for the current data point
        dists = all_dists(i, :);
        % dists = pdist2(data(i,:),sdata);
        % [~, idx_nearest] = min(dists);

        % Determine the radius for the nearest sphere
        % radius = min_dists(idx_nearest);

        % Check if the point is inside any sphere using logical indexing
        % inside_sphere_mask = dists <= radius;
        inside_sphere_mask = dists<=min_dists(I(i,:))';
        if any(inside_sphere_mask)
            % If inside spheres, find the closest sphere
            closest_sphere_dists = dists(inside_sphere_mask);
            [~, closest_sphere_idx] = min(closest_sphere_dists);
            ball_in = find(inside_sphere_mask, closest_sphere_idx, 'first');
            if length(ball_in)>1
                ball_in = ball_in(end);
            end
            % labels(i) = clusters(ball_in); % Assign the label of the closest sphere
            labels(i) = clusters(I(i,ball_in));
        else
            % If not inside any sphere, assign the label of the nearest point's cluster
            [~, idx_nearest] = min(dists);
            % labels(i) = clusters(idx_nearest);
            labels(i) = clusters(I(i,idx_nearest));
        end
        % clear dists
    end
    % clear all_dists
    % gscatter(data(:,1), data(:,2), labels);
    % hold on;
    % for i = 1:size(sdata,1)
    %     viscircles(sdata(i,:), min_dists(i), 'Color', 'k','LineWidth',1);
    % end
    % save("dataset\dbsc_plot_4C.mat","data","sdata","min_dists","labels")
    % plot(sdata(:,1), sdata(:,2), 'r.');
    
    else
        labels = spectralcluster(data, k, 'KernelScale', gama);
    end

end





%%


% function labels = DBSC(data, k, p, gama)
%     if nargin < 3
%         p = min(1000, size(data, 1));
%     end
% 
%     p = min(size(data, 1), p);
%     if p < size(data, 1)
%         % [~, sdata] = kmeans(data, p, 'MaxIter', 20);
%          CurtIndex = datasample(1:size(data, 1), p, 'Replace', false);
%         sdata = data(CurtIndex,:);
%     else
%         sdata = data;
%     end
% 
%     [~, W] = knnsearch(sdata, sdata, 'K', 2);
%     min_dists = W(:, 2); 
%     if nargin < 4 || gama==0
%         gama = mean(min_dists);
%     end
% 
%     % Perform spectral clustering on sdata
%     clusters = spectralcluster(sdata, k, 'KernelScale', gama);
% 
%     % Calculate distances to find the nearest points from the original data to sdata
%     labels = zeros(size(data, 1), 1);
% 
%     % Pre-compute distances from all data points to sdata
%     all_dists = pdist2(data, sdata);
% 
%     for i = 1:size(data, 1)
%         % Get distances for the current data point
%         dists = all_dists(i, :);
%         [~, idx_nearest] = min(dists);
% 
%         % Determine the radius for the nearest sphere
%         radius = min_dists(idx_nearest);
% 
%         % Check if the point is inside any sphere using logical indexing
%         inside_sphere_mask = dists <= radius;
% 
%         if any(inside_sphere_mask)
%             % If inside spheres, find the closest sphere
%             closest_sphere_dists = dists(inside_sphere_mask);
%             [~, closest_sphere_idx] = min(closest_sphere_dists);
%             ball_in = find(inside_sphere_mask, closest_sphere_idx, 'first');
%             if length(ball_in)>1
%                 ball_in = ball_in(end);
%             end
%             labels(i) = clusters(ball_in); % Assign the label of the closest sphere
%         else
%             % If not inside any sphere, assign the label of the nearest point's cluster
%             labels(i) = clusters(idx_nearest);
%         end
%     end
% end

























%% used
% function labels=DDSC(data, k, p, gama)
% if nargin < 3
%     p=min(1000,size(data,1));
% end
% 
% p = min(size(data,1), p);
% if p<size(data,1)
%     [~, sdata] = kmeans(data,p,'MaxIter',20);
% else
%     sdata = data;
% end
% [~, W] = knnsearch(sdata, sdata, 'K', 2);
% min_dists = W(:,2); 
% 
% if nargin < 4 || gama==0
%     gama = mean(min_dists);
% end
% 
% labels = zeros(size(data,1), 1);
% for i = 1:size(data,1)
%     dists = pdist2(data(i,:), sdata);
%     [~, idx] = min(dists);
%     if dists(idx) <= min_dists(idx)
%         labels(i) = idx;
%     else
%         labels(i) = -1; 
%     end
% end
% 
% 
% 
% % clusters = spectralcluster(sdata, k, "LaplacianNormalization",'symmetric','KernelScale',gama);
% clusters = spectralcluster(sdata, k,'KernelScale',gama);
% % gscatter(sdata(:,1), sdata(:,2), clusters);
% 
% for i = 1:size(data,1)
%     if labels(i) == -1
%         dists = pdist2(data(i,:), sdata);
%         [~, idx] = min(dists);
%         labels(i) = clusters(idx);
%      else
%         labels(i) = clusters(labels(i));
%     end
% end

%%
% hold on;
% for i = 1:size(sdata,1)
%     viscircles(sdata(i,:), min_dists(i), 'Color', 'r');
% end
% plot(data(:,1), data(:,2), 'b.');
% plot(sdata(:,1), sdata(:,2), 'ro');
% hold on;
% gscatter(data(:,1), data(:,2), labels);
% title('Spectral Clustering Result');


% end









% function labels = DDSC(data, k, p, gama)
%     if nargin < 3
%         p = min(1000, size(data, 1));
%     end
% 
%     p = min(size(data, 1), p);
% 
%     if p < size(data, 1)
%         [~, sdata] = kmeans(data, p, 'MaxIter', 20);
%     else
%         sdata = data;
%     end
% 
%     % 
%     [~, W] = knnsearch(sdata, sdata, 'K', gama);
%     min_dists = W(:, end); 
% 
%     if nargin < 4 || gama == 0
%         gama = mean(min_dists);
%     end
% 
%     % 
%     num_samples = size(sdata, 1);
%     similarity_matrix = zeros(num_samples);
% 
%     % 
%     radii = min_dists;
% 
%     % 
%     for i = 1:num_samples
%         for j = i+1:num_samples
%             dist = norm(sdata(i, :) - sdata(j, :));
%             if dist <= (radii(i) + radii(j))
%                 similarity_matrix(i, j) = 1;
%                 similarity_matrix(j, i) = similarity_matrix(i, j); % 对称矩阵
%             end
%         end
%     end
% 
%     % 
%     clusters = spectralcluster(similarity_matrix, k, 'LaplacianNormalization', 'symmetric','Distance','precomputed');
%     % gscatter(sdata(:,1), sdata(:,2), clusters);
%     % 
%     labels = zeros(size(data, 1), 1);
%     for i = 1:size(data, 1)
%         dists = pdist2(data(i, :), sdata);
%         [~, idx] = min(dists);
%         if dists(idx) <= radii(idx)
%             labels(i) = idx;
%         else
%             labels(i) = -1; 
%         end
%     end
% 
%     % 
%     for i = 1:size(data, 1)
%         if labels(i) == -1
%             dists = pdist2(data(i, :), sdata);
%             [~, idx] = min(dists);
%             labels(i) = clusters(idx);
%         else
%             labels(i) = clusters(labels(i));
%         end
%     end
% end
