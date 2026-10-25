function L = DBESC(fea, k, M, p)
if nargin < 7
    bcsUpK = 60;
end
if nargin < 6
    bcsLowK = 20;
end

if nargin < 5
    gama=0;
end
if nargin < 4
    p = 16;
end

if nargin < 3
    M = 20;
end

if bcsUpK < bcsLowK
    bcsUpK = bcsLowK;
end


baseCls = USENC_EnsembleGeneration(fea, M, p);
clear fea

L = USENC_ConsensusFunction(baseCls,k);


function members = USENC_EnsembleGeneration(fea, M, p)

if nargin < 3
    p = 16;
end

if nargin < 2
    M = 20;
end

N = size(fea,1);
if p>N
    p = N;
end
members = zeros(N,M);


warning('off');


parfor i = 1:M
    CurtIndex = datasample(1:N, p, 'Replace', false);
    sdata = fea(CurtIndex,:);
    [~, W] = knnsearch(sdata, sdata, 'K', 2);
    min_dists = W(:, 2); 
    labels = zeros(size(fea, 1), 1);
    all_dists = pdist2(fea, sdata);
    clusters = 1:1:p;
    for j = 1:size(fea, 1)
        % Get distances for the current data point
        dists = all_dists(j, :);
        [~, idx_nearest] = min(dists);

        % Determine the radius for the nearest sphere
        radius = min_dists(idx_nearest);

        % Check if the point is inside any sphere using logical indexing
        inside_sphere_mask = dists <= radius;

        if any(inside_sphere_mask)
            % If inside spheres, find the closest sphere
            closest_sphere_dists = dists(inside_sphere_mask);
            [~, closest_sphere_idx] = min(closest_sphere_dists);
            ball_in = find(inside_sphere_mask, closest_sphere_idx, 'first');
            if length(ball_in)>1
                ball_in = ball_in(end);
            end
            labels(j) = clusters(ball_in); % Assign the label of the closest sphere
        else
            % If not inside any sphere, assign the label of the nearest point's cluster
            labels(j) = clusters(idx_nearest);
        end
    end


    members(:,i) = labels;

end


function labels = USENC_ConsensusFunction(baseCls,k,maxTcutKmIters,cntTcutKmReps)
% Huang Dong. Mar. 20, 2019.
% Combine the M base clusterings in baseCls to obtain the final clustering
% result (with k clusters).




if nargin < 4
    cntTcutKmReps = 3; 
end
if nargin < 3
    maxTcutKmIters = 100; % maxTcutKmIters and cntTcutKmReps are used to limit the iterations of the k-means discretization in Tcut.
end

[N,M] = size(baseCls);

maxCls = max(baseCls);
for i = 1:numel(maxCls)-1
    maxCls(i+1) = maxCls(i+1)+maxCls(i);
end

cntCls = maxCls(end);
baseCls(:,2:end) = baseCls(:,2:end) + repmat(maxCls(1:end-1),N,1); clear maxCls

% Build the bipartite graph.
B=sparse(repmat([1:N]',1,M),baseCls(:),1,N,cntCls); clear baseCls
colB = sum(B);
B(:,colB==0) = [];

% Cut the bipartite graph.
labels = Tcut_for_bipartite_graph(B,k,maxTcutKmIters,cntTcutKmReps);



function labels = Tcut_for_bipartite_graph(B,Nseg,maxKmIters,cntReps)
% B - |X|-by-|Y|, cross-affinity-matrix

if nargin < 4
    cntReps = 3;
end
if nargin < 3
    maxKmIters = 100;
end

[Nx,Ny] = size(B);
if Ny < Nseg
    error('Need more columns!');
end

dx = sum(B,2);
dx(dx==0) = 1e-10; % Just to make 1./dx feasible.
Dx = sparse(1:Nx,1:Nx,1./dx); clear dx
Wy = B'*Dx*B;

%%% compute Ncut eigenvectors
% normalized affinity matrix
d = sum(Wy,2);
D = sparse(1:Ny,1:Ny,1./sqrt(d)); clear d
nWy = D*Wy*D; clear Wy
nWy = (nWy+nWy')/2;

% computer eigenvectors
[evec,eval] = eig(full(nWy)); clear nWy   
[~,idx] = sort(diag(eval),'descend');
Ncut_evec = D*evec(:,idx(1:Nseg)); clear D

%%% compute the Ncut eigenvectors on the entire bipartite graph (transfer!)
evec = Dx * B * Ncut_evec; clear B Dx Ncut_evec

% normalize each row to unit norm
evec = bsxfun( @rdivide, evec, sqrt(sum(evec.*evec,2)) + 1e-10 );

% k-means
% labels = litekmeans(evec,Nseg);
labels = kmeans(evec,Nseg,'MaxIter',maxKmIters,'Replicates',cntReps);
