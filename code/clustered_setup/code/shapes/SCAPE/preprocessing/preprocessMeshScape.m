function mesh = preprocessMeshScape(fname, suffix, force, params)

addpath(genpath('/home/haggaim/Work/Projects/GlobalHelpers/toolbox_graph'))
addpath(genpath('/home/haggaim/Work/Projects/GlobalHelpers/MeshProjectLocal'))
% default params
params.null = [];
n_farthest = getoptions(params, 'n_farthest', 500);  % number of farthest points for diameter calculation
diameter_prctile = getoptions(params, 'diameter_prctile', 75); % for calculation of diameter
smoothing_n_itersAGD = getoptions(params, 'smoothing_n_itersAGD', 200); % number of smoothing iterations for AGD
smoothing_n_itersGC = getoptions(params, 'smoothing_n_itersGC', 50); % number of smoothing iterations for GC
n_points = getoptions(params, 'n_points', 100); % number of points to extract
min_feature_dist = getoptions(params, 'min_feature_dist', 0.1); % minimal distance between features, relative to calculated diameter
save_figs = getoptions(params, 'save_figs', false); % save figures
visualize = ~isunix && getoptions(params, 'visualize', save_figs); % show results
agdExtremumType = getoptions(params, 'agdExtremumType', 'MinMax'); % Max-get only local maxima of the AGD, MinMax- max and min,Max2Min- max and only 2 min
sortCriticalFlag = getoptions(params, 'sortCriticalFlag', 'maxScore'); % sort by maxScore or distScore
saveGeoMat = getoptions(params,'saveGeoMat',false);

% init
[fnamePath, fnameName, fnameExt] = fileparts(fname);
fname_mat = fullfile(fnamePath,[fnameName suffix '.mat']); % preprocessed mat file name
fname_fig1 = fullfile(fnamePath,[fnameName suffix '_fig1']); % figures filenames
fname_fig2 = fullfile(fnamePath,[fnameName suffix '_fig2']); % figures filenames
fname_fig3 = fullfile(fnamePath,[fnameName suffix '_fig3']); % figures filenames
fname_fig4 = fullfile(fnamePath,[fnameName suffix '_fig4']); % figures filenames


% check if processed file already exists
if ~force && exist(fname_mat,'file')
    % load data
    data = load(fname_mat);
    % check if parameters match
    if isequal(params,data.params)
        fprintf('Found preprocessed data (%s) for %s.\n', suffix, fname)
        mesh = data.mesh; % compatible preprocessed data already exists
        return;
    else
        fprintf('Found INCOMPATIBLE preprocessed data (%s) for %s.\n', suffix, fname)
    end
end

% need to preprocess - do it
fprintf('Preprocessing %s...\n', fname);
[vertices, faces] = read_mesh(fname);
[vertices, faces] = check_face_vertex(vertices, faces); % fat and short matrices
mesh.m = struct('V',vertices,'F', faces,'Nv',max(size(vertices)));

% calculate all pairwise geodesic distances
fprintf('Calculating all pairwise distances (%d vertices)\n', mesh.m.Nv);
% calculate graph adjancency matrix
adj = triangulation2adjacency(mesh.m.F',mesh.m.V');
dist = graphallshortestpaths(adj,'directed',false);


% calculate diameter
% seed with the two farthest points
mesh.farthestPnts = zeros(n_farthest,1);
[~,indMax] = max(dist(:));
[mesh.farthestPnts(1), mesh.farthestPnts(2)] = ind2sub(size(dist),indMax);
for ii = 3:n_farthest
    indMin = min(dist(mesh.farthestPnts(1:ii-1),:),[],1);
    [~,mesh.farthestPnts(ii)] = max(indMin);
end

% calculate triangle areas
mesh.faceAreas = computeSurfAreas(mesh.m.V',mesh.m.F');

% calculate 1-ring areas
mesh.oneRingAreas = zeros(mesh.m.Nv,1);
for ii = 1:mesh.m.Nv
    ff = any(mesh.m.F==ii); % indices of faces in the ii'th vertex 1-ring
    mesh.oneRingAreas(ii) = (1/3)*sum(mesh.faceAreas(ff));
end

% calculate diameter- normalize by the areas of 1-ring
dist_farthest = dist(mesh.farthestPnts,mesh.farthestPnts);
mesh.diameter = prctile(dist_farthest(:), diameter_prctile);

% calculate AGD (average geodesic distance)
mesh.AGD_raw = dist*mesh.oneRingAreas;

% smoothing
mesh.AGD = perform_mesh_smoothing(mesh.m.F,mesh.m.V,mesh.AGD_raw,struct('niter_averaging',smoothing_n_itersAGD));
mesh.GaussCurvature = perform_mesh_smoothing(mesh.m.F,mesh.m.V,discrete_gaussian_curvature(mesh.m.V',mesh.m.F'),struct('niter_averaging',smoothing_n_itersGC));

% find critical points 
vring = compute_vertex_ring(mesh.m.F');
mesh.AGD_minima = false(mesh.m.Nv,1);
mesh.AGD_maxima = false(mesh.m.Nv,1);
mesh.GC_minima = false(mesh.m.Nv,1);
mesh.GC_maxima = false(mesh.m.Nv,1);
for ii = 1:mesh.m.Nv
    %currNeighborhood = vring{ii}; % one ring
    currNeighborhood = setdiff([vring{ii} vring{vring{ii}}], ii); % two ring
    mesh.AGD_maxima(ii) = all(mesh.AGD(ii)>=mesh.AGD(currNeighborhood));
    mesh.AGD_minima(ii) = all(mesh.AGD(ii)<=mesh.AGD(currNeighborhood));
    mesh.GC_minima(ii) = all(mesh.GaussCurvature(ii)>=mesh.GaussCurvature(currNeighborhood));
    mesh.GC_maxima(ii) = all(mesh.GaussCurvature(ii)<=mesh.GaussCurvature(currNeighborhood));
end
switch agdExtremumType
    case 'Max'
        mesh.AGD_critical = find(mesh.AGD_maxima);
        mesh.GC_critical = find(mesh.GC_maxima);
    case 'MinMax'
        mesh.AGD_critical = cat(1,find(mesh.AGD_minima),find(mesh.AGD_maxima));
        mesh.GC_critical = cat(1,find(mesh.GC_minima),find(mesh.GC_maxima));
    case 'Max2Min'
        mesh.AGD_critical = cat(1,find(mesh.AGD_minima),find(mesh.AGD_maxima));
        mesh.GC_critical = cat(1,find(mesh.GC_critical),find(mesh.GC_critical));
    otherwise
        error('Unknown agdExtremumType!');
end

% Concatanate the points- AGD and GC
mesh.criticalPoints = cat(1, mesh.AGD_critical, mesh.GC_critical);

% remove points closest in the combination of AGD and GC
distFarthest = dist(mesh.criticalPoints,mesh.criticalPoints);
validPnts = true(length(mesh.criticalPoints),1);
pruneMask = distFarthest(validPnts,validPnts)>=mesh.diameter*min_feature_dist;
pruneMaskSum = sum(~pruneMask);
while any(pruneMaskSum>1)
    [~,ff] = max(pruneMaskSum);
    validPnts(ff) = false;
    pruneMask(ff,:) = 1;
    pruneMask(:,ff) = 1;
    pruneMaskSum = sum(~pruneMask);
end

% valid AGD and GC points
validPntsAGD = validPnts(1:numel(mesh.AGD_critical));
validPntsGC = validPnts(numel(mesh.AGD_critical) + 1:end);
% take minimas of AGD
validMinimaAGD = find(validPntsAGD(1:length(find(mesh.AGD_minima))));
if strcmp(agdExtremumType, 'Max2Min')
    % take only 2 minimas
    if length(validMinimaAGD) > 2
        [ ~ , validMinimaIdx ]  = sort(mesh.AGD(mesh.AGD_critical(validMinimaAGD)));
        validPntsAGD(validMinimaAGD(validMinimaIdx(3:end))) = false;
        validMinimaAGD = validMinimaAGD(validMinimaIdx(1:2));
    end
end

% draw AGD before and after prunning
if visualize
    h = figure; axis equal; axis off;
    labels = cellstr( num2str([1:length(mesh.AGD_critical)]') );
    patch('vertices',mesh.m.V','faces',mesh.m.F','FaceVertexCData',mesh.AGD,'FaceColor','flat','EdgeColor','none','FaceAlpha',0.4);
    hold on;
    scatter3(mesh.m.V(1,mesh.AGD_critical), mesh.m.V(2,mesh.AGD_critical), mesh.m.V(3,mesh.AGD_critical),200,'r','filled');
    scatter3(mesh.m.V(1,mesh.AGD_critical(validPntsAGD)), mesh.m.V(2,mesh.AGD_critical(validPntsAGD)), mesh.m.V(3,mesh.AGD_critical(validPntsAGD)),100,'g','filled');
    text(mesh.m.V(1,mesh.AGD_critical), mesh.m.V(2,mesh.AGD_critical), mesh.m.V(3,mesh.AGD_critical), labels, 'VerticalAlignment','bottom',...
        'HorizontalAlignment','right','FontSize',20,'FontWeight','bold','color','k');
    cameratoolbar('show');
    cameratoolbar('SetMode','orbit')
    cameratoolbar('SetCoordSys','none')
    title('AGD points - before/after prunning');
    if save_figs
        saveas(h, fname_fig1, 'fig');
        saveas(h, fname_fig1, 'jpg');
        close(h);
    end
end

% draw GC before and after prunning
if visualize
    h = figure; axis equal; axis off;
    labels = cellstr( num2str([1:length(mesh.GC_critical)]') );
    patch('vertices',mesh.m.V','faces',mesh.m.F','FaceVertexCData',mesh.GaussCurvature,'FaceColor','flat','EdgeColor','none','FaceAlpha',0.4);
    hold on;
    scatter3(mesh.m.V(1,mesh.GC_critical), mesh.m.V(2,mesh.GC_critical), mesh.m.V(3,mesh.GC_critical),200,'r','filled');
    scatter3(mesh.m.V(1,mesh.GC_critical(validPntsGC)), mesh.m.V(2,mesh.GC_critical(validPntsGC)), mesh.m.V(3,mesh.GC_critical(validPntsGC)),100,'g','filled');
    text(mesh.m.V(1,mesh.GC_critical), mesh.m.V(2,mesh.GC_critical), mesh.m.V(3,mesh.GC_critical), labels, 'VerticalAlignment','bottom',...
        'HorizontalAlignment','right','FontSize',20,'FontWeight','bold','color','k');
    cameratoolbar('show');
    cameratoolbar('SetMode','orbit')
    cameratoolbar('SetCoordSys','none')
    title('GC points - before/after prunning');
    if save_figs
        saveas(h, fname_fig3, 'fig');
        saveas(h, fname_fig3, 'jpg');
        close(h);
    end
end

% draw Critical points before and after prunning
if visualize
    h = figure; axis equal; axis off;
    labels = cellstr( num2str([1:length(mesh.criticalPoints)]') );
    patch('vertices',mesh.m.V','faces',mesh.m.F','FaceVertexCData',mesh.GaussCurvature,'FaceColor','flat','EdgeColor','none','FaceAlpha',0.4);
    hold on;
    scatter3(mesh.m.V(1,mesh.criticalPoints), mesh.m.V(2,mesh.criticalPoints), mesh.m.V(3,mesh.criticalPoints),200,'r','filled');
    scatter3(mesh.m.V(1,mesh.criticalPoints(validPnts)), mesh.m.V(2,mesh.criticalPoints(validPnts)), mesh.m.V(3,mesh.criticalPoints(validPnts)),100,'g','filled');
    text(mesh.m.V(1,mesh.criticalPoints), mesh.m.V(2,mesh.criticalPoints), mesh.m.V(3,mesh.criticalPoints), labels, 'VerticalAlignment','bottom',...
        'HorizontalAlignment','right','FontSize',20,'FontWeight','bold','color','k');
    cameratoolbar('show');
    cameratoolbar('SetMode','orbit')
    cameratoolbar('SetCoordSys','none')
    title('Critical points - before/after prunning');
    if save_figs
        saveas(h, fname_fig4, 'fig');
        saveas(h, fname_fig4, 'jpg');
        close(h);
    end
end

% finally, remove the points
mesh.AGD_critical = mesh.AGD_critical(validPntsAGD);
mesh.GC_critical = mesh.GC_critical(validPntsGC);
mesh.criticalPoints = mesh.criticalPoints(validPnts);

switch sortCriticalFlag
    case 'distScore'
        % take farthest critical points
        v = zeros(size(mesh.criticalPoints));
        distFarthest = dist(mesh.criticalPoints,mesh.criticalPoints);
        [~,indMax] = max(distFarthest(:));
        [v(1), v(2)] = ind2sub(size(distFarthest),indMax);
        for ii = 3:length(v)
            indMin = min(distFarthest(v(1:ii-1),:),[],1);
            [~,v(ii)] = max(indMin);
        end
    case 'maxScore'
        % AGD critical scores
        mesh.AGD_score = [];
        mesh.GC_score = [];
        if strcmp(agdExtremumType, 'Max2Min')
            % in case it's Max2Min, sort only max, so minimas won't go away
            for ii = length(validMinimaAGD) + 1:length(mesh.AGD_critical)
                mesh.AGD_score(ii) = norm(mesh.AGD(mesh.AGD_critical(ii))-mesh.AGD(vring{mesh.AGD_critical(ii)}));
            end
            % sort scores
            [ ~, vAGD] = sort(mesh.AGD_score,'descend');
            vAGD = [(1:length(validMinimaAGD)) vAGD(1:end-length(validMinimaAGD))]; % the minimas would be zeros, it keeps them in
        else
            for ii = 1:length(mesh.AGD_critical)
                mesh.AGD_score(ii) = norm(mesh.AGD(mesh.AGD_critical(ii))-mesh.AGD(vring{mesh.AGD_critical(ii)}));
            end
            % sort scores
            [ ~, vAGD] = sort(mesh.AGD_score,'descend');
        end
        % GC scores
        for ii = 1:length(mesh.GC_critical)
            mesh.GC_score(ii) = norm(mesh.GaussCurvature(mesh.GC_critical(ii))-mesh.GaussCurvature(vring{mesh.GC_critical(ii)}));
        end
        % sort scores
        [ ~, vGC] = sort(mesh.GC_score,'descend');
        v = cat(2, vAGD, numel(vAGD) + vGC);
    otherwise
        error('Unsupported sortCriticalFlag!');
end
% sort by farthest points / best score
mesh.criticalPoints =  mesh.criticalPoints(v);

% gathering all points
mesh.featurePoints = zeros(n_points,1);
% more critical points than n_points- rejecting points
if length(mesh.criticalPoints)  > n_points
    mesh.featurePoints = mesh.criticalPoints(1:n_points);
else
    % less critical points than n_points- adding points
    mesh.featurePoints(1:length(mesh.criticalPoints),:) = mesh.criticalPoints;
    for ii = length(mesh.criticalPoints) + 1 : n_points
        indMin = min(dist(mesh.featurePoints(1:ii-1),:),[],1);
        [~,mesh.featurePoints(ii)] = max(indMin);
    end
end

% draw AGD & critical points
if visualize
    h = figure; axis equal; axis off;
    labels = cellstr( num2str([1:n_points]') );
    patch('vertices',mesh.m.V','faces',mesh.m.F','FaceVertexCData',mesh.AGD,'FaceColor','interp','EdgeColor','none','FaceAlpha',0.4);
    hold on;
    scatter3(mesh.m.V(1,mesh.featurePoints), mesh.m.V(2,mesh.featurePoints), mesh.m.V(3,mesh.featurePoints),200,'k','filled');
    scatter3(mesh.m.V(1,mesh.criticalPoints), mesh.m.V(2,mesh.criticalPoints), mesh.m.V(3,mesh.criticalPoints),200,'r','filled');
    text(mesh.m.V(1,mesh.featurePoints), mesh.m.V(2,mesh.featurePoints), mesh.m.V(3,mesh.featurePoints), labels, 'VerticalAlignment','bottom',...
        'HorizontalAlignment','right','FontSize',20,'FontWeight','bold','color','k');
    cameratoolbar('show');
    cameratoolbar('SetMode','orbit')
    cameratoolbar('SetCoordSys','none')
    title('All feature points');
    if save_figs
        saveas(h, fname_fig2, 'fig');
        saveas(h, fname_fig2, 'jpg');
        close(h);
    end
end
% save geodesic distances
if saveGeoMat
    mesh.m.GeoMat = dist;
else
    mesh.m.GeoMat = [];
end


% save preprocessed mesh
save(fname_mat,'mesh','params');
disp('Finished preprocessing')


