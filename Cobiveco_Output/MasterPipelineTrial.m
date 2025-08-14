%% === Master Pipeline: Validate + Visualize Effective Stress by Cobiveco Region ===
clear; clc;

% --- FILE PATHS ---
parcellationVTU = 'Y325_parcellation.vtu';     % Cobiveco VTU with parcellation labels
csvFile         = 'Y325output_region_stats.csv'; % Region-averaged stress stats

%% === STEP 0: Robust read of parcellation labels from VTU (no readVTU) ===
P = safeReadParcellation(parcellationVTU);   % returns P.labels (Nx1 double)
regionLabels = double(P.labels(:));

if isempty(regionLabels)
    error('Could not find pointData "parcellation" in %s', parcellationVTU);
end

% Basic report
fprintf('Regions found in VTU: %d total points\n', numel(regionLabels));
existing = unique(regionLabels(~isnan(regionLabels)));
fprintf('✅ Region IDs present: %s\n', mat2str(existing(:)'));

%% === STEP 1: Validate Region Mapping (counts per region id) ===
% NOTE: region IDs may be non-contiguous (e.g., [3..25, 28..50]).
maxID = max(existing);
counts = zeros(maxID+1,1);
valid = regionLabels(~isnan(regionLabels) & regionLabels>=0);
idx = floor(valid)+1;                 % assume integer IDs; +1 for 1-based indexing
counts = accumarray(idx, 1, [maxID+1, 1], @sum, 0);

figure;
bar(0:maxID, counts);                 % x-axis = actual region IDs (0..maxID)
title('Number of Points per Region'); xlabel('Region ID'); ylabel('Point Count'); grid on;

%% === STEP 2: Load Stress Data (rows = region IDs, cols = time steps) ===
fprintf('--- Loading Regional Stress CSV ---\n');
data = readmatrix(csvFile);
timeSteps    = data(1, 2:end);          % First row (excluding first col) = time steps
stressMatrix = data(2:end, 2:end);      % Rows = region IDs, Cols = time steps

% Sanity: keep only rows that exist in the VTU labels and fit within CSV
validRows = existing(existing <= size(stressMatrix,1)-1);  % region IDs present and in CSV
fprintf('CSV Shape: %d rows (region IDs 0..%d) x %d time steps\n', ...
    size(stressMatrix,1), size(stressMatrix,1)-1, size(stressMatrix,2));
fprintf('Plotting rows for region IDs present in VTU and CSV: %s\n', mat2str(validRows(:)'));

%% === STEP 3: Visualizations (lines, heatmap, geometry) ===
figure('Name','Regional Stress Visualization','Position',[60 60 1250 800]);

% --- Subplot 1: Line Trends ---
subplot(2,2,1); hold on;
colors = lines(numel(existing));
for k = 1:numel(existing)
    rID = existing(k);
    row = rID + 1;
    if row <= size(stressMatrix,1)
        plot(timeSteps, stressMatrix(row,:), 'LineWidth', 1.8, 'Color', colors(k,:), ...
            'DisplayName', sprintf('Region %d', rID));
    end
end
xlabel('Time Step'); ylabel('Effective Stress');
title('Regional Effective Stress Trends'); grid on;
legend('show','Location','northeastoutside');

% --- Subplot 2: Heatmap (only present regions) ---
subplot(2,2,3);
[sortedIDs, order] = sort(existing(:)');
rows = sortedIDs + 1;
imagesc(timeSteps, sortedIDs, stressMatrix(rows,:));
colormap(jet); hcb = colorbar; hcb.Label.String = 'Effective Stress';
title('Effective Stress by Region'); xlabel('Time Step'); ylabel('Region');

% --- Subplot 3: Geometry colored by region-average stress at a chosen time ---
% Choose a time index (change this to what you want)
[~, ti] = min(abs(timeSteps - timeSteps(1)));   % e.g., first time point
% [~, ti] = min(abs(timeSteps - 5));            % or: nearest to time=5

% Load VTU geometry + region labels (using vtkRead)
S = vtkRead(parcellationVTU);               % requires vtkToolbox on path
P = double(S.points);                       % Nx3
% Tetrahedra connectivity (cells)
T = double(S.cells);                        % Mx4 (should be tets)
% Region labels per point
if isfield(S,'pointData') && isfield(S.pointData,'parcellation')
    labels = double(S.pointData.parcellation(:));
else
    error('VTU does not contain pointData.parcellation');
end

% Build a per-region lookup for stress at time ti
regStress = nan(max([existing; 0])+1, 1);   % index by (regionID+1)
for k = 1:numel(existing)
    rID = existing(k);
    row = rID + 1;
    if row <= size(stressMatrix,1) && ti <= size(stressMatrix,2)
        regStress(row) = stressMatrix(row, ti);
    end
end

% Map region-average stress to each mesh point
% (points with labels not in CSV remain NaN and will show as gaps)
pointStress = nan(size(P,1),1);
validPts = labels >= 0 & (labels+1) <= numel(regStress);
pointStress(validPts) = regStress(labels(validPts)+1);

% Extract outer surface triangles from tetrahedra
TR = triangulation(T, P);
[Fb, Pb] = freeBoundary(TR);   % Fb: Kx3 faces; Pb: Kpts x 3 coords
% We need per-vertex colors for the surface vertices:
% map boundary vertices back to original point indices
bvi = TR.Points; %#ok<NASGU> % (not used, kept for clarity)
% Build a map from boundary vertex (row in Pb) to original point index:
% freeBoundary returns indices into TR.Points, which is P
[~, surfVertexIdx] = ismember(Pb, P, 'rows');
Csurf = pointStress(surfVertexIdx);

% Plot the colored surface
subplot(2,2,[2 4]); hold on;
hs = trisurf(Fb, Pb(:,1), Pb(:,2), Pb(:,3), Csurf, 'EdgeColor','none', 'FaceAlpha', 1.0);
axis equal tight vis3d; view(135, 20);
title(sprintf('Geometry colored by region-average stress (time index %d)', ti));
xlabel('X'); ylabel('Y'); zlabel('Z'); grid on;
colormap(jet); cb = colorbar; cb.Label.String = 'Avg Effective Stress (Region)';
lighting gouraud; camlight headlight;

sgtitle('Regional Effective Stress (Trends, Heatmap, Geometry)', 'FontWeight','bold');
%% ----------------- helper (robust VTU label reader) -----------------
function S = safeReadParcellation(vtuPath)
% Tries vtkRead first (if available), then regex parsing of ASCII VTU.
S = struct('labels', []);
useVTK = (exist('vtkRead','file')==3 || exist('vtkRead','file')==2);
if useVTK
    try
        V = vtkRead(vtuPath);
        % common field spellings
        cand = {'parcellation','Parcellation','PARCELLATION'};
        for k=1:numel(cand)
            if isfield(V,'pointData') && isfield(V.pointData, cand{k})
                S.labels = double(V.pointData.(cand{k})(:));
                return;
            end
        end
    catch
        % fall through
    end
end
% Fallback: text parse the DataArray named "parcellation"
txt = fileread(vtuPath);  % assumes ASCII VTU
patts = [ ...
    "Name=""parcellation""[^>]*>(.*?)</DataArray>"; ...
    "Name='parcellation'[^>]*>(.*?)</DataArray>"; ...
    "Name=""Parcellation""[^>]*>(.*?)</DataArray>"; ...
    "Name='Parcellation'[^>]*>(.*?)</DataArray>" ...
    ];
blk = [];
for p = 1:numel(patts)
    t = regexp(txt, patts(p), 'tokens', 'once');
    if ~isempty(t), blk = t{1}; break; end
end
if isempty(blk)
    warning('No DataArray named "parcellation" found in %s.', vtuPath);
    return;
end
S.labels = double(sscanf(blk, '%f'));
end