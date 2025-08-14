function plot_stress_by_region(parcellationVTU, csvFile)
% PLOT_STRESS_BY_REGION
% Robust pipeline to (1) read a VTU with Cobiveco parcellation labels and
% (2) plot regional effective stress time series + heatmap from a CSV.
%
% CSV layout expected:
%   Row 1: [ignored or region header]  t1  t2  t3 ...
%   Row 2+: [region id or label]       s11 s12 s13 ...
% i.e., timeSteps = data(1,2:end); stressMatrix = data(2:end,2:end);

fprintf('=== Regional Stress Pipeline ===\n');
fprintf('VTU: %s\n', parcellationVTU);
fprintf('CSV: %s\n', csvFile);

%% --- Ensure vtkRead is available (best-effort bootstrap) ---
if ~(exist('vtkRead','file')==2 || exist('vtkRead','file')==3)
    candidates = { ...
        '/Users/kristengarcia/Documents/Cobiveco_RatAddition/dependencies/vtkToolbox', ...
        '/Users/kristengarcia/Documents/Cobiveco_RatAddition/dependencies/vtkToolbox/MATLAB', ...
        fullfile(fileparts(mfilename('fullpath')), 'dependencies', 'vtkToolbox') ...
    };
    for p = candidates
        if isfolder(p{1}), addpath(genpath(p{1})); end
    end
    rehash toolboxcache;
end

%% --- Read VTU (Option B: try/catch the actual call) ---
S = tryReadVTU(parcellationVTU);

%% --- Extract parcellation labels (case-insensitive + formats) ---
pd = getPointDataStruct(S);  % handles pointData / point_data
parc = getArrayCaseInsensitive(pd, 'parcellation');

% Flatten to numeric column
parc = makeNumericVector(parc);
assert(~isempty(parc), 'Could not find/parse pointData.parcellation in VTU.');
fprintf('✅ Found parcellation array with %d points.\n', numel(parc));

% Basic sanity: VTU piece point count (if present)
nPtsVTU = getNumberOfPointsIfAvailable(S);
if ~isempty(nPtsVTU) && nPtsVTU ~= numel(parc)
    warning('Parcellation length (%d) differs from VTU point count (%d). Proceeding anyway...', ...
        numel(parc), nPtsVTU);
end

% Regions + counts
[uniqueRegions, ~, regionIdx] = unique(parc(:));
regionCounts = accumarray(regionIdx, 1);

fprintf('Regions present: %s\n', mat2str(uniqueRegions.'));
fprintf('Points per region (in same order):\n');
disp([uniqueRegions, regionCounts]); %#ok<DISPLAY>

figure('Name','Parcellation Point Counts','Position',[100 100 700 350]);
bar(uniqueRegions, regionCounts, 'FaceAlpha',0.8);
xlabel('Region'); ylabel('Point Count'); grid on; title('Number of Points per Region');

%% --- Load CSV stress data ---
assert(exist(csvFile,'file')==2, 'CSV not found: %s', csvFile);
data = readmatrix(csvFile);
assert(~isempty(data) && size(data,2) >= 2 && size(data,1) >= 2, ...
    'CSV must have at least 2 rows and 2 columns.');

% Flexible header handling:
% time in (1,2:end)
timeSteps = data(1,2:end);
% stressMatrix rows map to regions (in the same order your CSV uses)
stressMatrix = data(2:end, 2:end);

% Optional: If the first column (row headers) are region IDs, you can cross-check:
csvRegionIDs = data(2:end,1);
if ~any(isnan(csvRegionIDs))
    % Try to align plotting legend with CSV region IDs
    csvRegionIDs = csvRegionIDs(:).';
else
    csvRegionIDs = 1:size(stressMatrix,1);
end

fprintf('CSV shape: %d regions x %d time steps\n', size(stressMatrix,1), size(stressMatrix,2));

%% --- Plot: lines + heatmap ---
figure('Name', 'Regional Stress Visualization', 'Position', [80 80 950 750]);

% Subplot 1: Line plots
subplot(2,1,1); hold on;
C = lines(size(stressMatrix,1));
for r = 1:size(stressMatrix,1)
    plot(timeSteps, stressMatrix(r,:), 'LineWidth', 2, 'Color', C(r,:), ...
        'DisplayName', sprintf('Region %d', csvRegionIDs(r)));
end
xlabel('Time Step'); ylabel('Effective Stress');
title('Regional Effective Stress Trends'); grid on;
legend('show', 'Location', 'northoutside', 'NumColumns', 4);

% Subplot 2: Heatmap
subplot(2,1,2);
imagesc(timeSteps, 1:size(stressMatrix,1), stressMatrix);
set(gca, 'YTick', 1:size(stressMatrix,1), 'YTickLabel', arrayfun(@(x)sprintf('%d',x), csvRegionIDs, 'uni',0));
xlabel('Time Step'); ylabel('Region (CSV IDs)'); title('Effective Stress by Region (Absolute)');
colormap(jet); cb=colorbar; cb.Label.String = 'Effective Stress';

% Annotate values (optional; comment out if too dense)
for r = 1:size(stressMatrix,1)
    for t = 1:size(stressMatrix,2)
        text(timeSteps(t), r, sprintf('%.2f', stressMatrix(r,t)), ...
            'Color','w','FontSize',8, 'HorizontalAlignment','center');
    end
end

sgtitle('Regional Effective Stress Visualization','FontWeight','bold');

fprintf('✅ Visualization Complete!\n');

end

%% ==================== helpers ====================

function S = tryReadVTU(vtuPath)
% Try/catch reading VTU; give clear errors for missing vtkRead or bad path.
try
    S = vtkRead(vtuPath);  %#ok<NASGU>  % probe first to trigger a mex/m issue early
catch ME
    if ~(exist('vtkRead','file')==2 || exist('vtkRead','file')==3)
        error(['vtkRead not found. Please add vtkToolbox to your path, e.g.:' newline ...
            'addpath(''/Users/kristengarcia/Documents/Cobiveco_RatAddition/dependencies/vtkToolbox'');']);
    end
    % If vtkRead exists but path is wrong or file invalid, rethrow original
    rethrow(ME);
end
% actually read if the first probe didn't assign
S = vtkRead(vtuPath);
end

function pd = getPointDataStruct(S)
% Handle vtkToolbox structs that may use pointData or point_data
if isfield(S,'pointData')
    pd = S.pointData;
elseif isfield(S,'point_data')
    pd = S.point_data;
else
    error('VTU missing pointData/point_data section.');
end
end

function arr = getArrayCaseInsensitive(pd, name)
% Find array in pointData ignoring case, and return its numeric content.
fns = fieldnames(pd);
match = strcmpi(fns, name);
if ~any(match)
    % Sometimes DataArrays are stored in pd.scalars or pd.other arrays; try to search structs
    % Fallback: scan subfields for a match
    for k = 1:numel(fns)
        if isstruct(pd.(fns{k})) && strcmpi(fns{k}, name)
            arr = pd.(fns{k});
            return;
        end
    end
    error('Could not find pointData field (case-insensitive): %s', name);
end
arr = pd.(fns{match});
end

function v = makeNumericVector(x)
% Convert VTU DataArray (various shapes/cell formats) to a numeric column vector.
if isnumeric(x)
    v = x(:);
    return;
end
if iscell(x)
    try
        v = cellfun(@double, x(:));
        return;
    catch
        % flatten any nested and try str2double
        v = str2double(string(x(:)));
        return;
    end
end
if isstring(x) || ischar(x)
    v = str2double(splitlines(string(x)));
    v = v(~isnan(v));
    v = v(:);
    return;
end
try
    v = double(x(:));
catch
    v = [];
end
end

function n = getNumberOfPointsIfAvailable(S)
% Best-effort to retrieve the declared number of points in VTU (for a sanity check).
n = [];
try
    if isfield(S,'points')
        n = size(S.points,1);
    elseif isfield(S,'UnstructuredGrid') && isfield(S.UnstructuredGrid,'Piece') ...
            && isfield(S.UnstructuredGrid.Piece,'NumberOfPoints')
        n = S.UnstructuredGrid.Piece.NumberOfPoints;
    end
catch
    n = [];
end
end