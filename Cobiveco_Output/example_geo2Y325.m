addpath('..');

%% Create cobiveco object and compute coordinates

c = cobiveco(struct('inPrefix','Rat Y325 geo2 example/Y325W8', 'outPrefix','result_geo2Y325/'));
S = vtkRead(sprintf('%s.vtu', c.cfg.inPrefix));

% Fix 0-based connectivity (VTK) -> 1-based (MATLAB)
if isfield(S,'cells') && ~isempty(S.cells) && min(S.cells(:)) == 0
    S.cells = S.cells + 1;
end

% Sanity checks
assert(all(isfinite(S.points(:))), 'VTU has NaN/Inf in Points');
assert(all(S.cells(:) >= 1), 'Connectivity has 0/negative indices');
assert(all(S.cells(:) <= size(S.points,1)), 'Connectivity index exceeds point count');

TR = triangulation(double(S.cells), double(S.points));  % proceed
c.prepareMesh0;
if c.cfg.CobivecoX == true
    c.computeAllCobivecoX;
else 
    c.computeAllCobiveco;
end

