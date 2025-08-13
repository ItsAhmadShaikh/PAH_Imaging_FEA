function compute_effective_stress(stressFile, outputCSV)
% Computes von Mises effective stress from PostView export with extra columns

data = readmatrix(stressFile);
elemIDs = data(:,1); % first column is element ID
raw = data(:,8:end); % skip first 7 columns (based on your structure)
nElems = size(raw,1);
nCols = size(raw,2);

nSteps = nCols / 6;
if mod(nSteps,1) ~= 0
    error('Column count (%d) is not a multiple of 6 after skipping metadata.', nCols);
end

fprintf('Detected %d elements, %d time steps\n', nElems, nSteps);

% Preallocate
outTime = zeros(nElems*nSteps,1);
outElem = zeros(nElems*nSteps,1);
outEff  = zeros(nElems*nSteps,1);

rowIndex = 1;

for t = 1:nSteps
    cols = (1+(t-1)*6) : (t*6);
    sxx = raw(:,cols(1));
    syy = raw(:,cols(2));
    szz = raw(:,cols(3));
    sxy = raw(:,cols(4));
    syz = raw(:,cols(5));
    szx = raw(:,cols(6));

    % Von Mises
    effStress = sqrt(0.5*((sxx-syy).^2 + (syy-szz).^2 + (szz-sxx).^2 ...
                   + 6*(sxy.^2 + syz.^2 + szx.^2)));

    idxRange = rowIndex : (rowIndex+nElems-1);
    outTime(idxRange) = t;
    outElem(idxRange) = elemIDs;
    outEff(idxRange)  = effStress;
    rowIndex = rowIndex + nElems;
end

T = table(outTime, outElem, outEff, 'VariableNames', {'time','element_ID','effective_stress'});
writetable(T, outputCSV);
fprintf('✅ Effective stress saved: %s\n', outputCSV);

% Debug: show first few rows
disp(T(1:10,:));
end