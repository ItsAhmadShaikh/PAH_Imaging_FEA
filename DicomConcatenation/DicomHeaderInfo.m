% You should open "sorted dicoms" --> "click on a rat of choice" --> 
% open one of the "Cine_FLASH_flc_..." folders and will get result
folderPath = uigetdir([], 'Select Folder Containing DICOM Files');
files = dir(fullfile(folderPath, '*.dcm'));
numFiles = length(files);

z_positions = zeros(numFiles, 1);

% Loop 1: Extract Z positions
for i = 1:numFiles
    filePath = fullfile(folderPath, files(i).name);
    info = dicominfo(filePath);
    z_positions(i) = info.ImagePositionPatient(3);
end

% Sort Z positions and compute differences
[z_sorted, sortIdx] = sort(z_positions);
z_diff = diff(z_sorted);
[max_spacing, max_idx] = max(z_diff);

% Loop 2: Display file info
for i = 1:numFiles
    filePath = fullfile(folderPath, files(i).name);
    info = dicominfo(filePath);
    fprintf('File: %s\n', files(i).name);
    fprintf('  Slice Thickness: %.2f mm\n', info.SliceThickness);
    fprintf('  Pixel Spacing: [%.2f, %.2f] mm\n', info.PixelSpacing);
    fprintf('  Image Position: [%.2f, %.2f, %.2f]\n', info.ImagePositionPatient);
end

% Print max spacing between slices
fprintf('\n📏 Largest Z-spacing:\n');
fprintf('  Between slice %d and %d: %.3f mm\n', max_idx, max_idx+1, max_spacing);
