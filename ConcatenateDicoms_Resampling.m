clear; clc; close all;

%% --- Setup ---
fprintf('Choose data Folder \n');
dataPath = uigetdir([], 'Choose folder containing Cine DICOM series');
fprintf('Choose Save Folder \n');
savePath = uigetdir([], 'Choose where to save the .nrrd');
output_filename = input('Enter a name for the nrrd file: ',"s");
target_spacing_z = 1.0;  % mm — uniform Z-spacing for output, target slice thickness

fullVolume = [];
all_slice_positions = {};
num_slices_per_scan = [];
num_timeframes = [];

%% --- Read DICOMs from all Cine series ---
fprintf('Reading DICOM headers...\n');
for i = 1:5
    match = dir(fullfile(dataPath, sprintf('Cine_FLASH_flc_%02d_Series*', i)));
    if isempty(match)
        warning('No folder matched for series %d', i);
        continue;
    end

    currentPath = fullfile(match(1).folder, match(1).name);
    dcmList = dir(fullfile(currentPath, '*.dcm'));
    if isempty(dcmList)
        warning('No DICOMs found in %s', currentPath);
        continue;
    end

    [~, idx] = sort({dcmList.name});
    dcmList = dcmList(idx);

    pos = [];
    img_no = [];

    for t = 1:length(dcmList)
        currentdicom = fullfile(currentPath, dcmList(t).name);
        try
            cur = dicominfo(currentdicom);
        catch
            warning('Unreadable: %s', dcmList(t).name); continue;
        end

        if isfield(cur, 'ImagePositionPatient')
            pos = [pos; cur.ImagePositionPatient(3)];
        end

        if isfield(cur, 'ImageComments') && contains(cur.ImageComments, 'frame')
            img_no = [img_no; str2double(extractAfter(cur.ImageComments, 'frame '))];
        else
            img_no = [img_no; NaN];
        end

        if i == 1 && t == 1
            spacing_x = cur.PixelSpacing(1);
            spacing_y = cur.PixelSpacing(2);
        end
    end

    all_slice_positions{i} = pos;
    num_timeframes(i) = max(img_no);
    num_slices_per_scan(i) = nnz(img_no == num_timeframes(i));
end

valid_scans = find(num_timeframes > 0);
if any(diff(num_timeframes(valid_scans)))
    error('Unequal number of timeframes across valid series.');
end
num_timeframes = num_timeframes(valid_scans(1));

%% --- Slice Spacing and Z Resampling ---
all_z = unique(round(cell2mat(all_slice_positions), 4));
if isempty(all_z)
    error('No Z-positions found. Check DICOM headers.');
end

z_range = [min(all_z), max(all_z)];
new_z_positions = z_range(1):target_spacing_z:z_range(2);
num_new_slices = length(new_z_positions);

%% --- Build Full Volume (T, X, Y, Z) ---
fprintf('Building full volume...\n');
for t = 1:num_timeframes
    currentVolume = [];
    slice_count = 1;

    for i = 1:5
        match = dir(fullfile(dataPath, sprintf('Cine_FLASH_flc_%02d_Series*', i)));
        if isempty(match), continue; end
        currentPath = fullfile(match(1).folder, match(1).name);

        dcmList = dir(fullfile(currentPath, '*.dcm'));
        [~, idx] = sort({dcmList.name});
        dcmList = dcmList(idx);

        for j = 1:num_slices_per_scan(i)
            idx_dcm = 1 + (j - 1) * num_timeframes + (t - 1);
            if idx_dcm > length(dcmList), continue; end

            try
                img = dicomread(fullfile(currentPath, dcmList(idx_dcm).name));
                currentVolume(:, :, 1, slice_count) = img;
                slice_count = slice_count + 1;
            catch
                warning('Error reading image for slice %d', j);
            end
        end
    end

    fullVolume(t, :, :, :) = squeeze(currentVolume);
end

%% --- Interpolate to uniform Z spacing ---
fprintf('Resampling to %.2f mm slice spacing...\n', target_spacing_z);
[nT, nx, ny, nz] = size(fullVolume);
resampledVolume = zeros(nT, nx, ny, num_new_slices, 'like', fullVolume);

orig_z = sort(unique(round(cell2mat(all_slice_positions), 4)));

for t = 1:nT
    current = squeeze(fullVolume(t, :, :, :));   % [X Y Z]
    current = permute(current, [2 1 3]);         % [Y X Z]

    for y = 1:ny
        for x = 1:nx
            profile = squeeze(current(y, x, :));
            resampled_voxel = interp1(orig_z, profile, new_z_positions, 'linear', 'extrap');
            resampledVolume(t, x, y, :) = resampled_voxel;
        end
    end
end

%% --- Save as NRRD ---
fprintf('Saving NRRD file...\n');
img = struct();  % ← this fixes the dot indexing error
img.pixelData = resampledVolume;
img.ijkToLpsTransform = [ spacing_x 0 0 0;
                          0 spacing_y 0 0;
                          0 0 target_spacing_z 0;
                          0 0 0 1 ];
cd(savePath);
nrrdwrite(strcat(output_filename, '.nrrd'), img);
fprintf('✅ DONE: Saved %s\n', fullfile(savePath, [output_filename '.nrrd']));

disp('Volume size (T, X, Y, Z):');
disp(size(resampledVolume));
