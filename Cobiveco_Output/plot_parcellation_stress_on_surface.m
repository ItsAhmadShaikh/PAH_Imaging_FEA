function plot_parcellation_stress_on_surface(vtuPath,csvPath,varargin)
% Pretty AHA/Cobiveco surface colored by average regional stress per time step.
% Requires vtkRead from vtkToolbox.
% Optional name-value:
%   'Surface'    : 'epicardium' (default) | 'all' | 'endocardium'
%   'ABepiThresh': 0.90
%   'ShowText'   : true
%   'SaveFrames' : false
%   'OutDir'     : 'frames'

p = inputParser;
p.addParameter('Surface','all');
p.addParameter('ABepiThresh',0.90);
p.addParameter('ShowText',true);
p.addParameter('SaveFrames',false);
p.addParameter('OutDir','frames');
p.parse(varargin{:});
opt = p.Results;

%% --- Load VTU ---
S = vtkRead(vtuPath);              
V = double(S.points);
T = double(S.cells);
parc = double(S.pointData.parcellation(:));
hasAB = isfield(S.pointData,'ab');
if hasAB, ab = double(S.pointData.ab(:)); else, ab = []; end

%% --- Build boundary faces
F_raw = [T(:,[1 2 3]);
         T(:,[1 2 4]);
         T(:,[1 3 4]);
         T(:,[2 3 4])];
F_sort = sort(F_raw,2);
[~,~,ic] = unique(F_sort,'rows');
count = accumarray(ic,1);
keep = count==1;                             
F = F_raw(keep(ic),:);

switch lower(opt.Surface)
    case 'epicardium'
        assert(hasAB,'Need pointData.ab to extract epicardium.');
        faceAB = mean(ab(F),2);
        F = F(faceAB >= opt.ABepiThresh,:);
    case 'endocardium'
        assert(hasAB,'Need pointData.ab to extract endocardium.');
        faceAB = mean(ab(F),2);
        F = F(faceAB <= (1-opt.ABepiThresh),:);
    case 'all'
        % keep all boundary faces
    otherwise
        error('Surface must be epicardium|endocardium|all');
end

faceReg = mode(parc(F),2);

%% --- Load CSV stress data ---
M = readmatrix(csvPath);
regIDs = M(2:end,1);
stressByTime = M(2:end, 2:end);  % rows = regions, cols = time steps
timeSteps = M(1, 2:end);

% Map region IDs to row indices in stressByTime
maxReg = max(regIDs);
reg2row = nan(maxReg,1);
reg2row(regIDs) = 1:length(regIDs);

%% --- Prepare figure
if opt.SaveFrames && ~exist(opt.OutDir,'dir')
    mkdir(opt.OutDir);
end

for tIdx = 1:length(timeSteps)
    % Stress for this time step
    reg2stress = nan(maxReg,1);
    reg2stress(regIDs) = stressByTime(:,tIdx);
    faceColor = reg2stress(faceReg);

    figure('Color',[0.11 0.13 0.18],'Position',[100 100 800 700]);
    trisurf(F,V(:,1),V(:,2),V(:,3),faceColor,...
        'EdgeColor','none','FaceColor','flat');
    axis equal off
    colormap(turbo); colorbar; 
    caxis([nanmin(faceColor) nanmax(faceColor)]);
    title(sprintf('Time %.2f: AHA/Cobiveco Regions (Effective Stress)', timeSteps(tIdx)),...
        'Color','w','FontSize',14)

    camlight headlight; camlight('left'); lighting gouraud;
    material([0.3 0.8 0.2]);      
    set(gca,'Color',[0.11 0.13 0.18])
    set(findobj(gcf,'Type','ColorBar'),'Color','w')

    if opt.ShowText
        regsOnSurf = unique(faceReg);
        for r = regsOnSurf(:).'
            fIdx = find(faceReg==r);
            if isempty(fIdx) || isnan(reg2stress(r)), continue; end
            cTri = mean(reshape(V(F(fIdx,:).',:),3,[],3),1);
            c = squeeze(mean(cTri,2)); 
            text(c(1),c(2),c(3), sprintf('%.2f',reg2stress(r)),...
                'Color','w','FontSize',9,'FontWeight','bold',...
                'HorizontalAlignment','center');
        end
    end

    view(130,20);

    if opt.SaveFrames
        saveas(gcf, fullfile(opt.OutDir, sprintf('frame_%03d.png',tIdx)));
        close(gcf);
    end
end
end