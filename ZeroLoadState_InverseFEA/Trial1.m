clear all
clc
close all

% % Plot settings
fontSize=15;
faceAlpha1=0.8;
faceAlpha2=0.3;
markerSize=40;
lineWidth=3;

% % Control parameters
% Path names
% savePath=uigetdir();
savePath=fullfile(pwd,"FeBio");

% Defining file names
febioFebFileNamePart='tempModel';
febioFebFileName=fullfile(savePath,[febioFebFileNamePart,'.feb']); %FEB file name
febioLogFileName=[febioFebFileNamePart,'.txt']; %FEBio log file name
febioLogFileName_disp=[febioFebFileNamePart,'_disp_out.txt']; %Log file name for exporting displacement
febioLogFileName_force=[febioFebFileNamePart,'_force_out.txt']; %Log file name for exporting force
febioLogFileName_stress=[febioFebFileNamePart,'_stress_out.txt']; %Log file name for exporting stress

% %Material parameter set
% c1=1e-3; %Shear-modulus-like parameter
% m1=8; %Material parameter setting degree of non-linearity
% k_factor=1e2; %Bulk modulus factor
% k=c1*k_factor; %Bulk modulus


% FEA control settings
analysisType='DYNAMIC';
numTimeSteps=10; %Number of time steps desired
max_refs=25; %Max reforms
max_ups=10; %Set to zero to use full-Newton iterations
opt_iter=15; %Optimum number of iterations
max_retries=6; %Maximum number of retires
dtmin=(1/numTimeSteps)/100; %Minimum time step size
dtmax=1/numTimeSteps; %Maximum time step size
min_residual=1e-3; %1e-10;
symmetric_stiffness=0;
runMode='internal';% 'internal' or 'external'




% Load the stl file
[BaseStruct] = import_STL("Meshes\base.stl");
[LVStruct] = import_STL("Meshes\LV.stl");
[RVStruct] = import_STL("Meshes\RV.stl");
[epiStruct] = import_STL("Meshes\epi.stl");

F_base = BaseStruct.solidFaces{1};
V_base = BaseStruct.solidVertices{1};

F_LV = LVStruct.solidFaces{1};
V_LV = LVStruct.solidVertices{1};

F_RV = RVStruct.solidFaces{1};
V_RV = RVStruct.solidVertices{1};

F_epi = epiStruct.solidFaces{1};
V_epi = epiStruct.solidVertices{1};

C_base=1*ones(size(F_base,1),1);
C_LV=2*ones(size(F_LV,1),1);
C_RV=3*ones(size(F_RV,1),1);
C_epi=4*ones(size(F_epi,1),1);

[F,V,C]=joinElementSets({F_base,F_LV,F_RV,F_epi},{V_base,V_LV,V_RV,V_epi},{C_base,C_LV,C_RV,C_epi}); %joining sets together
%merging sets
[F,V]=mergeVertices(F,V);

% % % % % % % % % Visualize heart surface
% % % % % % % cFigure; hold on;
% % % % % % % gpatch(F,V,C,'k',1);
% % % % % % % axisGeom; camlight headlight;
% % % % % % % drawnow;
% % % % % % % 
% % % % % % % [N,P,NV]=patchNormal(F,V);
% % % % % % % s=mean(patchEdgeLengths(F,V));
% % % % % % % quiverVec(P,N,s,'k')


% % Mesh using tetgen
%Find interior point
V_inner=getInnerPoint(F,V);


% % % % % % % % % Visualize interior point
% % % % % % % cFigure; hold on;
% % % % % % % gpatch(F,V,'w','none',0.5);
% % % % % % % plotV(V_inner,'r.','MarkerSize',25)
% % % % % % % axisGeom; camlight headlight;
% % % % % % % drawnow;

%%

tetVolume=tetVolMeanEst(F,V); %Volume for regular tets

tetGenStruct.stringOpt='-pq1.2AaY';
tetGenStruct.Faces=F;
tetGenStruct.Nodes=V;
tetGenStruct.holePoints=[];
tetGenStruct.faceBoundaryMarker=C; %Face boundary markers
tetGenStruct.regionPoints=V_inner; %region points
tetGenStruct.regionA=tetVolume;

[meshOutput]=runTetGen(tetGenStruct); %Run tetGen

% Access elements, nodes, and boundary faces
E=meshOutput.elements;
V=meshOutput.nodes;
Fb=meshOutput.facesBoundary;
Cb=meshOutput.boundaryMarker;
CE=meshOutput.elementMaterialID;

% % % % % % % meshView(meshOutput);

% % % % % % % cFigure; hold on;
% % % % % % % gpatch(Fb,V,Cb,'k',1);
% % % % % % % axisGeom; camlight headlight;
% % % % % % % drawnow;


%%

%BCs
F_base_BC = Fb(Cb==1,:);
bcSupportList = unique(F_base_BC(:));

%PressureSurfaces
% F_LV_pressure=Fb(Cb==2,:);
F_LV_pressure= fliplr(Fb(Cb==2,:));

% F_RV_pressure=Fb(Cb==3,:);
F_RV_pressure= fliplr(Fb(Cb==3,:));

% Visualize BC and Load
hf=cFigure;
title('Boundary conditions and loads','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hold on;

gpatch(Fb,V,'kw','none',0.5);

hl(1)=plotV(V(bcSupportList,:),'k.','MarkerSize',markerSize/2);
hl(2)=gpatch(F_LV_pressure,V,'r','k',1);
hl(3)=gpatch(F_RV_pressure,V,'g','k',1);
patchNormPlot(F_LV_pressure,V);
patchNormPlot(F_RV_pressure,V);
legend(hl,{'Fixed','LV Pressure','RV Pressure'});

axisGeom(gca,fontSize);
camlight headlight;
drawnow;
%%
%Get a template with default settings
[febio_spec]=febioStructTemplate;

%febio_spec version
febio_spec.ATTR.version='4.0';

%Module section
febio_spec.Module.ATTR.type='solid';

%Control section
febio_spec.Control.analysis='DYNAMIC';
febio_spec.Control.time_steps=numTimeSteps;
febio_spec.Control.step_size=1/numTimeSteps;
febio_spec.Control.solver.max_refs=max_refs;
febio_spec.Control.solver.qn_method.ATTR.type='Broyden';
febio_spec.Control.solver.qn_method.max_ups=max_ups;
febio_spec.Control.solver.symmetric_stiffness=symmetric_stiffness;
febio_spec.Control.time_stepper.dtmin=dtmin;
febio_spec.Control.time_stepper.dtmax=dtmax;
febio_spec.Control.time_stepper.max_retries=max_retries;
febio_spec.Control.time_stepper.opt_iter=opt_iter;


%Material section
materialName1='Material1';
febio_spec.Material.material{1}.ATTR.name=materialName1;
febio_spec.Material.material{1}.ATTR.type='Holzapfel_Ogden';
febio_spec.Material.material{1}.ATTR.id=1;
febio_spec.Material.material{1}.a=7.75246;
febio_spec.Material.material{1}.b=22.827351;
febio_spec.Material.material{1}.af=2.306745;
febio_spec.Material.material{1}.bf=50.0;
febio_spec.Material.material{1}.as=0.114922;
febio_spec.Material.material{1}.bs=109.120469;
febio_spec.Material.material{1}.afs=0.0;
febio_spec.Material.material{1}.bfs=0.0;
febio_spec.Material.material{1}.asn=0.0;
febio_spec.Material.material{1}.bsn=0.0;
febio_spec.Material.material{1}.anf=0.0;
febio_spec.Material.material{1}.bnf=0.0;
febio_spec.Material.material{1}.k=10000.0;

% Mesh section
% -> Nodes
febio_spec.Mesh.Nodes{1}.ATTR.name='Object1'; %The node set name
febio_spec.Mesh.Nodes{1}.node.ATTR.id=(1:size(V,1))'; %The node id's
febio_spec.Mesh.Nodes{1}.node.VAL=V; %The nodel coordinates

% -> Elements
partName1='Part1';
febio_spec.Mesh.Elements{1}.ATTR.name=partName1; %Name of this part
febio_spec.Mesh.Elements{1}.ATTR.type='tet4'; %Element type
febio_spec.Mesh.Elements{1}.elem.ATTR.id=(1:1:size(E,1))'; %Element id's
febio_spec.Mesh.Elements{1}.elem.VAL=E; %The element matrix


%Mesh section
% -> Nodes
febio_spec.Mesh.Nodes{1}.ATTR.name='nodeSet_all'; %The node set name
febio_spec.Mesh.Nodes{1}.node.ATTR.id=(1:size(V,1))'; %The node id's
febio_spec.Mesh.Nodes{1}.node.VAL=V; %The nodel coordinates

% -> Elements
partName1='Part1';
febio_spec.Mesh.Elements{1}.ATTR.name=partName1; %Name of this part
febio_spec.Mesh.Elements{1}.ATTR.type='tet4'; %Element type
febio_spec.Mesh.Elements{1}.elem.ATTR.id=(1:1:size(E,1))'; %Element id's
febio_spec.Mesh.Elements{1}.elem.VAL=E; %The element matrix


% -> Surfaces
LV_surfaceName='LVPressure';
febio_spec.Mesh.Surface{1}.ATTR.name=LV_surfaceName;
febio_spec.Mesh.Surface{1}.tri3.ATTR.id=(1:1:size(F_LV_pressure,1))';
febio_spec.Mesh.Surface{1}.tri3.VAL=F_LV_pressure;

RV_surfaceName='RVPressure';
febio_spec.Mesh.Surface{2}.ATTR.name=RV_surfaceName;
febio_spec.Mesh.Surface{2}.tri3.ATTR.id=(1:1:size(F_RV_pressure,1))';
febio_spec.Mesh.Surface{2}.tri3.VAL=F_RV_pressure;

% -> NodeSets
nodeSetName1='bcSupportList';
febio_spec.Mesh.NodeSet{1}.ATTR.name=nodeSetName1;
febio_spec.Mesh.NodeSet{1}.VAL=mrow(bcSupportList);

%MeshDomains section
febio_spec.MeshDomains.SolidDomain.ATTR.name=partName1;
febio_spec.MeshDomains.SolidDomain.ATTR.mat=materialName1;

%Boundary condition section
% -> Fix boundary conditions
febio_spec.Boundary.bc{1}.ATTR.name='FixedDisplacement01';
febio_spec.Boundary.bc{1}.ATTR.type='zero displacement';
febio_spec.Boundary.bc{1}.ATTR.node_set=nodeSetName1;
febio_spec.Boundary.bc{1}.x_dof=1;
febio_spec.Boundary.bc{1}.y_dof=1;
febio_spec.Boundary.bc{1}.z_dof=1;


%Loads section
febio_spec.Loads.surface_load{1}.ATTR.type='pressure';
febio_spec.Loads.surface_load{1}.ATTR.surface=LV_surfaceName;
febio_spec.Loads.surface_load{1}.pressure.ATTR.lc=1;
febio_spec.Loads.surface_load{1}.pressure.VAL=6;
febio_spec.Loads.surface_load{1}.symmetric_stiffness=1;

febio_spec.Loads.surface_load{2}.ATTR.type='pressure';
febio_spec.Loads.surface_load{2}.ATTR.surface=RV_surfaceName;
febio_spec.Loads.surface_load{2}.pressure.ATTR.lc=1;
febio_spec.Loads.surface_load{2}.pressure.VAL=6;
febio_spec.Loads.surface_load{2}.symmetric_stiffness=1;


%LoadData section
% -> load_controller
febio_spec.LoadData.load_controller{1}.ATTR.name='LC_1';
febio_spec.LoadData.load_controller{1}.ATTR.id=1;
febio_spec.LoadData.load_controller{1}.ATTR.type='loadcurve';
febio_spec.LoadData.load_controller{1}.interpolate='LINEAR';
%febio_spec.LoadData.load_controller{1}.extend='CONSTANT';
febio_spec.LoadData.load_controller{1}.points.pt.VAL=[0 0; 1 1];



%Output section
% -> log file
febio_spec.Output.logfile.ATTR.file=febioLogFileName;
febio_spec.Output.logfile.node_data{1}.ATTR.file=febioLogFileName_disp;
febio_spec.Output.logfile.node_data{1}.ATTR.data='ux;uy;uz';
febio_spec.Output.logfile.node_data{1}.ATTR.delim=',';

febio_spec.Output.logfile.element_data{1}.ATTR.file=febioLogFileName_stress;
febio_spec.Output.logfile.element_data{1}.ATTR.data='s1';
febio_spec.Output.logfile.element_data{1}.ATTR.delim=',';

febio_spec.Output.plotfile.compression=0;

febioStruct2xml(febio_spec,febioFebFileName);

%% Running the FEBio analysis
% To run the analysis defined by the created FEBio input file the 
% |runMonitorFEBio| function is used. The input for this function is a
% structure defining job settings e.g. the FEBio input file name. The
% optional output runFlag informs the user if the analysis was run
% succesfully. 

febioAnalysis.run_filename=febioFebFileName; %The input file name
febioAnalysis.run_logname=febioLogFileName; %The name for the log file
febioAnalysis.disp_on=1; %Display information on the command window
febioAnalysis.runMode=runMode;

[runFlag]=runMonitorFEBio(febioAnalysis);%START FEBio NOW!!!!!!!!


%%

dataStruct=importFEBio_logfile(fullfile(savePath,febioLogFileName_disp),0,1);

%Access data
N_disp_mat=dataStruct.data; %Displacement
timeVec=dataStruct.time; %Time

%Create deformed coordinate set
V_DEF=N_disp_mat+repmat(V,[1 1 size(N_disp_mat,3)]);

% % Plotting the simulated results using anim8 to visualize and animate deformations
% DN_magnitude=sqrt(sum(N_disp_mat(:,:,end).^2,2)); %Current displacement magnitude
% 
% % Create basic view and store graphics handle to initiate animation
% hf=cFigure; %Open figure
% gtitle([febioFebFileNamePart,': Press play to animate']);
% title('Displacement magnitude [mm]','Interpreter','Latex')
% hp=gpatch(Fb,V_DEF(:,:,end),DN_magnitude,'k',1); %Add graphics object to animate
% hp.Marker='.';
% hp.MarkerSize=markerSize/10;
% hp.FaceColor='interp';
% 
% axisGeom(gca,fontSize);
% colormap(gjet(250)); colorbar;
% caxis([0 max(DN_magnitude)]);
% axis(axisLim(V_DEF)); %Set axis limits statically
% camlight headlight;
% 
% % Set up animation features
% animStruct.Time=timeVec; %The time vector
% for qt=1:1:size(N_disp_mat,3) %Loop over time increments
%     DN_magnitude=sqrt(sum(N_disp_mat(:,:,qt).^2,2)); %Current displacement magnitude
% 
%     %Set entries in animation structure
%     animStruct.Handles{qt}=[hp hp]; %Handles of objects to animate
%     animStruct.Props{qt}={'Vertices','CData'}; %Properties of objects to animate
%     animStruct.Set{qt}={V_DEF(:,:,qt),DN_magnitude}; %Property values for to set in order to animate
% end
% anim8(hf,animStruct); %Initiate animation feature
% drawnow;

%%

V_FinalDeformed = V_DEF(:,:,end);

% % Visualize interior point
cFigure; hold on;
gpatch(Fb,V,'w','k',0.75);
axisGeom; camlight headlight;
drawnow;

gpatch(Fb,V_FinalDeformed,'w','r',0.25);

% [F_new,V_new]=mergeVertices(Fb,V_FinalDeformed);
% V_inner_new = getInnerPoint(F_new,V_new,[],[],1);
% 
% tetVolume_new=tetVolMeanEst(F_new,V_new); %Volume for regular tets
% 
% tetGenStruct.stringOpt='-pq1.2AaY';
% tetGenStruct.Faces=F_new;
% tetGenStruct.Nodes=V_new;
% tetGenStruct.holePoints=[];
% tetGenStruct.faceBoundaryMarker=C; %Face boundary markers
% tetGenStruct.regionPoints=V_inner_new; %region points
% tetGenStruct.regionA=tetVolume_new;
% 
% [meshOutput_new]=runTetGen(tetGenStruct); %Run tetGen
% 
% % Access elements, nodes, and boundary faces
% E=meshOutput_new.elements;
% V_new=meshOutput_new.nodes;
% Fb_new=meshOutput_new.facesBoundary;
% Cb_new=meshOutput_new.boundaryMarker;
% CE=meshOutput_new.elementMaterialID;

% % % % % % meshView(meshOutput);

cFigure; hold on;
gpatch(Fb,V_FinalDeformed,Cb,'k',0.5);
axisGeom; camlight headlight;
drawnow;


%% Animations
if runFlag==1 %i.e. a succesful run

    dataStruct=importFEBio_logfile(fullfile(savePath,febioLogFileName_disp),0,1);

    %Access data
    N_disp_mat=dataStruct.data; %Displacement
    timeVec=dataStruct.time; %Time

    %Create deformed coordinate set
    V_DEF=N_disp_mat+repmat(V,[1 1 size(N_disp_mat,3)]);
    % Plotting the simulated results using anim8 to visualize and animate deformations
    DN_magnitude=sqrt(sum(N_disp_mat(:,:,end).^2,2)); %Current displacement magnitude

    % Create basic view and store graphics handle to initiate animation
    hf=cFigure; %Open figure
    gtitle([febioFebFileNamePart,': Press play to animate']);
    title('Displacement magnitude [mm]','Interpreter','Latex')
    hp=gpatch(Fb,V_DEF(:,:,end),DN_magnitude,'k',1); %Add graphics object to animate
    hp.Marker='.';
    hp.MarkerSize=markerSize/10;
    hp.FaceColor='interp';

    axisGeom(gca,fontSize);
    colormap(gjet(250)); colorbar;
    caxis([0 max(DN_magnitude)]);
    axis(axisLim(V_DEF)); %Set axis limits statically
    camlight headlight;

    % Set up animation features
    animStruct.Time=timeVec; %The time vector
    for qt=1:1:size(N_disp_mat,3) %Loop over time increments
        DN_magnitude=sqrt(sum(N_disp_mat(:,:,qt).^2,2)); %Current displacement magnitude

        %Set entries in animation structure
        animStruct.Handles{qt}=[hp hp]; %Handles of objects to animate
        animStruct.Props{qt}={'Vertices','CData'}; %Properties of objects to animate
        animStruct.Set{qt}={V_DEF(:,:,qt),DN_magnitude}; %Property values for to set in order to animate
    end
    anim8(hf,animStruct); %Initiate animation feature
    drawnow;
end


