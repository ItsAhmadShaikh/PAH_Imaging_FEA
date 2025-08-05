clear all
clc
close all

gm = fegeometry("C:\Users\ahmad\OneDrive - UC San Diego\DVJ_Lab\cMRI\Code_GitHub\DTI_testing\Meshes\N24_cut_full.stl");
pdegplot(gm,FaceLabels="on",FaceAlpha=0.3)

fiber = load("C:\Users\ahmad\OneDrive - UC San Diego\DVJ_Lab\cMRI\Code_GitHub\DTI_testing\Meshes\N24_fiber.mat");

x = 0:150;
y = 0:200;
z = 0:150;


