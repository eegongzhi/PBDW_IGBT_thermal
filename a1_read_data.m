close all
clear

%% load and normalization
% import the element information
elementCoord_back = importdata('./data/coord_back.txt');
load('./data/surface_area_back.mat');
elementCoord_copper = importdata('./data/coord_copper.txt');
load('./data/surface_area_copper.mat');

% import the snapshots
load('./data/T_back.mat'); 
load('./data/T_copper.mat');

D_full = [T_back; T_copper];
elementCoord = [elementCoord_back; elementCoord_copper];  % unit: m

%% Partition into 9:1
total_idx = 1:size(D_full, 2);
data_for_final = round(size(D_full, 2) * 0.1);
data_for_building = size(D_full, 2) - data_for_final;
idx_building = randperm(size(D_full, 2), data_for_building);
idx_final = setdiff(total_idx, idx_building);

D_building = D_full(:, idx_building);
D_final = D_full(:, idx_final);


%% clear and save
save('./results/D_full.mat', 'D_full')
save('./results/D_building.mat', 'D_building')
save('./results/D_final.mat', 'D_final')
% save the combined "elementCoord"
save('./results/elementCoord.mat', 'elementCoord')
save('./results/surface_area_back.mat','surface_area_back')
save('./results/surface_area_copper.mat','surface_area_copper')