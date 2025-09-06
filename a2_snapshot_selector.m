close all
clear

%% initialize
load('./results/D_building.mat')
D = D_building;
clear D_building

load('./results/surface_area_back.mat')
load('./results/surface_area_copper.mat')
ele_S_back = surface_area_back;  % unit: m^2
ele_S_copper = surface_area_copper;  % unit: m^2

N = 6; % The total number of selected snapshots 100

%% 1st snapshot
ZN= D(:,1); 
selected = ones(size(D,2),1);
[ZN_Q, ZN_R] = gson(ZN, ele_S_back, ele_S_copper);  % ZN = ZN_Q * ZN_R
% [ZN_Q, ZN_R] = qr(ZN,'econ');
jRecorder = 1;
selected(jRecorder) = 0;

%% loop to generate the Zn space
for i=2:N % loop until we have selected N snapshots
    % loop over the NOT-SELECTED snapshots
    norm_perpendicular = zeros(1,size(D,2));
    for idx_snap = 1:size(D,2)
        if selected(idx_snap) == 0
            continue;
        else
           %%%  hidden codes  %%%
           %%                                 %%
           %%%  hidden codes  %%%
        end
    end
    % find the maximum norm_perpendicular
    [~,recorder] = max(norm_perpendicular);
    selected(recorder) = 0;
    jRecorder = [jRecorder, recorder];

    % add to ZN and the jRecorder
    ZN = ZN_Q * ZN_R;
    ZN = [ZN, D(:,recorder)];

    [ZN_Q, ZN_R] = gson(ZN, ele_S_back, ele_S_copper);
%     [ZN_Q, ZN_R] = qr(ZN,'econ');
    i %#ok<NOPTS> 
end

save('./results/ZN_Q.mat', 'ZN_Q')  % we save the generated Z_N space
save('./results/ZN_R.mat', 'ZN_R')

%% visialize the jRecorder
distance = linspace(1, size(D,2), size(D,2));
scatter(jRecorder, ones(N,1));
grid on

% We also save the distance sequence number
chosen_snapshots = zeros(2,size(D,2));
chosen_snapshots(1,:) = distance;
chosen_snapshots(2,jRecorder) = 1;
save('./results/chosen_snapshots.mat', 'chosen_snapshots')

%% TEST
% innerProduct(ZN_Q(:,1),ZN_Q(:,7),alpha) % innerprodct expected to be zero
% innerProduct(ZN_Q(:,1),ZN_Q(:,1),alpha) % norm expected to be unity
% imagesc(reshape(ZN_Q(:,1), 200,200))
figure  % B field
true_DOF_back = length(ele_S_back);
n_basis = 5;
dataForPlot = ZN_Q(1:true_DOF_back,n_basis);
load('./results/elementCoord.mat')
result_plotter_back(elementCoord(1:true_DOF_back,:),dataForPlot);

figure
dataForPlot = ZN_Q(true_DOF_back+1:end, n_basis);
result_plotter_copper(elementCoord(true_DOF_back+1:end,:),dataForPlot);
