% The aim of the offline stage is to form matrices and vectors for the
% matrix equation.
clear
close all

%% import the spaces
load('./results/ZN_Q.mat');
load('./results/ZN_R.mat');
load('./results/UM_R.mat');
load('./results/UM_Q.mat');
load('./results/D_final.mat');

load('./results/sensorLocation.mat');

elementCoord_back = importdata('./data/coord_back.txt');
load('./data/surface_area_back.mat');
elementCoord_copper = importdata('./data/coord_copper.txt');
load('./data/surface_area_copper.mat');
elementCoord = [elementCoord_back; elementCoord_copper];

D = D_final;
clear D_final

%% control parameters
% UseArbitrary = 0;
% UseNoise = 0;
% SNR = 0.0;
rm = 1e-5;   % here we assume point-wise measurement

true_DOF_back = length(surface_area_back);
true_DOF_copper = length(surface_area_copper);

xi = 0;
M = size(UM_Q,2);
N = size(ZN_Q,2);

% generate L_eta
L_eta = zeros(M, M); 
for idm = 1:M
    for idm_prime = 1:M
        L_eta(idm, idm_prime) = L_obs(elementCoord, ...
                UM_Q(:,idm_prime),sensorLocation(idm,:), rm, true_DOF_back,...
                surface_area_back, surface_area_copper);
    end
end

% generate L_z
L_z = zeros(M, N);
for idm = 1 : M
    for idn = 1 : N
        L_z(idm, idn) = L_obs(elementCoord, ZN_Q(:,idn),...
                sensorLocation(idm,:), rm, true_DOF_back,...
                surface_area_back, surface_area_copper);
    end
end

K1 = xi * M * eye(M) + L_eta * L_eta'; 
K = [K1, L_z; L_z', zeros(N,N)];

%% online stage
%%% load the test set
err = zeros(1,size(D, 2));

% Large loop: loop all snapshots
for count = 1:size(D,2)
    groundTruth = D(:,count);

    % loop all sensors
    %%%  hidden codes  %%%
    %%                                 %%
    %%%  hidden codes  %%%

    RHS = [lo; zeros(N,1)];

    %%% solve the matrix equation
    sol = pinv(K) * RHS;
    eta_til = sol(1 : M);
    eta = L_eta' * eta_til;
    z = sol(M+1 : end);
    
    %%% re-construct
    u_ZN = zeros(size(ZN_Q,1),1);
    for idx = 1:N
        u_ZN = u_ZN + ZN_Q(:,idx) * z(idx);
    end
    
    u_UM = zeros(size(UM_Q,1),1);
    for idx = 1:M
        u_UM = u_UM + UM_Q(:,idx) * eta(idx);
    end
    
    % u_reconstruct = project(ZN_Q)*u_ZN + project(UM_Q)*u_UM;
    u_reconstruct = u_ZN + u_UM;
      
    %%% evaluate the error
    dif = groundTruth - u_reconstruct;
    
    err(count) = sqrt(innerProduct(dif, dif, surface_area_back, surface_area_copper)...
        /innerProduct(groundTruth,...
        groundTruth, surface_area_back, surface_area_copper));

    count %#ok<NOPTS> 
end

%%% error computation
mean_err = mean(err)  %#ok<NOPTS> 


%% visualize the reconstructed results: back 
figure
result_plotter_back(elementCoord(1:true_DOF_back,:),u_reconstruct(1:true_DOF_back));
colorbar

%%% plot the difference
difference = groundTruth - u_reconstruct;
figure
result_plotter_back(elementCoord(1:true_DOF_back,:),difference(1:true_DOF_back));
colorbar

%% visualize the reconstructed results: copper
figure
result_plotter_copper(elementCoord(true_DOF_back+1:end,:),...
    u_reconstruct(true_DOF_back+1:end));
colorbar

figure
result_plotter_back(elementCoord(true_DOF_back+1:end,:),difference(true_DOF_back+1:end));
colorbar
