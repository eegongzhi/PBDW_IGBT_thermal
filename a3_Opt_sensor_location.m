% input: the space ZN
% output: the sensor location coordinate, and the generated space UM

close all
clear

%% control parameters
rm = 1e-5;  % assume: point-wise measurement, i.e. small sensors  
MMax = 20;   

elementCoord_back = importdata('./data/coord_back.txt');
load('./data/surface_area_back.mat');
elementCoord_copper = importdata('./data/coord_copper.txt');
load('./data/surface_area_copper.mat');
elementCoord = [elementCoord_back; elementCoord_copper];

constraint = ones(length(surface_area_back)+length(surface_area_copper), 1);
true_DOF_back = length(surface_area_back);
true_DOF_copper = length(surface_area_copper);

%% Initialization
% recover the background space Zn (snapshots)
load('./results/ZN_Q.mat');
ZN_vector = ZN_Q;
clear ZN_Q
N = size(ZN_vector, 2); % number of background space bases

sensorLocation = [];
UM = [];

%% M=1
% first, find the location
Z_basis = ZN_vector(:,1);
node_idx = ...
    find(abs(Z_basis) == max(abs(Z_basis .* constraint)));
constraint(node_idx) = 0; % no longsonger available for sensor position choice
% if node_idx > true_DOF
%     constraint(node_idx - true_DOF) = 0; 
% else
%     constraint(node_idx + true_DOF) = 0;
% end
y_m = elementCoord(node_idx,2); 
x_m = elementCoord(node_idx,1); 

% find the true position (as coordinate)

sensorLocation = [sensorLocation; x_m, y_m, node_idx]; % add into the location set

% second, generate the 1st update basis. 
UM = [UM, phi_m(elementCoord, x_m, y_m, rm, node_idx, true_DOF_back,...
    surface_area_back, surface_area_copper)];
[UM_Q, UM_R] = gson(UM, surface_area_back, surface_area_copper);
% [UM_Q, UM_R] = qr(UM,'econ');

%% M>1. Currently, we do not consider the approximation stage
M = 1;
tol = 1;
beta = zeros(1,MMax-1);
beta_new = 0;

while M < MMax
    %%%  hidden codes  %%%
    %%                                 %%
    %%%  hidden codes  %%%


    % compute new sensor location
    %%% construct Lm
    L_m = zeros(M,1);
    for idm = 1:M
        L_m(idm, 1) = L_obs(elementCoord,z_min,...
            sensorLocation(idm,:), rm, true_DOF_back, surface_area_back, surface_area_copper);
    end

    %%% calculate Im
    term1 = pinv(L_eta)*L_m;
    Im = zeros(size(ZN_vector,1),1);
    for idm = 1:M
        Im = Im + term1(idm,1) * UM_Q(:,idm);
    end
    dif = abs(z_min - Im);
    [node_idx] = find(dif == max(dif .* constraint)); % indices
    node_idx = node_idx(1);
    y_new = elementCoord(node_idx,2);  % true position: 
    x_new = elementCoord(node_idx,1);


    isInArray = ismember([x_new,y_new,node_idx], sensorLocation,'rows');
    if (isInArray == 0) %&& beta_new < 0.5
        sensorLocation = [sensorLocation; x_new, y_new, node_idx];
        constraint(node_idx) = 0;

        %%% supplement update space
        UM = UM_Q * UM_R;
        phi_new = phi_m(elementCoord, x_new, y_new, rm, node_idx, true_DOF_back,...
    surface_area_back, surface_area_copper);
        UM = [UM, phi_new];
        [UM_Q, UM_R] = gson(UM, surface_area_back, surface_area_copper);
        M = M + 1;
        M %#ok<NOPTS> 
        beta_new %#ok<NOPTS> 
    else
        break
    end
end

%% visualization: here, we plot the selected sensors (not related to the real sensors)
% distinguish sensors
sensorLocation_back = [];
sensorLocation_copper = [];
for idx = 1:size(sensorLocation, 1)
    if sensorLocation(idx, 3) < true_DOF_back+1  
        sensorLocation_back = [sensorLocation_back;sensorLocation(idx,:)];
    else
        sensorLocation_copper = [sensorLocation_copper;sensorLocation(idx,:)];
    end
end

% plot the back layer
figure
x = elementCoord(1:true_DOF_back,1);
y = elementCoord(1:true_DOF_back,2);
grid_x = (-103/2 + 5 + 4):1:-(-103/2 + 5 + 4);
grid_y = (-78/2 + 5 + 4):1:-(-78/2 + 5 + 4);
[X, Y] = meshgrid(grid_x, grid_y);
Z = griddata(x,y,ZN_vector(1:true_DOF_back,1),...
        X,Y,'linear');
s = surface(X,Y,Z);
s.EdgeColor='none';
pointsX = sensorLocation_back(:,1);
pointsY = sensorLocation_back(:,2);
pointsZ = 1000*ones(length(pointsY),1);
hold on
scatter3(pointsX, pointsY, pointsZ,'filled', 'r')


% plot the copper layer
figure
x = elementCoord(true_DOF_back+1:end,1);
y = elementCoord(true_DOF_back+1:end,2);
grid_x = (-56/2 + 3):1:-(-56/2 + 3);
grid_y = (-58/2 + 3):1:-(-58/2 + 3);
[X, Y] = meshgrid(grid_x, grid_y);
Z = griddata(x,y,ZN_vector(true_DOF_back+1:end,1),...
        X,Y,'linear');
s = surface(X,Y,Z);
s.EdgeColor='none';
pointsX = sensorLocation_copper(:,1);
pointsY = sensorLocation_copper(:,2);
pointsZ = 1000*ones(length(pointsY),1);
hold on
scatter3(pointsX, pointsY, pointsZ,'filled', 'r')


%% plot beta
beta = [0, beta];

figure
plot(abs(beta),'-o')
grid on
xlabel('Number of Sensors')
ylabel('beta')

save('./results/sensorLocation.mat', 'sensorLocation')
save('./results/UM_Q.mat', 'UM_Q')
save('./results/UM_R.mat', 'UM_R')
save('./results/beta.mat', 'beta')


%% define the update generator
function phi_results = phi_m(elementCoord, x_m, y_m, rm, node_idx, true_DOF_B,...
    dS_B, dS_T)
    % This function returns 1D, not-orthornormal, phi
    if node_idx < true_DOF_B+1   % sensor for B field
        Rul = zeros(true_DOF_B,1);
        elementCoord_B = elementCoord(1:true_DOF_B,:);

        for idx = 1:length(Rul)
            disVec = [elementCoord_B(idx,1)-x_m, elementCoord_B(idx,2)-y_m];
            dis = disVec(1)^2 + disVec(2)^2;
            Rul(idx) = exp(-dis/(2*rm^2));
        end
        % calculate the constant
        C = 0;
        for idx = 1:length(Rul)
            disVec = [elementCoord_B(idx,1)-x_m, elementCoord_B(idx,2)-y_m];
            dis = disVec(1)^2 + disVec(2)^2;
            C = C + dS_B(idx) * 1* exp(-dis/(2*rm^2));
        end
        phi_results = [Rul*(1/C); zeros(length(elementCoord)-true_DOF_B,1)];

    else
        Rul = zeros(length(elementCoord)-true_DOF_B,1);
        elementCoord_T = elementCoord(true_DOF_B+1:end,:);

        for idx = 1:length(Rul)
            disVec = [elementCoord_T(idx,1)-x_m, elementCoord_T(idx,2)-y_m];
            dis = disVec(1)^2 + disVec(2)^2;
            Rul(idx) = exp(-dis/(2*rm^2));
        end
        % calculate the constant
        C = 0;
        for idx = 1:length(Rul)
            disVec = [elementCoord_T(idx,1)-x_m, elementCoord_T(idx,2)-y_m];
            dis = disVec(1)^2 + disVec(2)^2;
            C = C + dS_T(idx) * 1* exp(-dis/(2*rm^2));
        end
        phi_results = [zeros(true_DOF_B,1); Rul*(1/C)];
    end

end



