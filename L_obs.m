%% define the observation function
%%% Todo: add alpha here to see if the results are satisfactory
function L_obs = L_obs(elementCoord, u, center, rm, true_DOF_B, dS_B, dS_T)
    % In this project, u is 1D
    L_obs = 0;
    
    if center(3) < true_DOF_B+1
        u = u(1:true_DOF_B);
        elementCoord_B = elementCoord(1:true_DOF_B,:);
        
        for idx = 1:length(u)
            disVec = [elementCoord_B(idx,1)-center(1),...
                elementCoord_B(idx,2)-center(2)];
            dis = disVec(1)^2 + disVec(2)^2;
            L_obs = L_obs + dS_B(idx) * u(idx) * exp(-dis/(2*rm^2));
        end
    
        % calculate the constant
        C = 0;
        for idx = 1:length(u)
            disVec = [elementCoord_B(idx,1)-center(1),...
                elementCoord_B(idx,2)-center(2)];
            dis = disVec(1)^2 + disVec(2)^2;
            C = C + dS_B(idx) * 1 * exp(-dis/(2*rm^2));
        end
        L_obs = L_obs * (1/C);

    else
        u = u(true_DOF_B+1:end);
        elementCoord_T = elementCoord(true_DOF_B+1:end,:);

        for idx = 1:length(u)
            disVec = [elementCoord_T(idx,1)-center(1),...
                elementCoord_T(idx,2)-center(2)];
            dis = disVec(1)^2 + disVec(2)^2;
            L_obs = L_obs + dS_T(idx) * u(idx) * exp(-dis/(2*rm^2));
        end
    
        % calculate the constant
        C = 0;
        for idx = 1:length(u)
            disVec = [elementCoord_T(idx,1)-center(1),...
                elementCoord_T(idx,2)-center(2)];
            dis = disVec(1)^2 + disVec(2)^2;
            C = C + dS_T(idx) * 1 * exp(-dis/(2*rm^2));
        end
        L_obs = L_obs * (1/C);

    end
end

