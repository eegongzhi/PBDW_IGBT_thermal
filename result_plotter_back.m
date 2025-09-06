function s = result_plotter_back(elementCoord,dataForPlot)
% input: dataForPlot; elementCoord
% output: surface(X,Y,Z)

    grid_x = (-103/2 + 5 + 4):1:-(-103/2 + 5 + 4);
    grid_y = (-78/2 + 5 + 4):1:-(-78/2 + 5 + 4);

    [X, Y] = meshgrid(grid_x, grid_y);
    Z = griddata(elementCoord(:,1),elementCoord(:,2),dataForPlot,...
            X,Y,'linear');
    s = surface(X,Y,Z);
    s.EdgeColor='none';
    colorbar

end

