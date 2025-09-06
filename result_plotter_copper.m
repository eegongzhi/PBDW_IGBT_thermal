function s = result_plotter_copper(elementCoord,dataForPlot)
% input: dataForPlot; elementCoord
% output: surface(X,Y,Z)

    grid_x = (-56/2 + 3):1:-(-56/2 + 3);
    grid_y = (-58/2 + 3):1:-(-58/2 + 3);

    [X, Y] = meshgrid(grid_x, grid_y);
    Z = griddata(elementCoord(:,1),elementCoord(:,2),dataForPlot,...
            X,Y,'linear');
    s = surface(X,Y,Z);
    s.EdgeColor='none';


end

