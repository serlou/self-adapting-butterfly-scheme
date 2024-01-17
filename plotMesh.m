function plotMesh(X,Y,Z)
    s = surf(X,Y,Z);
    s.EdgeColor = 'none';
    s.FaceColor = 'Interp';
    colormap parula;
%     colormap gray;
    axis tight;
    axis square; 
    axis off;
    view(3);
    axis image;
    light('Position',[-0.1476 -1.4045 -0.0740]);
end