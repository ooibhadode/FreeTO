function display_3Dsmooth(xg,top)
    [nely,nelx,nelz] = size(xg);
    cla, hold on, view(30,30), rotate3d on, axis equal, axis([0 nelx 0 nely 0 nelz]), box
    set(gca,'YDir','reverse','ZDir','reverse','ZtickLabel',flipud(get(gca,'Ztick')'));
    [X,Y,Z] = meshgrid(1:nelx,1:nely,1:nelz);
    fcolor=[0 1 1];
    patch(isocaps(X,Y,Z,top,0),'FaceColor',fcolor,'EdgeColor','none');
    p = patch(isosurface(X,Y,Z,top,0),'FaceColor',fcolor,'EdgeColor','none');
    camlight('headlight')
    axis off; 
    box off; 
    set(gcf, 'color', [1 1 1]) %Background
    drawnow
end