function [xg,top,fnx,fny,fnz,nelx,nely,nelz,dx,dy,dz] = symmetry(xg,ls,nelx,nely,nelz,dx,dy,dz,ngrid,symm,dir)
if strcmp(symm,'x-y')
   if strcmp(dir,'left')
      [sx,sy,sz] = size(xg); xg1 = zeros(sx,sy,2*sz); xg_mir = flip(xg,3); xg1(:,:,1:sz) = xg_mir; xg1(:,:,sz+1:end) = xg;
      xg = xg1; 
      nelz = 2*nelz; dz = 2*dz; [fnx,fny,fnz] = meshgrid(0:1/ngrid:nelx+5,0:1/ngrid:nely+5,0:1/ngrid:nelz+5);
      fnx(:,:,size(xg,3)+1:end) = []; fny(:,:,size(xg,3)+1:end) = []; fnz(:,:,size(xg,3)+1:end) = [];
      fnx(:,size(xg,2)+1:end,:) = []; fny(:,size(xg,2)+1:end,:) = []; fnz(:,size(xg,2)+1:end,:) = []; 
      fnx(size(xg,1)+1:end,:,:) = []; fny(size(xg,1)+1:end,:,:) = []; fnz(size(xg,1)+1:end,:,:) = [];
   else
      [sx,sy,sz] = size(xg); xg1 = zeros(sx,sy,2*sz); xg_mir = flip(xg,3); xg1(:,:,1:sz) = xg; xg1(:,:,sz+1:end) = xg_mir;
      xg = xg1; 
      nelz = 2*nelz; dz = 2*dz; [fnx,fny,fnz] = meshgrid(0:1/ngrid:nelx+5,0:1/ngrid:nely+5,0:1/ngrid:nelz+5);
      fnx(:,:,size(xg,3)+1:end) = []; fny(:,:,size(xg,3)+1:end) = []; fnz(:,:,size(xg,3)+1:end) = [];
      fnx(:,size(xg,2)+1:end,:) = []; fny(:,size(xg,2)+1:end,:) = []; fnz(:,size(xg,2)+1:end,:) = []; 
      fnx(size(xg,1)+1:end,:,:) = []; fny(size(xg,1)+1:end,:,:) = []; fnz(size(xg,1)+1:end,:,:) = [];
   end 
elseif strcmp(symm,'y-z')
   if strcmp(dir,'left')
      [sx,sy,sz] = size(xg); xg1 = zeros(sx,2*sy,sz); xg_mir = flip(xg,2); xg1(:,1:sy,:) = xg_mir; xg1(:,sy+1:end,:) = xg;
      xg = xg1; 
      nelx = 2*nelx; dx = 2*dx; [fnx,fny,fnz] = meshgrid(0:1/ngrid:nelx+5,0:1/ngrid:nely+5,0:1/ngrid:nelz+5);
      fnx(:,:,size(xg,3)+1:end) = []; fny(:,:,size(xg,3)+1:end) = []; fnz(:,:,size(xg,3)+1:end) = [];
      fnx(:,size(xg,2)+1:end,:) = []; fny(:,size(xg,2)+1:end,:) = []; fnz(:,size(xg,2)+1:end,:) = []; 
      fnx(size(xg,1)+1:end,:,:) = []; fny(size(xg,1)+1:end,:,:) = []; fnz(size(xg,1)+1:end,:,:) = [];
   else 
      [sx,sy,sz] = size(xg); xg1 = zeros(sx,2*sy,sz); xg_mir = flip(xg,2); xg1(:,1:sy,:) = xg; xg1(:,sy+1:end,:) = xg_mir;
      xg = xg1; 
      nelx = 2*nelx; dx = 2*dx; [fnx,fny,fnz] = meshgrid(0:1/ngrid:nelx+5,0:1/ngrid:nely+5,0:1/ngrid:nelz+5);
      fnx(:,:,size(xg,3)+1:end) = []; fny(:,:,size(xg,3)+1:end) = []; fnz(:,:,size(xg,3)+1:end) = [];
      fnx(:,size(xg,2)+1:end,:) = []; fny(:,size(xg,2)+1:end,:) = []; fnz(:,size(xg,2)+1:end,:) = []; 
      fnx(size(xg,1)+1:end,:,:) = []; fny(size(xg,1)+1:end,:,:) = []; fnz(size(xg,1)+1:end,:,:) = [];
   end 
elseif strcmp(symm,'z-x')
   if strcmp(dir,'right')
      [sx,sy,sz] = size(xg); xg1 = zeros(2*sx,sy,sz); xg_mir = flip(xg,1); xg1(1:sx,:,:) = xg; xg1(sx+1:end,:,:) = xg_mir;
      xg = xg1; 
      nely = 2*nely; dy = 2*dy; [fnx,fny,fnz] = meshgrid(0:1/ngrid:nelx+5,0:1/ngrid:nely+5,0:1/ngrid:nelz+5);
      fnx(:,:,size(xg,3)+1:end) = []; fny(:,:,size(xg,3)+1:end) = []; fnz(:,:,size(xg,3)+1:end) = [];
      fnx(:,size(xg,2)+1:end,:) = []; fny(:,size(xg,2)+1:end,:) = []; fnz(:,size(xg,2)+1:end,:) = []; 
      fnx(size(xg,1)+1:end,:,:) = []; fny(size(xg,1)+1:end,:,:) = []; fnz(size(xg,1)+1:end,:,:) = [];
   else 
      [sx,sy,sz] = size(xg); xg1 = zeros(2*sx,sy,sz); xg_mir = flip(xg,1); xg1(1:sx,:,:) = xg_mir; xg1(sx+1:end,:,:) = xg;
      xg = xg1; 
      nely = 2*nely; dy = 2*dy; [fnx,fny,fnz] = meshgrid(0:1/ngrid:nelx+5,0:1/ngrid:nely+5,0:1/ngrid:nelz+5);
      fnx(:,:,size(xg,3)+1:end) = []; fny(:,:,size(xg,3)+1:end) = []; fnz(:,:,size(xg,3)+1:end) = [];
      fnx(:,size(xg,2)+1:end,:) = []; fny(:,size(xg,2)+1:end,:) = []; fnz(:,size(xg,2)+1:end,:) = []; 
      fnx(size(xg,1)+1:end,:,:) = []; fny(size(xg,1)+1:end,:,:) = []; fnz(size(xg,1)+1:end,:,:) = [];
   end
end
top = xg-ls;