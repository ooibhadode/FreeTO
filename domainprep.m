function [edofMat1,edofMat,ndn,nd,edofMatn,nnele,ele,n_vec,vol1,MusD]...
    = domainprep(oute,outeM,vol,nelx,nely,nelz,nele)

%% Node (edofMat1) and element (edofMat) connectivity matrices
nodenr = reshape(1:(1+nelx)*(1+nely)*(1+nelz),(1+nely),(1+nelx),(1+nelz));                   % node numbers in domain
nodenrs1 = nodenr-1; nodenrs1(end,:,:) = []; nodenrs1(:,end,:) = []; nodenrs1(:,:,end) = [];
edofMat1 = nodenrs1(:) + [2 2+(nely+1) 2+nely 1 (nelx+1)*(nely+1)+2 (nelx+1)*(nely+1)+2+(nely+1) ...
    (nelx+1)*(nely+1)+2+nely (nelx+1)*(nely+1)+1];                          % connectivity matrix for node numbers
nodenrs = 3*(nodenr-1); nodenrs(end,:,:) = []; nodenrs(:,end,:) = []; nodenrs(:,:,end) = [];
% element connectivity matrix
edofMat = nodenrs(:) + [4 5 6 4+3*(nely+1) 5+3*(nely+1) 6+3*(nely+1) 1+3*(nely+1) 2+3*(nely+1) 3+3*(nely+1) ...
    1 2 3 4+3*(nelx+1)*(nely+1) 5+3*(nelx+1)*(nely+1) 6+3*(nelx+1)*(nely+1) 4+3*(nelx+1)*(nely+1)+3*(nely+1)... 
    5+3*(nelx+1)*(nely+1)+3*(nely+1) 6+3*(nelx+1)*(nely+1)+3*(nely+1) 1+3*(nelx+1)*(nely+1)+3*(nely+1)...
     2+3*(nelx+1)*(nely+1)+3*(nely+1) 3+3*(nelx+1)*(nely+1)+3*(nely+1) 1+3*(nelx+1)*(nely+1)... 
     2+3*(nelx+1)*(nely+1) 3+3*(nelx+1)*(nely+1)];

%% List of important parameters needed for optimization
ndn = find(oute(:) == 0);                   % vector of passive domain elements
nd = length(ndn);                           % number of passive domain elements
edofMatn = edofMat; edofMatn(ndn,:) = [];   % connectivity matrix of active domain elements only
nnele = length(edofMatn(:,1));              % number of active domain elements
ele = setdiff(1:nele,ndn);                  % active domain element numbers from total domain element numbers
n_vec = unique(edofMatn(:));                % degrees of freedom for active domain elements only
vol1 = vol*(nnele/nele);                    % volume fraction based on active domains only
outeM = outeM(:); MusD = find(outeM == 1);
