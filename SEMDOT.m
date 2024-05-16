function [nelx,nely,nelz,ngrid,stp1,xg,top,lss,fnx,fny,fnz,del_x,del_y,del_z,S] = SEMDOT(MeshControl,domain,fixed,fixed_x,fixed_y,fixed_z,force_1,...
    force_2,force_3,force_4,force_5,force_6,force_7,force_8,force_9,force_10,keep_domain,...
    keep_BC,keep_BCxyz,vol,E00,nu,rmin,Fmagx,Fmagy,Fmagz,loadtype)
%% Arbitrary Domains
[~,oute,outeM,sup_all,sup_x,sup_y,sup_z,Fn,nelx,nely,nelz,nele,ndof,stp1,aa,del_x,del_y,del_z] = geomeshini(MeshControl,domain,fixed,fixed_x,...
    fixed_y,fixed_z,force_1,force_2,force_3,force_4,force_5,force_6,force_7,force_8,force_9,force_10,keep_domain,keep_BC,keep_BCxyz);
% node and DOF connectivity matrix
[~,~,~,~,edofMatn,nnele,ele,n_vec,~,MusD]...
    = domainprep(oute,outeM,vol,nelx,nely,nelz,nele);
% force nodes
[F,nf] = forcevec(Fn,ndof,Fmagx,Fmagy,Fmagz,loadtype);
% fixed DOFs
fixeddof = supportDOFs(sup_all,sup_x,sup_y,sup_z);
%% USER-DEFINED LOOP PARAMETERS
maxloop = 500;    % Maximum number of iterations
tolx = 0.003;      % Terminarion criterion
%% USER-DEFINED MATERIAL PROPERTIES
E0 = 1*E00;           % Young's modulus of solid material
Emin = 0.001;      % Young's modulus of void-like material
tol = 1; tol_thresh = 3e-3;
%% USER-DEFINED GRID POINTS
ngrid=4; rnmin=1;
%% INITIALIZE HEAVISIDE REGULARIZATION PARAMETER
beta=0.1; ER=0.05;
%% ELEMENTAL NODES AND COORDINATES
[nodex,nodey,nodez] = meshgrid(0:nelx,0:nely,0:nelz);
[fnx,fny,fnz] = meshgrid(0:1/ngrid:nelx,0:1/ngrid:nely,0:1/ngrid:nelz);
%% USER-DEFINED LOAD DOFs
%% Load domain 
U = zeros(ndof,nf);
freedofs = setdiff(n_vec,fixeddof);                                         % make a change here
KE = (aa/0.5)*lk_H8(nu);
iK = reshape(kron(edofMatn,ones(24,1))',24*24*nnele,1);
jK = reshape(kron(edofMatn,ones(1,24))',24*24*nnele,1);
%% PREPARE FILTER
[H,Hs]=HHs3D(nelx,nely,nelz,rmin,ele,nele);
%% PREPARE FILTER FOR NODELS
[Hn,Hns]=HnHns3D(nelx,nely,nelz,rnmin);
%% INITIALIZE ITERATION
vx = repmat(vol,[nely,nelx,nelz]); vx = vx(:); vx = vx(ele);                % make a change here
vxPhys = vx; 
loop = 0; 
change = 1;
% %% INITIALIZE MMA OPTIMIZER
% m     = 1;                % The number of general constraints.
% n     =  nnele;  % The number of design variables x_j.                      % make a change here
% xmin  = zeros(n,1);       % Column vector with the lower bounds for the variables x_j.
% xmax  = ones(n,1);        % Column vector with the upper bounds for the variables x_j.
% xold1 = vx(:);             % xval, one iteration ago (provided that iter>1).
% xold2 = vx(:);             % xval, two iterations ago (provided that iter>2).
% low   = ones(n,1);        % Column vector with the lower asymptotes from the previous iteration (provided that iter>1).
% upp   = ones(n,1);        % Column vector with the upper asymptotes from the previous iteration (provided that iter>1).
% a0    = 1;                % The constants a_0 in the term a_0*z.
% a     = zeros(m,1);       % Column vector with the constants a_i in the terms a_i*z.
% c_MMA = 10000*ones(m,1);  % Column vector with the constants c_i in the terms c_i*y_i.
% d     = zeros(m,1);       % Column vector with the constants d_i in the terms 0.5*d_i*(y_i)^2.
%% START ITERATION
while (change > tolx && tol>tol_thresh) && loop < maxloop 
%% UPDATE ITERATION
    loop = loop+1;
%% FE-ANALYSIS
    sK = KE(:)*(vxPhys(:)'*E0+(1-vxPhys(:))'*(Emin*E0));                    % make a change here
    K = sparse(iK(:),jK(:),sK(:)); K = (K+K')/2;
    U(freedofs,:)=K(freedofs, freedofs)\F(freedofs,:);
%% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
    ts_SA = tic;
    c = 0; dc = 0;
    for i = 1:size(F,2)
        Ui = U(:,i);
        ce = sum((Ui(edofMatn)*KE).*Ui(edofMatn),2);                        % make a change here
        c= c+sum((vxPhys.*E0+(1-vxPhys).*(Emin*E0)).*ce);                   % make a change here
        dc = dc-((1-vxPhys)*Emin+vxPhys).*E0.*ce;
    end
    dv = ones(nnele,1); cc(loop) = c;                                       % make a change here
%% FILTERING AND MODIFICATION OF SENSITIVITIES
    dc = H*(dc./Hs);                                                        % make a change here  
    dv = H*(dv./Hs);                                                        % make a change here
    %% OPTIMALITY CRITERIA UPDATE
    l1 = 0; l2 = 1e9; move = 0.1;
    while (l2-l1)/(l1+l2) > 1e-3
        lmid = 0.5*(l2+l1);
        vxnew = max(0,max(vxPhys-move,min(1,min(vxPhys+move,vxPhys.*sqrt(-dc./dv/lmid)))));
        if sum(vxnew) > vol*nnele, l1 = lmid; else l2 = lmid; end
    end
    vxPhys = (H*vxnew)./Hs;
%     %% METHOD OF MOVING ASYMPTOTES
%     xval  = vx;                                                             % make a change here
%     f0val = c;
%     df0dx = dc;                                                             % make a change here
%     fval  = sum(vxPhys)/(vol*nnele) - 1;                                    % make a change here
%     dfdx  = dv' / (vol*nnele);                                              % make a change here
%     [vxmma, ~, ~, ~, ~, ~, ~, ~, ~, low,upp] = ...
%     mmasub(m, n, loop, xval, xmin, xmax, xold1, xold2, ...
%     f0val,df0dx,fval,dfdx,low,upp,a0,a,c_MMA,d);
% %% Update MMA Variables
%     vxnew     = vxmma;                                                      % make a change here
%     vxPhys = (H*vxnew)./Hs;                                                 % make a change here 
%     xold2    = xold1(:);
%     xold1    = vx(:);
    %% Change to total number of elements
    vxPhys1 = zeros(nely,nelx,nelz); vxPhys2 = vxPhys1(:);                  % make a change here 
    vxPhys2(ele) = vxPhys; vxPhys = vxPhys2; vxPhys(MusD) = 1;              % make a change here
    %% Apply smooth edge algorithm
    [vxPhys,xg,lss,top,tol] = smoothedge3D(vxPhys,Hn,Hns,nelx,nely,nelz,nele,nnele,nodex,nodey,nodez,fnx,fny,fnz,beta,ngrid);
    %% CHECK CONVERGENCE
        change = sum(abs(vxnew(:)-vx(:)))/(vol*nnele);                      % make a change here 
        vx=vxnew; 
    if (change <= tolx || tol<=tol_thresh) || loop >= maxloop
       vxPhys(MusD) = 1;
       [vxPhys,xg,lss,top,tol] = smoothedge3D(vxPhys,Hn,Hns,nelx,nely,nelz,nele,nnele,nodex,nodey,nodez,fnx,fny,fnz,beta,ngrid); 
    end 
    %% PRINT RESULTS
        fvol(loop) = sum(vxPhys(:))/(nnele); 
        fprintf('It.:%5i Obj.:%11.3f Vol.:%7.3f ch.:%7.5f Topo.:%7.5f\n'...
            ,loop,cc(loop),fvol(loop),change,tol);              % make a change here
        figure(1)
        clf;
        display_3Dsmooth(flip(xg,3),flip(top,3));
        hold on;
        figure(2)
        clf;
        yyaxis left
        plot(cc,'b-','LineWidth',1.8);                                                      % make a change here 
        set(gca,'FontName','Times New Roman','Box','On','LineWidth',2,'FontSize',15,'defaultAxesColorOrder',[0 0 1; 1 0 0]);
        %plot(loop,cc(loop),'b*','LineWidth',2,'MarkerSize',8)               % make a change here
        xlabel('Iteration','Interpreter','Latex')
        ylabel('Compliance (Nm)', 'rotation', 90,'Interpreter','Latex')
        ax = gca; ax.YColor = 'b'; 
        
        yyaxis right
        plot(fvol,'r-','LineWidth',1.8);
        ylabel('Volume fraction', 'rotation', 270,'Interpreter','Latex'); ax = gca; ax.YColor = 'r'; 
        pause(1e-6); hold off
        %% UPDATE% HEAVISIDE REGULARIZATION PARAMETER
        if beta < 2
         beta=beta+ER;
        end
        fprintf('Parameter beta increased to %g.\n',beta);
        vxPhys_full = reshape(vxPhys,nely,nelx,nelz);
        vxPhys = vxPhys(ele);
end
cc_end = cc(end); S = struct('comp',cc_end,'finalvol',fvol(end),'eleden',vxPhys_full,'gridden',xg,'elenum1',nnele,'elenum2',nele);