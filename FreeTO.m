%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    FreeTO - 3D Topology Optimization of Freeform Structures      %%%
%%%%            Author: Osezua Ibhadode, March 2024               %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function S = FreeTO(domain,force1,MeshControl,volfrac,varargin)
%% Examples
%FreeTO('STLs\GE_domain.stl','STLs\GE_force.stl',70,0.3,'fixed','STLs\GE_fixed.stl','force2','STLs\GE_force.stl','Fmagy',[0 -2000],'Fmagz',[1500 0],'YoungsModulus',210e9)
%FreeTO('STLs\air_domain.stl','STLs\air_force.stl',90,0.2,'fixed','STLs\air_fixed.stl','zfixed','STLs\air_zfixed.stl','force2','STLs\air_force.stl','force3','STLs\air_force.stl','Fmagx',[1000 1324 0],'Fmagy',[0 -1324 -2500],'YoungsModulus',210e9,'Symmetry1',"x-y",'direction1',"right",'optimization','SEMDOT','modelName','STLs\air_brack_TO.stl')
%FreeTO('STLs\hand_domain.stl','STLs\hand_force1.stl',90,0.3,'fixed','STLs\hand_fixed.stl','force2','STLs\hand_force2.stl','force3','STLs\hand_force3.stl','force4','STLs\hand_force4.stl','force5','STLs\hand_force5.stl','Fmagz',2e3*ones(1,5),'YoungsModulus',210e9)
%FreeTO('STLs\quad_domain.stl','STLs\quad_force1.stl',70,0.3,'fixed','STLs\quad_fixed.stl','xfixed','STLs\quad_xfixed.stl','yfixed','STLs\quad_yfixed.stl','force2','STLs\quad_force2.stl','force3','STLs\quad_force3.stl','Fmagz',[-1500 -1500 -1000],'YoungsModulus',2e9,'Symmetry1','y-z','Symmetry2','z-x','direction1','right','direction2','left','modelName','STLs\quad_TO.stl','keep_BCz','yes')

%% Default values of parameters
default_fix = []; default_xfix = []; default_yfix = []; default_zfix = []; 
default_for2 = []; default_for3 = []; default_for4 = []; default_for5 = [];
default_for6 = []; default_for7 = []; default_for8 = []; default_for9 = [];
default_for10 = []; default_keepdom = [];
default_keep_BC = 'yes'; default_keep_BCx = 'no'; default_keep_BCy = 'no'; default_keep_BCz = 'no'; default_E0 = 1;
default_penal = 3; default_rmin = 1.5; default_Fmagx = 0; default_Fmagy = 0; default_Fmagz = 0;
default_text = []; default_load = 'distributed'; 
default_sym1 = []; default_dir1 = 'right'; default_sym2 = []; default_dir2 = 'right'; default_sym3 = []; default_dir3 = 'right'; 
default_v = 0.3; default_opt_type = 'SIMP';

%% Input Parser for required, optional and name-value parameters
pp = inputParser;
addRequired(pp,'domain'); addRequired(pp,'force1'); addRequired(pp,'MeshControl'); addRequired(pp,'volfrac');
addOptional(pp,'fixed',default_fix); addOptional(pp,'xfixed',default_xfix);
addOptional(pp,'yfixed',default_yfix); addOptional(pp,'zfixed',default_zfix);
addOptional(pp,'force2',default_for2); addOptional(pp,'force3',default_for3);
addOptional(pp,'force4',default_for4); addOptional(pp,'force5',default_for5);
addOptional(pp,'force6',default_for6); addOptional(pp,'force7',default_for7);
addOptional(pp,'force8',default_for8); addOptional(pp,'force9',default_for9); 
addOptional(pp,'force10',default_for10); addOptional(pp,'keepdom',default_keepdom);
addParameter(pp,'keep_BC',default_keep_BC); addParameter(pp,'keep_BCx',default_keep_BCx);
addParameter(pp,'keep_BCy',default_keep_BCy); addParameter(pp,'keep_BCz',default_keep_BCz);
addParameter(pp,'YoungsModulus',default_E0); addParameter(pp,'optimization',default_opt_type);
addParameter(pp,'penaltySIMP',default_penal); addParameter(pp,'filterRadius',default_rmin);
addParameter(pp,'Fmagx',default_Fmagx); addParameter(pp,'Fmagy',default_Fmagy);
addParameter(pp,'Fmagz',default_Fmagz); addParameter(pp,'modelName',default_text);
addParameter(pp,'loadtype',default_load); addParameter(pp,'PoissonRatio',default_v);
addParameter(pp,'Symmetry1',default_sym1); addParameter(pp,'direction1',default_dir1); 
addParameter(pp,'Symmetry2',default_sym2); addParameter(pp,'direction2',default_dir2);
addParameter(pp,'Symmetry3',default_sym3); addParameter(pp,'direction3',default_dir3);

parse(pp,domain,force1,MeshControl,volfrac,varargin{:});

%% Domain Initialization
domain = pp.Results.domain; MeshControl = pp.Results.MeshControl;

%% Optimization parameters
volfrac = pp.Results.volfrac; penal = pp.Results.penaltySIMP; 
rmin = pp.Results.filterRadius; E0 = pp.Results.YoungsModulus;v = pp.Results.PoissonRatio;
opt_type = pp.Results.optimization;

%% Load and boundary conditions
fixed = pp.Results.fixed; xfixed = pp.Results.xfixed;
yfixed = pp.Results.yfixed; zfixed = pp.Results.zfixed;
F1 = pp.Results.force1; F2 = pp.Results.force2; F3 = pp.Results.force3; F4 = pp.Results.force4;
F5 = pp.Results.force5; F6 = pp.Results.force6; F7 = pp.Results.force7; F8 = pp.Results.force8;
F9 = pp.Results.force9; F10 = pp.Results.force10; keep_domain = pp.Results.keepdom;
Fmagx = pp.Results.Fmagx; Fmagy = pp.Results.Fmagy; Fmagz = pp.Results.Fmagz; loadtype = pp.Results.loadtype;
keep_BC = pp.Results.keep_BC; keep_BCxyz = cellstr([string(pp.Results.keep_BCx) string(pp.Results.keep_BCy) string(pp.Results.keep_BCz)]);

%% Post-processing parameters
nm = pp.Results.modelName; 
symm = cellstr([string(pp.Results.Symmetry1) string(pp.Results.Symmetry2) string(pp.Results.Symmetry3)]); 
dir = cellstr([string(pp.Results.direction1) string(pp.Results.direction2) string(pp.Results.direction3)]);

%% Check input
if nargin == 4
    error('Please check that you included all required inputs in the right order and at least one support condition')
elseif isempty(fixed) && isempty(xfixed) && isempty(yfixed) && isempty(zfixed)
    error('Please include at least one support condition')
elseif all([Fmagx Fmagy Fmagz])
    error('Please include at least one non-zero load definition')
end 

%% Run topology optimization
if strcmpi(opt_type,'SIMP')
    % SIMP density-based
    [nelx,nely,nelz,ngrid,~,xg,top,lss,fnx,fny,fnz,dx,dy,dz,S] = SIMP(MeshControl,domain,fixed,xfixed,yfixed,...
        zfixed,F1,F2,F3,F4,F5,F6,F7,F8,F9,F10,keep_domain,keep_BC,keep_BCxyz,...
        volfrac,E0,v,penal,rmin,Fmagx,Fmagy,Fmagz,loadtype);
elseif strcmpi(opt_type,'SEMDOT')
    % SEMDOT density-based
    [nelx,nely,nelz,ngrid,~,xg,top,lss,fnx,fny,fnz,dx,dy,dz,S] = SEMDOT(MeshControl,domain,fixed,xfixed,yfixed,...
        zfixed,F1,F2,F3,F4,F5,F6,F7,F8,F9,F10,keep_domain,keep_BC,keep_BCxyz,...
        volfrac,E0,v,rmin,Fmagx,Fmagy,Fmagz,loadtype);
end

%% Post-processing
if ~isempty(symm)
    [xg,top,fnx,fny,fnz,nelx,nely,nelz,dx,dy,dz] = symmetry(xg,lss,nelx,nely,nelz,dx,dy,dz,ngrid,symm{1},dir{1});
    if length(symm) >= 2
        [xg,top,fnx,fny,fnz,nelx,nely,nelz,dx,dy,dz] = symmetry(xg,lss,nelx,nely,nelz,dx,dy,dz,ngrid,symm{2},dir{2});
        if length(symm) == 3
            [xg,top,fnx,fny,fnz,~,~,~,dx,dy,dz] = symmetry(xg,lss,nelx,nely,nelz,dx,dy,dz,ngrid,symm{3},dir{3});
        end
    end 
    figure (1)
    display_3Dsmooth(xg,top);
end
if ~isempty(nm)
    stlgen(top,fnx,fny,fnz,dx,dy,dz,nm)
end
end
% This Matlab code was written by Osezua Ibhadode                          %                                     %
% Please sent your comments to: ibhadode@ualberta.ca                       %
%                                                                          %
% The code is intended for educational purposes and theoretical details    %
% are discussed in a paper which is currently under review                 %  
%                                                                          %
% Disclaimer:                                                              %
% The author reserves all rights but do not guarantee that the code is     %
% free from errors. Furthermore, the author shall not be liable in any     %
% event caused by the use of the program.                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
