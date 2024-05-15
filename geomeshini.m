function [out_p,oute,outeM,sup_all,sup_x,sup_y,sup_z,Fn,nelx,nely,nelz,nele,ndof,stp1,a,...
    del_x,del_y,del_z] = geomeshini(MeshControl,inp1,inp2,inp3,...
    inp4,inp5,inp6,inp7,inp8,inp9,inp10,inp11,inp12,inp13,inp14,inp15,inp16,keep_BC,keep_BCxyz)
inp_f = {inp6,inp7,inp8,inp9,inp10,inp11,inp12,inp13,inp14,inp15}; inp_f(~cellfun(@ischar,inp_f)) = [];
% get stl vertices for inp1
st1 = stlread(inp1);                                                        % read stl
stp1 = st1.Points; stf1 = st1.ConnectivityList;                             % obtain vertices and faces
% get bounding box
stp_xmax = max(stp1(:,1)); stp_xmin = min(stp1(:,1));                       % get the minimum and maximum x for domain
del_x = stp_xmax-stp_xmin;
stp_ymax = max(stp1(:,2)); stp_ymin = min(stp1(:,2));                       % get the minimum and maximum y for domain
del_y = stp_ymax-stp_ymin;
stp_zmax = max(stp1(:,3)); stp_zmin = min(stp1(:,3));                       % get the minimum and maximum z for domain
del_z = stp_zmax-stp_zmin;
stt = abs([stp_xmin stp_xmax stp_ymin stp_ymax stp_zmin stp_zmax]);         % get the maximum of all maximums for domain
stm = max(stt); stff = find(stt == stm(1));                                 % get the position of the max of all maxs

if stff(1) == 1 || stff(1) == 2                                                             % if the maximum is in position 1
    gpx = MeshControl; 
    x_old = linspace(stp_xmin, stp_xmax, gpx); ssz_old = (x_old(2)-x_old(1))*1e-2;
    x = linspace(stp_xmin+ssz_old, stp_xmax-ssz_old, gpx);
    ssz = x(2)-x(1); gpx = length(x); a = ssz/2000;           % the MeshControl should be defined from the x-axis
    y_old = stp_ymin:ssz:stp_ymax; diff_y = stp_ymax-y_old(end);                    % obtain y MeshControl
    y = stp_ymin+diff_y/2:ssz:stp_ymax; gpy = length(y);
    z_old = stp_zmin:ssz:stp_zmax; diff_z = stp_zmax-z_old(end);                    % obtain z MeshControl
    z = stp_zmin+diff_z/2:ssz:stp_zmax; gpz = length(z);
    x_cen = x+ssz/2; x_cen(end) = []; y_cen = y+ssz/2; y_cen(end) = []; z_cen = z+ssz/2; z_cen(end) = [];
elseif stff(1) == 3 || stff(1) == 4
    gpy = MeshControl; 
    y_old = linspace(stp_ymin, stp_ymax, gpy); ssz_old = (y_old(2)-y_old(1))*1e-2;
    y = linspace(stp_ymin+ssz_old, stp_ymax-ssz_old, gpy);
    ssz = y(2)-y(1); gpy = length(y); a = ssz/2000; 
    x_old = stp_xmin:ssz:stp_xmax; diff_x = stp_xmax-x_old(end); 
    x = stp_xmin+diff_x/2:ssz:stp_xmax; gpx = length(x);
    z_old = stp_zmin:ssz:stp_zmax; diff_z = stp_zmax-z_old(end); 
    z = stp_zmin+diff_z/2:ssz:stp_zmax; gpz = length(z);
    x_cen = x+ssz/2; x_cen(end) = []; y_cen = y+ssz/2; y_cen(end) = []; z_cen = z+ssz/2; z_cen(end) = [];
else
    gpz = MeshControl; 
    z_old = linspace(stp_zmin, stp_zmax, gpz); ssz_old = (z_old(2)-z_old(1))*1e-2;
    z = linspace(stp_zmin+ssz_old,stp_zmax-ssz_old,gpz);
    ssz = z(2)-z(1); gpz = length(z); a = ssz/2000; 
    x_old = stp_xmin:ssz:stp_xmax; diff_x = stp_xmax-x_old(end); 
    x = stp_xmin+diff_x/2:ssz:stp_xmax; gpx = length(x);
    y_old = stp_ymin:ssz:stp_ymax; diff_y = stp_ymax-y_old(end);                    % obtain y MeshControl
    y = stp_ymin+diff_y/2:ssz:stp_ymax; gpy = length(y); 
    x_cen = x+ssz/2; x_cen(end) = []; y_cen = y+ssz/2; y_cen(end) = []; z_cen = z+ssz/2; z_cen(end) = [];
end

%create meshgrid for x, y, and z
[x1,y1,z1] = meshgrid(x,y,z);                                               % create a meshgrid from x,y,z MeshControl
[x1_cen,y1_cen,z1_cen] = meshgrid(x_cen,y_cen,z_cen);
% for inp1 - Domain
out = double(intriangulation(stp1,stf1,[x1(:),y1(:),z1(:)])); out(out<0) = 0; % check which points lie in and out of the domain
% for inp2 - Fixed all DOFs support
if ~isempty(inp2)
    st2 = stlread(inp2);
    stp2 = st2.Points; stf2 = st2.ConnectivityList;
    out2 = intriangulation(stp2,stf2,[x1(:),y1(:),z1(:)]); 
    out2 = flip(reshape(out2,gpy,gpx,gpz),1);
    sup_all = find(out2(:)>0);                                              % identify support nodes fixd in all directions in the domain
else
    sup_all = []; stp2 = []; stf2 = [];
end 
% for inp3 - Fixed in x
if ~isempty(inp3)
    st3 = stlread(inp3);
    stp3 = st3.Points; stf3 = st3.ConnectivityList;
    out3 = intriangulation(stp3,stf3,[x1(:),y1(:),z1(:)]); 
    out3 = flip(reshape(out3,gpy,gpx,gpz),1);
    sup_x = find(out3(:)>0);
else 
    sup_x = []; stp3 = []; stf3 = [];
end
% for inp4 - Fixed in y
if ~isempty(inp4)
    st4 = stlread(inp4);
    stp4 = st4.Points; stf4 = st4.ConnectivityList;
    out4 = intriangulation(stp4,stf4,[x1(:),y1(:),z1(:)]); 
    out4 = flip(reshape(out4,gpy,gpx,gpz),1);
    sup_y = find(out4(:)>0);
else
    sup_y = []; stp4 = []; stf4 = [];
end
% for inp5 - Fixed in z
if ~isempty(inp5)
    st5 = stlread(inp5);
    stp5 = st5.Points; stf5 = st5.ConnectivityList;
    out5 = intriangulation(stp5,stf5,[x1(:),y1(:),z1(:)]); 
    out5 = flip(reshape(out5,gpy,gpx,gpz),1);
    sup_z = find(out5(:)>0);
else
    sup_z = []; stp5 = []; stf5 = [];
end
% for inp6 - Force
Fn = sparse(length(out),length(inp_f)); stp6_1 = []; stf6_1 = [];
for i = 1:length(inp_f)
    st6 = stlread(inp_f{i});
    stp6 = st6.Points; stp6_1 = [stp6_1;stp6]; nnf1(i) = size(stp6,1);
    stf6 = st6.ConnectivityList; stf6_1 = [stf6_1;stf6]; nnf2(i) = size(stf6,1);
    out6 = intriangulation(stp6,stf6,[x1(:),y1(:),z1(:)]); 
    out6 = flip(reshape(out6,gpy,gpx,gpz),1);
    Fn(1:length(find(out6(:)>0)),i) = find(out6(:)>0);
end 
% for domains to keep
if ~isempty(inp16)
    st7 = stlread(inp16);
    stp7 = st7.Points; stf7 = st7.ConnectivityList;
    out7 = intriangulation(stp7,stf7,[x1(:),y1(:),z1(:)]); 
    out7 = flip(reshape(out7,gpy,gpx,gpz),1);
    dom_ke = find(out7(:)>0);
else
    dom_ke = []; stp7 = []; stf7 = [];
end
out_p = flip(reshape(out,gpy,gpx,gpz),1);                                   % flip the identified points
nelx = size(out_p,2)-1; nely = size(out_p,1)-1; nelz = size(out_p,3)-1; nele = nelx*nely*nelz; % obtain nelx, nely,nelz and nele
ndof = 3*(nelx+1)*(nely+1)*(nelz+1);                                        % obtain number of degrees of freedom
%% get elements in the domain
oute = double(intriangulation(stp1,stf1,[x1_cen(:),y1_cen(:),z1_cen(:)])); oute(oute<0) = 0; 
oute = flip(reshape(oute,nely,nelx,nelz),1);
%% get elements to keep during optimization
outeM = domainstokeep(oute,nelx,nely,nelz,keep_BC,keep_BCxyz,stp2,stp3,stp4,stp5,...
    stp6_1,stp7,stf2,stf3,stf4,stf5,stf6_1,stf7,nnf1,nnf2,x1_cen,y1_cen,z1_cen);