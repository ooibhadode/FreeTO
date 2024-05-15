function stlgen(top,fnx,fny,fnz,dx,dy,dz,nm)
%% check model name nm
nm1 = split(nm,'.');
if length(nm1) > 1
    nm = strcat(string(nm1{1}),'.stl');
else 
    nm = strcat(string(nm1),'.stl');
end
% Dimensions of initial phi
st_max = max([dx dy dz]);
phi = top; phi = flip(phi,1); 
gmax = max([max(fnx(:)) max(fny(:)) max(fnz(:))]);
gx = (fnx/gmax)*st_max; gy = (fny/gmax)*st_max; gz = (fnz/gmax)*st_max;

phi_stl = phi; % if positive/negative defines material (+/-phi).
limit = 0.00001; % the sign of the limit should match the sign of phi.
phi_stl(phi_stl<=limit) = -1;
phi_stl(phi_stl>limit) = 1;
phi_stl = smooth3(phi);

[faces1, vertices1] = isosurface(gx, gy, gz, phi_stl, 0);
[faces2, vertices2] = isocaps(gx, gy, gz, phi_stl, 0);
faces = [faces1; length(vertices1(:,1)) + faces2];
vertices = [vertices1; vertices2];
triangles = triangulation(faces, vertices);

normals = faceNormal(triangles);
vertices=vertices'; faces=faces'; normals=normals';
% Transform vertices into single precision.
datav = single(vertices);
% Put all vertices of the faces in a single 2D matrix (3 x m).
datavxyz = datav(:,faces);
% Transform data to [V1, V2, V3] for each face in a (3 x 3 x m/3) matrix.
dataxyz = reshape(datavxyz,3,3,numel(datavxyz)/9);
% Reshape normals from (3 x m) to (3 x 1 x m) matrix.
normals = reshape(normals,3,1,numel(normals)/3);
% Include normal vectors: [n, V1, V2, V3] in a (3 x 4 x m/3) matrix.
datanxyz = [normals(:,:,:), dataxyz(:,:,:)];
% Transform datanxyz to 16-bit unsigned integer.
datanxyz = typecast(datanxyz(:), 'uint16');
% Reshape data. Each column will contain the information of one face.
data = reshape(datanxyz(:), 24, numel(datanxyz)/24);
% The 25th is the attribute byte count.
data(25,:) = 0;

% Generating stl file in Binary mode
fileID = fopen(nm, 'w'); % name of stl file can be changed.
% Write a title up to 80 characters. (80 bytes)
fprintf(fileID, '%-80s',...
   'exported using phi2stl created by P. Vogiatzis (Advisor: S. Chen)');
% Write the number of faces in 32-bit unsigned integer.
fwrite(fileID, numel(data)/25, 'uint32');
% Write data.
fwrite(fileID, data, 'uint16');
fclose(fileID);
end
