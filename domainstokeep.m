function outeM = domainstokeep(oute,nelx,nely,nelz,keep_BC,keep_BCxyz,stp2,stp3,stp4,stp5,...
    stp6_1,stp7,stf2,stf3,stf4,stf5,stf6_1,stf7,nnf1,nnf2,x1_cen,y1_cen,z1_cen)
%% Collect the vertices and faces of the boundary condition regions to keep during optimization
if strcmp(keep_BC,'yes')
    % keep_BCxyz{1} = xfixed, keep_BCxyz{2} = yfixed, keep_BCxyz{3} = zfixed
    if  strcmp(keep_BCxyz{1},'no') && strcmp(keep_BCxyz{2},'no') && strcmp(keep_BCxyz{3},'no')
        stpM = [stp2;stp6_1;stp7]; stfM = [stf2;stf6_1;stf7];
        nn1 = [size(stp2,1) nnf1 size(stp7,1)];
        nn2 = [size(stf2,1) nnf2 size(stf7,1)];
        nn1(nn1 == 0) = []; nn2(nn2 == 0) = [];
    elseif strcmp(keep_BCxyz{1},'yes') && strcmp(keep_BCxyz{2},'no') && strcmp(keep_BCxyz{3},'no')
            stpM = [stp2;stp3;stp6_1;stp7]; stfM = [stf2;stf3;stf6_1;stf7];
            nn1 = [size(stp2,1) size(stp3,1) nnf1 size(stp7,1)];
            nn2 = [size(stf2,1) size(stf3,1) nnf2 size(stf7,1)];
            nn1(nn1 == 0) = []; nn2(nn2 == 0) = [];
    elseif strcmp(keep_BCxyz{1},'no') && strcmp(keep_BCxyz{2},'yes') && strcmp(keep_BCxyz{3},'no')
        stpM = [stp2;stp4;stp6_1;stp7]; stfM = [stf2;stf4;stf6_1;stf7];
        nn1 = [size(stp2,1) size(stp4,1) nnf1 size(stp7,1)];
        nn2 = [size(stf2,1) size(stf4,1) nnf2 size(stf7,1)];
        nn1(nn1 == 0) = []; nn2(nn2 == 0) = [];
    elseif strcmp(keep_BCxyz{1},'no') && strcmp(keep_BCxyz{2},'no') && strcmp(keep_BCxyz{3},'yes')
        stpM = [stp2;stp5;stp6_1;stp7]; stfM = [stf2;stf5;stf6_1;stf7];
        nn1 = [size(stp2,1) size(stp5,1) nnf1 size(stp7,1)];
        nn2 = [size(stf2,1) size(stf5,1) nnf2 size(stf7,1)];
        nn1(nn1 == 0) = []; nn2(nn2 == 0) = [];
    elseif strcmp(keep_BCxyz{1},'yes') && strcmp(keep_BCxyz{2},'yes') && strcmp(keep_BCxyz{3},'no')
        stpM = [stp2;stp3;stp4;stp6_1;stp7]; stfM = [stf2;stf3;stf4;stf6_1;stf7];
        nn1 = [size(stp2,1) size(stp3,1) size(stp4,1) nnf1 size(stp7,1)];
        nn2 = [size(stf2,1) size(stf3,1) size(stf4,1) nnf2 size(stf7,1)];
        nn1(nn1 == 0) = []; nn2(nn2 == 0) = [];
    elseif strcmp(keep_BCxyz{1},'yes') && strcmp(keep_BCxyz{2},'no') && strcmp(keep_BCxyz{3},'yes')
        stpM = [stp2;stp3;stp5;stp6_1;stp7]; stfM = [stf2;stf3;stf5;stf6_1;stf7];
        nn1 = [size(stp2,1) size(stp3,1) size(stp5,1) nnf1 size(stp7,1)];
        nn2 = [size(stf2,1) size(stf3,1) size(stf5,1) nnf2 size(stf7,1)];
        nn1(nn1 == 0) = []; nn2(nn2 == 0) = [];
    elseif strcmp(keep_BCxyz{1},'no') && strcmp(keep_BCxyz{2},'yes') && strcmp(keep_BCxyz{3},'yes')
        stpM = [stp2;stp4;stp5;stp6_1;stp7]; stfM = [stf2;stf4;stf5;stf6_1;stf7];
        nn1 = [size(stp2,1) size(stp4,1) size(stp5,1) nnf1 size(stp7,1)];
        nn2 = [size(stf2,1) size(stf4,1) size(stf5,1) nnf2 size(stf7,1)];
        nn1(nn1 == 0) = []; nn2(nn2 == 0) = [];
    elseif strcmp(keep_BCxyz{1},'yes') && strcmp(keep_BCxyz{2},'yes') && strcmp(keep_BCxyz{3},'yes')
        stpM = [stp2;stp3;stp4;stp5;stp6_1;stp7]; stfM = [stf2;stf3;stf4;stf5;stf6_1;stf7];
        nn1 = [size(stp2,1) size(stp3,1) size(stp4,1) size(stp5,1) nnf1 size(stp7,1)];
        nn2 = [size(stf2,1) size(stf3,1) size(stf4,1) size(stf5,1) nnf2 size(stf7,1)];
        nn1(nn1 == 0) = []; nn2(nn2 == 0) = [];
    else
        stpM = [stp2;stp6_1;stp7]; stfM = [stf2;stf6_1;stf7]; 
        nn1 = [size(stp2,1) nnf1 size(stp7,1)];
        nn2 = [size(stf2,1) nnf2 size(stf7,1)];
        nn1(nn1 == 0) = []; nn2(nn2 == 0) = [];
    end
else
    stpM = stp7; stfM = stf7; nn1 = size(stp7,1); nn2 = size(stf7,1);
end
% Identify elements according to the selected vertices and faces
outeM = zeros(size(oute(:)));
for i = 1:length(nn1)
    j1 = sum(nn1(1:i))-nn1(i)+1; k1 = j1+nn1(i)-1;
    j2 = sum(nn2(1:i))-nn2(i)+1; k2 = j2+nn2(i)-1;
    outeM = outeM+double(intriangulation(stpM(j1:k1,:),stfM(j2:k2,:),[x1_cen(:),y1_cen(:),z1_cen(:)])); 
    outeM(outeM<0) = 0; outeM(outeM>1) = 1;
end 
outeM = flip(reshape(outeM,nely,nelx,nelz),1);