function [F,nf] = forcevec(Fn,ndof,Fmagx,Fmagy,Fmagz,loadtype)

nf = max([length(Fmagx) length(Fmagy) length(Fmagz)]); 
Fx = zeros(ndof,nf); Fy = zeros(ndof,nf); Fz = zeros(ndof,nf);
for i = 1:nf
    frn = full(Fn(:,i)); frn(frn == 0) = [];
    if strcmp(loadtype,'point')
        frn = sort(frn); frn = frn(round(length(frn)/2));
    end
    frnm = [frn*3-2 frn*3-1 frn*3];
    if length(Fmagx) > 1
        Fx(frnm(:,1),i) = Fmagx(i)/length(frnm(:,1)); 
    else 
        Fx(frnm(:,1),1) = Fmagx/length(frnm(:,1));
    end 
    if length(Fmagy) > 1
        Fy(frnm(:,2),i) = Fmagy(i)/length(frnm(:,2)); 
    else 
        Fy(frnm(:,2),1) = Fmagy/length(frnm(:,2));
    end
    if length(Fmagz) > 1
        Fz(frnm(:,3),i) = Fmagz(i)/length(frnm(:,3));
    else 
        Fz(frnm(:,3),1) = Fmagz/length(frnm(:,3));
    end 
end 
F = Fx+Fy+Fz; F(:,all(F == 0)) = []; F = sparse(F);