function [vxPhys,xg,ls,top,tol] = smoothedge3D(vxPhys,Hn,Hns,nelx,nely,nelz,nele,nnele,nodex,nodey,nodez,fnx,fny,fnz,beta,ngrid)
xn = reshape((Hn*vxPhys(:)./Hns),nely+1,nelx+1,nelz+1);
%% UPDATE POINT DESNIGY BY A HEAVISIDE SMOOTH FUNCTION/ HEAVISIDE STEP FUNCTION
 xg = interp3(nodex,nodey,nodez,xn,fnx,fny,fnz,'linear');
 l1 =0; l2 = 1;
     while (l2-l1) > 1.0e-5
      ls = (l1+l2)/2.0;
    %% Heaviside smooth function
      xgnew = max(0.001,(tanh(beta*ls)+tanh(beta*(xg-ls)))/(tanh(beta*ls)+tanh(beta*(1-ls)))); 
          if sum(sum(sum(xgnew)))/((ngrid*nelx+1)*(ngrid*nely+1)*(ngrid*nelz+1)) - sum(vxPhys(:))/(nele) > 0
             l1 = ls;
          else
             l2 = ls;
          end
     end
    %% CONVERTING TO ELEMENTS
     vxPhys(:) = 0;
     Terr = 0;
     Tm=[];
     for nk = 1:nelz
        for ni = 1:nelx
            for nj = 1:nely
                ne = (nk-1)*nelx*nely+(ni-1)*nely+nj;
                for nk1 = ngrid*(nk-1)+1:ngrid*nk+1
                    for ni1 = ngrid*(ni-1)+1:ngrid*ni+1
                        for nj1 = ngrid*(nj-1)+1:ngrid*nj+1
                            Tm=[Tm;xgnew(nj1,ni1,nk1)]; 
                            vxPhys(ne) = vxPhys(ne)+xgnew(nj1,ni1,nk1);                      
                         end
                    end
                 end
                    if min(Tm)>0.001 && max(Tm)<1
                        Terr=Terr+1;
                    end
                    Tm=[];
            end
        end
     end
    vxPhys = vxPhys/((ngrid+1)^3);
    
    %% Topology Error
    tol = Terr/(nnele);  
    top = xg-ls;