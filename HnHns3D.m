function [Hn,Hns]=HnHns3D(nelx,nely,nelz,rnmin)
inH = ones((nelx+1)*(nely+1)*(nelz+1)*(2*(ceil(rnmin)+1))^2,1);
jnH = ones(size(inH)); snH = zeros(size(inH)); kn = 0;
[elex,eley,elez] = meshgrid(1.5:nelx+0.5,1.5:nely+0.5,1.5:nelz+0.5);
for kn1 = 1:nelz+1
    for in1 = 1:nelx+1
        for jn1 = 1:nely+1
            en1 = (kn1-1)*(nelx+1)*(nely+1)+(in1-1)*(nely+1)+jn1;
            for kn2 = max(kn1-(ceil(rnmin)),1):min(kn1+(ceil(rnmin)-1),nelz)
                for in2 = max(in1-(ceil(rnmin)),1):min(in1+(ceil(rnmin)-1),nelx)
                    for jn2 = max(jn1-(ceil(rnmin)),1):min(jn1+(ceil(rnmin)-1),nely)
                        en2 = (kn2-1)*nelx*nely + (in2-1)*nely+jn2;
                        kn = kn+1; inH(kn) = en1; jnH(kn) = en2;                        
                        snH(kn) = max(0,rnmin-sqrt((in1-elex(jn2,in2,kn2))^2+(jn1-eley(jn2,in2,kn2))^2+(kn1-elez(jn2,in2,kn2))^2));
                    end
                end
            end
        end
    end
end
Hn = sparse(inH,jnH,snH); Hns = sum(Hn,2);