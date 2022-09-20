% Vectorized assembly of mass matrix
% M = \int phi_i * phi_j dx

function M=vec_localmass_weighted(edet,FE,dofmap,func,weight)

nphi     = FE.nphi;
phi      = FE.phi;
mr       = FE.mr;
wr       = FE.wr;
nelement = size(edet,1);
M        = zeros(nelement,nphi,nphi);

for q=1:mr
    w = zeros(nelement,1);  
    for i=1:nphi
        w = w + weight(dofmap(:,i))*phi(i,q);
    end
    
    fac = func(w) .* edet(:,q) * wr(q);
    
    for i1=1:nphi
    for i2=1:nphi
        M(:,i1,i2) = M(:,i1,i2) + phi(i1,q)*phi(i2,q) * fac;
    end
    end
    
end
