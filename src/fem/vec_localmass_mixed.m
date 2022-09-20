% Vectorized assembly of mass matrix
% M = \int phi_i * phi_j dx

function M=vec_localmass_mixed(edet,FE1,FE2)


mr       = FE1.mr;
wr       = FE1.wr;
nelement = size(edet,1);

nphi1    = FE1.nphi;
phi1     = FE1.phi;
nphi2    = FE2.nphi;
phi2     = FE2.phi;

M        = zeros(nelement,nphi1,nphi2);

for q=1:mr
    
    fac = edet(:,q) * wr(q);
    
    for i1=1:nphi1
    for i2=1:nphi2
        M(:,i1,i2) = M(:,i1,i2) + phi1(i1,q)*phi2(i2,q) * fac;
    end
    end
    
end
