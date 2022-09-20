% Vectorized assembly of stiffness matrix
% S = \int grad(phi_i) * grad(phi_j) dx

function [Sx]=vec_localconv_mixed_1d(edet,dFinv,FE1,FE2)


mr       = FE1.mr;
wr       = FE1.wr;
ndim     = FE1.ndim;

nphi1     = FE1.nphi;
phi1      = FE1.phi;
gradphi1  = FE1.gradphi;

nphi2     = FE2.nphi;
phi2      = FE2.phi;
gradphi2  = FE2.gradphi;

nelement = size(edet,1);

Sx        = zeros(nelement,nphi1,nphi2);
Sy        = zeros(nelement,nphi1,nphi2);
dphi2     = zeros(nelement,nphi2,ndim,mr);

for q=1:mr
    for r1=1:ndim
    for r2=1:ndim
    for i=1:nphi2
        dphi2(:,i,r2,q) = dphi2(:,i,r2,q) + gradphi2(i,r1,q) * dFinv(:,r1,r2,q);
    end
    end
    end
    
    fac = edet(:,q) * wr(q);
    
    for i1=1:nphi1
    for i2=1:nphi2
        Sx(:,i1,i2) = Sx(:,i1,i2) + phi1(i1,q).*dphi2(:,i2,1,q).*fac;
        % Sy(:,i1,i2) = Sy(:,i1,i2) + phi1(i1,q).*dphi2(:,i2,2,q).*fac;
    end
    end
end
