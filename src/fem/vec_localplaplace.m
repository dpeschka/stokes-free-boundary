% Vectorized assembly of stiffness matrix
% S = \int grad(phi_i) * grad(phi_j) dx
function [S,F]=vec_localplaplace(edet,dFinv,FE,e2,u,p,kappa)


nphi     = FE.nphi;
phi      = FE.phi;
mr       = FE.mr;
wr       = FE.wr;
ndim     = FE.ndim;
gradphi  = FE.gradphi;

nelement = size(edet,1);

S        = zeros(nelement,nphi,nphi);
F        = zeros(nelement,nphi,nphi);
dphi     = zeros(nelement,nphi,ndim,mr);
gradu    = zeros(nelement,ndim,mr);
    
for q=1:mr
    for i=1:nphi
    for r1=1:ndim
    for r2=1:ndim
        dphi(:,i,r2,q) = dphi(:,i,r2,q) + gradphi(i,r1,q) * dFinv(:,r1,r2,q);
    end
    end
    end
    
    for i=1:nphi
        gradu(:,1,q) = gradu(:,1,q) + dphi(:,i,1,q).*u(e2(:,i));
        gradu(:,2,q) = gradu(:,2,q) + dphi(:,i,2,q).*u(e2(:,i));
    end
    
    gradu_n = (kappa + gradu(:,1,q).^2 + gradu(:,2,q).^2 ).^(1/2);
    
    fac = gradu_n.^(p-2) .* edet(:,q) * wr(q);
    
    for i1=1:nphi
    for i2=1:nphi
    for r=1:ndim
        S(:,i1,i2) = S(:,i1,i2) + dphi(:,i1,r,q).*dphi(:,i2,r,q).*fac;
    end
    end
    end
    
    fac = (p-2) * gradu_n.^(p-4) .* edet(:,q) * wr(q);
    
    for i1=1:nphi
    for i2=1:nphi
    for k=1:ndim
    for l=1:ndim
        F(:,i1,i2) = F(:,i1,i2) + dphi(:,i1,k,q).*gradu(:,k,q).*gradu(:,l,q).*dphi(:,i2,l,q).*fac;
    end
    end
    end
    end
    
end
