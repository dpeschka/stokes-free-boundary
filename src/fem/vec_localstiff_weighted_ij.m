% Vectorized assembly of stiffness matrix
% S = \int grad(phi_i) * grad(phi_j) dx

function S=vec_localstiff_weighted_ij(edet,dFinv,FE,dofmap,func,weight,k,l)

nphi     = FE.nphi;
phi      = FE.phi;
mr       = FE.mr;
wr       = FE.wr;
ndim     = FE.ndim;
gradphi  = FE.gradphi;

nelement = size(edet,1);

S        = zeros(nelement,nphi,nphi);
dphi     = zeros(nelement,nphi,ndim,mr);

for q=1:mr
    for i=1:nphi
    for r1=1:ndim
    for r2=1:ndim
        dphi(:,i,r2,q) = dphi(:,i,r2,q) + gradphi(i,r1,q) * dFinv(:,r1,r2,q);
    end
    end
    end
    
    w = zeros(nelement,1);  
    for i=1:nphi
        w = w + weight(dofmap(:,i))*phi(i,q);
    end
    
    fac = func(w) .* edet(:,q) * wr(q);
    
    for i1=1:nphi
    for i2=1:nphi
        S(:,i1,i2) = S(:,i1,i2) + dphi(:,i1,k,q).*dphi(:,i2,l,q).*fac;
    end
    end
end
