% Vectorized assembly of convection matrix
% S = \int grad(phi_i) * grad(phi_j) dx

function C=vec_localconv_cubic(edet,dFinv,FE)

nphi     = FE.nphi;
phi      = FE.phi;
mr       = FE.mr;
wr       = FE.wr;
ndim     = FE.ndim;
gradphi  = FE.gradphi;

nelement = size(edet,1);

C        = zeros(nelement,nphi,nphi,nphi,ndim);
dphi     = zeros(nelement,nphi,ndim,mr);

for q=1:mr
  
  % prepare dphi
  for i=1:nphi
    for r1=1:ndim
      for r2=1:ndim
        dphi(:,i,r2,q) = dphi(:,i,r2,q) + gradphi(i,r1,q) * dFinv(:,r1,r2,q);
      end
    end
  end
  
  % compute matrix elements
  fac = edet(:,q) * wr(q);
  for i1=1:nphi
    for i2=1:nphi
      for i3=1:nphi
        C(:,i1,i2,i3,1) = C(:,i1,i2,i3,1) + phi(i1,q)*phi(i2,q)*dphi(:,i3,1,q).*fac;
        C(:,i1,i2,i3,2) = C(:,i1,i2,i3,2) + phi(i1,q)*phi(i2,q)*dphi(:,i3,2,q).*fac;
      end
    end
  end
  
end