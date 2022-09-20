% Vectorized assembly of convection matrix
% S = \int grad(phi_i) * grad(phi_j) dx

function C=vec_localconv(edet,dFinv,FE,ux,uy,dofmap)

nphi     = FE.nphi;
phi      = FE.phi;
mr       = FE.mr;
wr       = FE.wr;
ndim     = FE.ndim;
gradphi  = FE.gradphi;

nelement = size(edet,1);

C        = zeros(nelement,nphi,nphi);
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
  
  % prepare wx,wy
  wx = zeros(nelement,1);
  wy = zeros(nelement,1);
  for i=1:nphi
    wx = wx + ux(dofmap(:,i))*phi(i,q);
    wy = wy + uy(dofmap(:,i))*phi(i,q);
  end
  
  % compute matrix elements
  fac = edet(:,q) * wr(q);
  for i1=1:nphi
    for i2=1:nphi
      C(:,i1,i2) = C(:,i1,i2) + phi(i1,q)*(wx.*dphi(:,i2,1,q)+wy.*dphi(:,i2,2,q)).*fac;
    end
  end  
end