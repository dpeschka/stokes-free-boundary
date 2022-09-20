% Vectorized assembly of convection matrix
% S = \int grad(phi_i) * grad(phi_j) dx

function [Dx,Dy]=vec_localderivatives_weighted(edet,dFinv,FE,dofmap,func,weight)

nphi     = FE.nphi;
phi      = FE.phi;
mr       = FE.mr;
wr       = FE.wr;
ndim     = FE.ndim;
gradphi  = FE.gradphi;

nelement = size(edet,1);

Dx       = zeros(nelement,nphi,nphi);
Dy       = zeros(nelement,nphi,nphi);
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
  
  w = zeros(nelement,1);  
  for i=1:nphi
    w = w + weight(dofmap(:,i))*phi(i,q);
  end
  % compute matrix elements 
  fac = func(w) .* edet(:,q) * wr(q);  
  % fac = edet(:,q) * wr(q);
  for i1=1:nphi
    for i2=1:nphi
      Dx(:,i1,i2) = Dx(:,i1,i2) + phi(i1,q).*dphi(:,i2,1,q).*fac;
      Dy(:,i1,i2) = Dy(:,i1,i2) + phi(i1,q).*dphi(:,i2,2,q).*fac;
    end
  end  
end
