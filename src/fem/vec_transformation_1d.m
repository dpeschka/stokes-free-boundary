% Compute the determinant of the Jacobian on every triangle.

function [edet,dFinv]=vec_transformation_1d(e2p,x,y,map_element)

mr      = map_element.mr;
gradphi = map_element.gradphi;
nphi    = map_element.nphi;
ndim    = map_element.ndim;

nelement = size(e2p,1);

dFinv = zeros(nelement,ndim,ndim,mr);
edet  = zeros(nelement,mr);
dF    = zeros(nelement,ndim,ndim,mr);

for q=1:mr   

    for r=1:ndim
        for i=1:nphi
            dF(:,1,r,q) = dF(:,1,r,q) + x(e2p(:,i))*gradphi(i,r,q);
            dF(:,2,r,q) = dF(:,2,r,q) + y(e2p(:,i))*gradphi(i,r,q);
        end
    end
    
    if dim==1
    elseif dim==2
      edet(:,q)  = dF(:,1,1,q).*dF(:,2,2,q)-dF(:,1,2,q).*dF(:,2,1,q);
      dFinv(:,1,1,q) = +dF(:,2,2,q) ./ edet(:,q);
      dFinv(:,2,2,q) = +dF(:,1,1,q) ./ edet(:,q);
      dFinv(:,1,2,q) = -dF(:,1,2,q) ./ edet(:,q);
      dFinv(:,2,1,q) = -dF(:,2,1,q) ./ edet(:,q);
    end
    
end