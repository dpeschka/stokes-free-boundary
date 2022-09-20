% Compute the determinant of the Jacobian on every triangle.

function [edet,dFinv]=vec_transformation_2d(e2p,x,y,map_element)

mr      = map_element.mr;
gradphi = map_element.gradphi;
nphi    = map_element.nphi;
ndim    = map_element.ndim;
nspace  = 2; % spatial dimension = 2 might differ from fe dimension 
% e.g. if we compute transformations in 2d for boundary curves then 
% nspace = 1 and ndim = 1

nelement = size(e2p,1);

dFinv = zeros(nelement,ndim,ndim,mr); % only tangential gradients, i.e. 
% ndim,dim instead of ndim,nspace or similar

edet  = zeros(nelement,mr);
dF    = zeros(nelement,nspace,ndim,mr);

for q=1:mr   

    for r=1:ndim
        for i=1:nphi
            dF(:,1,r,q) = dF(:,1,r,q) + x(e2p(:,i))*gradphi(i,r,q);
            dF(:,2,r,q) = dF(:,2,r,q) + y(e2p(:,i))*gradphi(i,r,q);
        end
    end
    
    if ndim==1
    
      edet(:,q)  = sqrt(dF(:,1,1,q).^2 + dF(:,2,1,q).^2);
      dFinv(:,1,1,q) = 1./ edet(:,q);
    
    elseif ndim==2
      
      edet(:,q)  = dF(:,1,1,q).*dF(:,2,2,q)-dF(:,1,2,q).*dF(:,2,1,q);
      dFinv(:,1,1,q) = +dF(:,2,2,q) ./ edet(:,q);
      dFinv(:,2,2,q) = +dF(:,1,1,q) ./ edet(:,q);
      dFinv(:,1,2,q) = -dF(:,1,2,q) ./ edet(:,q);
      dFinv(:,2,1,q) = -dF(:,2,1,q) ./ edet(:,q);
      
    else
      
      error('Dimension not supported in transformation')
      
    end
    
end