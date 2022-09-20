function [xe,ye,dofs] = dofs_p1_to_p2(x,y,e2)

e2 = e2(:,1:3);

ne = size(e2,1);
np = max(e2(:));

allfaces  = [e2(:,[1 2]);e2(:,[2 3]);e2(:,[3 1])];
face_sort = unique(sort(allfaces,2),'rows');
nf        = length(face_sort);

enumeration = sparse(face_sort(:,1),face_sort(:,2),np+(1:nf),np,np);
enumeration = enumeration + enumeration';

dofs = e2;
xe   = [x;mean(x(face_sort),2)];
ye   = [y;mean(y(face_sort),2)];

dofs = zeros(ne,6);

for i=1:ne
  dofs(i,1:3) = e2(i,1:3);
  dofs(i,4) = enumeration(e2(i,1),e2(i,2));
  dofs(i,5) = enumeration(e2(i,2),e2(i,3));
  dofs(i,6) = enumeration(e2(i,3),e2(i,1));
end
  
  

