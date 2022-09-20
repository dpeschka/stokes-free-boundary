% dof handler constructor
function [dofmap,ndof,dof_mapper]  = generate_dofs(e2p,FE,reduce)

%
% dof = 0+                (1..npoint) ---> vertex    0D--3D
% dof = ndim+             (1..nedge ) ---> edges     1D--3D
% dof = ndim+nedge+       (1..ntria ) ---> triangles 2D--3D
% dof = ndim+nedge+ntria  (1..ntetr ) ---> tetrahedrons  3D
%
% e2p size
%
% 0D  ---> nelement,1
% 1D  ---> nelement,2
% 2D  ---> nelement,6
% 3D  ---> nelement,14

% Implement: 
%    * Multiple dofs per elementary element (vertex, edge, face, ...)
%    * Region/boundary selector

e2p_sizes = [1 2 6 14]+1;

ndim = FE.ndim;
dof  = FE.dof;

% supported dimensions check
assert(ndim==0|ndim==1|ndim==2,'assert: ndim not yet supported.')

% supported e2p sizes check, for smaller sizes on could try to upgrade
% the size by computing the remaining dof indices
assert(size(e2p,2)==e2p_sizes(ndim+1),'assert: e2p size not yet supported.')

% check if dofs are assigned multiple times
assert(length(unique(dof))==length(dof),'assert: multiple assigned dofs not yet supported.')

nelement = size(e2p,1);
maxdof   = max(e2p(:));

% extend dofmap by element enumeration starting at latest dof
dofmap = e2p(:,dof);

% counts 1. points, 2. edges, 3. triangles, 4. tetrahedrons
switch ndim
     case 0
         dofmap_local = 1;
     case 1
         dofmap_local = [1 1 2];
     case 2
         dofmap_local = [1 1 1 2 2 2 3];
     case 3
         dofmap_local = [1 1 1 1 2 2 2 2 2 2 3 3 3 3 4];
     otherwise
         error('strange dimension');
end

for i=1:4
    countlist = unique(e2p(:,dofmap_local(1:end)==i));
    ncount(i) = length(countlist);
end
ndofs = sum(ncount.*FE.ndofs);

if reduce
    dof_reduce = sort(unique(dofmap(:)));
    ndof       = length(dof_reduce);
    dof_mapper = zeros(maxdof,1);
    dof_mapper(dof_reduce) = 1:ndof;  
    dofmap     = dof_mapper(dofmap);
else
    ndof = maxdof;
    dof_mapper=1:maxdof;
end

end