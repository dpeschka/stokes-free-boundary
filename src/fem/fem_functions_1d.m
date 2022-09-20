%
% Finite element function definition in 1D
%
% function element = fem_functions_1d(type,quadr)
%
% "dim"        reference dimension of finite element
% "type"       order(ing) of elements
% "quadr"      integration scheme to be used
%
% returns struct "element" which contains
%
% name       string name of the finite element
% type       type number of the finite element
% ndim       spatial dimension of reference element
% nphi       number of basis functions
% dof        assignment of basis functions of vertices,edges,...
% ndofs      assignment of #(basis functions) for global counting to
%            1=points,2=edges,3=triangles,4=tetrahedrons
%
% phi        basis functions at integration points
% gradphi    gradients of basis functions at integration points
% local_mass local mass matrix
%
% mr         number of integration points
% wr         integration weights
% quadr      integration scheme struct

function element = fem_functions_1d(type,quadr)

%
% Reference edge:
%
% 1---3---2
%

ndim = 1;
mr   = quadr.mr;
wr   = quadr.wr;
xq   = quadr.xq;

if quadr.ndim ~= 1
    error('Provided integration not compatible with element dimension!')
end


if type == 0
    
    name  = 'P0';
    nphi  = 1;
    dof   = 3;
    ndofs  = [0 1 0 0];
    
    coeff = [1];
    
    fphi  = @(x,c) c(1)+0*x;
    dphix = @(x,c) 0*x;
    
elseif type == 1
    
    name  = 'P1';
    nphi  = 2;
    dof   = [1 2];
    ndofs  = [1 0 0 0];
    
    coeff = [1 -1;...
             0 +1];
    
    fphi  = @(x,c) c(1)+c(2)*x;
    dphix = @(x,c) c(2)+0*x;
    
elseif type == 2
    
    name  = 'P2';
    nphi  = 3;
    dof   = [1 2 3];
    ndofs  = [1 1 0 0];
    
    coeff = [1 -3 +2;...
             0 -1 +2;...
             0 +4 -4;];
         
    fphi =@(x,c) c(1)+c(2)*x+c(3)*x.^2;
    dphix=@(x,c) c(2)+2*c(3).*x;
    
else
    
    error('Finite element function not found or implemented')
    
end

phi     = zeros(nphi,mr);
gradphi = zeros(nphi,ndim,mr);
for i=1:nphi
    for q=1:mr
        phi(i,q)=fphi(xq(q),coeff(i,:));
        gradphi(i,1,q)=dphix(xq(q),coeff(i,:));
    end
end

% computation of local mass matrix on reference triangle
local_mass = zeros(nphi,nphi);
for i=1:nphi
    for j=1:nphi
        for q=1:mr
            local_mass(i,j) = local_mass(i,j) + phi(i,q)*phi(j,q)*wr(q);
        end
    end
end

% store data in struct
element            = struct();

element.name       = name;
element.type       = type;
element.ndim       = ndim;
element.nphi       = nphi;
element.dof        = dof(:);
element.ndofs      = ndofs(:);

element.phi        = phi;
element.gradphi    = gradphi;
element.local_mass = local_mass;

element.mr         = quadr.mr;
element.wr         = quadr.wr;
element.quadr      = quadr;