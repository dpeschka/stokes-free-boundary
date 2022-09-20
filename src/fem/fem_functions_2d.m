%
% Finite element function definition in 2D
%
% function element = fem_functions_2d(type,quadr)
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

function element = fem_functions_2d(type,quadr)

%
% Reference triangle:
%  3
%  | \
%  |  \
%  6 7 5
%  |    \
%  1--4--2
%

ndim = 2;
mr   = quadr.mr;
wr   = quadr.wr;
xq   = quadr.xq;

if quadr.ndim ~= 2
    error('Provided integration not compatible with element dimension!')
end


% -------------------------------------------------------------------------
% definite finite element functions
if type == 0
    
    name  = 'P0(2D)';
    nphi  = 1;
    dof   = 7;
    
    ndofs  = [0 0 1 0];
    
    coeff = [1];
    
    fphi  =@(x,y,c) c(1)+0*x;
    dphix =@(x,y,c) 0*x;
    dphiy =@(x,y,c) 0*y;
    
elseif type == 1
    
    name  = 'P1(2D)';
    nphi  = 3;
    dof   = [1 2 3];
    ndofs  = [1 0 0 0];
    
    coeff = [1 -1 -1;...
        0  1  0;...
        0  0  1];
    
    fphi =@(x,y,c) c(1)+c(2)*x+c(3)*y;
    dphix=@(x,y,c) c(2)+0*x;
    dphiy=@(x,y,c) c(3)+0*x;
    
elseif type == 2
    
    name  = 'P2(2D)';
    nphi  = 6;
    dof   = [1 2 3 4 5 6];
    ndofs  = [1 1 0 0];
    
    coeff = [+1 -3 -3 +2 +4 +2;...
        +0 -1 +0 +2 +0 +0;...
        +0 +0 -1 +0 +0 +2;...
        +0 +4 +0 -4 -4 +0;...
        +0 +0 +0 +0 +4 +0;...
        +0 +0 +4 +0 -4 -4];
    
    fphi =@(x,y,c) c(1)+c(2)*x+c(3)*y+c(4)*x^2+c(5)*x*y+c(6)*y^2;
    dphix=@(x,y,c) c(2)+2*c(4).*x+c(5).*y;
    dphiy=@(x,y,c) c(3)+2*c(6).*y+c(5).*x;
    
elseif type == 3
    
    name  = 'P2-bubble(2D)';
    nphi  = 7;
    dof   = [1 2 3 4 5 6 7];
    ndofs  = [1 1 1 0];
       
    coeff = [1 -3 -3  2  4  2   3;...
             0 -1  0  2  0  0   3;...
             0  0 -1  0  0  2   3;...
             0  4  0 -4 -4  0 -12;...
             0  0  0  0  4  0 -12;...
             0  0  4  0 -4 -4 -12;...
             0  0  0  0  0  0  27];
    
    fphi =@(x,y,c) c(1)+c(2)*x+c(3)*y+c(4)*x^2+c(5)*x*y+c(6)*y^2+c(7)*x*y*(1-x-y);
    dphix=@(x,y,c) c(2)+2*c(4).*x+c(5).*y+c(7)*y*(1-y-2*x);
    dphiy=@(x,y,c) c(3)+2*c(6).*y+c(5).*x+c(7)*x*(1-x-2*y);
    
else
    
    error('Finite element function not found or implemented')
    
end

phi     = zeros(nphi,mr);
gradphi = zeros(nphi,ndim,mr);
for i=1:nphi
    for q=1:mr
        phi(i,q)=fphi(xq(1,q),xq(2,q),coeff(i,:));
        gradphi(i,1,q)=dphix(xq(1,q),xq(2,q),coeff(i,:));
        gradphi(i,2,q)=dphiy(xq(1,q),xq(2,q),coeff(i,:));
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