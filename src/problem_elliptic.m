% FEM code DEMO:
% Solves 2D isoparametric elliptic problem 
%
% 1) DIRICHLET = TRUE
% \int_\Omega grad(u)*grad(v) dx = \int_\Omega f*v dx
%
% or
%
% 2) DIRICHLET = FALSE
% \int_\Omega grad(u)*grad(v) dx + \int_\Gamma u*v ds = \int_\Omega f*v dx
%
% with f=1 and \Omega the unit disc and \Gamma its boundary.
%
% Author: Dirk Peschka
%
addpath('fem/')
addpath('io/')
addpath('dof/')

DIRICHLET = false;

qint2d    = 5;
qint1d    = 3;

qr1D   = fem_quadrature_rule( 1, qint1d);
qr2D   = fem_quadrature_rule( 2, qint2d);

mapFE  = fem_functions_2d(2, qr2D);
mapFEb = fem_functions_1d(2, qr1D);

FE     = fem_functions_2d(2, qr2D);
FEb    = fem_functions_1d(2, qr1D);

[x,y,npoint,nelement,e2p,idp,ide] = readtriamesh('mesh/domain_disc');
[e2pb,~] = extract_e2p_boundary(e2p);

[dof_FE ,ndof_FE ] = generate_dofs(e2p ,FE ,true );
[dof_FEb,ndof_FEb] = generate_dofs(e2pb,FEb,false);

%%
% project to manifold
bd    = (idp==1);
r     = sqrt(x.^2 + y.^2);
x(bd) = 3*x(bd)./r(bd);
y(bd) = 3*y(bd)./r(bd);

[edet , dFinv ]=vec_transformation_2d(e2p ,x,y,mapFE );
[edetb, dFinvb]=vec_transformation_2d(e2pb,x,y,mapFEb);

% assembly
aa  = vec_localstiff(edet ,dFinv,FE );
mm  = vec_localmass (edet       ,FE );
mmb = vec_localmass (edetb      ,FEb);

[ii ,jj ] = distribute_dofs(dof_FE ,FE );
[iib,jjb] = distribute_dofs(dof_FEb,FEb);

A  = sparse(ii(:) ,jj(:) ,aa(:) ,ndof_FE,ndof_FE);
M  = sparse(ii(:) ,jj(:) ,mm(:) ,ndof_FE,ndof_FE);
Mb = sparse(iib(:),jjb(:),mmb(:),ndof_FE,ndof_FE);

rhs = M*ones(ndof_FE,1);

if DIRICHLET
  u      = zeros(ndof_FE,1);
  sel    = ~bd(1:ndof_FE);
  u(sel) = A(sel,sel)\rhs(sel);
  u0     = 0;
else
  A  = A + Mb;
  u  = A\rhs;
  u0 = 1.5;
end

id = outvtk_triP2_scalar('output/elliptic.vtk',x,y,e2p,u,'sol');

xx = x(1:npoint);
yy = y(1:npoint);
uu = u(1:npoint);

uex = 0.25*(9-xx.^2-yy.^2) + u0;

trisurf(e2p(:,1:3),xx,yy,uu-0*uex)
fprintf('error=%e\n',max(abs(uu-uex)))