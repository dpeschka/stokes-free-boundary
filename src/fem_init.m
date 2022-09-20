% FE init for Taylor-Hood isoparametric elements

ndim     =    2;   % dimension
fe_deg_u =    2;   % degree of fe-space velocity
fe_deg_p =    1;   % degree of fe-space pressure
fname    = 'mesh/domain_droplet3x'; % filename

% load mesh and generate dofs
[x,y,npoint,nelement,e2p,idp,ide] = readtriamesh(fname);
[e2pb,~] = extract_e2p_boundary(e2p);

selfb = y(e2pb(:,3))>1d-2; 
selbb = ~selfb;         

e2pbb = e2pb(selbb,:);% free boundary
e2pfb = e2pb(selfb,:);% part dirichlet boundary

selfb = unique(e2pfb(:));
xfb = x(selfb);
yfb = y(selfb);
rfb = sqrt(xfb.^2+yfb.^2);

x(selfb) = xfb ./ rfb;
y(selfb) = yfb ./ rfb;

% construct quadrature rules & FE spaces
qr1D = fem_quadrature_rule( 1 , 3);
qr2D = fem_quadrature_rule( 2 , 5);

mapFE  = fem_functions_2d(2, qr2D);
mapFEb = fem_functions_1d(2, qr1D);

FE_u   = fem_functions_2d(fe_deg_u, qr2D);
FE_p   = fem_functions_2d(fe_deg_p, qr2D);
FE_ub  = fem_functions_1d(fe_deg_u, qr1D);

[dof_u    , ndof_u   ] = generate_dofs(e2p   , FE_u  ,  true);
[dof_p    , ndof_p   ] = generate_dofs(e2p   , FE_p  ,  true);
[dof_ubb  , ndof_ubb ] = generate_dofs(e2pbb , FE_ub , false);
[dof_ufb  , ndof_ufb ] = generate_dofs(e2pfb , FE_ub , false);

[dof_ufbb , ndof_ufbb, dm] = generate_dofs(e2pfb , FE_ub ,  true);
 
x = x(1:ndof_u);
y = y(1:ndof_u);

idp = 0*x;
idp(y<1d-4) = 1;