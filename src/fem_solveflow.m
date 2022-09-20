% SOLVES a SINGLE STEP of STOKES PROBLEM

% construct map
[edet  ,dFinv  ]=vec_transformation_2d(e2p  ,x,y,mapFE ); % bulk
[edetbb,dFinvbb]=vec_transformation_2d(e2pbb,x,y,mapFEb); % boundary 1
[edetfb,dFinvfb]=vec_transformation_2d(e2pfb,x,y,mapFEb); % boundary 2

%%  
% assembly for matrices, indices etc.

% stiffness
aaxx     = vec_localstiff_ij  (edet,dFinv,FE_u,1,1);
aaxy     = vec_localstiff_ij  (edet,dFinv,FE_u,1,2);
aayx     = vec_localstiff_ij  (edet,dFinv,FE_u,2,1);
aayy     = vec_localstiff_ij  (edet,dFinv,FE_u,2,2);

% mass for pressure and velocity
mm_p     = vec_localmass      (edet  ,FE_p );
mm_u     = vec_localmass      (edet  ,FE_u );

% boundary matrices
mm_ubb   = vec_localmass      (edetbb,FE_ub);
ss_ubb   = vec_localstiff     (edetbb,dFinvbb,FE_ub);
ss_ufb   = vec_localstiff     (edetfb,dFinvfb,FE_ub);

[stfbb]  = vec_localconv_mixed_1d(edetfb,dFinvfb,FE_ub,FE_ub);
[mtfbb]  = vec_localmass         (edetfb,FE_ub);

% convection bulk
[sx  ,sy  ] = vec_localconv_mixed(edet  ,dFinv  ,FE_p ,FE_u );

% bulk dof indices velocity and pressure
[ii_u  ,jj_u  ] = distribute_dofs (dof_u,FE_u);
[ii_p  ,jj_p  ] = distribute_dofs (dof_p,FE_p);

% boundary dofs indices
[ii_bb ,jj_bb ] = distribute_dofs (dof_ubb ,FE_ub);
[ii_fb ,jj_fb ] = distribute_dofs (dof_ufb ,FE_ub);
[ii_fbb,jj_fbb] = distribute_dofs (dof_ufbb,FE_ub);

% mixed dofs
[ii_pu ,jj_pu]  = distribute_dofs_mixed (dof_p,FE_p,dof_u,FE_u);

%%  
% build sparse matrices
Mbb  = sparse(ii_bb(:),jj_bb(:),mm_ubb(:),ndof_u,ndof_u);
Sbb  = sparse(ii_bb(:),jj_bb(:),ss_ubb(:),ndof_u,ndof_u);
Sfb  = sparse(ii_fb(:),jj_fb(:),ss_ufb(:),ndof_u,ndof_u);

Stfbb = sparse(ii_fbb(:),jj_fbb(:),stfbb(:),ndof_ufbb,ndof_ufbb);
Mtfbb = sparse(ii_fbb(:),jj_fbb(:),mtfbb(:),ndof_ufbb,ndof_ufbb);

Sxx = sparse(ii_u(:),jj_u(:),aaxx(:),ndof_u,ndof_u);
Sxy = sparse(ii_u(:),jj_u(:),aaxy(:),ndof_u,ndof_u);
Syx = sparse(ii_u(:),jj_u(:),aayx(:),ndof_u,ndof_u);
Syy = sparse(ii_u(:),jj_u(:),aayy(:),ndof_u,ndof_u);

M   = sparse(ii_u(:),jj_u(:),mm_u(:),ndof_u,ndof_u);
Mp  = sparse(ii_p(:),jj_p(:),mm_p(:),ndof_p,ndof_p);
Bx  = sparse(ii_pu(:),jj_pu(:),sx(:),ndof_p,ndof_u);
By  = sparse(ii_pu(:),jj_pu(:),sy(:),ndof_p,ndof_u);

Zp  = sparse(ndof_p,ndof_p);
Zup  = sparse(ndof_p,ndof_u);
SS  = sigma_fb*Sfb + sigma_bb*Sbb;

%%  
% build mixed problem from S including constraints encoded in B matrices
dof_cl = intersect(dof_ubb(:),dof_ufb(:));
Mcl = sparse(dof_cl(:),dof_cl(:),[1;1],ndof_u,ndof_u);
A   = [2*mu*Sxx+mu*Syy+beta*Mbb+gamma*Mcl+tau/2*SS       mu*Syx Bx' ;...
     mu*Sxy       mu*Sxx+2*mu*Syy+beta*Mbb+tau/2*SS By' ;...
     Bx        By        Zp  ];

% extra B matrices for Dirichlet bcs
jbd = ndof_u+find(idp>0);
Nbd = length(jbd);
ibd = (1:Nbd)';
Bbd = sparse(ibd(:),jbd(:),ones(Nbd,1),Nbd,2*ndof_u+ndof_p+1*0);
Zbd = sparse(Nbd,Nbd);
Abd = [A Bbd';Bbd Zbd];

%%  
% solve problem with forces fx,fy and disect the solution components
fx=+Bond+0*x;
fy=+0.0+0*x;
u = Abd\[M*fx+SS*x;M*fy+SS*y;zeros(ndof_p+Nbd+1*0,1)];

ux = u(1:ndof_u);
uy = u((ndof_u+1):2*ndof_u);
p  = u((2*ndof_u)+(1:ndof_p));