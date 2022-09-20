ix1 = find(dm>0);
ix2 = dm(ix1);

xt = x(ix1);
yt = y(ix1);


%%
ubx = ux(ix1);
uby = uy(ix1);

ux0 = mean(ubx);

tx = Mtfbb\(Stfbb*xt);
ty = Mtfbb\(Stfbb*yt);
tn = sqrt(tx.^2+ty.^2);

tx = tx ./ tn;
ty = ty ./ tn;
nbx =  ty;
nby = -tx;

nbx(yt<1d-3)=sign(nbx(yt<1d-3));
nby(yt<1d-3)=0;

ix2 = find(y<1d-3);

n1=length(ix1);
n2=length(ix2);

ii_bd1 = [1:n1 1:n1];
jj_bd1 = [ix1 ndof_u+ix1];
aa_bd1 = [nbx' nby'];

ii_bd2 = [1:n2];
jj_bd2 = [ndof_u+ix2];
aa_bd2 = ones(1,n2);

bd1 = sparse(ii_bd1,jj_bd1,aa_bd1,n1,2*ndof_u);
bd2 = sparse(ii_bd2,jj_bd2,aa_bd2,n2,2*ndof_u);


ZZ = sparse(ndof_u,ndof_u);
A = [Sxx+Syy      ZZ;...
     ZZ      Sxx+Syy;];

  
Aext = [ A    bd1' bd2';...
        bd1   sparse(n1,n1) sparse(n1,n2) ;...
        bd2   sparse(n2,n1) sparse(n2,n2)];
      
rhs = [zeros(2*ndof_u,1);nbx.*ubx+nby.*uby;zeros(n2,1)];

uext = Aext\rhs;
  
uext_x = uext(1:ndof_u);
uext_y = uext(ndof_u+(1:ndof_u));

x = x + tau*uext_x;
y = y + tau*uext_y;