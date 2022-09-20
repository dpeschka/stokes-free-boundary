% Computes global indices for finite element FE
% where dofs are counted as in dofmap

function [ ii,jj ] = distribute_dofs_mixed(dofmap1,FE1,dofmap2,FE2)

nphi1     = FE1.nphi;
nphi2     = FE2.nphi;

nelement = size(dofmap1,1);

ii = zeros(nelement,nphi1,nphi2);
jj = zeros(nelement,nphi1,nphi2);

for i1=1:nphi1
    for i2=1:nphi2
        ii(:,i1,i2) = dofmap1(:,i1);
        jj(:,i1,i2) = dofmap2(:,i2);
    end
end