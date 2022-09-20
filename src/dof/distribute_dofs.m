% Computes global indices for finite element FE
% where dofs are counted as in dofmap

function [ ii,jj ] = distribute_dofs(dofmap,FE)

nphi     = FE.nphi;
nelement = size(dofmap,1);

ii = zeros(nelement,nphi,nphi);
jj = zeros(nelement,nphi,nphi);

for i1=1:nphi
    for i2=1:nphi
        ii(:,i1,i2) = dofmap(:,i1);
        jj(:,i1,i2) = dofmap(:,i2);
    end
end