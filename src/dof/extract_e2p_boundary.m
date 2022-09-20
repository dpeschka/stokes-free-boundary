function [e2p_bd,nedge_bd] = extract_e2p_boundary(e2p)

nelement = size(e2p,1);

% create edge infrastructure
e0 = [e2p(:,1:2);e2p(:,2:3);e2p(:,[3 1])];
z0 = [ones(nelement,1)*1;ones(nelement,1)*2;ones(nelement,1)*3];
k0 = 1:nelement;k0=[k0';k0';k0'];

% sorting to tell interior and boundary edges
sel           = e0(:,1)>e0(:,2);
e0(sel,[1 2]) = e0(sel,[2 1]);

[e0,I] = sortrows(e0);
k0 = k0(I);
z0 = z0(I);

interior = false(3*nelement,1);

sel = find(all(e0(1:end-1,:) == e0(2:end,:),2));
interior(sel)   = true;
interior(sel+1) = true;

k00 = k0(~interior);
z00 = z0(~interior);

nedge_bd = size(k00,1);

e2p_bd = zeros(nedge_bd,3);

e2p_bd(z00==1,1:3) = e2p(k00(z00==1),[1 2 4]);
e2p_bd(z00==2,1:3) = e2p(k00(z00==2),[2 3 5]);
e2p_bd(z00==3,1:3) = e2p(k00(z00==3),[3 1 6]);