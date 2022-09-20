function [xnn,ynn,e2pnn]=uniformrefine(x,y,e2p)
npoint = max(max(e2p(:,1:3)));
[xn,yn,e2pn,~,~] = refinemesh(x(1:npoint),y(1:npoint),e2p(:,1:3));
[xnn,ynn,e2pnn] = dofs_p1_to_p2(xn,yn,e2pn);

maxdof = max(e2pnn(:));

e2pnn = [e2pnn maxdof+(1:size(e2pnn,1))'];

[e2pbnn,~] = extract_e2p_boundary(e2pnn);

xnn = [xnn;mean(xnn(e2p(:,1:3)),2)];
ynn = [ynn;mean(ynn(e2p(:,1:3)),2)];


selb = unique(e2pbnn(:));

r = sqrt(xnn(selb).^2+ynn(selb).^2);

xnn(selb) = xnn(selb) ./ r;
ynn(selb) = ynn(selb) ./ r;

