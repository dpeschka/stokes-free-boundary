function [xx1,yy1,e2p1,np,ne] = refinemesh(xx,yy,e2p) % ,idp)

np   = max(e2p(:));
eg   = sparse(np,np);

% recompute book-keeping 
e2 = zeros(size(e2p,1),3);
cyc  = [1 2;2 3;3 1];
edge = 0;
e1   = [];
%% idp1 = idp;
for k=1:size(e2p,1)
    for c=1:3
        
        i1 = e2p(k,cyc(c,1));
        i2 = e2p(k,cyc(c,2));
        
        if full(eg(i1,i2))==0 % edge not yet counted ?
            edge = edge + 1;
            eg(i1,i2)=edge;
            eg(i2,i1)=edge;
            e2(k,c)=edge;
            e1 = [e1;i1 i2];
        else
            e2(k,c)=eg(i1,i2);
        end
        
    end   
end

% generate new mesh
xx1 = [xx;mean(xx(e1),2)];
yy1 = [yy;mean(yy(e1),2)];
%% idp1 = [idp1;idp(e1(:,1)).*idp(e1(:,2))];

em = e2 + np;

e2p1 = [e2p(:,1) em(:,1) em(:,3);...
        em(:,1) e2p(:,2) em(:,2);...
        em(:,2) e2p(:,3) em(:,3);...
        em(:,1) em(:,2) em(:,3)];
    
ne = size(e2p1,1);
np = size(xx1,1);

end