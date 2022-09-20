%
% Mixed FE approach for an incompressible Stokes problem with free boundary
% and moving (dynamic) contact lines
%
clear all

addpath('fem/')
addpath('dof/')
addpath('io/')

mprint = 0;
t      = 0;
tau    = 0.005;

mu     = 1.0;
beta   = 0.1;
gamma  = 1.0;
Bond   = 0.5;

sigma_fb = -1; % surface tensions
sigma_bb =  0; % surface tensions

%%
fem_init
for it=1:900
    t = t + tau;
    fem_solveflow
    fem_ALE
    if mod(it,10)==0
        um = mean(ux(:));
        
        trisurf(e2p(:,1:3),x(1:npoint),y(1:npoint),0*p,ux(1:npoint))
        hold on
        quiver(x,y,ux-um,uy,'r')
        hold off
        
        shading faceted
        view(2)
        axis equal
        box on
        xlim([-1.5 3.5])
        ylim([0 1.3])
        drawnow
    end
end

