% MATLAB code that solves a Stokes free boundary problem with dynamic contact angle 
% using isoparametric finite elements in 2D. This packages provides a self-containted 
% finite element framework using P2/P1 Taylor finite elements and P2 mapping of mesh 
% description. A more detailed explanation of the model and its weak formulation can be found, 
% for example, in the paper
%
% Resolving the microscopic hydrodynamics at the moving contact line
% by Amal K. Giri, Paolo Malgaretti, Dirk Peschka, Marcello Sega
% published in Physical Review Fluids (2022).
%
% Author: Dirk Peschka
%

clear all

addpath('fem/')
addpath('dof/')
addpath('io/')

t      = 0;      % time
tau    = 0.005;  % time step size

mu     = 1.0;    % bulk viscosity
beta   = 0.1;    % navier slip
gamma  = 1.0;    % contact line dissipation
Bond   = 0.5;    % bond number (acceleration)

% surface tensions, where sigma_fb is the tension of the free boundary 
% and sigma_bb is the difference of solid/liquid and solid/gas surface
% tension and determines the equilibrium contact angle. If sigma_bb is 
% set zero, then the equilibrium contact angle is 90 degrees.
sigma_fb = -1;
sigma_bb =  0;

% main time loop
fem_init
for it=1:900
    t = t + tau;

    % solve stokes flow for pressure p and flow field u=(ux,uy)
    fem_solveflow

    % perform mesh update using ALE mesh motion strategy
    fem_ALE

    % plotting
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

