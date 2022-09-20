%
% Quadrature rules on simplices in 1D and 2D for different orders
%
% function quad_rule = fem_quadrature_rule(ndim,order)
%
% "ndim"    dimension of simplex integration rule
% "order"   order(ing) of integration rules
%
% returns struct "quad_rule" which contains
%
% name       string name of the quadrature rule
% ndim       spatial dimension of reference element
% order      order(ing) of quadrature rule
% mr         number of integration points
% wr         integration weights
% xq         location of integration points of reference simplex
% 

function quad_rule = fem_quadrature_rule(ndim,order)

scheme = [ndim order];

if all(scheme == [1 1])
    % https://en.wikipedia.org/wiki/Gaussian_quadrature
    name = 'Gauss, 1D, 1 point(s), 2nd order';
    mr   = 1;
    
    xq(1,1) = 1/2;
    wr(1)   = 1;
    
elseif all(scheme == [1 2])
    % https://en.wikipedia.org/wiki/Gaussian_quadrature
    name = 'Gauss, 1D, 2 point(s), 4th order';
    mr   = 2;
    ndim = 1;
    
    xq(1,1) = 1/2 - sqrt(1/12);
    xq(1,2) = 1/2 + sqrt(1/12);
    
    wr(1) = 1/2;
    wr(2) = 1/2;
    
elseif all(scheme == [1 3])
    % https://en.wikipedia.org/wiki/Gaussian_quadrature
    name = 'Gauss, 1D, 3 point(s), 6th order';
    mr   = 3;
    
    xq(1,1) = 1/2;
    xq(1,2) = 1/2 - sqrt(3/20);
    xq(1,3) = 1/2 + sqrt(3/20);
    
    wr(1) = 8/18;
    wr(2) = 5/18;
    wr(3) = 5/18;
    
elseif all(scheme == [1 4])
    % https://en.wikipedia.org/wiki/Gaussian_quadrature
    name = 'Gauss, 1D, 4 point(s), 8th order';
    mr   = 4;
    
    rt1=sqrt((15+2*sqrt(30))/35);
    rt2=sqrt((15-2*sqrt(30))/35);
    
    xq(1,1) = 1/2-rt1/2;
    xq(1,2) = 1/2+rt1/2;
    xq(1,3) = 1/2-rt2/2;
    xq(1,4) = 1/2+rt2/2;
    
    wr(1) = (1/2-sqrt(30)/36)/2;
    wr(2) = (1/2-sqrt(30)/36)/2;
    wr(3) = (1/2+sqrt(30)/36)/2;
    wr(4) = (1/2+sqrt(30)/36)/2;
    
 elseif all(scheme == [1 5])
    % https://en.wikipedia.org/wiki/Gaussian_quadrature
    name = 'Gauss, 1D, 5 point(s), 10th order';
    mr   = 5;
    
    rt1=1/3*sqrt(5-2*sqrt(10/7));
    rt2=1/3*sqrt(5+2*sqrt(10/7));
    
    xq(1,1) = 1/2;
    xq(1,2) = 1/2-rt1/2;
    xq(1,3) = 1/2+rt1/2;
    xq(1,4) = 1/2-rt2/2;
    xq(1,5) = 1/2+rt2/2;
    
    wr(1) = (128/225)/2;
    wr(2) = ((322+13*sqrt(70))/900)/2;
    wr(3) = ((322+13*sqrt(70))/900)/2;
    wr(4) = ((322-13*sqrt(70))/900)/2;
    wr(5) = ((322-13*sqrt(70))/900)/2;
    
elseif all(scheme == [2 1])
    
    % "Gaussian quadrature formulas for triangles" Cowper
    % "High degree efficient symmetrical Gaussian quadrature rules for the
    % ... triangle" Dunavant
    % "A compendium of FEM integration formulas for symbolic work" Felippa
    % "Quadrature on Simplices of Arbitrary Dimension" Walkington
    
    name = 'Gauss, 2D, 1 point(s), 1st order';
    mr   = 1;
    
    wr(1)  =1/2;
    
    xq = zeros(ndim,mr);
    xq(1,1)=1/3;
    xq(2,1)=1/3;
    
elseif all(scheme == [2 2])
    
    % "Gaussian quadrature formulas for triangles" Cowper
    % "High degree efficient symmetrical Gaussian quadrature rules for the
    % ... triangle" Dunavant
    % "A compendium of FEM integration formulas for symbolic work" Felippa
    % "Quadrature on Simplices of Arbitrary Dimension" Walkington
       
    name = 'Gauss, 2D, 3 point(s), 2nd order';
    mr   = 3;
    
    wr(1:3)  =1/6;
    
    xq = zeros(ndim,mr);
    xq(1,1)=1/6;
    xq(2,1)=1/6;

    xq(1,2)=2/3;
    xq(2,2)=1/6;
    
    xq(1,3)=1/6;
    xq(2,3)=2/3;
    
elseif all(scheme == [2 3])
    
    % "Gaussian quadrature formulas for triangles" Cowper
    % "High degree efficient symmetrical Gaussian quadrature rules for the
    % ... triangle" Dunavant
    % "A compendium of FEM integration formulas for symbolic work" Felippa
    % "Quadrature on Simplices of Arbitrary Dimension" Walkington
       
    name = 'Gauss, 2D, 4 point(s), 3rd order';
    mr   = 4;
    
    wr(1)    = -27/96;
    wr(2:4)  =  25/96;
    
    xq = zeros(ndim,mr);
    xq(1,1)=1/3;
    xq(2,1)=1/3;

    xq(1,2)=1/5;
    xq(2,2)=3/5;
    
    xq(1,3)=1/5;
    xq(2,3)=1/5;
    
    xq(1,4)=3/5;
    xq(2,4)=1/5;
        
elseif all(scheme == [2 4])
    
    % "Gaussian quadrature formulas for triangles" Cowper
    % "High degree efficient symmetrical Gaussian quadrature rules for the
    % ... triangle" Dunavant
    % "A compendium of FEM integration formulas for symbolic work" Felippa
    % "Quadrature on Simplices of Arbitrary Dimension" Walkington
       
    name = 'Gauss, 2D, 6 point(s), 4th order';
    mr   = 6;
    
    rt1 = (8-sqrt(10)+sqrt(38-44*sqrt(2/5)))/18;
    rt2 = (8-sqrt(10)-sqrt(38-44*sqrt(2/5)))/18;
    wt1 = (620+sqrt(213125-53320*sqrt(10)))/3720;
    wt2 = (620-sqrt(213125-53320*sqrt(10)))/3720;
    
    wr(1:3) = wt1/2;
    wr(4:6) = wt2/2;
    
    xq = zeros(ndim,mr);
    xq(1,1)=rt1;
    xq(2,1)=rt1;

    xq(1,2)=rt1;
    xq(2,2)=1-2*rt1;
    
    xq(1,3)=1-2*rt1;
    xq(2,3)=rt1;
    
    xq(1,4)=rt2;
    xq(2,4)=rt2;    
    
    xq(1,5)=rt2;
    xq(2,5)=1-2*rt2; 
    
    xq(1,6)=1-2*rt2;
    xq(2,6)=rt2; 
    
elseif all(scheme == [2 5])
    
    % "Gaussian quadrature formulas for triangles" Cowper
    % "High degree efficient symmetrical Gaussian quadrature rules for the
    % ... triangle" Dunavant
    % "A compendium of FEM integration formulas for symbolic work" Felippa
    % "Quadrature on Simplices of Arbitrary Dimension" Walkington
   
    name = 'Gauss, 2D, 7 point(s), 5th order';
    mr   = 7;
    
    rt=sqrt(15);
    wr(1)  =9/80;
    wr(2:4)=(155+rt)/2400;
    wr(5:7)=(155-rt)/2400;
    
    x1=1/3;
    x2=(6+  rt)/21;
    x3=(9-2*rt)/21;
    x4=(6-  rt)/21;
    x5=(9+2*rt)/21;
    
    xq = zeros(ndim,mr);
    
    xq(1,1)=x1;
    xq(2,1)=x1;
    
    
    xq(1,2)=x2;
    xq(2,2)=x3;
    
    xq(1,3)=x3;
    xq(2,3)=x2;
    
    
    xq(1,4)=x2;
    xq(2,4)=x2;
    
    
    xq(1,5)=x4;
    xq(2,5)=x5;
    
    xq(1,6)=x5;
    xq(2,6)=x4;
    
    xq(1,7)=x4;
    xq(2,7)=x4;
    
else
    
    error('Integration scheme not found or implemented')
    
end

quad_rule       = struct();
quad_rule.name  = name;
quad_rule.ndim  = ndim;
quad_rule.order = order;
quad_rule.mr    = mr;
quad_rule.wr    = wr;
quad_rule.xq    = xq;