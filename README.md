# MATLAB solver for Stokes free boundary problem

Description: MATLAB code that solves a Stokes free boundary problem with dynamic contact angle using isoparametric finite elements in 2D

It solves the following PDE problem:

Conservation of momentum and mass in $\Omega(t)$:

$$
-\nabla\cdot\boldsymbol{\sigma} = \boldsymbol{f}\quad\text{and}\quad \nabla\cdot\boldsymbol{u}=0
$$

Navier-slip on $\Gamma_{\mathrm{s}\ell}(t)$:

$$
\boldsymbol{t}\cdot\boldsymbol{\sigma}\boldsymbol{\nu}=-\mu_\Gamma\boldsymbol{t}\cdot\boldsymbol{u}
$$

Capillary forces on $\Gamma_{\ell}(t)$:

$$
\boldsymbol{\sigma}\boldsymbol{\nu}=\gamma\kappa\boldsymbol{\nu}
$$

Dynamic contact angle on $\Lambda(t)$:

$$
\mu_\Lambda\dot{\mathrm{x}}=f_\Lambda
$$

with the standard Cauchy stress $\boldsymbol{\sigma}=-p\mathbb{I}+\mu(\nabla\boldsymbol{u}+\nabla\boldsymbol{u}^\top)$, mean curvature $\kappa$ and outer normal vector field $\boldsymbol{\nu}$. Uncompensated Young stress is encoded in $f_\Lambda$ and dissipation parameters are $\mu,\mu_\Gamma,\mu_\Lambda$. Solutions are droplets as shown below, where these dissipation parameters impact bulk viscosity, Navier slip and contact line dissipation, respectively.

<img src="src/drop.png">
