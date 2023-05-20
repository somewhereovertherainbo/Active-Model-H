# Active-Model-H
AMH is simulated using Pseudo-spectral method 

The equations of motion for Active Model H are:

$$
\dot\phi = -\nabla\cdot(J+J^\Lambda)
$$

$J^\Lambda$ is the current arising from 'Thermal Noise'

J is the current arising from Stresses (forces) created in the fluid due to the concentration gradients in different kinds of substances that make up the fluid.


$$
\dot{\phi} =-\boldsymbol{\nabla}\cdot \left(\boldsymbol{J}^\phi + \bm J^{\bm v}
+\boldsymbol{J}^{\Lambda}\right).
$$

where the currents are
$$

\boldsymbol{J}^\phi=-M^{\phi}\boldsymbol{\nabla}\mu^{\phi},\qquad\mu^{\phi}=\frac{\delta\mathcal{F}}{\delta\phi},
\qquad \bm J^{\bm v} = \bm v \phi,
\qquad\boldsymbol{J}^{\Lambda}=\sqrt{2D^{\phi}M^{\phi}}\boldsymbol{\Lambda}

$$

$\boldsymbol{\Lambda}$ is a zero-mean, unit-variance Gaussian white noise. 


The equilibrium $\mathcal{F}$ is the Landau-Ginzburg
free energy functional in the previous chapter\ref{ch2}.

So, the final dynamical equation in the variable $\phi$ is 

$$
\dot{\phi} =-\nabla^2 \left[a\phi + b \phi^3 - \kappa \nabla^2 \phi \right] - \mathbf{v}\cdot\nabla\phi - \sqrt{2D^{\phi}M^{\phi}}\nabla\cdot\Lambda
$$
