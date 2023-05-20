# Active-Model-H
AMH is simulated using Pseudo-spectral method 

The equations of motion for Active Model H are:

$$
\dot{\phi} =-\boldsymbol{\nabla}\cdot \left(J^\phi +  J^{ v}+J^{\Lambda}\right).
$$

$J^\Lambda$ is the current arising from 'Thermal Noise'

$J^vv$ is the current arising from Stresses (forces) created in the fluid due to the concentration gradients in different kinds of substances that make up the fluid.

$J^\phi$ is the current that minimizes the Free energy of the system. 

the currents are written as follows:

$$
\boldsymbol{J}^\phi=-M^{\phi}\boldsymbol{\nabla}\mu^{\phi},\qquad\mu^{\phi}=\frac{\delta\mathcal{F}}{\delta\phi},\qquad  J^{ v} =  v \phi,
\qquad J^{\Lambda}=\sqrt{2D^{\phi}M^{\phi}}\boldsymbol{\Lambda}
$$

$\boldsymbol{\Lambda}$ is a zero-mean, unit-variance Gaussian white noise. 


The equilibrium $\mathcal{F}$ is the Landau-Ginzburg
free energy functional given by:

$$
\mathcal{F}[\phi]=\int\left(\frac{a}{2}\phi^{2}+\frac{b}{4}\phi^{4}+\frac{\kappa}{2}(\boldsymbol{\nabla}\phi)^{2}
%+\frac{D}{2}c^{2}-\beta\phi c
\right)d\boldsymbol{r} = \int f \ d\boldsymbol{r} 
$$

So, the final dynamical equation in the variable $\phi$ is 

$$
\dot{\phi} =-\nabla^2 \left[a\phi + b \phi^3 - \kappa \nabla^2 \phi \right] - \mathbf{v}\cdot\nabla\phi - \sqrt{2D^{\phi}M^{\phi}}\nabla\cdot\Lambda
$$

As we have fluid flow, we must look at the Navier Stokes equation to obtain a relation between the stresses and the flow in fluid. 

Navier-Stokes Equation:

$$
\rho(\dot{\mathbf{v}}+\mathbf{v} \cdot \nabla \mathbf{v})= \nabla\cdot\big(\sigma+\Sigma\big)
$$

where

$$
\sigma = -P I + \eta\left([\nabla\mathbf v]+[\nabla\mathbf v]^T\right);\quad \boldsymbol\Sigma = \tilde\kappa\textbf{S}; \quad \mathbf{S} \equiv(\boldsymbol{\nabla} \phi)(\boldsymbol{\nabla} \phi)-\frac{1}{d}|\boldsymbol{\nabla} \phi|^2 \mathbf{I}
$$

In the limit of low reynold's number, stokes equaiton becomes :

$$
\boldsymbol{\nabla} \cdot\big( \sigma + \Sigma + \hat{\mathbf{f}}\big) = 0
$$

Here, $\sigma$ is the viscous stress, $\Sigma$ is the stress due to the gradients in $\phi$ and $\hat{\mathbf{f}}$ is the thermal noise stress.

The simulation is done using a Class that simulated the above model, given the parameters like M, $\tilde\kappa$, $\eta$ etc.

# Numerical Solution using Pseudo-Spectral Method:

In Fourier space, the dynamical equations are of the form:

$$
\dot{\phi}_q=\alpha(q) \phi_q+\hat{N}_q
$$

Multiplying both sides by $\exp (-\alpha(k) t)$ gives

$$
\frac{d\left(\phi_q e^{-\alpha t}\right)}{dt}=\frac{d \psi}{dt}=e^{-\alpha t} \hat{N}_q
$$

Thus, we can solve for $\psi_q=\phi_q e^{-\alpha t}$ and then compute
  
$$
\phi_q(t+dt)=\psi_q(t+dt) e^{\alpha(t+dt)}=\left[\psi_q(t)+dt\left(e^{-\alpha t} \hat{N}_q\right)\right] e^{\alpha(t+dt)}
$$

This can be simplified to 
  
$$
\phi_q (t+dt) = \left[\phi_q (t) + dt\hat{N}_q   \right] e^{\alpha dt}
$$

The code for simulation and the results:

```
import numpy as np
import sys, time
import matplotlib.pyplot as plt
fft2  = np.fft.fft2
ifft2 = np.fft.ifft2
randn = np.random.randn

import scipy.fft
fft2 = scipy.fft.fft2
ifft2 = scipy.fft.ifft2


from matplotlib import rc
rc('text', usetex=True)
fSA=20
font = {'family' : 'normal',
         'weight' : 'bold',
         'size'   : 24}  
rc('font', **font)
```




