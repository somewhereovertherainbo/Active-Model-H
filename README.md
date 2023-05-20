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

```
class FPS():
    '''Class to simulate field theories using a PSEUDO-SPECTRAL TECHNIQUE WITH FILTERING'''
    def __init__(self, param):
        self.Nt = param['Nt'];    self.dt=param['dt'];    self.Nf = param['Nf']
        self.h  = param['h'];     self.a = param['a'];    self.b  = param['b']
        self.kp = param['kp'];    self.k_tilde= param['k_tilde'];  self.Ng= param['Ng']
        self.Df = 1j*param['Df']; Ng = self.Ng
        
        self.XX  = np.zeros((int(self.Nf+1), Ng*Ng)); 
        self.LL  = np.zeros((int(self.Nf+1), Ng*Ng)); 
        
        qf=np.fft.fftfreq(Ng)*(2*np.pi/self.h)
        self.qx, self.qy = np.meshgrid(qf, qf) 
        self.Dx, self.Dy = 1j*self.qx, 1j*self.qy

        self.q2 = self.qx*self.qx + self.qy*self.qy;    
        alpha_q = -self.q2*(self.a + self.kp*self.q2)*self.dt;   
        self.eA = np.exp(alpha_q);   self.q2b=self.q2*self.b 
        
        iq2=1/(self.q2+1e-16);  iq2[0,0]=0; self.iq2=iq2   
        self.iQ2x, self.iQ2y = self.qx*iq2, self.qy*iq2
        self.qmod = np.sqrt( (self.q2) )
         
        
    def integrate(self, u):
        '''  simulates the equation and plots it at different instants '''
        ii=0;  t=0;  dt=self.dt;    Df=self.Df; simC=(100/self.Nt)
        self.u=u
        for i in range(self.Nt):          
            self.rhs()

            if i%(int(self.Nt/self.Nf))==0:  
                self.XX[ii,:] = (np.real(ifft2(self.u))).flatten()
                #self.LL[ii,:] = self.getLength( )
                ii += 1   
                if ii%50==0:
                    print (int(simC*i), '% done...', end=' ')
        #print ('100% done.')  
            
        
    def rhs(self):
        '''
        integrator to simulate the dynamics of order paramter 
        using filtere-pseudo-spectral (FPS)
        '''       
        vx, vy = self.flow(self.u)
            
        ## evolve the order paramter using FPS
        uc  = ifft2(self.u);    
        N_u = -self.q2b*fft2(uc*uc*uc) - self.k_tilde*fft2(vx*self.u_x + vy*self.u_y)        

        self.u = (self.u + N_u*self.dt)*(self.eA) 
        return        
    
    
    def flow(self, u):
        ## compute all the stress tensors
        self.u_x = ifft2(self.Dx*u)        
        self.u_y = ifft2(self.Dy*u)
        sxx = 0.5*fft2( self.u_y*self.u_y-self.u_x*self.u_x )
        sxy =    -fft2( self.u_x*self.u_y )

        # compute flow
        fqx = self.Dx*sxx + self.Dy*sxy   
        fqy = self.Dx*sxy - self.Dy*sxx 
        fac = fqx*self.iQ2x + fqy*self.iQ2y
        
        vx = ifft2(self.iq2*fqx-fac*self.iQ2x) 
        vy = ifft2(self.iq2*fqy-fac*self.iQ2y)  
        return vx, vy
    
    
    def getLength(self):
        u1 = np.abs( (self.u) );  
        a2, b2 = avgFuncKspace(u1*u1/(Ng*Ng), self.qmod, int(self.Ng/4))
        
        LL = 2*np.pi*np.sum(b2)/np.sum(a2*b2)
        print (np.max(a2), np.max(b2), LL)
        return LL
            
    def configPlot(U, fig, n_, i):
        import matplotlib.pyplot as plt
        sp =  fig.add_subplot(2, 4, n_ )
        im=plt.pcolor(U, cmap=plt.cm.RdBu_r);  plt.clim(-1.1, 1.1); plt.axis('off'); plt.title('T = %1.2E'%(i))
        #cbar = plt.colorbar(im,fraction=0.04, pad=0.05, orientation="horizontal", ticks=[-1, 0, 1])
        plt.axis('off'); plt.title('T = %1.2E'%(i))     
  ```

```

 k_tilde, Ng   = 1, 64; 
Nt, dt, Nf = int(1e5), .02, 400
parameters = {'h':1, 'Ng':Ng, 'a':-0.25, 'b':0.25, 'kp':1, 'k_tilde':k_tilde,
         'Nt':Nt, 'dt':dt, 'Nf':Nf, 'Df':1}


am = FPS(parameters);    t1=time.perf_counter();     

phi0=0.0;  u0 = phi0 + (1-phi0)*(1-2*np.random.random((Ng,Ng)))          
am.integrate( fft2(u0) )

print ('total time taken: ', time.perf_counter()-t1)

fig = plt.figure(num=None, figsize=(18, 8), dpi=128);
ti=0;      configPlot(am.XX[ti,::].reshape(Ng, Ng), fig, 1, ti)
ti=10;     configPlot(am.XX[ti,::].reshape(Ng, Ng), fig, 2, ti)
ti=50;     configPlot(am.XX[ti,::].reshape(Ng, Ng), fig, 3, ti)
ti=80;     configPlot(am.XX[ti,::].reshape(Ng, Ng), fig, 4, ti)
ti=120;    configPlot(am.XX[ti,::].reshape(Ng, Ng), fig, 5, ti)
ti=200;    configPlot(am.XX[ti,::].reshape(Ng, Ng), fig, 6, ti)
ti=300;    configPlot(am.XX[ti,::].reshape(Ng, Ng), fig, 7, ti)
ti=399;    configPlot(am.XX[ti,::].reshape(Ng, Ng), fig, 8, ti);

X1 = am.XX

```

The value k_tilde can be made negative to see contractile stresses in action. The positive value of k_tilde gives rise to extensile stresses.
