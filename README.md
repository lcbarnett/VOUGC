# VOUGC
A small [Matlab](https://www.mathworks.com/products/matlab.html)&trade; toolbox<sup>*</sup> for calculating [Granger causality](https://en.wikipedia.org/wiki/Granger_causality), conditional and unconditional, for vector [Ornstein-Uhlenbeck (VOU) processes](https://en.wikipedia.org/wiki/Ornstein%E2%80%93Uhlenbeck_process).

<sup>*NOTE: some functionality requires the [Matlab Control System Toolbox](https://www.mathworks.com/products/control.html), in particular the [icare](https://www.mathworks.com/help/control/ref/icare.html) function for solving continuous-time algebraic Riccati equations (CAREs).</sup>

A VOU process is a continuous-time multivariate linear autoregressive process expressed as a stochastic differential equation ([SDE](https://en.wikipedia.org/wiki/Stochastic_differential_equation)) of the form:

$$d\boldsymbol{y}(t) = A\\,\boldsymbol{y}(t)\\,dt + d\boldsymbol{w}(t)$$

Here $A$ is the autoregressive (AR) coefficients matrix and $\boldsymbol{w}(t)$ a [Wiener process](https://en.wikipedia.org/wiki/Wiener_process) with $d\boldsymbol{w}(t) \sim \mathcal{N}(0,\Sigma\\!dt)$, where $\Sigma$ is a covariance matrix. The parameters of the VOU model are $(A,\Sigma)$. The AR matrix $A$ need not be stable—i.e., it may have eigenvalues with nonnegative real part—but the covariance matrix $\Sigma$ must be positive-definite (so that the process is purely-nondeterministic).

This code implements computation of the _zero-horizon Granger causality rate_ for continuous-time stochastic processes introduced in ref. [1], using a state-space method developed in ref. [2]. A particular application of the computations facilitated by this toolbox, is to construct the _Granger causality maps_ introduced in ref. [3].

### Developer and maintainer
[Lionel Barnett](https://users.sussex.ac.uk/~lionelb/) ([lionelb@sussex.ac.uk](mailto:lionelb@sussex.ac.uk)), Department of Informatics, University of Sussex, UK.

### References
[1]: [L. Barnett and A. K. Seth, Detectability of Granger causality for subsampled continuous-time neurophysiological processes, _J. Neurosci. Methods_ **275**, 93-121, 2017](http://www.sciencedirect.com/science/article/pii/S0165027016302564).

[2]: [L. Barnett and A. K. Seth, Granger causality for state-space models, _Phys. Rev. E **91**, 040101(R)_, 2015](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.91.040101).

[3]: [B. Wahl _et al_., Granger-causality maps of diffusion processes, _Phys. Rev. E_ **93**, 022213, 2016](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.93.022213).
