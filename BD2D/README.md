<script type="text/javascript"
  src="http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML">
</script>

# Brownian motion

## Brownian motion of a free particle
The fluid surronding the particle is made of small molecules, which collide with the particle, giving randomly fluctuating forces. This effect can be described by the following model equation
$$m\frac{dv}{dt}=-\zeta v + F_r(t),$$
where \\(\zeta\\) is called the friction constant, and \\(F _ r(t)\\) stands for the fluctuating force on the particle by the fluid, with zero mean value and obeys the **fluctuation-dissipation relation**:
$$\langle F _ r(t _ 1) F_ r(t _ 2)\rangle = 2 \zeta k_B T \delta(t _ 1 - t _ 2).$$

The mean square displacement can be calculated as:
$$\langle (x(t)-x(0))^2\rangle = 2Dt\quad \mathrm{for}\quad t\gg \tau _ v,$$
where *D* is called the self-diffusion constant, standing for the fluctuation of particle position. *D* is related to the friction constant by the **Einstein relaltion**:
$$D=\frac{k _ B T}{\zeta},$$
which is another example of the fluctuation-dissipation relation.

Ref: M. Doi *Soft Matter Physics*
## Multiple particles with external forces and torques

### 2D
The evolution of Brownian particles is governed by the coupled overdamped Langevin equations,
$$\mathbf{\dot{r}}_i =\beta D \mathbf{F} _ i + \sqrt{2D} \mathbf{\eta}^T _i,$$
$$\dot{\theta}=\beta D _ r T _ i + \sqrt{2D _ r}\eta^R_i$$

Here, \\(\mathbf{F}\\) is an excluded-volume repulsive force, \\(\beta=1/k _ B T\\). *D* and *D*<sub>r</sub> are translational and rotational diffusion constants, which in the low-Reynolds-number regime are:
$$D=\frac{k _ B T}{6\pi \eta a},$$
$$D _ r =\frac{k _ B T}{8\pi \eta a^3},$$
so \\(D_r=3D/\sigma^2\\), where \\(\sigma=2 a\\) is the diamiter of the particle. The \\(\eta\\) are Gaussian white noise variables with \\(\langle \eta _ i(t)\rangle=0\\) and \\(\langle\eta _ i(t) \eta _ j(t')\rangle=\delta _ {ij} \delta (t-t')\\).

### 3D
To update the orientation of the particle, Lowen et al. use the following equation:
$$f _ r \partial _t \mathbf{\hat{u}} _ i = - \mathbf{T} _ i \times \mathbf{\hat{u}} _ i,$$
where \\(f _ r = 1/\beta D _ r\\) denotes the rotational friction constant, and \\(\mathbf{T} _ i = \mathbf{\hat{u}} _ i\times \nabla _ {\mathbf{\hat{u}} _ i} U\\) is the torque  on particle *i* and \\(U=1/2\sum _ {i, j\neq i}U _ {ij}\\) the total interaction potential. (Ref: A. Kaiser, K. Popowa and H. Lowen 2015 PRE *Active dipole clusters: From helical motion to fission*.) Besides, the orientation is affected by fluctuation term:
$$\sqrt{2D _ r} \mathbf{\hat{u}} _ i \times \mathbf{\eta} _ i,$$
see B. van der Meer et al. 2016 Soft matter *Removing grain boundaries from three-dimensional colloidal crystals using active dopants*. Here I infer that \\(\mathbf{\eta} _ i\\) is a unit vector with random orientation, which need be verified in the future.

Another way is to use quaternion, such as Zhao-Yan Sun et al. 2015 Soft matter *A versatile model for soft patchy particles with various patch arrangements*.

# Numerical Integration

A set of stochastic differential equations, can be written in the Ito interpretation as:
$$\dot{x} _ i = f _ i (x) + \sigma _ {ij}(x) \eta _ i (t).$$

## The Euler Method
$$x _ i(h) = x _ i(0) + f _ i (x(0))h + \sigma _ {ij}(x(0)) W _ j (h),$$
where \\(W _ i(h) = \int^h _ 0 ds\ \eta _ i (s)\\) can be replaced by a simple Caussian random variable \\(\bar{W} _ i (h)\\) with variance *h*. Moreover, we can use a uniform random number *R* on the interval (0, 1) to generate \\(\bar{W}_i\\) through
$$\bar{W} _ i(h)=\sqrt{12 h}(R-0.5).$$

Ref: A. Greiner, W. Strittmatter and J. Honerkamp 1988 J. Stat. Phys. *Numerical Integration of Stochastic Differential Equations*

