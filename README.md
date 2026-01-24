# Diffusion Monte Carlo

This repository estimates the ground state of two well-known Bosonic systems: the Hydrogen ($`H`$), Hydrogen ion ($`H_2^+`$), and Hydrogen molecule ($`H_2`$) via Diffusion Monte Carlo (DMC).

![DMC simulation results for H, H2+, and H2](https://raw.githubusercontent.com/subinbg/Diff_MC/master/images/all.png)

**Figure 1.** Diffusion Monte Carlo simulations of $`H`$, $`H_2^+`$, and $`H_2`$. Simulations were executed with 10000 replicas and imaginary time interval of 0.01 atomic units.

---

## Table of Contents

- [Theoretical Background](#theoretical-background)
- [Integral Formulation](#integral-formulation-of-the-imaginary-time-schrödinger-equation)
- [Pure Diffusion Monte Carlo](#pure-diffusion-monte-carlo)
  - [Theory](#theory)
  - [Implementation](#implementation)
- [Importance Sampled Diffusion Monte Carlo](#importance-sampled-diffusion-monte-carlo)
  - [Motivation: Why Importance Sampling?](#motivation-why-importance-sampling)
  - [Theory](#theory-1)
    - [The Generalized Diffusion Equation](#the-generalized-diffusion-equation)
    - [The Local Energy](#the-local-energy)
    - [The Cusp Condition (Stability)](#the-cusp-condition-stability)
    - [The Drift-Diffusion Green's Function](#the-drift-diffusion-greens-function)
    - [Asymmetry of the Green's Function](#asymmetry-of-the-greens-function)
  - [Implementation](#implementation-1)

---

## Theoretical Background

### Time-Dependent Schrödinger Equation

We start with the time-dependent Schrödinger equation:

```math
i\frac{\partial}{\partial t} | \psi (\mathbf{x}, t) \rangle =\hat{H} | \psi (\mathbf{x}, t) \rangle
```

---

### Original Definition with the Auxiliary Constant Term $`E_T`$

We begin with the standard time-dependent Schrödinger equation for a single particle of mass $`m`$ in a potential $`\hat{V}`$. We introduce a constant energy shift $`E_T`$ (trial energy) into the Hamiltonian. The modified equation is:

```math
i\hbar \frac{\partial}{\partial t} | \Psi(t) \rangle = (\hat{H} - E_T) | \Psi(t) \rangle = (\hat{T} + \hat{V} - E_T) | \Psi(t) \rangle
```

> **Note on $`E_T`$:** Adding a constant scalar $`E_T`$ to the Hamiltonian shifts the energy spectrum ($`E_n \to E_n - E_T`$) but does not change the eigenstates. $`| \phi_n \rangle`$ remains an eigenstate of $`\hat{H} - E_T`$ if it is an eigenstate of $`\hat{H}`$. This constant is crucial for controlling the norm of the wavefunction in the imaginary time formalism.

---

### Explicit Expression of the Kinetic Operator

The kinetic energy operator $`\hat{T}`$ in the position representation is defined as:

```math
\hat{T} = -\frac{\hbar^2}{2m} \nabla_{\mathbf{x}}^2
```

where $`\nabla_{\mathbf{x}}^2`$ is the Laplacian operator acting on the spatial coordinates $`\mathbf{x}`$.

---

### Rewriting the Equation Using Imaginary Time

We perform a Wick rotation by substituting real time $`t`$ with imaginary time $`\tau`$ such that $`\tau = it`$, or $`t = -i\tau`$.

Using the chain rule $`\frac{\partial}{\partial t} = \frac{\partial \tau}{\partial t} \frac{\partial}{\partial \tau} = i \frac{\partial}{\partial \tau}`$, the equation transforms as follows:

```math
i\hbar \left( i \frac{\partial}{\partial \tau} \right) | \Psi(\tau) \rangle = (\hat{H} - E_T) | \Psi(\tau) \rangle
```

```math
-\hbar \frac{\partial}{\partial \tau} | \Psi(\tau) \rangle = (\hat{H} - E_T) | \Psi(\tau) \rangle
```

Rearranging for the time derivative gives the **Imaginary-Time Schrödinger Equation**:

```math
\frac{\partial}{\partial \tau} | \Psi(\tau) \rangle = -\frac{1}{\hbar} (\hat{H} - E_T) | \Psi(\tau) \rangle = -\frac{1}{\hbar} \left( -\frac{\hbar^2}{2m} \nabla_{\mathbf{x}}^2 + \hat{V} - E_T \right) | \Psi(\tau) \rangle
```

---

### Schrödinger Equation as the Diffusion Equation

#### Comparison to the Diffusion Equation

Projecting the imaginary-time equation onto the position basis $`| \mathbf{x} \rangle`$, where $`\Psi(\mathbf{x}, \tau) = \langle \mathbf{x} | \Psi(\tau) \rangle`$:

```math
\frac{\partial \Psi(\mathbf{x}, \tau)}{\partial \tau} = \frac{\hbar}{2m} \nabla_{\mathbf{x}}^2 \Psi(\mathbf{x}, \tau) - \frac{1}{\hbar}(V(\mathbf{x}) - E_T) \Psi(\mathbf{x}, \tau)
```

This equation is structurally isomorphic to the classical **Diffusion-Reaction Equation**:

```math
\frac{\partial \rho(\mathbf{x}, t)}{\partial t} = D \nabla_{\mathbf{x}}^2 \rho(\mathbf{x}, t) - k(\mathbf{x}) \rho(\mathbf{x}, t)
```

By comparing terms, we identify:

| **Quantity** | **Expression** |
|:-------------|:---------------|
| Diffusion coefficient | $`D = \frac{\hbar}{2m}`$ |
| Reaction rate (source/sink term) | $`k(\mathbf{x}) = \frac{1}{\hbar}(V(\mathbf{x}) - E_T)`$ |

---

### From Perron-Frobenius to Krein-Rutman Theorem

For the probabilistic interpretation of DMC to be valid, the wavefunction $`\Psi(\mathbf{x}, \tau)`$ must be interpreted as a probability density (or concentration), requiring it to be real and non-negative everywhere.

- **Perron-Frobenius Theorem:** Applies to finite-dimensional matrices with strictly positive entries. It states the eigenvector with the largest eigenvalue has strictly positive components.

- **Krein-Rutman Theorem:** Generalizes this to infinite-dimensional Banach spaces. The imaginary time propagator $`\hat{U}(\tau) = e^{-(\hat{H}-E_T)\tau/\hbar}`$ is a compact linear operator that leaves the cone of non-negative functions invariant.

Because the kinetic Green's function (Gaussian) is strictly positive and $`e^{-(V(\mathbf{x})-E_T)\tau/\hbar}`$ is strictly positive, the kernel of $`\hat{U}`$ is positive.

> **Conclusion:** The eigenstate with the maximal algebraic eigenvalue (which corresponds to the lowest physical energy ground state $`\Phi_0(\mathbf{x})`$) is unique, real, and strictly positive (nodeless).

---

### Mathematical vs. Quantum Mechanical Ground State

| **System** | **Description** |
|:-----------|:----------------|
| **Mathematical Ground State** | The Krein-Rutman theorem guarantees a nodeless ground state for any Hamiltonian of the form $`\hat{T} + \hat{V}`$. This is the Bosonic ground state. |
| **Bosonic Systems** | The physical ground state is symmetric and nodeless. DMC converges exactly to the physical ground state. |
| **Fermionic Systems** | The physical ground state must be antisymmetric ($`\Psi(\mathbf{x}_1, \mathbf{x}_2) = -\Psi(\mathbf{x}_2, \mathbf{x}_1)`$), implying the existence of nodal surfaces (where $`\Psi=0`$). |

The "Mathematical Ground State" (Bosonic) has lower energy than the Fermionic one. Without intervention, DMC collapses to the Bosonic state.

> **Fixed-Node Approximation:** We impose an artificial boundary condition $`\Psi(\mathbf{x}) = 0`$ at the nodes of a trial wavefunction $`\Psi_T`$. This effectively solves the Bosonic problem inside a restricted domain, providing an upper bound to the Fermionic energy.

---

### Proof of the Gaussian Green's Function (Free Diffusion)

Consider the case without potential ($`V(\mathbf{x}) = 0, E_T = 0`$). The equation is:

```math
\frac{\partial \Psi(\mathbf{x}, \tau)}{\partial \tau} = \frac{\hbar}{2m} \nabla_{\mathbf{x}}^2 \Psi(\mathbf{x}, \tau)
```

We solve via Fourier Transform. Let $`\tilde{\Psi}(\mathbf{k}, \tau) = \int d\mathbf{x} e^{-i\mathbf{k}\cdot\mathbf{x}} \Psi(\mathbf{x}, \tau)`$.

```math
\frac{\partial \tilde{\Psi}(\mathbf{k}, \tau)}{\partial \tau} = -\frac{\hbar k^2}{2m} \tilde{\Psi}(\mathbf{k}, \tau)
```

**Solution in k-space:**

```math
\tilde{\Psi}(\mathbf{k}, \tau) = \tilde{\Psi}(\mathbf{k}, 0) \exp\left( -\frac{\hbar k^2}{2m}\tau \right)
```

**Inverse Fourier transform** using the Convolution Theorem:

```math
\Psi(\mathbf{x}, \tau) = \int d\mathbf{y} G_0(\mathbf{x}-\mathbf{y}, \tau) \Psi(\mathbf{y}, 0)
```

Where the Green's function $`G_0`$ is the inverse transform of the Gaussian decay:

```math
G_0(\mathbf{x}, \tau) = \frac{1}{(2\pi)^3} \int d\mathbf{k} e^{i\mathbf{k}\cdot\mathbf{x}} e^{-\frac{\hbar k^2}{2m}\tau} = \left( \frac{m}{2\pi\hbar\tau} \right)^{3/2} \exp\left( -\frac{m |\mathbf{x}|^2}{2\hbar\tau} \right)
```

This is a normalized Gaussian distribution with variance $`\sigma^2 = \frac{\hbar\tau}{m}`$.

---

### Proof of Invariance of Total Number of Particles

Let $`N(\tau) = \int d\mathbf{x} \Psi(\mathbf{x}, \tau)`$.

```math
\frac{dN}{d\tau} = \int d\mathbf{x} \frac{\partial \Psi}{\partial \tau} = \frac{\hbar}{2m} \int d\mathbf{x} \nabla_{\mathbf{x}}^2 \Psi(\mathbf{x}, \tau) - \frac{1}{\hbar} \int d\mathbf{x} (V(\mathbf{x}) - E_T) \Psi(\mathbf{x}, \tau)
```

Using the Divergence Theorem on the kinetic term (assuming $`\Psi`$ vanishes at infinity):

```math
\int_{\text{vol}} \nabla \cdot (\nabla \Psi) d\mathbf{x} = \oint_{\text{surf}} (\nabla \Psi) \cdot d\mathbf{S} = 0
```

Thus, for pure diffusion ($`V=0`$), $`\frac{dN}{d\tau} = 0`$. The population is conserved.

With potential, the population change rate is:

```math
\frac{dN}{d\tau} = - \frac{1}{\hbar} \int d\mathbf{x} (V(\mathbf{x}) - E_T) \Psi(\mathbf{x}, \tau)
```

---

## Integral Formulation of the Imaginary-Time Schrödinger Equation

### Fourier Transforms and Momentum States

We define the momentum eigenstate $`| \mathbf{k} \rangle`$ such that $`\hat{\mathbf{p}} | \mathbf{k} \rangle = \hbar \mathbf{k} | \mathbf{k} \rangle`$.

The overlap with position states is:

```math
\langle \mathbf{x} | \mathbf{k} \rangle = \frac{1}{(2\pi)^{3/2}} e^{i\mathbf{k}\cdot\mathbf{x}}
```

Identity operator in momentum space:

```math
\hat{I} = \int d\mathbf{k} | \mathbf{k} \rangle \langle \mathbf{k} |
```

---

### Trotter-Suzuki Decomposition

The formal solution to $`\frac{\partial}{\partial \tau} | \Psi(\tau) \rangle = -\frac{1}{\hbar}(\hat{H}-E_T) | \Psi(\tau) \rangle`$ over a small time step $`\delta\tau`$ is:

```math
| \Psi(\tau + \delta\tau) \rangle = e^{-\frac{\delta\tau}{\hbar}(\hat{T} + \hat{V} - E_T)} | \Psi(\tau) \rangle
```

Since $`[\hat{T}, \hat{V}] \neq 0`$, we use the Trotter-Suzuki decomposition:

```math
e^{-\frac{\delta\tau}{\hbar}(\hat{T} + \hat{V} - E_T)} = e^{-\frac{\delta\tau}{\hbar}(\hat{V} - E_T)} e^{-\frac{\delta\tau}{\hbar}\hat{T}} + O(\delta\tau^2)
```

This splits the evolution into a kinetic step followed by a potential step.

---

### Expansion in Position Space

We project the evolution onto position space:

```math
\Psi(\mathbf{x}, \tau+\delta\tau) = \langle \mathbf{x} | e^{-\frac{\delta\tau}{\hbar}(\hat{V} - E_T)} e^{-\frac{\delta\tau}{\hbar}\hat{T}} | \Psi(\tau) \rangle + O(\delta\tau^2)
```

Since $`\hat{V}`$ is diagonal in $`| \mathbf{x} \rangle`$:

```math
\Psi(\mathbf{x}, \tau+\delta\tau) \approx e^{-\frac{\delta\tau}{\hbar}(V(\mathbf{x}) - E_T)} \langle \mathbf{x} | e^{-\frac{\delta\tau}{\hbar}\hat{T}} | \Psi(\tau) \rangle
```

Insert identity $`\int d\mathbf{y} | \mathbf{y} \rangle \langle \mathbf{y} | = \hat{I}`$:

```math
\Psi(\mathbf{x}, \tau+\delta\tau) \approx e^{-\frac{\delta\tau}{\hbar}(V(\mathbf{x}) - E_T)} \int d\mathbf{y} \langle \mathbf{x} | e^{-\frac{\delta\tau}{\hbar}\hat{T}} | \mathbf{y} \rangle \Psi(\mathbf{y}, \tau)
```

---

### Proof of Kinetic Operator as Green's Function

We evaluate the matrix element $`K(\mathbf{x}, \mathbf{y}) = \langle \mathbf{x} | e^{-\frac{\delta\tau}{\hbar}\hat{T}} | \mathbf{y} \rangle`$ by inserting momentum states $`\int d\mathbf{k} | \mathbf{k} \rangle \langle \mathbf{k} |`$:

```math
K(\mathbf{x}, \mathbf{y}) = \int d\mathbf{k} \langle \mathbf{x} | e^{-\frac{\delta\tau}{\hbar} \frac{\hat{p}^2}{2m}} | \mathbf{k} \rangle \langle \mathbf{k} | \mathbf{y} \rangle
```

```math
K(\mathbf{x}, \mathbf{y}) = \int d\mathbf{k} e^{-\frac{\delta\tau}{\hbar} \frac{\hbar^2 k^2}{2m}} \langle \mathbf{x} | \mathbf{k} \rangle \langle \mathbf{k} | \mathbf{y} \rangle
```

```math
K(\mathbf{x}, \mathbf{y}) = \frac{1}{(2\pi)^3} \int d\mathbf{k} e^{-\frac{\hbar \delta\tau k^2}{2m}} e^{i\mathbf{k}\cdot(\mathbf{x}-\mathbf{y})}
```

This is the standard Gaussian integral result:

```math
K(\mathbf{x}, \mathbf{y}) = \left( \frac{m}{2\pi\hbar\delta\tau} \right)^{3/2} \exp\left( -\frac{m |\mathbf{x}-\mathbf{y}|^2}{2\hbar\delta\tau} \right) \equiv G_{\text{diff}}(\mathbf{x}-\mathbf{y}, \delta\tau)
```

Thus, the integral update equation is:

```math
\Psi(\mathbf{x}, \tau+\delta\tau) = \int d\mathbf{y} \underbrace{e^{-\frac{\delta\tau}{\hbar}(V(\mathbf{x}) - E_T)}}_{\text{Branching Weight } W} \underbrace{G_{\text{diff}}(\mathbf{x}-\mathbf{y}, \delta\tau)}_{\text{Diffusion Probability}} \Psi(\mathbf{y}, \tau) + O(\delta\tau^2)
```

---

### Stationarity and Eigenfunction Correspondence

Assume stationarity, i.e., $`\Psi(\mathbf{x}, \tau+\delta\tau) = \Psi(\mathbf{x}, \tau) \equiv \Phi(\mathbf{x})`$.

```math
\Phi(\mathbf{x}) \approx \int d\mathbf{y} e^{-\frac{\delta\tau}{\hbar}(V(\mathbf{x}) - E_T)} G_{\text{diff}}(\mathbf{x}-\mathbf{y}, \delta\tau) \Phi(\mathbf{y})
```

Expanding the exponentials to first order in $`\delta\tau`$:

```math
\Phi(\mathbf{x}) \approx \left( 1 - \frac{\delta\tau}{\hbar}(V(\mathbf{x}) - E_T) \right) \int d\mathbf{y} G_{\text{diff}}(\mathbf{x}-\mathbf{y}, \delta\tau) \Phi(\mathbf{y})
```

The integral term represents diffusion over time $`\delta\tau`$: $`\Phi(\mathbf{y}) + \frac{\hbar \delta\tau}{2m} \nabla_{\mathbf{x}}^2 \Phi(\mathbf{x})`$.

```math
\Phi(\mathbf{x}) \approx \left( 1 - \frac{\delta\tau}{\hbar}(V(\mathbf{x}) - E_T) \right) \left( \Phi(\mathbf{x}) + \frac{\hbar \delta\tau}{2m} \nabla_{\mathbf{x}}^2 \Phi(\mathbf{x}) \right)
```

Dropping $`O(\delta\tau^2)`$ terms and canceling $`\Phi(\mathbf{x})`$:

```math
0 = \frac{\hbar \delta\tau}{2m} \nabla_{\mathbf{x}}^2 \Phi(\mathbf{x}) - \frac{\delta\tau}{\hbar}(V(\mathbf{x}) - E_T) \Phi(\mathbf{x})
```

```math
-\frac{\hbar^2}{2m} \nabla_{\mathbf{x}}^2 \Phi(\mathbf{x}) + V(\mathbf{x}) \Phi(\mathbf{x}) = E_T \Phi(\mathbf{x})
```

This recovers $`\hat{H} \Phi(\mathbf{x}) = E_T \Phi(\mathbf{x})`$. The stationary density is an eigenstate, and the trial energy $`E_T`$ is the eigenvalue.

---

### Convergence to the Ground State

The general solution in terms of eigenstates $`\{ | \phi_n \rangle \}`$ with energies $`E_n`$:

```math
| \Psi(\tau) \rangle = \sum_n c_n e^{-(E_n - E_T)\tau/\hbar} | \phi_n \rangle
```

As $`\tau \to \infty`$, the term with the minimum $`E_n`$ (ground state $`E_0`$) dominates, provided $`c_0 \neq 0`$:

```math
| \Psi(\tau \to \infty) \rangle \approx c_0 e^{-(E_0 - E_T)\tau/\hbar} | \phi_0 \rangle
```

This ensures convergence to the Bosonic ground state.

---

### Time Derivative of Number of Particles

If the system has converged to the ground state $`\Phi_0(\mathbf{x})`$, the population $`N(\tau)`$ behaves as:

```math
\frac{dN}{d\tau} = \frac{d}{d\tau} \langle \Phi_0 | \Psi(\tau) \rangle \approx \frac{d}{d\tau} \left( e^{-(E_0 - E_T)\tau/\hbar} \right) = -\frac{1}{\hbar}(E_0 - E_T) N(\tau)
```

---

### Population Control Bias and Feedback Law

To stabilize the population, we employ a feedback law for $`E_T(\tau)`$.

Let the estimator be the average potential energy $`\bar{V}(\tau)`$.

**Control Law:**

```math
E_T(\tau) = \bar{V}(\tau) - \alpha \ln\left( \frac{N(\tau)}{N_0} \right)
```

where $`N_0`$ is the target population.

#### Local Convergence Proof

Near equilibrium, let $`N(\tau) = N_0(1+\delta(\tau))`$ with $`\delta \ll 1`$. Assuming $`\bar{V} \approx E_0`$:

```math
\frac{dN}{d\tau} \approx -\frac{1}{\hbar}(E_0 - E_T) N_0
```

Substituting the control law (using $`\ln(1+\delta) \approx \delta`$):

```math
E_0 - E_T \approx E_0 - (E_0 - \alpha \delta) = \alpha \delta
```

```math
N_0 \frac{d\delta}{d\tau} \approx -\frac{1}{\hbar} (\alpha \delta) N_0 \implies \frac{d\delta}{d\tau} = -\frac{\alpha}{\hbar} \delta
```

This yields exponential decay $`\delta(\tau) \propto e^{-\frac{\alpha}{\hbar}\tau}`$, proving local stability of the population.

#### Bias

The estimator for Energy is $`E_{\text{DMC}} = \langle \frac{H\Psi}{\Psi} \rangle`$. Since $`E_T`$ depends on $`N`$, $`E_{\text{DMC}}`$ includes a bias term proportional to $`1/N`$ due to the correlation between the fluctuating population size and the energy estimate:

```math
E_{\text{DMC}} = E_{\text{exact}} + O\left( \frac{1}{N} \right)
```

---

## Pure Diffusion Monte Carlo

### Theory

Pure DMC directly simulates the imaginary-time Schrödinger equation by interpreting the wavefunction $`\Psi(\mathbf{x}, \tau)`$ as a probability density represented by a population of random walkers.

The update equation derived from Trotter-Suzuki decomposition is:

```math
\Psi(\mathbf{x}, \tau+\delta\tau) = \int d\mathbf{y} \, W(\mathbf{x}) \, G_{\text{diff}}(\mathbf{x}-\mathbf{y}, \delta\tau) \, \Psi(\mathbf{y}, \tau)
```

where:

| **Component** | **Expression** | **Role** |
|:--------------|:---------------|:---------|
| Diffusion kernel | $`G_{\text{diff}}(\mathbf{x}-\mathbf{y}) = \left(\frac{m}{2\pi\hbar\delta\tau}\right)^{3/2} e^{-\frac{m\|\mathbf{x}-\mathbf{y}\|^2}{2\hbar\delta\tau}}`$ | Gaussian random walk |
| Branching weight | $`W(\mathbf{x}) = e^{-\frac{\delta\tau}{\hbar}(V(\mathbf{x}) - E_T)}`$ | Birth/death of walkers |

**Key Properties:**
- The Green's function is **symmetric**: $`G_{\text{diff}}(\mathbf{x}-\mathbf{y}) = G_{\text{diff}}(\mathbf{y}-\mathbf{x})`$
- Walkers diffuse isotropically (no preferred direction)
- The weight $`W`$ depends directly on the potential $`V(\mathbf{x})`$

---

### Implementation

This section describes the algorithm used in the Fortran code provided in this repository. It solves for the distribution of walkers $`\rho(\mathbf{x}) \propto \Psi(\mathbf{x}, \tau)`$.

**Objective:** Evolve a population of $`N`$ walkers $`\{ \mathbf{x}_1, \dots, \mathbf{x}_N \}`$ in imaginary time to sample the ground state wavefunction $`\Phi_0(\mathbf{x})`$.

---

#### Step 1: Initialization

**Math:** Set $`\tau = 0`$. Initialize $`N`$ walkers at positions $`\mathbf{x}_i`$ according to an arbitrary initial probability distribution $`|\Psi(\mathbf{x}, 0)|^2`$ (often a delta function at the origin or a simple Gaussian). Initialize the Trial Energy $`E_T(\tau=0)`$ to a guess (e.g., the potential energy at the start).

**Implementation:**

```fortran
x(1:N) = 0.0          ! All walkers at origin
life(1:N) = 1         ! 1 = Alive, -1 = Dead
potential_avg = ...   ! Initial guess for E_T
```

---

#### Step 2: Diffusion (Kinetic Update)

**Equation:** The kinetic propagator $`e^{-\frac{\delta\tau}{\hbar}\hat{T}}`$ corresponds to sampling from the Gaussian Green's function:

```math
\mathbf{x}_i(\tau + \delta\tau) = \mathbf{x}_i(\tau) + \boldsymbol{\xi}
```

where $`\boldsymbol{\xi}`$ is a random vector distributed as:

```math
G_{\text{diff}}(\boldsymbol{\xi}) = \left( \frac{m}{2\pi\hbar\delta\tau} \right)^{3/2} \exp\left( -\frac{m |\boldsymbol{\xi}|^2}{2\hbar\delta\tau} \right)
```

**Implementation:**

```fortran
! Box-Muller generates standard normal random number (mean=0, var=1)
call box_muller(ranseed)
! Scale by sqrt(hbar * dt / m). In atomic units (hbar=1, m=1), scale = sqrt(dt)
x(i,j,k) = x(i,j,k) + sqrt(dt) * ranseed
```

---

#### Step 3: Weighting and Branching (Potential Update)

**Equation:** The potential propagator $`e^{-\frac{\delta\tau}{\hbar}(\hat{V} - E_T)}`$ acts as a multiplicative weight $`W`$.

```math
W_i = \exp\left( -\frac{\delta\tau}{\hbar} (V(\mathbf{x}_i) - E_T) \right)
```

> **Note:** The code uses the average potential as $`E_T`$, which is `potential_avg`.

**Implementation (Stochastic Branching):**

Instead of carrying decimal weights, we convert $`W_i`$ into an integer number of surviving walkers $`M_i`$.

```math
M_i = \lfloor W_i + u \rfloor
```

where $`u`$ is a uniform random number $`\in [0, 1)`$. This technique preserves the average weight $`\mathbb{E}[M_i] = W_i`$.

```fortran
W = exp(-(V(x(i)) - potential_avg) * dt)
replica_test = min(int(W + ranseed), 3)  ! Cap offspring at 3 to prevent explosion

if (replica_test == 0) then
   life(i) = -1              ! Walker dies
else if (replica_test == 2) then
   ! Clone: Find a dead slot and copy current walker x(i) into it
else if (replica_test == 3) then
   ! Clone twice
end if
```

---

#### Step 4: $`E_T`$ Computation (Population Control)

**Equation:** Update $`E_T`$ to stabilize the population $`N(\tau)`$ near the target $`N_0`$.

```math
E_T(\tau) = \langle V \rangle_{\text{walkers}} - \alpha \ln \left( \frac{N(\tau)}{N_0} \right)
```

**Implementation:**

The Fortran code uses a simplified estimator where $`E_T`$ is just the average potential of the current population.

```fortran
potential_avg = 0.0
do i = 1, N
   potential_avg = potential_avg + V(x(i))
end do
potential_avg = potential_avg / N_alive
! The code comments out the population correction term (alpha * ln(N/N0))
! implying it relies on the correlation between V and N for stability (meta-stable).
```

---

#### Step 5: Energy Estimation

**Equation:** The ground state energy $`E_0`$ is the average of the local potential (plus the kinetic contribution, which averages to the virial expectation). In pure DMC, we simply average the Trial Energy $`E_T`$ over the simulation steps once converged.

```math
E_0 \approx \frac{1}{M} \sum_{k=1}^M E_T(\tau_k)
```

---

## Importance Sampled Diffusion Monte Carlo

### Motivation: Why Importance Sampling?

Pure DMC, while mathematically elegant, suffers from a fundamental numerical instability when applied to atomic and molecular systems with Coulomb potentials.

**The Problem: Coulomb Singularity**

For a hydrogen-like atom, the potential is:

```math
V(r) = -\frac{Z}{r} \quad \xrightarrow{r \to 0} \quad -\infty
```

When a walker approaches the nucleus ($`r \to 0`$), the branching weight becomes:

```math
W = \exp\left( -\frac{\delta\tau}{\hbar}(V(r) - E_T) \right) = \exp\left( +\frac{Z\delta\tau}{\hbar r} \right) \xrightarrow{r \to 0} +\infty
```

This causes **numerical explosion**: walkers near the nucleus are replicated exponentially, leading to unbounded population growth, divergent energy estimates, and complete simulation breakdown.

**The Solution: Importance Sampling**

The key insight is to **replace the divergent potential $`V(\mathbf{x})`$ with a bounded "Local Energy" $`E_L(\mathbf{x})`$** by introducing a trial wavefunction $`\Psi_T(\mathbf{x})`$ that captures the correct behavior near singularities.

If $`\Psi_T`$ satisfies the **Cusp Condition** (discussed below), the kinetic energy singularity from $`\nabla^2 \Psi_T`$ exactly cancels the potential energy singularity from $`V(\mathbf{x})`$, leaving $`E_L(\mathbf{x})`$ finite everywhere.

> **Bottom Line:** Importance sampling transforms an inherently unstable algorithm into a robust, practical method for quantum Monte Carlo calculations.

---

### Theory

Importance sampling introduces a Trial Wavefunction $`\Psi_T(\mathbf{x})`$ to guide walkers into physically relevant regions. This transforms the equation for $`\Psi`$ into an equation for a new distribution $`f`$.

---

#### The Generalized Diffusion Equation

##### 1. Definition of the Target Distribution $`f`$

We define the distribution $`f(\mathbf{x}, \tau)`$ as the product of the true wavefunction and the trial wavefunction:

```math
f(\mathbf{x}, \tau) = \Psi(\mathbf{x}, \tau) \Psi_T(\mathbf{x})
```

> **Note:** If $`\Psi(\mathbf{x}, \tau) \to \Phi_0(\mathbf{x})`$ (the ground state), then $`f \to \Phi_0(\mathbf{x}) \Psi_T(\mathbf{x})`$. If $`\Psi_T \approx \Phi_0`$, then $`f \approx |\Phi_0|^2`$, which is the physical probability density.

---

##### 2. Time Derivative of $`f`$

Differentiate $`f`$ with respect to $`\tau`$:

```math
\frac{\partial f}{\partial \tau} = \frac{\partial \Psi}{\partial \tau} \Psi_T = -\frac{1}{\hbar} (\hat{H} - E_T) \Psi \cdot \Psi_T
```

Substitute $`\hat{H} = -\frac{\hbar^2}{2m} \nabla^2 + V`$:

```math
\frac{\partial f}{\partial \tau} = \Psi_T \left[ \frac{\hbar}{2m} \nabla^2 \Psi - \frac{1}{\hbar}(V - E_T)\Psi \right]
```

---

##### 3. Kinetic Term Transformation

We need to express $`\nabla^2 \Psi`$ in terms of $`f`$.

Since $`\Psi = f / \Psi_T`$, we apply the Laplacian:

```math
\nabla^2 \left( \frac{f}{\Psi_T} \right) = \frac{\nabla^2 f}{\Psi_T} - 2 \frac{\nabla f \cdot \nabla \Psi_T}{\Psi_T^2} - f \frac{\nabla^2 \Psi_T}{\Psi_T^2} + 2 f \frac{(\nabla \Psi_T)^2}{\Psi_T^3}
```

This algebra is cumbersome. Alternatively, use the identity:

```math
\Psi_T \nabla^2 \left( \frac{f}{\Psi_T} \right) = \nabla^2 f - \nabla \cdot \left( 2 f \frac{\nabla \Psi_T}{\Psi_T} \right) - f \frac{\nabla^2 \Psi_T}{\Psi_T}
```

Substitute this back into the time derivative equation:

```math
\frac{\partial f}{\partial \tau} = \frac{\hbar}{2m} \left[ \nabla^2 f - \nabla \cdot \left( 2 f \frac{\nabla \Psi_T}{\Psi_T} \right) - f \frac{\nabla^2 \Psi_T}{\Psi_T} \right] - \frac{1}{\hbar}(V - E_T)f
```

---

##### 4. Grouping Terms (The Local Energy)

We identify the **Quantum Drift Velocity** $`\mathbf{v}_D(\mathbf{x})`$:

```math
\mathbf{v}_D(\mathbf{x}) \equiv \frac{\hbar}{m} \frac{\nabla \Psi_T(\mathbf{x})}{\Psi_T(\mathbf{x})} = \frac{\hbar}{m} \nabla \ln |\Psi_T(\mathbf{x})|
```

And we group the potential $`V`$ and the kinetic remnant $`\frac{\nabla^2 \Psi_T}{\Psi_T}`$ into the **Local Energy** $`E_L(\mathbf{x})`$:

```math
E_L(\mathbf{x}) \equiv -\frac{\hbar^2}{2m} \frac{\nabla^2 \Psi_T(\mathbf{x})}{\Psi_T(\mathbf{x})} + V(\mathbf{x}) = \frac{\hat{H} \Psi_T(\mathbf{x})}{\Psi_T(\mathbf{x})}
```

---

##### 5. The Generalized Diffusion Equation

Substituting these definitions yields the final form:

```math
\frac{\partial f(\mathbf{x}, \tau)}{\partial \tau} = \underbrace{\frac{\hbar}{2m} \nabla^2 f(\mathbf{x}, \tau)}_{\text{Diffusion}} - \underbrace{\nabla \cdot (\mathbf{v}_D(\mathbf{x}) f(\mathbf{x}, \tau))}_{\text{Drift}} - \underbrace{\frac{1}{\hbar}(E_L(\mathbf{x}) - E_T) f(\mathbf{x}, \tau)}_{\text{Branching / Source-Sink}}
```

---

#### The Drift-Diffusion Green's Function

##### The Fokker-Planck Equation

The generalized diffusion equation (without the branching term) is a **Fokker-Planck equation**:

```math
\frac{\partial f}{\partial \tau} = D \nabla^2 f - \nabla \cdot (\mathbf{v}_D f)
```

where $`D = \frac{\hbar}{2m}`$ is the diffusion coefficient. This describes the probability density $`f(\mathbf{x}, \tau)`$ of particles undergoing both diffusion and drift.

---

##### Solving via the Langevin Equation

The Fokker-Planck equation is equivalent to the **Langevin stochastic differential equation** (SDE):

```math
d\mathbf{x} = \mathbf{v}_D(\mathbf{x}) d\tau + \sqrt{2D} \, d\mathbf{W}
```

where $`d\mathbf{W}`$ is a Wiener process (Brownian motion) with $`\langle d\mathbf{W} \rangle = 0`$ and $`\langle dW_i dW_j \rangle = \delta_{ij} d\tau`$.

**Discretizing** this SDE over a small time step $`\delta\tau`$:

```math
\mathbf{x}(\tau + \delta\tau) = \mathbf{x}(\tau) + \mathbf{v}_D(\mathbf{x}(\tau)) \delta\tau + \boldsymbol{\xi}
```

where $`\boldsymbol{\xi}`$ is a Gaussian random vector with:
- Mean: $`\langle \boldsymbol{\xi} \rangle = 0`$
- Variance: $`\langle \xi_i \xi_j \rangle = 2D \delta\tau \, \delta_{ij} = \frac{\hbar \delta\tau}{m} \delta_{ij}`$

---

##### The Green's Function (Transition Probability)

The probability of transitioning from $`\mathbf{y}`$ to $`\mathbf{x}`$ in time $`\delta\tau`$ is given by the Green's function. Since the displacement $`\boldsymbol{\xi}`$ is Gaussian, we have:

```math
G(\mathbf{x} \leftarrow \mathbf{y}, \delta\tau) = \frac{1}{(4\pi D \delta\tau)^{3/2}} \exp\left( -\frac{|\mathbf{x} - \mathbf{y} - \mathbf{v}_D(\mathbf{y})\delta\tau|^2}{4D\delta\tau} \right)
```

Substituting $`D = \frac{\hbar}{2m}`$ and $`\sigma^2 = 2D\delta\tau = \frac{\hbar\delta\tau}{m}`$:

```math
G(\mathbf{x} \leftarrow \mathbf{y}, \delta\tau) = \frac{1}{(2\pi \sigma^2)^{3/2}} \exp\left( -\frac{|\mathbf{x} - \mathbf{y} - \mathbf{v}_D(\mathbf{y})\delta\tau|^2}{2\sigma^2} \right)
```

This is a **Gaussian centered at** $`\mathbf{y} + \mathbf{v}_D(\mathbf{y})\delta\tau`$, not at $`\mathbf{y}`$.

---

#### Asymmetry of the Green's Function

##### Why Use Directional Notation $`G(\mathbf{x} \leftarrow \mathbf{y})`$?

In pure DMC, the Green's function is **symmetric**:

```math
G_{\text{pure}}(\mathbf{x}, \mathbf{y}) = G_{\text{pure}}(|\mathbf{x} - \mathbf{y}|) \implies G(\mathbf{x} \leftarrow \mathbf{y}) = G(\mathbf{y} \leftarrow \mathbf{x})
```

This is because the pure diffusion kernel is an isotropic Gaussian depending only on the distance $`|\mathbf{x} - \mathbf{y}|`$.

However, with drift, the Green's function becomes **asymmetric**:

```math
G(\mathbf{x} \leftarrow \mathbf{y}) \neq G(\mathbf{y} \leftarrow \mathbf{x})
```

**Physical Interpretation:** The drift velocity $`\mathbf{v}_D`$ depends on the *starting* position. A particle at $`\mathbf{y}`$ is pushed toward $`\mathbf{y} + \mathbf{v}_D(\mathbf{y})\delta\tau`$, while a particle at $`\mathbf{x}`$ is pushed toward $`\mathbf{x} + \mathbf{v}_D(\mathbf{x})\delta\tau`$. These are generally different directions.

**Example:** Consider a 1D hydrogen atom with $`\Psi_T = e^{-|x|}`$:
- At $`y = +1`$: $`v_D(y) = -\text{sign}(y) = -1`$ (drift toward nucleus)
- At $`x = -1`$: $`v_D(x) = -\text{sign}(x) = +1`$ (drift toward nucleus)

The transition $`y \to x`$ has drift pushing *left*, while $`x \to y`$ has drift pushing *right*. Thus $`G(x \leftarrow y) \neq G(y \leftarrow x)`$.

> **Notation Convention:** We use $`G(\mathbf{x} \leftarrow \mathbf{y})`$ to denote the transition probability *from* $`\mathbf{y}`$ *to* $`\mathbf{x}`$. The arrow indicates the direction of the transition, with the drift evaluated at the *source* point $`\mathbf{y}`$.

---

##### Consequence: Need for Detailed Balance Correction

The asymmetry of the Green's function means that naive sampling does **not** satisfy detailed balance:

```math
f(\mathbf{y}) G(\mathbf{x} \leftarrow \mathbf{y}) \neq f(\mathbf{x}) G(\mathbf{y} \leftarrow \mathbf{x})
```

This is why we need the **Metropolis acceptance/rejection step** to restore detailed balance (see [Step 3: Metropolis Acceptance/Rejection](#step-3-metropolis-acceptancerejection) in the Implementation section).

---

#### Summary: Pure DMC vs. Importance Sampled DMC

| **Feature** | **Pure DMC ($`\Psi`$)** | **Importance Sampled DMC ($`f = \Psi \Psi_T`$)** |
|:------------|:------------------------|:-------------------------------------------------|
| **Quantity Simulated** | Wavefunction $`\Psi(\mathbf{x})`$ | Mixed distribution $`f(\mathbf{x}) = \Psi \Psi_T`$ |
| **Propagator** | Isotropic Gaussian | Gaussian shifted by Drift |
| **Green's Function** | Symmetric: $`G(\mathbf{x}-\mathbf{y})`$ | **Asymmetric:** $`G(\mathbf{x} \leftarrow \mathbf{y}) \neq G(\mathbf{y} \leftarrow \mathbf{x})`$ |
| **Update Step** | $`\mathbf{x}' = \mathbf{x} + \boldsymbol{\xi}`$ | $`\mathbf{x}' = \mathbf{x} + \mathbf{v}_D(\mathbf{x})\delta\tau + \boldsymbol{\xi}`$ |
| **Weighting Term** | Potential Energy $`V(\mathbf{x})`$ | Local Energy $`E_L(\mathbf{x}) = \frac{\hat{H}\Psi_T}{\Psi_T}`$ |
| **Branching Factor** | $`W = \exp(-\frac{\delta\tau}{\hbar}(V - E_T))`$ | $`W = \exp(-\frac{\delta\tau}{\hbar}(E_L - E_T))`$ |
| **Detailed Balance** | Automatically satisfied | Requires Metropolis correction |
| **Singularity** | Unstable if $`V(\mathbf{x}) \to -\infty`$ | Stable if $`E_L(\mathbf{x})`$ is smooth (Cusp condition) |

---

#### Implementation Detail: Code Transformation

In the code, the pure DMC update:

```fortran
x(i,j,k) = x(i,j,k) + sqrt(dt) * ranseed
```

becomes for importance-sampled DMC:

```fortran
! Calculate drift based on Trial Wavefunction Gradient
drift = (hbar/m) * grad_ln_psi_T(x(i,j,:))
x(i,j,k) = x(i,j,k) + drift * dt + sqrt(dt) * ranseed
```

---

#### The Local Energy

The Local Energy is defined as the action of the Hamiltonian on the trial wavefunction, divided by the trial wavefunction itself:

```math
E_L(\mathbf{x}) \equiv \frac{\hat{H} \Psi_T(\mathbf{x})}{\Psi_T(\mathbf{x})} = -\frac{\hbar^2}{2m} \frac{\nabla^2 \Psi_T(\mathbf{x})}{\Psi_T(\mathbf{x})} + V(\mathbf{x})
```

---

##### Detailed Derivation

We derive this term naturally by transforming the imaginary-time Schrödinger equation for $`\Psi(\mathbf{x}, \tau)`$ into an equation for the product distribution $`f(\mathbf{x}, \tau) = \Psi(\mathbf{x}, \tau) \Psi_T(\mathbf{x})`$.

**Step 1: Start with the Schrödinger Equation**

```math
-\hbar \frac{\partial \Psi}{\partial \tau} = \hat{H} \Psi - E_T \Psi = -\frac{\hbar^2}{2m} \nabla^2 \Psi + (V - E_T)\Psi
```

**Step 2: Differentiate the Product Distribution $`f`$**

We take the time derivative of $`f = \Psi \Psi_T`$. Since $`\Psi_T`$ is time-independent:

```math
\frac{\partial f}{\partial \tau} = \Psi_T \frac{\partial \Psi}{\partial \tau} = -\frac{1}{\hbar} \Psi_T (\hat{H} - E_T) \Psi
```

**Step 3: Express $`\nabla^2 \Psi`$ in terms of $`f`$**

We need to replace $`\Psi`$ with $`f/\Psi_T`$ in the Hamiltonian.

Using the quotient rule for the Laplacian $`\nabla^2 (\frac{u}{v}) = \frac{\nabla^2 u}{v} - 2\frac{\nabla u \cdot \nabla v}{v^2} - u \frac{\nabla^2 v}{v^2} + 2 u \frac{(\nabla v)^2}{v^3}`$:

```math
\nabla^2 \Psi = \nabla^2 \left( \frac{f}{\Psi_T} \right) = \frac{\nabla^2 f}{\Psi_T} - 2\frac{\nabla f \cdot \nabla \Psi_T}{\Psi_T^2} + f \left( 2\frac{(\nabla \Psi_T)^2}{\Psi_T^3} - \frac{\nabla^2 \Psi_T}{\Psi_T^2} \right)
```

> **Note:** The term in the parenthesis can be rewritten using the identity $`\nabla^2 \ln \Psi_T = \frac{\nabla^2 \Psi_T}{\Psi_T} - \frac{(\nabla \Psi_T)^2}{\Psi_T^2}`$, but we stick to the expanded form for clarity.

**Step 4: Substitute and Group Terms**

Substitute $`\nabla^2 \Psi`$ back into the time derivative equation:

```math
\frac{\partial f}{\partial \tau} = -\frac{1}{\hbar} \Psi_T \left[ -\frac{\hbar^2}{2m} \left( \frac{\nabla^2 f}{\Psi_T} - 2\frac{\nabla f \cdot \nabla \Psi_T}{\Psi_T^2} + \dots \right) + (V - E_T)\frac{f}{\Psi_T} \right]
```

Multiply $`\Psi_T`$ through. The Kinetic Energy term $`-\frac{\hbar^2}{2m}`$ distributes over the expansion:

- **Diffusion Term:** $`-\frac{\hbar^2}{2m} \Psi_T (\frac{\nabla^2 f}{\Psi_T}) = -\frac{\hbar^2}{2m} \nabla^2 f`$

  > Note: The sign is flipped in the diffusion equation standard form on RHS, becoming positive coefficient.

- **Drift Term:** $`-\frac{\hbar^2}{2m} \Psi_T (-2\frac{\nabla f \cdot \nabla \Psi_T}{\Psi_T^2}) = \frac{\hbar^2}{m} \frac{\nabla \Psi_T}{\Psi_T} \cdot \nabla f = \nabla \cdot (\mathbf{v}_D f)`$

- **The Remainder (Local Energy):** The remaining terms involving only $`f`$ (no derivatives of $`f`$) are:

```math
\text{Remnant} = -\frac{\hbar^2}{2m} \Psi_T \left( f \left[ 2\frac{(\nabla \Psi_T)^2}{\Psi_T^3} - \frac{\nabla^2 \Psi_T}{\Psi_T^2} \right] \right) + V f
```

Through vector identities, the complex gradient terms simplify such that the equation becomes:

```math
\frac{\partial f}{\partial \tau} = \frac{\hbar}{2m} \nabla^2 f - \nabla \cdot (\mathbf{v}_D f) - \frac{1}{\hbar} \left[ \underbrace{-\frac{\hbar^2}{2m} \frac{\nabla^2 \Psi_T}{\Psi_T} + V(\mathbf{x})}_{E_L(\mathbf{x})} - E_T \right] f
```

Thus, $`E_L(\mathbf{x})`$ naturally emerges as the effective potential governing the growth/decay of the distribution $`f`$.

---

#### The Cusp Condition (Stability)

**Hypothesis:** The introduction of $`\Psi_T`$ leads to stability because $`E_L(\mathbf{x})`$ remains bounded (finite) even at Coulomb singularities (where $`V(\mathbf{x}) \to -\infty`$), provided $`\Psi_T`$ satisfies the Cusp Condition.

---

##### Proof Strategy

We examine the behavior of $`E_L(\mathbf{x})`$ for a Hydrogen-like atom (Nuclear charge $`Z`$) as the electron approaches the nucleus ($`r \to 0`$).

---

##### The Potential

The Coulomb potential is singular:

```math
V(r) = -\frac{Ze^2}{r} \quad (\text{In atomic units: } -\frac{Z}{r})
```

In pure DMC, the weight $`W = e^{-(V-E_T)\tau} \approx e^{+Z\tau/r} \to \infty`$ as $`r \to 0`$. This is the instability.

---

##### The Trial Wavefunction

Choose a trial function with the standard exponential cusp behavior:

```math
\Psi_T(r) = e^{-\alpha r}
```

where $`\alpha`$ is a variational parameter.

---

##### Calculating the Kinetic Term

We compute the Laplacian $`\nabla^2 \Psi_T`$ in spherical coordinates (radial part only, as angular derivatives are zero for s-orbitals):

```math
\nabla^2 = \frac{\partial^2}{\partial r^2} + \frac{2}{r} \frac{\partial}{\partial r}
```

Applying this to $`e^{-\alpha r}`$:

- $`\frac{\partial}{\partial r} e^{-\alpha r} = -\alpha e^{-\alpha r}`$
- $`\frac{\partial^2}{\partial r^2} e^{-\alpha r} = \alpha^2 e^{-\alpha r}`$

Substituting back:

```math
\nabla^2 \Psi_T = \left( \alpha^2 - \frac{2\alpha}{r} \right) e^{-\alpha r} = \left( \alpha^2 - \frac{2\alpha}{r} \right) \Psi_T
```

---

##### Constructing the Local Energy

Now compute $`E_L(r) = -\frac{1}{2}\frac{\nabla^2 \Psi_T}{\Psi_T} + V(r)`$ (in atomic units):

```math
E_L(r) = -\frac{1}{2} \left( \alpha^2 - \frac{2\alpha}{r} \right) - \frac{Z}{r}
```

```math
E_L(r) = -\frac{\alpha^2}{2} + \frac{\alpha}{r} - \frac{Z}{r}
```

Grouping the singular $`1/r`$ terms:

```math
E_L(r) = -\frac{\alpha^2}{2} + \frac{1}{r} (\alpha - Z)
```

---

##### The Stability (Cusp) Condition

To prevent divergence, the coefficient of the $`1/r`$ term must be zero.

```math
\alpha - Z = 0 \implies \alpha = Z
```

This condition ($`\alpha = Z`$) is the **Kato Cusp Condition**.

---

##### Resulting Local Energy

If we choose $`\Psi_T`$ satisfying this condition, the singularity vanishes exactly:

```math
E_L(r) = -\frac{Z^2}{2} \quad (\text{Constant finite value})
```

---

##### Conclusion on Stability

1. **Bounded Weights:** Because the kinetic singularity ($`+\frac{\alpha}{r}`$) cancels the potential singularity ($`-\frac{Z}{r}`$), $`E_L(r)`$ is finite everywhere. The branching weight $`W = e^{-(E_L - E_T)\tau}`$ is bounded, preventing numerical explosions.

2. **Reduced Variance:** Since $`E_L(\mathbf{x}) \approx \text{Constant}`$ (for a good $`\Psi_T`$), the fluctuations in weights are minimal. In the limit where $`\Psi_T`$ is the exact ground state, $`E_L(\mathbf{x}) = E_0`$ exactly, and the variance becomes zero. This is the **Zero-Variance Property** of Importance Sampled DMC.

---

### Implementation

This section outlines the step-by-step implementation of Importance Sampled DMC based on the theory developed above.

---

#### Step 1: Initialization

**Objective:** Create an initial population of $`N`$ "walkers" sampled from $`|\Psi_T(\mathbf{x})|^2`$.

**Math:** $`\mathbf{x}_i \sim |\Psi_T(\mathbf{x})|^2`$.

**Implementation:** Use Metropolis-Hastings or simple rejection sampling to distribute $`N`$ initial particle coordinates $`\{\mathbf{x}_1, \dots, \mathbf{x}_N\}`$ based on the trial probability density.

---

#### Step 2: Drift-Diffusion Update (The Kinetic Step)

**Objective:** Move walkers according to the kinetic energy (diffusion) and the trial function guidance (drift).

**Math (Green's Function with Drift):**

The propagator for the drift-diffusion part over time step $`\delta\tau`$ is:

```math
G(\mathbf{x} \leftarrow \mathbf{y}, \delta\tau) \approx \frac{1}{(2\pi \sigma^2)^{3/2}} \exp\left( -\frac{(\mathbf{x} - \mathbf{y} - \mathbf{v}_D(\mathbf{y})\delta\tau)^2}{2\sigma^2} \right)
```

where diffusion variance $`\sigma^2 = \frac{\hbar \delta\tau}{m}`$ and drift velocity $`\mathbf{v}_D(\mathbf{y}) = \frac{\hbar}{m} \nabla \ln \Psi_T(\mathbf{y})`$.

**Implementation:**

For each walker $`i`$:

1. **Calculate Drift:** Compute the gradient of the log trial function at the current position.
   ```python
   drift_vector = (hbar / mass) * grad_log_psi_T(x_old)
   ```

2. **Generate Random Displacement:** Sample a Gaussian random vector $`\boldsymbol{\eta}`$ with variance $`\frac{\hbar \delta\tau}{m}`$.
   ```python
   eta = sqrt(hbar * dt / mass) * random_normal_vector()
   ```

3. **Update Position:**
   ```python
   x_new = x_old + drift_vector * dt + eta
   ```

> **Note:** The drift pushes walkers away from regions where $`\Psi_T`$ is small (nodes) and toward regions where $`\Psi_T`$ is large, preventing them from exploring irrelevant space.

---

#### Step 3: Metropolis Acceptance/Rejection

**Objective:** Correct for the asymmetry of the drift-diffusion Green's function to satisfy detailed balance.

As discussed in the [Theory section](#asymmetry-of-the-greens-function), the Green's function with drift is **asymmetric**: $`G(\mathbf{x}' \leftarrow \mathbf{x}) \neq G(\mathbf{x} \leftarrow \mathbf{x}')`$. This means naive sampling violates detailed balance.

**Solution:** Accept proposed moves with probability:

```math
A(\mathbf{x}' \leftarrow \mathbf{x}) = \min \left( 1, \frac{|\Psi_T(\mathbf{x}')|^2 G(\mathbf{x} \leftarrow \mathbf{x}')}{|\Psi_T(\mathbf{x})|^2 G(\mathbf{x}' \leftarrow \mathbf{x})} \right)
```

**Explicit Form:**

```math
A = \min \left( 1, \frac{|\Psi_T(\mathbf{x}')|^2}{|\Psi_T(\mathbf{x})|^2} \exp\left[ -\frac{|\mathbf{x} - \mathbf{x}' - \mathbf{v}_D(\mathbf{x}')\delta\tau|^2 - |\mathbf{x}' - \mathbf{x} - \mathbf{v}_D(\mathbf{x})\delta\tau|^2}{2\sigma^2} \right] \right)
```

**Implementation:**

```python
# After proposing move in Step 2: x_old -> x_proposed
psi_ratio = (psi_T(x_proposed) / psi_T(x_old))**2
G_forward = exp(-|x_proposed - x_old - drift(x_old)*dt|**2 / (2*sigma**2))
G_backward = exp(-|x_old - x_proposed - drift(x_proposed)*dt|**2 / (2*sigma**2))
A = min(1, psi_ratio * G_backward / G_forward)

if random_uniform(0, 1) < A:
    x_new = x_proposed  # Accept
else:
    x_new = x_old       # Reject: keep old position
```

> **Physical Interpretation:** The Metropolis step ensures walkers sample the correct equilibrium distribution $`|\Psi_T|^2`$. It compensates for the directional bias introduced by the position-dependent drift.

---

#### Step 4: Branching (The Potential Step)

**Objective:** Reweight the population based on the "Local Energy" rather than the raw potential.

**Math (Weight Factor):**

The reaction term is governed by the Local Energy $`E_L(\mathbf{x})`$:

```math
E_L(\mathbf{x}) = \frac{\hat{H} \Psi_T(\mathbf{x})}{\Psi_T(\mathbf{x})} = -\frac{\hbar^2}{2m} \frac{\nabla^2 \Psi_T}{\Psi_T} + V(\mathbf{x})
```

The weight factor for a step $`\delta\tau`$ is:

```math
W = \exp\left( -\frac{\delta\tau}{\hbar} \left[ \frac{E_L(\mathbf{x}_{new}) + E_L(\mathbf{x}_{old})}{2} - E_T \right] \right)
```

> **Note:** Averaging $`E_L`$ at old and new positions reduces time-step error to $`O(\delta\tau^2)`$.

**Implementation:**

For each walker $`i`$:

1. **Evaluate Local Energy:**
   ```python
   E_L_old = H_psi(x_old) / psi(x_old)
   E_L_new = H_psi(x_new) / psi(x_new)
   E_L_avg = 0.5 * (E_L_old + E_L_new)
   ```

2. **Calculate Branching Factor:**
   ```python
   weight = exp(-(dt/hbar) * (E_L_avg - E_T))
   ```

3. **Stochastic Branching (Birth/Death):**
   Use the integer part of the weight plus a random number to determine the number of copies ($`m`$) for the next step.
   ```python
   m = int(weight + random_uniform(0, 1))
   ```
   - If `m == 0`: Delete the walker (Death).
   - If `m == 1`: Keep the walker (Survival).
   - If `m > 1`: Create `m-1` copies of the walker at `x_new` (Cloning).

---

#### Step 5: Updating $`E_T`$ (Population Control)

**Objective:** Adjust the reference energy to keep the population stable near target size $`N_0`$.

**Math (Feedback Law):**

```math
E_T(\tau) = \langle E_L \rangle_{\tau} - \alpha \ln\left( \frac{N(\tau)}{N_0} \right)
```

**Implementation:**

At the end of every block (or step):

1. **Calculate Average Energy:**
   ```python
   E_est = sum(E_L for all walkers) / current_population_size
   ```

2. **Update Trial Energy:**
   ```python
   E_T = E_est - alpha * log(current_population_size / N_target)
   ```
