# Diffusion Monte Carlo

This repository estimates the ground state of three well-known Bosonic systems: the Hydrogen ($`H`$), Hydrogen ion ($`H_2^+`$), and Hydrogen molecule ($`H_2`$) via Diffusion Monte Carlo (DMC).

---

## Table of Contents

- [Theoretical Background](#theoretical-background)
  - [Time-Dependent Schrödinger Equation](#time-dependent-schrödinger-equation)
  - [Original Definition with the Auxiliary Constant Term E_T](#original-definition-with-the-auxiliary-constant-term-e_t)
  - [Explicit Expression of the Kinetic Operator](#explicit-expression-of-the-kinetic-operator)
  - [Rewriting the Equation Using Imaginary Time](#rewriting-the-equation-using-imaginary-time)
  - [Schrödinger Equation as the Diffusion Equation](#schrödinger-equation-as-the-diffusion-equation)
  - [From Perron-Frobenius to Krein-Rutman Theorem](#from-perron-frobenius-to-krein-rutman-theorem)
  - [Mathematical vs. Quantum Mechanical Ground State](#mathematical-vs-quantum-mechanical-ground-state)
  - [Proof of the Gaussian Green's Function (Free Diffusion)](#proof-of-the-gaussian-greens-function-free-diffusion)
  - [Proof of Invariance of Total Number of Particles](#proof-of-invariance-of-total-number-of-particles)
- [Integral Formulation of the Imaginary-Time Schrödinger Equation](#integral-formulation-of-the-imaginary-time-schrödinger-equation)
  - [Fourier Transforms and Momentum States](#fourier-transforms-and-momentum-states)
  - [Trotter-Suzuki Decomposition](#trotter-suzuki-decomposition)
  - [Expansion in Position Space](#expansion-in-position-space)
  - [Proof of Kinetic Operator as Green's Function](#proof-of-kinetic-operator-as-greens-function)
  - [Stationarity and Eigenfunction Correspondence](#stationarity-and-eigenfunction-correspondence)
  - [Convergence to the Ground State](#convergence-to-the-ground-state)
  - [Time Derivative of Number of Particles](#time-derivative-of-number-of-particles)
  - [Population Control Bias and Feedback Law](#population-control-bias-and-feedback-law)
- [Pure Diffusion Monte Carlo](#pure-diffusion-monte-carlo)
  - [Theory](#theory)
  - [Implementation](#implementation)
- [Importance Sampled Diffusion Monte Carlo](#importance-sampled-diffusion-monte-carlo)
  - [Motivation: Why Importance Sampling?](#motivation-why-importance-sampling)
  - [Theory](#theory-1)
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

We derive the Green's function (fundamental solution) for the free diffusion equation, showing it is a Gaussian whose width grows with imaginary time.

---

#### Step 1: The Free Diffusion Equation

Consider the imaginary-time Schrödinger equation without potential ($`V(\mathbf{x}) = 0`$) and with $`E_T = 0`$:

```math
\frac{\partial \Psi(\mathbf{x}, \tau)}{\partial \tau} = D \nabla_{\mathbf{x}}^2 \Psi(\mathbf{x}, \tau)
```

where the diffusion coefficient is $`D = \frac{\hbar}{2m}`$.

**Goal:** Find the Green's function $`G_0(\mathbf{x}, \mathbf{y}; \tau)`$ satisfying:

```math
\Psi(\mathbf{x}, \tau) = \int_{\mathbb{R}^3} d\mathbf{y} \, G_0(\mathbf{x} - \mathbf{y}, \tau) \, \Psi(\mathbf{y}, 0)
```

with the initial condition $`G_0(\mathbf{x}, 0) = \delta^{(3)}(\mathbf{x})`$ (Dirac delta function).

---

#### Step 2: Fourier Transform Definition

Define the 3D Fourier transform and its inverse:

```math
\tilde{\Psi}(\mathbf{k}, \tau) = \mathcal{F}[\Psi] = \int_{\mathbb{R}^3} d\mathbf{x} \, e^{-i\mathbf{k}\cdot\mathbf{x}} \Psi(\mathbf{x}, \tau)
```

```math
\Psi(\mathbf{x}, \tau) = \mathcal{F}^{-1}[\tilde{\Psi}] = \frac{1}{(2\pi)^3} \int_{\mathbb{R}^3} d\mathbf{k} \, e^{i\mathbf{k}\cdot\mathbf{x}} \tilde{\Psi}(\mathbf{k}, \tau)
```

**Key Property:** The Fourier transform converts spatial derivatives to algebraic multiplication:

```math
\mathcal{F}[\nabla^2 \Psi] = -|\mathbf{k}|^2 \tilde{\Psi}(\mathbf{k}, \tau)
```

---

#### Step 3: Transforming the PDE to an ODE

Apply the Fourier transform to both sides of the diffusion equation:

```math
\mathcal{F}\left[ \frac{\partial \Psi}{\partial \tau} \right] = D \cdot \mathcal{F}\left[ \nabla^2 \Psi \right]
```

Since differentiation with respect to $`\tau`$ commutes with the spatial Fourier transform:

```math
\frac{\partial}{\partial \tau} \tilde{\Psi}(\mathbf{k}, \tau) = D \cdot (-|\mathbf{k}|^2) \tilde{\Psi}(\mathbf{k}, \tau)
```

This is now an ordinary differential equation (ODE) in $`\tau`$ for each fixed $`\mathbf{k}`$:

```math
\frac{\partial \tilde{\Psi}}{\partial \tau} = -D k^2 \tilde{\Psi}
```

where $`k^2 = |\mathbf{k}|^2 = k_x^2 + k_y^2 + k_z^2`$.

---

#### Step 4: Solving the ODE in Fourier Space

The ODE $`\frac{d\tilde{\Psi}}{d\tau} = -Dk^2 \tilde{\Psi}`$ is a first-order linear equation with solution:

```math
\tilde{\Psi}(\mathbf{k}, \tau) = \tilde{\Psi}(\mathbf{k}, 0) \exp\left( -D k^2 \tau \right)
```

Substituting $`D = \frac{\hbar}{2m}`$:

```math
\tilde{\Psi}(\mathbf{k}, \tau) = \tilde{\Psi}(\mathbf{k}, 0) \exp\left( -\frac{\hbar k^2}{2m} \tau \right)
```

**Interpretation:** Each Fourier mode $`\mathbf{k}`$ decays exponentially. High-frequency (large $`|\mathbf{k}|`$) components decay faster, corresponding to the smoothing effect of diffusion.

---

#### Step 5: Applying the Convolution Theorem

The solution in Fourier space is a product:

```math
\tilde{\Psi}(\mathbf{k}, \tau) = \tilde{G}_0(\mathbf{k}, \tau) \cdot \tilde{\Psi}(\mathbf{k}, 0)
```

where $`\tilde{G}_0(\mathbf{k}, \tau) = \exp\left( -\frac{\hbar k^2}{2m} \tau \right)`$ is the Fourier transform of the Green's function.

**Convolution Theorem:** Multiplication in Fourier space corresponds to convolution in real space:

```math
\mathcal{F}^{-1}[\tilde{f} \cdot \tilde{g}] = f * g = \int_{\mathbb{R}^3} d\mathbf{y} \, f(\mathbf{x} - \mathbf{y}) g(\mathbf{y})
```

Therefore:

```math
\Psi(\mathbf{x}, \tau) = \int_{\mathbb{R}^3} d\mathbf{y} \, G_0(\mathbf{x} - \mathbf{y}, \tau) \, \Psi(\mathbf{y}, 0)
```

---

#### Step 6: Computing the Inverse Fourier Transform of $`\tilde{G}_0`$

We must evaluate:

```math
G_0(\mathbf{x}, \tau) = \frac{1}{(2\pi)^3} \int_{\mathbb{R}^3} d\mathbf{k} \, e^{i\mathbf{k}\cdot\mathbf{x}} \exp\left( -\frac{\hbar \tau}{2m} k^2 \right)
```

Define $`\alpha = \frac{\hbar \tau}{2m}`$ for convenience. The integral factors into three independent 1D integrals:

```math
G_0(\mathbf{x}, \tau) = \frac{1}{(2\pi)^3} \prod_{j \in \{x,y,z\}} \int_{-\infty}^{\infty} dk_j \, e^{i k_j x_j} e^{-\alpha k_j^2}
```

Each 1D integral is a standard Gaussian integral of the form:

```math
I_j = \int_{-\infty}^{\infty} dk_j \, e^{-\alpha k_j^2 + i k_j x_j}
```

---

#### Step 7: Evaluating the Gaussian Integral (Completing the Square)

**Standard Result:** For $`\text{Re}(\alpha) > 0`$:

```math
\int_{-\infty}^{\infty} dk \, e^{-\alpha k^2 + \beta k} = \sqrt{\frac{\pi}{\alpha}} \exp\left( \frac{\beta^2}{4\alpha} \right)
```

**Derivation by Completing the Square:**

```math
-\alpha k^2 + \beta k = -\alpha \left( k^2 - \frac{\beta}{\alpha} k \right) = -\alpha \left( k - \frac{\beta}{2\alpha} \right)^2 + \frac{\beta^2}{4\alpha}
```

Substituting $`u = k - \frac{\beta}{2\alpha}`$ (shift of integration variable):

```math
\int_{-\infty}^{\infty} dk \, e^{-\alpha k^2 + \beta k} = e^{\frac{\beta^2}{4\alpha}} \int_{-\infty}^{\infty} du \, e^{-\alpha u^2} = e^{\frac{\beta^2}{4\alpha}} \sqrt{\frac{\pi}{\alpha}}
```

**Application:** With $`\beta = i x_j`$:

```math
I_j = \sqrt{\frac{\pi}{\alpha}} \exp\left( \frac{(ix_j)^2}{4\alpha} \right) = \sqrt{\frac{\pi}{\alpha}} \exp\left( -\frac{x_j^2}{4\alpha} \right)
```

---

#### Step 8: Assembling the Full Green's Function

Combining all three dimensions:

```math
G_0(\mathbf{x}, \tau) = \frac{1}{(2\pi)^3} \left( \sqrt{\frac{\pi}{\alpha}} \right)^3 \exp\left( -\frac{x^2 + y^2 + z^2}{4\alpha} \right)
```

```math
G_0(\mathbf{x}, \tau) = \frac{1}{(2\pi)^3} \left( \frac{\pi}{\alpha} \right)^{3/2} \exp\left( -\frac{|\mathbf{x}|^2}{4\alpha} \right)
```

Substituting back $`\alpha = \frac{\hbar \tau}{2m}`$:

```math
G_0(\mathbf{x}, \tau) = \frac{1}{(2\pi)^3} \left( \frac{2\pi m}{\hbar \tau} \right)^{3/2} \exp\left( -\frac{m |\mathbf{x}|^2}{2\hbar \tau} \right)
```

Simplifying:

```math
G_0(\mathbf{x}, \tau) = \left( \frac{m}{2\pi \hbar \tau} \right)^{3/2} \exp\left( -\frac{m |\mathbf{x}|^2}{2\hbar \tau} \right)
```

---

#### Step 9: Verification of Normalization

A valid probability density must integrate to unity. We verify:

```math
\int_{\mathbb{R}^3} d\mathbf{x} \, G_0(\mathbf{x}, \tau) = \left( \frac{m}{2\pi \hbar \tau} \right)^{3/2} \int_{\mathbb{R}^3} d\mathbf{x} \, \exp\left( -\frac{m |\mathbf{x}|^2}{2\hbar \tau} \right)
```

The integral factors into three 1D Gaussian integrals:

```math
\int_{-\infty}^{\infty} dx \, e^{-\frac{m x^2}{2\hbar\tau}} = \sqrt{\frac{2\pi\hbar\tau}{m}}
```

Therefore:

```math
\int_{\mathbb{R}^3} d\mathbf{x} \, G_0(\mathbf{x}, \tau) = \left( \frac{m}{2\pi \hbar \tau} \right)^{3/2} \cdot \left( \frac{2\pi\hbar\tau}{m} \right)^{3/2} = 1 \quad \checkmark
```

---

#### Step 10: Identification as a Gaussian Distribution

The Green's function is a multivariate Gaussian with:

- **Mean:** $`\boldsymbol{\mu} = \mathbf{0}`$ (centered at origin)
- **Covariance matrix:** $`\Sigma = \sigma^2 I_3`$ where $`\sigma^2 = \frac{\hbar \tau}{m}`$

In standard form:

```math
G_0(\mathbf{x}, \tau) = \frac{1}{(2\pi\sigma^2)^{3/2}} \exp\left( -\frac{|\mathbf{x}|^2}{2\sigma^2} \right)
```

**Physical Interpretation:**
- The variance $`\sigma^2 = \frac{\hbar \tau}{m}`$ grows linearly with imaginary time $`\tau`$
- This is characteristic of diffusive (Brownian) motion: $`\langle |\mathbf{x}|^2 \rangle = 3\sigma^2 = \frac{3\hbar\tau}{m}`$
- Heavier particles ($`m \uparrow`$) diffuse more slowly
- The Green's function represents the probability of a particle diffusing from $`\mathbf{y}`$ to $`\mathbf{x}`$ in imaginary time $`\tau`$

---

### Proof of Invariance of Total Number of Particles

We prove that diffusion alone conserves the total "particle number" (integral of $`\Psi`$), while the potential term acts as a source or sink.

---

#### Step 1: Definition of Total Particle Number

Define the total particle number as the integral of the wavefunction over all space:

```math
N(\tau) \equiv \int_{\mathbb{R}^3} d\mathbf{x} \, \Psi(\mathbf{x}, \tau)
```

---

#### Step 2: Differentiation Under the Integral Sign (Leibniz Rule)

To compute $`\frac{dN}{d\tau}`$, we interchange differentiation and integration:

```math
\frac{dN}{d\tau} = \frac{d}{d\tau} \int_{\mathbb{R}^3} d\mathbf{x} \, \Psi(\mathbf{x}, \tau) = \int_{\mathbb{R}^3} d\mathbf{x} \, \frac{\partial \Psi(\mathbf{x}, \tau)}{\partial \tau}
```

**Justification (Leibniz Integral Rule):** This interchange is valid when:
1. $`\Psi(\mathbf{x}, \tau)`$ and $`\frac{\partial \Psi}{\partial \tau}`$ are continuous on $`\mathbb{R}^3 \times [0, T]`$
2. There exists an integrable dominating function $`g(\mathbf{x})`$ such that $`\left| \frac{\partial \Psi}{\partial \tau} \right| \leq g(\mathbf{x})`$ for all $`\tau`$

For physically reasonable wavefunctions that decay exponentially at infinity (e.g., bound states), these conditions are satisfied.

---

#### Step 3: Substitution of the Imaginary-Time Schrödinger Equation

Recall the imaginary-time Schrödinger equation in position representation:

```math
\frac{\partial \Psi(\mathbf{x}, \tau)}{\partial \tau} = \frac{\hbar}{2m} \nabla_{\mathbf{x}}^2 \Psi(\mathbf{x}, \tau) - \frac{1}{\hbar}(V(\mathbf{x}) - E_T) \Psi(\mathbf{x}, \tau)
```

Substituting into the integral:

```math
\frac{dN}{d\tau} = \int_{\mathbb{R}^3} d\mathbf{x} \left[ \frac{\hbar}{2m} \nabla_{\mathbf{x}}^2 \Psi(\mathbf{x}, \tau) - \frac{1}{\hbar}(V(\mathbf{x}) - E_T) \Psi(\mathbf{x}, \tau) \right]
```

By linearity of integration, we split this into two terms:

```math
\frac{dN}{d\tau} = \underbrace{\frac{\hbar}{2m} \int_{\mathbb{R}^3} d\mathbf{x} \, \nabla_{\mathbf{x}}^2 \Psi(\mathbf{x}, \tau)}_{I_{\text{kinetic}}} - \underbrace{\frac{1}{\hbar} \int_{\mathbb{R}^3} d\mathbf{x} \, (V(\mathbf{x}) - E_T) \Psi(\mathbf{x}, \tau)}_{I_{\text{potential}}}
```

---

#### Step 4: Evaluation of the Kinetic Term via the Divergence Theorem

The Laplacian can be written as the divergence of a gradient:

```math
\nabla^2 \Psi = \nabla \cdot (\nabla \Psi)
```

Applying the Divergence Theorem (Gauss's Theorem) to convert the volume integral to a surface integral:

```math
I_{\text{kinetic}} = \frac{\hbar}{2m} \int_{\mathbb{R}^3} d\mathbf{x} \, \nabla \cdot (\nabla \Psi) = \frac{\hbar}{2m} \lim_{R \to \infty} \oint_{S_R} (\nabla \Psi) \cdot d\mathbf{S}
```

where $`S_R`$ is a sphere of radius $`R`$ centered at the origin.

**Boundary Condition:** For bound states, $`\Psi(\mathbf{x}, \tau)`$ decays exponentially as $`|\mathbf{x}| \to \infty`$:

```math
\Psi(\mathbf{x}, \tau) \sim e^{-\kappa |\mathbf{x}|} \quad \text{as } |\mathbf{x}| \to \infty
```

for some $`\kappa > 0`$. Consequently:

```math
|\nabla \Psi| \sim \kappa e^{-\kappa |\mathbf{x}|}
```

The surface integral over $`S_R`$ scales as:

```math
\left| \oint_{S_R} (\nabla \Psi) \cdot d\mathbf{S} \right| \leq \max_{S_R} |\nabla \Psi| \cdot 4\pi R^2 \sim \kappa e^{-\kappa R} \cdot R^2 \xrightarrow{R \to \infty} 0
```

Therefore:

```math
I_{\text{kinetic}} = 0
```

---

#### Step 5: Result for Pure Diffusion ($`V = 0`$)

In the absence of a potential ($`V(\mathbf{x}) = 0`$) and with $`E_T = 0`$:

```math
\frac{dN}{d\tau} = I_{\text{kinetic}} - I_{\text{potential}} = 0 - 0 = 0
```

**Conclusion:** Pure diffusion conserves the total particle number. This reflects the fact that diffusion merely redistributes probability density without creating or destroying it.

---

#### Step 6: Result with Potential (Source/Sink Term)

With a non-zero potential, only the kinetic term vanishes:

```math
\frac{dN}{d\tau} = 0 - \frac{1}{\hbar} \int_{\mathbb{R}^3} d\mathbf{x} \, (V(\mathbf{x}) - E_T) \Psi(\mathbf{x}, \tau)
```

Rearranging:

```math
\frac{dN}{d\tau} = -\frac{1}{\hbar} \left[ \int_{\mathbb{R}^3} d\mathbf{x} \, V(\mathbf{x}) \Psi(\mathbf{x}, \tau) - E_T \int_{\mathbb{R}^3} d\mathbf{x} \, \Psi(\mathbf{x}, \tau) \right]
```

```math
\frac{dN}{d\tau} = -\frac{1}{\hbar} \left[ \langle V \rangle_\Psi \cdot N(\tau) - E_T \cdot N(\tau) \right] = -\frac{N(\tau)}{\hbar} \left( \langle V \rangle_\Psi - E_T \right)
```

where $`\langle V \rangle_\Psi = \frac{\int d\mathbf{x} \, V(\mathbf{x}) \Psi(\mathbf{x}, \tau)}{\int d\mathbf{x} \, \Psi(\mathbf{x}, \tau)}`$ is the average potential weighted by $`\Psi`$.

**Interpretation:**
- If $`\langle V \rangle_\Psi > E_T`$: $`\frac{dN}{d\tau} < 0`$ (particle number decreases — net sink)
- If $`\langle V \rangle_\Psi < E_T`$: $`\frac{dN}{d\tau} > 0`$ (particle number increases — net source)
- If $`\langle V \rangle_\Psi = E_T`$: $`\frac{dN}{d\tau} = 0`$ (steady state)

---

## Integral Formulation of the Imaginary-Time Schrödinger Equation

This section develops the integral (path-integral-like) formulation of the imaginary-time Schrödinger equation, which forms the mathematical foundation for DMC algorithms.

---

### Fourier Transforms and Momentum States

We establish the mathematical machinery needed to evaluate the kinetic energy propagator in position space.

---

#### The Momentum Operator

In quantum mechanics, the momentum operator in position representation is:

```math
\hat{\mathbf{p}} = -i\hbar \nabla_{\mathbf{x}}
```

We seek eigenstates $`| \mathbf{k} \rangle`$ satisfying:

```math
\hat{\mathbf{p}} | \mathbf{k} \rangle = \hbar \mathbf{k} | \mathbf{k} \rangle
```

where $`\mathbf{k}`$ is the wavevector and $`\hbar \mathbf{k}`$ is the momentum eigenvalue.

---

#### Derivation of the Position-Momentum Overlap

To find $`\langle \mathbf{x} | \mathbf{k} \rangle`$, we solve the eigenvalue equation in position representation:

```math
\langle \mathbf{x} | \hat{\mathbf{p}} | \mathbf{k} \rangle = \hbar \mathbf{k} \langle \mathbf{x} | \mathbf{k} \rangle
```

Using $`\hat{\mathbf{p}} = -i\hbar \nabla`$ acting on the position basis:

```math
-i\hbar \nabla_{\mathbf{x}} \langle \mathbf{x} | \mathbf{k} \rangle = \hbar \mathbf{k} \langle \mathbf{x} | \mathbf{k} \rangle
```

This is a first-order differential equation:

```math
\nabla_{\mathbf{x}} \langle \mathbf{x} | \mathbf{k} \rangle = i \mathbf{k} \langle \mathbf{x} | \mathbf{k} \rangle
```

The solution is a plane wave:

```math
\langle \mathbf{x} | \mathbf{k} \rangle = C \, e^{i\mathbf{k}\cdot\mathbf{x}}
```

where $`C`$ is a normalization constant to be determined.

---

#### Normalization via Delta Function

We require orthonormality in the continuum sense:

```math
\langle \mathbf{k} | \mathbf{k}' \rangle = \delta^{(3)}(\mathbf{k} - \mathbf{k}')
```

Insert the position-space resolution of identity $`\hat{I} = \int d\mathbf{x} | \mathbf{x} \rangle \langle \mathbf{x} |`$:

```math
\langle \mathbf{k} | \mathbf{k}' \rangle = \int d\mathbf{x} \, \langle \mathbf{k} | \mathbf{x} \rangle \langle \mathbf{x} | \mathbf{k}' \rangle = |C|^2 \int d\mathbf{x} \, e^{-i\mathbf{k}\cdot\mathbf{x}} e^{i\mathbf{k}'\cdot\mathbf{x}}
```

```math
= |C|^2 \int d\mathbf{x} \, e^{i(\mathbf{k}'-\mathbf{k})\cdot\mathbf{x}}
```

Using the Fourier representation of the Dirac delta:

```math
\int d\mathbf{x} \, e^{i\mathbf{q}\cdot\mathbf{x}} = (2\pi)^3 \delta^{(3)}(\mathbf{q})
```

we obtain:

```math
\langle \mathbf{k} | \mathbf{k}' \rangle = |C|^2 (2\pi)^3 \delta^{(3)}(\mathbf{k}' - \mathbf{k})
```

For this to equal $`\delta^{(3)}(\mathbf{k} - \mathbf{k}')`$, we require:

```math
|C|^2 = \frac{1}{(2\pi)^3} \implies C = \frac{1}{(2\pi)^{3/2}}
```

**Result:**

```math
\langle \mathbf{x} | \mathbf{k} \rangle = \frac{1}{(2\pi)^{3/2}} e^{i\mathbf{k}\cdot\mathbf{x}}
```

---

#### Completeness Relation (Resolution of Identity)

The momentum states form a complete basis:

```math
\hat{I} = \int d\mathbf{k} \, | \mathbf{k} \rangle \langle \mathbf{k} |
```

**Verification:** Apply to an arbitrary state $`| \psi \rangle`$ and project onto position:

```math
\langle \mathbf{x} | \hat{I} | \psi \rangle = \int d\mathbf{k} \, \langle \mathbf{x} | \mathbf{k} \rangle \langle \mathbf{k} | \psi \rangle
```

Define the momentum-space wavefunction $`\tilde{\psi}(\mathbf{k}) = \langle \mathbf{k} | \psi \rangle`$:

```math
\psi(\mathbf{x}) = \int d\mathbf{k} \, \frac{e^{i\mathbf{k}\cdot\mathbf{x}}}{(2\pi)^{3/2}} \tilde{\psi}(\mathbf{k}) = \frac{1}{(2\pi)^{3/2}} \int d\mathbf{k} \, e^{i\mathbf{k}\cdot\mathbf{x}} \tilde{\psi}(\mathbf{k})
```

This is precisely the inverse Fourier transform, confirming completeness.

---

#### Connection to Fourier Transforms

The position-momentum overlap establishes the Fourier transform pair:

| **Transform** | **Definition** |
|:--------------|:---------------|
| Forward (position → momentum) | $`\tilde{\psi}(\mathbf{k}) = \langle \mathbf{k} \| \psi \rangle = \frac{1}{(2\pi)^{3/2}} \int d\mathbf{x} \, e^{-i\mathbf{k}\cdot\mathbf{x}} \psi(\mathbf{x})`$ |
| Inverse (momentum → position) | $`\psi(\mathbf{x}) = \langle \mathbf{x} \| \psi \rangle = \frac{1}{(2\pi)^{3/2}} \int d\mathbf{k} \, e^{i\mathbf{k}\cdot\mathbf{x}} \tilde{\psi}(\mathbf{k})`$ |

> **Note:** This symmetric convention places factors of $`(2\pi)^{-3/2}`$ in both transforms. An alternative convention places $`(2\pi)^{-3}`$ entirely in the inverse transform.

---

### Trotter-Suzuki Decomposition

The Trotter-Suzuki decomposition is essential for splitting the propagator into manageable kinetic and potential parts.

---

#### Formal Solution of the Imaginary-Time Equation

The imaginary-time Schrödinger equation:

```math
\frac{\partial}{\partial \tau} | \Psi(\tau) \rangle = -\frac{1}{\hbar}(\hat{H} - E_T) | \Psi(\tau) \rangle
```

has the formal solution over a time interval $`\delta\tau`$:

```math
| \Psi(\tau + \delta\tau) \rangle = e^{-\frac{\delta\tau}{\hbar}(\hat{H} - E_T)} | \Psi(\tau) \rangle = e^{-\frac{\delta\tau}{\hbar}(\hat{T} + \hat{V} - E_T)} | \Psi(\tau) \rangle
```

The propagator $`\hat{U}(\delta\tau) = e^{-\frac{\delta\tau}{\hbar}(\hat{T} + \hat{V} - E_T)}`$ involves the exponential of a sum of non-commuting operators.

---

#### The Non-Commutativity Problem

For operators $`\hat{A}`$ and $`\hat{B}`$ that do not commute ($`[\hat{A}, \hat{B}] \neq 0`$), in general:

```math
e^{\hat{A} + \hat{B}} \neq e^{\hat{A}} e^{\hat{B}}
```

The kinetic and potential operators do not commute:

```math
[\hat{T}, \hat{V}] = \left[ -\frac{\hbar^2}{2m}\nabla^2, V(\mathbf{x}) \right] \neq 0
```

This is because $`\hat{T}`$ involves derivatives that act on the position-dependent $`V(\mathbf{x})`$.

---

#### Baker-Campbell-Hausdorff Formula

The relationship between $`e^{\hat{A}+\hat{B}}`$ and $`e^{\hat{A}}e^{\hat{B}}`$ is given by the Baker-Campbell-Hausdorff (BCH) formula:

```math
e^{\hat{A}} e^{\hat{B}} = \exp\left( \hat{A} + \hat{B} + \frac{1}{2}[\hat{A}, \hat{B}] + \frac{1}{12}[\hat{A}, [\hat{A}, \hat{B}]] - \frac{1}{12}[\hat{B}, [\hat{A}, \hat{B}]] + \cdots \right)
```

Inverting this relationship:

```math
e^{\hat{A} + \hat{B}} = e^{\hat{A}} e^{\hat{B}} e^{-\frac{1}{2}[\hat{A}, \hat{B}]} \cdots
```

---

#### First-Order (Lie-Trotter) Decomposition

Let $`\hat{A} = -\frac{\delta\tau}{\hbar}\hat{T}`$ and $`\hat{B} = -\frac{\delta\tau}{\hbar}(\hat{V} - E_T)`$. Then:

```math
[\hat{A}, \hat{B}] = \frac{\delta\tau^2}{\hbar^2} [\hat{T}, \hat{V} - E_T] = \frac{\delta\tau^2}{\hbar^2} [\hat{T}, \hat{V}] = O(\delta\tau^2)
```

Since the commutator is already $`O(\delta\tau^2)`$, we can write:

```math
e^{-\frac{\delta\tau}{\hbar}(\hat{T} + \hat{V} - E_T)} = e^{-\frac{\delta\tau}{\hbar}(\hat{V} - E_T)} e^{-\frac{\delta\tau}{\hbar}\hat{T}} + O(\delta\tau^2)
```

This is the **first-order Trotter decomposition** (also called Lie-Trotter splitting).

**Physical Interpretation:** We approximate the simultaneous evolution under $`\hat{T}`$ and $`\hat{V}`$ by:
1. First evolving under $`\hat{T}`$ alone (diffusion step)
2. Then evolving under $`\hat{V} - E_T`$ alone (branching step)

The error from ignoring the non-commutativity is $`O(\delta\tau^2)`$ per step.

---

#### Second-Order (Strang) Decomposition

A more accurate **symmetric splitting** achieves $`O(\delta\tau^3)`$ error:

```math
e^{-\frac{\delta\tau}{\hbar}(\hat{T} + \hat{V} - E_T)} = e^{-\frac{\delta\tau}{2\hbar}(\hat{V} - E_T)} e^{-\frac{\delta\tau}{\hbar}\hat{T}} e^{-\frac{\delta\tau}{2\hbar}(\hat{V} - E_T)} + O(\delta\tau^3)
```

This "potential-kinetic-potential" or "VTV" scheme is sometimes used in production DMC codes for improved accuracy.

---

#### Error Accumulation

For a total imaginary time $`\tau = N \cdot \delta\tau`$:

| **Decomposition** | **Error per Step** | **Total Error** |
|:------------------|:-------------------|:----------------|
| First-order (TV) | $`O(\delta\tau^2)`$ | $`O(\tau \cdot \delta\tau) = O(\delta\tau)`$ |
| Second-order (VTV) | $`O(\delta\tau^3)`$ | $`O(\tau \cdot \delta\tau^2) = O(\delta\tau^2)`$ |

The time-step bias in DMC energy estimates originates from this Trotter error.

---

### Expansion in Position Space

We now project the decomposed propagator onto the position basis to obtain an integral equation.

---

#### Step 1: Projection onto Position Basis

Starting from the Trotter-decomposed evolution:

```math
| \Psi(\tau+\delta\tau) \rangle = e^{-\frac{\delta\tau}{\hbar}(\hat{V} - E_T)} e^{-\frac{\delta\tau}{\hbar}\hat{T}} | \Psi(\tau) \rangle + O(\delta\tau^2)
```

Project onto position eigenstate $`\langle \mathbf{x} |`$:

```math
\Psi(\mathbf{x}, \tau+\delta\tau) = \langle \mathbf{x} | e^{-\frac{\delta\tau}{\hbar}(\hat{V} - E_T)} e^{-\frac{\delta\tau}{\hbar}\hat{T}} | \Psi(\tau) \rangle + O(\delta\tau^2)
```

---

#### Step 2: Potential Operator is Diagonal

The potential operator $`\hat{V}`$ is diagonal in the position basis:

```math
\hat{V} | \mathbf{x} \rangle = V(\mathbf{x}) | \mathbf{x} \rangle
```

Therefore, any function of $`\hat{V}`$ is also diagonal:

```math
e^{-\frac{\delta\tau}{\hbar}(\hat{V} - E_T)} | \mathbf{x} \rangle = e^{-\frac{\delta\tau}{\hbar}(V(\mathbf{x}) - E_T)} | \mathbf{x} \rangle
```

Taking the adjoint (since $`\hat{V}`$ is Hermitian):

```math
\langle \mathbf{x} | e^{-\frac{\delta\tau}{\hbar}(\hat{V} - E_T)} = e^{-\frac{\delta\tau}{\hbar}(V(\mathbf{x}) - E_T)} \langle \mathbf{x} |
```

**Result:** The potential exponential becomes a multiplicative factor:

```math
\Psi(\mathbf{x}, \tau+\delta\tau) = e^{-\frac{\delta\tau}{\hbar}(V(\mathbf{x}) - E_T)} \langle \mathbf{x} | e^{-\frac{\delta\tau}{\hbar}\hat{T}} | \Psi(\tau) \rangle + O(\delta\tau^2)
```

---

#### Step 3: Insertion of Position-Space Identity

To evaluate $`\langle \mathbf{x} | e^{-\frac{\delta\tau}{\hbar}\hat{T}} | \Psi(\tau) \rangle`$, insert the resolution of identity in position space:

```math
\hat{I} = \int d\mathbf{y} \, | \mathbf{y} \rangle \langle \mathbf{y} |
```

This gives:

```math
\langle \mathbf{x} | e^{-\frac{\delta\tau}{\hbar}\hat{T}} | \Psi(\tau) \rangle = \int d\mathbf{y} \, \langle \mathbf{x} | e^{-\frac{\delta\tau}{\hbar}\hat{T}} | \mathbf{y} \rangle \langle \mathbf{y} | \Psi(\tau) \rangle
```

Recognizing $`\langle \mathbf{y} | \Psi(\tau) \rangle = \Psi(\mathbf{y}, \tau)`$:

```math
\langle \mathbf{x} | e^{-\frac{\delta\tau}{\hbar}\hat{T}} | \Psi(\tau) \rangle = \int d\mathbf{y} \, \langle \mathbf{x} | e^{-\frac{\delta\tau}{\hbar}\hat{T}} | \mathbf{y} \rangle \Psi(\mathbf{y}, \tau)
```

---

#### Step 4: The Integral Evolution Equation

Combining Steps 2 and 3:

```math
\Psi(\mathbf{x}, \tau+\delta\tau) = e^{-\frac{\delta\tau}{\hbar}(V(\mathbf{x}) - E_T)} \int d\mathbf{y} \, \langle \mathbf{x} | e^{-\frac{\delta\tau}{\hbar}\hat{T}} | \mathbf{y} \rangle \Psi(\mathbf{y}, \tau) + O(\delta\tau^2)
```

The matrix element $`K(\mathbf{x}, \mathbf{y}) \equiv \langle \mathbf{x} | e^{-\frac{\delta\tau}{\hbar}\hat{T}} | \mathbf{y} \rangle`$ is the **kinetic propagator** (free-particle Green's function), which we evaluate in the next section.

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
K(\mathbf{x}, \mathbf{y}) = \left( \frac{m}{2\pi\hbar\delta\tau} \right)^{3/2} \exp\left( -\frac{m |\mathbf{x}-\mathbf{y}|^2}{2\hbar\delta\tau} \right) \equiv G_0(\mathbf{x}-\mathbf{y}, \delta\tau)
```

> **Notation:** This is exactly the Green's function $`G_0`$ derived in [Proof of the Gaussian Green's Function](#proof-of-the-gaussian-greens-function-free-diffusion), evaluated at the discrete time step $`\delta\tau`$. In subsequent sections, we use $`G_0`$ consistently to denote the free diffusion kernel.

Thus, the integral update equation is:

```math
\Psi(\mathbf{x}, \tau+\delta\tau) = \int d\mathbf{y} \underbrace{e^{-\frac{\delta\tau}{\hbar}(V(\mathbf{x}) - E_T)}}_{\text{Branching Weight } W} \underbrace{G_0(\mathbf{x}-\mathbf{y}, \delta\tau)}_{\text{Diffusion Probability}} \Psi(\mathbf{y}, \tau) + O(\delta\tau^2)
```

---

### Stationarity and Eigenfunction Correspondence

We show that stationary solutions of the integral evolution equation are eigenfunctions of the Hamiltonian.

---

#### The Stationarity Condition

Assume the distribution reaches a stationary state where it no longer changes with imaginary time:

```math
\Psi(\mathbf{x}, \tau+\delta\tau) = \Psi(\mathbf{x}, \tau) \equiv \Phi(\mathbf{x})
```

Substituting into the integral evolution equation:

```math
\Phi(\mathbf{x}) = \int d\mathbf{y} \, e^{-\frac{\delta\tau}{\hbar}(V(\mathbf{x}) - E_T)} G_0(\mathbf{x}-\mathbf{y}, \delta\tau) \, \Phi(\mathbf{y}) + O(\delta\tau^2)
```

---

#### Step 1: Expand the Branching Factor

Taylor expand the exponential to first order in $`\delta\tau`$:

```math
e^{-\frac{\delta\tau}{\hbar}(V(\mathbf{x}) - E_T)} = 1 - \frac{\delta\tau}{\hbar}(V(\mathbf{x}) - E_T) + O(\delta\tau^2)
```

Substituting:

```math
\Phi(\mathbf{x}) = \left( 1 - \frac{\delta\tau}{\hbar}(V(\mathbf{x}) - E_T) \right) \int d\mathbf{y} \, G_0(\mathbf{x}-\mathbf{y}, \delta\tau) \, \Phi(\mathbf{y}) + O(\delta\tau^2)
```

---

#### Step 2: Evaluate the Diffusion Integral

The integral $`\int d\mathbf{y} \, G_0(\mathbf{x}-\mathbf{y}, \delta\tau) \, \Phi(\mathbf{y})`$ represents the effect of diffusion on $`\Phi`$.

**Key Identity:** For a normalized Gaussian kernel with variance $`\sigma^2 = \frac{\hbar\delta\tau}{m}`$:

```math
\int d\mathbf{y} \, G_0(\mathbf{x}-\mathbf{y}, \delta\tau) \, \Phi(\mathbf{y}) = \Phi(\mathbf{x}) + \frac{\sigma^2}{2} \nabla^2 \Phi(\mathbf{x}) + O(\sigma^4)
```

**Derivation:** Expand $`\Phi(\mathbf{y})`$ in a Taylor series around $`\mathbf{x}`$:

```math
\Phi(\mathbf{y}) = \Phi(\mathbf{x}) + (\mathbf{y} - \mathbf{x}) \cdot \nabla\Phi(\mathbf{x}) + \frac{1}{2} \sum_{ij} (y_i - x_i)(y_j - x_j) \frac{\partial^2 \Phi}{\partial x_i \partial x_j} + \cdots
```

Since $`G_0`$ is a symmetric Gaussian centered at $`\mathbf{x}`$:
- $`\int d\mathbf{y} \, G_0(\mathbf{x}-\mathbf{y}) = 1`$ (normalization)
- $`\int d\mathbf{y} \, G_0(\mathbf{x}-\mathbf{y}) (y_i - x_i) = 0`$ (odd moments vanish)
- $`\int d\mathbf{y} \, G_0(\mathbf{x}-\mathbf{y}) (y_i - x_i)(y_j - x_j) = \sigma^2 \delta_{ij}`$ (variance)

Substituting $`\sigma^2 = \frac{\hbar\delta\tau}{m}`$:

```math
\int d\mathbf{y} \, G_0(\mathbf{x}-\mathbf{y}, \delta\tau) \, \Phi(\mathbf{y}) = \Phi(\mathbf{x}) + \frac{\hbar\delta\tau}{2m} \nabla^2 \Phi(\mathbf{x}) + O(\delta\tau^2)
```

---

#### Step 3: Combine and Expand

Substituting the diffusion result:

```math
\Phi(\mathbf{x}) = \left( 1 - \frac{\delta\tau}{\hbar}(V(\mathbf{x}) - E_T) \right) \left( \Phi(\mathbf{x}) + \frac{\hbar\delta\tau}{2m} \nabla^2 \Phi(\mathbf{x}) \right) + O(\delta\tau^2)
```

Expanding the product and keeping only terms up to $`O(\delta\tau)`$:

```math
\Phi(\mathbf{x}) = \Phi(\mathbf{x}) + \frac{\hbar\delta\tau}{2m} \nabla^2 \Phi(\mathbf{x}) - \frac{\delta\tau}{\hbar}(V(\mathbf{x}) - E_T) \Phi(\mathbf{x}) + O(\delta\tau^2)
```

---

#### Step 4: Cancel and Rearrange

Subtract $`\Phi(\mathbf{x})`$ from both sides:

```math
0 = \frac{\hbar\delta\tau}{2m} \nabla^2 \Phi(\mathbf{x}) - \frac{\delta\tau}{\hbar}(V(\mathbf{x}) - E_T) \Phi(\mathbf{x}) + O(\delta\tau^2)
```

Divide by $`\delta\tau`$ and take the limit $`\delta\tau \to 0`$:

```math
0 = \frac{\hbar}{2m} \nabla^2 \Phi(\mathbf{x}) - \frac{1}{\hbar}(V(\mathbf{x}) - E_T) \Phi(\mathbf{x})
```

Multiply through by $`-\hbar`$:

```math
0 = -\frac{\hbar^2}{2m} \nabla^2 \Phi(\mathbf{x}) + (V(\mathbf{x}) - E_T) \Phi(\mathbf{x})
```

Rearranging:

```math
-\frac{\hbar^2}{2m} \nabla^2 \Phi(\mathbf{x}) + V(\mathbf{x}) \Phi(\mathbf{x}) = E_T \Phi(\mathbf{x})
```

---

#### Conclusion

This is exactly the time-independent Schrödinger equation:

```math
\hat{H} \Phi(\mathbf{x}) = E_T \Phi(\mathbf{x})
```

**Result:** The stationary distribution $`\Phi(\mathbf{x})`$ is an eigenfunction of the Hamiltonian with eigenvalue $`E_T`$. This proves that the integral formulation is equivalent to the original quantum mechanical problem.

---

### Convergence to the Ground State

We prove that imaginary-time evolution projects any initial state onto the ground state.

---

#### Spectral Decomposition

The Hamiltonian $`\hat{H}`$ has a complete set of orthonormal eigenstates $`\{ | \phi_n \rangle \}`$ with eigenvalues $`E_n`$:

```math
\hat{H} | \phi_n \rangle = E_n | \phi_n \rangle, \quad \langle \phi_m | \phi_n \rangle = \delta_{mn}
```

Order the eigenvalues: $`E_0 < E_1 \leq E_2 \leq \cdots`$ (assuming a non-degenerate ground state).

Any initial state can be expanded in this basis:

```math
| \Psi(0) \rangle = \sum_{n=0}^{\infty} c_n | \phi_n \rangle, \quad c_n = \langle \phi_n | \Psi(0) \rangle
```

---

#### Time Evolution in the Eigenstate Basis

The formal solution to the imaginary-time Schrödinger equation is:

```math
| \Psi(\tau) \rangle = e^{-(\hat{H} - E_T)\tau/\hbar} | \Psi(0) \rangle
```

Applying this to the eigenstate expansion:

```math
| \Psi(\tau) \rangle = \sum_{n=0}^{\infty} c_n e^{-(E_n - E_T)\tau/\hbar} | \phi_n \rangle
```

Each eigenstate evolves with its own exponential factor $`e^{-(E_n - E_T)\tau/\hbar}`$.

---

#### Asymptotic Behavior

Define the energy gaps $`\Delta_n = E_n - E_0 \geq 0`$. Then:

```math
| \Psi(\tau) \rangle = e^{-(E_0 - E_T)\tau/\hbar} \sum_{n=0}^{\infty} c_n e^{-\Delta_n \tau/\hbar} | \phi_n \rangle
```

```math
= e^{-(E_0 - E_T)\tau/\hbar} \left[ c_0 | \phi_0 \rangle + \sum_{n=1}^{\infty} c_n e^{-\Delta_n \tau/\hbar} | \phi_n \rangle \right]
```

Since $`\Delta_n > 0`$ for all $`n \geq 1`$, the excited state contributions decay exponentially:

```math
e^{-\Delta_n \tau/\hbar} \xrightarrow{\tau \to \infty} 0 \quad \text{for } n \geq 1
```

**Result:** As $`\tau \to \infty`$:

```math
| \Psi(\tau) \rangle \xrightarrow{\tau \to \infty} c_0 e^{-(E_0 - E_T)\tau/\hbar} | \phi_0 \rangle
```

---

#### Convergence Rate

The rate of convergence is determined by the **spectral gap** $`\Delta_1 = E_1 - E_0`$:

```math
\left| \frac{\langle \phi_1 | \Psi(\tau) \rangle}{\langle \phi_0 | \Psi(\tau) \rangle} \right| = \left| \frac{c_1}{c_0} \right| e^{-\Delta_1 \tau/\hbar}
```

The characteristic **equilibration time** is:

```math
\tau_{\text{eq}} \sim \frac{\hbar}{\Delta_1}
```

Systems with small spectral gaps (near-degeneracy) require longer equilibration.

---

#### Condition for Convergence

**Requirement:** The initial state must have non-zero overlap with the ground state:

```math
c_0 = \langle \phi_0 | \Psi(0) \rangle \neq 0
```

If $`c_0 = 0`$, the system converges to the lowest-energy state that has non-zero overlap.

**Practical Note:** For the Bosonic ground state, which is nodeless and strictly positive, almost any reasonable initial guess satisfies $`c_0 \neq 0`$.

---

### Time Derivative of Number of Particles

We derive how the total "particle number" (norm of $`\Psi`$) evolves after convergence to the ground state.

---

#### Definition of Particle Number

Define the total particle number as:

```math
N(\tau) = \int d\mathbf{x} \, \Psi(\mathbf{x}, \tau) = \langle \mathbf{1} | \Psi(\tau) \rangle
```

where $`| \mathbf{1} \rangle`$ represents integration over all space.

---

#### Evolution After Convergence

After equilibration, $`| \Psi(\tau) \rangle \approx c_0 e^{-(E_0 - E_T)\tau/\hbar} | \phi_0 \rangle`$. Therefore:

```math
N(\tau) \approx c_0 e^{-(E_0 - E_T)\tau/\hbar} \langle \mathbf{1} | \phi_0 \rangle
```

Taking the time derivative:

```math
\frac{dN}{d\tau} = -\frac{1}{\hbar}(E_0 - E_T) \cdot c_0 e^{-(E_0 - E_T)\tau/\hbar} \langle \mathbf{1} | \phi_0 \rangle = -\frac{1}{\hbar}(E_0 - E_T) N(\tau)
```

---

#### Interpretation

| **Condition** | **Population Behavior** |
|:--------------|:------------------------|
| $`E_T < E_0`$ | $`\frac{dN}{d\tau} < 0`$ → Population decays exponentially |
| $`E_T = E_0`$ | $`\frac{dN}{d\tau} = 0`$ → Population remains constant |
| $`E_T > E_0`$ | $`\frac{dN}{d\tau} > 0`$ → Population grows exponentially |

This motivates the need for **population control**: adjusting $`E_T`$ dynamically to maintain a stable population.

---

### Population Control Bias and Feedback Law

Without population control, the walker population either explodes or collapses. We introduce a feedback mechanism to stabilize $`N(\tau)`$.

---

#### The Population Control Problem

From the previous section:

```math
\frac{dN}{d\tau} = -\frac{1}{\hbar}(E_0 - E_T) N(\tau)
```

Since $`E_0`$ is unknown (it's what we're trying to compute!), we cannot simply set $`E_T = E_0`$.

---

#### Feedback Control Law

We introduce a feedback mechanism that adjusts $`E_T`$ based on the observed population:

```math
E_T(\tau) = \bar{E}(\tau) - \alpha \ln\left( \frac{N(\tau)}{N_0} \right)
```

where:
- $`\bar{E}(\tau)`$ is an estimator for the ground state energy (e.g., average potential or local energy)
- $`N_0`$ is the target population size
- $`\alpha > 0`$ is the feedback strength parameter (units of energy)

**Intuition:**
- If $`N > N_0`$: The $`\ln`$ term is positive, so $`E_T`$ decreases, increasing $`(E_0 - E_T)`$, causing population to decrease
- If $`N < N_0`$: The $`\ln`$ term is negative, so $`E_T`$ increases, decreasing $`(E_0 - E_T)`$, causing population to increase

---

#### Local Stability Analysis

**Setup:** Near equilibrium, let $`N(\tau) = N_0(1 + \delta(\tau))`$ with $`|\delta| \ll 1`$.

**Assumptions:**
1. The energy estimator is accurate: $`\bar{E} \approx E_0`$
2. Small population deviation: $`\ln(1 + \delta) \approx \delta`$

**Step 1:** Substitute the control law into the population equation:

```math
E_0 - E_T = E_0 - \left( \bar{E} - \alpha \ln\left(\frac{N}{N_0}\right) \right) \approx E_0 - E_0 + \alpha \delta = \alpha \delta
```

**Step 2:** Substitute into the population dynamics:

```math
\frac{dN}{d\tau} = -\frac{1}{\hbar}(E_0 - E_T) N \approx -\frac{\alpha}{\hbar} \delta \cdot N_0
```

**Step 3:** Since $`N = N_0(1 + \delta)`$, we have $`\frac{dN}{d\tau} = N_0 \frac{d\delta}{d\tau}`$:

```math
N_0 \frac{d\delta}{d\tau} = -\frac{\alpha}{\hbar} \delta N_0
```

```math
\frac{d\delta}{d\tau} = -\frac{\alpha}{\hbar} \delta
```

**Solution:** Exponential decay of population fluctuations:

```math
\delta(\tau) = \delta(0) e^{-\alpha\tau/\hbar}
```

The population stabilizes with time constant $`\tau_{\text{stab}} = \hbar/\alpha`$.

---

#### Choice of Feedback Parameter $`\alpha`$

The feedback strength $`\alpha`$ controls the trade-off between stability and bias:

| **$`\alpha`$ Value** | **Effect** |
|:---------------------|:-----------|
| Large $`\alpha`$ | Fast stabilization, but larger bias in energy estimate |
| Small $`\alpha`$ | Slower stabilization, smaller bias, larger population fluctuations |

A common choice is $`\alpha \sim 1`$ Hartree (in atomic units).

---

#### Population Control Bias

The feedback mechanism introduces a systematic bias in the energy estimate.

**Source of Bias:** The trial energy $`E_T`$ depends on the fluctuating population $`N(\tau)`$. Since $`E_T`$ and $`N`$ are correlated, the average $`\langle E_T \rangle`$ is not equal to the average over independent samples.

**Bias Formula:** The DMC energy estimator has a bias of order $`1/N`$:

```math
E_{\text{DMC}} = E_0 + O\left( \frac{1}{N} \right)
```

**Mitigation Strategies:**
1. **Large populations:** Use $`N \gg 1`$ to reduce the bias
2. **Extrapolation:** Run simulations at multiple population sizes and extrapolate to $`N \to \infty`$
3. **Correlated sampling:** Use variance reduction techniques to reduce the effective bias

---

## Pure Diffusion Monte Carlo

### Theory

Pure DMC directly simulates the imaginary-time Schrödinger equation by interpreting the wavefunction $`\Psi(\mathbf{x}, \tau)`$ as a probability density represented by a population of random walkers.

The update equation derived from Trotter-Suzuki decomposition is:

```math
\Psi(\mathbf{x}, \tau+\delta\tau) = \int d\mathbf{y} \, W(\mathbf{x}) \, G_0(\mathbf{x}-\mathbf{y}, \delta\tau) \, \Psi(\mathbf{y}, \tau)
```

where:

| **Component** | **Expression** | **Role** |
|:--------------|:---------------|:---------|
| Diffusion kernel | $`G_0(\mathbf{x}-\mathbf{y}) = \left(\frac{m}{2\pi\hbar\delta\tau}\right)^{3/2} e^{-\frac{m\|\mathbf{x}-\mathbf{y}\|^2}{2\hbar\delta\tau}}`$ | Gaussian random walk |
| Branching weight | $`W(\mathbf{x}) = e^{-\frac{\delta\tau}{\hbar}(V(\mathbf{x}) - E_T)}`$ | Birth/death of walkers |

**Key Properties:**
- The Green's function is **symmetric**: $`G_0(\mathbf{x}-\mathbf{y}) = G_0(\mathbf{y}-\mathbf{x})`$
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
G_0(\boldsymbol{\xi}) = \left( \frac{m}{2\pi\hbar\delta\tau} \right)^{3/2} \exp\left( -\frac{m |\boldsymbol{\xi}|^2}{2\hbar\delta\tau} \right)
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
G_0(\mathbf{x}, \mathbf{y}) = G_0(|\mathbf{x} - \mathbf{y}|) \implies G(\mathbf{x} \leftarrow \mathbf{y}) = G(\mathbf{y} \leftarrow \mathbf{x})
```

This is because $`G_0`$ (the free diffusion kernel) is an isotropic Gaussian depending only on the distance $`|\mathbf{x} - \mathbf{y}|`$.

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

This is why we need the **Metropolis acceptance/rejection step** to restore detailed balance (see [Step 3: Metropolis Acceptance/Rejection](#step-3-metropolis-acceptance-rejection) in the Implementation section).

---

#### Summary: Pure DMC vs. Importance Sampled DMC

| **Feature** | **Pure DMC ($`\Psi`$)** | **Importance Sampled DMC ($`f = \Psi \Psi_T`$)** |
|:------------|:------------------------|:-------------------------------------------------|
| **Quantity Simulated** | Wavefunction $`\Psi(\mathbf{x})`$ | Mixed distribution $`f(\mathbf{x}) = \Psi \Psi_T`$ |
| **Propagator** | Isotropic Gaussian | Gaussian shifted by Drift |
| **Green's Function** | Symmetric: $`G_0(\mathbf{x}-\mathbf{y})`$ | **Asymmetric:** $`G(\mathbf{x} \leftarrow \mathbf{y}) \neq G(\mathbf{y} \leftarrow \mathbf{x})`$ |
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
