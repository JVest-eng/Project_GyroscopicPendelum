# Gyroscopic Pendelum Project
3-DOF Gyroscopic pendelum stabilization project analyzing the behavior of LQR control against Disturbance Rejection.
<img width="1254" height="1254" alt="image" src="https://github.com/user-attachments/assets/b752818b-a949-42a5-992e-869638612e12" />

## System Description
**System Description**: A rigid, hollow rod with mass $M$, length $L$, and diameter $d$ fixed with a ball joint pivot and a probe of mass $m$ can undergo precession, nutation, and tilt. The purpose of this project is to examine whether reaction wheels could stabilize the rod such that it does not decrease its gravity potential. 
	
The wheels each have a spin inertia of $J_s$ each and the assembly is axis symmetric about the rod axis.
	
**Generalized Coordinates**: $\mathbf{q} = [\phi, \theta, \gamma]$ where $\phi$ is precession angle, $\theta$ is nutation angle, and $\gamma$ is spin angle.
	
**Control inputs**: Wheel relative spinrates $\boldsymbol{\Omega} = [\Omega_i, \Omega_j, \Omega_{\gamma}]$ driven by motor torques $\boldsymbol{\tau} = [\tau_i, \tau_j, \tau_{\gamma}]$.

This project uses a mix of Euler-Lagrangian formulations and Newtonian derivation to.

The lagrangian is given by $\mathcal{L} = T - V$ where $\mathcal{L}$ is the system energy descriptor, $T$ is the kinetic energy, and $V$ is the potential energy. The Lagrangian by itself is very powerful for defining and modelling *conservative* systems. The d'Alembert-lagrangian modification allows the system to be generalized for systems that do non-restorative work.

## Formulating the Lagrangian

### 1. Reference Frames
There is a **Fixed Frame** attached to the base of the hollow rod at the ball pivot joint and a **Body Frame** attached to the end of the rod at the probe.

This fixed frame $\mathbf{\hat{F}_O}$ = $\{ \hat{e}_x,\hat{e}_y, \hat{e}_z \}$ has $\hat{e}_z$ pointing vertically upward against the gravitational potential.

The body frame $\mathbf{\hat{F}\_B}$ = $\{ \hat{e}\_i, \hat{e}\_j, \hat{e}\_\gamma \}$ has $\hat{e}_\gamma$ pointing along the axis of the rod outward with $\hat{e}_i$ and $\hat{e}_j$ being tangential to the rod's rotation.

The Euler rotation that transforms any fixed frame coordinate to the body coordinate is given as:

$$R_{OB} = R_z(\phi) R_x(\theta) R_z(\gamma)$$

Expressing fixed frame basis vectors in the body frame basis vectors can be done by taking the elements of the rotation tensor $R_{OB}$. The results are collected below:

$$\hat{e}\_x = (c\_\phi c_\gamma - s_\phi c_\theta s_\gamma)\hat{e}\_i - (c\_\phi s\_\gamma + s_\phi c_\theta c_\gamma) \hat{e}\_j + (s_\phi s_\theta)\hat{e}_\gamma$$

$$\hat{e}\_y = (s_\phi c_\gamma + c_\phi c_\theta s_\gamma)\hat{e}_i + (c_\phi c_\theta c_\gamma - s_\phi s_\gamma) \hat{e}\_j - (c_\phi s_\theta)\hat{e}_\gamma$$

$$\hat{e}\_z = (s\_\phi s_\gamma)\hat{e}\_i + (s_\theta c_\gamma )\hat{e}\_j + c\_\theta \hat{e}\_\gamma$$

Do note that $c_x = \cos(x)$ and $s_x = \sin(x)$ for shorthand notation in this derivation.

### 2. Angular Velocity in the Body Frame

The Euler angles contribute to the total angular velocity vector $\vec{\omega}$ as:

$$\vec{\omega} = \dot{\phi} \hat{e}_z + \dot{\theta} (\cos\phi \hat{e}_x + \sin\phi \hat{e}_y)+ \dot{\gamma} \hat{e}_\gamma$$

Substituting the results for the fixed frame coordinates into $\vec{\omega}$ yields a vector of three components:

$$\omega_i = \dot{\theta}c_\gamma+\dot\phi s_\theta s_\gamma$$
$$\omega_j = -\dot{\theta}s_\gamma+\dot\phi s_\theta c_\gamma$$
$$\omega_\gamma = \dot{\gamma}+\dot\phi c_\theta$$

There is an axisymmetric identity from $\omega_i$ and $\omega_j$ as follows:

$${\omega_i}^2 + {\omega_j}^2 = \dot{\theta}^2+\dot\phi^2 {s_\theta}^2$$

### 3. Kinetic Energy

There are two mass moment of inertias that need to be defined before beginning: $I_p$ and $I_a$. 

Mass moment $I_{p}$ represents the mass moment of inertia of the rod and probe assembly *perpendicular and tangent to the pivot*. 

Likewise, $I_{a}$ represents the mass moment of inertia along the rod axis. Translation using the parallel axis theorem can be used to derive the specific geometry's inertia away from its center of mass.

The inertial frame is not moving, therefore the only kinetic energy comes from the rotation of the assembly about the pivot.

**Importantly**, each reaction wheel has its own moment $J_s$. This distinction is important because there will be coupled terms strictly dependent on $J_s$ and lumping the moment in with $I_{p}$ and $I_{a}$ formulates a fully rigid body with no wheels able to spin at the end of the rod.

Instead, an effective inertia tensor can be used which replaces the summands of lesser inertias defined as such:

$$I_T = I_{p} + J_s \quad \text{(transverse)}, \qquad I_A = I_{a} + J_s \quad \text{(axial)}$$

The **kinetic energy of the assembly** is given as:

$$T_{body} = \frac{1}{2}I_{p}\left(\omega_i^2 + \omega_j^2\right) + \frac{1}{2}I_{a}\left(\omega_\gamma\right)^2$$

The **kinetic energy of the reaction wheel** -- each wheel rotates at absolute rate $\omega_k + \Omega_k$:
 
$$T_{wheels} = \frac{1}{2}J_s\left[(\omega_i+\Omega_i)^2 + (\omega_j+\Omega_j)^2 + (\omega_\gamma+\Omega_\gamma)^2\right]$$

Expanding $T_{wheels}$ and grouping like terms gives the result:

$$T_{wheels} = \frac{1}{2}J_s\left(\omega_i^2 + \omega_j^2 + \omega_\gamma^2 + \Omega_i^2 + \Omega_j^2 + \Omega_\gamma^2 + 2\omega_i\Omega_i + 2\omega_j\Omega_j + 2\omega_\gamma\Omega_\gamma\right)$$

Combining $T_{wheels}$ with $T_{body}$:

$$T = \frac{1}{2}\left(I_p + J_s\right)(\omega_i^2 + \omega_j^2) + \frac{1}{2}\left(I_a + J_s\right)\omega_\gamma^2+\frac{1}{2}J_s\sum_{k}\left(\Omega_k^2 + 2\omega_k\Omega_k\right)$$

Using the axisymmetric identity of $\omega_i$ and $\omega_j$:

$$T = \frac{1}{2}I_T(\dot{\theta}^2+\dot\phi^2 {s_\theta}^2) + \frac{1}{2}I_A\left({\dot{\gamma}^2+2\dot{\gamma}\dot{\phi} c_\theta}+{\dot\phi}^2 {c_\theta}^2\right)+\frac{1}{2}J_s\sum_{k}\left(\Omega_k^2 + 2\omega_k\Omega_k\right)$$
