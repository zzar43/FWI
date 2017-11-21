---
title: "Note of Adjoint Equation"
author: Da Li
output:
  html_document:
    toc: true
    theme: flatly
    toc_depth: 2
    number_sections: true
    # toc_float:
    #   collapsed: false
    #   smooth_scroll: false
---


# Tensor Notion

- $s$ denotes a scalar
- $\mathbf{u}$ denotes a vector
- $\mathbf{F}$ denotes a tensor.

## Tensor Product of Two Second-Order Tensors

- Given two second-order tensors,
$$\mathbf{U} = U_{ij}\mathbf{e}_i\mathbf{e}_j, \quad \mathbf{V} = V_{ij}\mathbf{e}_i\mathbf{e}_j.$$

- Tensor Product (Outer Product):
$$\mathbf{U}\otimes\mathbf{V} = U_{ij}V_{kl}\mathbf{e}_i\mathbf{e}_j\mathbf{e}_k\mathbf{e}_l.$$

- Singal Dot Product (Contraction Operator):
$$\mathbf{U}\cdot\mathbf{V} = \mathbf{U}\mathbf{V} = (U_{ij}\mathbf{e}_i\mathbf{e}_j)\cdot(V_{kl}\mathbf{e}_k\mathbf{e}_l)
= U_{ij}\mathbf{e}_i\delta_{ij}V_{kl}\mathbf{e}_l = U_{ij}V_{jl}\mathbf{e}_i\mathbf{e}_l.$$

- Double Dot Product (Scalar or Inner Product):
$$\mathbf{U}: \mathbf{V} = (U_{ij}\mathbf{e}_i\mathbf{e}_j):(V_{kl}\mathbf{e}_k\mathbf{e}_l) = (U_{ij}\mathbf{e}_i)\delta_{jk}\cdot(V_{kl}\mathbf{e}_l) = U_{ij}V_{kl}\delta_{jk}\delta_{il} = U_{ij}V_{ji}.$$

## Gradient of a Field
A function of the position vector $\mathbf{x}$ is called a field.
Following discussion is in Cartesian Frame.

- Gradient of a Scalar:
$$\nabla\phi = \mathbf{e}_j\frac{\partial\phi}{x_j} = \mathbf{e}_1\frac{\partial\phi}{x_1} + \mathbf{e}_2\frac{\partial\phi}{x_2} + \mathbf{e}_3\frac{\partial\phi}{x_3}.$$

- Gradient of a Vector:
$$\nabla\mathbf{u} = \left(\mathbf{e}_j\frac{\partial}{x_j}\right)(u_j\mathbf{e}_j) = \mathbf{e}_i\mathbf{e}_j\frac{\partial \mathbf{u}_j}{\partial x_i}.$$
In short $(\nabla\mathbf{u})_{ij} = \partial \mathbf{u}_j/ \partial x_i$, in some books defination is different (the transpose of a gradient).

- Transpose of a Gradient:
$$\nabla\mathbf{u}^T = \mathbf{e}_i\mathbf{e}_j\frac{\partial \mathbf{u}_i}{\partial x_j}.$$

- Divergence of a Vector:
$$\nabla\cdot\mathbf{u} = \left(\mathbf{e}_j\frac{\partial}{x_j}\right)\cdot(u_j\mathbf{e}_j) = \mathbf{e}_i\cdot\mathbf{e}_j\frac{\partial \mathbf{u}_j}{\partial x_i} = \delta_{ij}\frac{\partial \mathbf{u}_j}{\partial x_i}.$$
That is,
$$\nabla\cdot\mathbf{u} = \frac{\partial u_i}{x_i} = \frac{\partial u_1}{x_1} + \frac{\partial u_2}{x_2} + \frac{\partial u_3}{x_3}.$$

- Divergence of a Tensor:
$$\nabla\cdot\mathbf{S} = \left(\mathbf{e}_j\frac{\partial}{x_j}\right) \cdot(S_{ij}\mathbf{e}_i\mathbf{e}_j) = \mathbf{e}_j\frac{\partial S_{ij}}{\partial x_i}.$$

## Prove Equality
$$\mathbf{u}^\dagger \cdot [\nabla \cdot (\dot{\mathbf{C}}:\nabla \mathbf{u})] - \mathbf{u} \cdot [\nabla \cdot (\dot{\mathbf{C}}:\nabla \mathbf{u}^\dagger)] = \nabla\cdot(\mathbf{u}^\dagger\cdot\dot{\mathbf{C}} : \nabla \mathbf{u}) - \nabla\cdot(\mathbf{u}\cdot\dot{\mathbf{C}} : \nabla \mathbf{u}^\dagger).$$

### Left Part
For $\mathbf{u}^\dagger \cdot [\nabla \cdot (\dot{\mathbf{C}}:\nabla \mathbf{u})]$:
$$\begin{align*}
\dot{\mathbf{C}}:\nabla \mathbf{u} &= (\dot C_{ijkl}\mathbf{e}_i\mathbf{e}_j\mathbf{e}_k\mathbf{e}_l):
\left(\left(\mathbf{e}_m\frac{\partial}{\partial x_m}\right)(u_n \mathbf{e}_n)\right) \\
&= (\dot C_{ijkl}\mathbf{e}_i\mathbf{e}_j\mathbf{e}_k\mathbf{e}_l):
\left(\mathbf{e}_m\mathbf{e}_n\frac{\partial u_n}{\partial x_m} \right)\\
&= \mathbf{e}_i\mathbf{e}_j \dot C_{ijkl} \frac{\partial u_n}{\partial x_m}\delta_{lm}\delta_{kn} = \mathbf{e}_i\mathbf{e}_j \dot C_{ijkl} \frac{\partial u_k}{\partial x_l}.
\end{align*}$$
Then,
$$\begin{align*}
\nabla \cdot (\dot{\mathbf{C}}:\nabla \mathbf{u})  &= \left(\mathbf{e}_m\frac{\partial}{\partial x_m}\right) \cdot \left(\mathbf{e}_i\mathbf{e}_j \dot C_{ijkl} \frac{\partial u_k}{\partial x_l}\right)\\
&=\mathbf{e}_j\frac{\partial}{\partial x_i} \left(\dot C_{ijkl}\frac{\partial u_k}{\partial x_l}\right).
\end{align*}$$
So,
$$\begin{align*}
\mathbf{u}^\dagger \cdot [\nabla \cdot (\dot{\mathbf{C}}:\nabla \mathbf{u})] &= (\mathbf{e}_m u_m^\dagger)\cdot \left(\mathbf{e}_j\frac{\partial}{\partial x_i} \left(\dot C_{ijkl}\frac{\partial u_k}{\partial x_l}\right)\right)\\
&= u_j^\dagger \frac{\partial}{\partial x_i} \left(\dot C_{ijkl}\frac{\partial u_k}{\partial x_l}\right)\\
&= u_j^\dagger \frac{\partial \dot C_{ijkl}}{\partial x_i} \frac{\partial u_k}{\partial x_l} + u_j^\dagger \dot C_{ijkl} \frac{\partial^2 u_k}{\partial x_i \partial x_l}.
\end{align*}$$

By the same way,
$$\mathbf{u} \cdot [\nabla \cdot (\dot{\mathbf{C}}:\nabla \mathbf{u}^\dagger)] = u_j \frac{\partial \dot C_{ijkl}}{\partial x_i} \frac{\partial u_k^\dagger}{\partial x_l} + u_j \dot C_{ijkl} \frac{\partial^2 u_k^\dagger}{\partial x_i \partial x_l}.$$

### Right Part
For the $\nabla\cdot(\mathbf{u}^\dagger\cdot\dot{\mathbf{C}} : \nabla \mathbf{u})$,
$$\mathbf{u}^\dagger\cdot\dot{\mathbf{C}} = (\mathbf{e}_m u_m^\dagger) \cdot (\dot C_{ijkl} \mathbf{e}_i\mathbf{e}_j\mathbf{e}_k\mathbf{e}_l) = \mathbf{e}_j\mathbf{e}_k\mathbf{e}_l u_i^\dagger \dot C_{ijkl}$$

Then,
$$\begin{align*}
\mathbf{u}^\dagger\cdot\dot{\mathbf{C}} : \nabla \mathbf{u} &= (\mathbf{e}_j\mathbf{e}_k\mathbf{e}_l u_i^\dagger \dot C_{ijkl}) : \left(\left(\mathbf{e}_m\frac{\partial}{\partial x_m}\right)(u_n\mathbf{e}_n)\right)\\
&= (\mathbf{e}_j\mathbf{e}_k\mathbf{e}_l u_i^\dagger \dot C_{ijkl}) : \left(\mathbf{e}_m\mathbf{e}_n\frac{\partial u_n}{\partial x_m}\right)\\
&= (\mathbf{e}_j\mathbf{e}_k u_i^\dagger \dot C_{ijkl})\delta_{lm}\cdot \left(\mathbf{e}_n\frac{\partial u_n}{\partial x_m}\right)\\
&= \mathbf{e}_j u_i^\dagger \dot C_{ijkl} \frac{\partial u_n}{\partial x_m} \delta_{lm}\delta_{kn}\\
&= \mathbf{e}_j u_i^\dagger \dot C_{ijkl} \frac{\partial u_k}{\partial x_l}.
\end{align*}$$
Then,
$$\begin{align*}
\nabla\cdot(\mathbf{u}^\dagger\cdot\dot{\mathbf{C}} : \nabla \mathbf{u}) &= \left(\mathbf{e}_m\frac{\partial}{\partial x_m}\right) \cdot \left(\mathbf{e}_j u_i^\dagger \dot C_{ijkl} \frac{\partial u_k}{\partial x_l}\right)\\
&= \frac{\partial}{\partial x_j}\left(u_i^\dagger \dot C_{ijkl} \frac{\partial u_k}{\partial x_l}\right)\\
&= \frac{\partial u_i^\dagger}{\partial x_j}\dot C_{ijkl}\frac{\partial u_k}{\partial x_l} + u_i^\dagger\frac{\partial \dot C_{ijkl}}{\partial x_j}\frac{\partial u_k}{\partial x_l} + u_i^\dagger \dot C_{ijkl}\frac{\partial^2 u_k}{\partial x_l\partial x_j}.
\end{align*}$$

By the same way,
$$\nabla\cdot(\mathbf{u}\cdot\dot{\mathbf{C}} : \nabla \mathbf{u}^\dagger)
= \frac{\partial u_i}{\partial x_j}\dot C_{ijkl}\frac{\partial u_k^\dagger}{\partial x_l} + u_i\frac{\partial \dot C_{ijkl}}{\partial x_j}\frac{\partial u_k^\dagger}{\partial x_l} + u_i \dot C_{ijkl}\frac{\partial^2 u_k^\dagger}{\partial x_l\partial x_j}$$

### Final

So,

$$\begin{align*}
\mathrm{LEFT} &=
u_j^\dagger \frac{\partial \dot C_{ijkl}}{\partial x_i} \frac{\partial u_k}{\partial x_l} + u_j^\dagger \dot C_{ijkl} \frac{\partial^2 u_k}{\partial x_i \partial x_l} - \left(u_j \frac{\partial \dot C_{ijkl}}{\partial x_i} \frac{\partial u_k^\dagger}{\partial x_l} + u_j \dot C_{ijkl} \frac{\partial^2 u_k^\dagger}{\partial x_i \partial x_l} \right).
\end{align*}$$

And,
$$\begin{align*}
\mathrm{RIGHT} &= \frac{\partial u_i^\dagger}{\partial x_j}\dot C_{ijkl}\frac{\partial u_k}{\partial x_l} + u_i^\dagger\frac{\partial \dot C_{ijkl}}{\partial x_j}\frac{\partial u_k}{\partial x_l} + u_i^\dagger \dot C_{ijkl}\frac{\partial^2 u_k}{\partial x_l\partial x_j}\\
&- \left(\frac{\partial u_i}{\partial x_j}\dot C_{ijkl}\frac{\partial u_k^\dagger}{\partial x_l} + u_i\frac{\partial \dot C_{ijkl}}{\partial x_j}\frac{\partial u_k^\dagger}{\partial x_l} + u_i \dot C_{ijkl}\frac{\partial^2 u_k^\dagger}{\partial x_l\partial x_j} \right).
\end{align*}$$

Since $\mathbf{C}$ is symmetric, that is $C_{ijkl} = C_{klij} = C_{jikl}$.
Then,
$$\mathrm{LEFT} = \mathrm{RIGHT}.$$

# General Discuss

Define an object function (misfit function) by,
$$\chi(\mathbf{m}) = \int_G\int_T \frac{1}{2} [\mathbf{u}^0(\mathbf{x},t) - \mathbf{u}(\mathbf{m};\mathbf{x},t)]^2 \delta(\mathbf{x} - \mathbf{x^r})\ \mathrm{d}t \ \mathrm{d}^3\mathbf{x}.$$

Let $\chi_1 = \frac{1}{2} [\mathbf{u}^0(\mathbf{x},t) - \mathbf{u}(\mathbf{m};\mathbf{x},t)]^2 \delta(\mathbf{x} - \mathbf{x^r})$, then,
$$\chi(\mathbf{m}) = \int_G\int_T \chi_1(\mathbf{u}(\mathbf{m};\mathbf{x},t))\ \mathrm{d}t \ \mathrm{d}^3\mathbf{x} = <\chi_1(\mathbf{m})>,$$
with $<\cdot>$ shows the integral over $G\times T$.

Our target is to compute $\nabla_m\chi\delta\mathbf{m}$.
$$\begin{align*}
\nabla_m\chi\delta\mathbf{m} = \nabla_u\chi\delta\mathbf{u} = \nabla_u\chi(\nabla_m\mathbf{u}\delta\mathbf{m}).
\end{align*}$$

And also,
$$\begin{align*}
\nabla_m\chi\delta\mathbf{m} &= (\delta\mathbf{m} \cdot \nabla_m)\chi \\
&= (\sum_i \delta m_i \partial_{m_i})\chi \\
&= (\sum_i \delta m_i \partial_{m_i}) \int_G\int_T \chi_1(\mathbf{u}(\mathbf{m};\mathbf{x},t))\ \mathrm{d}t \ \mathrm{d}^3\mathbf{x} \\
&= \int_G\int_T (\sum_i \delta m_i \partial_{m_i}) \chi_1(\mathbf{u}(\mathbf{m};\mathbf{x},t))\ \mathrm{d}t \ \mathrm{d}^3\mathbf{x} \\
&= \int_G\int_T \nabla_m\chi_1\delta\mathbf{m} \ \mathrm{d}t \ \mathrm{d}^3\mathbf{x} \\
&= \int_G\int_T \nabla_u\chi_1\delta\mathbf{u} \ \mathrm{d}t \ \mathrm{d}^3\mathbf{x} \\
&= \int_G\int_T (\delta\mathbf{u}\cdot \nabla_u)\chi_1 \ \mathrm{d}t \ \mathrm{d}^3\mathbf{x} \\
&= \int_G\int_T \delta\mathbf{u} \cdot \nabla_u\chi_1 \ \mathrm{d}t \ \mathrm{d}^3\mathbf{x} = <\delta\mathbf{u} \cdot \nabla_u\chi_1>.
\end{align*}$$

We can not compute $\nabla_m\chi\delta\mathbf{m}$ directly, so we need the adjoint method.

Define the wave equation be $\mathbf{L}(\mathbf{u}, \mathbf{m}) = \mathbf{f}$, we do the derivative with respect to $\mathbf{m}$,
$$\nabla_m\mathbf{L}\delta\mathbf{m} + \nabla_u\mathbf{L}\delta\mathbf{u} = \mathbf{0},$$

With $\delta\mathbf{u} = \nabla_m\mathbf{u}\delta\mathbf{m}$. Suppose there exist a $\mathbf{u}^\dagger$, then,
$$<\mathbf{u}^\dagger \cdot (\nabla_m\mathbf{L}\delta\mathbf{m} + \nabla_u\mathbf{L}\delta\mathbf{u})> = 0.$$

That is,
$$<\mathbf{u}^\dagger \cdot \nabla_m\mathbf{L}\delta\mathbf{m}> + <\mathbf{u}^\dagger \cdot \nabla_u\mathbf{L}\delta\mathbf{u}> = 0.$$

By introduce adjoint operator $\nabla_u\mathbf{L}^\dagger$, we have,
$$<\mathbf{u}^\dagger \cdot \nabla_m\mathbf{L}\delta\mathbf{m}> + <\delta\mathbf{u} \cdot \nabla_u\mathbf{L}^\dagger\mathbf{u}^\dagger> = 0.$$

And since,
$$\begin{align*}
\nabla_m\chi\delta\mathbf{m} &= <\delta\mathbf{u} \cdot \nabla_u\chi_1> \\
&= <\delta\mathbf{u} \cdot \nabla_u\chi_1> + <\mathbf{u}^\dagger \cdot \nabla_m\mathbf{L}\delta\mathbf{m}> + <\delta\mathbf{u} \cdot \nabla_u\mathbf{L}^\dagger\mathbf{u}^\dagger> \\
& = <\delta\mathbf{u} \cdot (\nabla_u\chi_1 + \nabla_u\mathbf{L}^\dagger\mathbf{u}^\dagger)> + <\mathbf{u}^\dagger \cdot \nabla_m\mathbf{L}\delta\mathbf{m}>.
\end{align*}$$

Then we have the adjoint equation,
$$\nabla_u\mathbf{L}^\dagger\mathbf{u}^\dagger = \nabla_u\chi_1.$$

By solving the adjoint equation, we can have the wave field $\mathbf{u}^\dagger$, then,
$$\nabla_m\chi\delta\mathbf{m} = <\mathbf{u}^\dagger \cdot \nabla_m\mathbf{L}\delta\mathbf{m}>.$$

The following problems are:

1. Determine what is $\nabla_u\mathbf{L}^\dagger$. Since $\mathbf{L}$ is linear in $\mathbf{u}$, and with suitable boundary condition, the differential operator $\partial^2/\partial x^2$ is self-adjoint operator. We have $\mathbf{L}^\dagger = \mathbf{L}$.
2. Determine adjoint source $\nabla_u\chi_1$.
3. Calculate the integral $\nabla_m\chi\delta\mathbf{m} = <\mathbf{u}^\dagger \cdot \nabla_m\mathbf{L}\delta\mathbf{m}>$.


# Adjoint Equation of Elastic Wave Equation

In this part, we prove the form of $\mathbf{L}^\dagger$.

For the elastic wave equation:
$$\begin{align*}
\rho(\mathbf{x})\ddot{\mathbf{u}}(\mathbf{x},t) - \nabla\cdot\boldsymbol{\sigma}(\mathbf{x},t) &= \mathbf{f}(\mathbf{x},t),
\quad \mathbf{x}\in G\subset\mathbb{R}^3, \quad t\in[t_0,\infty]\subset\mathbb{R},\\
\boldsymbol{\sigma}(\mathbf{x},t) &=  \int_{t_0}^\infty
\dot{\mathbf{C}}(\mathbf{x},t-\tau):
\nabla\mathbf{u}(\mathbf{x},\tau)\ \mathrm{d}\tau.
\end{align*}$$

Thus we have,
$$\mathbf{L}(\mathbf{u}, \rho, \mathbf{C}) = \rho(\mathbf{x})\ddot{\mathbf{u}}(\mathbf{x},t) - \nabla\cdot \int_{t_0}^\infty
\dot{\mathbf{C}}(\mathbf{x},t-\tau):
\nabla\mathbf{u}(\mathbf{x},\tau)\ \mathrm{d}\tau,$$
with the initial and boundary conditions,
$$\mathbf{u}|_{t\leq t_0} = \dot{\mathbf{u}}|_{t\leq t_0} = \mathbf{0},\quad \mathbf{n}\cdot\boldsymbol{\sigma} = 0.$$

For the adjoint operator $\mathbf{L}^\dagger$, we have defined,
$$\begin{align*}
<\mathbf{u}\cdot\mathbf{L}^\dagger\mathbf{u}^\dagger> &= <\mathbf{u}^\dagger\cdot \mathbf{L}\mathbf{u}> \\
&= \int_G\int_T \rho\mathbf{u}^\dagger\cdot\ddot{\mathbf{u}}\ \mathrm{d}t\ \mathbf{d}^3\mathbf{x} - \int_G\int_T \mathbf{u}^\dagger\cdot \left[\nabla\cdot \int_{t_0}^\infty
\dot{\mathbf{C}}(t-\tau):
\nabla\mathbf{u}(\tau)\ \mathrm{d}\tau\right]\ \mathrm{d}t\ \mathbf{d}^3\mathbf{x}
\end{align*}$$

For the first part, we use integration by parts,
$$\begin{align*}
<\mathbf{u}^\dagger \cdot \rho\ddot{\mathbf{u}}> &= \int_G\int_T \rho\mathbf{u}^\dagger\cdot\ddot{\mathbf{u}}\ \mathrm{d}t\ \mathbf{d}^3\mathbf{x} \\
&= \int_G\int \rho\mathbf{u}^\dagger\cdot\ \mathrm{d}\dot{\mathbf{u}}\ \mathbf{d}^3\mathbf{x} \\
&= \int_G \rho \left[\mathbf{u}^\dagger\cdot\dot{\mathbf{u}}|_{t_0}^{t_1} - \int_T \dot{\mathbf{u}} \cdot \dot{\mathbf{u}}^\dagger \ \mathrm{d}t \right]\ \mathbf{d}^3\mathbf{x} \\
&= \int_G \rho \left[\mathbf{u}^\dagger\cdot\dot{\mathbf{u}}|_{t_0}^{t_1} - \int \dot{\mathbf{u}}^\dagger \cdot\ \mathrm{d}\mathbf{u} \right]\ \mathbf{d}^3\mathbf{x} \\
&= \int_G \rho \left[\mathbf{u}^\dagger\cdot\dot{\mathbf{u}}|_{t_0}^{t_1} - \mathbf{u}\cdot\dot{\mathbf{u}}^\dagger|_{t_0}^{t_1} + \int_T \ddot{\mathbf{u}}^\dagger \cdot \mathbf{u}\ \mathrm{d}t \right]\ \mathbf{d}^3\mathbf{x}.
\end{align*}$$
We introduce the terminal conditions for adjoint equation, $\mathbf{u}^\dagger|_{t\geq t_0} = \dot{\mathbf{u}}^\dagger|_{t\geq t_0} = \mathbf{0}$. We have,
$$<\mathbf{u}^\dagger \cdot \rho\ddot{\mathbf{u}}> = \int_G\int_T\rho\ddot{\mathbf{u}}^\dagger \cdot \mathbf{u}\ \mathrm{d}t\ \mathrm{d}^3\mathbf{x} = <\mathbf{u}\cdot\rho\ddot{\mathbf{u}}^\dagger>.$$

Let $\gamma:= <\mathbf{u}^\dagger \cdot (\nabla\cdot\mathbf{\sigma})>$,
$$\begin{align*}\gamma &= \int_T\int_G \mathbf{u}^\dagger(t) \cdot \left[ \nabla\cdot \int_{t_0}^{t}\dot{\mathbf{C}}(t-\tau):\nabla\mathbf{u}(\tau)\ \mathrm{d}\tau\right]\ \mathrm{d}^3\mathbf{x}\ \mathrm{d}t \\
&= \int_T\int_G\int_{t_0}^{t} \mathbf{u}^\dagger(t) \cdot (\nabla\cdot(\dot{\mathbf{C}}(t-\tau):\nabla\mathbf{u}(\tau)))\ \mathrm{d}\tau\ \mathrm{d}^3\mathbf{x}\ \mathrm{d}t \\
&= \int_T\int_G\int_{t_0}^{t} \nabla\cdot(\mathbf{u}^\dagger(t)\cdot\dot{\mathbf{C}}(t-\tau):\nabla\mathbf{u}(\tau)) - \nabla\cdot(\mathbf{u}(\tau)\cdot\dot{\mathbf{C}}(t-\tau):\nabla\mathbf{u}^\dagger(t)) + \mathbf{u}(\tau)\cdot[\nabla\cdot(\dot{\mathbf{C}}(t-\tau):\nabla\mathbf{u}^\dagger(t))]\ \mathrm{d}\tau\ \mathrm{d}^3\mathbf{x}\ \mathrm{d}t \\
&= \int_T\int_G\int_{t_0}^{t} \nabla\cdot(\mathbf{u}^\dagger(t)\cdot\dot{\mathbf{C}}(t-\tau):\nabla\mathbf{u}(\tau))\ \mathrm{d}\tau\ \mathrm{d}^3\mathbf{x}\ \mathrm{d}t \\
&- \int_T\int_G\int_{t_0}^{t} \nabla\cdot(\mathbf{u}(\tau)\cdot\dot{\mathbf{C}}(t-\tau):\nabla\mathbf{u}^\dagger(t))\ \mathrm{d}\tau\ \mathrm{d}^3\mathbf{x}\ \mathrm{d}t \\
&+ \int_T\int_G\int_{t_0}^{t} \mathbf{u}(\tau)\cdot[\nabla\cdot(\dot{\mathbf{C}}(t-\tau):\nabla\mathbf{u}^\dagger(t))]\ \mathrm{d}\tau\ \mathrm{d}^3\mathbf{x}\ \mathrm{d}t.
\end{align*}$$

By Gauss' Theorem,
$$\begin{align*}
\gamma &= \int_{\partial G} \left\{ \int_{t_0}^{t_1}\int_{t_0}^{t} \mathbf{u}^\dagger(t)\cdot\dot{\mathbf{C}}(t-\tau):\nabla\mathbf{u}(\tau)\ \mathrm{d}\tau\  \mathrm{d}t\ \right\} \cdot \mathbf{n}\  \mathrm{d}^2\mathbf{x} \\
&- \int_{\partial G} \left\{ \int_{t_0}^{t_1}\int_{t_0}^{t} \mathbf{u}(\tau)\cdot\dot{\mathbf{C}}(t-\tau):\nabla\mathbf{u}^\dagger(t)\ \mathrm{d}\tau\  \mathrm{d}t\ \right\} \cdot \mathbf{n}\  \mathrm{d}^2\mathbf{x} \\
&+ \int_T\int_G\int_{t_0}^{t} \mathbf{u}(\tau)\cdot[\nabla\cdot(\dot{\mathbf{C}}(t-\tau):\nabla\mathbf{u}^\dagger(t))]\ \mathrm{d}\tau\ \mathrm{d}^3\mathbf{x}\ \mathrm{d}t.
\end{align*}$$

Since $\int_{t_0}^{t_1}\int_{t_0}^{t}\ \mathrm{d}\tau\ \mathrm{d}t = \int_{t_0}^{t_1}\int_{\tau}^{t_1}\ \mathrm{d}t\ \mathrm{d}\tau$,
$$\begin{align*}
\gamma &= \int_{\partial G} \left\{ \int_{t_0}^{t_1} \mathbf{u}^\dagger(t)\cdot \left[\int_{t_0}^{t} \dot{\mathbf{C}}(t-\tau):\nabla\mathbf{u}(\tau) \mathrm{d}\tau\right]\ \mathrm{d}t\ \right\} \cdot \mathbf{n}\  \mathrm{d}^2\mathbf{x} \\
&- \int_{\partial G} \left\{ \int_{t_0}^{t_1} \mathbf{u}(\tau)\cdot \left[\int_{\tau}^{t_1} \dot{\mathbf{C}}(t-\tau):\nabla\mathbf{u}^\dagger(t) \mathrm{d}t \right]\ \mathrm{d}\tau\ \right\} \cdot \mathbf{n}\  \mathrm{d}^2\mathbf{x} \\
&+ \int_T\int_G\int_{t_0}^{t} \mathbf{u}(\tau)\cdot[\nabla\cdot(\dot{\mathbf{C}}(t-\tau):\nabla\mathbf{u}^\dagger(t))]\ \mathrm{d}\tau\ \mathrm{d}^3\mathbf{x}\ \mathrm{d}t = A-B+C.
\end{align*}$$

- Since the boundary condition $\mathbf{\sigma}\cdot\mathbf{n} = 0$ on $\partial G$, then $A=0$.
- Since so far we only have terminal conditions for adjoint condition, $\mathbf{u}^\dagger |_{t\geq t_1} = \dot{\mathbf{u}}^\dagger |_{t\geq t_1} = \mathbf{0}$, we are free to impose a boundary condition $\mathbf{n}\cdot\mathbf{\sigma}^\dagger |_{\mathbf{x}\in\partial G} = 0$. Where $\mathbf{\sigma}^\dagger(t) = \int_t^{t_1} \dot{\mathbf{C}}(t-\tau):\nabla\mathbf{u}^\dagger(t)\ \mathrm{d}\tau$. Then $B=0$.

Thus,

$$\begin{align*}
<\mathbf{u}\cdot\nabla_uL^\dagger \mathbf{u}^\dagger> &= <\mathbf{u}^\dagger\cdot \nabla_uL\mathbf{u}>\\
&= \int_T\int_G \mathbf{\rho}\mathbf{u}^\dagger(t)\cdot \ddot{\mathbf{u}}(t)\ \mathrm{d}^3 \mathbf{x} \mathrm{d}t - \int_T\int_G \mathbf{u}^\dagger(t) \cdot [ \nabla\cdot \int_{t_0}^{t}\dot{\mathbf{C}}(t-\tau):\nabla\mathbf{u}(\tau)\ \mathrm{d}\tau]\ \mathrm{d}^3\mathbf{x}\ \mathrm{d}t \\
&= <\mathbf{u}^\dagger\cdot\rho\ddot{\mathbf{u}}> - <\mathbf{u}\cdot(\nabla\cdot\mathbf{\sigma}^\dagger)> \\
&= <\mathbf{u}\cdot\mathbf{\rho}\ddot{\mathbf{u}}^\dagger> - <\mathbf{u}\cdot(\nabla\cdot\mathbf{\sigma}^\dagger)> \\
&= <\mathbf{u}\cdot(\mathbf{\rho}\ddot{\mathbf{u}}^\dagger - (\nabla\cdot\mathbf{\sigma}^\dagger))>.
\end{align*}$$

Then,

$$\nabla_uL^\dagger\mathbf{u}^\dagger = \mathbf{\rho}\ddot{\mathbf{u}}^\dagger - (\nabla\cdot\mathbf{\sigma}^\dagger) = -\nabla_u\chi_1^\dagger.$$

# Adjoint Source Term

Suppose we are working on the $L_2$ waveform misfit function,
$$\begin{align*}
\chi(\mathbf{m}) &= \int_G\int_T \frac{1}{2} [\mathbf{u}^0(\mathbf{x},t) - \mathbf{u}(\mathbf{m};\mathbf{x},t)]^2 \delta(\mathbf{x} - \mathbf{x^r})\ \mathrm{d}t \ \mathrm{d}^3\mathbf{x} \\
&= \int_G\int_T \chi_1(\mathbf{m};\mathbf{x},t)\ \mathrm{d}t \ \mathrm{d}^3\mathbf{x}
\end{align*}$$
with $\chi_1(\mathbf{m};\mathbf{x},t) = \frac{1}{2} [\mathbf{u}^0(\mathbf{x},t) - \mathbf{u}(\mathbf{m};\mathbf{x},t)]^2 \delta(\mathbf{x} - \mathbf{x^r})$.

We let the adjoint source $\mathbf{f}^\dagger(\mathbf{x},t)$ be,
$$\mathbf{f}^\dagger(\mathbf{x},t) = -\nabla_u\chi_1 = [\mathbf{u}^0(\mathbf{x},t) - \mathbf{u}(\mathbf{m};\mathbf{x},t)]\delta(\mathbf{x} - \mathbf{x^r}).$$

Then the adjoint equation is,
$$\begin{align*}
\nabla_uL^\dagger\mathbf{u}^\dagger &= -\nabla_u\chi_1^\dagger, \\
\mathbf{\rho}\ddot{\mathbf{u}}^\dagger - (\nabla\cdot\mathbf{\sigma}^\dagger) &= [\mathbf{u}^0(\mathbf{x},t) - \mathbf{u}(\mathbf{m};\mathbf{x},t)]\delta(\mathbf{x} - \mathbf{x^r}).
\end{align*}$$
With the initial and boundary conditions,
$$\mathbf{u}^\dagger |_{t\geq t_1} = \dot{\mathbf{u}}^\dagger |_{t\geq t_1} = \mathbf{0},\quad \mathbf{n}\cdot\mathbf{\sigma}^\dagger |_{\mathbf{x}\in\partial G} = 0.$$
By solving adjoint equation with time reverse, we can have $\mathbf{u}^\dagger$.



# Kernels

## Gerenal Case
For the derivative of misfit function,
$$\begin{align*}
\nabla_m\chi\delta\mathbf{m} &= <\mathbf{u}^\dagger \cdot \nabla_mL\delta\mathbf{m}> \\
&= \int_G\int_T \mathbf{u}^\dagger(t) \cdot \left[\delta\boldsymbol{\rho}\ddot{\mathbf{u}}(t) - \nabla\cdot\int_{t_0}^t \delta\dot{\mathbf{C}}(t-\tau):\nabla\mathbf{u}(\tau)\ \mathrm{d}\tau\right]\ \mathrm{d}t\ \mathrm{d}^3\mathbf{x} \\
&= \int_G\int_T \delta\boldsymbol{\rho}\mathbf{u}^\dagger(t)\cdot\ddot{\mathbf{u}}(t)\ \mathrm{d}t\ \mathrm{d}^3\mathbf{x} - \int_G\int_T \mathbf{u}^\dagger(t) \cdot \left[\nabla\cdot\int_{t_0}^t \delta\dot{\mathbf{C}}(t-\tau):\nabla\mathbf{u}(\tau)\ \mathrm{d}\tau\right]\ \mathrm{d}t\ \mathrm{d}^3\mathbf{x} \\
&= A - B.
\end{align*}$$
With the initial conditions of original equation and adjoint equation, we have,
$$\begin{align*}
A &=
\int_G\int \delta\boldsymbol{\rho}\mathbf{u}^\dagger(t) \cdot \mathrm{d}\dot{\mathbf{u}}(t)\ \mathrm{d}^3\mathbf{x} \\
&= \int_G \delta\boldsymbol{\rho} \left[\mathbf{u}^\dagger(t)\cdot\dot{\mathbf{u}}(t)|_{t_0}^{t_1} - \int_T \dot{\mathbf{u}}^\dagger(t)\cdot\dot{\mathbf{u}}(t)\ \mathrm{d}t\right]\ \mathrm{d}^3\mathbf{x} \\
&= -\int_G\int_T \delta\boldsymbol{\rho}\dot{\mathbf{u}}^\dagger(t)\cdot\dot{\mathbf{u}}(t)\ \mathrm{d}t\ \mathrm{d}^3\mathbf{x}.
\end{align*}$$
Since,
$$\begin{align*}
\mathbf{u}^\dagger \cdot (\nabla \cdot (\delta \dot{\mathbf{C}} : \nabla \mathbf{u})) &= \mathbf{u}^\dagger \cdot \left(\nabla \cdot \left((\delta C_{ijkl}\mathbf{e}_i\mathbf{e}_j\mathbf{e}_k\mathbf{e}_l) : \left(\mathbf{e}_m\mathbf{e}_n\frac{\partial u_n}{\partial x_m}\right)\right)\right) \\
&=  \mathbf{u}^\dagger \cdot \left(\nabla \cdot \left( \mathbf{e}_i\mathbf{e}_j \delta C_{ijkl} \frac{\partial u_k}{\partial x_l}\right)\right) \\
&= \mathbf{u}^\dagger \cdot \left(\left(\mathbf{e}_m \frac{\partial}{\partial x_m}\right) \cdot \left( \mathbf{e}_i\mathbf{e}_j \delta C_{ijkl} \frac{\partial u_k}{\partial x_l}\right)\right) \\
&= (\mathbf{e}_m u_m^\dagger) \cdot \left(\mathbf{e}_j \frac{\partial}{\partial x_i}\left(\delta C_{ijkl} \frac{\partial u_k}{\partial x_l} \right)\right) \\
&= u_j^\dagger \frac{\partial}{\partial x_i} \left(\delta C_{ijkl}\frac{\partial u_k}{\partial x_l}\right).
\end{align*}$$
And also,
$$\begin{align*}
\nabla\mathbf{u}^\dagger : \delta\dot{\mathbf{C}} : \nabla\mathbf{u} &= \left(\mathbf{e}_m\mathbf{e}_n\frac{\partial u_n^\dagger}{\partial x_m}\right) : (\mathbf{e}_i\mathbf{e}_j\mathbf{e}_k\mathbf{e}_l\delta C_{ijkl}) : \nabla\mathbf{u}\\
&= \left(\mathbf{e}_k\mathbf{e}_l\frac{\partial u_i^\dagger}{\partial x_j} \delta C_{ijkl}\right) : \left(\mathbf{e}_m\mathbf{e}_n\frac{\partial u_n}{\partial x_m}\right) \\
&= \frac{\partial u_i^\dagger}{\partial x_j} \delta C_{ijkl} \frac{\partial u_k}{\partial x_l}.
\end{align*}$$
Since in previous step, we have proved $\mathbf{u}^\dagger \cdot (\nabla \cdot (\delta \dot{\mathbf{C}} : \nabla \mathbf{u})) + \nabla\mathbf{u}^\dagger : \delta\dot{\mathbf{C}} : \nabla\mathbf{u} = \nabla\cdot (\mathbf{u}^\dagger \cdot \dot{\mathbf{C}} : \nabla\mathbf{u})$. And the $\iiint \nabla\cdot (\mathbf{u}^\dagger \cdot \dot{\mathbf{C}} : \nabla\mathbf{u}) \ \mathrm{d}\tau\ \mathrm{d}t \ \mathrm{d}^3\mathbf{x} = 0$.
Then,
$$B = -\int_T\int_G\int_{t_0}^{t_1} \nabla\mathbf{u}^\dagger : \delta\dot{\mathbf{C}} : \nabla\mathbf{u} \ \mathrm{d}\tau\ \mathrm{d}t \ \mathrm{d}^3\mathbf{x}.$$
So,
$$\nabla_m\chi\delta\mathbf{m} = A - B = -\int_G\int_T \delta\boldsymbol{\rho}\dot{\mathbf{u}}^\dagger(t)\cdot\dot{\mathbf{u}}(t)\ \mathrm{d}t\ \mathrm{d}^3\mathbf{x} + \int_T\int_G\int_{t_0}^{t_1} \nabla\mathbf{u}^\dagger : \delta\dot{\mathbf{C}} : \nabla\mathbf{u} \ \mathrm{d}\tau\ \mathrm{d}t \ \mathrm{d}^3\mathbf{x}$$

## Perfectly Elastic and Isotropic Medium
For perfectly elastic and isotropic medium, we have the following simplification, $\mathbf{C}(\mathbf{x},t) = \mathbf{C}(\mathbf{x})H(t)$ and $\delta\mathbf{C}(\mathbf{x},t) = \delta\mathbf{C}(\mathbf{x})H(t)$. In this case,
$$\sigma(\mathbf{x},t) = \int_{-\infty}^{\infty} \dot{\mathbf{C}}(\mathbf{x},t-t') : \nabla\mathbf{u}(\mathbf{x},t')\ \mathrm{d}t' = \mathbf{C} : \nabla\mathbf{u}(\mathbf{x},t).$$
Then,
$$\nabla_m\chi\delta\mathbf{m} = -\int_G\int_T \delta\boldsymbol{\rho}\dot{\mathbf{u}}^\dagger(t)\cdot\dot{\mathbf{u}}(t)\ \mathrm{d}t\ \mathrm{d}^3\mathbf{x} + \int_T\int_G \nabla\mathbf{u}^\dagger(t) : \delta\mathbf{C} : \nabla\mathbf{u}(t) \ \mathrm{d}t \ \mathrm{d}^3\mathbf{x}.$$

Since $C_{ijkl} = \lambda\delta_{ij}\delta_{kl} + \mu\delta_{ik}\delta_{jl} + \mu\delta_{il}\delta_{jk}$, $\lambda$ and $\mu$ are Lamé parameters. We have proved,
$$\begin{align*}
\nabla\mathbf{u}^\dagger : \delta\mathbf{C} : \nabla\mathbf{u} &= \frac{\partial u_i^\dagger}{\partial x_j} \delta C_{ijkl} \frac{\partial u_k}{\partial x_l} \\
&= \frac{\partial u_i^\dagger}{\partial x_j} [\delta\lambda\delta_{ij}\delta_{kl} + \delta\mu\delta_{ik}\delta_{jl} + \delta\mu\delta_{il}\delta_{jk}]  \frac{\partial u_k}{\partial x_l} \\
&= \delta\lambda\frac{\partial u_i^\dagger}{\partial x_j}\frac{\partial u_k}{\partial x_k} + \delta\mu\frac{\partial u_i^\dagger}{\partial x_j}\frac{\partial u_i}{\partial x_j} + \delta\mu\frac{\partial u_i^\dagger}{\partial x_j}\frac{\partial u_j}{\partial x_i} \\
&= \delta\lambda(\nabla\cdot\mathbf{u}^\dagger)(\nabla\cdot\mathbf{u}) + \delta\mu[(\nabla\mathbf{u}^\dagger):(\nabla\mathbf{u}) + (\nabla\mathbf{u}^\dagger) : (\nabla\mathbf{u})^T].
\end{align*}$$
Then we have,
$$\begin{align*}
\nabla_m\chi\delta\mathbf{m} &= - \int_G\int_T\delta\boldsymbol{\rho}\dot{\mathbf{u}}^\dagger(t)\cdot\dot{\mathbf{u}}(t)\ \mathrm{d}t\ \mathrm{d}^3\mathbf{x} + \int_G\int_T \delta\boldsymbol{\lambda}(\nabla\cdot\mathbf{u}^\dagger)(\nabla\cdot\mathbf{u})\ \mathrm{d}t\ \mathrm{d}^3\mathbf{x}
+ \int_G\int_T \delta\boldsymbol{\mu}[(\nabla\mathbf{u}^\dagger):(\nabla\mathbf{u}) + (\nabla\mathbf{u}^\dagger) : (\nabla\mathbf{u})^T]\ \mathrm{d}t\ \mathrm{d}^3\mathbf{x} \\
&= \int_G K_\rho^0 \delta\boldsymbol{\rho}\ \mathrm{d}^3\mathbf{x} + \int_G K_\lambda^0 \delta\boldsymbol{\lambda}\ \mathrm{d}^3\mathbf{x} + \int_G K_\mu^0 \delta\boldsymbol{\mu}\ \mathrm{d}^3\mathbf{x} \\
&= \nabla_\rho\chi\delta\boldsymbol{\rho} +
\nabla_\lambda\chi\delta\boldsymbol{\lambda} +
\nabla_\mu\chi\delta\boldsymbol{\mu}.
\end{align*}$$
With sensitivity kernal:
$$\begin{align*}
K_\rho^0 &= -\int_T \dot{\mathbf{u}}^\dagger(t)\cdot\dot{\mathbf{u}}(t)\ \mathrm{d}t, \\
K_\lambda^0 &= \int_T (\nabla\cdot\mathbf{u}^\dagger)(\nabla\cdot\mathbf{u})\ \mathrm{d}t, \\
K_\mu^0 &= \int_T [(\nabla\mathbf{u}^\dagger):(\nabla\mathbf{u}) + (\nabla\mathbf{u}^\dagger) : (\nabla\mathbf{u})^T]\ \mathrm{d}t.
\end{align*}$$
