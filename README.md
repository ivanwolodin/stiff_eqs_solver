# Numerical Solution of Stiff Differential Equations with `py-pde`

A scientific computing project focused on the numerical solution of stiff differential equations using the [`py-pde`](https://py-pde.readthedocs.io/en/stable/) Python package.

The repository contains numerical implementations of evolution equations governing the dynamics of a thin liquid film under vertical vibrations, including coupled hydrodynamic and thermal effects.

The computational framework was developed as part of my PhD research and supports the numerical investigations presented in my doctoral dissertation and related scientific publications.

## Scientific Background

The project addresses the numerical solution of stiff differential equations arising in mathematical models of fluid dynamics.

In particular, it considers the averaged dynamics of a thin liquid film subject to vertical vibrations and thermocapillary (Marangoni) effects.

The governing equations form a coupled nonlinear system of evolution equations. Their numerical solution can be challenging due to the presence of different characteristic time and spatial scales.

The [`py-pde`](https://py-pde.readthedocs.io/en/stable/) package provides a convenient framework for implementing and solving such systems using Python and established numerical integration methods.

## 1. A Stiff ODE

The repository includes a tool to solve system of stiff PDE.

Simple example of stiff DE is:
$$
\frac{dy}{dx}=-50(y-\cos x).
$$

This equation is discussed in:

> Hairer, E. & Wanner, G. (1996).  
> *Solving Ordinary Differential Equations II: Stiff and Differential-Algebraic Problems.*  
> Springer.

The example provides a simple benchmark for investigating numerical integration methods for stiff DE.

## 2. Governing Equations: Thin-Film Dynamics under Vertical Vibrations

The main application of the numerical framework is the solution of a coupled system of nonlinear evolution equations describing thermocapillary convection in a thin liquid film under vertical vibrations.

The model is based on the averaged equations derived in:

> Volodin, I. & Alabuzhev, A. (2025).  
> *Linear Stability of Marangoni Convection in a Thin Film under Vertical Vibrations.*  
> Microgravity Science and Technology, 37, 48.  
> https://doi.org/10.1007/s12217-025-10201-8

### Mathematical Model

The system describes the evolution of the dimensionless film thickness \(\bar{\xi}\) and the depth-averaged temperature \(\bar{T}\):


$$
\begin{align}
    \frac{\partial \bar{\xi}}{\partial \bar{t}}
    &= \frac{\partial}{\partial X} \left(
    \frac{\bar{\xi}^3}{3} \frac{\partial \bar{\Pi}}{\partial X}
    + \frac{\text{Ma}}{2} \bar{\xi}^2 \frac{\partial f}{\partial X} \right. \notag \\
    &\quad \left. - \frac{(b \Omega)^2}{2} \left(
    Q_1 \bar{\xi}^2 \left( \frac{\partial \bar{\xi}}{\partial X} \right)^3
    + Q_2 \bar{\xi}^2 \frac{\partial \bar{\xi}}{\partial X} \frac{\partial^2 \bar{\xi}}{\partial X^2}
    \right)
    \right),
\end{align}
$$

$$
\begin{align}
    \frac{\bar{\xi} \, \partial \bar{T}}{\partial \bar{t}}
    &= \frac{1}{\text{Pr}} \frac{\partial}{\partial X} \left( \bar{\xi} \frac{\partial \bar{T}}{\partial X} \right)
    - \frac{1}{2 \, \text{Pr}} \left( \frac{\partial \bar{\xi}}{\partial X} \right)^2
    - \frac{\beta}{\text{Pr}} f \notag \\
    &\quad + \left(
    \frac{\bar{\xi}^3}{3} \frac{\partial \bar{\Pi}}{\partial X}
    + \frac{\text{Ma}}{2} \bar{\xi}^2 \frac{\partial f}{\partial X}
    \right) \frac{\partial f}{\partial X} \notag \\
    &\quad + \frac{\partial}{\partial X} \left(
    \frac{\bar{\xi}^4}{8} \frac{\partial \bar{\Pi}}{\partial X}
    + \frac{\text{Ma}}{6} \bar{\xi}^3 \frac{\partial f}{\partial X}
    \right)
\end{align}
$$

where $\bar{\Pi} = G \bar{\xi} - C \frac{\partial^2 \bar{\xi}}{\partial X^2} + \frac{(b \Omega)^2}{2} \Re\{H\} - \varphi $, $f=\bar{T} - \bar{\xi}$ is the temperature perturbation at the film surface, and

$$
\begin{align*}
	H = - \frac{\partial}{\partial X} \left[ \bar{\xi} F(\theta \bar{\xi}) \frac{\partial \bar{\xi}}{\partial X} \right], F (\theta \bar{\xi}) = 1 - \frac{\tan(\theta \bar{\xi})}{\theta \bar{\xi}},
\end{align*}
$$
$$
\begin{align*}
	Q_1 = 3 \frac{4 \sinh \gamma \sin \gamma - \gamma (\psi_+ \varphi_+ - \psi_- \varphi_-)}{2 \gamma^2 \psi_+^2},
\end{align*}
$$
$$
\begin{align*}
	Q_2 = -\frac{1}{3} + \frac{11 \varphi_- - 3 \gamma \psi_-}{\gamma^3 \psi_+}, \gamma=\bar{\xi} \sqrt{2 \omega},
\end{align*}
$$
$$
\begin{align*}
	\psi_+ = \cosh \gamma + \cos \gamma,
	\psi_- = \cosh \gamma - \cos \gamma, \\
	\varphi_+ = \sinh \gamma + \sin \gamma,
	\varphi_- = \sinh \gamma - \sin \gamma.
\end{align*}
$$

The model incorporates gravitational, capillary, thermocapillary, and vibration-induced contributions to the film dynamics.

## 3. Numerical Implementation

The system is implemented in Python using the `py-pde` package.

The numerical framework is designed to handle the coupled evolution equations and investigate their long-time dynamics.

Particular attention is paid to the numerical challenges associated with stiff differential equations, including:

- Different characteristic time scales.
- Nonlinear coupling between hydrodynamic and thermal fields.
- Higher-order spatial derivatives.
- Numerical stability and accuracy of time integration.
- Dependence of computational efficiency on the choice of integration method.

The implementation provides a basis for investigating the averaged dynamics of thin liquid films and reproducing the numerical results reported in the associated scientific studies.

## 4. Scientific Publications and Dissertation

The numerical investigations carried out using this computational framework constitute part of my PhD dissertation.

The results obtained with the implemented model have also been published in the following paper:

> Volodin, I. & Alabuzhev, A. (2026).  
> *Averaged Dynamics of Marangoni Convection in a Thin Liquid Film Under Vertical Vibrations.*  
> Microgravity Science and Technology, 38, 70.  
> https://doi.org/10.1007/s12217-026-10289-6

The repository provides a computational implementation of the mathematical model used in this research.

The dissertation, including the theoretical background and numerical results obtained using this solver, is available in a separate repository:

**[My PhD thesis — Dynamics of Free Deformable Surface in Vertical Vibrational Field](https://github.com/ivanwolodin/thesis)**

## 5. References

1. Hairer, E. & Wanner, G. (1996). *Solving Ordinary Differential Equations II: Stiff and Differential-Algebraic Problems.* Springer.
2. Volodin, I. & Alabuzhev, A. (2025). *Linear Stability of Marangoni Convection in a Thin Film under Vertical Vibrations.* Microgravity Science and Technology, 37, 48. https://doi.org/10.1007/s12217-025-10201-8
3. Volodin, I. & Alabuzhev, A. (2026). *Averaged Dynamics of Marangoni Convection in a Thin Liquid Film Under Vertical Vibrations.* Microgravity Science and Technology, 38, 70. https://doi.org/10.1007/s12217-026-10289-6

## Software and Computational Framework

All numerical computations in this project are based on [`py-pde`](https://py-pde.readthedocs.io/en/stable/), an open-source Python package for solving partial differential equations (PDEs).

The package provides a flexible framework for defining differential operators, implementing PDE systems, and performing numerical time integration. It serves as the primary computational tool used to implement and solve the coupled nonlinear evolution equations considered in this research.

**Official documentation:** https://py-pde.readthedocs.io/en/stable/

**Software publication:**

> Zwicker, D. (2020).  
> *py-pde: A Python package for solving partial differential equations.*  
> Journal of Open Source Software, 5(48), 2158.  
> https://doi.org/10.21105/joss.02158

The numerical implementation presented in this repository builds upon the functionality provided by this package, extending its application to the study of stiff coupled evolution equations arising in thin-film hydrodynamics.

## Acknowledgements

This project relies on the `py-pde` package developed by David Zwicker. I gratefully acknowledge the author and contributors for providing an open-source computational framework that made this research possible.