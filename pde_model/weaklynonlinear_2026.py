from pde import PDEBase, FieldCollection
from numpy import sqrt, sin, sinh, cos, cosh
import numpy as np

from pde.solvers.scipy import ScipySolverError
from exceptions.film_vanishes_exc import FilmVanishedError


class WeaklyNonlinear2026(PDEBase):
    """Solution to the Nonlinear PDE system Eq.(20)-Eq.(21):
    https://doi.org/10.1007/s12217-026-10289-6
    """

    def __init__(
        self,
        Ca,
        G0,
        phi,
        Pr,
        Bi,
        Ma,
        omega,
        b,
        bc,
    ):
        super().__init__()
        self.Ca = Ca
        self.G0 = G0
        self.phi = phi
        self.Pr = Pr
        self.Bi = Bi
        self.Ma = Ma
        self.omega = omega
        self.b = b
        self.bc = bc or {'x': 'periodic'}
        self.g = self._gamma0(omega)
        self.m1 = self._m1()
        self.m2 = self._m2()
        self.R1, self.R2 = self._R1_R2()
        self.output_params()

    def output_params(self):
        print()
        print(
            f'Ca={self.Ca}, G0={self.G0}, phi={self.phi}, Pr={self.Pr}, Bi={self.Bi}, '
            f'Ma={self.Ma}, omega={self.omega}, b={self.b}, bc={self.bc}, '
            f'g={self.g}, m1={self.m1}, m2={self.m2}, R1={self.R1}, R2={self.R2}'
        )
        print()

    def _gamma0(self, omega: float) -> float:
        return sqrt(2.0 * omega)

    def _psi_phi(self):
        psi_plus = cosh(self.g) + cos(self.g)
        psi_minus = cosh(self.g) - cos(self.g)
        phi_plus = sinh(self.g) + sin(self.g)
        phi_minus = sinh(self.g) - sin(self.g)
        return psi_plus, psi_minus, phi_plus, phi_minus

    def _R1_R2(self):
        psi_p, psi_m, phi_p, phi_m = self._psi_phi()
        R1 = (
            3.0
            * (
                4.0 * sinh(self.g) * sin(self.g)
                - self.g * (psi_p * phi_p - psi_m * phi_m)
            )
            / (2.0 * self.g**2 * psi_p**2)
        )
        R2 = -1.0 / 3.0 + (11.0 * phi_m - 3.0 * self.g * psi_m) / (
            self.g**3 * psi_p
        )
        return R1, R2

    def _m1(self) -> float:
        return 1.0 - (2.0 / self.g) * (sin(self.g) + sinh(self.g)) / (
            cos(self.g) + cosh(self.g)
        )

    def _m2(self) -> float:
        return (sinh(self.g) ** 2 - sin(self.g) ** 2) / (
            cos(self.g) + cosh(self.g)
        ) ** 2

    def _n1(self) -> float:
        g2 = self.g / 2.0
        return (
            sqrt(2.0)
            * self.g
            * (cos(self.g) + 1.0)
            * (
                sinh(self.g) * (1.0 - cos(self.g) + cosh(self.g))
                + sin(self.g) * (1.0 - cosh(self.g))
            )
            + 1.5 * self.g * sinh(self.g) ** 2 * sin(self.g)
            - sinh(self.g) ** 3
            + 4.0 * sqrt(2.0) * self.g * (cosh(g2) - sin(g2))
        )

    def _n2(self) -> float:
        g2 = self.g / 2.0
        return (
            sqrt(2.0)
            * self.g
            * sin(self.g)
            * (3.0 * (cosh(self.g) - 1.0) * cosh(g2) - cos(self.g) - 1.0)
        )

    def _m3(self, xi1, dxi1_dx) -> float:
        dxi = xi1**2 * dxi1_dx  # ξ₁² · ∂ξ̄₁/∂X

        return (
            sqrt(2.0)
            / (2.0 * (cosh(self.g) + cos(self.g)) ** 3)
            * (dxi * self._n1() + self._n2())
        )

    def evolution_rate(self, state: FieldCollection, t: float = 0.0):
        h, T = state

        if np.any(np.isnan(T.data)) or np.any(np.isinf(T.data)):
            raise ScipySolverError(f'NaN/Inf in T at t={t}')
        if np.any(np.isnan(h.data)) or np.any(np.isinf(h.data)):
            raise ScipySolverError(f'NaN/Inf in h at t={t}')
        if np.any(h.data <= 1e-6):
            raise FilmVanishedError(f'Film vanished at t={t}')

        f = T - h

        dx_f = f.gradient(bc=self.bc)[0]    # ∂f₁/∂x

        dx_h = h.gradient(bc=self.bc)[0]    # ∂ξ₁/∂x
        d2x_h = h.laplace(bc=self.bc)       # ∂²ξ₁/∂x²
        dx_T = T.gradient(bc=self.bc)[0]    # ∂T₁/∂x

        # ---- Π⁽¹⁾ field and its gradient ----
        Pi1 = (
            self.G0 * h
            - self.Ca * d2x_h
            - ((self.omega**2 * self.b**2) / 2.0) * self.m1 * dx_h
        )

        dx_Pi1 = Pi1.gradient(bc=self.bc)[0]

        # ---- Π⁽²⁾ ----
        Pi2 = (
            # self.G * xi2
            # - self.Ca * d2x_xi2
            -(self.omega**2 * self.b**2)
            / 2.0
            * ((1.0 - self.m1) * h * dx_h + h * dx_h * self.m2)
        )

        dx_Pi2 = Pi2.gradient(bc=self.bc)[0]

        # ---- Π⁽³⁾ = B²Ω²/2 · m₃ ----
        m3 = self._m3(h, dx_h)
        Pi3 = (self.omega**2 * self.b**2) / 2.0 * m3

        rhs1 = (
            dx_Pi1 * h**2
            + dx_Pi2 * h
            + 1.0 / 3.0 * Pi3
            + self.Ma / 2.0 * h**2 * (dx_T - dx_h)
        )

        h_t = rhs1.gradient(bc=self.bc)[0] - (
            self.omega**2 * self.b**2
        ) / 2.0 * (self.R1 * dx_h**3 + self.R2 * dx_h * h * d2x_h)

        term1 = (
            -self.Bi / self.Pr * (T - h)
            + (h * dx_Pi1 + self.Ma * h * dx_f) * dx_f
        )

        term2 = (
            3.0 / 4.0 * h**2 * dx_Pi1
            + h / 2.0 * dx_Pi2
            + 1.0 / 3.0 * Pi3
            + self.Ma / 2.0 * h**2 * dx_f
        ).gradient(bc=self.bc)[0]

        T_t = term1 + term2

        return FieldCollection([h_t, T_t])
