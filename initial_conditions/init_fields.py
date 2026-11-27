import numpy as np
from pde import CartesianGrid, ScalarField


def prepare_initial_conditions(perturbation_amplitude, wave_k, nx=100):
    """Начальные поля h и T: 1 + a * cos(k * x).

    Область — одна длина волны 2π/k, границы периодические, nx — число узлов.
    """
    wavelength = 2 * np.pi / wave_k
    grid = CartesianGrid([[0, wavelength]], [nx], periodic=True)
    x = grid.axes_coords[0]

    h0 = ScalarField(grid, 1 + perturbation_amplitude * np.cos(wave_k * x))
    T0 = ScalarField(grid, 1 + perturbation_amplitude * np.cos(wave_k * x))
    return h0, T0
