"""Параметры расчётов по умолчанию для каждой модели."""
from pathlib import Path

from copy import deepcopy


MODEL_PARAMS = {
    'volodin_2025': {
        't_range': 500,
        'Ca': 2.5,
        'G0': 3,
        'phi': 0,
        'Pr': 7,
        'Bi': 1,
        'Ma': 5.5,
        'omega': 1.5,
        'b': 4.5,
        'bc': {'x': 'periodic'},
    },
    'weaklynonlinear_2026': {
        't_range': 500,
        'Ca': 2.5,
        'G0': 3,
        'phi': 0,
        'Pr': 7,
        'Bi': 1,
        'Ma': 5.5,
        'omega': 1.5,
        'b': 4.5,
        'bc': {'x': 'periodic'},
    },
    'shklyaev_2008_viii_b': {
        't_range': 500,
        'V': 0,
        'k': 1.1,
        'G0': 3.33e-4,
        'omega': 0.02,
        'bc': {'x': 'periodic'},
    },
}


def get_simulation_params(pde_model='volodin_2025'):
    """Копия параметров модели (можно менять, оригинал не пострадает)."""
    return deepcopy(MODEL_PARAMS[pde_model])


def decimal_range(start, stop, step, ndigits=1):
    """Числа от start до stop включительно. Считает через индекс, чтобы не
    накапливалась погрешность (0.3, а не 0.30000000000000004 в имени папки)."""
    count = round((stop - start) / step) + 1
    return [round(start + i * step, ndigits) for i in range(count)]


PDE_MODEL = 'volodin_2025'
WEAKLY_NONLINEAR_MODEL = 'weaklynonlinear_2026'
RESULTS_DIR = Path('results')

# Начальное возмущение по умолчанию
DEFAULT_PERTURBATION_AMPLITUDE = 0.01
DEFAULT_WAVE_K = 2.0


# volodin_2025: серия по начальным условиям.
INITIAL_AMPLITUDES = (
    0.001,
    0.025,
    0.005,
    0.075,
    0.1,
    0.125,
    0.15,
    0.175,
    0.2,
    0.225,
)
WAVE_K_VALUES = decimal_range(2.0, 3.1, 0.1)

# volodin_2025: серия по внешнему воздействию
EXTERNAL_AMPLITUDES = decimal_range(0.0, 4.0, 0.1)
EXTERNAL_OMEGAS = decimal_range(0.5, 4.0, 0.1)

# weaklynonlinear_2026: базовая точка; варьируется один параметр,
# остальные берутся отсюда
WEAKLY_NONLINEAR_BASE_POINT = {'b': 1, 'Ma': 5.5, 'omega': 1.5}

K_WINDOW_081 = [0.66, 0.71, 0.76, 0.81, 0.86, 0.91, 0.96]  # k ≈ 0.81
K_WINDOW_091 = [0.76, 0.81, 0.86, 0.91, 0.96, 1.05, 1.10]  # k ≈ 0.91
K_WINDOW_228 = [2.15, 2.20, 2.25, 2.28, 2.34, 2.39, 2.42]  # k ≈ 2.28
K_WINDOW_289 = [2.75, 2.80, 2.85, 2.89, 2.94, 2.99, 3.02]  # k ≈ 2.89
K_WINDOW_312 = [2.95, 3.00, 3.05, 3.09, 3.12, 3.15, 3.20, 3.24]  # k ≈ 3.12

K_WINDOW_086 = [0.71, 0.76, 0.81, 0.86, 0.91, 0.96, 1.05]  # k ≈ 0.86
K_WINDOW_255 = [2.41, 2.46, 2.50, 2.55, 2.57, 2.60, 2.62]  # k ≈ 2.55

K_WINDOW_407 = [3.95, 4.00, 4.05, 4.07, 4.12, 4.15, 4.20]  # k ≈ 4.07

K_WINDOW_095 = [0.81, 0.86, 0.91, 0.95, 1.00, 1.05, 1.10]  # k ≈ 0.95
K_WINDOW_209 = [2.075, 2.08, 2.085, 2.09, 2.14, 2.19, 2.23]  # k ≈ 2.09

WEAKLY_NONLINEAR_SCANS = {
    # амплитуда
    'b': {
        0.1: [K_WINDOW_312],
        1.0: [K_WINDOW_081, K_WINDOW_289],
        2.0: [K_WINDOW_091, K_WINDOW_228],
    },
    # частота
    'omega': {
        0.5: [K_WINDOW_312],
        1.5: [K_WINDOW_081, K_WINDOW_289],
        2.0: [K_WINDOW_086, K_WINDOW_255],
    },
    # число Марангони
    'Ma': {
        6.0: [K_WINDOW_407],
        5.5: [K_WINDOW_081, K_WINDOW_289],
        5.0: [K_WINDOW_095, K_WINDOW_209],
        12.0: [K_WINDOW_407],
    },
}
