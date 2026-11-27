"""Запуск серий расчётов.

Использование:
    python main.py [эксперимент]

Эксперименты:
    initial        сетка по амплитуде и k начального возмущения, volodin_2025 (по умолчанию)
    external       сетка по частоте и амплитуде внешнего воздействия, volodin_2025
    external-fast  то же, но перебор амплитуд останавливается на первой is_solved=True
    single         один расчёт, volodin_2025
    weak-b         скан по k при разных амплитудах b, weaklynonlinear_2026
    weak-omega     скан по k при разных частотах omega, weaklynonlinear_2026
    weak-ma        скан по k при разных числах Марангони Ma, weaklynonlinear_2026
"""

import argparse
from itertools import chain, product

from initial_conditions.init_fields import prepare_initial_conditions
from io_utils.saver import params_to_folder_name, save_results


from params import (
    DEFAULT_PERTURBATION_AMPLITUDE,
    DEFAULT_WAVE_K,
    EXTERNAL_AMPLITUDES,
    EXTERNAL_OMEGAS,
    INITIAL_AMPLITUDES,
    PDE_MODEL,
    RESULTS_DIR,
    WAVE_K_VALUES,
    WEAKLY_NONLINEAR_BASE_POINT,
    WEAKLY_NONLINEAR_MODEL,
    WEAKLY_NONLINEAR_SCANS,
    get_simulation_params,
)
from runner.run import run_simulation


def run_with_different_initial_conditions():
    """Сетка по амплитуде и волновому числу начального возмущения."""
    params = get_simulation_params(PDE_MODEL)
    root_dir = RESULTS_DIR / params_to_folder_name(params)

    for amplitude, wave_k in product(INITIAL_AMPLITUDES, WAVE_K_VALUES):
        initial_h, initial_T = prepare_initial_conditions(amplitude, wave_k)
        _, storage, is_solved = run_simulation(
            PDE_MODEL, initial_h, initial_T, **params
        )
        folder = f'k={wave_k:.2f}_a={amplitude}_is_solved={is_solved}'
        save_results(storage, root_dir / folder, is_solved)


def run_with_different_external_freqs_and_amplitudes(stop_at_boundary=False):
    """Сетка по частоте omega и амплитуде b внешнего воздействия.

    stop_at_boundary=True: для каждой omega перебор амплитуд прекращается на
    первой, где is_solved=True (для больших амплитуд результат уже не меняется).
    """
    params = get_simulation_params(PDE_MODEL)
    initial_h, initial_T = prepare_initial_conditions(
        DEFAULT_PERTURBATION_AMPLITUDE, DEFAULT_WAVE_K
    )

    for omega in EXTERNAL_OMEGAS:
        print(f'\n=== omega = {omega} ===')
        boundary_found = False

        for amplitude in EXTERNAL_AMPLITUDES:
            run_params = {**params, 'b': amplitude, 'omega': omega}
            _, storage, is_solved = run_simulation(
                PDE_MODEL, initial_h, initial_T, **run_params
            )
            folder = f"Ma={params['Ma']}_ampl={amplitude}_omfr={omega}_is_solved={is_solved}"
            save_results(storage, RESULTS_DIR / folder, is_solved)

            if stop_at_boundary and is_solved:
                print(
                    f'>>> Граница найдена: omega={omega}, b={amplitude}. Большие b пропускаем.'
                )
                boundary_found = True
                break

        if stop_at_boundary and not boundary_found:
            print(
                f'>>> Для omega={omega} граница не найдена (всё неустойчиво).'
            )


def run_single():
    """Один расчёт с параметрами по умолчанию."""
    params = get_simulation_params(PDE_MODEL)
    initial_h, initial_T = prepare_initial_conditions(
        DEFAULT_PERTURBATION_AMPLITUDE, DEFAULT_WAVE_K
    )
    _, storage, is_solved = run_simulation(
        PDE_MODEL, initial_h, initial_T, **params
    )
    save_results(storage, RESULTS_DIR, is_solved)


def run_weaklynonlinear_k_range(scan):
    """Скан по k вблизи крит. чисел. scan — какой параметр меняется вместе с k:
    'b', 'omega' или 'Ma' (ключ WEAKLY_NONLINEAR_SCANS)."""
    default_params = get_simulation_params(WEAKLY_NONLINEAR_MODEL)

    for scan_value, k_windows in WEAKLY_NONLINEAR_SCANS[scan].items():
        params = {
            **default_params,
            **WEAKLY_NONLINEAR_BASE_POINT,
            scan: scan_value,
        }

        for wave_k in chain.from_iterable(k_windows):
            print(f'\n=== {scan}={scan_value}, wave_k={wave_k} ===')
            initial_h, initial_T = prepare_initial_conditions(
                DEFAULT_PERTURBATION_AMPLITUDE, wave_k
            )
            _, storage, is_solved = run_simulation(
                WEAKLY_NONLINEAR_MODEL, initial_h, initial_T, **params
            )
            folder = (
                f"wave_k={wave_k}_ext_amplitude={params['b']}_"
                f"Ma={params['Ma']}_Omega={params['omega']}_is_solved={is_solved}"
            )
            save_results(
                storage,
                RESULTS_DIR / WEAKLY_NONLINEAR_MODEL / folder,
                is_solved,
            )


EXPERIMENTS = {
    'initial': run_with_different_initial_conditions,
    'external': run_with_different_external_freqs_and_amplitudes,
    'external-fast': lambda: run_with_different_external_freqs_and_amplitudes(
        stop_at_boundary=True
    ),
    'single': run_single,
    'weak-b': lambda: run_weaklynonlinear_k_range('b'),
    'weak-omega': lambda: run_weaklynonlinear_k_range('omega'),
    'weak-ma': lambda: run_weaklynonlinear_k_range('Ma'),
}


def main():
    parser = argparse.ArgumentParser(description='Серии численных расчётов')
    parser.add_argument(
        'experiment',
        nargs='?',
        default='initial',
        choices=EXPERIMENTS,
        help='какой эксперимент запустить (по умолчанию: initial)',
    )
    args = parser.parse_args()
    EXPERIMENTS[args.experiment]()


if __name__ == '__main__':
    main()
