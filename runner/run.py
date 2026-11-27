"""Запуск одного расчёта выбранной PDE-модели."""
from pde import FieldCollection, MemoryStorage
from pde.solvers.scipy import ScipySolverError

from exceptions.film_vanishes_exc import FilmVanishedError
from pde_model.shklyaev_2008 import Shklyaev2008, Shklyaev2008_VIII_B
from pde_model.volodin_2025 import Volodin2025
from pde_model.weaklynonlinear_2026 import WeaklyNonlinear2026

# имя модели -> класс уравнения
MODELS = {
    'volodin_2025': Volodin2025,
    'shklyaev_2008': Shklyaev2008,
    'shklyaev_2008_viii_b': Shklyaev2008_VIII_B,
    'weaklynonlinear_2026': WeaklyNonlinear2026,
}

# модели с двумя полями (h и T); у остальных только h
MODELS_WITH_TEMPERATURE = {'volodin_2025', 'weaklynonlinear_2026'}


def run_simulation(pde_model, initial_h, initial_T, t_range, **params):
    """Решает уравнение и возвращает (результат, storage, is_solved).

    params — параметры модели (Ca, G0, b, ...), ровно те, что лежат в params.py.
    is_solved=False, если расчёт оборвался: плёнка исчезла или упал решатель.
    """
    equation = MODELS[pde_model](**params)

    if pde_model in MODELS_WITH_TEMPERATURE:
        state = FieldCollection([initial_h, initial_T])
    else:
        state = initial_h

    print(f'Модель: {pde_model}, t_range={t_range}')
    print(f'Параметры: {params}')

    storage = MemoryStorage()
    try:
        result = equation.solve(
            state,
            t_range=t_range,
            solver='scipy',
            method='Radau',
            rtol=1e-6,
            atol=1e-9,
            tracker=['progress', storage.tracker(interrupts=1)],
        )
    except FilmVanishedError as error:
        print(f'Расчёт остановлен: {error!r}')
        last_state = storage[-2] if len(storage) >= 2 else state
        return last_state, storage, False

    except ScipySolverError as error:
        print(f'Расчёт остановлен: {error!r}')
        last_state = storage[-2] if len(storage) >= 2 else state
        return last_state, storage, False

    return result, storage, True
