from io_utils.plotter import plot_results_to_file

from os import makedirs
from numpy import savetxt


def save_res_to_txt(storage, output_dir=None, is_solved=None):
    # TODO: fix if storage is 1D instead of 2D

    makedirs(output_dir, exist_ok=True)

    # создаем поддиректории
    h_dir = f'{output_dir}/h'
    T_dir = f'{output_dir}/T'
    makedirs(h_dir, exist_ok=True)
    makedirs(T_dir, exist_ok=True)

    # сохраняем поля
    for i, step in enumerate(storage.data):
        h_field = step[0].data
        T_field = step[1].data

        savetxt(f'{h_dir}/h_{i}.txt', h_field)
        savetxt(f'{T_dir}/T_{i}.txt', T_field)

    print(f'Сохранено {len(storage.data)} шагов в {output_dir}/')

    # создаем текстовый файл статуса
    status_path = f'{output_dir}/status.txt'
    with open(status_path, 'w') as f:
        if is_solved is True:
            f.write('solved')
        elif is_solved is False:
            f.write('not_solved')
        else:
            f.write('unknown')

    print(f'Статус записан в {status_path}')


def params_to_folder_name(params):
    """Имя папки из параметров расчёта (t_range и bc в имя не попадают)."""
    return '_'.join(
        f'{key}={value}'
        for key, value in params.items()
        if key not in ('t_range', 'bc')
    )


def save_results(storage, output_dir, is_solved):
    """Создаёт папку и сохраняет в неё данные и графики."""
    output_dir.mkdir(parents=True, exist_ok=True)
    save_res_to_txt(storage, str(output_dir), is_solved)
    plot_results_to_file(storage, str(output_dir), is_solved)
