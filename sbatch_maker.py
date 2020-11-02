import numpy as np

def make_sbatch_file(file_path, job_name, code_folder, output_folder_name, params, nthreads, mkdirs=True):

    exec_command = f'../exec {params["nr"]} {params["nth"]} {params["nph"]} {params["cfl"]} {params["chi"]} {params["A"]} {params["r0"]} {params["w"]} {params["mu_s"]}'
    dir_name = f'{output_folder_name}_{params["nr"]}_{params["nth"]}_{params["nph"]}_{params["chi"]}'
    if mkdirs:
        mkdir_line = f'mkdir output/{dir_name} ; cd output/{dir_name} ; mkdir checkpoint_data'
    else:
        mkdir_line = f'cd output/{dir_name}'

    with open(file_path, "w") as file:
        file.write(f'''#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --output=/home/goncaloa/{code_folder}/output/out/{job_name}%j.out
#SBATCH --error=/home/goncaloa/{code_folder}/output/err/{job_name}%j.err
#SBATCH --nodes=1
#SBATCH --ntasks={nthreads}
#SBATCH --mail-user=goncalo.j.c.andrade@tecnico.ulisboa.pt
#SBATCH --mail-type=ALL
#SBATCH --time=72:00:00
#SBATCH --mem=32G

RUNPATH=/home/goncaloa/{code_folder} ; cd $RUNPATH
{mkdir_line}
export OMP_NUM_THREADS={nthreads}
{exec_command}
''')


def xx0_to_r(xx0s, rmax, sinhw):
    return rmax * np.sinh(xx0s / sinhw) / np.sinh(1 / sinhw)


def points_inside_horizon(nr, rmax, sinhw, M, chi):

    a = M * chi
    rp = M + np.sqrt(M**2 - a**2)
    rh = rp / 4

    i = np.arange(nr)
    xx0s = (i + 1 / 2) / nr
    rs = xx0_to_r(xx0s, rmax, sinhw)

    return rs[np.where(rs < rh)].size


def find_nr_from_horizon_pts(horizon_pts, rmax, sinhw, M, chi, increment=2):

    nr = 0
    n_horizon = 0

    while n_horizon < horizon_pts:

        n_horizon = 0
        nr += increment
        n_horizon = points_inside_horizon(nr, rmax, sinhw, M, chi)

    return nr


def find_nr_from_gaussian_width(w, rmax, sinhw, increment=2):

    min_distance = w / 2
    nr = 0
    distance = 1e10 # something ridiculous just to meet the while criterium

    while distance > min_distance:

        nr += increment
        i = np.arange(nr)
        xx0s = (i + 1 / 2) / nr
        rs = xx0_to_r(xx0s, rmax, sinhw)
        distance = rs[-1] - rs[-2]

    return nr
        

if __name__ == "__main__":
    
    M = 1
    chi = 0
    w = 0.5
    rmax = 500

    sinhw = 0.15
    sinhws = []
    nrs = []
    pts = []

    while sinhw < 1:

        sinhw += 0.01
        nr = find_nr_from_gaussian_width(w, rmax, sinhw)
        horizon_pts = points_inside_horizon(nr, rmax, sinhw, M, chi)
        if horizon_pts <= 20:
            break

    print(sinhw, nr, horizon_pts)

