import numpy as np
from scipy.special import erf
import multiprocessing as mp

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
#SBATCH --output=/home/goncaloa/{code_folder}/output/out/{job_name}_%j.out
#SBATCH --error=/home/goncaloa/{code_folder}/output/err/{job_name}_%j.err
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


def find_nr_from_gaussian_width_sinh(w, rmax, sinhw, increment=2):

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


def erfspherical(x0, rmax, erfw, offset):
    return - 1 / 2 * rmax * x0 * (-2 + erf(offset - erfw * x0) + erf(offset + erfw * x0))


def find_timestep_erf(rmax, erfw, offset, nr, nth, nph, cfl=0.5):
    
    ds_min = 1e38

    xx0min = 0
    xx0max = 1
    dxx0 = (xx0max - xx0min) / nr
    i0s = np.arange(nr)
    xx0s = xx0min + (i0s + 1 / 2) * dxx0

    xx1min = 0
    xx1max = np.pi
    dxx1 = (xx1max - xx1min) / nth
    i1s = np.arange(nth)
    xx1s = xx1min + (i1s + 1 / 2) * dxx1

    xx2min = -np.pi
    xx2max = np.pi
    dxx2 = (xx2max - xx2min) / nph
    i2s = np.arange(nph)
    xx2s = xx2min + (i2s + 1 / 2) * dxx2

    for i0 in i0s:
        xx0 = xx0s[i0]
        for i1 in i1s:
            xx1 = xx1s[i1]
            for i2 in i2s:
                xx2 = xx2s[i2]

                ds_dirn0 = dxx0*(-rmax*xx0*(2*erfw*np.exp(-(erfw*xx0 + offset)**2)/np.sqrt(np.pi) - \
                        2*erfw*np.exp(-(erfw * xx0 - offset)**2)/np.sqrt(np.pi))/2 - rmax*(-erf(erfw*xx0 - offset) + erf(erfw*xx0 + offset) - 2)/2)
                ds_dirn1 = -rmax*dxx1*xx0 * \
                    (-erf(erfw*xx0 - offset) + erf(erfw*xx0 + offset) - 2)/2
                ds_dirn2 = -rmax*dxx2*xx0*(-erf(erfw*xx0 - offset) + erf(erfw*xx0 + offset) - 2)*np.sin(xx1)/2
                ds_min = min(ds_min, ds_dirn0, ds_dirn1, ds_dirn2)

    return ds_min * cfl


def find_nr_for_erfw_and_off(params):

    # params = (erfw, off, nth, nph, nr_step, rmax, min_distance, horizon_pts)
    erfw, off, nth, nph, nr_step, rmax, min_distance, horizon_pts, rh = params

    print(f'erfw = {erfw}, off = {off}')

    nr = 0
    distance = 1e38  # something ridiculous
    while distance > min_distance:

        nr += nr_step
        i = np.arange(nr)
        xx0s = (i + 1 / 2) / nr
        rs = erfspherical(xx0s, rmax, erfw, off)

        hpts = rs[np.where(rs < rh)].size
        if hpts < horizon_pts:
            continue

        distances = np.zeros(nr - 1)
        for i in range(nr - 1):
            distances[i] = rs[i + 1] - rs[i]
        distance = np.amax(distances)

        rmax_correspondence = bool(
            np.absolute(rs[-1] - rmax) < rmax * 0.05)

    timestep = find_timestep_erf(rmax, erfw, off, nr, nth, nph)
    time_param = timestep / (nr * nth * nph)

    return erfw, off, nr, distance, hpts, time_param, rmax_correspondence



def find_erfspherical_params(nth, nph, rmax, M, chi, horizon_pts, gaussian_width):

    min_distance = gaussian_width / 2
    a = M * chi
    rp = M + np.sqrt(M**2 - a**2)
    rh = rp / 4

    nr_step   = 10
    min_erfw  = 1.0
    max_erfw  = 10.1
    erfw_step = 0.2
    min_off   = 0.1
    max_off   = 5.1
    off_step  = 0.1

    best_nr   = rmax * 1e10 # something ridiculous
    best_erfw = 0
    best_off  = 0
    best_dist = 0
    best_hpts = 0
    best_tp   = 0

    all_params = []
    erfw = min_erfw
    while erfw <= max_erfw:
        off = min_off
        while off <= max_off:
            all_params.append((erfw, off, nth, nph, nr_step,
                               rmax, min_distance, horizon_pts, rh))
            # if rmax_correspondence and time_param > best_tp:
            #     best_nr = nr
            #     best_erfw = erfw
            #     best_off = off
            #     best_dist = distance
            #     best_hpts = hpts
            #     best_tp = time_param
            off += off_step
        erfw += erfw_step

    p = mp.Pool(mp.cpu_count())
    erfws, offs, nrs, distances, hptss, time_params, rmax_correspondences = np.asarray(
        p.map(find_nr_for_erfw_and_off, all_params)).T
    p.close()
    p.join()

    rmax_correspondences = np.asarray(rmax_correspondences, dtype=bool)

    erfws = erfws[rmax_correspondences]
    offs = offs[rmax_correspondences]
    nrs = nrs[rmax_correspondences]
    distances = distances[rmax_correspondences]
    hptss = hptss[rmax_correspondences]
    time_params = time_params[rmax_correspondences]

    best_idx = np.argmax(time_params)
    best_erfw = erfws[best_idx]
    best_off  = offs[best_idx]
    best_nr   = nrs[best_idx]
    best_dist = distances[best_idx]
    best_hpts = hptss[best_idx]
    best_tp   = time_params[best_idx]

    return best_nr, best_erfw, best_off, best_dist, best_hpts, best_tp


def gaussianspherical(x0, rmax, A, B, sigma):
    return 1 / 2 * rmax * x0 * (2 - A * np.exp(- (B - x0)**2 / (2 * sigma**2)) - A * np.exp(- (B + x0)**2 / (2 * sigma**2)))


def find_timestep_gaussian(rmax, A, B, sigma, nr, nth, nph, cfl=0.5):

    ds_min = 1e38

    xx0min = 0
    xx0max = 1
    dxx0 = (xx0max - xx0min) / nr
    i0s = np.arange(nr)
    xx0s = xx0min + (i0s + 1 / 2) * dxx0

    xx1min = 0
    xx1max = np.pi
    dxx1 = (xx1max - xx1min) / nth
    i1s = np.arange(nth)
    xx1s = xx1min + (i1s + 1 / 2) * dxx1

    xx2min = -np.pi
    xx2max = np.pi
    dxx2 = (xx2max - xx2min) / nph
    i2s = np.arange(nph)
    xx2s = xx2min + (i2s + 1 / 2) * dxx2

    for i0 in i0s:
        xx0 = xx0s[i0]
        for i1 in i1s:
            xx1 = xx1s[i1]
            for i2 in i2s:
                xx2 = xx2s[i2]

                ds_dirn0 = dxx0*(rmax*xx0*(A*(-2*B + 2*xx0)*np.exp(-(B - xx0)**2/(2*sigma**2))/(2*sigma**2) + A*(2*B + 2*xx0)*np.exp(-(B + xx0)**2/(2*sigma**2))/(2*sigma**2))/2 + rmax*(-A*np.exp(-(B + xx0)**2/(2*sigma**2)) - A*np.exp(-(B - xx0)**2/(2*sigma**2)) + 2)/2)
                ds_dirn1 = rmax*dxx1*xx0*(-A*np.exp(-(B + xx0)**2/(2*sigma**2)) - A*np.exp(-(B - xx0)**2/(2*sigma**2)) + 2)/2
                ds_dirn2 = rmax*dxx2*xx0*(-A*np.exp(-(B + xx0)**2/(2*sigma**2)) - A*np.exp(-(B - xx0)**2/(2*sigma**2)) + 2)*np.sin(xx1)/2
                ds_min = min(ds_min, ds_dirn0, ds_dirn1, ds_dirn2)
            
    return ds_min * cfl


def find_nr_for_gauss_params(params):

    A, B, sig, nth, nph, nr_step, rmax, min_distance, horizon_pts, rh = params
    print(f'A = {A}, B = {B}, sig = {sig}')

    nr = 0
    distance = 1e38
    while distance > min_distance:

        nr += nr_step
        i = np.arange(nr)
        xx0s = (i + 1 / 2) / nr
        rs = gaussianspherical(xx0s, rmax, A, B, sig)

        hpts = rs[np.where(rs < rh)].size
        if hpts < horizon_pts:
            continue

        distances = np.zeros(nr - 1)
        for i in range(nr - 1):
            distances[i] = rs[i + 1] - rs[i]
        distance = np.amax(distances)

        rmax_correspondence = bool(np.absolute(rs[-1] - rmax) < rmax * 0.05)

    timestep = find_timestep_gaussian(rmax, A, B, sig, nr, nth, nph)
    time_param = timestep / (nr * nth * nph)

    return A, B, sig, nr, distance, hpts, time_param, rmax_correspondence


def find_gaussianspherical_params(nth, nph, rmax, M, chi, horizon_pts, gaussian_width):

    min_distance = gaussian_width / 2
    a = M * chi
    rp = M + np.sqrt(M**2 - a**2)
    rh = rp / 4

    nr_step = 10
    min_A = 0.5
    max_A = 1.5
    A_step = 0.2
    min_B = 0.01
    max_B = 0.21
    B_step = 0.02
    min_sig = 0.01
    max_sig = 0.41
    sig_step = 0.04

    best_nr = rmax * 1e10  # something ridiculous
    best_A = 0
    best_B = 0
    best_sig = 0
    best_dist = 0
    best_hpts = 0
    best_tp = 0

    all_params = [(A, B, sig, nth, nph, nr_step, rmax, min_distance, horizon_pts, rh) for A in np.arange(min_A, max_A, A_step) for B in np.arange(min_B, max_B, B_step) for sig in np.arange(min_sig, max_sig, sig_step)]
    
    p = mp.Pool(mp.cpu_count())
    As, Bs, sigs, nrs, distances, hptss, time_params, rmax_correspondences = np.asarray(p.map(
        find_nr_for_gauss_params, all_params)).T
    p.close()
    p.join()

    rmax_correspondences = np.asarray(rmax_correspondences, dtype=bool)

    As = As[rmax_correspondences]
    Bs = Bs[rmax_correspondences]
    sigs = sigs[rmax_correspondences]
    nrs = nrs[rmax_correspondences]
    distances = distances[rmax_correspondences]
    hptss = hptss[rmax_correspondences]
    time_params = time_params[rmax_correspondences]

    best_idx = np.argmax(time_params)
    best_A = As[best_idx]
    best_B = Bs[best_idx]
    best_sig = sigs[best_idx]
    best_nr = nrs[best_idx]
    best_dist = distances[best_idx]
    best_hpts = hptss[best_idx]
    best_tp = time_params[best_idx]

    return best_nr, best_A, best_B, best_sig, best_dist, best_hpts, best_tp


if __name__ == "__main__":

    # Best GaussianSpherical config
    M = 1
    chi = 0.95
    rmax = 1000
    nth = 16
    nph = 32
    horizon_pts = 15
    width = 10

#     best_nr, best_A, best_B, best_sig, best_dist, best_hpts, best_tp = find_gaussianspherical_params(
#         nth, nph, rmax, M, chi, horizon_pts, width)
#     with open('best_configs.txt', 'w') as file:
#         file.write(f'''Found best gaussian config!
#     nr             = {best_nr}
#     A              = {best_A}
#     B              = {best_B}
#     sig            = {best_sig}
#     time_parameter = {best_tp}
# These parameters result in a max distance between points of {best_dist} and {best_hpts} inside the horizon.
# ''')
    
    # Best ErfSpherical config
    # M = 1
    # chi = 0.0
    # rmax = 100
    # nth = 8
    # nph = 16
    # horizon_pts = 15
    # width = 0.5

    # find_erfspherical_params(nth, nph, rmax, M, chi, horizon_pts, width)
#     best_nr, best_erfw, best_off, best_dist, best_hpts, best_tp = find_erfspherical_params(nth, nph, rmax, M, chi, horizon_pts, width)
#     with open('best_configs.txt', 'a') as file:
#         file.write(f'''Found best erf config!
#     nr = {best_nr}
#     ERFW = {best_erfw}
#     OFFSET = {best_off}
#     time_parameter = {best_tp} gp M/s
# These parameters result in a max distance between points of {best_dist} and {best_hpts} points inside the horizon.
# ''')

#     old_erfw = 2.7
#     old_off = 1
#     old_nr = 580
#     old_timestep = find_timestep(rmax, old_erfw, old_off, old_nr, nth, nph)
#     old_time_param = old_timestep / (old_nr * nth * nph)
#     print(f'''Previous best parameters were:
#     nr = {old_nr}
#     ERFW = {old_erfw}
#     OFFSET = {old_off}
#     time_parameter = {old_time_param}
# These parameters result in a max distance between points of 0.24675 and 18 points inside the horizon.''')

    # Best SinhSpherical config
    # M = 1
    # chi = 0
    # w = 0.5
    # rmax = 100

    sinhw = 0.1
    sinhws = []
    nrs = []
    pts = []

    while sinhw < 1:

        sinhw += 0.01
        nr = find_nr_from_gaussian_width_sinh(width, rmax, sinhw)
        hpts = points_inside_horizon(nr, rmax, sinhw, M, chi)
        if hpts <= horizon_pts:
            break

    print(f'''Best SinhSpherical config:
        sinhw = {sinhw}
        nr = {nr}
        horizon_pts = {horizon_pts}\n''')

    # with open('best_configs.txt', 'a') as file:
    #     file.write(f'''Best SinhSpherical config:
    #     sinhw = {sinhw}
    #     nr = {nr}
    #     horizon_pts = {horizon_pts}\n''')

