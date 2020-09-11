from pathlib import Path


def make_sbatch_file(file_path, job_name, code_folder, output_folder_name, params, nthreads):

    exec_command = f'../exec {params["nr"]} {params["nth"]} {params["nph"]} {params["cfl"]}'
    dir_name = f'{output_folder_name}_{params["nr"]}_{params["nth"]}_{params["nph"]}'

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
mkdir output/{dir_name} ; cd output/{dir_name} ; mkdir checkpoint_data
export OMP_NUM_THREADS={nthreads}
{exec_command}
        ''')


nrs = [300, 400]
nths = [8]
nphs = [16]
cfl = 0.5

job_name = 'gaussian'
code_folder = 'Gaussian/FD_Order8'
output_folder_name = 'fields'
nthreads = 24

sbatch_folder = Path.cwd().parent / 'nrpytutorial' / \
    'BSSN_LinearScalarFieldEvolution_Ccodes' / 'sbatch'
if not sbatch_folder.exists():
    sbatch_folder.mkdir()

count = 0
for nr in nrs:
    for nth in nths:
        for nph in nphs:

            params = {}
            params['nr'] = nr
            params['nth'] = nth
            params['nph'] = nph
            params['cfl'] = cfl

            file_path = sbatch_folder / f'job_{count}.sbatch'
            make_sbatch_file(file_path, job_name, code_folder, output_folder_name, params, nthreads)

            count += 1

for i in range(count):
    print(f'sbatch sbatch/job_{i}.sbatch')