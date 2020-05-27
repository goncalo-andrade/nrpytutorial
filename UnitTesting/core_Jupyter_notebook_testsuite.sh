#!/bin/bash

set -e # Error out if any commands complete with an error.

# Skip Baikal and all Start-to-Finish notebooks, except ScalarWave for now
# Tutorial-Start_to_Finish-ScalarWave*.ipynb
# Let's try all but Psi4 Start-to-Finish and Baikal notebooks
for i in Tutorial-[A]*.ipynb Tutorial-[C-RT-Z]*.ipynb Tutorial-B[B-Z]*.ipynb Tutorial-S[A-SU-Z]*.ipynb Tutorial-Start_to_Finish-*[^4].ipynb NRPyPN/PN*.ipynb; do
    ./run_Jupyter_notebook.sh $i notimer
    cat $i | sed "s/\\\r\\\n/\\\n/g" > $i-new && mv $i-new $i
    git diff $i |grep -v "image/png"|grep -E "^\-|^\+"|grep -v  '^\-\-\-'|cdiff |cat
#    git diff $i | grep -v "image/png" | cdiff | cat
#    echo Number of lines different in the git diff: `git diff|grep -v image/png|wc -l`
done
# ./run_Jupyter_notebook.sh Tutorial-Finite_Difference_Derivatives.ipynb && git diff Tutorial-Finite_Difference_Derivatives.ipynb && \
# ./run_Jupyter_notebook.sh Tutorial-Numerical_Grids.ipynb && git diff Tutorial-Numerical_Grids.ipynb && \
# ./run_Jupyter_notebook.sh Tutorial-Coutput__Parameter_Interface.ipynb && git diff Tutorial-Coutput__Parameter_Interface.ipynb && \
# ./run_Jupyter_notebook.sh Tutorial-cmdline_helper.ipynb && git diff Tutorial-cmdline_helper.ipynb && \
# ./run_Jupyter_notebook.sh Tutorial-Symbolic_Tensor_Rotation.ipynb && git diff Tutorial-Symbolic_Tensor_Rotation.ipynb
