This folder contains the following files/directories:

- 0.setup.R: This script loads the required packages and functions for the simulation code.
- 1.structure_learning_w_SS.R: This script runs TreeMHN and genotype MHN of Schill et al. (2020) to do structure learning with and without stability selection.
- 2.pathKL_runtime_memory.R: This script gets the pathway KL divergence, runtime, and memory for TreeMHN (MLE/MCEM), genotype MHN, REVOLVER, and HINTRA.
- 3.plot_script.R: This script contains code to generate plots.
- structure_learning_w_SS: Sample output from 1.structure_learning_w_SS.R.
- pathKL_runtime_memory: Sample output from 2.pathKL_runtime_memory.R.

To run the simulation code, one could for example type the following in the command line:

	Rscript 1.structure_learning_w_SS.R -n 10 -N 500 -i 1 --ncores=125

	Rscript 2.pathKL_runtime_memory.R -n 10 -N 200 --method=TreeMHN_MLE

We highly recommend running the scripts using many nodes in a cluster, because the simulation runs and the stability selection procedure can be parallelized.

