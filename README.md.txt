## Purpose
Scoring pipeline to prioritize putative Type III Secretion System (T3SS) effectors by combining SWISS-MODEL, AlphaFold, HHpred, genomic context, and predicted interactions.

## Quick usage
# install packages (one-time)
Rscript packages.R

# create a sample input 
Rscript data_sample_input_generation.R

# run scoring
Rscript run_scoring_combined.R path/to/input.xlsx path/Rscript run_scoring_combined.R /mnt/c/Users/Reviewer/Downloads/test2.xlsx /home/reviewer/results/scoring_combined.xlsxto/output.xlsx