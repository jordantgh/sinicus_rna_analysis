# Rhinolophus sinicus RNA sequencing analysis pipelines

- reference_pipeline.sh gets the rna seq and wgs data, and does genome indexing and alignment

- denovo.sh is work in progress to get the rna seq data aligned against a reference de novo transcriptome rather than a reference genome

- The R_files directory has a (very) basic plotting function, can be improved
- The tabulate_tx.py script is used in denovo.sh to create a .csv with 
transcript sequences matching our target sequence (currently using the Sinicus IFITM3-like variant x2)

- I use poetry and pyenv for specifying the local python environment (hence the pyproject.toml and .python-version files). I am not sure what options are there on the Eddie cluster for this purpose, maybe you can't use them there. However, I think you may be able to run denovo.sh locally as it is much less demanding than reference_pipeline.sh. If this is something you want to try I can provide instructions about poetry and pyenv

- I will be running denovo.sh for *all* the sra files soon and will provide the .csvs once this is completed.