# Methylr2py
Python wrapped R packages for methylation array analysis.

This repo contains code to perform methylation array analysis. 

We have taken some freely available R packages as
'minfi'
'ChAMP'
'DNAmArray' and some other packages and wrapped them with r2py in order to get them connected to Python.

One approach we adopted for this is PyMethylProcess and MethylNet which have already taken that path. We adjusted many features and put some new features into our pipeline.
You can see an example notebook on how to use this code and we hope it helps in analysing MethylationData in the context especially for Python Programmers.

One major advantage is that it (Python programming language) facilitates the use of Machine and Deep Learning (as also shown in MehtylNet) and is also shown in an extra example notebook.

All software needed to be installed is shown in our "https://github.com/FAU-DLM/GPU-Jupyterhub" as we do all analysis via a docker hosted jupyterhub deep learning platform.


Root:.

├───Script (Ipython Notebook i.e.)

│    └───IDAT_Folder

│        └───Subfolder (*) 

│            └───Processed Phenosheet

└─── Samplesheet       

* It is important that there is no other pheno- or sample-sheet in the path (recursively also!!) specified where the module looks up the phenosheet!!!
