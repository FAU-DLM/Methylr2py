# Methylr2py
Python wrapped R packages for methylation array analysis.

This repo contains code to perform methylation array analysis as used for the following publications:


Kobow, K., Jabari, S., Pieper, T. et al. Mosaic trisomy of chromosome 1q in human brain tissue associates with unilateral polymicrogyria, very early-onset focal epilepsy, and severe developmental delay. Acta Neuropathol 140, 881–891 (2020). 
https://doi.org/10.1007/s00401-020-02228-5

Jabari, S., Kobow, K., Pieper, T. et al. DNA methylation-based classification of malformations of cortical development in the human brain. Acta Neuropathol 143, 93–104 (2022). 
https://doi.org/10.1007/s00401-021-02386-0

Hoffmann, L., Coras, R., Kobow, K. et al. Ganglioglioma with adverse clinical outcome and atypical histopathological features were defined by alterations in PTPN11/KRAS/NF1 and other RAS-/MAP-Kinase pathway genes. Acta Neuropathol 145, 815–827 (2023). 
https://doi.org/10.1007/s00401-023-02561-5
      
       

 We have taken some freely available R packages such as
'minfi'
'ChAMP'
'DNAmArray' and some other packages and wrapped them with r2py in order to get them connected to Python.

One approach we adopted for this is PyMethylProcess and MethylNet which have already taken that path. We adjusted many features and put some new features into our pipeline.
You can see an example notebook on how to use this code and we hope it helps in analysing MethylationData in the context especially for Python Programmers.

One major advantage is that it (Python programming language) facilitates the use of Machine and Deep Learning (as also shown in MehtylNet) and is also shown in an extra example notebook.

All software needed to be installed is shown in our "https://github.com/FAU-DLM/GPU-Jupyterhub" as we do all analysis via a docker hosted jupyterhub deep learning platform.

### *It is important that there is no other pheno- or sample-sheet in the path (recursively also!!) specified where the module looks up the phenosheet!!!

    .
    ├── Script                   # (Ipython Notebook i.e.)
    │   └───IDAT_Folder
    │       └───Subfolder (*) 
    │           └───Processed Phenosheet
    └─── Samplesheet       
