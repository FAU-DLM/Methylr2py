import pandas as pd
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr 
from rpy2.robjects import r
import rpy2.robjects.packages as rpackages
from pymethylprocess.MethylationDataTypes import *
from rpy2.robjects import pandas2ri, numpy2ri
import pickle
import sqlite3
import os, glob, subprocess
from os.path import basename
from collections import Counter
import seaborn as sns
import matplotlib.pyplot as plt 
from bokeh.models import ColumnDataSource
from bokeh.plotting import figure, show
from bokeh.palettes import d3
import bokeh.models as bmo
from bokeh.io import output_notebook
from bokeh.layouts import gridplot
from bokeh.transform import factor_cmap
from bokeh.models import NumeralTickFormatter
import methylcheck
from scipy.special import comb
PATH="/home/Deep_Learner/private/third_party_repos/Methylr2py/manifest/EPIC.hg19.manifest.tsv.gz" ###--->can be downloaded here :"http://zwdzwd.io/InfiniumAnnotation/current/EPIC/"
pandas2ri.activate()
numpy2ri.activate()


class PreProcessIDATs:
    '''Class that will preprocess IDATs using R pipelines.
    Many things have been adapted from the original file-->"https://github.com/Christensen-Lab-Dartmouth/PyMethylProcess/blob/master/pymethylprocess/PreProcessDataTypes.py"
    
    
    idat_dir
        Location of idats or samplesheet csv.
    minfi
        Rpy2 importr minfi library, default to None will load through rpy2
    enmix
        Rpy2 importr enmix library, default to None will load through rpy2
    base
        Rpy2 importr base library, default to None will load through rpy2
    meffil
        Rpy2 importr meffil library, default to None will load through rpy2
        '''
    
    
    def __init__(self, idat_dir, minfi=None, enmix=None, base=None, meffil=None):
        self.idat_dir = idat_dir
        if minfi == None:
            self.minfi = importr('minfi')
        else:
            self.minfi = minfi
        if enmix == None:
            self.enmix = importr("ENmix")
        else:
            self.enmix = enmix
        if base == None:
            self.base = importr('base')
        else:
            self.base = base
        try:
            if meffil==None:
                self.meffil = importr('meffil')
            else:
                self.meffil = meffil
        except:
            self.meffil=None
        self.qcinfo=robjects.r('NULL')
        
        
        
    def load_idats(self,savedir=None, phenosheetdir=None, use_cache=False, rename_samples=True, parallel=True, nworkers=2, verbose=True, extended=True):
        """For minfi pipeline, load IDATs from specified idat_dir."""
        
        if savedir is not None:
            cache_storage_path = os.path.join(savedir,'RGSet.rds')
        else:
            cache_storage_path = os.path.join(self.idat_dir,'RGSet.rds')
        
        if phenosheetdir is not None:           
            self.pheno = self.minfi.read_metharray_sheet(phenosheetdir)
        else:
            self.pheno = self.minfi.read_metharray_sheet(self.idat_dir) 
        
        if use_cache:
            self.RGset=robjects.r('readRDS')(cache_storage_path)            
        else:
            if parallel:
                self.RGset=self.__load_idats_parallel(targets=self.pheno,verbose=verbose,extended=extended,nworkers=nworkers)
                
            else:    
                self.RGset = self.minfi.read_metharray_exp(targets=self.pheno, extended=extended)
                
        
        if rename_samples:
            self.rename_samples()            
        
        if not os.path.exists(cache_storage_path):
            robjects.r('saveRDS')(self.RGset,cache_storage_path)  
            
        if not use_cache:
            robjects.r('saveRDS')(self.RGset,cache_storage_path) 
        self.RGset_orig=self.RGset    
        self.pheno_orig=self.pheno
        self.pheno_orig_py=self.ri2py_dataframe(self.pheno_orig, matrix=True)
            
    
    def __load_idats_parallel(self, targets, verbose=True, extended=True, nworkers=2 ):       
        targets=targets
        verbose=verbose
        nworkers=nworkers
        extended=extended
        RGset= robjects.r("""function (targets, verbose, extended, nworkers) {
            read.metharray.exp.par <- function(targets=targets, verbose = verbose, extended=extended, n_workers=nworkers, ...) {
                library(BiocParallel)
                nworkers <- bpworkers(bpparam())
                
                if (n_workers<=nworkers) {
                        nworkers=n_workers
                      }
                else {
                        nworkers=nworkers
                     }   
                param<-MulticoreParam(workers=nworkers)
                #cat(nworkers)               
                
                if (nworkers <= 1)
                    stop("Did you registered a biocparallel back-end?")
                y <- rep(1, ceiling(nrow(targets)/nworkers))
                for (i in 2:nworkers) y <- c(y, rep(i, ceiling(nrow(targets)/nworkers)))
                y <- y[1:nrow(targets)]
                jobs <- split(targets, y)

                fun <- function(x, ...) {
                    ##these need to be loaded on the worker nodes explicitly for BatchJobs!
                    requireNamespace("minfi")
                    requireNamespace("Biobase")
                    read.metharray.exp(targets = x, extended=extended)
                }

                if(verbose)cat("Reading multiple idat-files in parallel")
                res <- bplapply(jobs, FUN = fun, BPPARAM = param)
                combine <- minfi::combine
                #if(verbose)
                #    message(str(res))
                if(verbose)cat("Combining the RGsets to one big RGset")
                rgSet <- res[[1]]
                
                for (i in 2:length(res)) rgSet <- combine(rgSet, res[[i]])
                return(rgSet)
            }
        
            rgSet<-read.metharray.exp.par(targets=targets, verbose = verbose, extended=extended, n_workers=nworkers)
        
        
            
            return(rgSet)
            }""")(targets, verbose, extended, nworkers)
        return RGset
        
            
    
    def rename_samples(self):       
        
        self.RGset= robjects.r("""function (rgSet,pheno) {
            #pheno$ID <- paste(disease, sample, sep=".")
            sampleNames(rgSet) <- pheno$ID
            result=list(rg=rgSet)
            return(result$rg)
            }""")(self.RGset, self.pheno)        
        
            
    
    def getQC(self, addQC=False, phenotype=None, RGset=None):        
        
        RGset=RGset if RGset else self.RGset
        pheno= robjects.r("pData")(RGset)
        if phenotype is not None:            
            if phenotype not in pandas2ri.ri2py(robjects.r['as'](pheno,'data.frame')).columns.tolist():
                print('The pheno sheet does not contain a '+phenotype+' column you specified \n'
                         'These are the available column names:')
                print(pandas2ri.ri2py(robjects.r['as'](pheno,'data.frame')).columns.tolist())
                return None, None, None, None
        
        MSet = self.minfi.preprocessRaw(RGset)
        qc,cols,rows=robjects.r("""function ( MSet ) {
        qc = getQC(MSet)         
        
        result=list(qc,colnames(qc), rownames(qc))
        return(result)
             }""")(MSet)          
        
        data_py=self.ri2py_dataframe(qc, matrix=True)
        
        rows_py=pd.Series(pandas2ri.ri2py(rows))
        cols_py=pd.Series(pandas2ri.ri2py(cols))
        datfr=pd.DataFrame(data_py.to_numpy(),index=rows_py,columns=cols_py)
        if addQC:            
            pheno=pd.concat([pandas2ri.ri2py(robjects.r['as'](pheno,'data.frame')), data_py], axis=1, sort=False)
            self.pheno=pandas2ri.py2ri(pheno)
        
        if phenotype is not None:
            try:       
                datfr[phenotype]=pd.DataFrame(pandas2ri.ri2py(pheno))[phenotype].to_numpy()
            except:
                datfr[phenotype]=pd.DataFrame(pandas2ri.ri2py(robjects.r['as'](pheno,'data.frame'))[phenotype]).to_numpy()
        
        self.QC_df=datfr
        if phenotype is not None:
            return MSet, qc, datfr, phenotype
        else:
            return MSet, qc, datfr
                             
              
    
    def getM(self,Object=None):
        if Object:
            Object=Object
        else:
            Object=self.RGset
            
        mvals_py, mvals=self.compute_mvals(Object)    
        return mvals_py, mvals        
    
    
    def getBeta(self,Object=None):
        if Object:
            Object=Object
        else:
            Object=self.RGset
            
        betas_py,betas=self.compute_betas(Object) 
        return betas_py, betas       
    
    
    def compute_betas(self, Object):
        """Get beta value matrix from minfi after finding RSet."""       
                
        beta = self.minfi.getBeta(Object)
        beta_py = self.ri2py_dataframe(beta)
        
        return beta_py , beta   
    
    
    def compute_mvals(self, Object):
        """Get mvalue matrix from minfi after finding RSet."""    
        
        rgset = robjects.r("""function (obj) {             
             if(class(obj) == "RGChannelSet" | class(obj) == "RGChannelSetExtended"  ){                    
                    return(TRUE)  
                    }
                    else{
                    return(FALSE)
                    }             
            }""")(Object) 

        if pandas2ri.ri2py(rgset):
            mval = self.minfi.getM(self.minfi.preprocessRaw(Object))        
        
        else:
            mval = self.minfi.getM(Object)
            
        mval_py = self.ri2py_dataframe(mval)
        
        return mval_py, mval 
          
    def detectionP(self, RGset=None):
        if RGset is None:
            RGset = self.RGset
            
        # calculate the detection p-values
        detP =robjects.r("detectionP")(RGset) 
       
        detP_py=self.ri2py_dataframe(detP, matrix=False)
        return detP,detP_py   
    
    def ri2py_dataframe(self, r_dat, matrix=False):
        if matrix:
            py_dat=pd.DataFrame(pandas2ri.ri2py(robjects.r['as'](r_dat,'data.frame')))
        elif not matrix:
            py_dat=pd.DataFrame(pandas2ri.ri2py(r_dat),
                index=numpy2ri.ri2py(robjects.r("rownames")(r_dat)),
                columns=numpy2ri.ri2py(robjects.r("colnames")(r_dat)))
        return py_dat  
    
     
    
    def remove_badprobes(self, obj=None, RGset=None, detP=None, ProbeCutoff=0, SampleCutoff=0.1, detPcut=0.01, verbose=True):        
        
        if obj:
            obj=obj
        else:
            obj=self.GRset
            
        RGset=RGset if RGset else self.RGset    
        grset = robjects.r("""function (obj) {             
            if(class(obj) == "GenomicRatioSet")
                    return(TRUE)  
             else return(FALSE)                
            }""")(obj)   
        
        if not pandas2ri.ri2py(grset):          
            if (verbose):    
                print('sorry obj needs to be of type "GenomicRatioSet"...')             
            return None
        
        if detP:
            detP=detP  
        else:
            detP,_=self.detectionP(RGset=RGset)
        
        
        GRset = robjects.r("""function (obj, cutprobes, detP, detPcut, cutsamples, verbose)
        
        {
               
               if(identical(rownames(detP),rownames(obj)) & identical(colnames(detP),colnames(obj)))
                    {
                        if(verbose)cat("    detP check success.")
                    }
                    else {
                        
                        if(verbose)cat("Your detP matrix has been aligned to match the EXACT same rowname and colname as Data Matrix.")
                        
                        keep <- ( rownames(detP) %in% rownames(obj))
                        detP <- detP[keep,]  
                        keep <- (colnames(detP) %in% colnames(obj))
                        detP <- detP[,keep]
                        
                    }  
               
                numfail <- matrix(colMeans(detP >= detPcut))
                rownames(numfail) <- colnames(detP)
                colnames(numfail) <- "Failed CpG Fraction."

                RemainSample <- which(numfail <= cutsamples)
               
                if(any(numfail > cutsamples))
                {
                    cat("\n    The detSamplecut parameter is : ",cutsamples)
                    cat("    Samples : ", paste(rownames(numfail)[which(numfail >= cutsamples)],collapse=",")," will be deleted.")
                    cat("    There are ",length(RemainSample)," samples remained for analysis.")
                    obj <- obj[,RemainSample]

                }
               
               
               RemainProbe <- rowSums(detP > detPcut) <= cutprobes * length(RemainSample)
               
               if(cutprobes==0)
               
                   {   
                       if(verbose)cat("\n    Filtering probes with a detection p-value above ",detPcut,".")
                       if(verbose)cat("    Removing ",sum(RemainProbe==FALSE)," probes.")
                       if(verbose)cat("    If a large number of probes have been removed, we suggests you to identify potentially bad samples")
                   } else {
                       if(verbose)cat("\n    Filtering probes with a detection p-value above ",detPcut," in at least", cutprobes*100,"% Samples.")
                       if(verbose)cat("    Removing ",sum(RemainProbe==FALSE)," probes.")
                       if(verbose)cat("    If a large number of probes have been removed, we suggests you to identify potentially bad samples")
                   }
              
               obj <-  obj[RemainProbe,]
               
               if(sum(detP > detPcut) > 0) cat("    There are still ",sum(detP > detPcut), " failed probes exist in your data set, imputation is recommended.")                
                return(obj)
                
                 }""")(obj, ProbeCutoff, detP, detPcut, SampleCutoff,verbose ) 
        pheno= robjects.r("pData")(GRset)
        return GRset, pheno
    
    
    def remove_badsamples(self, badSampleCutoff=10, rm_badsamples=True, detPFilter=False, detPcut=0.01, SampleCutoff=0.1, addQC=False, verbose=True, RGset=None):
        ''' 
        excluding samples if avarage of unmeth and meth values is below 10
        and if detectionP probe wise threshold is below 0,01(detPcut) 
        and fraction of affected probes for a given sample is above 0.1 (cutsamples)
             
        '''
        
        if RGset:
            RGset=RGset
        else:
            RGset=self.RGset
            
        MSet,qc, datfr = self.getQC(addQC=addQC, phenotype=None, RGset=RGset)
        #print(qc)    
        RGset, pheno = robjects.r("""function (rgset, badSampleCutoff, detPcut, cutsamples, verbose, rm_badsamples, detPFilter) {
            
            #### meth and unmeth filtering ##########
            MSet<-preprocessRaw(rgset)
            qc <- getQC(MSet)
            if (rm_badsamples)
            {
                meds <- (qc$uMed + qc$mMed) / 2            
                keepIndex <- which(meds > badSampleCutoff)            
                if (length(keepIndex) == 0) {
                    stop("All samples found to be bad")
                }
                if (length(keepIndex) < ncol(rgset)) {               
                    if(verbose){
                        cat(
                            sprintf("Found and removed %s bad samples",
                                    ncol(rgset) - length(keepIndex)))
                    }
                    rgset <- rgset[, keepIndex]
                }
            }
            
            ####  detP filtering  ########                
            if(detPFilter)
            {
                if(verbose){   
                    cat("\n  Filtering Detect P value Start\n")

                }

                detP <- detectionP(rgset)
                numfail <- matrix(colMeans(detP >= detPcut))
                rownames(numfail) <- colnames(detP)
                colnames(numfail) <- "Failed CpG Fraction."

                RemainSample <- which(numfail <= cutsamples)
                print(numfail)
                if(any(numfail > cutsamples))
                {
                    cat("\n    The detSamplecut parameter is : ",cutsamples)
                    cat("    Samples : ", paste(rownames(numfail)[which(numfail >= cutsamples)],collapse=",")," will be deleted.")
                    cat("    There are ",length(RemainSample)," samples remained for analysis.")
                    rgset <- rgset[,RemainSample]

                } 

                else{
                    cat("\ There are no failed samples in your datset")
                }

               if(sum(detP > detPcut) > 0) cat("    There are still ",sum(detP > detPcut), " failed probes exist in your data set, imputation is recommended.") 
            }
           
           if('filenames' %in% names(colData(rgset)))
               {
               colData(rgset) <- subset(colData(rgset), select = -filenames )
               }
           pheno = pData(rgset)
           result=list(rgset, pheno)
           return(result)                       
            
            }""")(RGset, badSampleCutoff, detPcut, SampleCutoff, verbose, rm_badsamples, detPFilter) 
        
        return RGset, pheno
               
        
    #### probe wise QC filtering ######
    
    def probeFiltering(self, RGset=None, cutbead=3, zeropoint=True, verbose=True):
        
        if RGset:
            RGset=RGset
        else:
            RGset=self.RGset

        self.RGset_filt, pheno = robjects.r("""function (RGset, cutbead, zeropoint, verbose) { 


            ##' filter on probe-level: number of beads and zero intensities
            ##'
            ##'
            ##' @title filter on probe-level: number of beads and zero intensities
            ##' @param RGset a RGChannelSetExtended
            ##' @param cutbead threshold number of beads
            ##' @param zeroint filter out zero intenseties default is TRUE
            ##' @param verbose Default is TRUE
            ##' @return RGChannelSet
            ##' @author mvaniterson
            ##' @export
            ##' @import minfi
            ##' @importFrom Biobase varMetadata AnnotatedDataFrame phenoData featureData experimentData annotation protocolData assayDataElement
            probeFiltering <- function(RGset, cutbead, zeroint, verbose){

                if(class(RGset) != "RGChannelSetExtended")
                    stop("RGset should be of class 'RGChannelSetExtended' in order to perform filtering on number of beads!")

                ##Filter on number of beads
                if(verbose)
                    cat("Filtering on number of beads... \n")

                beadmat <- getNBeads(RGset)
                #RedSD <- RGset$RedSD
                #GreenSD <- RGset$GreenSD

                idBeadmat <- beadmat < cutbead
                ##beadmat[idBeadmat] <- NA
                
                
                Grn <- getGreen(RGset)
                Red <- getRed(RGset)
                
                if(verbose)
                    cat("On average", round(100*sum(idBeadmat)/prod(dim(idBeadmat)), 2),"% of the probes (",nrow(idBeadmat),") have number of beads below", cutbead, "\n")

                ##Filter on Red and Green intensity <1
                if(zeroint) {
                    if(verbose)
                        cat("Filtering on zero intensities... \n")

                    

                    ##determine if Grn and/or Red intensities of type II probes are <1
                    idT2 <- Grn[getProbeInfo(RGset, type = "II")$AddressA,] < 1 | Red[getProbeInfo(RGset, type = "II")$AddressA,] < 1

                    ##determine if either Grn or Red intensities of Type I probes are <1
                    idT1Grn <- Grn[c(getProbeInfo(RGset, type = "I-Green")$AddressA,
                                     getProbeInfo(RGset, type = "I-Green")$AddressB),] < 1

                    idT1Red <- Red[c(getProbeInfo(RGset, type = "I-Red")$AddressA,
                                     getProbeInfo(RGset, type = "I-Red")$AddressB),] < 1

                    if(verbose) {
                        cat("On average", round(100*sum(idT2)/prod(dim(idT2)), 3),"% of the Type II probes (",nrow(idT2),") have Red and/or Green intensity below 1 \n")
                        cat("On average", round(100*sum(idT1Grn)/prod(dim(idT1Grn)), 3),"% of the Type I probes (",nrow(idT1Grn),"), measured in Green channel, have intensity below 1 \n")
                        cat("On average", round(100*sum(idT1Red)/prod(dim(idT1Red)), 3),"% of the Type I probes (",nrow(idT1Red),"), measured in Red channel, have intensity below 1 \n")
                    }
                }

                ##combine all filtered results and set NA in Red and/or Green channels
                Red[idBeadmat] <- Grn[idBeadmat] <- NA

                if(zeroint) {
                    if(verbose){
                        cat("Set filtered probes in Red and/or Green channels to NA... \n")
                    }

                    for(i in 1:ncol(RGset)) {
                        if(verbose & i%%100 == 0)
                            cat("... done ",i," out of ",ncol(RGset)," ... \n")
                        idRed <- c(names(which(idT2[,i])), names(which(idT1Red[,i])))
                        midRed <- match(idRed, rownames(Red))
                        Red[midRed, i] <- NA
                        idGrn <- c(names(which(idT2[,i])), names(which(idT1Grn[,i])))
                        midGrn <- match(idGrn, rownames(Grn))
                        Grn[midGrn, i] <- NA
                    }
                }
                #RedSD <- RedSD[idBeadmat]
                #GreenSD <- GreenSD[idBeadmat]

                RGset<-RGChannelSet(Green = Grn, Red = Red,
                             colData = colData(RGset),
                             annotation = annotation(RGset)                        
                             )

                

            }


            RGset<-probeFiltering(RGset=RGset, cutbead=cutbead, zeroint=zeropoint, verbose=verbose)
            
            pheno = pData(RGset)
            result=list(RGset, pheno)
            return(result)  
            



        }""")(RGset, cutbead, zeropoint, verbose)
    
    
        return self.RGset_filt, pheno
        
    
    def preprocessFunnorm(self, celltype_adoption=False, use_cell_count2=False, nPCs=2, RGset=None):   
            RGset=RGset if RGset else self.RGset
            
            self.GRset = robjects.r("""function (RGset, nPCs) {
                 mSetSq <- preprocessFunnorm(RGset, nPCs=nPCs, sex = NULL, bgCorr = TRUE, dyeCorr = TRUE, verbose = TRUE)
                return(mSetSq)
                }""")(RGset, nPCs)
            
            if celltype_adoption:
                if use_cell_count2:
                    self.NeunC,self.NeunC_table,_  = self.est_cell_counts2(RGset=RGset, processMethod="preprocessFunnorm", nPCs=nPCs)

                else:    
                    self.NeunC,self.NeunC_table, _ = self.est_cell_counts(RGset=RGset, processMethod="preprocessFunnorm", nPCs=nPCs)    
            
                self.insert_cell_types()
            pheno=robjects.r("pData")(self.GRset)
            return self.GRset, pheno


    def preprocessQuantile(self,celltype_adoption=False, use_cell_count2=False,nPCs=2,  RGset=None):   
            RGset=RGset if RGset else self.RGset
            
            self.GRset = robjects.r("""function (RGset) {
                    mSetSq <- preprocessQuantile(RGset, fixOutliers = TRUE, removeBadSamples = FALSE,
                               badSampleCutoff = 10.5, quantileNormalize = TRUE,
                               stratified = TRUE, mergeManifest = FALSE, sex = NULL,
                               verbose = TRUE)
                    return(mSetSq)
                    }""")(RGset)
            if celltype_adoption:
                if use_cell_count2:
                    self.NeunC,self.NeunC_table, _ = self.est_cell_counts2(RGset=RGset,processMethod="preprocessQuantile", nPCs=nPCs)

                else:    
                    self.NeunC,self.NeunC_table, _ = self.est_cell_counts(RGset=RGset,processMethod="preprocessQuantile", nPCs=nPCs)       
            
                self.insert_cell_types()
            pheno=robjects.r("pData")(self.GRset)
            return self.GRset, pheno   
    
    def preprocessNoob(self,celltype_adoption=False, use_cell_count2=False, nPCs=2, RGset=None):
        RGset=RGset if RGset else self.RGset
        
        self.GRset = robjects.r("""function (RGset) {
           mSet <- preprocessNoob(rgSet=RGset,offset = 15, dyeCorr = TRUE, verbose = TRUE,
               dyeMethod=c("single")) 
           GMset <- mapToGenome(mSet)
           GRset<-ratioConvert(GMset)
           
           return(GRset)
            }""")(RGset)
        if celltype_adoption:
            if use_cell_count2:
                self.NeunC,self.NeunC_table, _ = self.est_cell_counts2(RGset=RGset,processMethod="preprocessNoob", nPCs=nPCs)

            else:    
                self.NeunC,self.NeunC_table, _ = self.est_cell_counts(RGset=RGset,processMethod="preprocessNoob", nPCs=nPCs)       
            
            self.insert_cell_types()
        pheno=robjects.r("pData")(self.GRset)
        return self.GRset, pheno 

    def est_cell_counts(self,RGset=None,**kwargs):
            
            processMethod=list(kwargs.keys())[0]     
            processMethodname=list(kwargs.values())[0]
            PCs=list(kwargs.keys())[1]     
            nPCs=list(kwargs.values())[1]   
            rgset=RGset if RGset else self.RGset
            robjects.r('library(FlowSorted.DLPFC.450k)')
            cell_count_estimates = robjects.r("""function (RGset, PCs, nPCs, processMethodname) {
                cellCounts <- estimateCellCounts(RGset, compositeCellType = "DLPFC",
                           processMethod = processMethodname, probeSelect = "both",
                           cellTypes = c("NeuN_neg", "NeuN_pos"),
                           referencePlatform = c("IlluminaHumanMethylation450k",
                                                 "IlluminaHumanMethylationEPIC",
                                                 "IlluminaHumanMethylation27k"),
                           returnAll = TRUE, meanPlot = FALSE, verbose = TRUE, PCs=nPCs)

                return(cellCounts)
                }""")(rgset, PCs, nPCs, processMethodname)        
            return cell_count_estimates

    def est_cell_counts2(self,RGset=None, **kwargs):
            """Given RGSet object, estimate cell counts using reference approach via FlowSorted.Blood.EPIC 
            estimatecellcounts2 method.

            Parameters
            ----------
            rgset
                RGSet object stored in python via rpy2

           """
            rgset=RGset if RGset else self.RGset
            

            processMethod=list(kwargs.keys())[0]     
            processMethodname=list(kwargs.values())[0]
            PCs=list(kwargs.keys())[1]     
            nPCs=list(kwargs.values())[1]        
            robjects.r('library(FlowSorted.Blood.EPIC)')
            cell_count_estimates = robjects.r("""function (RGset, PCs, nPCs, processMethodname) {
                cellCounts <- estimateCellCounts2(RGset, compositeCellType = "DLPFC",
                           processMethod = processMethodname, probeSelect = "both",
                           cellTypes = c("NeuN_neg", "NeuN_pos"),
                           referencePlatform = c("IlluminaHumanMethylation450k",
                                                 "IlluminaHumanMethylationEPIC",
                                                 "IlluminaHumanMethylation27k"),
                           referenceset = NULL, IDOLOptimizedCpGs = NULL, returnAll = TRUE,
                  meanPlot = FALSE, verbose = TRUE, lessThanOne = FALSE, fixOutliers=TRUE, removeBadSamples = FALSE, PCs=nPCs )

                return(cellCounts)
                }""")(rgset, PCs, nPCs, processMethodname)        
            return cell_count_estimates    


    def insert_cell_types(self):
            celltypes=pd.DataFrame(pandas2ri.ri2py(self.NeunC), columns=['Glia','Neurons'])
            #pheno=pd.concat([pandas2ri.ri2py(robjects.r['as'](self.pheno,'data.frame')), celltypes], axis=1, sort=False)
            celltypes_ri=pandas2ri.py2ri(celltypes)  
            #self.pheno=pandas2ri.py2ri(pheno)    
            self.RGset,self.GRset, self.pheno=robjects.r("""function (RGset,GRset, celltypes) {        
                .pDataAdd <- function(object, df) {
                    stopifnot(is(df, "data.frame") || is(df, "DataFrame"))
                    pD <- colData(object)
                    if (any(names(df) %in% names(pD))) {
                        alreadyPresent <- intersect(names(df), names(pD))
                        warning(sprintf(
                            "replacing the following columns in colData(object): %s",
                            paste(alreadyPresent, collapse = ", ")))
                        pD[, alreadyPresent] <- df[, alreadyPresent]
                        df <- df[, !names(df) %in% alreadyPresent]
                    }
                    if (ncol(df) > 0) {
                        # NOTE: Work around for bug in cbind(DataFrame, DataFrame)
                        rownam <- rownames(pD)
                        pD <- cbind(pD, df)
                        rownames(pD) <- rownam
                    }
                    colData(object) <- pD
                     return(object)
                }
                RGset<-.pDataAdd(RGset, celltypes)
                GRset<-.pDataAdd(GRset, celltypes)
                pheno<- pData(GRset)
                result=list(RGset, GRset, pheno)
                return(result)
                }""")(self.RGset, self.GRset, robjects.r.DataFrame(celltypes_ri)) 

        
        
    
    def DNAmArray_processing(self, GRset=None, RGset=None, filterXY=True, filterNoCG=True, excludeXreactiveprobes=True, dropSnPs=True, mask_probes=True, cutbead=3, zeropoint=True, what="both", ProbeCutoff=0.05, SampleCutoff=0.05, array_type='EPIC', badSampleCutoff=10, rm_badsamples=False, rm_badprobes=False, detPFilter=False, detPcut=0.01, addQC=False, verbose=True, autoimpute=True, imputation_method="imputePCA"):
        
        
               
        if GRset:
            GRset=GRset
        else:
            GRset = self.GRset
        if RGset:
            RGset=RGset
        else:
            RGset = self.RGset 
        
        if imputation_method!="methyLImp" and imputation_method!="imputePCA" and imputation_method!="knn":
            print('You did not specify a valid imputation method!!\n Choose "imputePCA" or "methyLImp" or "knn"...\n now using fallbackmethod "imputePCA"')
            imputation_method="imputePCA"   
        
        if detPFilter or rm_badsamples:
            if verbose:
                print('\n Now performing badsample removal')
            RGset, pheno = self.remove_badsamples(badSampleCutoff=badSampleCutoff, rm_badsamples=rm_badsamples, detPFilter=detPFilter, detPcut=detPcut, SampleCutoff=SampleCutoff, addQC=addQC, verbose=verbose, RGset=RGset)
        
        if verbose:
            print('\n Now performing probefiltering on beadcount')
        RGset_filt, pheno = self.probeFiltering(cutbead=cutbead, RGset=RGset,zeropoint=zeropoint, verbose=verbose) 
                 
        if rm_badprobes:
            if verbose:
                    print('\n Now removing bad probes ')  
            detP, _ = self.detectionP(RGset=RGset)        
            GRset, pheno=self.remove_badprobes(obj=GRset, RGset=RGset, detP=detP, ProbeCutoff=ProbeCutoff,SampleCutoff=SampleCutoff, detPcut=detPcut, verbose=True)      
            
            
        if rm_badprobes or rm_badsamples:
            if verbose:
                    print('\n Aligning RGset and GRset ') 
            GRset,self.pheno = robjects.r("""function (RGset, GRset) {  
                      #colnames(GRset) %in% colnames(RGset)
                      GRset <- GRset[,colnames(GRset) %in% colnames(RGset)]  
                      pheno<-pData(GRset)
                      result=list(GRset,pheno)
                      return(result)
             }""")(RGset_filt,GRset)  
            
        
        if verbose:
                print('\n Now performing reduce function')  
                
        if what=='M':            
            self.mval_py, self.pheno_py=self.reduce(GRset=GRset, RGset=RGset_filt, what=what, detPcut=detPcut, SampleCutoff=SampleCutoff, ProbeCutoff=ProbeCutoff, verbose=verbose,autoimpute=autoimpute, imputation_method=imputation_method)
            
            if filterNoCG or excludeXreactiveprobes or dropSnPs or filterXY or mask_probes is not False: 
                if verbose:
                    print('\n Now removing specific probes ')
                self.mval, self.pheno=self.filterCpGs(obj=self.mval , 
                                    GRset=GRset if GRset else None,         
                                    dropSnPs=dropSnPs,                                     
                                    filterXY=filterXY, 
                                    filterNoCG=filterNoCG, 
                                    excludeXreactiveprobes=excludeXreactiveprobes,
                                    mask_probes=mask_probes,                   
                                    array_type=array_type, 
                                    verbose=verbose)
            
            
                try:
                    self.pheno_py=self.ri2py_dataframe(self.pheno, matrix=True)
                except:
                    self.pheno_py=self.ri2py_dataframe(self.pheno, matrix=False)  
        
        
                #self.mval=obj[0]
                self.mval_py = self.ri2py_dataframe(r_dat=self.mval, matrix=False)
                #self.beta=obj[1]
                #self.beta_py = self.ri2py_dataframe(r_dat=self.beta, matrix=False)     

            
            return self.mval_py,self.pheno_py      
        
        
        
        if what=='beta':    
            
            self.beta, self.pheno_py=self.reduce(GRset=GRset, RGset=RGset_filt, what=what, detPcut=detPcut, SampleCutoff=SampleCutoff, ProbeCutoff=ProbeCutoff, verbose=verbose, autoimpute=autoimpute, imputation_method=imputation_method)
            
            if filterNoCG or excludeXreactiveprobes or dropSnPs or filterXY or mask_probes is not False: 
                if verbose:
                    print('\n Now removing specific probes ')
                self.beta, self.pheno=self.filterCpGs(obj=self.beta , 
                                    GRset=GRset if GRset else None,         
                                    dropSnPs=dropSnPs,                                     
                                    filterXY=filterXY, 
                                    filterNoCG=filterNoCG, 
                                    excludeXreactiveprobes=excludeXreactiveprobes, 
                                    mask_probes=mask_probes,                  
                                    array_type=array_type, 
                                    verbose=verbose)
            
            
                try:
                    self.pheno_py=self.ri2py_dataframe(self.pheno, matrix=True)
                except:
                    self.pheno_py=self.ri2py_dataframe(self.pheno, matrix=False)  
        
        
                #self.mval=obj[0]
                self.beta_py = self.ri2py_dataframe(r_dat=self.mval, matrix=False)
                #self.beta=obj[1]
                #self.beta_py = self.ri2py_dataframe(r_dat=self.beta, matrix=False)    
            
            
            return self.beta_py,self.pheno_py
        
        if what=='both':            
            self.beta_py, self.mval_py, self.pheno_py =self.reduce(GRset=GRset, RGset=RGset_filt, what=what, detPcut=detPcut, SampleCutoff=SampleCutoff, ProbeCutoff=ProbeCutoff, verbose=verbose,autoimpute=autoimpute, imputation_method=imputation_method) 
            
                
            if filterNoCG or excludeXreactiveprobes or dropSnPs or filterXY or mask_probes is not False: 
                if verbose:
                    print('\n Now removing specific probes for m-values')
                self.mval, self.pheno=self.filterCpGs(obj=self.mval , 
                                    GRset=GRset if GRset else None,         
                                    dropSnPs=dropSnPs,                                     
                                    filterXY=filterXY, 
                                    filterNoCG=filterNoCG, 
                                    excludeXreactiveprobes=excludeXreactiveprobes,
                                    mask_probes=mask_probes,                   
                                    array_type=array_type, 
                                    verbose=verbose)
            
            
                try:
                    self.pheno_py=self.ri2py_dataframe(self.pheno, matrix=True)
                except:
                    self.pheno_py=self.ri2py_dataframe(self.pheno, matrix=False)        
                #self.mval=obj[0]
                self.mval_py = self.ri2py_dataframe(r_dat=self.mval, matrix=False)
                if verbose:
                    print('\n Now removing specific probes for beta-values ')
                self.beta, _ =self.filterCpGs(obj=self.beta , 
                                    GRset=GRset if GRset else None,         
                                    dropSnPs=dropSnPs,                                     
                                    filterXY=filterXY, 
                                    filterNoCG=filterNoCG, 
                                    excludeXreactiveprobes=excludeXreactiveprobes,
                                    mask_probes=mask_probes,           
                                    array_type=array_type, 
                                    verbose=verbose)
                
                self.beta_py = self.ri2py_dataframe(r_dat=self.beta, matrix=False)    
            return self.beta_py, self.mval_py, self.pheno_py
              
        
    def mask_probes(self, obj, path=None):

        if path is None:
            path=PATH

        probelist=robjects.r("""function(path){                                      

            xReactiveProbes <- read.csv(path,
              sep="\t", header=T, stringsAsFactors=FALSE)  

            return(xReactiveProbes)

            }""")(path)

        probe_df=self.ri2py_dataframe(probelist, matrix=True)
        probe_pyarray=probe_df['probeID'][probe_df['MASK_general']==True].to_numpy()
        probe_rarray=numpy2ri.py2ri(probe_pyarray)

        obj=robjects.r("""function(obj, array){

            keep <- !(rownames(obj) %in% array)
            objkeep <- obj[keep,]

            #cat(nrow(obj))
            #cat(nrow(objkeep))

            return(objkeep)

            }""")(obj,probe_rarray)   
        
        return obj
        
        
    def dropLociWithSnps(self,GRset=None, obj=None):
        if GRset:
            GRset=GRset
        else:
            GRset=self.GRset
        if obj==None:
            obj=False
        result = robjects.r("""function (grset, obj) {             
             if(class(grset) != "GenomicRatioSet")
                    stop("The object should be of class 'GenomicRatioSet'!")
             #cat(obj)       
             #print(nrow(grset))
             
             grsetFlt <- dropLociWithSnps(grset)
             #print(nrow(grsetFlt))
             if (obj != FALSE)
             {
             #cat('here')
             snps<- grset[!(rownames(grset) %in% rownames(grsetFlt)),]
             #print(head(snps, 100))
             objflt<-obj[!(rownames(obj) %in% rownames(snps)),]     
             #print(nrow(obj))
             #print(nrow(objflt))
             result=list(grsetFlt, objflt)
             }
             else
             {
             result=list(grsetFlt)
             }
             
             return(result)
             
            }""")(GRset, obj)
        
       
        if len(pandas2ri.ri2py(result))>1:            
            GRset=result[0]
            obj=result[1]
            return GRset, obj
        else:
            GRset=result[0]
            return GRset, None
    
    def filterXY(self, obj=None):
        # if your data includes males and females, remove probes related to the sex chromosomes 
        self.get_annotation()
        obj = robjects.r("""function (obj,ann) {               
            keep <- !(rownames(obj) %in% ann$Name[ann$chr %in% c("chrX","chrY")])
            obj <- obj[keep,]    
            return(obj)
            }""")(obj, self.annotation)
        return obj
    
    def filterNoCG(self, obj=None, verbose=True):
        obj = robjects.r("""function (obj, verbose) {             
            if(verbose)cat("\n  Filtering NoCG Start")
            RemainProbe <- substr(rownames(obj),1,2) == "cg"
            if(verbose)cat("    Only Keep CpGs, removing ", sum(RemainProbe == FALSE) ," probes from the analysis.")
            obj <- obj[RemainProbe,]               
            return(obj)
            }""")(obj, verbose)
        return obj
    
    def excludeXreactiveprobes(self, obj=None, array_type='EPIC'):  
        
                    
        k450_criteria = ['Chen2013', 'Price2013', 'Naeem2014', 'DacaRoszak2015',
                    'Polymorphism', 'CrossHybridization', 'BaseColorChange', 'RepeatSequenceElements']
        EPIC_criteria = ['McCartney2016', 'Zhou2016', 'Polymorphism', 'CrossHybridization', 'BaseColorChange', 'RepeatSequenceElements']

        if array_type=='450k':
            #print('450k')
            for crit in k450_criteria:
                #print(crit, len(methylcheck.list_problem_probes('450k', [crit])))
                criteria= k450_criteria         
        elif array_type=='EPIC':            
            #print('EPIC')
            for crit in EPIC_criteria:
                #print(crit, len(methylcheck.list_problem_probes('EPIC', [crit])))
                criteria= EPIC_criteria          

        sketchy_probes_list = methylcheck.list_problem_probes(array=array_type, criteria=criteria)
            #df2 = methylcheck.exclude_probes(betaquant_5000, sketchy_probes_list)
        sketchy_probes_array=pd.DataFrame(np.array(sketchy_probes_list), columns=['TargetID'])    
        xReactiveProbes=pandas2ri.py2ri(sketchy_probes_array)

        obj = robjects.r("""function (obj, xReactiveProbes) {  
                    #xReactiveProbes <- read.csv(file=paste(dataDirectory,
                                                          # "48639-non-specific-probes-Illumina450k.csv",
                                                           #sep="/"), stringsAsFactors=FALSE)
                    keep <- !(rownames(obj) %in% xReactiveProbes$TargetID)
                    obj <- obj[keep,]                     
                    return(obj)
                        }""")(obj, xReactiveProbes)
        return obj  
    
    
    def filterCpGs(self, obj=None, GRset=None, dropSnPs=True, filterXY=True, filterNoCG=True, excludeXreactiveprobes=True, mask_probes=True, array_type='EPIC', verbose=True):
        orig_obj=obj        
        if GRset is None:
            GRset=self.GRset
        rgset = robjects.r("""function (obj) {             
             if(class(obj) == "RGChannelSet" | class(obj) == "RGChannelSetExtended"  ){                    
                    return(TRUE)  
                    }
                    else{
                    return(FALSE)
                    }             
            }""")(obj) 

        if pandas2ri.ri2py(rgset):
            print('Sorry this does not work with RGChannelSets, please use "GenomicRatioSet", "MethySets" , "RatioSets", "Mvals" or "betas"! ')
            return obj                 
        
                
        if excludeXreactiveprobes:
            if verbose:
                print('Dropping cross- reactive probes')
            obj = self.excludeXreactiveprobes(obj, array_type=array_type)      
        if filterXY:
            if verbose:
                print('Dropping XY-Chromosome-related probes')
            obj = self.filterXY(obj)    
        if filterNoCG:            
            obj = self.filterNoCG(obj=obj, verbose=verbose) 
            
        grset = robjects.r("""function (obj) {             
            if(class(obj) == "GenomicRatioSet")
                    return(TRUE)  
             else return(FALSE)                
            }""")(obj)   
        
        if dropSnPs and not pandas2ri.ri2py(grset):
            _,obj=self.dropLociWithSnps(GRset=GRset, obj=obj)
            
            
        elif dropSnPs and pandas2ri.ri2py(grset):
            obj,_=self.dropLociWithSnps(GRset=obj)
            #self.GRset=obj           
        
        if mask_probes:
            if verbose:
                print('\n Now performing mask_probes function')                    
            obj = self.mask_probes(obj, path=None)
        
        mset = robjects.r("""function (obj) {             
             if(class(obj) == "MethylSet")
                    return(TRUE)    
             else return(FALSE)                
            }""")(obj)
        #if pandas2ri.ri2py(mset):   
                  #self.Mset=obj
        
        rset = robjects.r("""function (obj) {             
             if(class(obj) == "RatioSet")
                    return(TRUE)                         
             else return(FALSE)
            }""")(obj)
        #if pandas2ri.ri2py(rset):   
                 # self.Rset=obj           
                  
        
        gmset = robjects.r("""function (obj) {             
             if(class(obj) == "GenomicMethylSet")
                    return(TRUE)                         
             else return(FALSE)
            }""")(obj)      
        #if pandas2ri.ri2py(gmset):   
                 # self.GMset=obj 
        if pandas2ri.ri2py(gmset) or pandas2ri.ri2py(rset) or pandas2ri.ri2py(mset) or pandas2ri.ri2py(grset):
            pheno=robjects.r("pData")(obj)
        
        else:
            
            pheno, GRset = robjects.r("""function (obj, GRset) {                   
                  
                  GRset <- GRset[,colnames(GRset) %in% colnames(obj)]
                  pheno=pData(GRset)
                  result=list(pheno,GRset)
                  return(result)

                }""")(obj, GRset) 
            
            
        if verbose:    
            robjects.r("""function (obj, orig_obj) { 
                start = length(rownames(orig_obj))
                left = length(rownames(obj))
                
                startcol = length(colnames(orig_obj))
                leftcol = length(colnames(obj))

                cat("\n In total there were ",start," probes for the analysis before filtering.")

                cat("\n",start-left," probes have been removed from further analysis.")

                cat("\n In total there are",left," probes left for the analysis.")
                
                
                cat("\n In total there were ",startcol," samples for the analysis before filtering.")

                cat("\n",startcol-leftcol," samples have been removed from further analysis.")

                cat("\n In total there are",leftcol," samples left for the analysis.")   

                }""")(obj, orig_obj)            
        
        return obj, pheno
       
    def beadCount(self, RGset=None):
        
        if RGset is None:
            RGset=self.RGset
        else:
            RGset=RGset

        bc=robjects.r("""function(RGset)                                               
                        {                                             
         #beadcount function, creates matrix with NAs representing probes with beadcount <3 from Extended RG Channel Set


        #' Creates matrix of beacounts from minfi data.
        #' 
        #' Creates matrix of beacounts from data read in using the minfi package.  NAs
        #' represent probes with beadcount <3.  An Extended RG Channel Set is required
        #' for this function to work.
        #' 
        #' 
        #' @param x 450K methylation data read in using minfi to create an Extended RG
        #' Channel Set
        #' @return A matrix of bead counts with bead counts <3 represented by NA for
        #' use in the pfilter function for quality control
        #' @note The beadcount function is internal to the pfilter function
        #' @author Ruth.Pidsley@@kcl.ac.uk
        #' @references [1] Pidsley R, Wong CCY, Volta M, Lunnon K, Mill J, Schalkwyk
        #' LC: A data-driven approach to preprocessing Illumina 450K methylation array
        #' data (submitted)
        #' @export beadcount
        beadcount<-function(x){
            #select out bead count dataframe
            getNBeads(x) -> nb
            #match rownames of beadcount dataframe to addresses
            getProbeInfo(x,type="I")->typeIadd
            match(typeIadd$AddressA,rownames(nb))->typeImatchA
            match(typeIadd$AddressB,rownames(nb))->typeImatchB

            #match rownames of beadcount dataframe to addresses
            getProbeInfo(x,type="II")->typeIIadd
            match(typeIIadd$Address,rownames(nb))->typeIImatch

            nb->nbcg

                locusNames <- getManifestInfo(x, "locusNames")
                bc_temp <- matrix(NA_real_, ncol = ncol(x), nrow = length(locusNames),
                dimnames = list(locusNames, sampleNames(x)))

                TypeII.Name <- getProbeInfo(x, type = "II")$Name
                bc_temp[TypeII.Name, ] <- nbcg[getProbeInfo(x, type = "II")$AddressA,]

            TypeI <- getProbeInfo(x, type = "I")

                bc_temp->bcB
                bc_temp->bcA

                bcB[TypeI$Name, ] <- nbcg[TypeI$AddressB,]
                bcA[TypeI$Name, ] <- nbcg[TypeI$AddressA,]

            which(bcB<3)->bcB3
            which(bcA<3)->bcA3
                bcA->bcA2
                bcB->bcB2

                bcA2[bcA3]<-NA
                bcB2[bcB3]<-NA
                bcA2[bcB3]<-NA
                bcB2[bcA3]<-NA

                data.frame(bcA2)->bcM
                data.frame(bcA2)->bcU
                data.frame(bcA2)->bc
                #row.names = rownames(bc)
                bc<-as.matrix(bc)
        }   


            bc<-beadcount(RGset)        

            return(bc)

            }""")(RGset)

        bc_py=self.ri2py_dataframe(bc, matrix=False)
        return bc, bc_py



    def getNBeads(self,RGset=None):
        
        if RGset is None:
            RGset=self.RGset
        else:
            RGset=RGset

        nBeads = robjects.r("""function(RGset)                                               
                {                                             
                 nbeads <- getNBeads(RGset)                         


                 return (nbeads)         

                }""")(RGset)
        return nBeads   
        




    def champ_processing(self, pheno=None, GRset=None, RGset=None, beta=None, M=None, 
                         autoimpute=True, filterDetP=True, ProbeCutoff=0, SampleCutoff=0.1,
                         filterBeads=True,beadCutoff=0.05,fixOutlier = True, dropSnPs=True, 
                         filterXY=True, filterNoCG=True, excludeXreactiveprobes=True, mask_probes=True, array_type='EPIC',                                  verbose=True, badSampleCutoff=10, rm_badsamples=True, detPFilter=False, detPcut=0.01, addQC=False, imputation_method="imputePCA"):
        
        if imputation_method!="methyLImp" and imputation_method!="imputePCA" and imputation_method!="knn":
            print('You did not specify a valid imputation method!!\n Choose "imputePCA" or "methyLImp" or "knn"...\n now using fallbackmethod "imputePCA"')
            imputation_method="imputePCA"   
       
        if pheno is None:
            pheno = self.pheno        
        else:
            pheno = pheno
        
        if RGset is None:
            RGset = self.RGset        
        else:
            RGset = RGset 
            
        if GRset is None:
            try:
                GRset = self.GRset
            except:
                GRset=None    
        else:
            GRset = GRset     
        if rm_badsamples: 
            if verbose:
                print('\n Now performing badsample removal')           
            RGset, pheno = self.remove_badsamples(badSampleCutoff=badSampleCutoff, rm_badsamples=rm_badsamples, 
                                              detPFilter=detPFilter, detPcut=detPcut, SampleCutoff=SampleCutoff,                                                                            addQC=addQC, verbose=verbose, RGset=RGset) 
            
            if verbose:
                print('\n Aligning RGset and GRset')
            GRset,pheno = robjects.r("""function (RGset, GRset) {  
                          #colnames(GRset) %in% colnames(RGset)
                          GRset <- GRset[,colnames(GRset) %in% colnames(RGset)]  
                          pheno<-pData(GRset)
                          result=list(GRset,pheno)
                          return(result)
                 }""")(RGset,GRset)
        
            
        
        if beta is None:
            _,beta = self.getBeta(Object = GRset if GRset else RGset)    
        else:
            beta = beta  

        if M is None:
            _,M = self.getM(Object = GRset if GRset else RGset)    
        else:
            M = M    

        detP, _ = self.detectionP(RGset=RGset)

        #self.getNBeads(RGset=RGset)     
        if filterBeads:
            nBeads,_ = self.beadCount(RGset=RGset)
        else:
            nBeads=None
            
        
        
        if verbose:
            print('\n Now performing champ_filter function') 

        self.beta_py, self.mval_py,self.pheno_py = self.champ_filter(
                                                                  beta=beta,
                                                                  M=M, 
                                                                  pheno=pheno, 
                                                                  detP=detP,
                                                                  beadcount=nBeads,
                                                                  autoimpute=autoimpute,
                                                                  filterDetP=filterDetP,
                                                                  ProbeCutoff=ProbeCutoff,
                                                                  SampleCutoff = SampleCutoff,
                                                                  detPcut= detPcut,
                                                                  filterBeads=filterBeads,
                                                                  beadCutoff= beadCutoff,
                                                                  fixOutlier = fixOutlier,
                                                                  imputation_method=imputation_method
                                                                 )
        
        GRset,self.pheno = robjects.r("""function (mval, GRset) {  
                          #colnames(GRset) %in% colnames(RGset)
                          GRset <- GRset[,colnames(GRset) %in% colnames(mval)]  
                          pheno<-pData(GRset)
                          result=list(GRset,pheno)
                          return(result)
                 }""")(self.mval,GRset)
        
        
        if filterNoCG or excludeXreactiveprobes or dropSnPs or filterXY or mask_probes is not False: 
                if verbose:
                    print('\n Now removing specific probes for m-values')
                self.mval, self.pheno=self.filterCpGs(obj=self.mval , 
                                    GRset=GRset if GRset else None,         
                                    dropSnPs=dropSnPs,                                     
                                    filterXY=filterXY, 
                                    filterNoCG=filterNoCG, 
                                    excludeXreactiveprobes=excludeXreactiveprobes,
                                    mask_probes=mask_probes,                   
                                    array_type=array_type, 
                                    verbose=verbose)
                
                try:
                    self.pheno_py=self.ri2py_dataframe(self.pheno, matrix=True)
                except:
                    self.pheno_py=self.ri2py_dataframe(self.pheno, matrix=False)        
                #self.mval=obj[0]
                self.mval_py = self.ri2py_dataframe(r_dat=self.mval, matrix=False)
                
                if verbose:
                    print('\n Now removing specific probes for beta-values ')
                self.beta, _ =self.filterCpGs(obj=self.beta , 
                                    GRset=GRset if GRset else None,         
                                    dropSnPs=dropSnPs,                                     
                                    filterXY=filterXY, 
                                    filterNoCG=filterNoCG, 
                                    excludeXreactiveprobes=excludeXreactiveprobes,
                                    mask_probes=mask_probes,           
                                    array_type=array_type, 
                                    verbose=verbose)
              
                self.beta_py = self.ri2py_dataframe(r_dat=self.beta, matrix=False)           
        
        
        
        
            
        return self.beta_py, self.mval_py, self.pheno_py




    def champ_filter(self,
                      beta=None,
                      M=None, 
                      pheno=None, 
                      detP=None,
                      beadcount=None,
                      autoimpute=True,
                      filterDetP=True,
                      ProbeCutoff=0,
                      SampleCutoff=0.1,
                      detPcut=0.01,
                      filterBeads=True,
                      beadCutoff=0.05,
                      fixOutlier = True,
                      imputation_method="imputePCA"):
            
            if pheno==None:
                pheno=self.pheno
            
            
            if imputation_method!="methyLImp" and imputation_method!="imputePCA" and imputation_method!="knn":
                print('You did not specify a valid imputation method!!\n Choose "imputePCA" or "methyLImp" or "knn"...\n now using fallbackmethod "imputePCA"')
                imputation_method="imputePCA"             


            if (imputation_method=="methyLImp"): 
                methyLImp=self.methyLImp()            


            if (imputation_method!="methyLImp"):
                methyLImp=False

            if imputation_method=="imputePCA":
                imputePCA=self.imputePCA()


            if (imputation_method!="imputePCA"):
                imputePCA=False
            
            
            self.beta, self.mval, self.pheno = robjects.r("""function( beta, 
                                                                      M,
                                                                      pd,
                                                                      detP,
                                                                      beadcount,
                                                                      autoimpute,
                                                                      filterDetP,
                                                                      ProbeCutoff,
                                                                      SampleCutoff,
                                                                      detPcut,
                                                                      filterBeads,
                                                                      beadCutoff,
                                                                      fixOutlier,
                                                                      imputation_method,
                                                                      imputePCA,
                                                                      methyLImp
                                                                      )    
         { 
            
            champ.filter <- function(beta,
                                     M,
                                     pd,
                                     detP,
                                     beadcount,
                                     autoimpute,
                                     filterDetP,
                                     ProbeCutoff,
                                     SampleCutoff,
                                     detPcut,
                                     filterBeads,
                                     beadCutoff,
                                     fixOutlier,
                                     imputation_method,
                                     imputePCA,
                                     methyLImp
                                     )
           {
            cat("[===========================]")
            cat("<<<FILTER START>>> ")
            cat("-----------------------------")

            
           
            cat(" Section 1:  Check Input Start ")
            R_MAX_MEM_SIZE=memory.limit(size = NA)
            Objects <- list("beta"=beta,
                            "M"=M)

            Objects <- Objects[which(lapply(Objects,FUN=is.null)==FALSE)]

            if(length(Objects)==0) stop("  At least one Data Set needed.")
            cat("  You have inputed ",paste(names(Objects),collapse=",")," for Analysis.")
            if(length(unique(lapply(Objects,FUN=rownames))) != 1) stop("  !!!  You objects have different rownames. Please make sure your Matrix are in accord on Rows.")
            if(length(unique(lapply(Objects,FUN=colnames))) != 1) stop("  !!!  You objects have different colnames. Please make sure your Matrix are in accord on Rows.")
            
            Accessory <- list("detP"=detP,
                              "beadcount"=beadcount)

            FilterOption <- list("filterDetP"=filterDetP,
                                 "autoimpute"=autoimpute,
                                 "filterBeads"=filterBeads
                                 )

            ### Checking pd file
            if(!is.null(pd)) 
           {
                cat(" pd file provided, checking if it's in accord with Data Matrix...")
                if(nrow(pd) == ncol(Objects[[1]])) {
                    if("ID" %in% names(pd)){
                        if(identical(as.character(pd$ID),colnames(Objects[[1]]))) 
                            cat("    pd file check success.")
                        else
                            stop("    Your pd file's ID is different from your Data Matrix colnames.")
                    } else {
                        cat("    !!! Your pd file does not have ID column, we can not check your Sample_Name, please make sure the pd file is correct.")
                    }
                } else {
                    stop("    pd file and Data matrix have different dimensions.")
                }
            }



           
           
           

            ### Checking Detect P value
            if(FilterOption$filterDetP == TRUE)
            {   



                cat("  Parameter filterDetP is TRUE, checking if detP in accord with Data Matrix...")
                if(!is.null(Accessory$detP)) {
                    if(identical(rownames(Accessory$detP),rownames(Objects[[1]])) & identical(colnames(Accessory$detP),colnames(Objects[[1]])))
                    {
                        cat("    detP check success.")
                    }
                    else {
                        
                        cat("    !!! Your detP matrix has been aligned to match the EXACT same rowname and colname as Data Matrix.")
                        
                        #keep <-match(colnames(nbeads), colnames(m))
                        
                        if (nrow(Accessory$detP) > nrow(Objects[[1]]))
                        {
                        
                            keep <- ( rownames(Accessory$detP) %in% rownames(Objects[[1]]))
                            Accessory$detP <- Accessory$detP[keep,]                       
                        }
                        
                        else if (nrow(Accessory$detP) < nrow(Objects[[1]]))
                        {
                        
                            keep <- ( rownames(Objects[[1]]) %in% rownames(Accessory$detP))
                            Accessory$detP <- Accessory$detP[keep,] 
                        
                        }
                        
                        else if (nrow(Accessory$detP) == nrow(Objects[[1]]))
                        {
                            keep <- match(rownames(Objects[[1]]), rownames(Accessory$detP))

                            Accessory$detP <- Accessory$detP[keep,]
                        
                        }
                        
                        
                        if (ncol(Accessory$detP) <= ncol(Objects[[1]]))
                        {
                        
                            keep <- match( colnames(Accessory$detP), colnames(Objects[[1]]))
                            Accessory$detP <- Accessory$detP[,keep]                       
                        }
                        
                        else if (ncol(Accessory$detP) > ncol(Objects[[1]]))
                        {
                        
                            keep <- match(colnames(Objects[[1]]), colnames(Accessory$detP))
                            Accessory$detP <- Accessory$detP[,keep] 
                        
                        }
                       
                        
                                         
                        
                        
                        
                    }  
                } else {
                    cat("    !!! Parameter detP is not found, filterDetP is reset FALSE now.")
                    FilterOption$filterDetP <- FALSE
                    Accessory$detP <- NULL
                }
            }






            ### Checking Beadcount value
            if(FilterOption$filterBeads == TRUE)
            {
                cat("  Parameter filterBeads is TRUE, checking if beadcount in accord with Data Matrix...")
                if(!is.null(Accessory$beadcount)) {
                    if(identical(rownames(Accessory$beadcount),rownames(Objects[[1]])) & identical(colnames(Accessory$beadcount),colnames(Objects[[1]])))
                    {
                        cat("    beadcount check success.")
                    }
                    else {
                                                                       
                        cat("    !!! Your  beadcount matrix has been aligned to match the EXACT same rowname and colname as Data Matrix.")
                        
                                 
                        
                        if (nrow(Accessory$beadcount) > nrow(Objects[[1]]))
                        {
                        
                            keep <- ( rownames(Accessory$beadcount) %in% rownames(Objects[[1]]))
                            Accessory$beadcount <- Accessory$beadcount[keep,]                       
                        }
                        
                        else if (nrow(Accessory$beadcount) < nrow(Objects[[1]]))
                        {
                        
                            keep <- ( rownames(Objects[[1]]) %in% rownames(Accessory$beadcount))
                            Accessory$beadcount <- Accessory$beadcount[keep,] 
                        
                        }
                        
                        else if (nrow(Accessory$beadcount) == nrow(Objects[[1]]))
                        {
                            keep <- match(rownames(Objects[[1]]), rownames(Accessory$beadcount))

                            Accessory$beadcount <- Accessory$beadcount[keep,]
                        
                        }
                        
                        
                        if (ncol(Accessory$beadcount) <= ncol(Objects[[1]]))
                        {
                        
                            keep <- match( colnames(Accessory$beadcount), colnames(Objects[[1]]))
                            Accessory$beadcount <- Accessory$beadcount[,keep]                       
                        }
                        
                        else if (ncol(Accessory$beadcount) > ncol(Objects[[1]]))
                        {
                        
                            keep <- match(colnames(Objects[[1]]), colnames(Accessory$beadcount))
                            Accessory$beadcount <- Accessory$beadcount[,keep] 
                        
                        }
                       
                        
                        
                        
                        
                        
                        
                    }  
                } else {
                    cat("    !!! Parameter beadcount is not found, filterBeads is reset FALSE now.")
                    FilterOption$filterBeads <- FALSE
                    Accessory$beadcount <- NULL
                }
            }


            if(FilterOption$autoimpute == TRUE)
            {
                cat("  parameter autoimpute is TRUE. Checking if the conditions are fulfilled...")
                if("beta" %in% names(Objects) | "M" %in% names(Objects)){
                    if(!is.null(detP)) {
                        cat("    autoimpute check success.")
                    } else {
                        cat("    You need to provide a detP value matrix for this option to work. autoimpute has been reset FALSE.")
                        FilterOption$autoimpute <- FALSE
                    }
                } else {
                    cat("    !!! beta matrix or M matrix are required for impute. autoimpute has been reset FALSE.")
                    FilterOption$autoimpute <- FALSE
                }
            }
            
            #cat(nrow(Objects$M)) 
            #cat(nrow(Objects$beta) )
            #cat(nrow( Accessory$beadcount) )
            #cat(nrow( Accessory$detP) )
            
            ### Start Filtering Here
            cat(" Section 2: Filtering Start >>")  


            if(FilterOption$filterDetP == TRUE)
            {
                cat("  Filtering Detect P value Start")
                cat("    The fraction of failed positions per sample")
                cat("    You may need to delete samples with high proportion of failed probes:\n")

                numfail <- matrix(colMeans(Accessory$detP >= detPcut))
                rownames(numfail) <- colnames(Accessory$detP)
                colnames(numfail) <- "Failed CpG Fraction."
                
                RemainSample <- which(numfail < SampleCutoff)
                
                
                if(any(numfail >= SampleCutoff))
                {  
                    cat("    The detSamplecut parameter is : ",SampleCutoff)
                    cat("    Samples : ", paste(rownames(numfail)[which(numfail >= SampleCutoff)],collapse=",")," will be deleted.")
                    cat("    There are ",length(RemainSample)," samples remained for analysis.")
                    Objects <- lapply(Objects,function(x) x[,RemainSample])
                    Accessory <- lapply(Accessory,function(x) x[,RemainSample])
                   
                }
               
               
               
               RemainProbe <- rowSums(Accessory$detP > detPcut) <= ProbeCutoff * length(RemainSample)
               if(ProbeCutoff==0)
               {
                   cat("\n    Filtering probes with a detection p-value above ",detPcut,".")
                   cat("    Removing ",sum(RemainProbe==FALSE)," probes.")
                  cat("    If a large number of probes have been removed, we suggests you to identify potentially bad samples")
               } else {
                   cat("\n    Filtering probes with a detection p-value above ",detPcut," in at least", ProbeCutoff*100,"% Samples.")
                   cat("    Removing ",sum(RemainProbe==FALSE)," probes.")
                   cat("    If a large number of probes have been removed, we suggests you to identify potentially bad samples")
               }

               Objects <- lapply(Objects,function(x) x[RemainProbe,])
               Accessory <- lapply(Accessory,function(x) x[RemainProbe,])
               if(sum(Accessory$detP > detPcut) > 0) cat(" There are still ",sum(Accessory$detP > detPcut), " failed probes exist in your data set, imputation is recommended.")

                
            }


           
           
           
           
           
            if(FilterOption$autoimpute == TRUE)
            {
               cat("  Autoimpute Start") 
               if(sum(Accessory$detP > detPcut) == 0) 
               {
                   cat("    No NAs (failed probes) exist in your data set any more, you don't need to do any imputation.")
               } else
               {
                   cat("    There are ",sum(Accessory$detP > detPcut), " NAs (failed probes) exists in your data set.")
                   cat("    Imputing will be conducted for remain NAs. \n") 

                                     
                   if("beta" %in% names(Objects))
                   {
                      cat("    Doing imputation on beta matrix.")
                      Objects$beta[Accessory$detP > detPcut] <- NA
                      if (imputation_method=="knn")
                                     {  
                                     cat("   Performing impute.knn.")
                                     ### Using sink function to remove messages
                                     zz <- file("ImputeMessage.Rout", open="wt")
                                     sink(zz)
                                     sink(zz, type="message")                                     
                                      library(impute)
                                     
                                     Objects$beta <- impute.knn(Objects$beta,k=5)$data
                                     
                                     sink(type="message")
                                     sink()  
                                     
                                     }
                      
                      if(imputation_method=="methyLImp")
                                     { cat("   Performing methyLImp")
                                         methyLImps<-methyLImp$methyLImp
                                         plogit<-methyLImp$plogit
                                         inv.plogit<-methyLImp$inv.plogit 
                                         pinvr<-methyLImp$pinvr                                         
                                         
                                         Objects$beta<-methyLImps(t(Objects$beta), min = 0, max = 1, max.sv = NULL, col.list = NULL)
                                         Objects$beta<-t(Objects$beta)
                                        
                                     }
                      
                      if(imputation_method=="imputePCA")
                                     { cat("   Performing imputePCA")
                                       imputePCAs<-imputePCA$imputePCA
                                       FactoMineR.svd.triplet<-imputePCA$FactoMineR.svd.triplet
                                                                   
                                         
                                         Objects$beta <- imputePCAs(Objects$beta)
                                         if(class(Objects$beta)!="matrix")
                                         {
                                             Objects$beta <- Objects$beta$completeObs
                                         } 
                                         
                                         Objects$beta[Objects$beta<0] <- 0
                                         Objects$beta[Objects$beta>1] <- 1  
                                         
                                     }
                      
                      
                      
                      
                       
                   }
                   if("M" %in% names(Objects))
                   {
                       cat("    Doing imputation on M matrix.")
                       Objects$M[Accessory$detP > detPcut] <- NA
                       
                       if (imputation_method=="knn")
                                     {  
                                     cat("   Performing impute.knn.")
                                     ### Using sink function to remove messages
                                     zz <- file("ImputeMessage.Rout", open="wt")
                                     sink(zz)
                                     sink(zz, type="message")                                     
                                     library(impute)
                                     
                                     Objects$M <- impute.knn(Objects$M,k=5)$data
                                     
                                     sink(type="message")
                                     sink()  
                                     
                                     }
                      
                      if(imputation_method=="methyLImp")
                                     { cat("   Performing methyLImp")
                                         methyLImps<-methyLImp$methyLImp
                                         plogit<-methyLImp$plogit
                                         inv.plogit<-methyLImp$inv.plogit 
                                         pinvr<-methyLImp$pinvr                                         
                                         
                                         Objects$M<-methyLImps(t(Objects$M), min = -16, max = 16, max.sv = NULL, col.list = NULL)
                                         Objects$M<-t(Objects$M)
                                        
                                     }
                      
                      if(imputation_method=="imputePCA")
                                     { cat("   Performing imputePCA")
                                       imputePCAs<-imputePCA$imputePCA
                                       FactoMineR.svd.triplet<-imputePCA$FactoMineR.svd.triplet
                                                                   
                                         
                                         Objects$M <- imputePCAs(Objects$M)
                                         
                                         if(class(Objects$M)!="matrix")
                                         {
                                             Objects$M <- Objects$M$completeObs
                                         } 
                                                                                  
                                     }                
                       
                       
                       
                   }
                  
               }
            }
             
           #keep <- (rownames(Accessory$detP) %in% rownames(Objects[[1]]))
           #Accessory$detP <- Accessory$detP[keep,]  
           #keep <- (rownames(Accessory$beadcount) %in% rownames(Objects[[1]]))
           #Accessory$beadcount <- Accessory$beadcount[keep,]
           
           #cat(nrow(Objects$M)) 
           #cat(nrow(Objects$beta) )
           #cat(nrow( Accessory$beadcount) )
           #cat(nrow( Accessory$detP) )
           
           if(FilterOption$filterBeads == TRUE)
            {   
                
                cat("  Filtering BeadCount Start")
                RemainProbe <- rowSums(is.na(Accessory$beadcount)) < beadCutoff*(ncol(Accessory$beadcount))
                
                cat("    Filtering probes with a beadcount <3 in at least ",beadCutoff*100,"% of samples.")
                cat("    Removing ",sum(RemainProbe == FALSE)," probes")
                
                Objects <- lapply(Objects,function(x) x[RemainProbe,])
                Accessory <- lapply(Accessory,function(x) x[RemainProbe,])
            } 

            
            #rownames(Accessory$detP) %in% rownames(Objects[[1]]))
            #mval <- mval[,colnames(mval) %in% rownames(pheno)] 
            if(!is.null(pd))
            {
                cat("  Updating PD file")
                if(FilterOption$filterDetP==FALSE)
                {
                    cat("    filterDetP parameter is FALSE, so no Sample Would be removed.")
                    RemainSample <- 1:nrow(pd)
                }
                pd <- pd[RemainSample,]
                Objects <- append(Objects,list(pd=pd))
            }




            

            if(fixOutlier & "beta" %in% names(Objects))
            {
                cat("  Fixing Outliers Start")
                cat("    Replacing all value smaller/equal to 0 with smallest positive value.")
                Objects$beta[Objects$beta <= 0] <- min(Objects$beta[which(Objects$beta > 0)])
                cat("    Replacing all value greater/equal to 1 with largest value below 1..")
                Objects$beta[Objects$beta >= 1] <- max(Objects$beta[which(Objects$beta < 1)])
            }



            cat("[ Section 2: Filtering Done ]")

            cat(" All filterings are Done, now you have ", nrow(Objects[[1]]), " probes and ",ncol(Objects[[1]]), " samples.\n")

            cat("[<<<<< ChAMP.FILTER END >>>>>>]")
            cat("[===========================]")
            cat("[You may want to process champ.QC() next.]\n")
            return(Objects)
      }
    


            Objects<-champ.filter(beta,
                                     M,
                                     pd,
                                     detP,
                                     beadcount,
                                     autoimpute,
                                     filterDetP,
                                     ProbeCutoff,
                                     SampleCutoff,
                                     detPcut,
                                     filterBeads,
                                     beadCutoff,
                                     fixOutlier,
                                     imputation_method,
                                     imputePCA,
                                     methyLImp
                                     )

            result=list(Objects$beta, Objects$M, Objects$pd)

            return (result)         

                                            
                                            
                                            }""")(beta, 
                                                  M,
                                                  pheno,                                                                                  
                                                  detP,
                                                  beadcount,
                                                  autoimpute,
                                                  filterDetP,
                                                  ProbeCutoff,
                                                  SampleCutoff,
                                                  detPcut,
                                                  filterBeads,
                                                  beadCutoff,                                            
                                                  fixOutlier,
                                                  imputation_method,
                                                  imputePCA,
                                                  methyLImp
                                                  )       
            
            
            
            self.mval_py = self.ri2py_dataframe(r_dat=self.mval, matrix=False)

            self.beta_py = self.ri2py_dataframe(r_dat=self.beta, matrix=False)
            try:
                self.pheno_py = self.ri2py_dataframe(r_dat=self.pheno, matrix=False)
            except:
                self.pheno_py = self.ri2py_dataframe(r_dat=self.pheno, matrix=True)

            
                
            return self.beta_py, self.mval_py, self.pheno_py

    def methyLImp(self):
        methyLImp=robjects.r("""function(){
        
        
            ## methyLImp - linear imputation model for continuous variables
            ## 
            ## This R script contains the methyLImp function and other
            ## subroutines needed for calculation.
            ## 
            ## Author: P. Di Lena, pietro.dilena@unibo.it


            # Pseudo-logit function with range restricted to [-X,X] where
            # X depends on the double machine precision.
            plogit <- function(x, min=0, max=1)
            {
                p <- (x-min)/(max-min)
                # fix -Inf
                p <- ifelse(p <   .Machine$double.neg.eps,  .Machine$double.neg.eps,p)
                # fix +Inf
                p <- ifelse(p > 1-.Machine$double.neg.eps,1-.Machine$double.neg.eps,p)
                log(p/(1-p))
            }

            # Inverse of the pseudo-logit function.
            inv.plogit <- function(x, min=0, max=1)
            {
                p <- exp(x)/(1+exp(x))
                # fix problems with +Inf
                p <- ifelse(is.na(p) & !is.na(x), 1, p )               
                # fix 0 rounding
                p <- ifelse(p <= exp(plogit(0))/(1+exp(plogit(0))), 0, p)
                p * (max-min) + min
            }

            # Computes the Moore-Penrose generalized inverse of a matrix. Allows rank
            # reduction of the generalized inverse.
            # 
            # This function is directly taken from MASS package (code on GPLv3 license)
            # and modified in order to include the rank reduction option. The added code
            # for rank reduction is commented in the implementation.
            #
            pinvr <- function(X, max.sv = min(dim(X)), tol = sqrt(.Machine$double.eps))
            {
                #
                # based on suggestions of R. M. Heiberger, T. M. Hesterberg and WNV
                #
                if(length(dim(X)) > 2L || !(is.numeric(X) || is.complex(X)))
                    stop("'X' must be a numeric or complex matrix")
                if(!is.matrix(X)) X <- as.matrix(X)
                Xsvd <- svd(X)
                if(is.complex(X)) Xsvd$u <- Conj(Xsvd$u)
                Positive <- Xsvd$d > max(tol * Xsvd$d[1L], 0)

                # Rank reduction extension: START
                max.sv      <- min(ifelse(max.sv < 0, 1, max.sv),min(dim(X)))
                L           <- logical(length(Positive))
                L[seq_len(max.sv)] <- TRUE
                Positive    <- Positive & L
                # Rank reduction extension: END

                if (all(Positive)) Xsvd$v %*% (1/Xsvd$d * t(Xsvd$u))
                else if(!any(Positive)) array(0, dim(X)[2L:1L])
                else Xsvd$v[, Positive, drop=FALSE] %*% ((1/Xsvd$d[Positive]) * t(Xsvd$u[, Positive, drop=FALSE]))
            }

            # methylation data Linear Imputation model
            #
            # Arguments:
            # dat      = data matrix with missing values, where samples are on the rows 
            #            and variables on the columns
            # min      = minimum value for bounded-range variables. Default: 0 (we assume
            #            beta-value representation of methylation data). Unrestricted range
            #            if min or max = NULL        
            # max      = maximum value for bounded-range variables. Default: 1 (we assume
            #            beta-value representation of methylation data). Unrestricted range
            #            if min or max = NULL.
            # max.sv   = maximum number of singuar value to be used in the pseudoinverse.
            #            If NULL use all singular values.
            # col.list = restricts the imputation on the specified columns.
            #      
        
        
        
            methyLImp <- function(dat, min = 0, max = 1, max.sv = NULL, col.list = NULL)
        {
            out    <- dat
            NAcols <- colSums(is.na(dat)) > 0; NAcols <- which(NAcols)

            # Convert col.list, if any, from names to numbers
            if(is.character(col.list)) {
                if(is.null(colnames(dat)))
                    col.list <- NULL
                else
                    col.list <- which(colnames(dat) %in% col.list)
            }

            # If all the colums have a missing value we cannot do anything
            if(length(NAcols) < ncol(dat)) {
                # Columns with all NAs or a single not NA value (to be excluded:
            # not enough information for imputation)
                NAall  <- colSums(is.na(dat))<(nrow(dat)-1); NAall <- which(NAall)
                # List of columns to impute
                NAlist <- intersect(NAcols, NAall)
                # Filter the columns to impute according to col.list
                if(!is.null(col.list))
                    NAlist <- intersect(NAlist,col.list)

                while(length(NAlist) != 0) {
                    col_id <- NAlist[1]

                    # List of rows for which col_id is NA
                    row_id <- which(is.na(dat[,col_id])==TRUE)

                    # Colum indexes of NA columns for all the row_id(s)
                    if(length(row_id) == 1)
                        tmp1 <- which(is.na(dat[row_id,])==TRUE)
                    else
                        tmp1 <- which(colSums(is.na(dat[row_id,])) == length(row_id))

                    # Column indexes: no colum element is NA for the rows not in row_id
                    tmp2   <- which(colSums(is.na(dat[-row_id,])) == 0)

                    # List of colums in NAlist that are NA only for all the row_id(s)
                    NAcols_rowid <- intersect(intersect(tmp1,tmp2),NAlist)

                    # Extract submatrices for regression
                    A <- dat[-row_id,-NAcols]
                    B <- dat[-row_id,NAcols_rowid]
                    C <- dat[row_id,-NAcols]

                    # Updates or computes max.sv from A. Negative or zero value not allowed
                    max.sv <- max(ifelse(is.null(max.sv),min(dim(A)),max.sv),1)

                    if(is.null(min) || is.null(max)) {
                        # Unrestricted-range imputation
                        # X <- pinvr(A,rank)%*%B (X = A^-1*B)
                        # O <- C%*%X             (O = C*X)
                        out[row_id,NAcols_rowid] <- C%*%(pinvr(A,max.sv)%*%B)
                    } else { 
                        # Bounde-range imputation
                        # X <- pinvr(A,rank)%*%logit(B,min,max) (X = A^-1*logit(B))
                        # P <- inv.logit(C%*%X,min,max)         (P = logit^-1(C*X))
                        out[row_id,NAcols_rowid] <- inv.plogit(C%*%(pinvr(A,max.sv)%*%plogit(B,min,max)),min,max)
                    }

                    # Update NA column list
                    NAlist <- setdiff(NAlist,NAcols_rowid)
                }
            }
            return(out)
        }
        
        result=list(plogit=plogit, inv.plogit=inv.plogit, pinvr=pinvr, methyLImp= methyLImp)
        
        return(result)
        }""")()
        
        return methyLImp

    def imputePCA(self):
        
        imputePCA=robjects.r("""function(){
             
            FactoMineR.svd.triplet = function (X, row.w = NULL, col.w = NULL,ncp=Inf) {
 
 
            tryCatch.W.E <- function(expr){  ## function proposed by Maechler
                W <- NULL
                w.handler <- function(w){ # warning handler
                    W <<- w
                    invokeRestart("muffleWarning")
                }
                list(value = withCallingHandlers(tryCatch(expr, error = function(e) e),
                                                 warning = w.handler),
                     warning = W)
            }
               if (is.null(row.w)) row.w <- rep(1/nrow(X), nrow(X))
               if (is.null(col.w)) col.w <- rep(1, ncol(X))
               ncp <- min(ncp,nrow(X)-1,ncol(X))
               row.w <- row.w / sum(row.w)
                X <- t(t(X)*sqrt(col.w))*sqrt(row.w)
            if (ncol(X)==1){
                  mult <- sign(as.vector(crossprod(rep(1,nrow(V)),as.matrix(V))))
                  mult[mult==0] <- 1
                  U <- t(t(U)*mult)
                  V <- t(t(V)*mult)
                
                U <- U/sqrt(row.w)
                V <- V/sqrt(col.w)
            }
            else{
                svd.usuelle <- tryCatch.W.E(svd(t(X),nu=ncp,nv=ncp))$val
                if (names(svd.usuelle)[[1]]=="message"){
                  svd.usuelle <- tryCatch.W.E(svd(X,nu=ncp,nv=ncp))$val
                  if (names(svd.usuelle)[[1]]=="d"){
                    aux <- svd.usuelle$u
                    svd.usuelle$u <- svd.usuelle$v
                    svd.usuelle$v <- aux
                  } else{
                      bb <- eigen(crossprod(t(X),t(X)),symmetric=TRUE)
                      svd.usuelle <- vector(mode = "list", length = 3)
                      svd.usuelle$d[svd.usuelle$d<0]=0
                      svd.usuelle$d <- sqrt(svd.usuelle$d)
                      svd.usuelle$v <- bb$vec[,1:ncp]
                      svd.usuelle$u <- t(t(crossprod(X,svd.usuelle$v))/svd.usuelle$d[1:ncp])
                  }
                }
                U <-  svd.usuelle$v
                V <- svd.usuelle$u
                mult <- sign(as.vector(crossprod(rep(1,nrow(V)),as.matrix(V))))
                mult[mult==0] <- 1
                V <- t(t(V)*mult)/sqrt(col.w)
                U <- t(t(U)*mult)/sqrt(row.w)
            }
                vs <- svd.usuelle$d[1:min(ncol(X),nrow(X)-1)]
                num <- which(vs[1:ncp]<1e-15)
                if (length(num)==1){
                  U[,num] <- U[,num,drop=FALSE]*vs[num]
                  V[,num] <- V[,num,drop=FALSE]*vs[num]
                } 
                if (length(num)>1){
                  U[,num] <- t(t(U[,num])*vs[num])
                  V[,num] <- t(t(V[,num])*vs[num])
                }
                res <- list(vs = vs, U = U, V = V)
                return(res)
            }
        
             
            
            imputePCA <- function (X, ncp = 2, scale=TRUE, method=c("Regularized","EM"),row.w=NULL,coeff.ridge=1,threshold = 1e-6,seed = NULL,nb.init=1,maxiter=1000,...){
        impute <- function (X, ncp = 4, scale=TRUE, method=NULL,threshold = 1e-6,seed = NULL,init=1,maxiter=1000,row.w=NULL,coeff.ridge=1,...){
            moy.p <- function(V, poids) {
                res <- sum(V * poids,na.rm=TRUE)/sum(poids[!is.na(V)])
            }
            ec <- function(V, poids) {
                res <- sqrt(sum(V^2 * poids,na.rm=TRUE)/sum(poids[!is.na(V)]))
            }
           nb.iter <- 1
           old <- Inf
           objective <- 0
           if (!is.null(seed)){set.seed(seed)}
           X <- as.matrix(X)
           ncp <- min(ncp,ncol(X),nrow(X)-1)
           missing <- which(is.na(X))
           mean.p <- apply(X, 2, moy.p,row.w)
           Xhat <- t(t(X)-mean.p)
           et <- apply(Xhat, 2, ec,row.w)
           if (scale) Xhat <- t(t(Xhat)/et)
           if (any(is.na(X))) Xhat[missing] <- 0
           if (init>1) Xhat[missing] <- rnorm(length(missing)) ## random initialization
           fittedX <- Xhat
           if (ncp==0) nb.iter=0
           while (nb.iter > 0) {
               Xhat[missing] <- fittedX[missing]
               if (scale) Xhat=t(t(Xhat)*et)
               Xhat <- t(t(Xhat)+mean.p)
               mean.p <- apply(Xhat, 2, moy.p,row.w)
               Xhat <- t(t(Xhat)-mean.p)
               et <- apply(Xhat, 2, ec,row.w)
               if (scale) Xhat <- t(t(Xhat)/et)
               svd.res <- FactoMineR.svd.triplet(Xhat,row.w=row.w,ncp=ncp)
        #       sigma2 <- mean(svd.res$vs[-(1:ncp)]^2)
               sigma2  <- nrow(X)*ncol(X)/min(ncol(X),nrow(X)-1)* sum((svd.res$vs[-c(1:ncp)]^2)/((nrow(X)-1) * ncol(X) - (nrow(X)-1) * ncp - ncol(X) * ncp + ncp^2))
               sigma2 <- min(sigma2*coeff.ridge,svd.res$vs[ncp+1]^2)
               if (method=="em") sigma2 <-0
               lambda.shrinked=(svd.res$vs[1:ncp]^2-sigma2)/svd.res$vs[1:ncp]
               fittedX = tcrossprod(t(t(svd.res$U[,1:ncp,drop=FALSE]*row.w)*lambda.shrinked),svd.res$V[,1:ncp,drop=FALSE])
               fittedX <- fittedX/row.w
               diff <- Xhat-fittedX
               diff[missing] <- 0
               objective <- sum(diff^2*row.w)
        #       objective <- mean((Xhat[-missing]-fittedX[-missing])^2)
               criterion <- abs(1 - objective/old)
               old <- objective
               nb.iter <- nb.iter + 1
               if (!is.nan(criterion)) {
                 if ((criterion < threshold) && (nb.iter > 5))  nb.iter <- 0
                 if ((objective < threshold) && (nb.iter > 5))  nb.iter <- 0
               }
               if (nb.iter > maxiter) {
                 nb.iter <- 0
                 warning(paste("Stopped after ",maxiter," iterations"))
               }
           }
           if (scale) Xhat <- t(t(Xhat)*et)
           Xhat <- t(t(Xhat)+mean.p)
           completeObs <- X
           completeObs[missing] <- Xhat[missing]
           if (scale) fittedX <- t(t(fittedX)*et)
           fittedX <- t(t(fittedX)+mean.p)
           result <- list()
           result$completeObs <- completeObs
           result$fittedX <- fittedX
           return(result) 
        }
        #### Main program
         method <- match.arg(method,c("Regularized","regularized","EM","em"),several.ok=T)[1]
         obj=Inf
         method <- tolower(method)
         if (ncp>min(nrow(X)-2,ncol(X)-1)) stop("ncp is too large")
         if (is.null(row.w)) row.w = rep(1,nrow(X))/nrow(X)
         for (i in 1:nb.init){
          if (!any(is.na(X))) return(X)
          res.impute=impute(X, ncp=ncp, scale=scale, method=method, threshold = threshold,seed=if(!is.null(seed)){(seed*(i-1))}else{NULL},init=i,maxiter=maxiter,row.w=row.w,coeff.ridge=coeff.ridge)
          if (mean((res.impute$fittedX[!is.na(X)]-X[!is.na(X)])^2) < obj){
            res <- res.impute
            obj <- mean((res.impute$fittedX[!is.na(X)]-X[!is.na(X)])^2)
          }
         }
        return(res)
        }
        
        result=list(FactoMineR.svd.triplet=FactoMineR.svd.triplet, imputePCA=imputePCA)
        return(result) 
        
        }""")()
        
        return imputePCA
        
    def detectionP_NA(self):
        
        detectionP=robjects.r("""function(){
        library(matrixStats)
        detectionP <- function(rgSet, type = "m+u", na.rm=FALSE) {
            
            locusNames <- getManifestInfo(rgSet, "locusNames")
            detP <- matrix(NA_real_, ncol = ncol(rgSet), nrow = length(locusNames),
                           dimnames = list(locusNames, sampleNames(rgSet)))

            controlIdx <- getControlAddress(rgSet, controlType = "NEGATIVE")
            r <- getRed(rgSet)
            rBg <- r[controlIdx,]
            rMu <- matrixStats::colMedians(rBg, na.rm=na.rm)
            rSd <- matrixStats::colMads(rBg, na.rm=na.rm)

            g <- getGreen(rgSet)
            gBg <- g[controlIdx,]
            gMu <- matrixStats::colMedians(gBg, na.rm=na.rm)
            gSd <- matrixStats::colMads(gBg, na.rm=na.rm)

            TypeII <- getProbeInfo(rgSet, type = "II")
            TypeI.Red <- getProbeInfo(rgSet, type = "I-Red")
            TypeI.Green <- getProbeInfo(rgSet, type = "I-Green")
            for (i in 1:ncol(rgSet)) {
                ## Type I Red
                intensity <- r[TypeI.Red$AddressA, i] + r[TypeI.Red$AddressB, i]
                detP[TypeI.Red$Name, i] <- 1-pnorm(intensity, mean=rMu[i]*2, sd=rSd[i]*2)
                ## Type I Green
                intensity <- g[TypeI.Green$AddressA, i] + g[TypeI.Green$AddressB, i]
                detP[TypeI.Green$Name, i] <- 1-pnorm(intensity, mean=gMu[i]*2, sd=gSd[i]*2)
                ## Type II
                intensity <- r[TypeII$AddressA, i] + g[TypeII$AddressA, i]
                detP[TypeII$Name, i] <- 1-pnorm(intensity, mean=rMu[i]+gMu[i], sd=rSd[i]+gSd[i])
            }
            detP
        }
        return(detectionP)
        
        }""")()
        
        return detectionP
        
    def reduce(self, GRset=None, RGset=None, what="both", detPcut=0.01, SampleCutoff=0.1, ProbeCutoff=0.1, verbose=True, autoimpute=True, imputation_method="imputePCA"):
       
        if GRset:
            GRset=GRset
        else:
            GRset = self.GRset
        if RGset:
            RGset = RGset
        else:
            RGset = self.RGset_filt 
        
        if imputation_method!="methyLImp" and imputation_method!="imputePCA" and imputation_method!="knn":
            print('You did not specify a valid imputation method!!\n Choose "imputePCA" or "methyLImp" or "knn"...\n now using fallbackmethod "imputePCA"')
            imputation_method="imputePCA"                
                 
            
        if (imputation_method=="methyLImp"): 
            methyLImp=self.methyLImp()            
                        
            
        if (imputation_method!="methyLImp"):
            methyLImp=False
            
        if imputation_method=="imputePCA":
            imputePCA=self.imputePCA()
            
        
        if (imputation_method!="imputePCA"):
            imputePCA=False
           
        detectionP=self.detectionP_NA()
        
        obj, self.pheno = robjects.r("""function (GRset, RGset, what=c("beta", "M","both"), cutp, cutsamples, cutcpgs, verbose, autoimpute, methyLImp, imputePCA, imputation_method,detectionPs,...) { 

            ##' Extract functional normalized data according to filter data
            ##'
            ##' Since functional normalization doesn't accept NA e.g. after probe filter the
            ##' normalized 'GenomicRatioSet' and filtered 'RGChannelSet' need to be merged.
            ##' Additionally if M-values are calculated we set -/+Inf values to +/-16
            ##' more details
            ##' @title Extract functional normalized data according to filter data
            ##' @param GRset 'GenomicRatioSet' output from functional normalization
            ##' @param RGset 'RGset' output from filtering
            ##' @param what return type either M- or beta-values
            ##' @param cutp detection p-value threshold
            ##' @param cutsamples threshold removing samples
            ##' @param cutcpgs threshold removing probes
            ##' @param verbose Default is TRUE
            ##' @param ... optional arguments getBeta
            ##' @return matrix with  M- or beta-values
            ##' @author mvaniterson
            ##' @export
            ##' @importFrom utils data
            reduce <- function(GRset, RGset, what=c("beta", "M", "both"), cutp=0.01, cutsamples=0.1, cutcpgs=0, verbose=TRUE, autoimpute,methyLImp, imputePCA, imputation_method,detectionPs, ...) {

                what <- match.arg(what)

                if(class(GRset) != "GenomicRatioSet")
                    stop("First argument should be an object of class 'GenomicRatioSet' from preprocessFunnorm!")
                if(!(class(RGset) == "RGChannelSet" | class(RGset) == "RGChannelSetExtended"))
                    stop("Second argument should be an object of class RGChannelSet or RGChannelSetExtended from probeFiltering!")

                ##Filter on detection P-value using minfi's detectionP
                if(verbose)
                    cat("Calculate and filter on detection P-value... \n")
                
                pvalmat <- detectionPs(RGset, na.rm=TRUE)
                idPvalmat <- pvalmat > cutp
                idPvalmat[is.na(idPvalmat)] <- TRUE ##set those for which a detection P-value could not be calculate TRUE to be filtered out
               

                if(verbose)
                    cat("On average", round(100*sum(idPvalmat)/prod(dim(idPvalmat)), 2),"% of the CpGs (",nrow(idPvalmat),") have detection P-value above the threshold ",cutp, "\n")

                if(verbose)
                    cat("Transform to ",what,"-values... \n")

                if(what == "M"){
                    matfilt <- minfi::getM(preprocessRaw(RGset), ...)
                    matnorm <- minfi::getM(GRset, ...)
                }
                if (what =="beta"){
                    matfilt <- minfi::getBeta(RGset, ...)
                    matnorm <- minfi::getBeta(GRset, ...)
                }
                else 
                {  if(what=="both") 
                    {
                        matfiltbeta <- minfi::getBeta(RGset, ...)
                        matfiltM <- minfi::getM(preprocessRaw(RGset), ...)
                        matnormbeta <- minfi::getBeta(GRset, ...)
                        matnormM <- minfi::getM(GRset, ...)
                    }
                }    

                ##set max/min M-values to +/- 16
                if(verbose & what=="M")
                    cat("Set +/-Inf to +/-16... \n")

                if(what == "M")
                    matnorm[!is.finite(matnorm)] <-  sign(matnorm[!is.finite(matnorm)])*16
                    
                if(verbose & what=="both")
                    cat("Set +/-Inf to +/-16... \n")    
                if(what == "both")
                    matnormM[!is.finite(matnormM)] <-  sign(matnormM[!is.finite(matnormM)])*16    

                ##set NA from probeFiltering
                ##!!!NOTE orders are not the same
                
                if(verbose & what!="both")
                {
                         cat("On average", round(100*sum(is.na(matfilt))/prod(dim(matfilt)), 2),"% of the probes (",nrow(matfilt),") were set to NA in the probe filtering step! \n")
                    
                    mid <- match(rownames(matfilt), rownames(matnorm))
                  
                    matnorm <- matnorm[mid,]
                    matnorm[is.na(matfilt)] <- NA

                    ##set NA from detectionP
                    ##order seems OK just to be sure
                    mid <- match(rownames(matnorm), rownames(idPvalmat))
                    idPvalmat <- idPvalmat[mid,]
                    matnorm[idPvalmat] <- NA

                    
                    ##calculate failure rates and reduce
                    if(verbose)
                        cat("Calculate failure/success rates and reduce... \n")                    
                    
                    srRows<-rowSums(is.na(matnorm))/ncol(matnorm)
                    
                    srCols<-colSums(is.na(matnorm))/nrow(matnorm)
                    
                    #srCols <- apply(matnorm, 2, function(x) sum(is.na(x))/(length(x)))                    
                    #srRows <- apply(matnorm, 1, function(x) sum(is.na(x))/length(x))
                    

                    if(verbose){
                        cat("Percentage of samples having failure rates below", cutsamples, "is", round(100*sum(srCols < cutsamples)/length(srCols),2),"% \n")
                        cat("Samples having failure rates above", cutsamples, paste(colnames(matnorm[,srCols > cutsamples]),collapse=",")," and will be deleted. ")                        
                        
                        
                        cat("Percentage of CpGs having failure rates below", cutcpgs, "is", round(100*sum(srRows < cutcpgs)/length(srRows),2),"% \n")
                        
                        
                        
                    }
                    
                    matnorms<-matnorm[srRows < cutcpgs,  srCols < cutsamples]
                }
                
                if(verbose & what=="both")
                {    
                   
                      cat("On average", round(100*sum(is.na(matfiltbeta))/prod(dim(matfiltbeta)), 2),"% of the beta probes (",nrow(matfiltbeta),") were set to NA in the probe filtering step! \n")
                     
                    mid <- match(rownames(matfiltbeta), rownames(matnormbeta))
                    matnormbeta <- matnormbeta[mid,]
                    matnormbeta[is.na(matfiltbeta)] <- NA
                                        
                   
                    cat("On average", round(100*sum(is.na(matfiltM))/prod(dim(matfiltM)), 2),"% of the M probes (",nrow(matfiltM),") were set to NA in the probe filtering step! \n")
                    
                    mid <- match(rownames(matfiltM), rownames(matnormM))
                    matnormM <- matnormM[mid,]
                    matnormM[is.na(matfiltM)] <- NA

                    ##set NA from detectionP
                    ##order seems OK just to be sure
                    mid <- match(rownames(matnormbeta), rownames(idPvalmat))
                    idPvalmat <- idPvalmat[mid,]
                    matnormbeta[idPvalmat] <- NA
                    
                    ##set NA from detectionP
                    ##order seems OK just to be sure
                    mid <- match(rownames(matnormM), rownames(idPvalmat))
                    idPvalmat <- idPvalmat[mid,]
                    matnormM[idPvalmat] <- NA

                    
                    ##calculate success rates and reduce
                    if(verbose)
                         cat("Calculate failure/success rates and reduce... \n")

                    #srColsbeta <- apply(matnormbeta, 2, function(x) sum(is.na(x))/(length(x)))
                    srColsbeta<-colSums(is.na(matnormbeta))/nrow(matnormbeta)
                    #srRowsbeta <- apply(matnormbeta, 1, function(x) sum(is.na(x))/length(x))
                    srRowsbeta<-rowSums(is.na(matnormbeta))/ncol(matnormbeta)
                    
                    
                                        
                    #srColsM <- apply(matnormM, 2, function(x) sum(is.na(x))/(length(x)))
                    srColsM<-colSums(is.na(matnormM))/nrow(matnormM)
                    #srRowsM <- apply(matnormM, 1, function(x) sum(is.na(x))/length(x))
                    srRowsM<-rowSums(is.na(matnormM))/ncol(matnormM)

                    if(verbose){
                        cat("Percentage of beta samples having failure rates below", cutsamples, "is", round(100*sum(srColsbeta < cutsamples)/length(srColsbeta),2),"% \n")
                        cat("Percentage of M samples having failure rates below", cutsamples, "is", round(100*sum(srColsM < cutsamples)/length(srColsM),2),"% \n")
                        cat("Percentage of beta CpGs having failure rates below", cutcpgs, "is", round(100*sum(srRowsbeta < cutcpgs)/length(srRowsbeta),2),"% \n")
                        cat("Percentage of M CpGs having failure rates below", cutcpgs, "is", round(100*sum(srRowsM < cutcpgs)/length(srRowsM),2),"% \n")
                        
                        cat("Samples having failure rates above", cutsamples, paste(colnames(matnormM[,srColsM > cutsamples]),collapse=",")," and will be deleted. ")                        
                        
                    }
                    matnormsbeta<-matnormbeta[srRowsbeta < cutcpgs,  srColsbeta < cutsamples]
                    matnormsM<-matnormM[srRowsM < cutcpgs,  srColsM < cutsamples]
                }
                
                
                
                if(autoimpute)
                {  
                       if(verbose){
                       cat("\n  Autoimpute Start") 
                       }
                      

                           if(sum(is.na(idPvalmat)) == 0) 
                          {

                               if(verbose){
                               cat("    No NAs (failed probes) exist in your data set any more, you don't need to do any imputation.")
                               }
                           } 

                           else
                           {   if(verbose){
                               cat("    There are ",sum(is.na(idPvalmat)), " NAs (failed probes) exists in your data set.")
                               cat("    Imputing will be conducted for remain NAs. \n")                               
                               
                               }

                                            

                                if(verbose){

                                    cat("    Doing imputation on matrix.")
                                    }
                                 library(impute)
                                 if (what!="both")
                                 {    
                                     if (imputation_method=="knn")
                                     {  
                                     cat("   Performing impute.knn.")
                                     ### Using sink function to remove messages
                                     zz <- file("ImputeMessage.Rout", open="wt")
                                     sink(zz)
                                     sink(zz, type="message")                                     
                                     
                                     matnorms <- impute.knn(matnorms,k=5)$data 
                                     
                                     sink(type="message")
                                     sink()  
                                     
                                     }
                                      if(imputation_method=="methyLImp")
                                     { cat("   Performing methyLImp")
                                         methyLImps<-methyLImp$methyLImp
                                         plogit<-methyLImp$plogit
                                         inv.plogit<-methyLImp$inv.plogit 
                                         pinvr<-methyLImp$pinvr
                                         
                                         if  (what=="M")
                                         {
                                         matnorms<-methyLImps(t(matnorms), min = -16, max = 16, max.sv = NULL, col.list = NULL)
                                         matnorms<-t(matnorms)
                                         }
                                         else
                                         {
                                         matnorms<-methyLImps(t(matnorms), min = 0, max = 1, max.sv = NULL, col.list = NULL)
                                         matnorms<-t(matnorms)
                                         }
                                     }
                                      if(imputation_method=="imputePCA")
                                     { cat("   Performing imputePCA")
                                       imputePCAs<-imputePCA$imputePCA
                                       FactoMineR.svd.triplet<-imputePCA$FactoMineR.svd.triplet
                                         if  (what=="M")
                                         {                             
                                         matnorms <- imputePCAs(matnorms)

                                             if(class(matnorms)!="matrix")
                                             {
                                             matnorms <- matnorms$completeObs
                                             }
                                         
                                         }
                                         else
                                         {
                                         matnorms <- imputePCAs(matnorms)
                                         
                                         if(class(matnorms)!="matrix")
                                         {
                                             matnorms <- matnorms$completeObs
                                         }    

                                         matnorms[matnorms<0] <- 0
                                         matnorms[matnorms>1] <- 1  
                                         }
                                     }
                                     
                                     
                                 }
                                 else
                                 {                                     
                                     if (imputation_method=="knn")
                                     {                                         
                                     cat("   Performing impute.knn.")  
                                     ### Using sink function to remove messages
                                     zz <- file("ImputeMessage.Rout", open="wt")
                                     sink(zz)
                                     sink(zz, type="message")
                                     
                                     matnormsM <- impute.knn(matnormsM,k=5)$data
                                     matnormsbeta <- impute.knn(matnormsbeta,k=5)$data
                                     
                                     sink(type="message")
                                     sink()  
                                     
                                     }
                                     if(imputation_method=="methyLImp")
                                     {cat("   Performing methyLImp")
                                     methyLImps<-methyLImp$methyLImp
                                     plogit<-methyLImp$plogit
                                     inv.plogit<-methyLImp$inv.plogit 
                                     pinvr<-methyLImp$pinvr                                     
                                     
                                     matnormsM<-methyLImps(t(matnormsM), min = -16, max = 16, max.sv = NULL, col.list = NULL)
                                     matnormsM<-t(matnormsM)
                                     matnormsbeta<-methyLImps(t(matnormsbeta), min = 0, max = 1, max.sv = NULL, col.list = NULL)
                                     matnormsbeta<-t(matnormsbeta)
                                     }
                                     if(imputation_method=="imputePCA")
                                     {cat("   Performing imputePCA")
                                     
                                     imputePCAs<-imputePCA$imputePCA
                                     FactoMineR.svd.triplet<-imputePCA$FactoMineR.svd.triplet
                                     
                                     matnormsbeta <- imputePCAs(matnormsbeta)
                                     
                                     if(class(matnormsbeta)!="matrix")
                                         {
                                         matnormsbeta <- matnormsbeta$completeObs
                                         }
                                     matnormsbeta[matnormsbeta<0] <- 0
                                     matnormsbeta[matnormsbeta>1] <- 1                                     
                                     
                                     matnormsM <- imputePCAs(matnormsM)
                                     if(class(matnormsM)!="matrix")
                                         {
                                         matnormsM <- matnormsM$completeObs
                                         }
                                     }
                                     
                                 }                                                      

                            }                         
                        
                }
                         
                
                if (what!="both")
                {
                  result=list(matr=matnorms)
                }
                else
                {
                    result=list(matrM=matnormsM, matrbeta=matnormsbeta)
                }
                
                return(result)
            }

            R_MAX_MEM_SIZE=memory.limit(size = NA)
            object<-reduce(GRset, RGset, what, cutp, cutsamples, cutcpgs, verbose, autoimpute,methyLImp, imputePCA, imputation_method, detectionPs)
            pheno = pData(GRset)
            if (what!="both")
            {
                keep <- match(colnames(object$matr), rownames(pheno))
                pheno <- pheno[keep,]
                result=list(object$matr, pheno)            
            }
            else
            {
                keep <- match(colnames(object$matrM), rownames(pheno))
                pheno <- pheno[keep,]           
                
                result=list(object, pheno)
            }
            
            return (result) 


     }""")(GRset, RGset, what, detPcut, SampleCutoff, ProbeCutoff, verbose, autoimpute, methyLImp, imputePCA, imputation_method, detectionP)
        
        
        try:
            self.pheno_py=self.ri2py_dataframe(self.pheno, matrix=True)
        except:
            self.pheno_py=self.ri2py_dataframe(self.pheno, matrix=False)  
        
        if (what=="both"):
            self.mval=obj[0]
            self.mval_py = self.ri2py_dataframe(r_dat=self.mval, matrix=False)
            self.beta=obj[1]
            self.beta_py = self.ri2py_dataframe(r_dat=self.beta, matrix=False)                      
            
            return self.beta_py, self.mval_py, self.pheno_py
        
        elif(what=="M"):   
            self.mval=obj
            self.mval_py = self.ri2py_dataframe(r_dat=self.mval, matrix=False)
            return self.mval_py, self.pheno_py
        
        elif(what=="beta"):   
            self.beta=obj
            self.beta_py = self.ri2py_dataframe(r_dat=self.beta, matrix=False)  
            return self.beta_py, self.pheno_py      
    
    def subsetter(self, matrix, pheno_py, disease_list=None, sample_list=None, include=True, phenotype=None, verbose=True ):
           
            if sample_list is not None or disease_list is not None:
                
                mval=self.mval if self.mval else False
                beta=self.beta if self.beta else False 
                
                if sample_list is not None:
                    if include:                        
                        pheno_py = pheno_py[pheno_py['ID'].isin(sample_list)]
                    else:
                        pheno_py = pheno_py[~pheno_py['ID'].isin(sample_list)]
                        
                    
                    self.sub_pheno_py=pheno_py

                    pheno=pandas2ri.py2ri(self.sub_pheno_py)
                    self.sub_pheno=pheno                             

                
                if disease_list is not None:
                    if include: 
                        pheno_py = pheno_py[pheno_py[phenotype].isin(disease_list)]
                    else:
                        pheno_py = pheno_py[~pheno_py[phenotype].isin(disease_list)]
                        
                    self.sub_pheno_py=pheno_py

                    pheno=pandas2ri.py2ri(self.sub_pheno_py)
                    self.sub_pheno=pheno               
                
            
                matrix, self.sub_mval, self.sub_beta = robjects.r("""function(matrix, mval, beta, pheno){

                if (mval!=FALSE)mval<- mval[,colnames(mval) %in% pheno[,c('ID')] ]
                if (beta!=FALSE)beta<- beta[,colnames(beta) %in% pheno[,c('ID')]]
                matrix<- matrix[,colnames(matrix) %in% pheno[,c('ID')]]        
                result=list(matrix, mval,beta)
                return(result)
                }""")(matrix, mval, beta, self.sub_pheno) 
                
                self.sub_mval_py=self.ri2py_dataframe(self.sub_mval, matrix=False)
                self.sub_beta_py=self.ri2py_dataframe(self.sub_beta , matrix=False)   
            
            
            
            elif disease_list is None and sample_list is None:
                if verbose:
                    print('Using the subsetting function without excluding any disease or sample does not make any sense\n...')
                print('Setting and returning original values')
                self.sub_mval_py=self.mval_py
                self.sub_beta_py=self.beta_py
                self.sub_pheno=self.pheno
                self.sub_pheno_py=self.pheno_py
                self.sub_mval=self.mval
                self.sub_beta=self.beta        
   
                
            return self.sub_mval_py, self.sub_beta_py, self.sub_pheno_py, self.sub_mval, self.sub_beta, self.sub_pheno, matrix          
    
    
    def dmp_finder(self,matrix, pheno, disease_list=None, sample_list=None,include=True, phenotype=None, adjust_vars=None, correction_vars=None, 
               useCombat=False, sva=False, number=10000, pvalue=0.05, adjpval=1, save_csv=False, path=None, 
               adjust_method='BH', method='separate', top=None ):
            
            from scipy.special import comb
            if not phenotype:
                print('Please specify a target of interest \n'
                    'These are the available column names:')
                print(pandas2ri.ri2py(robjects.r['as'](pheno,'data.frame')).columns.tolist())
                return None
            try:
                pheno_py=self.ri2py_dataframe(pheno, matrix=True)
            except:
                pheno_py=self.ri2py_dataframe(pheno, matrix=False)       
            
            self.sub_mval_py, self.sub_beta_py, self.sub_pheno_py, self.sub_mval, self.sub_beta, self.sub_pheno, matrix =self.subsetter(matrix, 
              pheno_py, 
              disease_list=disease_list, 
              sample_list=sample_list, 
              phenotype=phenotype,          
              include=include,
              verbose=False )       
                
                
            pheno_py=self.sub_pheno_py  
            pheno=pandas2ri.py2ri(self.sub_pheno_py)
            
            py_array=pheno_py[phenotype].unique()
            r_array=numpy2ri.py2ri(py_array)
            numeric=[]
            categorical=[]
            if adjust_vars:

                if not type(adjust_vars)==list:
                    print('Please provide a list specifying your adjustment variables obtained from the pheno_sheet \n')
                    print(pandas2ri.ri2py(robjects.r['as'](pheno,'data.frame')).columns.tolist())
                    return None, None, None


                for var in adjust_vars:

                    if (pheno_py[var].dtype=='float64') or (pheno_py[var].dtype=='float32') or (pheno_py[var].dtype=='int'):

                        numeric.append(var)
                    else:    

                        categorical.append(var)



            adjust_vars=numeric+categorical  

            numeric=[]
            categorical=[]
            if correction_vars:

                if not type(correction_vars)==list:
                    print('Please provide a list specifying your correction variables obtained from the pheno_sheet \n')
                    print(pandas2ri.ri2py(robjects.r['as'](pheno,'data.frame')).columns.tolist())
                    return None, None, None


                for var in correction_vars:

                    if (pheno_py[var].dtype=='float64') or (pheno_py[var].dtype=='float32') or (pheno_py[var].dtype=='int'):

                        numeric.append(var)
                    else:    

                        categorical.append(var)


            correction_vars=numeric+categorical

            contrasts=dict()
            n = len(py_array)
            k=0
            for i in range(0,n):
                #print(k)
                for j in range(i+1,n):
                    k+=1

                    contrasts[k]= py_array[i]+"-"+py_array[j]


            print('These are the possible pairwise groups for your comparisons \n %s' %(contrasts))
            print('Enter your comparison of choice;- to get all comparisons type "all":')
            coef = input()
            if coef=='all':
                print('you selected all')
            else:
                coef=int(coef)
                max_coef=comb(len(py_array),2)
                #print(int(max_coef))
                if (coef>int(max_coef)):
                    print('The maximum coefficient number for comparisons is %s' %(int(max_coef)))
                    return 0,0   
                else:
                    print('You have chosen this comparison group: %s:%s' %(coef,contrasts[coef]))
                    
                    
            if top:  
                print('Taking only the %s most variable cpgs into account for analysis ' %(top))
                to_do=[self.sub_mval_py, self.sub_beta_py]
                done=[]
                for element in to_do:               
                    if len(element.columns)<=len(element):
                                print('Dataframe needed to be transposed')
                                element=element.transpose()        
                    sorter=element.var(axis=0).sort_values(ascending=False)[:top]
                    element=element.transpose().reindex(sorter.index)#.transpose()

                #sorter=preproidat.mval_py.var(axis=0).sort_values(ascending=False)[:1000]
                #preproidat.mval_py.transpose().reindex(sorter.index).transpose()
                    done.append(element)
                self.sub_mval_py=done[0]
                self.sub_beta_py=done[1]
                
                self.sub_beta=pandas2ri.py2ri(self.sub_beta_py)
                self.sub_mval=pandas2ri.py2ri(self.sub_mval_py)
                
                matrix = robjects.r("""function(matrix, mval ){

                matrix <- matrix[intersect(rownames(mval),rownames(matrix)), ]              
                       
               
                return(matrix)
                }""")(matrix, self.sub_mval) 
                
            elif top==None:
                print('Taking all cpgs into account for analysis')        
                    

            if sva:
                print('You have chosen to include surrogate variable analysis')  
            if len(adjust_vars)!=0:    
                print('You are adjusting for these variables: %s' % (adjust_vars))
            if len(correction_vars)!=0:
                print('You are correcting for these variables: %s' % (correction_vars))

            self.get_annotation()     
            #mod, mod0

            self.dmps, data, dectest, self.vennc, self.fit2= robjects.r("""function ( pheno, phenotype, M,  ann, number, array, coef, pvalue, 
            adjpval, adj_vars,corr_vars, sva,useCombat, adjust_method, method) 
            {

            phen <- as.factor(pheno[,c(phenotype)])


            if (length(adj_vars>0))
            {
                form_mod <- as.formula(paste(" ~0+",phenotype,'+',paste(c(adj_vars),collapse=" + "),sep=""))
                form_mod0 <- as.formula(paste(" ~",paste(c(adj_vars),collapse=" + "),sep="")) 
            }
            else
            {
                form_mod <- as.formula(paste(" ~0+",phenotype))
                form_mod0 <- as.formula(paste(" ~1") )
            }


            mod<-model.matrix(form_mod,data=pheno)
            colnames(mod)<- make.names(colnames(mod))
            rownames(mod)<-make.names(rownames(mod))

            mod0<-model.matrix(form_mod0,data=pheno)
            colnames(mod0)<- make.names(colnames(mod0))
            rownames(mod0)<-make.names(rownames(mod0))        

            modSv=mod

            print(form_mod)
            #print(form_mod0)


            if(length(corr_vars)>0)
            {   cat('Start correcting for batches')
                if(useCombat)
                {
                  cat('Using ComBat function')
                }
                else
                {
                  cat('Using removeBatchEffect function')
                }
                M_2 <- M           

                for(i in 1:length(corr_vars))
                    {   
                        batchname=toString(c(corr_vars[i]))
                        cat("\n<< Start Correcting for ",batchname," >>")           

                        batch=as.factor(pheno[,batchname])
                        #print(batch)
                        if(useCombat)
                        {   
                            M <- ComBat(dat=M,batch=batch,mod=modSv,par.prior=TRUE)
                        }
                        else
                        {                           
                            M = removeBatchEffect(M, batch=batch, design=modSv)
                        }                    
                    }                 
             } 

            ####SVA#####

            if (sva)
            {    cat('creating model for SVA')
                 n.sv = num.sv(M,modSv,method="leek")

                 svobj = sva(M,modSv,mod0,n.sv=n.sv)

                 modSv = cbind(modSv,svobj$sv)
                 colnames(modSv)<- make.names(colnames(modSv))
                 rownames(modSv)<-make.names(rownames(modSv))
            }

            else
            {
                 cat('creating model')
                 #modSv = mod

            }  

            fit <- lmFit(M,modSv)        



            # create a contrast matrix for specific comparisons      
            cat('Creating contrast matrix for experiment')            
            design.pairs <- function(levels) {
                                n <- length(levels)
                                design <- matrix(0,n,choose(n,2))
                                rownames(design) <- levels
                                colnames(design) <- 1:choose(n,2)

                                k <- 0
                                for (i in 1:(n-1))
                                for (j in (i+1):n) {
                                    k <- k+1
                                    design[i,k] <- 1
                                    design[j,k] <- -1
                                    colnames(design)[k] <- paste(phenotype,levels[i],"-",phenotype,levels[j],sep="")
                                    }
                                result<-list(des=design,cols=colnames(design), rows=rownames(design) )    
                                return(result)    
                                #design
                                }         

            matr <-design.pairs(make.names(array))


            x <- matr$cols


            contrast.matrix <-makeContrasts(contrasts=x, levels=modSv)    
            
            fit2 <- contrasts.fit(fit, contrast.matrix)
            fit2 <- eBayes(fit2)
            
           
            # look at the numbers of DM CpGs at FDR < 0.01
            cat('Computing statistics for experiment')   
            
            
            results = decideTests(fit2, method=method, adjust.method=adjust_method,p.value=pvalue)
            dectest<- summary(results)  
            
            
            #ii <- apply(abs(results)>0, 1, all)
            #vennc=vennCounts(results)
            #venndia=vennDiagram(results)
            cat('Aligning annotation')
            #annSub <- ann[match(rownames(M),ann$Name),]
            annSub <- ann[ann$Name %in% rownames(M),]
            #annSub <- ann[match(rownames(fit2$coefficients),ann$Name),]

            if (coef == 'all'){

                            cat('Computing contrasts for experiment')           
                            topresults <- c()             
                            for (i in 1:choose(length(array),2)) {        
                                #print(i)
                                ii <- which(abs(results[,i])==1)
                               
                                top<-topTable(fit2, genelist=annSub, coef=i, number=number,  adjust.method=adjust_method, sort.by="P", 
                                p.value=adjpval)
                               
                                
                                #top<-top[ii,]
                                top<-top[which(top$adj.P.Val < adjpval) ,]
                                #top<-top[which( top$adj.P.Val < adjpval & abs(top$logFC) > 1),]                                
                                
                                topresults[[i]] <- top    

                            }

                            result=list(topresults,matr$des, dectest,results, fit2)
                        }

                        else{

                            cat('Computing contrast for experiment') 
                            ii <- which(abs(results[,coef])==1)
                            top<-topTable(fit2, genelist=annSub, coef=coef,number=number,  adjust.method=adjust_method, sort.by="P", 
                            p.value=adjpval )
                            
                            #top<-top[which(top$adj.P.Val < adjpval) ,]  
                           
                            #top<-top[ii,]
                            #top<-top[which(top$adj.P.Val < adjpval  & abs(top$logFC) > 1) ,]
                            top<-top[which(top$adj.P.Val < adjpval) ,]
                             
                            top<-list(top)
                            result=list(top,matr$des,dectest, results, fit2 )

                        }     


            return(result)                     




           }""")  (pheno, 
                   phenotype,                
                   matrix,                
                   self.annotation, 
                   number, 
                   r_array, 
                   coef,
                   pvalue,
                   adjpval,
                   adjust_vars,               
                   correction_vars,
                   sva,
                   useCombat,
                   adjust_method,
                   method )     

            self.dmps_list=[]
            for elem in self.dmps:
                #print(i)
                if len(elem)>0:
                    #print(elem)
                    #print(self.ri2py_dataframe(elem).index)
                    self.dmps_list.append(pd.DataFrame(pandas2ri.ri2py(elem)))
            self.dmps_df = pd.concat(self.dmps_list).reset_index().drop_duplicates(subset='Name',keep = False)
            
            self.dmps = self.dmps_df['Name'].tolist()
            #self.dmps = list(dict.fromkeys(lists))            
            self.contrastmatrix=self.ri2py_dataframe(data)

            self.dectest=self.ri2py_dataframe(dectest)
            ###to-do
            if save_csv:
                for i,elem in enumerate(self.dmps_list):
                    if path:
                        elem.to_csv(path+'{}.csv'.format(i,contrasts[i]), index=False)
                    else:
                        elem.to_csv(self.idat_dir+'{}.csv'.format(i,contrasts[i]), index=False)

            print('done')      
   

    def get_annotation(self):
        self.annotation=robjects.r("""function (RGset) {
            ann450k <- getAnnotation(RGset)
            return(ann450k)
            
            }""")(self.RGset)
        self.annotation_py=pandas2ri.ri2py(robjects.r['as'](self.annotation,'data.frame'))
        return self.annotation_py  
        
        
        

    def extract_pheno_data(self, methylset=False, grset=False):
        """Extract pheno data from MSet or RGSet, minfi.
        Parameters
        ----------
        methylset
            If MSet has beenn created, set to True, else extract from original RGSet.
        grset
            If GRset has beenn created, set to True, else extract from original RGSet.    """
        
        if grset and not methylset:
            self.pheno = robjects.r("pData")(self.GRset)
        
        elif methylset and not grset:
            self.pheno = robjects.r("pData")(self.MSet) 
        
        else:
            self.pheno = robjects.r("pData")(self.RGset)
        
        return self.pheno

    def extract_manifest(self):
        """Get manifest from RGSet."""
        self.manifest = self.minfi.getManifest(self.RGset)
        return self.manifest     
    

    def output_pheno_beta(self, meffil=False, rset=False):
        """Get pheno and beta dataframe objects stored as attributes for input to MethylationArray object.
        Parameters
        ----------
        meffil
            True if ran meffil pipeline."""
        indexsetter=self.RSet if rset else self.GRset 
            
        self.pheno_py=pandas2ri.ri2py(robjects.r['as'](self.pheno,'data.frame'))
                   
        if not meffil:
            self.beta_py=pd.DataFrame(pandas2ri.ri2py(self.beta_final),index=numpy2ri.ri2py(robjects.r("featureNames")(indexsetter)),columns=numpy2ri.ri2py(robjects.r("sampleNames")(indexsetter))).transpose()
            self.pheno_py['Sample_Name']=np.vectorize(lambda x: x.split('/')[-1])(self.pheno_py['Basename'])
            self.pheno_py = self.pheno_py.set_index('Sample_Name').loc[self.beta_py.index,:]
        else:
            self.beta_py=pd.DataFrame(pandas2ri.ri2py(self.beta_final),index=robjects.r("rownames")(self.beta_final),columns=robjects.r("colnames")(self.beta_final)).transpose()
            print(self.beta_py)
            print(self.beta_py.index)
            print(self.pheno_py)
            self.pheno_py = self.pheno_py.set_index('Sample_Name').loc[self.beta_py.index,:]

    def export_pickle(self, output_pickle, disease=''):
        """Export pheno and beta dataframes to pickle, stored in python dict that can be loaded into MethylationArray
        Parameters
        ----------
        output_pickle
            Where to store MethylationArray.
        disease
            Custom naming scheme for data."""
        output_dict = {}
        if os.path.exists(output_pickle):
            output_dict = pickle.load(open(output_pickle,'rb'))
        output_dict['pheno' if not disease else 'pheno_{}'.format(disease)] = self.pheno_py
        output_dict['beta' if not disease else 'beta_{}'.format(disease)] = self.beta_py
        pickle.dump(output_dict, open(output_pickle,'wb'),protocol=4)

    def export_sql(self, output_db, disease=''):
        """Export pheno and beta dataframes to SQL
        Parameters
        ----------
        output_db
            Where to store data, sqlite db.
        disease
            Custom naming scheme for data."""
        conn = sqlite3.connect(output_db)
        self.pheno_py.to_sql('pheno' if not disease else 'pheno_{}'.format(disease), con=conn, if_exists='replace')
        self.beta_py.to_sql('beta' if not disease else 'beta_{}'.format(disease), con=conn, if_exists='replace')
        conn.close()

    def export_csv(self, output_dir):
        """Export pheno and beta dataframes to CSVs
        Parameters
        ----------
        output_dir
            Where to store csvs."""
        self.pheno_py.to_csv('{}/pheno.csv'.format(output_dir))
        self.beta_py.to_csv('{}/beta.csv'.format(output_dir))

    def to_methyl_array(self,disease=''):
        """Convert results from preprocessing into MethylationArray, and directly return MethylationArray object.
        Parameters
        ----------
        disease
            Custom naming scheme for data."""
        return MethylationArray(self.pheno_py,self.beta_py, disease)
    
        
        
  
