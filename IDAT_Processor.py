import pandas as pd
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr #from utils import importr_ as importr#
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

pandas2ri.activate()
numpy2ri.activate()

class TCGADownloader:
    """Downloads TCGA and GEO IDAT and clinical data"""
    def __init__(self):
        pass

    def download_tcga(self, output_dir):
        """Download TCGA IDATs.
        Parameters
        ----------
        output_dir
            Where to output idat files."""
        tcga = importr("TCGAbiolinks")
        print(tcga)
        robjects.r("""
                   library(TCGAbiolinks)
                   projects <- TCGAbiolinks:::getGDCprojects()$project_id
                   projects <- projects[grepl('^TCGA',projects,perl=T)]
                   match.file.cases.all <- NULL
                   for(proj in projects){
                        print(proj)
                        query <- GDCquery(project = proj,
                                          data.category = "Raw microarray data",
                                          data.type = "Raw intensities",
                                          experimental.strategy = "Methylation array",
                                          legacy = TRUE,
                                          file.type = ".idat",
                                          platform = "Illumina Human Methylation 450")
                        match.file.cases <- getResults(query,cols=c("cases","file_name"))
                        match.file.cases$project <- proj
                        match.file.cases.all <- rbind(match.file.cases.all,match.file.cases)
                        tryCatch(GDCdownload(query, method = "api", files.per.chunk = 20),
                                 error = function(e) GDCdownload(query, method = "client"))
                    }
                    # This will create a map between idat file name, cases (barcode) and project
                    readr::write_tsv(match.file.cases.all, path = "idat_filename_case.txt")
                    # code to move all files to local folder
                    for(file in dir(".",pattern = ".idat", recursive = T)){
                        TCGAbiolinks:::move(file,file.path('%s',basename(file)))
                    }
                   """%output_dir)

    def download_clinical(self, output_dir):
        """Download TCGA Clinical Data.
        Parameters
        ----------
        output_dir
            Where to output clinical data csv."""
        robjects.r("""
                   library(TCGAbiolinks)
                   library(data.table)
                   projects <- TCGAbiolinks:::getGDCprojects()$project_id
                   projects <- projects[grepl('^TCGA',projects,perl=T)]
                   match.file.cases.all <- NULL
                   data <- list()
                   for(n in 1:length(projects)){
                        proj <- projects[n]
                        clin.query <- GDCquery_clinic(project = proj,
                                          type='clinical', save.csv=F)
                        data[[length(data)+1]] = clin.query
                    }
                    df <- rbindlist(data)
                    write.csv(df, file=file.path('%s','clinical_info.csv'))
                   """%output_dir)

    def download_geo(self, query, output_dir):
        """Download GEO IDATs.
        Parameters
        ----------
        query
            GEO accession number to query, must be 450k/850k.
        output_dir
            Output directory to store idats and clinical information csv"""
        base=importr('base')
        geo = importr("GEOquery")
        geo.getGEOSuppFiles(query)
        raw_tar_file = "{0}/{0}_RAW.tar".format(query)
        if not os.path.exists(raw_tar_file):
            print("Warning: GEO file {} not downloaded. Check accession on GEO and make sure there is this file available.".format(raw_tar_file))
        robjects.r["untar"](raw_tar_file, exdir = "{}/idat".format(query), tar='internal')
        idatFiles = robjects.r('list.files("{}/idat", pattern = "idat.gz$", full = TRUE)'.format(query))
        robjects.r["sapply"](idatFiles, robjects.r["gunzip"], overwrite = True)
        subprocess.call('mv {}/idat/*.idat {}/'.format(query, output_dir),shell=True)
        self.download_geo_clinical_info(query,output_dir)

    def download_geo_clinical_info(self, query, output_dir, other_output_fname=''):
        geo = importr("GEOquery")
        robjects.r("write.csv(as(getGEO('{}')[[1]]@phenoData@data,'data.frame'),'{}')".format(query,'{}/{}_clinical_info.csv'.format(output_dir,query) if not other_output_fname else other_output_fname))
        #pandas2ri.ri2py(robjects.r('as')(robjects.r("getGEO('{}')[[1]]@phenoData@data".format(query)),'data.frame')).to_csv('{}/{}_clinical_info.csv'.format(output_dir,query) if not other_output_fname else other_output_fname)# ,GSEMatrix = FALSE

class PreProcessPhenoData:
    """Class that will manipute phenotype samplesheet before preprocessing of IDATs.
    pheno_sheet
        Location of clinical info csv.
    idat_dir
        Location of idats
    header_line
        Where to start reading clinical csv"""
    def __init__(self, pheno_sheet, idat_dir, header_line=0):
        self.xlsx = True if pheno_sheet.endswith('.xlsx') or pheno_sheet.endswith('.xls') else False
        if self.xlsx:
            self.pheno_sheet = pd.read_excel(pheno_sheet,header=header_line)
        else:
            self.pheno_sheet = pd.read_csv(pheno_sheet, header=header_line)
        self.idat_dir = idat_dir

    def format_geo(self, disease_class_column="methylation class:ch1", include_columns={}):
        """Format clinical sheets if downloaded geo idats.
        Parameters
        ----------
        disease_class_column
            Disease column of clinical info csv.
        include_columns
            Dictionary specifying other columns to include, and new names to assign them to."""
        idats = glob.glob("{}/*.idat".format(self.idat_dir))
        idat_basenames = np.unique(np.vectorize(lambda x: basename(x).replace('_Grn.idat','').replace('_Red.idat',''))(idats)) # '_'.join(x.split('/')[-1].split('_')[:3])
        idat_geo_map = dict(zip(np.vectorize(lambda x: basename(x).split('_')[0])(idat_basenames),np.array(idat_basenames)))
        self.pheno_sheet['Basename'] = self.pheno_sheet['geo_accession'].map(idat_geo_map)
        self.pheno_sheet = self.pheno_sheet[self.pheno_sheet['Basename'].isin(idat_basenames)]
        self.pheno_sheet.loc[:,'Basename'] = self.pheno_sheet['Basename'].map(lambda x: self.idat_dir+x)
        col_dict = {'geo_accession':'AccNum',disease_class_column:'disease'}
        col_dict.update(include_columns)
        self.pheno_sheet = self.pheno_sheet[['Basename', 'geo_accession',disease_class_column]+(list(include_columns.keys()) if include_columns else [])].rename(columns=col_dict)

    def format_tcga(self, mapping_file="idat_filename_case.txt"):
        """Format clinical sheets if downloaded tcga idats.
        Parameters
        ----------
        mapping_file
            Maps uuids to proper tcga sample names, should be downloaded with tcga clinical information."""
        def decide_case_control(barcode):
            case_control_num = int(barcode.split('-')[3][:2])
            if case_control_num < 10:
                return 'case'
            elif case_control_num < 20:
                return 'normal'
            else:
                return 'control'
            return 0
        idats = glob.glob("{}/*.idat".format(self.idat_dir))
        barcode_mappings = pd.read_csv(mapping_file,sep='\t')
        barcode_mappings['barcodes'] = np.vectorize(lambda x: '-'.join(x.split('-')[:3]))(barcode_mappings['cases'])
        barcode_mappings['idats'] = barcode_mappings['file_name'].map(lambda x: x[:x.rfind('_')])
        barcode_mappings_d1 = dict(barcode_mappings[['barcodes','idats']].values.tolist())
        barcode_mappings['case_controls']= barcode_mappings['cases'].map(decide_case_control)
        barcode_mappings_d2 = dict(barcode_mappings[['barcodes','case_controls']].values.tolist())
        self.pheno_sheet['Basename'] = self.pheno_sheet['bcr_patient_barcode'].map(barcode_mappings_d1)
        self.pheno_sheet['case_control'] = self.pheno_sheet['bcr_patient_barcode'].map(barcode_mappings_d2)
        idat_basenames = np.unique(np.vectorize(lambda x: '_'.join(x.split('/')[-1].split('_')[:2]))(idats))
        self.pheno_sheet = self.pheno_sheet[self.pheno_sheet['Basename'].isin(idat_basenames)]
        self.pheno_sheet.loc[:,['Basename']] = self.pheno_sheet['Basename'].map(lambda x: self.idat_dir+x)
        self.pheno_sheet = self.pheno_sheet[['Basename', 'bcr_patient_barcode', 'disease', 'tumor_stage', 'vital_status', 'age_at_diagnosis', 'gender', 'race', 'ethnicity','case_control']].rename(columns={'tumor_stage':'stage','bcr_patient_barcode':'PatientID','vital_status':'vital','gender':'Sex','age_at_diagnosis':'age'})

   
    def format_custom(self, basename_col, disease_class_column, include_columns={}):
        """Custom format clinical sheet if user supplied idats.
        Parameters
        ----------
        basename_col
            Column name of sample names.
        disease_class_column
            Disease column of clinical info csv.
        include_columns
            Dictionary specifying other columns to include, and new names to assign them to.
        """
        #idats = glob.glob("{}/*.idat".format(self.idat_dir))
        idats = [os.path.join(root, name) for root, dirs, files in os.walk(self.idat_dir) for name in files if name.endswith((".idat"))]   
        idat_basenames = np.unique(np.vectorize(lambda x: '_'.join(x.strip(self.idat_dir).split('_')[:-1]))(idats))        
        idat_batchnames = np.vectorize(lambda x: 'b'+'_'.join(x.split('_')[:-1]).split('/')[-1])(idat_basenames.tolist())
        #np.unique(np.vectorize(lambda x: '_'.join(x.split('/')[-1].split('_')[:-1]))(idats))
        idat_count_underscores = np.vectorize(lambda x: x.count('_'))(idat_basenames)
        #print(basename_col)
        self.pheno_sheet['Basename'] = self.pheno_sheet[basename_col]
        basename_count_underscores = np.vectorize(lambda x: x.count('_'))(self.pheno_sheet['Basename'])
        min_underscores=min(np.hstack([idat_count_underscores,basename_count_underscores]))
        basic_basename_fn = np.vectorize(lambda x: '_'.join(x.split('/')[-1].split('_')[-min_underscores-1:]))
        basic_basename=dict(zip(basic_basename_fn(self.pheno_sheet['Basename']),self.pheno_sheet['Basename'].values))
        batchnames= dict(zip(idat_basenames, idat_batchnames))
        
        #print(batchnames)
        basic_idat=dict(zip(basic_basename_fn(idat_basenames),idat_basenames))
        basic_batchnum=dict(zip(basic_basename_fn(idat_basenames),idat_batchnames))
        complete_mapping={basic_basename[basename]:basic_idat[basename] for basename in basic_basename}
        complete_mapp={basic_basename[basename]:basic_batchnum[basename] for basename in basic_basename}
        self.pheno_sheet.loc[:,'Batchnum']=self.pheno_sheet['Basename'].map(complete_mapp)
        #print(complete_mapping)
        self.pheno_sheet.loc[:,'Basename']=self.pheno_sheet['Basename'].map(complete_mapping).map(lambda x: self.idat_dir+x)
        self.pheno_sheet['disease'] = self.pheno_sheet[disease_class_column.replace("'",'')]
        #self.pheno_sheet['Batchnum'] = idat_batchnames
        self.pheno_sheet = self.pheno_sheet[np.unique(['Basename', 'disease', 'Batchnum']+list(include_columns.keys()))].rename(columns=include_columns)
        if ('sex' in include_columns.keys()) or ('sex' in include_columns.values()) :
            
            if any(item not in self.pheno_sheet['sex'].unique() for item in ['M', 'F']):
                print('Check whether your input data at the "sex"-column contains M for male and/ or F for Female')
            
        return self.pheno_sheet
    
    
    def format_sex_column(self, male_name=None, female_name=None):
        if any(item ==None for item in [male_name, female_name]):
            print('Please specifiy your sexes right')
            return
        self.pheno_sheet['sex']=self.pheno_sheet['sex'].map({male_name:'M', female_name:'F'})
        return self.pheno_sheet
    
    
    def merge(self, other_formatted_sheet, use_second_sheet_disease=True, no_disease_merge=False):
        """Merge multiple PreProcessPhenoData objects, merge their dataframes to accept more than one saplesheet/dataset or add more pheno info.
        Parameters
        ----------
        other_formatted_sheet
            Other PreProcessPhenoData to merge.
        use_second_sheet_disease
            Change disease column to that of second sheet instead of first.
        no_disease_merge
            Keep both disease columns from both sheets.
        """
        disease_dict = {False:'disease_x',True:'disease_y'}
        self.pheno_sheet = self.pheno_sheet.merge(other_formatted_sheet.pheno_sheet,how='inner', on='Basename')
        if not no_disease_merge:
            self.pheno_sheet['disease'] = self.pheno_sheet[disease_dict[use_second_sheet_disease]]
        cols=list(self.pheno_sheet)
        self.pheno_sheet = self.pheno_sheet[[col for col in cols if col!='Unnamed: 0_x' and col!='Unnamed: 0_y' and col!=disease_dict[use_second_sheet_disease]]]

    def concat(self, other_formatted_sheet):
        """Concat multiple PreProcessPhenoData objects, concat their dataframes to accept more than one smaplesheet/dataset.
        Parameters
        ----------
        other_formatted_sheet
            Other PreProcessPhenoData to concat.
        """
        self.pheno_sheet=pd.concat([self.pheno_sheet,other_formatted_sheet.pheno_sheet],join='inner').reset_index(drop=True)
        self.pheno_sheet=self.pheno_sheet[[col for col in list(self.pheno_sheet) if not col.startswith('Unnamed:')]]

    def export(self, output_sheet_name):
        """Export pheno data to csv after done with manipulation.
        Parameters
        ----------
        output_sheet_name
            Output csv name.
        """
        self.pheno_sheet.to_csv(output_sheet_name, index=False)
        print("Next Step: Please move all other sample sheets out of this directory.")

    def split_key(self, key, subtype_delimiter):
        """Split pheno column by key, with subtype delimiter, eg. entry S1,s2 -> S1 with delimiter ",".
        Parameters
        ----------
        key
            Pheno column name.
        subtype_delimiter
            Subtype delimiter to split on.
        """
        new_key = '{}_only'.format(key)
        self.pheno_sheet[new_key] = self.pheno_sheet[key].map(lambda x: x.split(subtype_delimiter)[0])
        return new_key

    def get_categorical_distribution(self, key, disease_only=False, subtype_delimiter=','):
        """Print categorical distribution, counts for each unique value in phenotype column.
        Parameters
        ----------
        key
            Phenotype Column.
        disease_only
            Whether to split phenotype column entries by delimiter.
        subtype_delimiter
            Subtype delimiter to split on.
        """
        if type(key) == type('string'):
            if disease_only:
                key=self.split_key(key,subtype_delimiter)
            return Counter(self.pheno_sheet[key])
        else:
            cols=self.pheno_sheet[list(key)].astype(str)
            cols=reduce(lambda a,b: a+'_'+b,[cols.iloc[:,i] for i in range(cols.shape[1])])
            return Counter(cols)

    def remove_diseases(self,exclude_disease_list, low_count=False, disease_only=False,subtype_delimiter=False):
        """Remove samples with certain diseases from disease column.
        Parameters
        ----------
        exclude_disease_list
            List containing diseases to remove.
        low_count
            Remove samples that have less than x disease occurances in column.
        disease_only
            Whether to split phenotype column entries by delimiter.
        subtype_delimiter
            Subtype delimiter to split on.
        """
        if low_count:
            low_count = int(low_count)
            cat_dist = self.get_categorical_distribution('disease', disease_only,subtype_delimiter).items()
            count_diseases=pd.DataFrame(list(cat_dist),columns=['disease','count'])
            count_diseases.loc[:,'count'] = count_diseases.loc[:,'count'].astype(int)
            exclude_diseases_more=count_diseases.loc[count_diseases['count'].values<low_count,'disease']
            exclude_diseases_more=exclude_diseases_more.unique().tolist()
            if disease_only:
                exclude_diseases_more=self.pheno_sheet.loc[self.pheno_sheet['disease_only'].isin(exclude_diseases_more),'disease'].unique().tolist()
        else:
            exclude_diseases_more=[]
        self.pheno_sheet = self.pheno_sheet[~self.pheno_sheet['disease'].isin(exclude_disease_list+exclude_diseases_more)]
        
    def create_ID(self, disease='disease', sample='Sample'):
        # give the samples descriptive names
        check_list=[disease, sample]
        for val in check_list: 
            if val not in self.pheno_sheet.columns.tolist():
                print('The pheno sheet does not contain a '+disease+' or '+sample+' column you specified \n'
                     'These are the available column names:')
                print(self.pheno_sheet.columns.tolist())
                return
            
        self.pheno_sheet['ID']=self.pheno_sheet[disease]+'.'+self.pheno_sheet[sample]
        
        return self.pheno_sheet   



class PreProcessIDATs:
    """Class that will preprocess IDATs using R pipelines.
    idat_dir
        Location of idats or samplesheet csv.
    minfi
        Rpy2 importr minfi library, default to None will load through rpy2
    enmix
        Rpy2 importr enmix library, default to None will load through rpy2
    base
        Rpy2 importr base library, default to None will load through rpy2
    meffil
        Rpy2 importr meffil library, default to None will load through rpy2"""
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
        
        
        
    def load_idats(self, use_cache=False, rename_samples=True, parallel=True, nworkers=2, verbose=True, extended=True):
        """For minfi pipeline, load IDATs from specified idat_dir."""
        
        #from methylprep.files import get_sample_sheet, find_sample_sheet
        #sample_sheet = get_sample_sheet(data_dir, filepath=sample_sheet_filepath)
        
        cache_storage_path = os.path.join(self.idat_dir,'RGSet.rds')
        self.pheno = self.minfi.read_metharray_sheet(self.idat_dir)
        
        if use_cache:
            self.RGset=robjects.r('readRDS')(cache_storage_path)            
        else:
            if parallel:
                self.RGset=self.load_idats_parallel(targets=self.pheno,verbose=verbose,extended=extended,nworkers=nworkers)
                
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
        #return pd.DataFrame(pandas2ri.ri2py(self.pheno))    
    
    
    def load_idats_parallel(self, targets, verbose=True, extended=True, nworkers=2 ):
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
        # give the samples descriptive names
        #check_list=[disease, sample]
        #for val in check_list: 
        #    if val not in pandas2ri.ri2py(robjects.r['as'](self.pheno,'data.frame')).columns.tolist():
        #        print('The pheno sheet does not contain a '+disease+' or '+sample+' column you specified \n'
        #             'These are the available column names:')
        #        print(pandas2ri.ri2py(robjects.r['as'](self.pheno,'data.frame')).columns.tolist())
        #        return
        
        self.RGset= robjects.r("""function (rgSet,pheno) {
            #pheno$ID <- paste(disease, sample, sep=".")
            sampleNames(rgSet) <- pheno$ID
            result=list(rg=rgSet)
            return(result$rg)
            }""")(self.RGset, self.pheno)
        #pandas2ri.ri2py(self.pheno)[disease], pandas2ri.ri2py(self.pheno)[sample] 
        #pd.DataFrame(pandas2ri.ri2py(self.pheno))
        #return self.RGset 
        
        
    #### QC ##################################################  
    
    
    
    def plt_mds(self, dataframe=None, pheno=None,  top=1000, n_components=2, group='disease', components=(0,1)): 
        if dataframe is None: 
            beta=self.minfi.getBeta(self.RGset)
            dataframe=pd.DataFrame(pandas2ri.ri2py(beta),
                       index=numpy2ri.ri2py(robjects.r("rownames")(beta)),
                       columns=numpy2ri.ri2py(robjects.r("colnames")(beta))).transpose()
        if pheno is None:
            pheno=pandas2ri.ri2py(self.pheno)
            
        import matplotlib.pyplot as plt
        import matplotlib.cm as cm
        from sklearn.manifold import MDS
        import seaborn as sns
        from seaborn import cubehelix_palette
        sns.set()
        #print(len(dataframe.columns))
        #print(len(dataframe))
        if len(dataframe.columns)<=len(dataframe):
            #print('Dataframe needed to be transposed')
            dataframe=dataframe.transpose()        
        sorter=dataframe.var(axis=0).sort_values(ascending=False)[:top]
        dataframe=dataframe.transpose().reindex(sorter.index).transpose()
        from sklearn.impute import SimpleImputer
        imp = SimpleImputer(missing_values=np.nan, strategy='mean')
        imp.fit(dataframe)
        dataframe=pd.DataFrame(imp.transform(dataframe),index = dataframe.index.to_numpy(),
                  columns=dataframe.columns)
        embedding = MDS(n_components=n_components)
        X_transformed = embedding.fit_transform(X=dataframe.to_numpy(), y=pheno[group])
        X_transformed=pd.DataFrame(X_transformed)
        X_transformed[group]=pheno[group].reset_index(drop=True)           
        
        fig, ax = plt.subplots()    
        sns.scatterplot(components[0],components[1],hue=group, cmap=cubehelix_palette(as_cmap=True),data=X_transformed, ax=ax )
            
        #ax.legend(categories.unique())
        ax.set_xlabel('Principal Component %s'  % (components[0]+1))
        ax.set_ylabel('Principal Component %s'  % (components[1]+1))
        ax.set_title('MDS Plot')
        plt.show()       
    
    
    def getQC(self, addQC=False, phenotype=None):
        # give the samples descriptive names
        #check_list=[disease, sample]
        #for val in check_list: 
        if phenotype is not None:
            if phenotype not in pandas2ri.ri2py(robjects.r['as'](self.pheno,'data.frame')).columns.tolist():
                print('The pheno sheet does not contain a '+phenotype+' column you specified \n'
                         'These are the available column names:')
                print(pandas2ri.ri2py(robjects.r['as'](self.pheno,'data.frame')).columns.tolist())
                return None, None, None, None
        
        MSet = self.minfi.preprocessRaw(self.RGset)
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
            pheno=pd.concat([pandas2ri.ri2py(robjects.r['as'](self.pheno,'data.frame')), data_py], axis=1, sort=False)
            self.pheno=pandas2ri.py2ri(pheno)
        
        if phenotype is not None:
            try:       
                datfr[phenotype]=pd.DataFrame(pandas2ri.ri2py(self.pheno))[phenotype].to_numpy()
            except:
                datfr[phenotype]=pd.DataFrame(pandas2ri.ri2py(robjects.r['as'](self.pheno,'data.frame'))[phenotype]).to_numpy()
        
        self.QC_df=datfr
        if phenotype is not None:
            return MSet, qc, datfr, phenotype
        else:
            return MSet, qc, datfr
                   
       
    def plt_mu(self, hue=None, thresh=None):              
        
        if hue is not None:
            _,_, datfr, phenotype =self.getQC(addQC=False,phenotype=hue)   
            if phenotype is not None:
                hue=datfr[phenotype]
            else:
                return
           
        else:
            _,_, datfr =self.getQC(addQC=False)   
         
        
        hymin=datfr['mMed'].min()
        vymin=datfr['uMed'].min()
        
        import seaborn as sns
        fig, ax = plt.subplots()  
        sns.scatterplot(x=datfr['mMed'],y=datfr['uMed'], hue=hue, ax=ax) 
        if thresh is not None:
            ax.vlines(x=thresh, ymin=vymin, ymax=thresh, ls='--')
            ax.hlines(y=thresh, xmin=hymin, xmax=thresh, ls='--') 
        
    def plt_meandetP(self, detPcut=0.01, cutsamples=0.1, cutcpgs=0, log_scale=True, plot='all' ):
        
        #dataframe=self.detectionP()        
               
        #if len(dataframe)<=len(dataframe.columns):
        #    print('Dataframe needed to be transposed')
        #    dataframe=dataframe.transpose()  
        ### all, allsamples, badsamples, goodsamples    
            
        detP, detP_keep, detP_fails = robjects.r("""function (rgset, detPcut, cutsamples, cutcpgs) { 
            detP <- detectionP(rgset)
            numfail <- matrix(colMeans(detP >= detPcut))
            rownames(numfail) <- colnames(detP)
            colnames(numfail) <- "Failed CpG Fraction."
            #print(numfail)
            RemainSample <- which(numfail < cutsamples)
            FailSample <- which(numfail > cutsamples)

            if(any(numfail >= cutsamples))
            {   
                rgset_fails <- rgset[,FailSample]
                rgset_keep <- rgset[,RemainSample]
                
            }

           #RemainProbe <- rowSums(detP > detPcut) <= cutcpgs * length(RemainSample)    
           #FailProbe <- rowSums(detP < detPcut) <= cutcpgs * length(FailSample)

           #rgset_fails <- rgset_fails[FailProbe,]    
           detP_fails <- detectionP(rgset_fails) 
           
           #rgset_keep <- rgset_keep[RemainProbe,] 
           detP_keep <- detectionP(rgset_keep) 
           
           result=list(detP,detP_keep, detP_fails)
           return(result)
           
        }""")(self.RGset,  detPcut, cutsamples, cutcpgs)       
        
        
        detP_py=self.ri2py_dataframe(detP, matrix=False)
        detP_keep_py=self.ri2py_dataframe(detP_keep, matrix=False)
        detP_fails_py=self.ri2py_dataframe(detP_fails, matrix=False)
        
        if plot== 'all':    
            
            
            fig, (ax1, ax2, ax3) = plt.subplots(nrows=3, ncols=1, figsize=(20,5))    
            ax1 = sns.barplot(x=detP_py.columns, y=detP_py.mean(axis=0).to_numpy(), ax=ax1)
            ax2 = sns.barplot(x=detP_keep_py.columns, y=detP_keep_py.mean(axis=0).to_numpy(), ax=ax2)
            ax3 = sns.barplot(x=detP_fails_py.columns, y=detP_fails_py.mean(axis=0).to_numpy(), ax=ax3)
            if log_scale:
                ax1.set_yscale('log')
                ax2.set_yscale('log')
                ax3.set_yscale('log')
            
            ax1.axhline(detPcut, ls='--')
            ax2.axhline(detPcut, ls='--')
            ax3.axhline(detPcut, ls='--')
            ax1.set_title('All Samples')   
            ax2.set_title('Good Samples')   
            ax3.set_title('Bad Samples') 
            ax1.set_xticklabels(ax1.get_xticklabels(), rotation=45, horizontalalignment='right', fontsize='small')
            ax2.set_xticklabels(ax2.get_xticklabels(), rotation=45, horizontalalignment='right', fontsize='small')
            ax3.set_xticklabels(ax3.get_xticklabels(), rotation=45, horizontalalignment='right', fontsize='small')
            plt.tight_layout()
            plt.show()
            return detP_py
        
        if plot== 'allsamples':
            
            dataframe=detP_py
            if len(dataframe)<=len(dataframe.columns):
                print('Dataframe needed to be transposed')
                dataframe=dataframe.transpose()  
            fig, ax = plt.subplots(figsize=(20, 7))
            fig.subplots_adjust=0.85
            ax = sns.barplot(x=dataframe.columns, y=dataframe.mean(axis=0).to_numpy())
            if log_scale:
                ax.set_yscale('log')
            ax.axhline(detPcut, ls='--')    
            ax.set_title('All Samples') 
            ax.set_xticklabels(ax.get_xticklabels(), rotation=45, horizontalalignment='right', fontsize='small')          
            plt.tight_layout()
            plt.show()  
            return
            
        if plot== 'badsamples': 
            dataframe=detP_fails_py
            if len(dataframe)<=len(dataframe.columns):
                print('Dataframe needed to be transposed')
                dataframe=dataframe.transpose()  
            fig, ax = plt.subplots()     
            ax = sns.barplot(x=dataframe.columns, y=dataframe.mean(axis=0).to_numpy())
            if log_scale:
                ax.set_yscale('log')
            ax.axhline(detPcut, ls='--')    
            ax.set_title('Bad Samples') 
            ax.set_xticklabels(ax.get_xticklabels(), rotation=45, horizontalalignment='right', fontsize='small')
            plt.tight_layout()
            plt.show()   
            return
            
        if plot== 'goodsamples':   
            dataframe=detP_keep_py
            if len(dataframe)<=len(dataframe.columns):
                print('Dataframe needed to be transposed')
                dataframe=dataframe.transpose()  
            fig, ax = plt.subplots()     
            ax = sns.barplot(x=dataframe.columns, y=dataframe.mean(axis=0).to_numpy())
            if log_scale:
                ax.set_yscale('log')
            ax.axhline(detPcut, ls='--')    
            ax.set_title('Good Samples')
            ax.set_xticklabels(ax.get_xticklabels(), rotation=45, horizontalalignment='right', fontsize='small')
            plt.tight_layout()
            plt.show()   
            return
            
        else:
            print('Please specify which plots you would want: \n'
            'either "all" "allsamples" "badsamples" or "goodsamples" may be specified')
            return
    
    def plt_curves(self, matrix=None, pheno=None, variant="grouped", group_by=None, nrows=1):
        
        if pheno is None:
            pheno=self.pheno_orig_py
        import matplotlib.pyplot as plt
        import seaborn as sns
        import copy
        from seaborn import cubehelix_palette
        fig, ax = plt.subplots()         
        from collections import OrderedDict
        import matplotlib.pyplot as plt

        cmap=sns.color_palette('muted') 
        
        if type(matrix)==dict:
            #print('dict')
            if (variant=="grouped"):
                
                for i, (key, matr) in enumerate(matrix.items()):
                    #print(len(betas))  
                    pheno=pheno[pheno['ID'].isin(matr.columns.to_numpy())]
                    group_list=pheno[group_by].unique().tolist()
                    matr=copy.deepcopy(matr.transpose())
                    c=len(matrix)/nrows
                    if nrows >= 1:
                        ind=i+1                        
                    plt.subplot(nrows, c, ind)
                    if nrows > 1:
                        if ind==nrows:
                            ind=0                        
                    for val in group_list:
                        #print(val)
                        sns.distplot(matr[np.array(pheno[group_by]==val)].mean(), label=val, hist=False)
                    plt.xlabel('curve-values')
                    plt.ylabel('density')
                    plt.title(key)
                    handles, labels = plt.gca().get_legend_handles_labels()
                    by_label = OrderedDict(zip(labels, handles))
                    plt.legend(by_label.values(), by_label.keys())
                plt.tight_layout()    
                plt.show()
            elif(variant=="single"):
                for i, (key, matr) in enumerate(matrix.items()):
                    #print(len(betas))
                    pheno=pheno[pheno['ID'].isin(matr.columns.to_numpy())]
                    matr=copy.deepcopy(matr)
                    c=len(matrix)/nrows
                    if nrows >= 1:
                        ind=i+1                        
                    plt.subplot(nrows, c, ind)
                    if nrows > 1:
                        if ind==nrows:
                            ind=0                      
                    #print(matr)
                    #label=val,    
                    #print(len(matr.columns))
                    colors = ['#2300A8', '#00A658','#1f77b4', '#ff7f0e', 
                       '#d62728', '#9467bd', '#8c564b', 
                      '#e377c2', '#7f7f7f', '#bcbd22', '#17becf','yellow','#070d0d', '#8ffe09']

                    colordict={}
                    for items, values in zip(pheno[group_by].unique().tolist(), colors):
                        colordict[items]=values
                                        
                    for col in matr.columns:
                        #print(pheno[group_by][pheno['ID']==col].to_numpy())
                        label=pheno[group_by][pheno['ID']==col].to_numpy()[0]
                        #print(matr[col])
                        if len(matr.columns) <= 30:
                            sns.distplot(
                                matr[col], hist=False, rug=False, label=label,
                                color=colordict[label])
                        else:
                            sns.distplot(
                                matr[col], hist=False, rug=False,label=label, color=colordict[label])
                    #sns.distplot(matr.to_numpy(),  hist=False)
                    plt.xlabel('curve-values')
                    plt.ylabel('density')
                    plt.title(key)
                    handles, labels = plt.gca().get_legend_handles_labels()
                    by_label = OrderedDict(zip(labels, handles))
                    plt.legend(by_label.values(), by_label.keys())
                #plt.legend()   
                plt.tight_layout()
                plt.show()
                
                
        else:
            print('Please provide dataframes in key value dict manner, where the key is the label of the plot and the values is the ' 
                  'dataframe')
              
    
    def getM(self,Object=None):
        if Object:
            Object=Object
        else:
            Object=self.RGset
            
        mvals_py, mvals=self.compute_mvals(Object)    
        return mvals_py, mvals
        #mat <- getM(preprocessRaw(RGset), ...)
        
    
    
    def getBeta(self,Object=None):
        if Object:
            Object=Object
        else:
            Object=self.RGset
            
        betas_py,betas=self.compute_betas(Object) 
        return betas_py, betas
        #mat <- getBeta(RGset)
    
    
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
        
        return mval_py ,mval 
          
    def detectionP(self, RGset=None):
        # calculate the detection p-values
        self.detP =robjects.r("detectionP")(RGset) 
        #detP_py=pd.DataFrame(pandas2ri.ri2py(self.detP),
        #        index=numpy2ri.ri2py(robjects.r("rownames")(self.detP)),
        #        columns=numpy2ri.ri2py(robjects.r("colnames")(self.detP))).transpose()
        detP_py=self.ri2py_dataframe(self.detP, matrix=False)
        return detP_py   
    
    def ri2py_dataframe(self, r_dat, matrix=False):
        if matrix:
            py_dat=pd.DataFrame(pandas2ri.ri2py(robjects.r['as'](r_dat,'data.frame')))
        elif not matrix:
            py_dat=pd.DataFrame(pandas2ri.ri2py(r_dat),
                index=numpy2ri.ri2py(robjects.r("rownames")(r_dat)),
                columns=numpy2ri.ri2py(robjects.r("colnames")(r_dat)))
        return py_dat
    #index=numpy2ri.ri2py(robjects.r("rownames")(r_dat)),
    #                        columns=numpy2ri.ri2py(robjects.r("colnames")(r_dat))
    
    def plot_qc_minfi_report(self, output_dir, diseasegroup='disease', samplenames='ID'):
        check_list=[diseasegroup, samplenames]
        for val in check_list: 
            if val not in pandas2ri.ri2py(robjects.r['as'](self.pheno,'data.frame')).columns.tolist():
                print('The pheno sheet does not contain a '+diseasegroup+' or '+samplenames+' column you specified \n'
                     'These are the available column names:')
                print(pandas2ri.ri2py(robjects.r['as'](self.pheno,'data.frame')).columns.tolist())
                return
        
        self.minfi.qcReport(self.RGset, 
                            sampNames=pandas2ri.ri2py(self.pheno)[samplenames], 
                            sampGroups=pandas2ri.ri2py(self.pheno)[diseasegroup],                            
                            pdf = "{}/qcReport.pdf".format(output_dir))
        
        
    def remove_badsamples(self, badSampleCutoff=10,rm_badsamples=False, detPcut=0.01, cutsamples=0.1, addQC=False, verbose=True, RGset=None):
        ''' 
        excluding samples if avarage of unmeth and meth values is below 10
        and if detectionP probe wise threshold is below 0,01(detPcut) 
        and fraction of affected probes for a given sample is above 0.1 (cutsamples)
             
        '''
        
        if RGset:
            RGset=RGset
        else:
            RGset=self.RGset
            
        MSet,qc, datfr = self.getQC(addQC=addQC)
        #print(qc)    
        self.RGset, self.pheno = robjects.r("""function (rgset, qc, badSampleCutoff, detPcut, cutsamples, verbose, rm_badsamples) {
            
            #### meth and unmeth filtering ##########
            #MSet<-preprocessRaw(rgset)
            #qc <- getQC(MSet)
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
            
            if(verbose){   
                cat("\n  Filtering Detect P value Start\n")
                
            }
           
            detP <- detectionP(rgset)
            numfail <- matrix(colMeans(detP >= detPcut))
            rownames(numfail) <- colnames(detP)
            colnames(numfail) <- "Failed CpG Fraction."
            
            RemainSample <- which(numfail < cutsamples)
            print(numfail)
            if(any(numfail >= cutsamples))
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
           
           if('filenames' %in% names(colData(rgset)))
               {
               colData(rgset) <- subset(colData(rgset), select = -filenames )
               }
           pheno = pData(rgset)
           result=list(rgset, pheno)
           return(result)                       
            
            }""")(RGset, qc, badSampleCutoff, detPcut, cutsamples, verbose, rm_badsamples) 
               
        
    #### probe wise QC filtering ######
    
    def probeFiltering(self, RGset=None, cutbead=3, zeropoint=True, verbose=True):
        
        if RGset:
            RGset=RGset
        else:
            RGset=self.RGset

        self.RGset_filt, self.pheno = robjects.r("""function (RGset, cutbead, zeropoint, verbose) { 


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
                RedSD <- RGset$RedSD
                GreenSD <- RGset$GreenSD

                idBeadmat <- beadmat < cutbead
                ##beadmat[idBeadmat] <- NA

                if(verbose)
                    cat("On average", round(100*sum(idBeadmat)/prod(dim(idBeadmat)), 2),"% of the probes (",nrow(idBeadmat),") have number of beads below", cutbead, "\n")

                ##Filter on Red and Green intensity <1
                if(zeroint) {
                    if(verbose)
                        cat("Filtering on zero intensities... \n")

                    Grn <- getGreen(RGset)
                    Red <- getRed(RGset)

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
                RedSD <- RedSD[idBeadmat]
                GreenSD <- GreenSD[idBeadmat]

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
    
    
        return self.RGset_filt
    
    
    def screeplot(self, RGset=None, nmax=10):
        if RGset:
            RGset=RGset
        else:
            RGset=self.RGset
        
        pcs,nmax = robjects.r("""function (RGset, nmax) {  
            library(matrixStats)
            screeplot <- function(RGset=RGset, nmax=nmax){
                
                extr<-minfi:::.extractFromRGSet450k(RGset)
                
                controlMatrix <- minfi:::.buildControlMatrix450k(extr)
               
                pc <- prcomp(controlMatrix)

                nmax <- ifelse(nmax > nrow(controlMatrix), nrow(controlMatrix), nmax)

                #barplot(summary(pc)$importance[2,1:nmax], ylab="Proportion of Variance", main="Scree Plot", col="#7cb4c9")
                #print(summary(pc)$importance[2,1:nmax])
                result=list(pc=summary(pc)$importance[2,1:nmax], nmax=nmax)
                return(result)
            }            
            
             pcs<-screeplot(RGset=RGset, nmax=nmax)
             
            
             result=list(pcs$pc,pcs$nmax)
             return(result)
            }""")(RGset, nmax)
      
        m=numpy2ri.ri2py(pcs)
        df=pd.DataFrame(np.array(m)).transpose()
        
        import seaborn as sns
        import matplotlib.pyplot as plt
        fig,ax=plt.subplots()
        #for i in df.transpose().to_numpy().tolist():
            #for j in i:
                #print(j)
                #sns.barplot(j)
        x = np.array(df.columns.tolist())
        y1 = df.to_numpy()[0]
        sns.barplot(x=x, y=y1, ax=ax)
        ax.set_title("Screeplot")
        ax.set_ylabel("Proportion of Variance")
        ax.set_xlabel("Principle components")
        plt.show()            
        
    
    def plt_covariates(self, pheno=None, matrix=None, pcs=2):
        
        if not matrix or not pheno:
            print('Please input data')
            return 
        
        cxy = robjects.r("""function (pcs, matrix, pheno) {               
             pheno <- pheno[rownames(pheno) %in% colnames(matrix) ,]
             
             df <- apply(pheno, 2, function(x) as.numeric(factor(x)))
             
             keep <- apply(df, 2, sd) > 0
             
             df <- df[ , keep]
             
             #library(irlba)
             pc <- prcomp(t(matrix), rank=pcs)
             #pc <- prcomp_irlba(t(matrix), n=pcs)
             
             cxy <- round(cor(pc$x, scale(df)),2)
             
             return(cxy)
            }""")(pcs, matrix, pheno)         
        
        
        corr=self.ri2py_dataframe(cxy, matrix=False)
        
        f, ax = plt.subplots(figsize=(10, 7))

        # Generate a custom diverging colormap
        cmap = sns.diverging_palette(220, 10, as_cmap=True)

        # Draw the heatmap with the mask and correct aspect ratio
        ax=sns.heatmap(corr.transpose(),  cmap=cmap, vmax=1, vmin=-1, center=0, annot=True,
                     linewidths=.5, cbar_kws={"shrink": .5, 'label': "Pearson's r"})
        ax.set_title("Heatmap of correlations")
        bottom, top = ax.get_ylim()
        ax.set_ylim(bottom + 0.5, top - 0.5) 
        plt.show() 
    
    
    
    def preprocessFunnorm(self, celltype_adoption=False, use_cell_count2=False, nPCs=2):   
            if celltype_adoption:
                if use_cell_count2:
                    self.NeunC,self.NeunC_table,_  = self.est_cell_counts2(processMethod="preprocessFunnorm", nPCs=nPCs)

                else:    
                    self.NeunC,self.NeunC_table, _ = self.est_cell_counts(processMethod="preprocessFunnorm", nPCs=nPCs)

                 

            self.GRset = robjects.r("""function (RGset, nPCs) {
                 mSetSq <- preprocessFunnorm(RGset, nPCs=nPCs, sex = NULL, bgCorr = TRUE, dyeCorr = TRUE, verbose = TRUE)
                return(mSetSq)
                }""")(self.RGset, nPCs)
            self.insert_cell_types() 
            return self.GRset


    def preprocessQuantile(self,celltype_adoption=False, use_cell_count2=False):   
            if celltype_adoption:
                if use_cell_count2:
                    self.NeunC,self.NeunC_table, _ = self.est_cell_counts2(processMethod="preprocessQuantile", nPCs=2)

                else:    
                    self.NeunC,self.NeunC_table, _ = self.est_cell_counts(processMethod="preprocessQuantile", nPCs=2)

                   

            self.GRset = robjects.r("""function (RGset) {
                    mSetSq <- preprocessQuantile(RGset, fixOutliers = TRUE, removeBadSamples = FALSE,
                               badSampleCutoff = 10.5, quantileNormalize = TRUE,
                               stratified = TRUE, mergeManifest = FALSE, sex = NULL,
                               verbose = TRUE)
                    return(mSetSq)
                    }""")(self.RGset)
            self.insert_cell_types()
            return self.GRset    


    def est_cell_counts(self,**kwargs):
            
            processMethod=list(kwargs.keys())[0]     
            processMethodname=list(kwargs.values())[0]
            PCs=list(kwargs.keys())[1]     
            nPCs=list(kwargs.values())[1]   
            rgset=self.RGset
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

    def est_cell_counts2(self, **kwargs):
            """Given RGSet object, estimate cell counts using reference approach via FlowSorted.Blood.EPIC 
            estimatecellcounts2 method.

            Parameters
            ----------
            rgset
                RGSet object stored in python via rpy2

           """

            rgset=self.RGset

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

        
        
    
    def rm_bad_probes_DNAmArray(self, GRset=None, RGset=None, remove_sex=True, removeXreactiveProbes=True, removeSNPs=True, cutbead=3, zeropoint=True, verbose=True, what="M", cutp=0.01, cutsamples=0.90, cutcpgs=0.90, normalization='funnorm', celltype_adoption=False, use_cell_count2=False):
        
        if GRset:
            GRset=GRset
        else:
            GRset = self.GRset
        if RGset:
            RGset=RGset
        else:
            RGset = self.RGset    
        
        self.RGset_filt=self.probeFiltering(cutbead=cutbead, zeropoint=zeropoint, verbose=zeropoint)  
        
        if removeSNPs:
            self.dropLociWithSnps(GRset=None, obj=None)
            
        if removeXreactiveProbes:
            self.excludeXreactiveprobes(GRset=None, obj=None, array_type='EPIC')
        
        if remove_sex:
            self.dropSexRelatedProbes( GRset=None, obj=None)
            
        self.mval_py=self.reduce(GRset=None, RGset=None, what=what, cutp=cutp, cutsamples=cutsamples, cutcpgs=cutcpgs, verbose=verbose)
        
        self.pheno=pheno[pheno['ID'].isin(self.mval_py.columns.to_list())]
        
        return self.mval_py
        
        
    def dropLociWithSnps(self, GRset=None, obj=None):
        if GRset:
            GRset=GRset
        else:
            GRset=self.GRset
        GRset = robjects.r("""function (grset) {             
             if(class(grset) != "GenomicRatioSet")
                    stop("The object should be of class 'GenomicRatioSet'!")
             grsetFlt <- dropLociWithSnps(grset)             
             return(grsetFlt)
            }""")(GRset)
        return GRset
    
    def filterXY(self, obj=None):
        # if your data includes males and females, remove probes on the sex chromosomes 
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
            bj <- obj[RemainProbe,]               
            return(obj)
            }""")(obj, verbose)
        return obj
    
    def excludeXreactiveprobes(self, obj=None, array_type='EPIC'):  
        import methylcheck
                    
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
    
    
    def filterCpGs(self, obj=None, dropSnPs=True, GRset=True, filterXY=True, filterNoCG=True, excludeXreactiveprobes=True, array_type='EPIC', verbose=True):
        orig_obj=obj        
        
        rgset = robjects.r("""function (obj) {             
             if(class(obj) == "RGChannelSet" | class(obj) == "RGChannelSetExtended"  ){                    
                    return(TRUE)  
                    }
                    else{
                    return(FALSE)
                    }             
            }""")(obj) 

        if pandas2ri.ri2py(rgset):
            print('Sorry this does not work with RGChannelSets, please use "GenomicRatioSet", "MethySets" or "RatioSets"! ')
            return obj                 
        
        
        if excludeXreactiveprobes:
            if verbose:
                print('Dropping cross- reactive probes')
            obj=self.excludeXreactiveprobes(obj,array_type=array_type)      
        if filterXY:
            if verbose:
                print('Dropping XY-Chromosome-related probes')
            obj=self.filterXY(obj)    
        if filterNoCG:            
            obj=self.filterNoCG(obj=obj, verbose=verbose) 
            
        grset = robjects.r("""function (obj) {             
            if(class(obj) == "GenomicRatioSet")
                    return(TRUE)  
             else return(FALSE)                
            }""")(obj)   
        
        if dropSnPs and not pandas2ri.ri2py(grset):
            if verbose:
                print('sorry can not remove snps object need to be of type "GenomicRatioSet"...\n processing without performing snp removal')             
        if dropSnPs and pandas2ri.ri2py(grset):
            obj=self.dropLociWithSnps(GRset=obj)
            self.GRset=obj                               
                      
        
        mset = robjects.r("""function (obj) {             
             if(class(obj) == "MethylSet")
                    return(TRUE)    
             else return(FALSE)                
            }""")(obj)
        if pandas2ri.ri2py(mset):   
                  self.Mset=obj
        
        rset = robjects.r("""function (obj) {             
             if(class(obj) == "RatioSet")
                    return(TRUE)                         
             else return(FALSE)
            }""")(obj)
        if pandas2ri.ri2py(rset):   
                  self.Rset=obj           
                  
        
        gmset = robjects.r("""function (obj) {             
             if(class(obj) == "GenomicMethylSet")
                    return(TRUE)                         
             else return(FALSE)
            }""")(obj)      
        if pandas2ri.ri2py(gmset):   
                  self.GMset=obj           
        
            
        if verbose:    
            robjects.r("""function (obj, orig_obj) { 
                start = length(rownames(orig_obj))
                left = length(rownames(obj))

                cat("\n In total there were ",start," probes for the analysis before filtering.")

                cat("\n",start-left," probes have been removed from further analysis.")

                cat("\n In total there are",left," probes left for the analysis.")    

                }""")(obj, orig_obj)            
        
        return obj
       
   
    
    def reduce(self, GRset=None, RGset=None, what="both", cutp=0.01, cutsamples=0.95, cutcpgs=0.95, verbose=True, autoimpute=True):
        
        if GRset:
            GRset=GRset
        else:
            GRset = self.GRset
        if RGset:
            RGset = RGset
        else:
            RGset = self.RGset_filt 

        obj, self.pheno = robjects.r("""function (GRset, RGset, what=c("beta", "M"), cutp, cutsamples, cutcpgs, verbose, autoimpute, ...) { 

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
            reduce <- function(GRset, RGset, what=c("beta", "M", "both"), cutp=0.01, cutsamples=0.95, cutcpgs=0.95, verbose=TRUE, autoimpute, ...) {

                what <- match.arg(what)

                if(class(GRset) != "GenomicRatioSet")
                    stop("First argument should be an object of class 'GenomicRatioSet' from preprocessFunnorm!")
                if(!(class(RGset) == "RGChannelSet" | class(RGset) == "RGChannelSetExtended"))
                    stop("Second argument should be an object of class RGChannelSet or RGChannelSetExtended from probeFiltering!")

                ##Filter on detection P-value using minfi's detectionP
                if(verbose)
                    cat("Calculate and filter on detection P-value... \n")

                pvalmat <- detectionP(RGset)
                idPvalmat <- pvalmat > cutp
                idPvalmat[is.na(idPvalmat)] <- TRUE ##set those for which a detection P-value could not be calculate TRUE to be filtered out
                ##pvalmat[idPvalmat] <- NA

                if(verbose)
                    cat("On average", round(100*sum(idPvalmat)/prod(dim(idPvalmat)), 2),"% of the CpGs (",nrow(idPvalmat),") have detection P-value above the threshold ",cutp, "\n")

                if(verbose)
                    cat("Transform to ",what,"-values... \n")

                if(what == "M"){
                    matfilt <- getM(preprocessRaw(RGset), ...)
                    matnorm <- getM(GRset, ...)
                }
                if (what =="beta"){
                    matfilt <- getBeta(RGset, ...)
                    matnorm <- getBeta(GRset, ...)
                }
                else 
                {  if(what=="both") 
                    {
                        matfiltbeta <- getBeta(RGset, ...)
                        matfiltM <- getM(preprocessRaw(RGset), ...)
                        matnormbeta <- getBeta(GRset, ...)
                        matnormM <- getM(GRset, ...)
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

                    ##Replaced by gap_hunting
                    ##set chen CpGs/probes NA
                    ##if(verbose)
                    ##    message("Removing cross-reactive or polymorphic probes...")

                    ##data("chen", package="Leiden450K")
                    ##matnorm[rownames(matnorm) %in% names(chenProbes),] <- NA

                    ##calculate success rates and reduce
                    if(verbose)
                        cat("Calculate success rates and reduce... \n")

                    ##srCols <- apply(matnorm, 2, function(x) sum(!is.na(x))/(length(x) - 30969)) ##chen CpGs excluded
                    srCols <- apply(matnorm, 2, function(x) sum(!is.na(x))/(length(x)))
                    srRows <- apply(matnorm, 1, function(x) sum(!is.na(x))/length(x))

                    if(verbose){
                        cat("Percentage of samples having success rate above", cutsamples, "is", round(100*sum(srCols > cutsamples)/length(srCols),2),"% \n")
                        cat("Percentage of CpGs having success rate above", cutcpgs, "is", round(100*sum(srRows > cutcpgs)/length(srRows),2),"% \n")
                    }
                    matnorms<-matnorm[srRows > cutcpgs,  srCols > cutsamples]
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
                        cat("Calculate success rates and reduce... \n")

                    srColsbeta <- apply(matnormbeta, 2, function(x) sum(!is.na(x))/(length(x)))
                    srRowsbeta <- apply(matnormbeta, 1, function(x) sum(!is.na(x))/length(x))
                    
                    srColsM <- apply(matnormM, 2, function(x) sum(!is.na(x))/(length(x)))
                    srRowsM <- apply(matnormM, 1, function(x) sum(!is.na(x))/length(x))

                    if(verbose){
                        cat("Percentage of beta samples having success rate above", cutsamples, "is", round(100*sum(srColsbeta > cutsamples)/length(srColsbeta),2),"% \n")
                        cat("Percentage of beta CpGs having success rate above", cutcpgs, "is", round(100*sum(srRowsbeta > cutcpgs)/length(srRowsbeta),2),"% \n")
                        cat("Percentage of M samples having success rate above", cutsamples, "is", round(100*sum(srColsM > cutsamples)/length(srColsM),2),"% \n")
                        cat("Percentage of M CpGs having success rate above", cutcpgs, "is", round(100*sum(srRowsM > cutcpgs)/length(srRowsM),2),"% \n")
                    }
                    matnormsbeta<-matnormbeta[srRowsbeta > cutcpgs,  srColsbeta > cutsamples]
                    matnormsM<-matnormM[srRowsM > cutcpgs,  srColsM > cutsamples]
                }
                
                
                
                if(autoimpute)
                {  
                   if(verbose){
                   cat("\n  Autoimpute Start") 
                   }
                   if(sum(is.na(idPvalmat)) == 0) {
                      
                       if(verbose){
                       cat("    No NAs (failed probes) exist in your data set any more, you don't need to do any imputation.")
                       }
                   } 
                   
                   else
                   {   if(verbose){
                       cat("    There are ",sum(is.na(idPvalmat)), " NAs (failed probes) exists in your data set.")
                       cat("    Impute.knn will be conducted for remain NAs. (NOT suitable for small data sets)\n")
                       }
                      
                      ### Using sink function to remove messages
                        zz <- file("ImputeMessage.Rout", open="wt")
                        sink(zz)
                        sink(zz, type="message")               

                        if(verbose){
                            
                            cat("    Doing imputation on matrix.")
                            }
                         library(impute)
                         if (what!="both")
                         {
                             matnorms <- impute.knn(matnorms,k=5)$data 
                         }
                         else
                         {
                             matnormsM <- impute.knn(matnormsM,k=5)$data
                             matnormsbeta <- impute.knn(matnormsbeta,k=5)$data
                         }
                         sink(type="message")
                         sink()                       

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


            object<-reduce(GRset, RGset, what, cutp, cutsamples, cutcpgs, verbose, autoimpute)
            pheno = pData(RGset)
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


     }""")(GRset, RGset, what, cutp, cutsamples, cutcpgs, verbose, autoimpute)
        
        
        if (what=="both"):
            self.mval=obj[0]
            self.mval_py = self.ri2py_dataframe(r_dat=self.mval, matrix=False)
            self.beta=obj[1]
            self.beta_py = self.ri2py_dataframe(r_dat=self.beta, matrix=False)
            return self.beta_py, self.mval_py
        
        elif(what=="M"):   
            self.mval=obj
            self.mval_py = self.ri2py_dataframe(r_dat=self.mval, matrix=False)
            return self.mval_py
        
        elif(what=="beta"):   
            self.beta=obj
            self.beta_py = self.ri2py_dataframe(r_dat=self.beta, matrix=False)  
            return self.mval_py      
    
        
    
    
    def dmp_finder(self,matrix, pheno, phenotype=None, adjust_vars=None, correction_vars=None, useCombat=False, sva=False, number=10000, pvalue=0.05, adjpval=1, save_csv=False, path=None ):  
            from scipy.special import comb
            if not phenotype:
                print('Please specify a target of interest \n'
                    'These are the available column names:')
                print(pandas2ri.ri2py(robjects.r['as'](pheno,'data.frame')).columns.tolist())
                return None
            try:
                pheno_py=pd.DataFrame(pandas2ri.ri2py(pheno))
            except:
                pheno_py=pandas2ri.ri2py(robjects.r['as'](pheno,'data.frame'))            

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

            if sva:
                print('You have chosen to include surrogate variable analysis')  
            print('You are adjusting for these variables: %s' % (adjust_vars))
            print('You are correcting for these variables: %s' % (correction_vars))


            self.get_annotation()     
            #mod, mod0

            self.dmps, data, dectest= robjects.r("""function ( pheno, phenotype, M,  ann, number, array, coef, pvalue, adjpval, adj_vars,corr_vars, sva,useCombat) 
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

            ####SVA#####

            if (sva)
            {    cat('creating model for SVA')
                 n.sv = num.sv(M,mod,method="leek")

                 svobj = sva(M,mod,mod0,n.sv=n.sv)

                 modSv = cbind(mod,svobj$sv)
                 colnames(modSv)<- make.names(colnames(modSv))
                 rownames(modSv)<-make.names(rownames(modSv))
            }

            else
            {
                 cat('creating model')
                 modSv = mod

            }  

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



            # look at the numbers of DM CpGs at FDR < 0.05
            cat('Computing statistics for experiment')    
            dectest<- summary(decideTests(fit2, p.value=pvalue))           

            cat('Aligning annotation')
            #annSub <- ann[match(rownames(M),ann$Name),]
            annSub <- ann[ann$Name %in% rownames(M),]
            #annSub <- ann[match(rownames(fit2$coefficients),ann$Name),]

            if (coef == 'all'){

                            cat('Computing contrasts for experiment')           
                            topresults <- c()             
                            for (i in 1:choose(length(array),2)) {        
                                #print(i)

                                top<-topTable(fit2, genelist=annSub, coef=i,number=number, p.value=adjpval )

                                topresults[[i]] <- top    

                            }

                            result=list(topresults,matr$des, dectest)
                        }

                        else{

                            cat('Computing contrast for experiment') 
                            top<-topTable(fit2, genelist=annSub, coef=coef, number=number, p.value=adjpval )   
                            top<-list(top)
                            result=list(top,matr$des,dectest)

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
                   useCombat)     

            self.dmps_list=[]
            for elem in self.dmps:

                self.dmps_list.append(pd.DataFrame(pandas2ri.ri2py(elem)))
            lists=[]
            for elements in self.dmps_list:            
                for elem in elements['Name'].tolist():         
                    lists.append(elem)  
            self.dmps = list(dict.fromkeys(lists))          
            self.contrastmatrix=self.ri2py_dataframe(data)

            self.dectest=self.ri2py_dataframe(dectest)
            ###to-do
            if save_csv:
                if path:
                    for elem in self.dmps_list:


                        elem.to_csv(path+'{}.csv'.format(i,contrasts[i]), index=False)
                else:
                    for i,elem in enumerate(self.dmps_list):


                        elem.to_csv(self.idat_dir+'{}.csv'.format(i,contrasts[i]), index=False)

            print('done') 
    
    
    
    
    
    


    
    

    def move_jpg(self):
        """Move jpeg files from current working directory to the idat directory.
        """
        subprocess.call('mv *.jpg {}'.format(self.idat_dir),shell=True)
    
    
    def preprocess_cns_data(self):
        self.qcinfo = self.enmix.QCinfo(self.RGset, detPthre=1e-7)
        self.move_jpg()      
        print('raw processing')
        self.MSet = self.minfi.preprocessRaw(self.RGset)
        print('filtering MSet')
        self.MSet = self.enmix.QCfilter(self.MSet,qcinfo=self.qcinfo,outlier=True)
        self.move_jpg()                 
        return self.MSet
    
    def preprocessNoob(self):
        """Run minfi preprocessing with Noob normalization"""
        self.qcinfo = self.enmix.QCinfo(self.RGset, detPthre=1e-7)
        self.move_jpg()
        self.MSet = self.minfi.preprocessNoob(self.RGset)
        self.move_jpg()
        self.MSet = self.enmix.QCfilter(self.MSet,qcinfo=self.qcinfo,outlier=True)
        self.move_jpg()
        return self.MSet

    def preprocessRAW(self):
        """Run minfi preprocessing with RAW normalization"""
        self.qcinfo = self.enmix.QCinfo(self.RGset, detPthre=1e-7)
        self.move_jpg()
        self.MSet = self.minfi.preprocessRaw(self.RGset)
        self.move_jpg()
        self.MSet = self.enmix.QCfilter(self.MSet,qcinfo=self.qcinfo,outlier=True)
        self.move_jpg()
        return self.MSet

    def preprocessENmix(self, n_cores=6):
        """Run ENmix preprocessing pipeline.
        Parameters
        ----------
        n_cores
            Number of CPUs to use."""
        self.qcinfo = self.enmix.QCinfo(self.RGset, detPthre=1e-7)
        self.move_jpg()
        self.MSet = self.enmix.preprocessENmix(self.RGset, QCinfo=self.qcinfo, nCores=n_cores)
        self.move_jpg()
        self.MSet = self.enmix.QCfilter(self.MSet,qcinfo=self.qcinfo,outlier=True)
        self.move_jpg()
        return self.MSet

    def preprocess_quant(self):
        self.qcinfo = self.enmix.QCinfo(self.RGset, detPthre=1e-7)
        self.move_jpg()
        rgset=self.RGset
        self.GRset = robjects.r("""function (RGset) {
            cellCounts <- preprocessQuantile(RGset, fixOutliers = TRUE, removeBadSamples = TRUE,
                       badSampleCutoff = 10.5, quantileNormalize = TRUE,
                       stratified = TRUE, mergeManifest = FALSE, sex = NULL,
                       verbose = TRUE)
            return(cellCounts)
            }""")(rgset)
        return self.GRset
    
    
    
    def preprocessNOOB(self):
        self.MSet = self.minfi.preprocessNoob(self.RGSet)
        return self.MSet
    
    def preprocessRaw(self):
        self.MSet = robjects.r("""function (RGset) {
           mSetRaw <- preprocessRaw(RGset) 
           return(mSetRaw)
            }""")(self.RGset)
        return self.MSet
    
        
           
    
        
    
    
    
    
    def dmp_find(self, M_val, GRset, path, phenotype=None, save_csv=False, number=10000):  
        from scipy.special import comb
        if not phenotype:
            print('Please specify a target of interest \n'
                'These are the available column names:')
            print(pandas2ri.ri2py(robjects.r['as'](self.pheno,'data.frame')).columns.tolist())
            return
        try:
            pheno_py=pd.DataFrame(pandas2ri.ri2py(self.pheno))
        except:
            pheno_py=pandas2ri.ri2py(robjects.r['as'](self.pheno,'data.frame'))            
        
            
        py_array=pheno_py[phenotype].unique()
        #prtink=0
        contrasts=dict()
        n = len(py_array)
        k=0
        for i in range(0,n):
            #print(k)
            for j in range(i+1,n):
                k+=1
                #print(k)
                contrasts[k]= py_array[i]+"-"+py_array[j]
                #print(py_array[i]+"-"+py_array[j])
                
                
        
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
            #return 0,0       
              
          
        r_array=numpy2ri.py2ri(py_array)
        #self.compute_mvals(self.GRset)
        self.get_annotation()
        
        
        self.dmps, data,cols,rows, self.dectest = robjects.r("""function ( pheno, phenotype, M, GRset,ann, number, array, coef) {
    
            library(limma)                  

            rownames(M) <- featureNames(GRset)           

            design <- model.matrix(~0+phenotype,data=pheno)
            colnames(design) <- array
            colnames(design)<- make.names(colnames(design))
            samplenames <- sampleNames(GRset)

            fit <- lmFit(M,design)
            
            # create a contrast matrix for specific comparisons
            
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
                        colnames(design)[k] <- paste(levels[i],"-",levels[j],sep="")
                        }
                    result<-list(des=design,cols=colnames(design), rows=rownames(design) )    
                    return(result)    
                    #design
                    }         
           
            matr <-design.pairs(array)                                          
                                        
            #print(contMatrix)
            #print(matr$des)
            fit2 <- contrasts.fit(fit, matr$des)
            fit2 <- eBayes(fit2)

            # look at the numbers of DM CpGs at FDR < 0.05
            dectest<- summary(decideTests(fit2))
            
           
            annSub <- ann[match(rownames(M),ann$Name),]
                             # c(1:4,12:19,24:ncol(ann))]
                             
            if (coef == 'all'){
                
                           
                topresults <- c()             
                for (i in 1:choose(length(array),2)) {                                  

                    top<-topTable(fit2, coef=i, genelist=annSub, number=number )  
                    topresults[[i]] <- top                    
                }
                result=list(topresults,matr$des,matr$cols,matr$rows, dectest)
            }
            
            else{
                
                top<-topTable(fit2, coef=coef, genelist=annSub, number=number )
                result=list(top,matr$des,matr$cols,matr$rows, dectest)
            
            }       

            return(result)
         }""")(self.pheno, 
               pheno_py[phenotype], 
               M_val, 
               GRset, 
               self.annotation, 
               number, 
               r_array, 
               coef)
        
        self.dmps_list=[]
        for elem in self.dmps:
            #print(elem)
            self.dmps_list.append(pd.DataFrame(pandas2ri.ri2py(elem)))
        lists=[]
        for elements in self.dmps_list:            
            for elem in elements['Name'].tolist():         
                lists.append(elem)  
        self.dmps = list(dict.fromkeys(lists))          
        data_py=pd.DataFrame(pandas2ri.ri2py(data))
        rows_py=pd.Series(pandas2ri.ri2py(rows))
        cols_py=pd.Series(pandas2ri.ri2py(cols))
        self.contrastmatrix=pd.DataFrame(data_py.to_numpy(),index=rows_py,columns=cols_py)
        ###to-do
        if save_csv:
            if path:
                for elem in self.dmps_list:
                    
                    #dmps_list.append(pd.DataFrame(pandas2ri.ri2py(elem)))
                    #elem.to_csv(, index=False)
                    elem.to_csv(path+'{}.csv'.format(i,contrasts[i]), index=False)
            else:
                for i,elem in enumerate(self.dmps_list):
                    
                    #dmps_list.append(pd.DataFrame(pandas2ri.ri2py(elem)))
                    elem.to_csv(self.idat_dir+'{}.csv'.format(i,contrasts[i]), index=False)
                    #dmps_py.to_csv(self.idat_dir+'dmps.csv',index=False)
                    #print("{:3d} {:4d} {:5d}".format(i, i*i, i*i*i))
                
        
    
    def dmr_find_bumphunter(self, path, n_cores=4, cutoff=0.2, iterations=0, phenotype='disease', save_csv=False):
                
        if phenotype not in pandas2ri.ri2py(robjects.r['as'](self.pheno,'data.frame')).columns.tolist():
            print('The pheno sheet does not contain a '+phenotype+' column you specified \n'
                     'These are the available column names:')
            print(pandas2ri.ri2py(robjects.r['as'](self.pheno,'data.frame')).columns.tolist())
            return
        
        setter=self.GRset            
        
                    
        self.bump_dmrs = robjects.r("""function (Setter, pheno, Phenotype, N_cores,Cutoff, Iterations) {

                designMatrix <- model.matrix(~0+Phenotype, data=pheno)
                colnames(designMatrix) <- c('FCD IIa', 'FCD IIb', 'FCD 1a', 'TSC', 'PMG 1q', 'PMG')
                library(doParallel)
                registerDoParallel(cores = N_cores)
                dmrs <- bumphunter(Setter, design = designMatrix, 
                    cutoff = Cutoff, B=Iterations, type="Beta")
                return(dmrs)
                }""")(setter, self.pheno, pandas2ri.ri2py(self.pheno)[phenotype], n_cores, cutoff, iterations)
        
        bump_dmrs_py = pd.DataFrame(pandas2ri.ri2py(self.bump_dmrs[0]))
        if save_csv:
            if path:
                bump_dmrs_py.to_csv(path, index=False)
            else:
                bump_dmrs_py.to_csv(self.idat_dir+'bump_dmrs.csv',index=False)
        return bump_dmrs_py

    def get_annotation(self):
        self.annotation=robjects.r("""function (RGset) {
            ann450k <- getAnnotation(RGset)
            return(ann450k)
            
            }""")(self.RGset)
        self.annotation_py=pandas2ri.ri2py(robjects.r['as'](self.annotation,'data.frame'))
        return self.annotation_py
    
       
            
    def ensure_ordering(self, cutoff=0.05):
        # ensure probes are in the same order in the mSetSq and detP objects
        # remove any probes that have failed in one or more samples
        self.detectionP(RGset=self.RGset)
        self.GRset, self.detP, self.RGset = robjects.r("""function (grset, detP, cutoff, pheno, rgSet) {        
             detP <- detP[match(featureNames(grset),rownames(detP)),]
             keep <- rowSums(detP < cutoff) == ncol(grset) 
             rgSet <- rgSet[keep,]
             grsetFlt <- grset[keep,]              
             result = list(grsetFlt, detP, rgSet)
             return(result)
            }""")(self.GRset, self.detP, cutoff, self.pheno, self.RGset)
        return self.GRset
    
    def convert_sets(self,mset=True, rset=True):
        if mset:
            self.RSet = self.minfi.ratioConvert(self.MSet, what = "both", keepCN = True)
        if rset:
            self.GRset = self.minfi.mapToGenome(self.RSet)
            return self.GRset
        return self.RSet
    
    
    
    
        
    

    def plot_qc_metrics(self, output_dir):
        """Plot QC results from ENmix pipeline and possible minfi. Still experimental.
        Parameters
        ----------
        output_dir
            Where to store plots."""
        self.enmix.plotCtrl(self.RGset)
        grdevice = importr("grDevices")
        geneplotter = importr("geneplotter")
        base = importr('base')
        anno=self.minfi.getAnnotation(self.RGset)
        anno_py = pandas2ri.ri2py(robjects.r['as'](anno,'data.frame'))
        beta_py = pandas2ri.ri2py(self.beta)
        beta1=numpy2ri.py2ri(beta_py[anno_py["Type"]=="I"])
        beta2=numpy2ri.py2ri(beta_py[anno_py["Type"]=="II"])
        grdevice.jpeg(output_dir+'/dist.jpg',height=900,width=600)
        base.par(mfrow=robjects.vectors.IntVector([3,2]))
        self.enmix.multidensity(self.beta, main="Multidensity")
        self.enmix.multifreqpoly(self.beta, xlab="Beta value")
        self.enmix.multidensity(beta1, main="Multidensity: Infinium I")
        self.enmix.multifreqpoly(beta1, main="Multidensity: Infinium I", xlab="Beta value")
        self.enmix.multidensity(beta2, main="Multidensity: Infinium II")
        self.enmix.multifreqpoly(beta2, main="Multidensity: Infinium II", xlab="Beta value")
        grdevice.dev_off()
        self.minfi.qcReport(self.RGset, pdf = "{}/qcReport.pdf".format(output_dir))
        self.minfi.mdsPlot(self.RGset)
        self.minfi.densityPlot(self.RGset, main='Beta', xlab='Beta')

    def return_beta(self):
        """Return minfi RSet after having created MSet."""
        self.RSet = self.minfi.ratioConvert(self.MSet, what = "both", keepCN = True)
        return self.RSet

    def get_beta(self, bmiq = False, n_cores=6):
        """Get beta value matrix from minfi after finding RSet."""
        from pymethylprocess.meffil_functions import bmiq_mc
        self.beta = self.minfi.getBeta(self.RSet)
        if bmiq:
            self.beta = bmiq_mc(self.beta, nCores=n_cores, nfit=10000)
        return self.beta

    def filter_beta(self):
        """After creating beta, filter out outliers."""
        self.beta_final=self.enmix.rm_outlier(self.beta,qcscore=self.qcinfo)
        return self.beta_final

    def get_meth(self):
        """Get methylation intensity matrix from MSet"""
        return self.minfi.getMeth(self.MSet)

    def get_unmeth(self):
        """Get unmethylated intensity matrix from MSet"""
        return self.minfi.getUnmeth(self.MSet)

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

    def preprocess_enmix_pipeline(self, n_cores=6, pipeline='enmix', noob=False, qc_only=False, use_cache=False, bmiq=False):
        """Run complete ENmix or minfi preprocessing pipeline.
        Parameters
        ----------
        n_cores
            Number CPUs.
        pipeline
            Run enmix or minfi
        noob
            Noob norm or RAW if minfi running.
        qc_only
            Save and quit after only running QC?
        use_cache
            Load preexisting RGSet instead of running QC again."""
        cache_storage_path = os.path.join(self.idat_dir,'RGSet.rds')
        self.load_idats(use_cache=use_cache)
        
        if qc_only:
            robjects.r('saveRDS')(self.RGset,cache_storage_path)
            exit()
        if not os.path.exists(cache_storage_path):
            robjects.r('saveRDS')(self.RGset,cache_storage_path)
        if pipeline =='enmix':
            self.preprocessENmix(n_cores)
        else:
            if noob:
                self.preprocessNoob()
            else:
                self.preprocessRAW()
        self.return_beta()
        if bmiq:
            self.get_beta(bmiq = True, n_cores=n_cores)
        else:
            self.get_beta()
            
          
        self.filter_beta()       
        
        
        self.extract_pheno_data(methylset=True, grset=False)
        return self.pheno, self.beta_final
    
    def preprocess_cns_pipeline(self, n_cores=6, pipeline='enmix',  noob=False, qc_only=False, use_cache=False, bmiq=False, array_type='epic', rm_sex=False, use_cell_count2=False):
        """Run complete ENmix or minfi preprocessing pipeline.
        Parameters
        ----------
        n_cores
            Number CPUs.
        pipeline
            Run enmix or minfi
        noob
            Noob norm or RAW if minfi running.
        qc_only
            Save and quit after only running QC?
        use_cache
            Load preexisting RGSet instead of running QC again."""
        print('starting')
        cache_storage_path = os.path.join(self.idat_dir,'RGSet.rds')
        if use_cache:
            self.RGset=robjects.r('readRDS')(cache_storage_path)
        else:
            self.load_idats()
        if qc_only:
            robjects.r('saveRDS')(self.RGset,cache_storage_path)
            exit()
        if not os.path.exists(cache_storage_path):
            robjects.r('saveRDS')(self.RGset,cache_storage_path)
        
                
        if pipeline =='enmix':
            self.preprocessENmix(n_cores)
        else:
            if noob:
                self.preprocessNoob()
            else:
                self.preprocessRAW()
        
        if use_cell_count2:
            self.NeunC,self.NeunC_table,self.GRset = self.est_cell_counts2(use_mset=False)
        else:    
            self.NeunC,self.NeunC_table,self.GRset = self.est_cell_counts()
        
        self.return_beta()
        if bmiq:
            self.get_beta(bmiq = True, n_cores=n_cores)
        else:
            self.get_beta()            
          
        self.filter_beta()   
                
        self.extract_pheno_data(methylset=True, grset=False)
        print('finished')
        return self.pheno, self.beta_final
    
    def preprocess_cns_pipeline2(self, n_cores=6, qc_only=False, use_cache=False, bmiq=False, array_type='epic', rm_sex=False, use_cell_count2=False):
        """Run complete ENmix or minfi preprocessing pipeline.
        Parameters
        ----------
        n_cores
            Number CPUs.
        pipeline
            Run enmix or minfi
        noob
            Noob norm or RAW if minfi running.
        qc_only
            Save and quit after only running QC?
        use_cache
            Load preexisting RGSet instead of running QC again."""
        print('starting')
        cache_storage_path = os.path.join(self.idat_dir,'RGSet.rds')
        if use_cache:
            self.RGset=robjects.r('readRDS')(cache_storage_path)
        else:
            self.load_idats()
        if qc_only:
            robjects.r('saveRDS')(self.RGset,cache_storage_path)
            exit()
        if not os.path.exists(cache_storage_path):
            robjects.r('saveRDS')(self.RGset,cache_storage_path)
        
                
        self.preprocess_quant()
        
        print('estimating cell counts')   
        if use_cell_count2:
            self.NeunC,self.NeunC_table,_ = self.est_cell_counts2(use_mset=False)
        else:    
            self.NeunC,self.NeunC_table,_ = self.est_cell_counts()
        
        self.beta = self.minfi.getBeta(self.GRset)
        if bmiq:
            from pymethylprocess.meffil_functions import bmiq_mc
            self.beta = bmiq_mc(self.beta, nCores=n_cores, nfit=10000)          
          
        self.filter_beta()   
        
        
        self.extract_pheno_data(methylset=False, grset=True)
        print('finished')
        return self.pheno, self.beta_final
    
    def est_cell_counts2_dlpfc(self, use_mset=True):
        """Given RGSet object, estimate cell counts using reference approach via FlowSorted.Blood.EPIC 
        estimatecellcounts2 method.
        
        Parameters
        ----------
        rgset
            RGSet object stored in python via rpy2
            
        or an
        
        mset    
            MSet object stored in python via rpy2
            """
        if use_mset:
            mset=self.MSet
        else:
            mset=self.GRset
        robjects.r('library(FlowSorted.Blood.EPIC)')
        cell_count_estimates = robjects.r("""function (Mset) {
            cellCounts <- estimateCellCounts2(Mset, compositeCellType = "DLPFC",
                       processMethod = "preprocessQuantile", probeSelect = "auto",
                       cellTypes = c("NeuN_neg", "NeuN_pos"),
                       referencePlatform = c("IlluminaHumanMethylation450k",
                                             "IlluminaHumanMethylationEPIC",
                                             "IlluminaHumanMethylation27k"),
                       referenceset = NULL, IDOLOptimizedCpGs = NULL, returnAll = TRUE,
              meanPlot = FALSE, verbose = TRUE, lessThanOne = FALSE, fixOutliers=TRUE, removeBadSamples = TRUE )
              
            return(cellCounts)
            }""")(mset)        
        return cell_count_estimates
        ##beta<-getBeta(cellCounts$normalizedData)
        
    
    
    
    def plot_original_qc(self, output_dir):
        """Plot QC results from ENmix pipeline and possible minfi. Still experimental.
        Parameters
        ----------
        output_dir
            Where to store plots."""
        self.preprocessRAW()
        self.return_beta()
        self.get_beta()
        self.plot_qc_metrics(output_dir)

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
    
    
    def findDMRS(self, n_cores=4, cutoff=0.2, iterations=0, methylset=True):
        """Given RGSet object, estimate cell counts using reference approach via minfi.
            Parameters
            ----------
            phen_group
               The phenotype to be tested for association with methylation.
            typ
                Is the phenotype '‘continuous’ or ‘categorical’?
                """
       
        setter=self.MSet if methylset else self.GRset
        phenotype=phen_group        
        n_cores=n_cores
        cutoff=cutoff
        iterations=iterations

            #self.pheno = robjects.r("pData")(self.MSet) #if methylset else robjects.r("pData")(self.RGset)
            #phenotype=self.pheno[phen_group]
        dmrs = robjects.r("""function (Setter, Phenotype, N_cores,Cutoff, Iterations) {
            grp  <- pData(Setter)$disease
            designMatrix <- model.matrix(~ grp)
            library(doParallel)
            registerDoParallel(cores = N_cores)
            dmrs <- bumphunter(Setter, design = designMatrix, 
                cutoff = Cutoff, B=Iterations, type="Beta")
            return(dmrs)
            }""")(setter, phenotype, n_cores, cutoff, iterations)
        
        #dmrs=pd.DataFrame(pandas2ri.ri2py(dmrs))
        #dmrs=pandas2ri.ri2py(robjects.r['as'](dmrs,'data.frame'))
        return dmrs  
    
     

        #pheno <- pData(GRset)$status
        #designMatrix <- model.matrix(~ pheno)

        #library(doParallel)
        #registerDoParallel(cores = 3)

        #dmrs <- bumphunter(self.GRset, design = designMatrix, 
        #         cutoff = 0.2, B=0, type="Beta")


        #dmrs <- bumphunter(self.GRset, design = designMatrix, 
          #       cutoff = 0.2, B=1000, type="Beta")

    def findDMP(self, phen_group='disease', typ='categorical', methylset=True ):
            """Given RGSet object, estimate cell counts using reference approach via minfi.
            Parameters
            ----------
            phen_group
               The phenotype to be tested for association with methylation.
            typ
                Is the phenotype '‘continuous’ or ‘categorical’?
                """
            
            beta=self.beta_final        
            setter=self.MSet if methylset else self.GRset 
            phenotype=phen_group
            typ=typ
            #self.pheno = robjects.r("pData")(self.MSet) #if methylset else robjects.r("pData")(self.RGset)
            #phenotype=self.pheno[phen_group]
            dmp = robjects.r("""function (Beta,Setter, Phenotype, Typ) {
                beta <- getBeta(Setter)
                grp  <- pData(Setter)$disease
                dmp <- dmpFinder(Beta, pheno=grp, type=Typ)
                return(dmp)
                }""")(beta, setter, phenotype, typ)
            #dmp=pd.DataFrame(pandas2ri.ri2py(dmp))
            #dmp=pandas2ri.ri2py(robjects.r['as'](dmp,'data.frame'))
            return dmp   
        
  