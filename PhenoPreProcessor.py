import pandas as pd
import numpy as np
import os, glob
from os.path import basename
from collections import Counter
'''Most if this is originally from "https://github.com/Christensen-Lab-Dartmouth/PyMethylProcess/blob/master/pymethylprocess/PreProcessDataTypes.py" '''

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
            Maps uuids to proper tcga sample names, should be downloaded with tcga clinical information.
            
            
            This function has been edited by our group!!!
            
            
            """
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
        idat_basenames = np.unique(np.vectorize(lambda x: '_'.join(x.replace(self.idat_dir, '').split('_')[:-1]))(idats))        
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
        
        '''Function to fromat the sex column name into "M" for male and "F" for female '''
        
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
        
    def combine_columns(self, first_column=None, second_column=None, column_name=None, seperator='.'):
        '''New function to fuse two column data into one with a separation indicated by the "seperator" '''
        if column_name is None:
            print('You did not specify a column name... using default value and creating "ID" column')
            column_name='ID'            
        
        if not isinstance(column_name, str) :
            print('Please specify the name of type string')
            return
        
        # give the samples descriptive names
        check_list=[first_column, second_column]
        for val in check_list: 
            if val not in self.pheno_sheet.columns.tolist():
                print('The pheno sheet does not contain a '+first_column+' or '+second_column+' column you specified \n'
                     'These are the available column names:')
                print(self.pheno_sheet.columns.tolist())
                return
            
        self.pheno_sheet[column_name]=self.pheno_sheet[first_column].apply(str)+seperator+self.pheno_sheet[second_column].apply(str)
        
        return self.pheno_sheet 
    
    
    def remove_columns(self, column_name=None):
        '''Function to remove a specific column from the pheno sheet indicated with the "column_name" '''
        
        if column_name is None:
            print('Please specify the column by name you want to remove \n'
                  'These are the available column names:\n')
            print(self.pheno_sheet.columns.tolist())
            return
        column_name_list=[]
        if isinstance(column_name, list):
            
            for name in column_name:
                print(name)
                if isinstance(name, str):
                    
                    column_name_list.append(name)
                if not isinstance(name, str):
                    try:
                        column_name_list.append(str(name))
                    except:
                        print('ommiting '+name+'from columns to remove because it could not be converted to type string')                           
        if isinstance(column_name, str):
            
            column_name_list=[column_name]   
        elif not isinstance(column_name, str) and not isinstance(column_name, list):
            
            try:
                column_name_list=[str(column_name)]
            except:
                print('ommiting '+column_name+'from columns to remove because it could not be converted to type string')  
                
        if column_name_list==[]:
            print('No valid columns to remove where given or found')
            return        
        
        self.pheno_sheet=self.pheno_sheet.drop(column_name_list, axis=1)            
       
        return self.pheno_sheet