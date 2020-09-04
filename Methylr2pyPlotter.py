import pandas as pd
import numpy as np
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri, numpy2ri
from rpy2.robjects.packages import importr 
from rpy2.robjects import r
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

class Plotter:
    """Class that will take an Idat_Processor and plot all associated QC and related values.
    """
    def __init__(self, IDAT_Processor):
        self.IDAT_Processor = IDAT_Processor
        
 
    def plt_mds(self, dataframe=None, pheno=None,  top=1000, n_components=2, group='disease', components=(0,1)):

        if dataframe is None:
            beta=self.IDAT_Processor.minfi.getBeta(self.IDAT_Processor.RGset)
            dataframe=self.IDAT_Processor.ri2py_dataframe(beta, matrix=False).transpose()

        if pheno is None:
            try:
                pheno=self.IDAT_Processor.ri2py_dataframe(self.IDAT_Processor.pheno, matrix=True)
            except:
                pheno=self.IDAT_Processor.ri2py_dataframe(self.IDAT_Processor.pheno, matrix=False)

        import matplotlib.pyplot as plt
        import matplotlib.cm as cm
        from sklearn.manifold import MDS

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

        X_transformed['ID']=pheno['ID'].to_numpy()
        X_transformed['case_ID']=pheno['case_ID'].to_numpy()

        X_transformed=X_transformed.rename(columns={components[0]: "a", components[1]: "b"})

        output_notebook()

        p = figure(
            tools="hover,pan,wheel_zoom,save",
            toolbar_location="above",
            title='MDS Plot',
            plot_width=900,
            plot_height=700,
        )

        if len(X_transformed[group].unique())!=2:
            palette = d3['Category20'][len(X_transformed[group].unique())]
            color_map = bmo.CategoricalColorMapper(factors=X_transformed[group].unique(),
                                               palette=palette)
        elif len(X_transformed[group].unique())==2:
            colors=['red', 'green']
            #from bokeh.palettes import brewer
            #colors = brewer["Spectral"][len(X_transformed[group].unique())]

            # Create a map between factor and color.
            colormap = {k: colors[i] for i,k in enumerate(X_transformed[group].unique())}

            # Create a list of colors for each value that we will be looking at.
            X_transformed['colors'] = [colormap[x] for x in X_transformed[group]]


        for key in X_transformed[group].unique():

            keys=ColumnDataSource(X_transformed[X_transformed[group]==key])

            if len(X_transformed[group].unique())!=2:
                p.scatter(x='a', y='b', size=10, source=keys ,legend_label=key, color={'field': group, 'transform': color_map})
            else:
                p.scatter(x='a', y='b', size=10, source=keys ,legend_label=key, color='colors')


        p.legend.location = "top_left"
        p.legend.click_policy="hide"
        p.xaxis.axis_label = 'Principal Component %s'  % (components[0]+1)
        p.yaxis.axis_label = 'Principal Component %s'  % (components[1]+1)
        p.hover.tooltips = [("ID", '@ID'), ("category", '@'+group), ("Name", '@case_ID')]
        p.title.text = 'MDS Plot'
        show(p)

    def plt_mu(self, hue=None, thresh=None, RGset=None):
        RGset=RGset if RGset else self.IDAT_Processor.RGset
        if hue is not None:
            _,_, datfr, phenotype =self.IDAT_Processor.getQC(addQC=False,phenotype=hue)
            if phenotype is not None:

                hue=datfr[phenotype]
            else:
                return
        else:
            _,_, datfr =self.IDAT_Processor.getQC(addQC=False, phenotype=None, RGset=None)

        hymin=datfr['mMed'].min()
        vymin=datfr['uMed'].min()
        pheno= robjects.r("pData")(RGset)
        try:
            datfr['ID']=pd.DataFrame(pandas2ri.ri2py(pheno))['ID'].to_numpy()
        except:
            datfr['ID']=pd.DataFrame(pandas2ri.ri2py(robjects.r['as'](pheno,'data.frame'))['ID']).to_numpy()


        X_transformed=datfr


        output_notebook()

        p = figure(
            tools="hover,pan,wheel_zoom,save",
            toolbar_location="above",
            title='MU Plot',
            plot_width=900,
            plot_height=700,
        )

        if (len(X_transformed[phenotype].unique())>2):
            palette = d3['Category20'][len(X_transformed[phenotype].unique())]
            color_map = bmo.CategoricalColorMapper(factors=X_transformed[phenotype].unique(),
                                               palette=palette)

        elif len(X_transformed[phenotype].unique())==2:
            colors=['red', 'green']
            #from bokeh.palettes import brewer
            #colors = brewer["Spectral"][len(X_transformed[group].unique())]

            # Create a map between factor and color.
            colormap = {k: colors[i] for i,k in enumerate(X_transformed[phenotype].unique())}

            # Create a list of colors for each value that we will be looking at.
            X_transformed['colors'] = [colormap[x] for x in X_transformed[phenotype]]

        elif len(X_transformed[phenotype].unique())==1:
            colors=['red']
            #from bokeh.palettes import brewer
            #colors = brewer["Spectral"][len(X_transformed[group].unique())]

            # Create a map between factor and color.
            colormap = {k: colors[i] for i,k in enumerate(X_transformed[phenotype].unique())}

            # Create a list of colors for each value that we will be looking at.
            X_transformed['colors'] = [colormap[x] for x in X_transformed[phenotype]]

        for key in X_transformed[phenotype].unique():

            keys=ColumnDataSource(X_transformed[X_transformed[phenotype]==key])
            if len(X_transformed[phenotype].unique())>2:
                p.scatter(x='mMed', y='uMed', size=10, source=keys ,legend_label=key, color={'field': phenotype, 'transform': color_map})
            else:
                p.scatter(x='mMed', y='uMed', size=10, source=keys ,legend_label=key, color='colors')
        #p.multi_line(xs=[[1, 2, 3], [2, 3, 4]], ys=[[6, 7, 2], [4, 5, 7]],
        #     color=['red','green'])
        p.multi_line(xs=[[thresh,thresh],[hymin,thresh]], ys=[[vymin,thresh],[thresh,thresh]],color='red')
        p.legend.location = "top_left"
        p.legend.click_policy="hide"
        p.hover.tooltips = [ ("category", '@'+phenotype), ("ID", '@ID')]
        p.title.text = 'MU Plot'
        show(p)

    def plt_failedbeads(self,RGset=None, percent=True):

            if RGset is None:
                RGset=self.IDAT_Processor.RGset

            _,nbeads_df=self.IDAT_Processor.beadCount(RGset=RGset)
            nbeads_df=nbeads_df.transpose()
            if percent:
                #nbeads_df['failure_df']=nbeads_df.isnull().sum(axis=1)/len(nbeads_df)*100
                nbeads_df=pd.DataFrame(nbeads_df.isnull().sum(axis=1)/len(nbeads_df.transpose())*100,columns=['failure_df'], index=nbeads_df.index)
            else:
                #nbeads_df['failure_df']=nbeads_df.isnull().sum(axis=1)
                nbeads_df=pd.DataFrame(nbeads_df.isnull().sum(axis=1),columns=['failure_df'], index=nbeads_df.index)


            pheno= robjects.r("pData")(RGset)

            nbeads_df['ID']=self.IDAT_Processor.ri2py_dataframe(pheno, matrix=True)['ID'].to_numpy()

            nbeads_df['disease']=self.IDAT_Processor.ri2py_dataframe(pheno, matrix=True)['disease'].to_numpy()

            nbeads_df['case_ID']=self.IDAT_Processor.ri2py_dataframe(pheno, matrix=True)['case_ID'].to_numpy()


            p = figure(
                    tools="hover,pan,wheel_zoom,save",
                    toolbar_location="above",

                    plot_width=950,
                    plot_height=700,
                    x_range=nbeads_df['ID'].to_numpy(),

                    )
            palette = d3['Category20'][len(nbeads_df['disease'].unique())]
            color_map = bmo.CategoricalColorMapper(factors=nbeads_df['disease'].unique(),
                                           palette=palette)
            legend_it = []
            for key in nbeads_df['disease'].unique():


                keys=ColumnDataSource(nbeads_df[nbeads_df['disease']==key])
                p.vbar(x='ID', top='failure_df', width=1, line_color="white", legend_label=key,
                           fill_color={'field': 'disease', 'transform': color_map}, source=keys)
                #legend_it.append((key, [p]))


            p.xgrid.grid_line_color = None
            p.xaxis.axis_label = "Samples"
            if percent:
                p.yaxis.axis_label = "Percentage"
            else:
                p.yaxis.axis_label = 'Numbers'

            p.x_range.range_padding = 0
            p.x_range.group_padding=0
            p.xaxis.major_label_text_font_size ='3pt'
            p.xaxis.major_label_orientation = 1
            p.legend.location = "top_left"
            #legend = Legend(items=legend_it, location=(0, -60))
            p.legend.click_policy="hide"
            #p.add_layout(legend, 'right')
            p.hover.tooltips = [ ("category", '@disease'), ("ID", '@ID'), ('case_ID','@case_ID')]

            if percent:
                p.title.text = 'Percentage of beadcounts below 3 per sample'

            else:
                p.title.text = 'Number of beadcounts below 3 per sample'

            output_notebook()
            #p = gridplot([s1, s2,s3], ncols=1, toolbar_location="above")
            show(p)

    def plt_meandetP(self, detPcut=0.01, SampleCutoff=0.1, log_scale=True, plot='all', RGset=None ):
        
        RGset= RGset if RGset else self.IDAT_Processor.RGset
          
        detP = robjects.r("""function (rgset, detPcut, cutsamples) {
            detP <- detectionP(rgset)
            numfail <- matrix(colMeans(detP >= detPcut))
            rownames(numfail) <- colnames(detP)
            colnames(numfail) <- "Failed CpG Fraction."
            #print(numfail)
            RemainSample <- which(numfail < cutsamples)

            if(any(numfail >= cutsamples))
            {
                rgset_keeps <- rgset[,RemainSample]

                detP_keeps <- detectionP(rgset_keeps)
                result=list(detP,detP_keeps)
            }
            else{
            result=list(detP)
            }
           return(result)

        }""")(self.IDAT_Processor.RGset,  detPcut, SampleCutoff)

        pheno= robjects.r("pData")(RGset)
        pheno_py=pd.DataFrame(self.IDAT_Processor.ri2py_dataframe(pheno, matrix=True))
        detP_bad_py=None
        if len(pandas2ri.ri2py(detP))==2:
            detP_py=self.IDAT_Processor.ri2py_dataframe(detP[0], matrix=False)

            detP_keep_py=self.IDAT_Processor.ri2py_dataframe(detP[1], matrix=False)

            detP_py=detP_py.transpose()
            detP_py=pd.DataFrame(detP_py.mean(axis=1),columns=['mean'], index=detP_py.index)
            detP_keep_py=detP_keep_py.transpose()
            detP_keep_py=pd.DataFrame(detP_keep_py.mean(axis=1),columns=['mean'], index=detP_keep_py.index)

            if len(np.transpose(pandas2ri.ri2py(detP[0])))!=len(np.transpose(pandas2ri.ri2py(detP[1]))):
                bad_ind=[x for x in detP_py.index.tolist() if x not in detP_keep_py.index.tolist()]
                detP_bad_py=detP_py.loc[bad_ind]
                detP_bad_py=pd.DataFrame(detP_bad_py.mean(axis=1),columns=['mean'], index=detP_bad_py.index)

            detP_py['ID']=pheno_py['ID'].to_numpy()


            detP_py['disease']=pheno_py['disease'].to_numpy()

            detP_py['case_ID']=pheno_py['case_ID'].to_numpy()


            detP_keep_py['ID']=pheno_py[pheno_py['ID'].isin(detP_keep_py.index)]['ID'].to_numpy()

            detP_keep_py['disease']=pheno_py[pheno_py['ID'].isin(detP_keep_py.index)]['disease'].to_numpy()

            detP_keep_py['case_ID']=pheno_py[pheno_py['ID'].isin(detP_keep_py.index)]['case_ID'].to_numpy()


            if detP_bad_py is not None:
                detP_bad_py['ID']=pheno_py[pheno_py['ID'].isin(detP_bad_py.index)]['ID'].to_numpy()

                detP_bad_py['disease']=pheno_py[pheno_py['ID'].isin(detP_bad_py.index)]['disease'].to_numpy()

                detP_bad_py['case_ID']=pheno_py[pheno_py['ID'].isin(detP_bad_py.index)]['case_ID'].to_numpy()


        elif len(pandas2ri.ri2py(detP))==1:
            detP_py=self.IDAT_Processor.ri2py_dataframe(detP[0], matrix=False)
            detP_py['ID']=pheno_py['ID'].to_numpy()


            detP_py['disease']=pheno_py['disease'].to_numpy()

            detP_py['case_ID']=pheno_py['case_ID'].to_numpy()
        
        if plot== 'all' and len(pandas2ri.ri2py(detP))==2:


            s1 = figure(
                tools="hover,pan,wheel_zoom,save",
                toolbar_location="above",

                plot_width=950,
                plot_height=700,
                x_range=detP_py['ID'].to_numpy(),
                #y_range=(test_py['mean'].min(), test_py['mean'].max()),
                y_axis_type="log"
                )
            palette = d3['Category20'][len(detP_py['disease'].unique())]
            color_map = bmo.CategoricalColorMapper(factors=detP_py['disease'].unique(),
                                           palette=palette)

            for key in detP_py['disease'].unique():


                keys=ColumnDataSource(detP_py[detP_py['disease']==key])
                s1.vbar(x='ID', top='mean', bottom=detP_py['mean'].min(), width=1, line_color="white", legend_label=key,
                       fill_color={'field': 'disease', 'transform': color_map}, source=keys)
                s1.line(x=[0,len(detP_py['ID'].to_numpy())], y=[detPcut,detPcut],color='red')

            #p.vbar(x='ID', top='mean', width=1, line_color="white", legend_field='ID',
            #          fill_color={'field': 'disease', 'transform': color_map}, source=source)
            #caption = Label(text='cutoff')
            #s2.add_layout(caption, 'above')
            s1.xgrid.grid_line_color = None
            s1.xaxis.axis_label = "Samplenames"
            s1.yaxis.axis_label = "Mean detP"
            #s1.y_range.start = 0
            s1.x_range.range_padding = 0
            s1.x_range.group_padding=0
            s1.xaxis.major_label_text_font_size ='3pt'
            s1.xaxis.major_label_orientation = 1
            s1.legend.location = "top_left"
            s1.legend.click_policy="hide"
            s1.hover.tooltips = [ ("category", '@disease'), ("ID", '@ID'), ('case_ID','@case_ID')]
            s1.title.text = 'All Samples'


            s2 = figure(
                tools="hover,pan,wheel_zoom,save",
                toolbar_location="above",

                plot_width=950,
                plot_height=700,
                x_range=detP_keep_py['ID'].to_numpy(),
                #y_range=(test_py['mean'].min(), test_py['mean'].max()),
                y_axis_type="log"
                )

            palette = d3['Category20'][len(detP_keep_py['disease'].unique())]
            color_map = bmo.CategoricalColorMapper(factors=detP_keep_py['disease'].unique(),
                                           palette=palette)

            for key in detP_keep_py['disease'].unique():


                keys=ColumnDataSource(detP_keep_py[detP_keep_py['disease']==key])
                s2.vbar(x='ID', top='mean', bottom=detP_keep_py['mean'].min(), width=1, line_color="white", legend_label=key,
                       fill_color={'field': 'disease', 'transform': color_map}, source=keys)
                s2.line(x=[0,len(detP_keep_py['ID'].to_numpy())], y=[detPcut,detPcut],color='red')
            
            s2.xgrid.grid_line_color = None
            s2.xaxis.axis_label = "Samplenames"
            s2.yaxis.axis_label = "Mean detP"            
            s2.x_range.range_padding = 0
            s2.x_range.group_padding=0
            s2.xaxis.major_label_text_font_size ='3pt'
            s2.xaxis.major_label_orientation = 1
            s2.legend.location = "top_left"
            s2.legend.click_policy="hide"
            s2.hover.tooltips = [ ("category", '@disease'), ("ID", '@ID'), ('case_ID','@case_ID')]
            s2.title.text = 'Good Samples'

            if detP_bad_py is not None:
                s3 = figure(
                    tools="hover,pan,wheel_zoom,save",
                    toolbar_location="above",

                    plot_width=950,
                    plot_height=700,
                    x_range=detP_bad_py['ID'].to_numpy(),                    
                    y_axis_type="log"
                    )
                palette = d3['Category20'][len(detP_bad_py['disease'].unique())]
                color_map = bmo.CategoricalColorMapper(factors=detP_bad_py['disease'].unique(),
                                           palette=palette)

                for key in detP_bad_py['disease'].unique():


                    keys=ColumnDataSource(detP_bad_py[detP_bad_py['disease']==key])
                    s3.vbar(x='ID', top='mean', bottom=detP_bad_py['mean'].min(), width=1, line_color="white", legend_label=key,
                           fill_color={'field': 'disease', 'transform': color_map}, source=keys)
                    s3.line(x=[0,len(detP_bad_py['ID'].to_numpy())], y=[detPcut,detPcut],color='red')

                s3.xgrid.grid_line_color = None
                s3.xaxis.axis_label = "Samplenames"
                s3.yaxis.axis_label = "Mean detP"                
                s3.x_range.range_padding = 0
                s3.x_range.group_padding=0
                s3.xaxis.major_label_text_font_size ='3pt'
                s3.xaxis.major_label_orientation = 1
                s3.legend.location = "top_left"
                s3.legend.click_policy="hide"
                s3.hover.tooltips = [ ("category", '@disease'), ("ID", '@ID'), ('case_ID','@case_ID')]
                s3.title.text = 'Good Samples'

                output_notebook()
                p = gridplot([s1, s2,s3], ncols=1, toolbar_location="above")
                show(p)

            else:
                output_notebook()
                p = gridplot([s1, s2], ncols=1, toolbar_location="above")
                show(p)

            return



        if plot== 'all' and len(pandas2ri.ri2py(detP))==1:


            p = figure(
                tools="hover,pan,wheel_zoom,save",
                toolbar_location="above",

                plot_width=950,
                plot_height=700,
                x_range=detP_py['ID'].to_numpy(),                
                y_axis_type="log"
                )
            palette = d3['Category20'][len(detP_py['disease'].unique())]
            color_map = bmo.CategoricalColorMapper(factors=detP_py['disease'].unique(),
                                           palette=palette)

            for key in detP_py['disease'].unique():


                keys=ColumnDataSource(detP_py[detP_py['disease']==key])
                p.vbar(x='ID', top='mean', bottom=detP_py['mean'].min(), width=1, line_color="white", legend_label=key,
                       fill_color={'field': 'disease', 'transform': color_map}, source=keys)
                p.line(x=[0,len(detP_py['ID'].to_numpy())], y=[detPcut,detPcut],color='red')

            p.xgrid.grid_line_color = None
            p.xaxis.axis_label = "Samplenames"
            p.yaxis.axis_label = "Mean detP"            
            p.x_range.range_padding = 0
            p.x_range.group_padding=0
            p.xaxis.major_label_text_font_size ='3pt'
            p.xaxis.major_label_orientation = 1
            p.legend.location = "top_left"
            p.legend.click_policy="hide"
            p.hover.tooltips = [ ("category", '@disease'), ("ID", '@ID'), ('case_ID','@case_ID')]
            p.title.text = 'All Samples'
            output_notebook()

            show(p)
            return



        if plot== 'allsamples':

            p = figure(
                tools="hover,pan,wheel_zoom,save",
                toolbar_location="above",
                plot_width=950,
                plot_height=700,
                x_range=detP_py['ID'].to_numpy(),                
                y_axis_type="log"
                )
            palette = d3['Category20'][len(detP_py['disease'].unique())]
            color_map = bmo.CategoricalColorMapper(factors=detP_py['disease'].unique(),
                                           palette=palette)

            for key in detP_py['disease'].unique():

                keys=ColumnDataSource(detP_py[detP_py['disease']==key])
                p.vbar(x='ID', top='mean', bottom=detP_py['mean'].min(), width=1, line_color="white", legend_label=key,
                       fill_color={'field': 'disease', 'transform': color_map}, source=keys)
                p.line(x=[0,len(detP_py['ID'].to_numpy())], y=[detPcut,detPcut],color='red')

            p.xgrid.grid_line_color = None
            p.xaxis.axis_label = "Samplenames"
            p.yaxis.axis_label = "Mean detP"
            #s1.y_range.start = 0
            p.x_range.range_padding = 0
            p.x_range.group_padding=0
            p.xaxis.major_label_text_font_size ='3pt'
            p.xaxis.major_label_orientation = 1
            p.legend.location = "top_left"
            p.legend.click_policy="hide"
            p.hover.tooltips = [ ("category", '@disease'), ("ID", '@ID'), ('case_ID','@case_ID')]
            p.title.text = 'All Samples'
            output_notebook()

            show(p)
            return


        if plot== 'goodsamples' and len(pandas2ri.ri2py(detP))==2:
            p = figure(
                tools="hover,pan,wheel_zoom,save",
                toolbar_location="above",

                plot_width=950,
                plot_height=700,
                x_range=detP_keep_py['ID'].to_numpy(),                
                y_axis_type="log"
                )
            palette = d3['Category20'][len(detP_keep_py['disease'].unique())]
            color_map = bmo.CategoricalColorMapper(factors=detP_keep_py['disease'].unique(),
                                           palette=palette)

            for key in detP_keep_py['disease'].unique():


                keys=ColumnDataSource(detP_keep_py[detP_keep_py['disease']==key])
                p.vbar(x='ID', top='mean', bottom=detP_keep_py['mean'].min(), width=1, line_color="white", legend_label=key,
                       fill_color={'field': 'disease', 'transform': color_map}, source=keys)
                p.line(x=[0,len(detP_keep_py['ID'].to_numpy())], y=[detPcut,detPcut],color='red')

            p.xgrid.grid_line_color = None
            p.xaxis.axis_label = "Samplenames"
            p.yaxis.axis_label = "Mean detP"            
            p.x_range.range_padding = 0
            p.x_range.group_padding=0
            p.xaxis.major_label_text_font_size ='3pt'
            p.xaxis.major_label_orientation = 1
            p.legend.location = "top_left"
            p.legend.click_policy="hide"
            p.hover.tooltips = [ ("category", '@disease'), ("ID", '@ID'), ('case_ID','@case_ID')]
            p.title.text = 'Good Samples'
            output_notebook()

            show(p)
            return

        if plot== 'goodsamples' and len(pandas2ri.ri2py(detP))==1:
            print('no goodsamples present')
            return


        if plot== 'badsamples' and detP_bad_py is not None:
            p = figure(
                tools="hover,pan,wheel_zoom,save",
                toolbar_location="above",

                plot_width=950,
                plot_height=700,
                x_range=detP_bad_py['ID'].to_numpy(),                
                y_axis_type="log"
                )
            palette = d3['Category20'][len(detP_bad_py['disease'].unique())]
            color_map = bmo.CategoricalColorMapper(factors=detP_bad_py['disease'].unique(),
                                           palette=palette)

            for key in detP_bad_py['disease'].unique():


                keys=ColumnDataSource(detP_bad_py[detP_bad_py['disease']==key])
                p.vbar(x='ID', top='mean', bottom=detP_bad_py['mean'].min(), width=1, line_color="white", legend_label=key,
                       fill_color={'field': 'disease', 'transform': color_map}, source=keys)
                p.line(x=[0,len(detP_bad_py['ID'].to_numpy())], y=[detPcut,detPcut],color='red')

            p.xgrid.grid_line_color = None
            p.xaxis.axis_label = "Samplenames"
            p.yaxis.axis_label = "Mean detP"            
            p.x_range.range_padding = 0
            p.x_range.group_padding=0
            p.xaxis.major_label_text_font_size ='3pt'
            p.xaxis.major_label_orientation = 1
            p.legend.location = "top_left"
            p.legend.click_policy="hide"
            p.hover.tooltips = [ ("category", '@disease'), ("ID", '@ID'), ('case_ID','@case_ID')]
            p.title.text = 'Bad Samples'
            output_notebook()
            show(p)
            return

        if plot== 'badsamples' and detP_bad_py is None:
            print('no badsamples present')
            return

        else:
            print('Please specify which plots you would want: \n'
            'either "all", "allsamples", "goodsamples" or "badsamples" may be specified')
            return



    def plt_curves(self, matrix=None, pheno=None, variant="grouped", group_by=None, nrows=1):

        if pheno is None:
            pheno=self.IDAT_Processor.pheno_orig_py
        import matplotlib.pyplot as plt
        import seaborn as sns
        import copy
        from seaborn import cubehelix_palette
        fig, ax = plt.subplots()
        from collections import OrderedDict
        import matplotlib.pyplot as plt

        cmap=sns.color_palette('muted')

        if type(matrix)==dict:            
            if (variant=="grouped"):

                for i, (key, matr) in enumerate(matrix.items()):                    
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
                    
                    pheno=pheno[pheno['ID'].isin(matr.columns.to_numpy())]
                    matr=copy.deepcopy(matr)
                    c=len(matrix)/nrows
                    if nrows >= 1:
                        ind=i+1
                    plt.subplot(nrows, c, ind)
                    if nrows > 1:
                        if ind==nrows:
                            ind=0                    
                    colors = ['#2300A8', '#00A658','#1f77b4', '#ff7f0e',
                       '#d62728', '#9467bd', '#8c564b',
                      '#e377c2', '#7f7f7f', '#bcbd22', '#17becf','yellow','#070d0d', '#8ffe09']

                    colordict={}
                    for items, values in zip(pheno[group_by].unique().tolist(), colors):
                        colordict[items]=values

                    for col in matr.columns:                        
                        label=pheno[group_by][pheno['ID']==col].to_numpy()[0]                        
                        if len(matr.columns) <= 30:
                            sns.distplot(
                                matr[col], hist=False, rug=False, label=label,
                                color=colordict[label])
                        else:
                            sns.distplot(
                                matr[col], hist=False, rug=False,label=label, color=colordict[label])
                    
                    plt.xlabel('curve-values')
                    plt.ylabel('density')
                    plt.title(key)
                    handles, labels = plt.gca().get_legend_handles_labels()
                    by_label = OrderedDict(zip(labels, handles))
                    plt.legend(by_label.values(), by_label.keys())
                
                plt.tight_layout()
                plt.show()


        else:
            print('Please provide dataframes in key value dict manner, where the key is the label of the plot and the values is the '
                  'dataframe')


    def plot_qc_minfi_report(self, output_dir, diseasegroup='disease', samplenames='ID'):
        check_list=[diseasegroup, samplenames]
        for val in check_list:
            if val not in pandas2ri.ri2py(robjects.r['as'](self.IDAT_Processor.pheno,'data.frame')).columns.tolist():
                print('The pheno sheet does not contain a '+diseasegroup+' or '+samplenames+' column you specified \n'
                     'These are the available column names:')
                print(pandas2ri.ri2py(robjects.r['as'](self.IDAT_Processor.pheno,'data.frame')).columns.tolist())
                return

        self.IDAT_Processor.minfi.qcReport(self.IDAT_Processor.RGset,
                            sampNames=pandas2ri.ri2py(self.IDAT_Processor.pheno)[samplenames],
                            sampGroups=pandas2ri.ri2py(self.IDAT_Processor.pheno)[diseasegroup],
                            pdf = "{}/qcReport.pdf".format(output_dir))


    def screeplot(self, RGset=None, nmax=10):

        if RGset:
            RGset=RGset
        else:
            RGset=self.IDAT_Processor.RGset

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


        x = df.columns.tolist()
        x=[str(element+1) for element in x]
        y1 = df.to_numpy()[0].tolist()

        source = ColumnDataSource(data=dict(x=x, y=y1))
        p = figure(
                tools="hover,pan,wheel_zoom,save",
                toolbar_location="above",
                plot_width=950,
                plot_height=700,
                x_range=x
                )

        palette = d3['Category20'][len(x)]
        p.vbar(x='x', top='y', width=1, line_color="white", source=source, legend_field="x",
              fill_color=factor_cmap('x', palette=palette, factors=x))

        p.xgrid.grid_line_color = None
        p.xaxis.axis_label = "Principile components"
        p.yaxis.axis_label = "Proportion of explained variance"
        p.yaxis.formatter = NumeralTickFormatter(format='0 %')
        p.legend.orientation = "horizontal"
        p.legend.location = "top_center"        
        p.legend.title = 'Principle components'
        p.hover.tooltips = [ ("Explained variance", '@y')]
        p.title.text = 'Screeplot'
        output_notebook()

        show(p)


    def plt_covariates(self, pheno=None, matrix=None, pcs=2, ids=False):

        if not matrix or not pheno:
            print('Please input data')
            return

        cxy = robjects.r("""function (pcs, matrix, pheno,ids) {
             if (ids!=FALSE)rownames(pheno)<-pheno[,ids]

             keep<-intersect(colnames(matrix),rownames(pheno))
             pheno<- pheno[keep,]
             matrix<- matrix[,keep]
             #pheno <- pheno[rownames(pheno) %in% colnames(matrix) ,]

             df <- apply(pheno, 2, function(x) as.numeric(factor(x)))

             keep <- apply(df, 2, sd) > 0

             df <- df[ , keep]

             library(irlba)
             #pc <- prcomp(t(matrix), rank=pcs)
             pc <- prcomp_irlba(t(matrix), n=pcs)

             cxy <- round(cor(pc$x, scale(df)),2)

             return(cxy)
            }""")(pcs, matrix, pheno,ids)


        corr=self.IDAT_Processor.ri2py_dataframe(cxy, matrix=False)

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


    def plot_qc_metrics(self, output_dir):
        """Plot QC results from ENmix pipeline and possible minfi. Still experimental.
        Parameters
        ----------
        output_dir
            Where to store plots."""
        self.IDAT_Processor.enmix.plotCtrl(self.IDAT_Processor.RGset)
        grdevice = importr("grDevices")
        geneplotter = importr("geneplotter")
        base = importr('base')
        anno=self.IDAT_Processor.minfi.getAnnotation(self.IDAT_Processor.RGset)
        anno_py = pandas2ri.ri2py(robjects.r['as'](anno,'data.frame'))
        beta_py = pandas2ri.ri2py(self.IDAT_Processor.beta)
        beta1=numpy2ri.py2ri(beta_py[anno_py["Type"]=="I"])
        beta2=numpy2ri.py2ri(beta_py[anno_py["Type"]=="II"])
        grdevice.jpeg(output_dir+'/dist.jpg',height=900,width=600)
        base.par(mfrow=robjects.vectors.IntVector([3,2]))
        self.IDAT_Processor.enmix.multidensity(self.IDAT_Processor.beta, main="Multidensity")
        self.IDAT_Processor.enmix.multifreqpoly(self.IDAT_Processor.beta, xlab="Beta value")
        self.IDAT_Processor.enmix.multidensity(beta1, main="Multidensity: Infinium I")
        self.IDAT_Processor.enmix.multifreqpoly(beta1, main="Multidensity: Infinium I", xlab="Beta value")
        self.IDAT_Processor.enmix.multidensity(beta2, main="Multidensity: Infinium II")
        self.IDAT_Processor.enmix.multifreqpoly(beta2, main="Multidensity: Infinium II", xlab="Beta value")
        grdevice.dev_off()
        self.IDAT_Processor.minfi.qcReport(self.IDAT_Processor.RGset, pdf = "{}/qcReport.pdf".format(output_dir))
        self.IDAT_Processor.minfi.mdsPlot(self.IDAT_Processor.RGset)
        self.IDAT_Processor.minfi.densityPlot(self.IDAT_Processor.RGset, main='Beta', xlab='Beta')


    def plot_original_qc(self, output_dir):
        """Plot QC results from ENmix pipeline and possible minfi. Still experimental.
        Parameters
        ----------
        output_dir
            Where to store plots."""
        self.IDAT_Processor.preprocessRAW()
        self.IDAT_Processor.return_beta()
        self.IDAT_Processor.get_beta()
        self.IDAT_Processor.plot_qc_metrics(output_dir)