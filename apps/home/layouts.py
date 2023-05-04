from dash import dcc
from dash import html

import dash_bootstrap_components as dbc
menu = None

body = [
    dbc.Row(dbc.Col(
        dbc.Card([
            dbc.CardHeader(html.H2("GeneScopeX")),
            dbc.CardBody(
                dbc.Row([
                    dbc.Col(dcc.Markdown('''
GeneScopeX is a web application that provides a comprehensive and user-friendly interface for analyzing single-cell gene expression data using the scanpy library. To explore the data click on one of the tabs above.

1. **Expression**: Explore expression of genes projected on PCA or UMAP coordinates, perform clustering, identify differentially expressed genes and visualize the representation on various plots.
2. **Download**: Provides various formats of sample files to download and run the expression.
3. **Merge**: Merges multiple 10x Genomic files into a single file and download the resulting h5ad file to run the expression.
''')),
                dbc.Col(
                    dbc.Carousel(items=[
                        {"key": "1", "src": "assets/images/newplot-scatter.png"},
                        {"key": "2", "src": "assets/images/newplot-clustering.png"},
                        {"key": "3", "src": "assets/images/newplot-heatmap.png"},
                    ], interval=5000, className="carousel-fade", indicators=False)
                )
                ])),
        ],),
    ), class_name="mb-4 mt-4"),
    dbc.Row([
        dbc.Col(dbc.Card([
            dbc.CardHeader(html.H4("Citation")),
            dbc.CardBody(dcc.Markdown('''
Jialu Hu1, Mengjie Chen2, and Xiang Zhou, "Effective and scalable single-cell data alignment with non-linear canonical correlation analysis", 2021

F. Alexander Wolf, Philipp Angerer & Fabian J. Theis, "SCANPY: large-scale single-cell gene expression data analysis",2018
''')),
        ],),),
        dbc.Col(dbc.Card([
            dbc.CardHeader(html.H4("Contact")),
            dbc.CardBody([dbc.Row(dcc.Markdown('''
For questions or comments please contact:

Divya.Thota at umb dot edu

Divya.Thota@umb.edu
''')), 
dbc.Row(dbc.CardLink(
    html.Img(src='assets/images/github-logo.png', className="github_logo"), 
    href="https://github.com/umbibio/Gene-Scope-X", class_name="githublink"),
    class_name = "githubrow")]),
        ],),),
    ],class_name="mb-4 mt-4"),
    dbc.Row(dbc.Col(
        dbc.Card([
            dbc.CardHeader(html.H2("Features")),
            dbc.CardBody(
dcc.Markdown(''' The main features of the GeneScopeX application are:

1.	File Upload: The application allows the user to upload 10x Genomic files in different formats such as .h5ad, .mtx and .tsv compressed into a .zip,  and .h5
2.	Preprocessing: The uploaded files undergo various preprocessing steps such as normalization, Logarithmization, Identify highly variable genes, Perform principle component analysis, Compute neighbors, Uniform Manifold Approximation and Projection for Dimension Reduction, etc. to ensure high-quality analysis.
3.	Clustering: GeneScopeX allows the user to perform unsupervised clustering analysis on the preprocessed data using Leiden clustering algorithm or manual clustering.
4.	Differential Gene Analysis: The application allows the user to perform differential gene analysis on the clustered data to identify genes that are differentially expressed between the different conditions.
5.	Data Visualization: GeneScopeX provides various visualization tools to help the user explore and analyze the data. The application provides different types of plots such as scatterplots, heatmaps, violin plots, etc.
6.	Data Download: The application allows the user to download the differential gene alaysis data and generated plots.
7.	Merge Multiple Files: GeneScopeX also provides an option to merge multiple 10x Genomic files into a single file and download it, making it easier to analyze and compare data from different experiments.
'''),
            ),
        ],),
    ),class_name="mb-4 mt-4"),
    dbc.Row(dbc.Col(
        dbc.Card([
            dbc.CardHeader(html.H2("Usage")),
            dbc.CardBody(dcc.Markdown(''' To use the application, follow these steps:

1.	If you have a 10x Genomic file, skip to step 4.
2.	Click on “Download” Tab.
3.	Download any sample file.
4.	If you want to view the data from different experiments, click on “Merge” tab. Otherwise, skip to step 8.
5.	Drag and drop or select the various 10x Genomic files to merge.
6.	Once the upload is complete, click on the "Merge Files" button to start merging the uploaded files. The application will show the progress of the merging process in a progress bar.
7.	Once the merging process is complete, you can download the merged file by clicking on the "Download Merged h5ad File" button.
8.	Click on “Expression” Tab.
9.	Drag and drop or select the downloaded 10x Genomic file to upload in the upload area.
10.	Select the preprocessing parameters and click on “OK”. The application will show the progress of the preprocessing steps.
11.	Once the Preprocessing is complete, Choose the type of plot PCA/UMAP and the gene and click on “Plot Current Attributes” button. The scatter plot for the selected attributes will be displayed. 
12.	Click on one of the clustering buttons:
a.	Manual Clustering: Click on “Add Cluster” until the required number of clusters are visible and choose Lasso Select tool on the scatterplot.
i.	Manually select the points for the first cluster, name the cluster and click on “OK” below the name. Wait until the loading message disappears.
ii.	Repeat the previous step for all the other clusters.
iii.	Click on “Compare” button.
b.	Automatic Clustering
13.	The differential gene analysis for the clustered data will be processed and the scatter plot will be updated to represent the clusters. 
14.	Click on “Download Differential Gene Expression” button to download the .csv file containing the scores, logfoldchanges, pvals etc.
15.	Choose various visualization options and the number of Genes and click on “Plot Graph” button to view the graphs.
16.	Click on the download icon on the plot to download the plots.
'''),
                    
            ),
        ],),
    ),class_name="mb-4 mt-4"),  
]
