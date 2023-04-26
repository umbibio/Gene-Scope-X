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
                        {"key": "1", "src": "assets/newplot-scatter.png"},
                        {"key": "2", "src": "assets/newplot-clustering.png"},
                        {"key": "3", "src": "assets/newplot-heatmap.png"},
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
            dbc.CardBody(dcc.Markdown('''
For questions or comments please contact:

Divya.Thota at umb dot edu

Divya.Thota@umb.edu
''')),
        ],),),
    ]),
]
