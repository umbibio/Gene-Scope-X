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
GeneScopeX is a web application that provides a comprehensive and user-friendly interface for analyzing single-cell gene expression data using the scanpy library. The application takes anndata, performs preprocessing, displays PCA and map data, facilitates manual and automatic clustering, and performs differential gene expression analysis. Additionally, GeneScopeX provides a wide range of visualization options like dendrogram, violin plot, heat map, matrix plot, and more, enabling users to gain deep insights into their data.

The application is designed to be highly flexible and adaptable, allowing users to customize their analysis and visualization workflows to fit their specific research needs. With GeneScopeX, researchers can quickly and easily visualize and analyze gene expression patterns across different cell types, identify differentially expressed genes, perform clustering and dimensionality reduction, and explore complex relationships within their data.

Overall, GeneScopeX is an all-in-one solution for single-cell analysis that combines powerful data processing and analysis tools with intuitive visualization options, making it an ideal choice for researchers and analysts working with large-scale single-cell gene expression datasets.
''')),
                # dbc.Col(
                #     dbc.Carousel(items=[
                #         {"key": "1", "src": "/babesiasc/assets/newplot-expression.png"},
                #         {"key": "2", "src": "/babesiasc/assets/newplot-pstime.png"},
                #         {"key": "3", "src": "/babesiasc/assets/network-graph.png"},
                #     ], interval=5000, className="carousel-fade", indicators=False)
                # )
                ])),
        ],),
    ), class_name="mb-4 mt-4"),
    dbc.Row([
#         dbc.Col(dbc.Card([
#             dbc.CardHeader(html.H4("Citation")),
#             dbc.CardBody(dcc.Markdown('''
# Yasaman Rezvani*, Caroline D Keroack*, Brendan Elsworth, Argenis Arriojas, Marc-Jan Gubbels, Manoj T Duraisingh, Kourosh Zarringhalam, "Single cell transcriptional atlas of Babesia species reveals coordinated progression of transcriptomes and identifies conserved and species-specific transcriptional profiles", 2022, BioRxiv
# ''')),
#         ],),),
        dbc.Col(dbc.Card([
            dbc.CardHeader(html.H4("Contact")),
            dbc.CardBody(dcc.Markdown('''
For questions or comments please contact:

Kourosh.Zarringhalam at umb dot edu

Kourosh.Zarringhalam@umb.edu
''')),
        ],),),
    ]),
]
