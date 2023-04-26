from dash import dcc
from dash import html

import dash_bootstrap_components as dbc
import dash_uploader as du
import uuid
def get_upload_component(id):
    return du.Upload(
        id=id,
        text="Drag and Drop or Select Files",
        cancel_button=True,
        pause_button=True,
        filetypes=['h5ad', 'zip', 'h5'],
        upload_id=uuid.uuid4(),  # Unique session id
        max_files=5,
        default_style = {'min-height': '150px'}
    )

menu = None
body = [
    dcc.Store(id='uploaded-files'),
    dcc.Store(id='merge-file-path'),
    dbc.Row(dbc.Col(
        dbc.Card([
            dbc.CardHeader(html.H2("Merge Multiple 10x Genomic Files")),
            dbc.CardBody(
                dbc.Row([
                    dcc.Markdown('''
The Merge feature is used to merge multiple 10x Genomic files into a single file. The script uses the scanpy and VIPCCA packages for preprocessing and cell type annotation. The merged file can be used for downstream analysis, such as clustering, differential expression analysis, and visualization. 
                    
The feature takes as input a list of 10x Genomic files, which can be in either .h5ad, .zip, or .h5 format. It preprocesses the input files by removing low-quality cells, filtering genes, and normalizing gene expression values. It then uses VIPCCA to perform cell type annotation and RNA velocity estimation. Finally, it saves the merged file in .h5ad format, which can be downloaded and used for further analysis.''')
                ])
            ),
        ],),
    ),class_name="mb-4 mt-4"),
    dbc.Row( dbc.Col(
        dbc.Card([
            dbc.CardHeader(html.H2("Merge files")),
            dbc.CardBody(
                dbc.Row([
                    dbc.Col([get_upload_component(id='merge-uploader'),
                             html.Br(),
                             html.Button("Merge files", id="merge-button", disabled=True),
                             dbc.Progress(id="merge-progress-bar", value="0", label="", style={'display': 'none'}),]),
                    dbc.Col([html.Button("Download merged h5ad file", id="download-merge-button", disabled=True),
                             dcc.Download(id="download-merge")]),
                ])
            ),
        ], ),
    ),class_name="mb-4 mt-4"),
]
