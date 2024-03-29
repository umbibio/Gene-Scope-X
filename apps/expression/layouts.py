from dash import dcc
from dash import html

import plotly.graph_objects as go
import dash_bootstrap_components as dbc
import uuid
import dash_uploader as du



def get_upload_component(id):
    return du.Upload(
        id=id,
        text="Drag and Drop or Select Files",
        cancel_button=True,
        pause_button=True,
        filetypes=['h5ad', 'zip', 'h5'],
        upload_id=uuid.uuid4(),  # Unique session id
        max_files=1,
        default_style = {'min-height': '150px'}
    )

menu = [
    html.Div([dbc.Spinner([html.Div(id='adata-div', children=[get_upload_component(id='dash-uploader'), 
              dcc.Markdown('''
The supported file formats are:
* h5ad file
* h5 file
* zip file compressing the folder containing matrix, features and barcode files''', id='supported-file'),
              html.Div(id='preprocessing-div', children=[
                  html.Br(),
                  html.P("PCA"),
                  html.P("SVD solver to use:"),
                  dcc.Dropdown(id='svd_solver',
                               options=['arpack', 'randomized', 'auto', 'lobpcg'], value='arpack'),
                  html.Hr(),
                  html.P("UMAP"),
                  html.P("Number of Neighbors:"),
                  dcc.Input(id="n_neighbors", type="number"),
                  html.Br(), html.Br(),
                  html.Button("OK", id="preprocessing-button"),
              ], style={'display': 'none'}),])],id='adata-div-spinner', spinner_style={'display': 'display'}, type='border', fullscreen=False, color='primary', delay_hide=100,),
              dcc.Store(id='adata-path'),
              dcc.Store(id='session-id'),
              dcc.Store(id='genes'),
              dcc.Store(id='obsmdf'),
              dcc.Store(id='clustertype'),
              dcc.Store(id='grouplen'),
              dcc.Store(id='manual-cluster-dict'),
              html.Div(id = 'scatter-plot', children = [
                  html.Div(id='scatter-plot-div', children = [html.P("Type: "),
                    dcc.Dropdown(id='type',options=['PCA', 'UMAP'], value='PCA'),
                    dbc.RadioItems(id='2d3d-radio',
                                   options=[{'label': x, 'value': x} for x in ['2D', '3D']],
                                   value='2D',
                                   inline=True),
                    html.P("Gene: "),
                    dcc.Dropdown(id='gene'),
                    html.Br(),]),
                  html.Button("Plot Current Attributes", id="plot-button"),
                  html.Hr(),], style={'display':'none'}),
              html.Div(id='cluster-option-div', children =[
                  html.Button("Automatic Clustering", id="auto-cluster-button"),
                  html.Br(), html.Br(),
                  html.Button("Manual Clustering", id="manual-cluster-button"),], style={'display': 'none'}),
              html.Div(id='manual-cluster-div', style={'display': 'none'}, children =[
                            html.Button("Add Cluster", id="add-cluster", n_clicks=0),
                            html.Div(id='cluster-container', children=[]),
                            html.Br(),
                            html.Div(id = 'manual-progress-text'),
                            html.Button("Compare", id="compare-button")]),
              html.Div(id='visualization-div', style={'display': 'none'}, children=[
                        html.Button("Download Differential Gene Expression", id="download-button"),
                        dcc.Download(id="download"),
                        html.Br(),
                        html.P("Choose graph:"),
                        dcc.Dropdown(id='plotComboBox',
                                       options=['Dendrogram', 'Gene Ranking', 'Dot Plot', 'Violin', 'Stacked Violin',
                                                'Matrix Plot', 'Heatmap', 'Tracksplot']),
                        html.P("Number of Genes:"),
                        dcc.Input(id="n_genes", type="number"),
                        html.Br(),
                        html.Div(id = "violinplot-radio-buttons", style={'display': 'none'}),
                        html.Button("Plot Graph", id="visualization-plot-button"),
              ]),
            ],),
    html.Div(id = 'progress-text'),
]


body = [
    dbc.Spinner([dcc.Graph(
        id='expression-graph',
        figure={'layout': {'height': 700}},
    )]),
]