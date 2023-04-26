from dash.dependencies import Input, Output, State
from app import app
import dash
from dash import dcc

@app.callback(Output("download-h5ad", "data"),
              Input('download-h5ad-button','n_clicks'),
              prevent_initial_call=True)
def anndata_download(n):
    if n is None:
        return dash.no_update
    else:
        return dcc.send_file("assets/adata.h5ad")

@app.callback(Output("download-10x", "data"),
              Input('download-10x-button','n_clicks'),
              prevent_initial_call=True)
def anndata_download(n):
    if n is None:
        return dash.no_update
    else:
        return dcc.send_file("assets/filtered_feature_bc_matrix.zip")

@app.callback(Output("download-h5", "data"),
              Input('download-h5-button','n_clicks'),
              prevent_initial_call=True)
def anndata_download(n):
    if n is None:
        return dash.no_update
    else:
        return dcc.send_file("assets/filtered_feature_bc_matrix.h5")