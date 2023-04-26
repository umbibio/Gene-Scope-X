import dash
import dash_uploader as du
from dash import html, dcc, ctx, ALL, DiskcacheManager, CeleryManager
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output, State

from app import app
from skimage import io
from pathlib import Path
from zipfile import ZipFile

import os
import json
import matplotlib
import scanpy as sc
import pandas as pd
import plotly.express as px
import matplotlib.pyplot as plt
import plotly.graph_objects as go

"""
    Configuration of background callbacks

    A DiskCache backend that runs callback logic in a separate process and stores the results to disk 
    using the diskcache library.This is the easiest backend to use for local development, but is not 
    recommended for production.

    A Celery backend that runs callback logic in a Celery worker and returns results to the Dash app 
    through a Celery broker like Redis. This is recommended for production as, unlike Disk Cache, 
    it queues the background callbacks, running them one-by-one in the order that they were received 
    by dedicated Celery worker(s). Celery is a widely adopted, production-ready job queue library.
"""
if 'REDIS_URL' in os.environ:
    # Use Redis & Celery if REDIS_URL set as an env variable
    from celery import Celery

    celery_app = Celery(__name__, broker=os.environ['REDIS_URL'], backend=os.environ['REDIS_URL'])
    background_callback_manager = CeleryManager(celery_app)

else:
    # Diskcache for non-production apps when developing locally
    import diskcache

    cache = diskcache.Cache(r"./cache")
    background_callback_manager = DiskcacheManager(cache)

matplotlib.pyplot.switch_backend('Agg')
UPLOAD_FOLDER_ROOT = r"Uploads"
du.configure_upload(app, UPLOAD_FOLDER_ROOT)


def write_adata(adata, adata_path):
    """
    Write the AnnData object to a file.

    Args:
        adata (AnnData): AnnData object to write.
        adata_path (Path): Path of the file to write to.
    """
    adata.write(adata_path)


@du.callback(
    [Output('adata-path', 'data'),
     Output('session-id', 'data'),
     Output('preprocessing-div', 'style'),
     Output('supported-file', 'style'),
     Output('dash-uploader', 'disabled'),
     Output('dash-uploader', 'disabledMessage')],
    id='dash-uploader',
)
def get_adata(status: du.UploadStatus):
    """
        Callback function to get the path of the uploaded file, save it in H5AD format and return necessary outputs for the Dash app.

        Args:
            filename (str): The path of the uploaded file.

        Returns:
        tuple: A tuple of the following:
            - (str) Path of the generated AnnData file.
            - (str) session id.
            - (dict) Display style for the preprocessing div.
            - (bool) Boolean flag to disable the dash-uploader component.
            - (str) Message to display when the dash-uploader component is disabled.
    """
    try:
        filename = str(status.uploaded_files[0])
        print('File uploaded in :', filename)
        if '\\' in filename:
            filename = filename.replace('\\', '/')
        folderpath = filename.rsplit('/', 1)[0]
        if filename.endswith('h5ad'):
            adata = sc.read(filename)
        elif filename.endswith('zip'):
            with ZipFile(filename, 'r') as zip:
                name = zip.infolist()[0].filename
                zip.extractall(folderpath)
            adata = sc.read_10x_mtx(
                Path(folderpath + '/' + name),
                var_names='gene_symbols',
                cache=False)
        elif filename.endswith('h5'):
            adata = sc.read_10x_h5(filename)
        write_adata(adata, folderpath + '/adata.h5ad')
        return folderpath + '/adata.h5ad', filename.split('/')[1], {'display': 'block'},{'display': 'none'}, True, 'Uploaded: ' + \
               filename.rsplit('/', 1)[1]
    except Exception as e:
        print(f"Error occurred while processing file: {e}")
        return None, None, {'display': 'none'},{'display': 'block'}, False, 'Error: Please upload a valid file'


def make_unique(adata):
    """
        Modifies the variable names in the input AnnData object to ensure uniqueness.

        This function modifies the variable names in the input AnnData object to ensure that they are unique.
        If any of the variable names are duplicates, a unique suffix will be appended to each name.

        Parameters:
            adata (anndata.AnnData): The input AnnData object.

        Returns:
            anndata.AnnData: The modified AnnData object with unique variable names.
    """
    print("Preprocessing annData file...")
    adata.var_names_make_unique()
    return adata


def normalize_adata(adata):
    """
        Normalize the gene expression data in the input AnnData object.

        This function normalizes the gene expression data in the input AnnData object by scaling the expression
        levels of each cell so that they sum to a fixed value (1e4 by default). This scaling factor is applied
        to all genes in each cell to correct for differences in sequencing depth between cells.

        Parameters:
            adata (anndata.AnnData): The input AnnData object.

        Returns:
            anndata.AnnData: The modified AnnData object with normalized gene expression data.
    """
    print("Normalize annData file...")
    sc.pp.normalize_total(adata, target_sum=1e4)
    return adata


def logarithmize_total(adata):
    """
        Transform the gene expression data in the input AnnData object to logarithmic scale.

        This function applies a natural logarithm transformation to the gene expression data in the input AnnData object.
        This transformation is typically applied to gene expression data to convert it to a logarithmic scale, which can
        make the data easier to work with and interpret.

        Parameters:
            adata (anndata.AnnData): The input AnnData object.

        Returns:
            anndata.AnnData: The modified AnnData object with logarithmically transformed gene expression data.
    """
    print("Logarithmize annData file...")
    sc.pp.log1p(adata)
    return adata


def highly_variable_genes(adata):
    """
        Identify highly variable genes in the input AnnData object.

        This function identifies the highly variable genes in the input AnnData object by calculating the mean expression
        and the dispersion (variance/mean) for each gene across all cells. The genes that fall outside the specified
        thresholds for mean expression and dispersion are considered highly variable and are marked as such in the
        `var.highly_variable` attribute of the AnnData object.

        Parameters:
            adata (anndata.AnnData): The input AnnData object.

        Returns:
            anndata.AnnData: The modified AnnData object with identified highly variable genes.
    """
    print("Identify highly variable genes...")
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    return adata


def pca_adata(adata, svd_solver):
    """
        Perform principal component analysis (PCA) on the input AnnData object.

        This function performs principal component analysis (PCA) on the gene expression data in the input AnnData object,
        reducing the dimensionality of the data by identifying the major sources of variation across cells. The PCA
        components are stored in the `obsm['X_pca']` attribute of the AnnData object.

        Parameters:
            adata (anndata.AnnData): The input AnnData object.
            svd_solver (str, optional): The solver to use for calculating the PCA. Must be one of 'arpack' or 'svd'.
                Defaults to 'arpack'.

        Returns:
            anndata.AnnData: The modified AnnData object with PCA components stored in `obsm['X_pca']`.
    """
    print("Perform principle component analysis...")
    sc.tl.pca(adata, svd_solver=svd_solver)
    return adata


def compute_neighbors(adata, n_neighbors):
    """
        Compute nearest neighbors for the input AnnData object.

        This function computes the nearest neighbors for each cell in the input AnnData object, based on the gene
        expression similarity between cells. The nearest neighbor information is stored in the `uns['neighbors']`
        attribute of the AnnData object.

        Parameters:
            adata (anndata.AnnData): The input AnnData object.
            n_neighbors (int, optional): The number of neighbors to compute for each cell. Defaults to 10.

        Returns:
            anndata.AnnData: The modified AnnData object with computed nearest neighbor information.
    """
    print("Compute neighbors...")
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, knn=True, n_pcs=40)
    return adata


def umap_adata(adata, adata_path):
    """
        Compute UMAP embedding for the input AnnData object.

        This function computes the Uniform Manifold Approximation and Projection (UMAP) embedding for the input AnnData
        object, reducing the dimensionality of the data to two dimensions for visualization purposes. The UMAP
        coordinates for each cell are stored in the `obsm['X_umap']` attribute of the AnnData object. Additionally, the
        modified AnnData object is saved to the specified `adata_path`, and the UMAP coordinates are returned as a Pandas
        DataFrame.

        Parameters:
            adata (anndata.AnnData): The input AnnData object.
            adata_path (str): The file path to save the modified AnnData object.

        Returns:
            pandas.DataFrame: The UMAP coordinates for each cell in the input AnnData object.
    """
    print("Uniform Manifold Approximation and Projection for Dimension Reduction...")
    sc.tl.umap(adata)
    adata.write(adata_path)
    obsmdf = adata.obsm.to_df()
    return obsmdf


@app.callback(
    Output('genes', 'data'),
    Output('obsmdf', 'data'),
    Output('adata-div', 'style'),
    Output('gene', 'options'),
    Output('gene', 'value'),
    Output('scatter-plot', 'style'),
    State('adata-path', 'data'),
    Input('preprocessing-button', 'n_clicks'),
    State('svd_solver', 'value'),
    State('n_neighbors', 'value'),
    State('scatter-plot-div', 'children'),
    background=True,
    manager=background_callback_manager,
    running=[
        (Output("progress-text", "style"), {'display': 'block'}, {'display': 'none'}),
        (Output("adata-div-spinner", "style"), {'display': 'block'}, {'display': 'none'}),
    ],
    progress=Output("progress-text", "children"),
    prevent_initial_call=True)
def preprocess_adata(set_progress, adata_path, n, svd_solver, n_neighbors, children):
    """
        Preprocess the annData file and update the relevant components on the UI.

        Parameters:
            set_progress (callable): Function used to update the progress on the UI.
            adata_path (str): Path to the annData file.
            n (int): Number of times the preprocessing button has been clicked.
            svd_solver (str): The string value for principal component analysis.
            n_neighbors (int): Number of neighbors to consider for computing the UMAP.
            children (list): Children of the scatter plot div.

        Returns:
            tuple: A tuple of outputs for the callback.
        """
    try:
        if n is None:
            return dash.no_update
        else:
            adata = sc.read(adata_path)
            genes = pd.DataFrame(adata.var.gene_ids)
            gene_ids = genes['gene_ids'].to_numpy()
            set_progress("Preprocessing annData file...")
            adata = make_unique(adata)
            set_progress("Normalize annData file...")
            adata = normalize_adata(adata)
            set_progress("Logarithmize annData file...")
            adata = logarithmize_total(adata)
            set_progress("Identify highly variable genes...")
            adata = highly_variable_genes(adata)
            set_progress("Perform principle component analysis...")
            adata = pca_adata(adata, svd_solver)
            set_progress("Compute neighbors...")
            if n_neighbors:
                adata = compute_neighbors(adata, n_neighbors)
            else:
                n_neighbors = 10
                adata = compute_neighbors(adata, n_neighbors)
            set_progress("Uniform Manifold Approximation and Projection for Dimension Reduction...")
            obsmdf = umap_adata(adata, adata_path)
            set_progress("Preprocessing Complete...")
        return genes.to_json(), obsmdf.to_json(), {'display': 'none'}, gene_ids, gene_ids[0], {'display': 'block'}
    except Exception as e:
        print(f"Error occurred while preprocessing annData file: {e}")
        set_progress(f"Error occurred while preprocessing annData file: {e}")
        return None, None, {'display': 'none'}, [], '', {'display': 'none'}


@app.callback(
    Output('2d3d-radio', 'style'),
    Input('type', 'value'),
    prevent_initial_call=True
)
def update_2d3d_visibility(value):
    """
        Updates the visibility style of the 2D/3D radio button based on the selected value of the "type" dropdown.

        Args:
            value (str): The selected value of the "type" dropdown.

        Returns:
            dict: A dictionary with the style to be applied to the 2D/3D radio button.
    """
    try:
        if value == 'PCA':
            style = {'display': 'block'}
        else:
            style = {'display': 'none'}
        return style
    except Exception as e:
        print(f"Error in update_2d3d_visibility: {e}")
        return {'display': 'none'}


@app.callback (Output('expression-graph', 'figure', allow_duplicate=True),
               Output('cluster-option-div', 'style', allow_duplicate=True),
               Output('visualization-div', 'style', allow_duplicate=True),
               Output('manual-cluster-div', 'style', allow_duplicate=True),
               Output('clustertype', 'data', allow_duplicate=True),
               Output('manual-cluster-dict', 'data', allow_duplicate=True),
               Input('plot-button', 'n_clicks'),
               State('adata-path', 'data'),
               State('genes', 'data'),
               State('2d3d-radio', 'value'),
               State('obsmdf', 'data'),
               State('type', 'value'),
               State('gene', 'value'),
               prevent_initial_call=True)
def initial_scatter_plot(n, adata_path, genes, pca_radio, obsmdf, type, gene):
    """
        Update the scatter plot when the user selects a gene and clicks the "Plot" button.

        Parameters:
            n_clicks (int): Number of times the "Plot" button has been clicked.
            adata_path (str): Path to the AnnData file.
            genes_data (str): JSON-encoded pandas DataFrame containing gene IDs.
            pca_radio_value (str): Value of the 2D/3D radio button.
            obsmdf_data (str): JSON-encoded pandas DataFrame containing UMAP/PCA coordinates.
            type_value (str): Type of dimensionality reduction (PCA/UMAP).
            gene_value (str): Selected gene ID.

        Returns:
            tuple: A tuple containing the updated figure, the style of the cluster options, and the cluster type.
        """
    try:
        if n is None:
            return dash.no_update
        else:
            clustertype = ""
            adata = sc.read(adata_path)
            obsmdf = pd.read_json(obsmdf)
            dfadata = adata.to_df()
            genes = pd.read_json(genes)
            genename = genes.index[genes['gene_ids'] == gene].tolist()[0]
            if type == 'PCA' and pca_radio == '2D':
                fig = px.scatter(x=obsmdf.X_pca1, y=obsmdf.X_pca2, labels={'x': 'PC1', 'y': 'PC2'}, template="simple_white",
                                 title=genename, color=dfadata[genename])
            elif type == 'PCA' and pca_radio == '3D':
                fig = px.scatter_3d(x=obsmdf.X_pca1, y=obsmdf.X_pca2, z=obsmdf.X_pca3,
                                    labels={'x': 'PC1', 'y': 'PC2', 'z': 'PC3'}, template="simple_white",
                                    title=genename, color=dfadata[genename])
                fig.update_traces(marker_size=2)
                fig.update_layout(scene=dict(xaxis=dict(showgrid=True),
                                             yaxis=dict(showgrid=True),
                                             zaxis=dict(showgrid=True)
                                             ))
            elif type == 'UMAP':
                fig = px.scatter(x=obsmdf.X_umap1, y=obsmdf.X_umap2, labels={'x': 'UMAP1', 'y': 'UMAP2'},
                                 template="simple_white", title=genename, color=dfadata[genename])
            figure = go.Figure(data=fig)
            return figure, {'display': 'block'},{'display': 'none'}, {'display': 'none'}, clustertype, {}
    except Exception as e:
        print(f"Error in initial_scatter_plot: {e}")
        return {}, {'display': 'none'},{'display': 'none'}, {'display': 'none'}, clustertype, {}


@app.callback(Output('expression-graph', 'figure', allow_duplicate=True),
               Output('cluster-option-div', 'style', allow_duplicate=True),
               Output('visualization-div', 'style', allow_duplicate=True),
               Output('clustertype', 'data', allow_duplicate=True),
               Output('grouplen', 'data', allow_duplicate=True),
               Input('auto-cluster-button', 'n_clicks'),
               State('adata-path', 'data'),
               State('session-id', 'data'),
               State('type', 'value'),
               prevent_initial_call=True)
def auto_scatter_plot(n, adata_path, session_id, type):
    """
        Generates an auto-clustered scatter plot and saves the plot image file.

        Args:
            n (int): The number of times the auto-cluster button was clicked.
            adata_path (str): The path to the AnnData object.
            session_id (str): The ID of the session.
            type (str): The type of plot ('PCA' or 'UMAP').

        Returns:
            tuple: A tuple containing the plotly figure, cluster option div style, visualization div style,
                cluster type data, and group length data.

    """
    try:
        if n is None:
            return dash.no_update
        else:
            adata = sc.read(adata_path)
            clustertype = "auto"
            adata.uns['log1p']["base"] = None
            adata.obs['leiden'] = None
            sc.tl.leiden(adata)
            grouplen = len(adata.obs['leiden'].value_counts())
            sc.set_figure_params(dpi_save=200, figsize=(10, 7), fontsize=10)
            if type == 'PCA':
                sc.pl.pca(adata, color='leiden', show=False, save=session_id + '.png', size=100)
                image_path = str(sc.settings.figdir) + '/pca' + session_id + '.png'
            elif type == 'UMAP':
                sc.pl.umap(adata, color='leiden', show=False, save=session_id + '.png', size=100)
                image_path = str(sc.settings.figdir) + '/umap' + session_id + '.png'
            sc.tl.rank_genes_groups(adata, groupby='leiden', method='wilcoxon', key_added='wilcoxon')
            sc.tl.dendrogram(adata, 'leiden')
            adata.write(adata_path)
            plotly_fig = px.imshow(io.imread(image_path))
            plotly_fig.update_xaxes(visible=False)
            plotly_fig.update_yaxes(visible=False)
            plotly_fig.update_traces(hovertemplate=None, hoverinfo='skip')
            plotly_fig.update_layout(width=1100, height=700, template="simple_white")
            return plotly_fig, {'display': 'none'}, {'display': 'block'}, clustertype, grouplen
    except Exception as e:
        print(f"Error in auto_scatter_plot: {e}")
        return {}, {'display': 'none'}, {'display': 'none'}, clustertype, grouplen

@app.callback(
    Output('manual-cluster-div', 'style', allow_duplicate=True),
    Output('cluster-option-div', 'style', allow_duplicate=True),
    Output('clustertype', 'data', allow_duplicate=True),
    Input('manual-cluster-button', 'n_clicks'),
    prevent_initial_call=True)
def manual_clustering(n):
    """
        Callback function for manual clustering button.

        Args:
            n (int): Number of times the button has been clicked.

        Returns:
            tuple: A tuple containing the following:

                - A dictionary containing CSS styles for the manual cluster div.
                - A dictionary containing CSS styles for the cluster option div.
                - A string representing the cluster type.
    """
    try:
        if n is None:
            return dash.no_update
        else:
            clustertype = "manual"
            return {'display': 'block'}, {'display': 'none'}, clustertype
    except Exception as e:
        print(f"Error in manual_clustering: {e}")
        # Return default values in case of an error
        return {'display': 'none'}, {'display': 'block'}, ""

@app.callback(
    Output('cluster-container', 'children'),
    Input('add-cluster', 'n_clicks'),
    State('cluster-container', 'children'),
    prevent_initial_call=True)
def display_added_cluster(n_clicks, children):
    """
        Display a ok button and input field when the "Add Cluster" button is clicked.

        Args:
            n_clicks (int): The number of times the "Add Cluster" button has been clicked.
            children (list): A list of the children elements of the cluster container.

        Returns:
            list: The updated list of children elements of the cluster container.
    """
    try:
        children.append(html.Br())
        children.append(html.Br())
        children.append(dcc.Input(id={
            'type': 'cluster-name',
            'index': n_clicks
        }, type="text"))
        children.append(html.Button("OK", id={
            'type': 'ok-button',
            'index': n_clicks
        }))
        return children
    except Exception as e:
        print(f'Error in display_added_cluster: {e}')
        return children

@app.callback(
    Output('manual-cluster-dict', 'data'),
    Input({'type': 'ok-button', 'index': ALL}, 'n_clicks'),
    State('manual-cluster-dict', 'data'),
    State('obsmdf', 'data'),
    State('expression-graph', 'selectedData'),
    State({'type': 'cluster-name', 'index': ALL}, 'value'),
    State('type', 'value'),
    background=True,
    manager=background_callback_manager,
    running=[
        (Output("manual-progress-text", "style"), {'display': 'block'}, {'display': 'none'})
    ],
    progress=Output("manual-progress-text", "children"),
    prevent_initial_call=True
)
def add_data_to_dict(set_progress,n, manual_cluster_dict, obsmdf, selectedData, clustername, type):
    """
    Updates the manual_cluster_dict with new data based on the selected genes and clusters.

    Args:
        set_progress (function): A function to update the progress status.
        n (list): The number of times each OK button was clicked.
        manual_cluster_dict (dict): A dictionary of manual clusters and associated genes.
        obsmdf (str): A string of JSON-encoded Pandas DataFrame containing the gene expression data.
        selectedData (dict): A dictionary of selected points in the expression plot.
        clustername (list): A list of strings representing the cluster names assigned by the user.
        type (str): The type of dimensionality reduction method (PCA or UMAP).

    Returns:
        dict: A dictionary of manual clusters and associated genes.
    """
    try:
        set_progress("Loading...")
        if not any(None in trigg.values() for trigg in ctx.triggered):
            clusterCategory = pd.Series(data=(str(x) for x in clustername), dtype='category')
            obsmdf = pd.read_json(obsmdf)
            if manual_cluster_dict == None:
                manual_cluster_dict = {}
            if type == 'PCA':
                xaxis = 'X_pca1'
                yaxis = 'X_pca2'
            elif type == 'UMAP':
                xaxis = 'X_umap1'
                yaxis = 'X_umap2'
            if selectedData['points']:
                for point in selectedData['points']:
                    selectedGene = obsmdf.index[(obsmdf[xaxis] == point['x']) & (obsmdf[yaxis] == point['y'])][0]
                    manual_cluster_dict[selectedGene] = clusterCategory[
                        json.loads(ctx.triggered[0]['prop_id'].replace('.n_clicks', ''))['index']]
        return manual_cluster_dict
    except Exception as e:
        print(f'Error in add_data_to_dict: {e}')
        return {}

@app.callback (Output('expression-graph', 'figure', allow_duplicate=True),
               Output('cluster-option-div', 'style', allow_duplicate=True),
               Output('visualization-div', 'style', allow_duplicate=True),
               Output('manual-cluster-div', 'style', allow_duplicate=True),
               Output('grouplen', 'data', allow_duplicate=True),
               Input('compare-button', 'n_clicks'),
               State('adata-path', 'data'),
               State('session-id', 'data'),
               State('type', 'value'),
               State('manual-cluster-dict', 'data'),
               prevent_initial_call=True)
def manual_scatter_plot(n, adata_path, session_id, type, manual_cluster_dict):
    """
    This callback creates a scatter plot of the PCA/UMAP embedding colored by the leiden clustering
    with manually annotated clusters superimposed, calculates differential expression analysis and
    creates a dendrogram. The callback returns the scatter plot, cluster option div style,
    visualization div style, manual cluster div style, and the length of unique clusters.
    """
    try:
        if n is None:
            return dash.no_update
        else:
            adata = sc.read(adata_path)
            adata.uns['log1p']["base"] = None
            adata.obs['leiden'] = None
            for key in manual_cluster_dict.keys():
                adata.obs.loc[key, 'leiden'] = manual_cluster_dict[key]
            grouplen = len(adata.obs['leiden'].value_counts())
            sc.set_figure_params(dpi_save=200, figsize=(10, 7), fontsize=10)
            if type == 'PCA':
                sc.pl.pca(adata, color='leiden', show=False, save=session_id + '.png', size=100)
                image_path = str(sc.settings.figdir) + '/pca' + session_id + '.png'
            elif type == 'UMAP':
                sc.pl.umap(adata, color='leiden', show=False, save=session_id + '.png', size=100)
                image_path = str(sc.settings.figdir) + '/umap' + session_id + '.png'
            print('Executing differential gene analysis...')
            sc.tl.rank_genes_groups(adata, groupby='leiden', method='wilcoxon', key_added='wilcoxon')
            sc.tl.dendrogram(adata, 'leiden')
            adata.write(adata_path)
            plotly_fig = px.imshow(io.imread(image_path))
            plotly_fig.update_xaxes(visible=False)
            plotly_fig.update_yaxes(visible=False)
            plotly_fig.update_traces(hovertemplate=None, hoverinfo='skip')
            plotly_fig.update_layout(width=1100, height=700, template="simple_white")
            return plotly_fig, {'display': 'none'}, {'display': 'block'},{'display': 'none'}, grouplen
    except Exception as e:
        print(f"Error in manual_scatter_plot: {e}")
        return {}, {'display': 'none'}, {'display': 'none'},{'display': 'none'}, grouplen

@app.callback(Output("download", "data"),
              Input('download-button', 'n_clicks'),
              State('adata-path', 'data'),
              State('clustertype', 'data'),
              State({'type': 'cluster-name', 'index': ALL}, 'value'),
              State('grouplen', 'data'), 
              prevent_initial_call=True)
def differential_gene_download(n,adata_path, clustertype, clustername, grouplen):
    """
    Callback function for downloading differential gene expression data as a CSV file.

    Parameters:
        n (int): number of clicks on the download button
        adata_path (str): path to AnnData object file
        clustertype (str): type of clustering (manual/auto)
        clustername (list): list of selected cluster names
        grouplen (int): number of auto-generated clusters

    Returns:
        Dash data object with a CSV file containing the differential gene expression data
    """
    try:
        if n is None:
            return dash.no_update
        else:
            adata = sc.read(adata_path)
            df = pd.DataFrame()
            if clustertype == "manual":
                for group in clustername:
                    # Fetching the differential gene expression data to a dataframe
                    if df.empty:
                        df = sc.get.rank_genes_groups_df(adata, group=str(group), key='wilcoxon')
                    else:
                        df = pd.concat([df, sc.get.rank_genes_groups_df(adata, group=str(group), key='wilcoxon')])
            elif clustertype == "auto":
                for group in range(grouplen):
                    # Fetching the differential gene expression data to a dataframe
                    if df.empty:
                        df = sc.get.rank_genes_groups_df(adata, group=str(group), key='wilcoxon')
                    else:
                        df = pd.concat([df, sc.get.rank_genes_groups_df(adata, group=str(group), key='wilcoxon')])
            return dcc.send_data_frame(df.to_csv, filename="Differential_Gene_Expression_Information.csv")
    except Exception as e:
        print(f"Error in differential_gene_download: {e}")
        return dash.no_update

@app.callback (Output('expression-graph', 'figure', allow_duplicate=True),
               Input('visualization-plot-button', 'n_clicks'),
               State('adata-path', 'data'),
               State('session-id', 'data'),
               State('plotComboBox', 'value'),
               State('n_genes', 'value'),
               State('clustertype', 'data'),
               prevent_initial_call=True)
def visualization_plot(n, adata_path, session_id, plotComboBox, n_genes, clustertype):
    """
    Generate visualization plot based on user input and return a Plotly figure.

    Parameters:
        n (int): Number of times the plot button has been clicked.
        adata_path (str): Path to the input AnnData file.
        session_id (str): Unique ID for the current session.
        plotComboBox (str): Type of plot to generate.
        n_genes (int): Number of genes to include in the plot.
        cluster_type (str): Type of clustering to use.

    Returns:
        Plotly figure object.
    """
    try:
        if n is None:
            return dash.no_update
        else:
            adata = sc.read(adata_path)
            if n_genes is None:
                n_genes = 5
            if plotComboBox == 'Dendrogram':
                sc.pl.dendrogram(adata, 'leiden', show=False, save=session_id + '.png')
                image_path = str(sc.settings.figdir) + '/dendrogram' + session_id + '.png'
            elif plotComboBox == 'Gene Ranking':
                sc.pl.rank_genes_groups(adata, n_genes=n_genes, sharey=False, key='wilcoxon', show=False,
                                        save=session_id + '.png')
                image_path = str(sc.settings.figdir) + '/rank_genes_groups_leiden' + session_id + '.png'
            elif plotComboBox == 'Dot Plot':
                sc.pl.rank_genes_groups_dotplot(adata, n_genes=n_genes, key='wilcoxon', groupby='leiden', dendrogram=False,
                                                show=False, save=session_id + '.png')
                image_path = str(sc.settings.figdir) + '/dotplot_' + session_id + '.png'
            elif plotComboBox == 'Violin':
                return {}
            elif plotComboBox == 'Stacked Violin':
                sc.pl.rank_genes_groups_stacked_violin(adata, n_genes=n_genes, key='wilcoxon', groupby='leiden', show=False,
                                                       save=session_id + '.png', dendrogram=False)
                image_path = str(sc.settings.figdir) + '/stacked_violin_' + session_id + '.png'
            elif plotComboBox == 'Matrix Plot':
                sc.pl.rank_genes_groups_matrixplot(adata, n_genes=n_genes, key='wilcoxon', groupby='leiden', show=False,
                                                   save=session_id + '.png', dendrogram=False)
                image_path = str(sc.settings.figdir) + '/matrixplot_' + session_id + '.png'
            elif plotComboBox == 'Heatmap':
                if clustertype == "auto":
                    sc.pl.rank_genes_groups_heatmap(adata, n_genes=n_genes, key='wilcoxon', groupby='leiden',
                                                    show_gene_labels=True, show=False, save=session_id + '.png',
                                                    dendrogram=False)
                elif clustertype == "manual":
                    adata_sub = adata[adata.obs['leiden'].isin(adata.obs['leiden'].dropna().sort_values().unique()), :]
                    sc.pl.rank_genes_groups_heatmap(adata_sub, n_genes=n_genes, key='wilcoxon', groupby='leiden',
                                                    show_gene_labels=True, show=False, save=session_id + '.png',
                                                    dendrogram=False)
                image_path = str(sc.settings.figdir) + '/heatmap' + session_id + '.png'
            elif plotComboBox == 'Tracksplot':
                if clustertype == "auto":
                    sc.pl.rank_genes_groups_tracksplot(adata, n_genes=n_genes, key='wilcoxon', groupby='leiden', show=False,
                                                       save=session_id + '.png', dendrogram=False)
                elif clustertype == "manual":
                    adata_sub = adata[adata.obs['leiden'].isin(adata.obs['leiden'].dropna().sort_values().unique()), :]
                    sc.pl.rank_genes_groups_tracksplot(adata_sub, n_genes=n_genes, key='wilcoxon', groupby='leiden',
                                                       show_gene_labels=True, show=False, save=session_id + '.png',
                                                       dendrogram=False)
                image_path = str(sc.settings.figdir) + '/tracksplot' + session_id + '.png'
            if image_path:
                plotly_fig = px.imshow(io.imread(image_path))
                plotly_fig.update_xaxes(visible=False)
                plotly_fig.update_yaxes(visible=False)
                plotly_fig.update_traces(hovertemplate=None, hoverinfo='skip')
                plotly_fig.update_layout(width=1100, height=700, template="simple_white")
                return plotly_fig
    except Exception as e:
        print(f"Error in visualization_plot: {e}")
        return {}

@app.callback(Output('violinplot-radio-buttons', 'children'),
              Output('violinplot-radio-buttons', 'style'),
              Output('expression-graph', 'style'),
              Input('visualization-plot-button', 'n_clicks'),
              State('plotComboBox', 'value'),
              State('adata-path', 'data'),
              State('session-id', 'data'),
              State('n_genes', 'value'),
              prevent_initial_call=True)
def violin_plot_radio(n, plotComboBox, adata_path, session_id, n_genes):
    """
    Generate radio items and display styles for violin plot visualization.

    Parameters:
        n (int): The number of clicks on the "visualization-plot-button" button.
        plotComboBox (str): The type of plot selected from the combo box.
        adata_path (str): The path to the AnnData object.
        session_id (str): The session ID.
        n_genes (int): The number of genes to use for the violin plot.

    Returns:
        tuple: A tuple of three elements containing the radio items for the violin plot, the style for the radio
        items, and the style for the expression graph.
    """
    if ctx.triggered_id == 'visualization-plot-button' and plotComboBox == 'Violin':
        adata = sc.read(adata_path)
        sc.pl.rank_genes_groups_violin(adata, n_genes=n_genes, key='wilcoxon', show=False,
                                       save=session_id + '.png')
        return [dbc.RadioItems(id='dimred-radio',
                               options=[{'label': x, 'value': x} for x in
                                        adata.obs['leiden'].dropna().sort_values().unique()],
                               value=adata.obs['leiden'].dropna().sort_values().unique()[0],
                               inline=True)], {'display': 'block'}, {'display': 'block'}
    else:
        return [], {'display': 'none'}, {'display': 'block'}


@app.callback(Output('expression-graph', 'figure', allow_duplicate=True),
              Input('dimred-radio', 'value'),
              State('session-id', 'data'),
              prevent_initial_call=True)
def violin_plot(dimred_radio, session_id):
    """
    Updates the expression-graph figure with a violin plot generated based on the selected `dimred_radio`.

    Args:
        dimred_radio (str): The selected value of the `dimred-radio` component.
        session_id (str): The id of the current session.

    Returns:
        go.Figure: A plotly Figure object.
    """
    image_path = str(sc.settings.figdir) + '/rank_genes_groups_leiden_' + str(dimred_radio) + session_id + '.png'
    plotly_fig = px.imshow(io.imread(image_path))
    plotly_fig.update_xaxes(visible=False)
    plotly_fig.update_yaxes(visible=False)
    plotly_fig.update_traces(hovertemplate=None, hoverinfo='skip')
    plotly_fig.update_layout(width=1100, height=700, template="simple_white")
    return plotly_fig












