import dash_uploader as du
from pathlib import Path
import scanpy as sc
import dash
import pandas as pd
import anndata as anndata
from dash.dependencies import Input, Output, State
import dash_bootstrap_components as dbc
from app import app
from dash import dcc, ctx
from dash import html, ALL
import plotly.express as px
import plotly.graph_objects as go
from skimage import io
import json
import matplotlib.pyplot as plt
from zipfile import ZipFile
import matplotlib
matplotlib.pyplot.switch_backend('Agg')
UPLOAD_FOLDER_ROOT = r"Uploads"
du.configure_upload(app, UPLOAD_FOLDER_ROOT)
adata = anndata.AnnData(
    X=None,
    obs=None,
    var=None)
genes = pd.DataFrame()
obsmdf = pd.DataFrame()
dfadata = pd.DataFrame()
grouplen = 0
selectedArray = []
clustertype = ""

@du.callback(
    [Output('preprocessing-div', 'style'),
    Output('dash-uploader','disabled'),
    Output('dash-uploader','disabledMessage')],
    id='dash-uploader',
)
def get_adata(filename):
    print(filename)
    global adata
    filename = filename[0]
    sc.settings.figdir = filename.rsplit('/', 1)[0] + '/figures/'
    if filename.endswith('h5ad'):
        adata = sc.read(filename)
    elif filename.endswith('zip'):
        with ZipFile(filename, 'r') as zip:
            zip.extractall(filename.rsplit('/', 1)[0])
        adata = sc.read_10x_mtx(
            Path(filename.replace('.zip','')),
            var_names='gene_symbols',
            cache=False)
    elif filename.endswith('h5'):
        adata = sc.read_10x_h5(filename)
    return {'display': 'block'}, True, 'Uploaded: '+filename.rsplit('/', 1)[1]

@app.callback(
    Output('adata-div','style'),
    Output('gene', 'options'),
    Output('gene', 'value'),
    Output('scatter-plot', 'style'),
    Input('preprocessing-button','n_clicks'),
    State('svd_solver','value'),
    State('n_neighbors','value'),
    State('scatter-plot-div', 'children'))
def preprocess_adata(n, svd_solver, n_neighbors, children):
    if n is None:
        return dash.no_update
    else:
        global genes, obsmdf, dfadata
        print("Preprocessing annData file...")
        adata.var_names_make_unique()
        print("Normalize annData file...")
        sc.pp.normalize_total(adata, target_sum=1e4)
        print("Logarithmize annData file...")
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
        print("Perform principle component analysis...")
        sc.tl.pca(adata, svd_solver=svd_solver)
        print("Compute neighbors...")
        sc.pp.neighbors(adata, n_neighbors=n_neighbors, knn=True, n_pcs=40)
        print("Uniform Manifold Approximation and Projection for Dimension Reduction...")
        sc.tl.umap(adata)
        adata.uns['log1p']["base"] = None
        obsmdf = adata.obsm.to_df()
        dfadata = adata.to_df()
        genes = pd.DataFrame(adata.var.gene_ids)
        gene_ids = genes['gene_ids'].to_numpy()
        # children.append(dcc.Dropdown(id='gene',options=gene_ids,value=gene_ids[0]))
    return {'display': 'none'}, gene_ids, gene_ids[0],{'display': 'block'}

@app.callback(
    Output('2d3d-radio', 'style'),
    Input('type', 'value'),)
def update_2d3d_visibility(type):
    if type == 'PCA':
        style = {'display': 'block'}
    else:
        style = {'display': 'none'}
    return style

@app.callback([Output('expression-graph', 'figure'),
              Output('cluster-option-div', 'style'),
              Output('visualization-div', 'style'),
              Output('violinplot-radio-buttons', 'style'),],
              [Input('plot-button', 'n_clicks_timestamp'),
              Input('auto-cluster-button','n_clicks_timestamp'),
              Input('manual-cluster-button','n_clicks_timestamp'),
              Input('compare-button','n_clicks_timestamp'),
              Input('visualization-plot-button','n_clicks_timestamp'),
              State('type','value'),
              State('gene', 'value'),
              State('plotComboBox', 'value'),
              State('n_genes', 'value'),
              State('2d3d-radio', 'value')])
def scatter_plot(scatter_plot, autocluster_plot,manualcluster_plot, comparebutton, visualization_plot,  type, gene, plotComboBox, n_genes, pca_radio):
    if ctx.triggered_id == None:
        return dash.no_update, dash.no_update, dash.no_update, dash.no_update
    global xaxis, yaxis, clustertype, grouplen
    listedTimestamps = [scatter_plot, autocluster_plot, manualcluster_plot, comparebutton, visualization_plot]
    listedTimestamps = [0 if v is None else v for v in listedTimestamps]
    sortedTimestamps = sorted(listedTimestamps)
    pickedButton=""
    if scatter_plot == sortedTimestamps[-1]:
        pickedButton = "scatter_plot"
    if autocluster_plot == sortedTimestamps[-1]:
        pickedButton = "autocluster_plot"
    if manualcluster_plot == sortedTimestamps[-1]:
        pickedButton = "manualcluster_plot"
    if comparebutton == sortedTimestamps[-1]:
        pickedButton = "comparebutton"
    if visualization_plot == sortedTimestamps[-1]:
        pickedButton = "visualization_plot"
    if pickedButton == "scatter_plot" or pickedButton == "manualcluster_plot":
        adata.obs['leiden'] = None
        genename = genes.index[genes['gene_ids'] == gene].tolist()[0]
        if type == 'PCA' and pca_radio == '2D':
            fig = px.scatter(x = obsmdf.X_pca1, y = obsmdf.X_pca2,labels={'x':'PC1', 'y':'PC2'}, template="simple_white", title=genename, color=dfadata[genename])
            xaxis = 'X_pca1'
            yaxis = 'X_pca2'
        elif type == 'PCA' and pca_radio == '3D':
            fig = px.scatter_3d(x=obsmdf.X_pca1, y=obsmdf.X_pca2,z=obsmdf.X_pca3, labels={'x': 'PC1', 'y': 'PC2', 'z': 'PC3'}, template="simple_white",
                             title=genename, color=dfadata[genename])
            fig.update_traces(marker_size=2)
            fig.update_layout(scene=dict(xaxis=dict(showgrid=True),
                                         yaxis=dict(showgrid=True),
                                         zaxis=dict(showgrid=True)
                                         ))
            xaxis = 'X_pca1'
            yaxis = 'X_pca2'
        elif type == 'UMAP':
            fig = px.scatter(x = obsmdf.X_umap1, y = obsmdf.X_umap2,labels={'x':'UMAP1', 'y':'UMAP2'}, template="simple_white", title=genename, color=dfadata[genename])
            xaxis = 'X_umap1'
            yaxis = 'X_umap2'
        figure = go.Figure(data = fig)
        if pickedButton == "manualcluster_plot":
            clustertype = "manual"
            if type == 'PCA' and pca_radio == '3D':
                fig = px.scatter(x=obsmdf.X_pca1, y=obsmdf.X_pca2, labels={'x': 'PC1', 'y': 'PC2'},
                                 template="simple_white", title=genename, color=dfadata[genename])
                figure = go.Figure(data=fig)
                return figure, {'display': 'none'}, {'display': 'none'}, {'display': 'none'}
            else:
                return figure, {'display': 'none'}, {'display': 'none'}, {'display': 'none'}
        return figure, {'display': 'block'}, {'display': 'none'}, {'display': 'none'}
    elif pickedButton == "autocluster_plot":
        clustertype = "auto"
        adata.obs['leiden'] = None
        sc.tl.leiden(adata)
        grouplen = len(adata.obs['leiden'].value_counts())
        sc.set_figure_params(dpi_save=200 ,figsize=(10,7), fontsize=10)
        if type == 'PCA':
            sc.pl.pca(adata, color='leiden', show=False, save='.png', size=100)
            image_path = str(sc.settings.figdir)+'/pca.png'
        elif type == 'UMAP':
            sc.pl.umap(adata, color='leiden', show=False, save='.png', size=100)
            image_path = str(sc.settings.figdir)+'/umap.png'
        sc.tl.rank_genes_groups(adata, groupby='leiden', method='wilcoxon', key_added='wilcoxon')
        sc.tl.dendrogram(adata, 'leiden')
        plotly_fig = px.imshow(io.imread(image_path))
        plotly_fig.update_xaxes(visible=False)
        plotly_fig.update_yaxes(visible=False)
        plotly_fig.update_traces(hovertemplate = None,hoverinfo = 'skip')
        plotly_fig.update_layout(width=1100,height=700, template="simple_white")
        return plotly_fig, {'display': 'none'}, {'display': 'block'}, {'display': 'none'}
    elif pickedButton == "comparebutton":
        grouplen = len(adata.obs['leiden'].value_counts())
        sc.set_figure_params(dpi_save=200, figsize=(10, 7), fontsize=10)
        if type == 'PCA':
            sc.pl.pca(adata, color='leiden', show=False, save='.png', size=100)
            image_path = str(sc.settings.figdir) + '/pca.png'
        elif type == 'UMAP':
            sc.pl.umap(adata, color='leiden', show=False, save='.png', size=100)
            image_path = str(sc.settings.figdir) + '/umap.png'
        sc.tl.rank_genes_groups(adata, groupby='leiden', method='wilcoxon', key_added='wilcoxon')
        sc.tl.dendrogram(adata, 'leiden')
        plotly_fig = px.imshow(io.imread(image_path))
        plotly_fig.update_xaxes(visible=False)
        plotly_fig.update_yaxes(visible=False)
        plotly_fig.update_traces(hovertemplate=None, hoverinfo='skip')
        plotly_fig.update_layout(width=1100, height=700, template="simple_white")
        return plotly_fig, {'display': 'none'}, {'display': 'block'}, {'display': 'none'}
    elif pickedButton == "visualization_plot":
        if n_genes is None:
            n_genes = 5
        if plotComboBox == 'Dendrogram':
            sc.pl.dendrogram(adata, 'leiden', show=False, save='.png')
            image_path = str(sc.settings.figdir)+'/dendrogram.png'
        elif plotComboBox == 'Gene Ranking':
            sc.pl.rank_genes_groups(adata, n_genes=n_genes, sharey=False, key='wilcoxon', show=False, save='.png')
            image_path = str(sc.settings.figdir)+'/rank_genes_groups_leiden.png'
        elif plotComboBox == 'Dot Plot':
            sc.pl.rank_genes_groups_dotplot(adata, n_genes=n_genes, key='wilcoxon', groupby='leiden', dendrogram=False,
                                            show=False, save='.png')
            image_path = str(sc.settings.figdir)+'/dotplot_.png'
        elif plotComboBox == 'Violin':
            return {}, {'display': 'none'}, {'display': 'block'}, {'display': 'block'}
        elif plotComboBox == 'Stacked Violin':
            sc.pl.rank_genes_groups_stacked_violin(adata, n_genes=n_genes, key='wilcoxon', groupby='leiden', show=False,
                                                   save='.png', dendrogram=False)
            image_path = str(sc.settings.figdir)+'/stacked_violin_.png'
        elif plotComboBox == 'Matrix Plot':
            sc.pl.rank_genes_groups_matrixplot(adata, n_genes=n_genes, key='wilcoxon', groupby='leiden', show=False,
                                               save='.png', dendrogram=False)
            image_path = str(sc.settings.figdir)+'/matrixplot_.png'
        elif plotComboBox == 'Heatmap':
            if clustertype == "auto":
                sc.pl.rank_genes_groups_heatmap(adata, n_genes=n_genes, key='wilcoxon', groupby='leiden',
                                            show_gene_labels=True, show=False, save='.png', dendrogram=False)
            elif clustertype == "manual":
                adata_sub = adata[adata.obs['leiden'].isin(adata.obs['leiden'].dropna().sort_values().unique()),:]
                sc.pl.rank_genes_groups_heatmap(adata_sub, n_genes=n_genes, key='wilcoxon', groupby='leiden',
                                                show_gene_labels=True, show=False, save='.png', dendrogram=False)
            image_path = str(sc.settings.figdir)+'/heatmap.png'
        elif plotComboBox == 'Tracksplot':
            if clustertype == "auto":
                sc.pl.rank_genes_groups_tracksplot(adata, n_genes=n_genes, key='wilcoxon', groupby='leiden', show=False,
                                               save='.png', dendrogram=False)
            elif clustertype == "manual":
                adata_sub = adata[adata.obs['leiden'].isin(adata.obs['leiden'].dropna().sort_values().unique()), :]
                sc.pl.rank_genes_groups_tracksplot(adata_sub, n_genes=n_genes, key='wilcoxon', groupby='leiden',
                                                show_gene_labels=True, show=False, save='.png', dendrogram=False)
            image_path = str(sc.settings.figdir)+'/tracksplot.png'
        plotly_fig = px.imshow(io.imread(image_path))
        plotly_fig.update_xaxes(visible=False)
        plotly_fig.update_yaxes(visible=False)
        plotly_fig.update_traces(hovertemplate=None, hoverinfo='skip')
        plotly_fig.update_layout(width=1100, height=700, template="simple_white")
        return plotly_fig, {'display': 'none'}, {'display': 'block'}, {'display': 'none'}

@app.callback(Output('violinplot-radio-buttons', 'children'),
              Output('violin-graph', 'style'),
              Output('expression-graph', 'style'),
              Input('plotComboBox', 'value'),
              Input('visualization-plot-button','n_clicks'),
              State('n_genes', 'value'))
def violin_plot(plotComboBox, visualization_plot_button, n_genes):
    if ctx.triggered_id == 'visualization-plot-button' and plotComboBox == 'Violin':
        if n_genes is None:
            n_genes = 5
        sc.pl.rank_genes_groups_violin(adata, n_genes=n_genes, key='wilcoxon', show=False, save='.png')
        return [dbc.RadioItems(id='dimred-radio',
            options=[{'label': x, 'value': x} for x in adata.obs['leiden'].dropna().sort_values().unique()],
            value=adata.obs['leiden'].dropna().sort_values().unique()[0],
            inline=True)], {'display': 'block'}, {'display': 'none'}
    else:
        return [], {'display': 'none'}, {'display': 'block'}

@app.callback(Output('violin-graph', 'figure'),
              Input('dimred-radio', 'value'))
def violin_plot(dimred_radio):
    image_path = str(sc.settings.figdir)+'/rank_genes_groups_leiden_'+str(dimred_radio)+'.png'
    plotly_fig = px.imshow(io.imread(image_path))
    plotly_fig.update_xaxes(visible=False)
    plotly_fig.update_yaxes(visible=False)
    plotly_fig.update_traces(hovertemplate=None, hoverinfo='skip')
    plotly_fig.update_layout(width=1100, height=700, template="simple_white")
    return plotly_fig

@app.callback(Output("download", "data"),
              Input('download-button','n_clicks'),
              State({'type': 'cluster-name','index': ALL}, 'value'))
def differential_gene_download(n, clustername):
    if n is None:
        return dash.no_update
    else:
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

@app.callback(
    Output('manual-cluster-div','style'),
    Input('manual-cluster-button','n_clicks'),
    Input('compare-button','n_clicks'),)
def manual_clustering(n, compare):
    if ctx.triggered_id == 'compare-button':
        return {'display': 'none'}
    elif ctx.triggered_id == 'manual-cluster-button':
        return {'display': 'block'}
    else:
        return dash.no_update

@app.callback(
    Output('cluster-container', 'children'),
    Input('add-cluster', 'n_clicks'),
    State('cluster-container', 'children'))
def display_dropdowns(n_clicks, children):
    children.append(html.Br())
    children.append(html.Br())
    children.append(dcc.Input(id={
            'type': 'cluster-name',
            'index': n_clicks
        }, type="text"))
    children.append( html.Button("OK", id={
            'type': 'ok-button',
            'index': n_clicks
        }))
    return children

@app.callback(
    Output('dropdown-container-output', 'children'),
    State('expression-graph', 'selectedData'),
    Input({'type': 'ok-button','index': ALL}, 'n_clicks'),
    State({'type': 'cluster-name','index': ALL}, 'value')
)
def display_output(selectedData,n_clicks, clustername):
    if not any(None in trigg.values() for trigg in ctx.triggered):
        clusterCategory = pd.Series(data=(str(x) for x in clustername), dtype='category')
        if selectedData['points']:
            for point in selectedData['points']:
                selectedGene = obsmdf.index[(obsmdf[xaxis] == point['x']) & (obsmdf[yaxis] == point['y'])][0]
                adata.obs.loc[selectedGene, 'leiden'] = clusterCategory[json.loads(ctx.triggered[0]['prop_id'].replace('.n_clicks',''))['index']]
    return dash.no_update
