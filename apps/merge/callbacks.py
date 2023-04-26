import os
from dash import html
from dash import dcc
from app import app
from dash import html, dcc, DiskcacheManager, CeleryManager
from dash.dependencies import Input, Output, State
import dash_uploader as du
from pathlib import Path
from zipfile import ZipFile
import scbean.model.vipcca as vip
import scbean.tools.utils as tl
import dash
import json


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

@du.callback(
    [Output('uploaded-files', 'data'),
    Output('merge-button','disabled')],
    id="merge-uploader",
)
def upload_files(status: du.UploadStatus):
    """
    Callback function to handle completion of file upload.

    Parameters:
    -----------
    status: UploadStatus
        The status of the file upload.

    Returns:
    --------
    A dictionary containing a list of the uploaded file paths.
    """
    try:
        if status.is_completed:
            uploaded_files={'uploaded_files':[]}
            for filename in status.uploaded_files:
                uploaded_files['uploaded_files'].append(str(filename))
            return uploaded_files, False
    except Exception as e:
        print(f"Error occurred while uploading file: {e}")
        return uploaded_files, True


@app.callback(
    Output('merge-file-path', 'data'),
    Output('download-merge-button','disabled'),
    Input('merge-button', 'n_clicks'),
    State('uploaded-files','data'),
    background=True,
    manager=background_callback_manager,
    running=[
        (Output("merge-progress-bar", "style"),{"visibility": "visible"},{"visibility": "hidden"},),
    ],
    progress=[Output("merge-progress-bar", "value"), Output("merge-progress-bar", "max"), Output("merge-progress-bar", "label")],
    prevent_initial_call=True
)
def merge_adata(set_progress,n, uploaded_files):
    """
    Merge multiple single-cell RNA sequencing data sets using VIPCCA algorithm.

    Parameters:
    set_progress (Tuple[int, int, str]): Tuple containing the progress bar value, maximum value, and progress message.
    n (Optional[int]): The number of times the 'merge-button' has been clicked.
    uploaded_files (dict): The uploaded files from the user.

    Returns:
    Optional[str]: The file path to the merged AnnData object in h5ad format.

    """
    if n is None:
        return dash.no_update
    else:
        adatas = []
        for filename in uploaded_files['uploaded_files']:
            if '\\' in filename:
                filename = filename.replace('\\', '/')
            foldarpath = filename.rsplit('/', 1)[0]
            if filename.endswith('h5ad'):
                adata = tl.read_sc_data(filename, batch_name=filename[filename.rfind('/')+1:filename.rfind('.')])
                adatas.append(adata)
            elif filename.endswith('.zip'):
                with ZipFile(filename, 'r') as zip:
                    name = zip.infolist()[0].filename
                    zip.extractall(foldarpath)
                    adata = tl.read_sc_data(Path(foldarpath+'/'+name), fmt='10x_mtx', batch_name=filename[filename.rfind('/')+1:filename.rfind('.')])
                    adatas.append(adata)
            elif filename.endswith('.h5'):
                # print(filename[filename.rfind('/')+1:filename.rfind('.')])
                adata = tl.read_sc_data(filename,fmt='10x_h5', batch_name=filename[filename.rfind('/')+1:filename.rfind('.')])
                adatas.append(adata)
        print(adatas)
        try:
            os.environ["TF_CPP_MIN_LOG_LEVEL"] = "2"
            set_progress(((str(2), str(10), str("Preprocessing..."))))
            print("Preprocessing...")
            adata_all = tl.preprocessing(adatas, index_unique="-")
            set_progress(((str(4), str(10), str("Generate Model..."))))
            print("Generate Model...")
            handle = vip.VIPCCA(
                adata_all=adata_all,
                res_path=foldarpath+'/data/vipcca/mixed_cell_lines/',
                mode='CVAE',
                split_by="_batch",
                epochs=20,
                lambda_regulizer=5,
                batch_input_size=128,
                batch_input_size2=16
            )
            set_progress(((str(6), str(10), str("Fit and integrate model..."))))
            print("Fit and integrate model...")
            # do integration and return an AnnData object
            adata_integrate = handle.fit_integrate()
            print("Saving merged anndata file...")
            set_progress(((str(8), str(10), str("Saving merged anndata file..."))))
            print(adata_integrate)
            adata_integrate.write(foldarpath+'/integrated_adata.h5ad')
            set_progress(((str(10), str(10), str("Completed..."))))
            return foldarpath+'/integrated_adata.h5ad', False
        except Exception as e:
            print(f"Error in merge_adata: {e}")
            return '', True

@app.callback(Output("download-merge", "data"),
              Input('download-merge-button','n_clicks'),
              State('merge-file-path', 'data'),
              prevent_initial_call=True)
def download_merged_anndata(n,merge_file_path):
    """
    Downloads the merged AnnData file when the download button is clicked.

    Parameters:
    n (int): The number of times the download button has been clicked.
    merge_file_path (str): The file path of the merged AnnData file.

    Returns:
    The merged AnnData file for download, or `dash.no_update` if `n_clicks` is `None` or `merge_file_path` is `None`.
    """
    try:
        if n is None:
            return dash.no_update
        else:
            if merge_file_path:
                return dcc.send_file(merge_file_path)
            else:
                return dash.no_update
    except Exception as e:
        print(f"Error in download_merged_anndata: {e}")
        return dash.no_update
