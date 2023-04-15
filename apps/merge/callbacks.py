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
import tensorflow as tf


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
    output=Output('merge-file-path', 'data'),
    id="merge-uploader",
)
def callback_on_completion(status: du.UploadStatus):
    adatas = []
    for filename in status.uploaded_files:
        filename = str(filename)
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
            print(filename[filename.rfind('/')+1:filename.rfind('.')])
            adata = tl.read_sc_data(filename,fmt='10x_h5', batch_name=filename[filename.rfind('/')+1:filename.rfind('.')])
            adatas.append(adata)
    print(adatas)
    os.environ["TF_CPP_MIN_LOG_LEVEL"] = "2"
    adata_all = tl.preprocessing(adatas, index_unique="-")
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
    # do integration and return an AnnData object
    adata_integrate = handle.fit_integrate()
    print(adata_integrate)
    adata_integrate.write(foldarpath+'/integrated_adata.h5ad')
    print(adata_integrate.obs)
    return foldarpath+'/integrated_adata.h5ad'



@app.callback(Output("download-merge", "data"),
              Input('download-merge-button','n_clicks'),
              State('merge-file-path', 'data'),)
def anndata_download(n,merge_file_path):
    if n is None:
        return dash.no_update
    else:
        return dcc.send_file(merge_file_path)


# class MyCallback(tf.keras.callbacks.Callback):
#     def on_epoch_begin(epoch):
#         print(epoch)