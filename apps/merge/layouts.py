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
    )

menu = None
body = [
    dcc.Store(id='merge-file-path'),
    dbc.Row(dbc.Col(
        dbc.Card([
            dbc.CardHeader(html.H2("Merge Multiple 10x Genomic Files")),
            dbc.CardBody(
                dbc.Row([
                    dcc.Markdown(''' Lorem ipsum dolor sit amet, consectetur adipiscing elit. Nunc bibendum augue vel massa tincidunt, at rhoncus orci fermentum. Aliquam consequat justo eu risus dictum, ac eleifend felis vestibulum. Aenean viverra metus at est maximus, sed molestie nisi ultricies. In nec odio eu purus iaculis laoreet eu a felis. Sed ac nulla id quam aliquet sodales. Nulla facilisi. Integer vestibulum dolor in diam pulvinar luctus. Sed eu dolor malesuada, pellentesque nunc eu, lobortis sapien. Nunc eget aliquam lacus, sit amet tristique nisl. Sed dictum, est vel sagittis eleifend, arcu leo auctor elit, in vestibulum elit sapien eu arcu. Ut ornare fermentum nibh, non vestibulum nulla tincidunt a. Sed luctus ut justo in ullamcorper. Sed ut elit urna. Curabitur finibus vel mi eu iaculis.''')
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
                             html.Div(id='callback-output'),]),
                    dbc.Col([html.Button("Download merged h5ad file", id="download-merge-button"),
                             dcc.Download(id="download-merge")]),
                ])
            ),
        ], ),
    ),class_name="mb-4 mt-4"),
]
