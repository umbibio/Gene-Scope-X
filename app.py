import dash

import dash_bootstrap_components as dbc
from dash_bootstrap_templates import load_figure_template

# import data.csv_db as db
# import data.sql_db as db

url_base_pathname = '/genescopex/'
# url_base_pathname = '/'

###################
### Theme Setup ###
###################

# theme="slate"
theme="spacelab"

# This loads the themed figure template from dash-bootstrap-templates library,
# adds it to plotly.io and makes it the default figure template.
load_figure_template(theme)

app = dash.Dash(
    __name__,
    external_stylesheets=[getattr(dbc.themes, theme.upper()), dbc.icons.BOOTSTRAP],
    url_base_pathname=url_base_pathname,
    suppress_callback_exceptions=True,
    title='Single Cell RNA Sequencing', update_title='Loading...')

server = app.server