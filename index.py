from importlib import import_module

from dash import dcc
from dash import html
from dash.dependencies import Input, Output, State, MATCH

from app import app, server, url_base_pathname
import apps

import dash_bootstrap_components as dbc


pages = [
    {'title': 'Home', 'app': 'home', 'href': f'{url_base_pathname}'},
    {'title': 'Expression', 'app': 'expression', 'href': f'{url_base_pathname}expression'},
    {'title': 'Download', 'app': 'download', 'href': f'{url_base_pathname}download'},
]


for page in pages:
    import_module(f'apps.{page["app"]}')


nav = dbc.Nav([
    dbc.NavItem(dbc.NavLink(page['title'], href=page['href'], id={'type': 'navlink', 'page': page['app']}))
    for page in pages ], pills=True)

logo = html.Img(src=app.get_asset_url('logo.png'), className="img-fluid")

app.layout = dbc.Container(
    dbc.Row([
        dcc.Location(id='url', refresh=False),
        dbc.Col([
            dbc.Row(dbc.Col(logo)),
            dbc.Row(dbc.Col(id='left-menu')),
        ], width=2),
        dbc.Col([
            dbc.Row(dbc.Col(nav), class_name="p-2"),
            dbc.Row(dbc.Col(id='page-content')),
        ], width=10),
    ])
)


@app.callback(
    Output({'type': 'navlink', 'page': MATCH}, 'active'),
    Input('url', 'pathname'),
    State({'type': 'navlink', 'page': MATCH}, 'href'),)
def update_active_menu(pathname, page_href):
    return pathname == page_href


@app.callback(Output('left-menu', 'children'),
              Output('page-content', 'children'),
              Input('url', 'pathname'))
def display_page(pathname):
    print(f"display page: {pathname}")

    for page in pages:
        if pathname == page['href']:

            page_app = getattr(apps, page["app"])
            menu = getattr(page_app, 'menu')
            body = getattr(page_app, 'body')

            return menu, body

    print(f"page not found: {pathname}")
    return None, None


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('--host', type=str, default='127.0.0.1')
    parser.add_argument('--port', type=int, default=8050)
    parser.add_argument('--debug', action='store_true')
    
    kvargs = vars(parser.parse_args())

    app.run_server(**kvargs)

