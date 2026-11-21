import sqlite3
import sys
from pathlib import Path
from urllib.parse import parse_qs, quote

import dash_ag_grid as dag
from dash import Dash, dcc, html, Input, Output, no_update, ctx

import pandas as pd
import plotly.express as px


DB = sys.argv[1]
PORT = int(sys.argv[2]) if len(sys.argv) > 2 else 8050

QUERY_FILE = Path(__file__).resolve().parent.parent / "db" / "queries" / "q1.sql"


def get_runs():
    with sqlite3.connect(DB) as con:
        con.row_factory = sqlite3.Row

        with open(QUERY_FILE) as f:
            query = f.read()

        return [dict(row) for row in con.execute(query).fetchall()]

def get_samples(run_id):
    query = """
        SELECT
            rs.sample_id,
            rs.condition,
            s.phenotype,
            s.is_synthetic_phenotype
        FROM run_samples AS rs
        LEFT JOIN samples AS s
            ON rs.sample_id = s.sample_id
        WHERE rs.run_id = ?
        ORDER BY rs.sample_id
    """

    with sqlite3.connect(DB) as con:
        con.row_factory = sqlite3.Row
        return [
            dict(row)
            for row in con.execute(query, (run_id,)).fetchall()
        ]

app = Dash(__name__, suppress_callback_exceptions=True)

app.layout = html.Div([
    dcc.Location(id="url"),
    html.Div(id="page-content")
])

def get_sites(run_id):
    query = """
        SELECT
            chrom,
            pos,
            ref,
            alt,
            qual,
            filter
        FROM variants
        WHERE run_id = ?
        ORDER BY chrom, pos
    """

    with sqlite3.connect(DB) as con:
        con.row_factory = sqlite3.Row

        return [
            dict(row)
            for row in con.execute(query, (run_id,)).fetchall()
        ]
def create_filter_plot(sites):
    df = pd.DataFrame(sites)

    if df.empty:
        return px.bar(title="How many variants have each FILTER value?")

    filter_counts = (
        df["filter"]
        .fillna("Missing")
        .value_counts()
        .reset_index()
    )

    filter_counts.columns = ["filter", "count"]

    fig = px.bar(
        filter_counts,
        x="filter",
        y="count",
        title="How many variants have each FILTER value?",
        labels={
            "filter": "FILTER Value",
            "count": "Number of Variants"
        }
    )

    return fig

def classify_variant(ref, alt):
    if "," in alt:
        return "multiallelic"

    if len(ref) == 1 and len(alt) == 1:
        return "SNP"

    return "indel"

def classify_titv(ref, alt):
    transitions = {
        ("A", "G"),
        ("G", "A"),
        ("C", "T"),
        ("T", "C")
    }

    if (ref.upper(), alt.upper()) in transitions:
        return "transition"

    return "transversion"

def create_type_plot(sites):
    df = pd.DataFrame(sites)

    if df.empty:
        return px.pie(title="How many variants are SNPs, indels, or multiallelic?")

    df["variant_type"] = df.apply(
        lambda row: classify_variant(row["ref"], row["alt"]),
        axis=1
    )

    type_counts = (
        df["variant_type"]
        .value_counts()
        .reset_index()
    )

    type_counts.columns = ["variant_type", "count"]

    fig = px.pie(
        type_counts,
        names="variant_type",
        values="count",
        title="How many variants are SNPs, indels, or multiallelic?"
    )

    return fig



def create_titv_plot(sites):
    df = pd.DataFrame(sites)

    if df.empty:
        return px.bar(
            title="How many SNPs are transitions versus transversions?"
        )

    # Keep only SNPs
    df = df[
        df.apply(
            lambda row: classify_variant(row["ref"], row["alt"]) == "SNP",
            axis=1
        )
    ].copy()

    if df.empty:
        return px.bar(
            title="How many SNPs are transitions versus transversions?"
        )

    # Classify each SNP
    df["titv"] = df.apply(
        lambda row: classify_titv(row["ref"], row["alt"]),
        axis=1
    )

    # Count transitions and transversions
    titv_counts = (
        df["titv"]
        .value_counts()
        .reset_index()
    )

    titv_counts.columns = ["titv", "count"]

    fig = px.bar(
        titv_counts,
        x="titv",
        y="count",
        title="How many SNPs are transitions versus transversions?",
        labels={
            "titv": "Mutation Type",
            "count": "Number of SNPs"
        }
    )

    return fig




def runs_page():
    runs = get_runs()

    return html.Div([
        html.H1("Variant Calling Dashboard"),
        html.H2("Pipeline Runs"),

        html.H3("Select a Run"),

        html.Div([
            html.Div(
                dcc.Link(
                    run["run_id"],
                    href=f"/results?run={quote(str(run['run_id']), safe='')}"
                ),
                style={"marginBottom": "8px"}
            )
            for run in runs
        ]),

        dag.AgGrid(
            id="runs-table",
            rowData=runs,
            columnDefs=[
                {"field": "run_id", "headerName": "Run ID"},
                {"field": "implementation", "headerName": "Implementation"},
                {"field": "platform", "headerName": "Platform"},
                {"field": "n_samples", "headerName": "Samples"},
                {"field": "minutes", "headerName": "Minutes"}
            ],
            defaultColDef={
                "sortable": True,
                "filter": True,
                "resizable": True
            },
            dashGridOptions={
                "pagination": True,
                "paginationPageSize": 10
            },
            style={"height": "350px", "width": "100%"}
        )
    ], style={"padding": "30px"})




def results_page(run_id):
    samples = get_samples(run_id)
    sites = get_sites(run_id)
    filter_fig = create_filter_plot(sites)
    type_fig = create_type_plot(sites)
    titv_fig = create_titv_plot(sites)

    return html.Div([
        dcc.Link("← Back to Runs", href="/"),

        html.H1("Run Results"),
        html.H2(f"Run: {run_id}"),

        html.H3("Sample Metadata"),

        dag.AgGrid(
            id="samples-table",
            rowData=samples,
            columnDefs=[
                {"field": "sample_id", "headerName": "Sample ID"},
                {"field": "condition", "headerName": "Condition"},
                {"field": "phenotype", "headerName": "Phenotype"},
                {
                    "field": "is_synthetic_phenotype",
                    "headerName": "Synthetic Phenotype"
                }
            ],
            defaultColDef={
                "sortable": True,
                "filter": True,
                "resizable": True
            },
            columnSize="sizeToFit",
            style={"height": "300px", "width": "100%"}
        ),

        html.H3("Variant Analysis"),

        dcc.Graph(
            id="plot-filter",
            figure=filter_fig
        ),

        dcc.Graph(
            id="plot-type",
            figure=type_fig
        ),

        dcc.Graph(
            id="plot-titv",
            figure=titv_fig
        ),

        html.H3("Variant Sites"),

        dag.AgGrid(
            id="sites-table",
            rowData=sites,
            columnDefs=[
                {"field": "chrom", "headerName": "Chromosome"},
                {"field": "pos", "headerName": "Position"},
                {"field": "ref", "headerName": "Reference"},
                {"field": "alt", "headerName": "Alternate"},
                {"field": "qual", "headerName": "Quality"},
                {"field": "filter", "headerName": "Filter"}
            ],
            defaultColDef={
                "sortable": True,
                "filter": True,
                "resizable": True
            },
            dashGridOptions={
                "pagination": True,
                "paginationPageSize": 20
            },
            columnSize="sizeToFit",
            style={"height": "450px", "width": "100%"}
        )

    ], style={"padding": "30px"})


@app.callback(
    Output("page-content", "children"),
    Input("url", "pathname"),
    Input("url", "search")
)
def display_page(pathname, search):
    if pathname == "/":
        return runs_page()

    if pathname == "/results":
        params = parse_qs((search or "").lstrip("?"))
        run_id = params.get("run", [None])[0]

        if not run_id:
            return html.H2("No run selected")

        return results_page(run_id)



@app.callback(
    Output("url", "href"),
    Input("runs-table", "cellClicked"),
    prevent_initial_call=True
)
def navigate_to_results(cell):
    if cell and cell.get("colId") == "run_id":
        run_id = cell.get("value")

        if run_id is not None:
            return f"/results?run={quote(str(run_id), safe='')}"

    return no_update


@app.callback(
    Output("sites-table", "rowData"),
    Input("plot-filter", "clickData"),
    Input("plot-type", "clickData"),
    Input("plot-titv", "clickData"),
    Input("url", "search"),
    prevent_initial_call=True
)
def filter_sites(filter_click, type_click, titv_click, search):
    params = parse_qs((search or "").lstrip("?"))
    run_id = params.get("run", [None])[0]

    if not run_id:
        return []

    sites = get_sites(run_id)

    # Reset the table when navigating to another run
    if ctx.triggered_id == "url":
        return sites

    if ctx.triggered_id == "plot-filter":
        click = filter_click
    elif ctx.triggered_id == "plot-type":
        click = type_click
    elif ctx.triggered_id == "plot-titv":
        click = titv_click
    else:
        return no_update

    if not click or not click.get("points"):
        return no_update

    point = click["points"][0]

    # Pie charts use "label"; bar charts usually use "x"
    category = point.get("label", point.get("x"))

    if category is None:
        return no_update

    if ctx.triggered_id == "plot-filter":
        return [
            site for site in sites
            if (site["filter"] if site["filter"] is not None else "Missing")
            == category
        ]

    if ctx.triggered_id == "plot-type":
        return [
            site for site in sites
            if classify_variant(site["ref"], site["alt"]) == category
        ]

    if ctx.triggered_id == "plot-titv":
        return [
            site for site in sites
            if classify_variant(site["ref"], site["alt"]) == "SNP"
            and classify_titv(site["ref"], site["alt"]) == category
        ]

    return no_update


if __name__ == "__main__":
    app.run(port=PORT, debug=True)