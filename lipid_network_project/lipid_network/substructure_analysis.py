from .utils import to_json_series, NpEncoder
from .models import (
    LipidomeSummary, SubstructureAnalysis
)
from linex2.lipid import Lipid
from linex2 import plot_pca
from linex2._network_lipidome_summary import _ceiling_division_, _get_colours_
from linex2 import LipidNetwork as LINEX2Network
import plotly.figure_factory as ff
import plotly.graph_objects as go
from typing import Union
import pandas as pd
import numpy as np


SUBSTRUCTURE_FUNCTIONS = [
    Lipid.get_lipid_class, Lipid.sum_length, Lipid.sum_dbs,
    Lipid.get_headgroup
]
SUBSTRUCTURE_COMBINATIONS = 2
N_FEATURES = 20


def _add_subplots_axis_labels_(layout: dict, xlabel: str, ylablel: str):
    layout['annotations'].append(
        {
            'font': {'size': 14},
            'showarrow': False,
            'text': xlabel,
            'x': 0.5,
            'xanchor': 'center',
            'xref': 'paper',
            'y': 0,
            'yanchor': 'top',
            'yref': 'paper',
            'yshift': -30
        }
    )
    layout['annotations'].append(
        {
            'font': {'size': 14},
            'showarrow': False,
            'x': 0,
            'xshift': -50,
            'xanchor': 'top',
            'xref': 'paper',
            'y': 0.5,
            'yanchor': 'center',
            'yref': 'paper',
            'text': ylablel,
            'textangle': -90
        }
    )


def add_summary_plots(
    summary: LipidomeSummary, network: LINEX2Network
) -> Union[None, dict]:
    if not summary.data:
        summ_data = {}
        # lipidome summary plots
        class_data, len_data, db_data, c_db_scatter = \
            network.lipidome_summary(
                as_dict=True, alpha=.5,
                plot_args={"linewidth": 2, "linecolour": "rgba(0, 0, 0, 1.0)"}
            )
        # lipid class compositions
        summ_data['lipid_class_plot'] = to_json_series(class_data)
        summ_data['lipid_class_layout'] = to_json_series(
            {'xaxis': {'title': 'Lipid Class'},
             'yaxis': {'title': 'Relative Abundance'},
             'boxmode': 'group',
             'title': 'Lipid Class composition'}
        )
        # C vs. DB scatter
        c_db_scatter[1]['title'] = 'Chain Length - Double Bond Abundance'
        _add_subplots_axis_labels_(c_db_scatter[1], 'Chain Length', 'Double Bond')
        summ_data['c_db_plot'] = to_json_series(c_db_scatter[0], encoder=NpEncoder)
        summ_data['c_db_layout'] = to_json_series(c_db_scatter[1], encoder=NpEncoder)
        # chain length
        len_data[1]['title'] = 'Total Chain Lengths'
        len_data[1]['boxmode'] = 'group'
        len_data[1]['showlegend'] = True
        _add_subplots_axis_labels_(len_data[1], 'Chain Length', 'Relative Abundance')
        summ_data['sum_length_plot'] = to_json_series(len_data[0])
        summ_data['sum_length_layout'] = to_json_series(len_data[1])
        # double bonds
        db_data[1]['title'] = 'Total Number of Double Bonds'
        db_data[1]['boxmode'] = 'group'
        _add_subplots_axis_labels_(db_data[1], 'Double Bond', 'Relative Abundance')
        summ_data['sum_db_plot'] = to_json_series(db_data[0])
        summ_data['sum_db_layout'] = to_json_series(db_data[1])
        return summ_data
    return summary.data


def plotly_clustermap(data: pd.DataFrame, labels: pd.Series) -> dict:
    scaled_data = pd.DataFrame(
        ((data.values.T - np.nanmean(data, axis=1)) / np.nanstd(data, axis=1)).T,
        index=data.index, columns=data.columns
    )
    # Adapted from: https://plotly.com/python/dendrogram/
    # Initialize figure by creating upper dendrogram
    dendro_data = scaled_data.replace([np.inf, -np.inf], np.nan).fillna(0).values
    dendro_top = ff.create_dendrogram(
        dendro_data.T,
        orientation='bottom'
    )

    # Create Side Dendrogram
    dendro_side = ff.create_dendrogram(
        dendro_data,
        orientation='right'
    )

    # Create Heatmap (ordered by dendrogram)
    fig = go.Figure()
    row_leaves = dendro_side['layout']['yaxis']['ticktext']
    row_leaves = list(map(int, row_leaves))
    column_leaves = dendro_top['layout']['xaxis']['ticktext']
    column_leaves = list(map(int, column_leaves))
    fig.add_trace(
        go.Heatmap(
            x=scaled_data.columns[column_leaves],
            y=scaled_data.index[row_leaves],
            z=dendro_data[row_leaves, :][:, column_leaves],
            text=[
                dendro_data.shape[0] * [[f'Group: {label}' for label in labels[scaled_data.columns[column_leaves]]]]
            ][0]
        )
    )

    # TODO: add group colouring
    # Edit Layout
    unique_labels = labels.unique()
    colours = _get_colours_(unique_labels.size)
    label_colors = [[i/unique_labels.size, colours[i]] for i in range(unique_labels.size)]
    label_map = dict(zip(sorted(unique_labels), np.arange(unique_labels.size)))
    fig.add_trace(
        go.Heatmap(
            x=scaled_data.columns[column_leaves],
            y=scaled_data.shape[1] * [1],
            z=list(labels[scaled_data.columns[column_leaves]].map(label_map)),
            text=[f'Group: {label}' for label in labels[scaled_data.columns[column_leaves]]],
            colorscale=label_colors,
            xaxis='x2',
            yaxis='y2',
            showscale=False
        )
    )
    # Heatmap x-axis
    fig.update_layout(xaxis={'domain': [0, 1],
                             'mirror': False,
                             'showgrid': False,
                             'showline': False,
                             'zeroline': False,
                             'ticks': ""})
    # Group colour x-axis
    fig.update_layout(xaxis2={'domain': [0, 1],
                              'mirror': False,
                              'showgrid': False,
                              'showline': False,
                              'zeroline': False,
                              'showticklabels': False,
                              'ticks': ""})

    # Heatmap y-axis
    fig.update_layout(yaxis={'domain': [0, .95],
                             'mirror': False,
                             'showgrid': False,
                             'showline': False,
                             'zeroline': False,
                             'ticks': ""
                             })
    # Group colour y-axis
    fig.update_layout(yaxis2={'domain': [.95, 1],
                              'mirror': False,
                              'showgrid': False,
                              'showline': False,
                              'zeroline': False,
                              'showticklabels': False,
                              'ticks': ""})
    return fig.to_plotly_json()


def _heatmap_to_json_(figure: dict, as_string: bool = False) -> Union[dict, str]:
    for i, trace in enumerate(figure['data']):
        if trace['type'] == 'scatter':
            for arr in ['x', 'y']:
                if isinstance(trace[arr], np.ndarray):
                    if trace[arr].ndim == 1:
                        trace[arr] = list(trace[arr])
                    else:
                        # NOTE: this assumes that ndim > 2 is not occuring
                        trace[arr] = [list(trace[arr][i, :])
                                      for i in range(trace[arr].shape[0] - 1, -1, -1)]
            figure['data'][i] = trace
        if trace['type'] == 'heatmap':
            if not isinstance(trace['z'], list):
                trace['z'] = [list(trace['z'][i, :])
                              for i in range(trace['z'].shape[0] - 1, -1, -1)]
            if not isinstance(trace['x'], list):
                trace['x'] = list(trace['x'])
            if not isinstance(trace['y'], list):
                trace['y'] = list(trace['y'])
            figure['data'][i] = trace
    figure['layout'].pop('height', '')
    figure['layout'].pop('width', '')
    if as_string:
        return {
            'data': to_json_series(figure['data']),
            'layout': to_json_series(figure['layout'])
        }
    return figure


def add_substructure_analysis(
    substruct: SubstructureAnalysis, network: LINEX2Network
):
    if not substruct.data:
        plotting = {}
        substructure_data, coefficients = network.substructure_selection(
            functions=SUBSTRUCTURE_FUNCTIONS,
            feature_combinations=SUBSTRUCTURE_COMBINATIONS
        )
        if coefficients.shape[0] == 2:
            coefficients = coefficients.iloc[0, :]
        substructure_results = {'data': substructure_data,
                                'coefficients': coefficients}
        # plotting entire substructure data
        plotting['full_clustermap'] = _heatmap_to_json_(
            plotly_clustermap(substructure_data, network.groups),
            as_string=True
        )
        full_pca = plot_pca(substructure_data, network.groups,
                            to_plotly=True, title='')
        plotting['full_pca_data'] = to_json_series(full_pca[0])
        plotting['full_pca_layout'] = to_json_series(full_pca[1])
        if coefficients.ndim == 1 or coefficients.shape[0] == 1:
            # finding most discriminative features
            sorted_coeffs = np.argsort(-abs(coefficients))
            format_features = [[(coefficients.index[sorted_coeffs[i]],
                                 round(coefficients.values[sorted_coeffs[i]], ndigits=3),
                                 round(abs(coefficients.values[sorted_coeffs[i]]), ndigits=3))
                                for i in range(coefficients.size)
                                ]]
            plotting['substructure_groups'] = [coefficients.name]
            best_features = coefficients.index[sorted_coeffs][:N_FEATURES]
        else:
            format_features = {
                coefficients.index[i]: coefficients.iloc[i, :].sort_values(ascending=True)
                for i in range(coefficients.shape[0])
            }
            n_features = _ceiling_division_(N_FEATURES, len(format_features))
            best_features = []
            for coeff_data in format_features.values():
                best_features += list(coeff_data.index[:n_features])
            format_features = [
                [(format_features[group].index[i], round(format_features[group][i], ndigits=3),
                  round(abs(format_features[group][i]), ndigits=3))
                 for group in coefficients.index]
                for i in range(coefficients.shape[1])
            ]
            plotting['substructure_groups'] = list(coefficients.index)
        plotting['substructures'] = format_features
        # plotting most discriminative features
        plotting['best_clustermap'] = _heatmap_to_json_(
            plotly_clustermap(substructure_data.loc[best_features, :],
                              network.groups),
            as_string=True
        )
        best_pca = plot_pca(substructure_data.loc[best_features, :],
                            network.groups,
                            to_plotly=True, title='')
        plotting['best_pca_data'] = to_json_series(best_pca[0])
        plotting['best_pca_layout'] = to_json_series(best_pca[1])
        return substructure_results, plotting
    return substruct.data, substruct.plot
