from .utils import to_json_series, NpEncoder
from .models import ChainLengthAnalysis, ComputedNetwork
from linex2 import LipidNetwork
from linex2._network_lipidome_summary import (
    _subplot_layout_, _get_colours_,
    DEFAULT_COLOUR
)
from typing import Tuple, Dict, List
import pandas as pd
import numpy as np
import re


def _lineplot_(data: pd.DataFrame, plot_list: List[dict], naxis: int,
               groups: bool = False, **kwargs):
    hovertemplate = '<b>Chain Length:</b> %{x}<br><b>Fold Change: %{y:.3f}</b>'
    if not groups:
        if not np.all(data['FC'].isna()):
            if 'name' in kwargs.keys():
                hovertemplate += f'<br><b>Lipid Class:</b> {kwargs["name"]}'
            data = data.sort_values(['variable'])
            plot_list.append({
                'x': list(data['variable']),
                'y': list(data['FC'].replace([np.inf, -np.inf], np.nan).fillna(0)),
                'type': 'scatter',
                'mode': 'lines+markers',
                'marker': kwargs.pop('marker', {'color': DEFAULT_COLOUR}),
                'hovertemplate': hovertemplate,
                **kwargs
            })
    else:
        unique_groups = data['groups'].unique()
        colours = _get_colours_(unique_groups.size)
        for i, group in enumerate(sorted(unique_groups)):
            tmp_data = data.loc[data['groups'] == group, :]
            tmp_data = tmp_data.sort_values(['variable'])
            if not np.all(tmp_data['FC'].isna()):
                hovertemplate_group = f"{hovertemplate}<br></b>Group:</b> {group}"
                plot_list.append({
                    'x': list(tmp_data['variable']),
                    'y': list(tmp_data['FC'].replace([np.inf, -np.inf], np.nan).fillna(0)),
                    'name': group,
                    'legendgroup': group,
                    'type': 'scatter',
                    'mode': 'lines+markers',
                    'marker': {
                        'color': colours[i]
                    },
                    'xaxis': f'x{naxis + 1}',
                    'yaxis': f'y{naxis + 1}',
                    'hovertemplate': hovertemplate_group,
                    **kwargs
                })


def _data_to_plotly_(data: Dict[str, pd.DataFrame], multi_group: bool) -> Tuple[Dict[str, pd.DataFrame], list, dict]:
    plots = []
    if multi_group:
        for i, (lclass, class_data) in enumerate(sorted(data.items())):
            if i == 0:
                _lineplot_(class_data, plots, i, title=lclass, showlegend=True,
                           groups=multi_group)
            else:
                _lineplot_(class_data, plots, i, showlegend=False, groups=multi_group)
        layout = _subplot_layout_(data, xaxis={'title': 'Chain Length'},
                                  yaxis={'title': 'Fold Change'})
    else:
        colours = _get_colours_(len(data))
        for i, (lclass, class_data) in enumerate(sorted(data.items())):
            if i == 0:
                _lineplot_(class_data, plots, i, title=lclass, showlegend=True,
                           groups=multi_group, name=lclass, legendgroup=lclass,
                           marker={'color': colours[i]})
            else:
                _lineplot_(class_data, plots, i, showlegend=True, groups=multi_group,
                           name=lclass, legendgroup=lclass,
                           marker={'color': colours[i]})
        layout = {
            'xaxis': {'title': 'Chain Length'},
            'yaxis': {'title': 'Fold Change'}
        }
    return data, plots, layout


def _get_chain_length_data_(
    network: LipidNetwork, reference_group: str,
    multi_group: bool = False, data_is_log: bool = False,
    **kwargs
) -> Tuple[Dict[str, pd.DataFrame], list, dict]:
    cl_data = network.analysis_chain_length(reference_group, return_data=True,
                                            data_is_log=data_is_log, **kwargs)
    if not data_is_log:
        # this is to have the data < 0 on the same scale
        # in order to not underestimate 'negative' fold changes
        for lclass, data in cl_data.items():
            mask = data['FC'].values < 1
            data['FC'].values[mask] = -(1 / data['FC'].values[mask])
    return _data_to_plotly_(cl_data, multi_group=multi_group)


def chainlength_analysis(
    chainlength_model: ChainLengthAnalysis, network: LipidNetwork,
    reference_group: str, is_log: bool, **kwargs
):
    if not chainlength_model.data:
        if reference_group == "":
            cl_data = {}
            cl_analysis = {}
            seen_groups = set()
            groups = set(network.unique_groups)
            for group in sorted(groups):
                seen_groups.add(group)
                if len(seen_groups) == len(groups):
                    break
                data, plot, layout = _get_chain_length_data_(
                    network, group,
                    other_groups=groups.difference(seen_groups),
                    data_is_log=is_log, multi_group=True, **kwargs
                )
                cl_data[group] = data
                cl_analysis[group] = {
                    'plot': to_json_series(plot, encoder=NpEncoder),
                    'layout': to_json_series(layout),
                    'name': re.sub('[^0-9a-zA-Z_]', '', group)
                }
            return cl_data, cl_analysis
        else:
            data, plot, layout = _get_chain_length_data_(
                network, reference_group,
                multi_group=False, data_is_log=is_log,
                **kwargs
            )
            cl_analysis = {
                'plot': to_json_series(plot),
                'layout': to_json_series(layout)
            }
            return data, cl_analysis
    return chainlength_model.data, chainlength_model.plot
