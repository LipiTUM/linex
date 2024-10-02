from .LINEX1.Network.vis_utils import ATTR_NAMES as LINEX1_ATTR_NAMES
from linex2.vis_utils import ATTR_NAMES as LINEX2_ATTR_NAMES
from linex2.vis_utils import DynamicVisParser
from .utils import to_json_series, NpEncoder
from .models import ComputedNetwork
from django.http import HttpResponse, JsonResponse
from .utils import TIMEOUT_DELTA_MIN
import pandas as pd
import os
from typing import List
import subprocess


def get_network_view_as_dict(
    network_entry, network, version, network_type='native'
):
    mol_net = False
    if version == 2:
        if network_type != 'native':
            net_split = network_type.split('_')
            if len(net_split) == 1:
                network.compute_native_network(
                    bipartite=network_type == 'bipartite',
                    bipartite_type="hyper"
                )
            elif len(net_split) == 2:
                mol_net = True
                network.compute_molecular_network(
                    bipartite=net_split[0] == 'bipartite',
                    bipartite_type="hyper"
                )
        vis_net: DynamicVisParser = network.dynamic_network(
            node_colour_attributes=network_entry.node_colours,
            edge_colour_attributes=network_entry.edge_colours,
            node_size_attributes=network_entry.node_sizes,
            random_seed=42
        )
        vis_opts = {
            "node_colours": network_entry.node_colours,
            "node_sizes": network_entry.node_sizes,
            "edge_colours": network_entry.edge_colours
        }
        if hasattr(network, "comparisons"):
            comps = network.comparisons
        else:
            comps = None
        if hasattr(network, "unique_groups"):
            unique_groups = network.unique_groups
        else:
            unique_groups = None
        vis_context = vis_net.write_vis(
            "", vis_opts, force_hierarchical=True,
            groups=unique_groups,
            as_dict=True, group_combinations=comps,
            use_molecular_strings=mol_net
        )
        # converting to be able to loop through for selects
        vis_context["node_colours"] = {attr: LINEX2_ATTR_NAMES[attr] for attr in
                                       network_entry.node_colours}
        vis_context["edge_colours"] = {attr: LINEX2_ATTR_NAMES[attr] for attr in
                                       network_entry.edge_colours}
        vis_context["node_sizes"] = {attr: LINEX2_ATTR_NAMES[attr] for attr in
                                     network_entry.node_sizes}

    else:
        network = network
        # updating to have all attributes at download
        vis_net = network.dynamic_network(
            node_colour_attributes=network_entry.node_colours,
            edge_colour_attributes=network_entry.edge_colours,
            node_size_attributes=network_entry.node_sizes,
            directed=network_entry.directed, random_seed=42
        )
        network = network

        vis_opts = {
            "node_colours": network_entry.node_colours,
            "node_sizes": network_entry.node_sizes,
            "edge_colours": network_entry.edge_colours
        }
        if hasattr(network_entry.network, "comparisons"):
            comps = network_entry.network.comparisons
        else:
            comps = None
        vis_context = vis_net.write_vis(
            "", vis_opts, force_hierarchical=True,
            groups=network.unique_groups,
            as_dict=True, group_combinations=comps
        )
        # converting to be able to loop through for selects
        vis_context["node_colours"] = {attr: LINEX1_ATTR_NAMES[attr] for attr in
                                       network_entry.node_colours}
        vis_context["edge_colours"] = {attr: LINEX1_ATTR_NAMES[attr] for attr in
                                       network_entry.edge_colours}
        vis_context["node_sizes"] = {attr: LINEX1_ATTR_NAMES[attr] for attr in
                                     network_entry.node_sizes}

    vis_context["nodes"] = to_json_series(vis_context["nodes"], encoder=NpEncoder)
    vis_context["edges"] = to_json_series(vis_context["edges"], encoder=NpEncoder).replace(': NaN', ': "nan"')
    vis_context["legend_colour_nodes"] = to_json_series(vis_context["legend_nodes"], encoder=NpEncoder)
    vis_context["legend_colour_edges"] = to_json_series(vis_context["legend_edges"], encoder=NpEncoder)
    vis_context["legend_size_nodes"].append(
        {
            'x': 2,
            'y': 3,
            'shape': 'triangle',
            'label': 'Enzyme',
            'fixed': True,
            'id': len(vis_context['legend_size_nodes']),
            'size': 30,
            'color': "#0065bd"
        }
    )
    vis_context['legend_size_nodes'].append(
        {
            'x': 2,
            'y': 1,
            'shape': 'dot',
            'label': 'Lipid Species',
            'fixed': True,
            'id': len(vis_context['legend_size_nodes']),
            'size': 30,
            'color': "#0065bd"
        }
    )
    vis_context["legend_size_nodes"] = to_json_series(vis_context["legend_size_nodes"], encoder=NpEncoder)
    vis_context["legend_size_edges"] = to_json_series(vis_context["legend_size_edges"], encoder=NpEncoder)

    # TODO: add edge sizes
    vis_context["version"] = version
    return vis_context, network


def check_version2_network(request):
    if ComputedNetwork.objects.filter(userid=request.session["session_id"]).count() != 1:
        message = f'No network found. Either you did not compute one or your data has been' \
                  f'deleted after {TIMEOUT_DELTA_MIN} minutes.'
        return JsonResponse(
            {'error': True, 'message': message}, safe=False)
    network_entry = ComputedNetwork.objects.get(userid=request.session["session_id"])
    if not network_entry.network:
        message = f'No network found. Either you did not compute one or your data has been' \
                  f'deleted after {TIMEOUT_DELTA_MIN} minutes.'
        return JsonResponse(
            {'error': True, 'message': message}, safe=False)
    return network_entry


def build_function_call(parameters):
    return f'{parameters["group1"]}_{parameters["group2"]}_{parameters["ratio_diff"]}_{parameters["fa_penalty"]}_' \
           f'{parameters["min_size"]}_{parameters["max_size"]}_{parameters["max_iter"]}'


def enrichment_key(request_type):
    return ','.join(sorted([request_type['group1'], request_type['group2']]))


def zip_files(out_file: str, files: List[str]) -> int:
    process = subprocess.Popen(
        ["zip", "--junk-paths"] + [out_file] + files, stdout=subprocess.PIPE)
    _ = process.communicate()[0]
    return process.returncode


def bool_from_url_query(query):
    if isinstance(query, str):
        return query == 'true'
    if len(query) > 1 or len(query) == 0:
        return False
    if query[0] == 'true':
        return True
    else:
        return False


def check_get_summary_model(userid, model_type, data_attr="data"):
    if model_type.objects.filter(userid=userid).count() != 1:
        # TODO: error message
        return HttpResponse()
    entry = model_type.objects.get(userid=userid)
    if not hasattr(entry, data_attr):
        return HttpResponse()
    return entry


def save_summary_files(sum_dict, base_folder, base_path):
    files = []
    for key, summary in sum_dict.items():
        if isinstance(summary, dict):
            dataframes = []
            for group, data in summary.items():
                data['C'] = data.index
                data['Group'] = key
                data = data.melt(id_vars=["C", "Group"], var_name="DB")
                dataframes.append(data)
            summary = pd.concat(dataframes)
        file_path = os.path.join(base_folder, f"{base_path}_{key}.csv")
        summary.to_csv(file_path)
        files.append(file_path)
    return files
