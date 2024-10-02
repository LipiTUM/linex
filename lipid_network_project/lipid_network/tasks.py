import warnings
import pickle
from .models import (
    ComputedNetwork, NetworkEnrichment, UserIds,
    ChainLengthAnalysis, SubstructureAnalysis, LipidomeSummary
)
from background_task import background
from .LINEX1.Network import LipidNetwork as LINEX1Network
from .LINEX1.exceptions import (
    LipidDataError, FaSettingError,
    LipidSettingError, NameConversionError,
    MissingClassError, SignificanceTestNAError
)
from linex2 import LipidNetwork as LINEX2Network
from linex2.exceptions import (
    MolecularSpeciesError,
    ParsingError,
    ReferenceLipidError
)
from typing import Union, Tuple
from .utils import PATH, EXAMPLE_PATH, set_running_task
from .task_utils import (
    reset_progress, save_failed_network,
    load_data, load_groups,
    check_group_data_names,
    check_lynx_conversion,
    compute_metrics,
    unknown_error,
    scores_to_plotly,
    network_to_json_dict
)
from .view_utils import get_network_view_as_dict, build_function_call, enrichment_key
from .substructure_analysis import add_summary_plots, add_substructure_analysis
from .chainlength_analysis import chainlength_analysis
import os
from os.path import exists
from operator import not_
import pandas as pd
import matplotlib
matplotlib.use('Agg')


def get_timestamp(userid):
    return UserIds.objects.get(userid=userid).timestamp


@background(schedule=0)
def compute_network_1(
    userid,
    has_groups: bool,
    lipid_level: str = None,
    lipid_settings: str = None,  # optional file for lipid classes + rules
    fa_settings: str = None,  # optional file for fatty acids + rules
    sample_axis: str = "index",  # other option: 'column'
    is_in_lynx_format: bool = True,
    compute_correlations: bool = True,
    compute_partial_correlations: bool = True,
    compute_fold_changes: bool = True,
    compute_pvalues: bool = True,
    test_method: str = "ttest",
    directed: bool = False,
    log_data: bool = False,
    to_log: bool = True,
    reference_group: str = None,
    significance_threshold: float = 0.05,
    forbidden_group_chars: bool = False
):
    # setting progress to be displayed at page
    reset_progress(userid, reading="In Progress")
    # starting actual computations
    network = None
    if forbidden_group_chars:
        user_warnings = [
            'Forbidden characters were found in the group names and removed\n'
        ]
    else:
        user_warnings = []
    unconverted = []
    try:
        group_file, data, nas = load_data(
            userid, PATH, has_groups, user_warnings, 1)
        if has_groups:
            try:
                groups, reference_group = load_groups(
                    userid, group_file, PATH, has_groups, user_warnings,
                    data_samples=getattr(data, sample_axis),
                    reference_group=reference_group,
                    version=1
                )
            except ValueError:
                set_running_task(userid, value=False)
                return
        else:
            groups = None
            compute_fold_changes = False
            compute_pvalues = False

        reset_progress(userid, lynx="In Progress", class_init="In Progress",
                       groups=has_groups, warning=user_warnings, reading="Done")
        lipid_settings = os.path.join(PATH, lipid_settings) if lipid_settings is not None else None
        fa_settings = os.path.join(PATH, fa_settings) if fa_settings is not None else None

        with warnings.catch_warnings():
            warnings.simplefilter(action='ignore', category=FutureWarning)
            lipid_network = LINEX1Network(
                lipid_data=data, class_file_path=lipid_settings,
                fa_file_path=fa_settings, names_in_convention=is_in_lynx_format,
                sample_axis=sample_axis, sample_groups=groups,
                return_unconverted=True, lipid_level=lipid_level,
                reference_group=reference_group)
        network = lipid_network

        # checks for matching and complete sample names
        user_warnings, removals, error = check_group_data_names(
            network.data, network.groups, user_warnings)
        if error is not None:
            reset_progress(
                userid, message=error, error=True, warning=user_warnings,
                groups=True, reading="Failed", version=1
            )
            set_running_task(userid, value=False)
            return
        if removals:
            if "groups" in removals.keys():
                network.groups.drop(index=removals["groups"], inplace=True)
                uni_groups = network.groups.unique
                mis_groups = set(network.unique_groups).difference(uni_groups)
                if mis_groups:
                    network.unique_groups = uni_groups
                    network.comparisons = [
                        comp for comp in network.comparisons
                        if not uni_groups.intersect(comp)
                    ]
            network.data.drop(index=removals.get("data", []), inplace=True)

        network.directed = directed
        unconverted = network.unconverted
        # checking whether naming conversion worked
        check_lynx_conversion(userid, unconverted, network, has_groups, user_warnings, version=1)
        network.compute_network()
        reset_progress(userid, lynx="Done", class_init="Done",
                       compute_network="Done",
                       compute_statistics="In Progress",
                       groups=has_groups, warning=user_warnings,
                       unconverted=network.unconverted, reading="Done")

        reset_progress(userid, lynx="Done", class_init="Done",
                       compute_network="Done",
                       compute_statistics="Done",
                       groups=has_groups,
                       warning=user_warnings,
                       unconverted=network.unconverted, reading="Done")

        node_colours, node_sizes, edge_colours, edge_sizes, network = compute_metrics(
            userid, network, groups, has_groups,
            nas, log_data, to_log, significance_threshold,
            test_method, compute_correlations,
            compute_partial_correlations, compute_pvalues,
            compute_fold_changes, user_warnings,
            network.unconverted, version=1
        )

        network_file = os.path.join(PATH, f"lipid_network_{userid}.pickle")

        com_net = ComputedNetwork(
            message="Network Computation Successful!",
            userid=userid,
            timestamp=get_timestamp(userid),
            network=network_file,
            node_colours=node_colours,
            node_sizes=node_sizes,
            edge_colours=edge_colours,
            edge_sizes=edge_sizes,
            directed=directed,
            converted=not_(is_in_lynx_format),
            version=1
        )
        # Computing view
        vis_network_file = os.path.join(PATH, f"vis_network_{userid}.pickle")

        native_context, network = get_network_view_as_dict(com_net, network, 1)
        pickle.dump(native_context, open(vis_network_file, "wb"))
        com_net.vis_network = vis_network_file
        com_net.bip_net = ""
        com_net.bip_mol_net = ""
        com_net.native_mol_net = ""
        com_net.reference_group = reference_group if reference_group is not None else ""
        com_net.is_log = log_data

        pickle.dump(network, open(network_file, "wb"))
        print(str(com_net))
        com_net.save()

        print("Done")
        reset_progress(
            userid, lynx="Done", class_init="Done",
            compute_network="Done",
            compute_statistics="Done",
            done=True, warning=user_warnings,
            groups=has_groups,
            unconverted=network.unconverted, reading="Done"
        )

    except NameConversionError as nce:
        # NOTE: this is the case when is_lynx_format is True but the name
        # is corrupted
        reset_progress(
            userid, lynx="Failed",
            message=f"Error while reading lipid species. Please check your "
                    f"input species names for LINEX/LipidLynxX "
                    f"compatibility.\n"
                    f"Row/Column at which name reading failed: '{nce.lipid}'\n"
                    f"If '{nce.lipid}' is a sample and not a lipid, transpose "
                    f"your input data so samples are in "
                    "rows and lipids in columns.\n",
            error=True, groups=has_groups, warning=user_warnings,
            unconverted=unconverted, reading="Done"
        )
        save_failed_network(userid, version=1)
    except LipidDataError as lde:
        reset_progress(
            userid, lynx="Done", class_init="Failed",
            message=f"The following error occured during data reading: '{lde}'.\n"
                    "Please make sure samples are in rows and check the tutorial for data format requirements.\n",
            error=True, groups=has_groups, warning=user_warnings,
            unconverted=unconverted, reading="Done"
        )
        save_failed_network(userid, version=1)
    except LipidSettingError as lse:
        reset_progress(
            userid, lynx="Done", class_init="Failed",
            message=f"The following error occured during loading of lipid class settings: '{lse.message}'\n"
                    " Please check the tutorial for data format requirements and make sure your input complies.",
            error=True, groups=has_groups, warning=user_warnings,
            unconverted=unconverted, reading="Done"
        )
        save_failed_network(userid, version=1)
    except FaSettingError as fse:
        reset_progress(
            userid, lynx="Done", class_init="Failed",
            message=f"The following error occured during loading of fatty acid settings: \n'{fse.message}'"
                    " Please check the tutorial for data format requirements and make sure your input complies.",
            error=True, groups=has_groups, warning=user_warnings,
            unconverted=unconverted, reading="Done"
        )
        save_failed_network(userid, version=1)
        raise fse
    except MissingClassError as me:
        reset_progress(
            userid, lynx="Done", class_init="Failed",
            message=f"Unknown lipid class '{me.class_}' found in lipid '{me.lipid}'.\n"
                    f"To fix it add '{me.class_}' to the lipid class settings file ('{me.location}' section).",
            error=True, groups=has_groups, warning=user_warnings,
            unconverted=unconverted, reading="Done"
        )
        save_failed_network(userid, version=1)
    except SignificanceTestNAError as snae:
        reset_progress(
            userid, lynx="Done", class_init="Done",
            message=snae.message + "\nPlease choose a different test or impute all missing values.",
            error=True, groups=has_groups, warning=user_warnings,
            unconverted=unconverted, reading="Done", compute_network="Done",
            compute_statistics="Failed"
        )
    except Exception as exception:
        unknown_error(userid, has_groups, exception, user_warnings, unconverted, 1)

    # Set running_task value for user back to False
    set_running_task(userid, value=False)


@background(schedule=0)
def compute_network_2(
    userid,
    has_groups: bool,
    fa_settings: str = None,  # optional file for fatty acids + rules
    lipids_are_columns: bool = True,  # other option: 'column'
    is_in_lynx_format: bool = True,
    compute_correlations: bool = True,
    compute_partial_correlations: bool = True,
    compute_fold_changes: bool = True,
    compute_pvalues: bool = True,
    test_method: str = "ttest",
    log_data: bool = False,
    to_log: bool = True,
    reference_group: str = None,
    significance_threshold: float = 0.05,
    organism: str = 'HSA',
    database: Union[str, Tuple[str, str]] = ('Rhea', 'Reactome'),
    ether_conversions: bool = True,
    forbidden_group_chars: bool = False,
    confirmed_species_file: str = None,
    use_example_data: bool = False
):
    # setting progress to be displayed at page
    reset_progress(userid, reading="In Progress", version=2)
    # starting actual computations
    if forbidden_group_chars:
        user_warnings = [
            'Forbidden characters were found in the group names and removed\n'
        ]
    else:
        user_warnings = []
    unconverted = []
    try:
        group_file = None
        if use_example_data:
            data = pd.read_csv(
                os.path.join(EXAMPLE_PATH, "example_data.csv"), index_col=0)
            nas = data.isna().values.any()
            if nas:
                user_warnings.append(
                    f"NaNs found in data. "
                    f"Some statistics might not be computed for all lipids!")
        else:
            group_file, data, nas = load_data(
                userid, PATH, has_groups, user_warnings, version=2)

        if has_groups:
            if use_example_data:
                groups = pd.read_csv(
                    os.path.join(EXAMPLE_PATH, f"example_groups.csv"),
                    index_col=0).iloc[:, 0]
                if isinstance(reference_group, str):
                    if reference_group.strip() == "":
                        reference_group = None
            else:
                try:
                    groups, reference_group = load_groups(
                        userid, group_file, PATH, has_groups, user_warnings,
                        data_samples=getattr(data, "index" if lipids_are_columns else "columns"),
                        reference_group=reference_group, version=2
                    )
                except ValueError:
                    set_running_task(userid, value=False)
                    return
        else:
            groups = None
            compute_fold_changes = False
            compute_pvalues = False

        reset_progress(userid, lynx="In Progress", class_init="In Progress",
                       groups=has_groups, warning=user_warnings, reading="Done",
                       version=2)

        if use_example_data:
            if exists(os.path.join(EXAMPLE_PATH, "example_fatty_acids.txt")):
                fa_settings = os.path.join(EXAMPLE_PATH, "example_fatty_acids.txt")
            else:
                fa_settings = ""

            if exists(os.path.join(PATH, "example_confirmed_species.txt")):
                confirmed_species_file = os.path.join(PATH, "example_confirmed_species.txt")
            else:
                confirmed_species_file = ""
        else:
            fa_settings = os.path.join(PATH, fa_settings) if fa_settings is not None else ""
            confirmed_species_file = os.path.join(PATH, confirmed_species_file) if confirmed_species_file is not None else ""
        print("starting network")
        lipid_network = LINEX2Network(
            data=data, fatty_acids=fa_settings,
            lipids_lm_compatible=is_in_lynx_format,
            groups=groups, lipids_are_columns=lipids_are_columns,
            organism=organism, database=database,
            ether_conversions=ether_conversions,
            confirmed_lipid_species=confirmed_species_file,
            allow_molspec_fails=True, reference_group=reference_group
        )
        reset_progress(
            userid, lynx="Done", class_init="Done",
            compute_network="In Progress", groups=has_groups,
            warning=user_warnings, reading="Done", version=2
        )

        # checking whether molecular species inference worked
        if lipid_network.gln.has_failed_molspec():
            user_warnings.append(
                'For some sum species no molecular species could be inferred '
                'with the given fatty acids. If you wish to include these '
                'lipids consider adding more possible fatty acids. For how '
                'to modify the default fatty acid file see the '
                '<a href="/tutorial">tutorial</a>. To download the list of '
                'removed sum species got to the '
                '<a href="/download">download page</a>.'
            )
            reset_progress(
                userid, lynx="Done", class_init="Done",
                compute_network="In Progress", groups=has_groups,
                warning=user_warnings,
                failed_molspec=lipid_network.gln.failed_molspec()
            )

        # checking whether name conversion worked
        unconverted = lipid_network.get_incompatible_lipids()
        check_lynx_conversion(
            userid, unconverted, lipid_network, has_groups, user_warnings,
            version=2)

        # computing native network
        lipid_network.compute_native_network()
        reset_progress(
            userid, lynx="Done", class_init="Done", compute_network="Done",
            compute_statistics="In Progress", groups=has_groups,
            warning=user_warnings, reading="Done", version=2)

        node_colours, node_sizes, edge_colours, edge_sizes, lipid_network = \
            compute_metrics(
                userid, lipid_network, groups, has_groups,
                nas, log_data, to_log, significance_threshold,
                test_method, compute_correlations,
                compute_partial_correlations, compute_pvalues,
                compute_fold_changes, user_warnings, unconverted,
                version=2)

        reset_progress(
            userid, lynx="Done", class_init="Done", compute_network="Done",
            compute_statistics="Done", groups=has_groups,
            warning=user_warnings, reading="Done", compute_views="In Progress",
            version=2)

        network_file = os.path.join(PATH, f"lipid_network_{userid}.pickle")

        com_net = ComputedNetwork(
            message="Network Computation Successful!",
            userid=userid,
            timestamp=get_timestamp(userid),
            network=network_file,
            node_colours=node_colours,
            node_sizes=node_sizes,
            edge_colours=edge_colours,
            edge_sizes=edge_sizes,
            converted=not_(is_in_lynx_format),
            version=2,
            reference_group=reference_group
        )
        # computing network views
        network = lipid_network
        # pre-computing all four 'modes'
        vis_network_file = os.path.join(PATH, f"vis_network_{userid}.pickle")

        # native network
        vis_network, network = get_network_view_as_dict(com_net, network, 2)
        vis_network['networks'] = {'native': network_to_json_dict(vis_network)}
        # optional networks
        # TODO: remove overhead => only store nodes and edges
        nat_mol, _ = get_network_view_as_dict(
            com_net, network, 2, 'native_mol'
        )
        vis_network['networks']['native_mol'] = network_to_json_dict(nat_mol)

        bip_net, _ = get_network_view_as_dict(
            com_net, network, 2, 'bipartite'
        )
        vis_network['networks']['bipartite'] = network_to_json_dict(bip_net)

        bip_mol, _ = get_network_view_as_dict(
            com_net, network, 2, 'bipartite_mol'
        )
        vis_network['networks']['bipartite_mol'] = network_to_json_dict(bip_mol)

        # restoring native network as default
        com_net.reference_group = "" if reference_group is None else reference_group
        com_net.is_log = log_data

        pickle.dump(network, open(network_file, "wb"))

        reset_progress(
            userid, lynx="Done", class_init="Done",
            compute_network="Done",
            compute_statistics="Done",
            compute_views="Done",
            compute_summaries="In Progress",
            done=False, warning=user_warnings,
            groups=has_groups,
            unconverted=unconverted, reading="Done",
            version=2
        )

        # chain length analysis
        # setting reference if only 2 groups and no reference given
        if network.unique_groups.size == 2 and com_net.reference_group == "":
            com_net.reference_group = network.unique_groups[0]
        vis_network['global_reference'] = com_net.reference_group

        com_net.vis_network = vis_network_file
        pickle.dump(vis_network, open(vis_network_file, "wb"))

        com_net.save()

        chainlength_model = ChainLengthAnalysis(
            userid=userid, timestamp=get_timestamp(userid),
            reference_group=com_net.reference_group,
            data={}, plot={}
        )
        with warnings.catch_warnings():
            warnings.simplefilter(action='ignore', category=FutureWarning)
            cl_data, plotting = chainlength_analysis(
                chainlength_model, network,
                com_net.reference_group, com_net.is_log
            )
        chainlength_model.data = cl_data
        chainlength_model.plot = plotting
        chainlength_model.save()

        # computing lipidome summary, substructure and chain-length analysis
        # lipidome summary
        summ = LipidomeSummary(
            userid=userid, timestamp=get_timestamp(userid),
            data={}
        )
        with warnings.catch_warnings():
            warnings.simplefilter(action='ignore', category=FutureWarning)
            summary = add_summary_plots(summ, network)
        summ.data = summary
        summ.save()

        # substructure analysis
        substruct_model = SubstructureAnalysis(
            userid=userid, timestamp=get_timestamp(userid),
            data={}
        )
        with warnings.catch_warnings():
            warnings.simplefilter(action='ignore', category=FutureWarning)
            substruct_data = add_substructure_analysis(
                substruct_model, network)
        substruct_model.data = substruct_data[0]
        substruct_model.plot = substruct_data[1]
        substruct_model.save()

        reset_progress(
            userid, lynx="Done", class_init="Done",
            compute_network="Done",
            compute_statistics="Done",
            compute_views="Done",
            compute_summaries="Done",
            done=True, warning=user_warnings,
            groups=has_groups,
            unconverted=unconverted, reading="Done",
            version=2
        )

    except ParsingError as pe:
        reset_progress(
            userid, lynx="Done", class_init="Failed",
            message=str(pe),
            error=True, groups=has_groups, warning=user_warnings,
            unconverted=unconverted, reading="Done", version=2
        )
    except ReferenceLipidError as rle:
        reset_progress(
            userid, lynx="Done", class_init="Failed",
            message=str(rle),
            error=True, groups=has_groups, warning=user_warnings,
            unconverted=unconverted, reading="Done", version=2
        )
    except MolecularSpeciesError as mse:
        reset_progress(
            userid, lynx="Done", class_init="Failed",
            message=str(mse),
            error=True, groups=has_groups, warning=user_warnings,
            unconverted=unconverted, reading="Done", version=2
        )
    except SignificanceTestNAError as snae:
        reset_progress(
            userid, lynx="Done", class_init="Done",
            message=snae.message + "\nPlease choose a different test or impute all missing values.",
            error=True, groups=has_groups, warning=user_warnings,
            unconverted=unconverted, reading="Done", compute_network="Done",
            compute_statistics="Failed", version=2
        )
    except Exception as exception:
        unknown_error(userid, has_groups, exception, user_warnings, unconverted, 2)

    # Set running_task value for user back to False
    set_running_task(userid, value=False)


@background(schedule=0)
def run_enrichment(request):
    print('Background task running')
    try:
        if NetworkEnrichment.objects.filter(userid=request["session_id"]).count() == 1:
            enrichment = NetworkEnrichment.objects.get(userid=request["session_id"])
        else:
            raise ValueError("Not enrichment entry found in DB.")
        network_entry = ComputedNetwork.objects.get(userid=request['session_id'])
        network = pickle.load(
            open(os.path.join(PATH, network_entry.network), "rb")
        )
        # formatting parameters
        difference_as_change = False if request['ratio_diff'][0] == 'ratio' else True
        request['group1'] = request['group1'][0]
        request['group2'] = request['group2'][0]
        graphs, scores, p_values, _ = network.analysis_reaction_network_enrichment(
            group1=request['group1'],
            group2=request['group2'],
            difference_as_change=difference_as_change,
            penalty_fa_reaction=float(request['fa_penalty'][0]),
            penalty_n_reactions=float(request['nreac_penalty'][0]),
            min_size=int(request['min_size'][0]),
            max_size=int(request['max_size'][0]),
            max_iter=int(request['max_iter'][0]),
            temp=float(request['temp'][0]),
            repetitions=5,
            plot_scores=False
        )

        vis_graphs = {}
        for comp, graph in graphs.items():
            nodes = [
                {'label': graph.nodes[node]['representation_id'] if graph.nodes[node]['node_molecule_type'] != 'Lipid'
                    else graph.nodes[node]['data_name'],
                 'color': graph.nodes[node]['color'],
                 'shape': 'triangle' if graph.nodes[node]['color'] == 'grey' else 'dot',
                 'id': node, 'physics': False}
                for node in graph.nodes
            ]
            edges = [
                {'from': edge[0], 'to': edge[1], 'id': i,
                 'dashes': graph.edges[edge].get('enzyme_edge', False),
                 'color': 'lightgrey'}
                for i, edge in enumerate(graph.edges)
            ]
            vis_graphs[comp] = {'nodes': nodes,
                                'edges': edges,
                                'max_score': round(max(scores[comp]), 2),
                                'p_value': p_values[comp]}
        # Write result in DB
        enrichment.computed_enrichments[enrichment_key(request)] = {
            'subnetworks': vis_graphs,
            'scores': scores_to_plotly(scores),
            'function_call': build_function_call(request),
            'groups': enrichment_key(request)
        }
        enrichment.status = 0
        enrichment.started = False
        enrichment.save()
        print('Enrichment task done')

    except Exception as e:
        enrichment.status = 1
        enrichment.started = False
        enrichment.save()
        print(e)

    # Set running_task value for user back to False
    set_running_task(request["session_id"], value=False)
