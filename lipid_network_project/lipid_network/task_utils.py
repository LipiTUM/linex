import os
import sys
import pandas as pd
import numpy as np
from django.conf import settings
from .models import (
    UploadedData, UploadedGroups,
    ComputationProgress, ComputedNetwork
)
from .LINEX1.exceptions import (
    PartialCorrelationError,
    SignificanceTestError,
    NameConversionError
)
from .LINEX1.Lipids.LipidSpecies import unify_name
import traceback


def save_failed_network(userid, **kwargs):
    com_net = ComputedNetwork(
        message="Network Computation Failed!",
        userid=userid, **kwargs
    )
    # com_net.message = "Network Computation Failed!"
    # com_net.userid = userid
    com_net.save()


def reset_progress(userid, **kwargs):
    ComputationProgress.objects.filter(userid=userid).delete()
    com_progress = ComputationProgress(userid=userid, **kwargs)
    com_progress.save()


def unknown_error(userid, has_groups, exception, user_warnings, unconverted, version):
    e = sys.exc_info()[0]
    print(e)
    print(exception)
    print(traceback.print_exc())
    progress = ComputationProgress.objects.get(userid=userid)
    reset_progress(
        userid=userid, error=True, groups=has_groups,
        message=f'{exception}<br><b>Unknown internal error.<b> '
                'If you cannot resolve it please contact the developers on <a '
                'href="https://github.com/lipitum/linex">GitHub</a>.',
        warning=user_warnings, unconverted=unconverted, reading=progress.reading,
        lynx=progress.lynx, compute_network=progress.compute_network,
        compute_statistics=progress.compute_statistics, class_init=progress.class_init
    )
    save_failed_network(userid, version=version)


def load_data(userid, path, has_groups, user_warnings, version):
    data_model = UploadedData.objects.get(userid=userid)
    data_file = data_model.docfile.path

    file_ending = data_file.split('.')[-1].lower()
    try:
        if file_ending == "tsv":
            data = pd.read_csv(os.path.join(path, data_file),
                               index_col=0, sep="\t")
        elif file_ending == "csv":
            data = pd.read_csv(os.path.join(path, data_file),
                               index_col=0)
        else:
            reset_progress(
                userid,
                message=f'Lipid data must be provided as a .csv or a .tsv file.\n'
                        f'You provided a {file_ending} file.',
                error=True, warning=user_warnings,
                groups=has_groups, reading="Failed",
                version=version
            )
            return
    except FileNotFoundError:
        log_file = f"{userid}_upload_error.log"
        with open(os.path.join(settings.USER_DATA, log_file)) as file:
            file.write(
                f"Expected filename: {data_file}\n\n=====\n\n"
                "Found files:\n-------------\n"
            )
            for userfile in sorted(os.listdir(settings.USER_DATA)):
                file.write(f"{userfile}\n")
        reset_progress(
            userid,
            message="Uploaded file was not found! Please try to remove "
                    "special characters from the file and filename. If this "
                    "does not resolve the error please contact the developers "
                    f"to check the log file {log_file}",
            error=True, warning=user_warnings,
            groups=has_groups, reading="Failed",
            version=version
        )
        return
    try:
        _ = unify_name(data.index.name, {})
        reset_progress(
            userid,
            message='The first column of your data file was recognized as a '
                    'lipid column based on its name.\n Please make sure that '
                    'the first column in your data file contains the sample '
                    'names and all subsequent columns represent lipid species.',
            error=True, warning=user_warnings,
            groups=has_groups, reading="Failed",
            version=version
        )
        return
    except NameConversionError:
        pass

    data.index = data.index.astype(str)
    nas = data.isna().values.any()
    if nas:
        user_warnings.append(
            f"NaNs found in data. Some statistics might not be "
            f"computed for all lipids!")

    if UploadedGroups.objects.filter(userid=userid).count() == 1:
        group_model = UploadedGroups.objects.get(userid=userid)
        group_file = group_model.group_docfile.path
    else:
        group_file = None

    return group_file, data, nas


def load_groups(userid, group_file, path, has_groups, user_warnings,
                data_samples, reference_group, version):
    file_ending = group_file.split('.')[-1].lower()
    if file_ending == "tsv":
        groups = pd.read_csv(os.path.join(path, f"uploaded_group_file_{userid}.tsv"),
                             index_col=0, sep="\t")
    elif file_ending == "csv":
        groups = pd.read_csv(os.path.join(path, f"uploaded_group_file_{userid}.csv"),
                             index_col=0)
    else:
        reset_progress(
            userid,
            message=f'Sample groups must be provided as a .csv or a .tsv file.\n'
                    f'You provided a {file_ending} file.',
            error=True, warning=user_warnings,
            groups=has_groups, reading="Failed",
            version=version
        )
        save_failed_network(userid, version=version)
        raise ValueError(
            f'Sample groups must be provided as a .csv or a .tsv file.\n'
            f'You provided a {file_ending} file.'
        )
    if groups.shape[1] == 0:
        reset_progress(
            userid,
            message=f'Sample group contains only a single column, but two are required.\n'
                    f'Please make sure the first column are sample names and the second column group annotations.',
            error=True, warning=user_warnings,
            groups=has_groups, reading="Failed",
            version=version
        )
        save_failed_network(userid, version=version)
        raise ValueError(
            f'Sample group contains only a single column, but two are required.\n'
            f'Please make sure the first column are sample names and the second column group annotations.'
        )
    else:
        groups = groups.iloc[:, 0]
    groups.index = groups.index.astype(str)
    group_nas = groups.isna()
    if group_nas.any():
        user_warnings.append(
            f"{group_nas.sum()} NaNs found in sample group annotations. "
            "This can lead to complications in computations and visualisation!"
        )

    l1_template = 'Sample(s) {0} from uploaded {1} was/were not found ' \
                  'in {2} file.\n Please ensure all sample names are correct and samples are in the first ' \
                  'column of your data file. If all names are correct '\
                  'make sure you select/de-select "Lipids are in columns" depending on your data table ' \
                  'formatting.'
    l2_template = 'Sample(s) {0} from uploaded {1} was/were not found ' \
                  'in {2} file.\n Please ensure all sample names are correct and samples are in the first ' \
                  'column of your data file. If all names are correct '\
                  'make sure samples are in rows and lipids are in columns.'

    sample_matches = np.isin(data_samples, groups.index)
    if not sample_matches.all():
        if version == 2:
            #           'formatting.'
            message = l2_template.format(list(data_samples.values[~sample_matches]), 'data', 'groups')
        else:
            message = l1_template.format(list(data_samples.values[~sample_matches]), 'data', 'groups')
        reset_progress(
            userid, message=message,
            error=True, warning=user_warnings,
            groups=has_groups, reading="Failed",
            version=version
        )
        save_failed_network(userid, version=version)
        raise ValueError(message)

    group_index_matches = np.isin(groups.index, data_samples)
    if not group_index_matches.all():
        if version == 2:
            message = l2_template.format(list(groups.index[~group_index_matches]), 'groups', 'data')
        else:
            message = l1_template.format(list(groups.index[~group_index_matches]), 'groups', 'data')
        reset_progress(
            userid, message=message,
            error=True, warning=user_warnings,
            groups=has_groups, reading="Failed",
            version=version
        )
        save_failed_network(userid, version=version)
        raise ValueError(message)

    if isinstance(reference_group, str):
        if reference_group.strip() == "":
            reference_group = None
        else:
            if not np.isin(reference_group, groups).all():
                user_warnings.append(
                    f"Reference group '{reference_group}' not found in "
                    "group data. Computations will run with no reference. "
                    "Please check for spelling errors if you wish to compute "
                    "with a reference."
                )

    group_counts = groups.value_counts()
    suff_group_sizes = np.sum(group_counts > 3)
    if suff_group_sizes < 2:
        error_message = "At least 2 samples groups with at least 3 samples " \
                        "per group are required to run LINEX. Only found " \
                        f"{suff_group_sizes} groups with at least 3 samples."
        reset_progress(
            userid=userid, message=error_message, error=True,
            warning=user_warnings, groups=has_groups, reading="Failed",
            version=version
        )
    for group in group_counts.index:
        gsum = group_counts[group]
        if gsum < 3:
            user_warnings.append(
                f"{gsum} samples of group '{group}' found, but > 2"
                f"required! '{group}' will be removed for the analysis."
            )
            groups = groups[groups.values != group]
            if group == reference_group:
                reference_group = groups.values[0]
    return groups, reference_group


def check_lynx_conversion(
    userid, unconverted, network, has_groups, user_warnings,
    version
):
    if len(unconverted) == network.data.shape[1]:
        reset_progress(
            userid, lynx="Failed",
            message='Lipids could not be converted to LipidLynxX format. '
                    'Please make sure your nomenclature is '
                    'compatible (see <a href="/tutorial"><i '
                    'class="fa fa-chalkboard-teacher fa-fw"aria-hidden="true">'
                    '</i>Tutorial</a> and <a '
                    'href="https://lipidmaps.org/lipidlynxx/">'
                    'LipidLynxX@LipidMaps</a>) '
                    'and make sure your data contains samples in rows and '
                    'lipids in columns.<br>'
                    'If you need further support with adapting lipid names, '
                    'please contact us via '
                    'email (see <a href="/about"><i '
                    'class="fa fa-info-circle fa-fw" aria-hidden='
                    '"true"></i>About</a> page)',
            error=True, warning=user_warnings, unconverted=unconverted,
            groups=has_groups, reading="Done", version=version
        )
        save_failed_network(userid, version=version)
        return
    elif len(unconverted) > 0:
        user_warnings.append(
            'Some Lipids could not be converted to LipidLynxX format. '
            'Please check the unconverted list and adapt nomenclature if '
            'these lipids should be included. For a more in-depth explanation '
            'of the LipidLynxX format and the convertable nomenclature please '
            'consult '
            '<a href="https://lipidmaps.org/lipidlynxx/">its homepage</a>. '
            'The list of unconverted as well as the uploaded data with names '
            'in the LipidLynxX format can be downloaded '
            '<a href="/download">here</a>.'
        )
        reset_progress(
            userid, lynx="Done", class_init="Done",
            compute_network="In Progress",
            groups=has_groups,
            warning=user_warnings,
            unconverted=unconverted, reading="Done",
            version=version
        )
    else:
        reset_progress(userid, lynx="Done", class_init="Done",
                       compute_network="In Progress",
                       groups=has_groups, warning=user_warnings,
                       unconverted=unconverted, reading="Done",
                       version=version)


def check_group_data_names(data, groups, user_warnings):
    remove = {}
    if data.index.duplicated().any():
        dup_names = ", ".join(list(data.index[data.index.duplicated()]))
        error = "Input data contains duplicated sample names: " \
                f"'{dup_names}'. Please remove all duplicates and make " \
                "sure that the sample groups also contain no duplicates " \
                "(if given)\n"

        return [], [], error
    if groups is not None:
        if groups.index.duplicated().any():
            dup_names = ", ".join(
                list(groups.index[groups.index.duplicated()]))
            error = "Sample groups contain duplicated names: " \
                    f"'{dup_names}'. Please remove all duplicates."
            return [], [], error
        elif data.shape[0] < groups.size:
            missing = set(groups.index).difference(data.index)
            user_warnings.append(
                "Sample groups contain samples that are missing in data file: "
                f"'{', '.join(missing)}'.\n"
            )
            remove["groups"] = list(missing)
        elif groups.size < data.shape[0]:
            missing = set(groups.index).difference(data.index)
            user_warnings.append(
                "Lipid data contain samples that are missing in group file: "
                f"'{', '.join(missing)}'. These will be removed.\n"
            )
            remove["data"] = list(data)
    return user_warnings, remove, None


def add_pcorrs(network, edge_colours, user_warnings,
               groups, estimator, significance_threshold):
    network.compute_partial_correlations(
        estimator=estimator,
        significance=significance_threshold
    )
    if groups is not None:
        try:
            network.compute_correlation_changes(partial_corrs=True)
            edge_colours.append("partial_correlations")
            edge_colours.append("partial_correlation_changes")
        except SignificanceTestError as ste:
            # fisher's z was not possible
            user_warnings.append(
                "Partial Correlation Significance could not be computed "
                "because of too small sample to feature ratio. No partial "
                "correlation changes will be computed."
            )
            print(ste)
            edge_colours.append("partial_correlations")
    else:
        edge_colours.append("partial_correlations")
    return edge_colours, user_warnings


def compute_pcors(
    userid, unconverted, network, groups, has_groups, nas,
    edge_colours, significance_threshold, user_warnings,
    version
):
    if nas:
        try:
            edge_colours, user_warnings = add_pcorrs(
                network, edge_colours, user_warnings,
                groups, estimator="empirical",
                significance_threshold=significance_threshold
            )
        except FloatingPointError:
            user_warnings.append(
                "Partial Correlations with empirical covariances failed. No partial correlations could be "
                "computed!"
            )
            reset_progress(userid, lynx="Done", class_init="Done",
                           compute_network="Done",
                           compute_statistics="In Progress",
                           groups=has_groups,
                           warning=user_warnings,
                           unconverted=unconverted, reading="Done",
                           version=version)
        except PartialCorrelationError:
            user_warnings.append(
                "Invalid values in partial correlations calculation. Skipping..."
            )
            reset_progress(userid, lynx="Done", class_init="Done",
                           compute_network="Done",
                           compute_statistics="In Progress",
                           groups=has_groups,
                           warning=user_warnings,
                           unconverted=unconverted, reading="Done",
                           version=version)
    else:
        try:
            network.compute_partial_correlations(significance=significance_threshold)
            edge_colours, user_warnings = add_pcorrs(network, edge_colours, user_warnings,
                                                     groups, estimator="GraphLasso",
                                                     significance_threshold=significance_threshold)
        except FloatingPointError:
            user_warnings.append(f"\n Partial Correlations with GraphLasso failed due to ill-conditioned solver")
            reset_progress(userid, lynx="Done", class_init="Done",
                           compute_network="Done",
                           compute_statistics="In Progress",
                           groups=has_groups,
                           warning=user_warnings,
                           unconverted=unconverted, reading="Done",
                           version=version)
            try:
                edge_colours, user_warnings = add_pcorrs(network, edge_colours, user_warnings,
                                                         groups, estimator="LedoitWolf",
                                                         significance_threshold=significance_threshold)
            except FloatingPointError:
                user_warnings.append(f"\n Partial Correlations with Ledoit-Wolf failed due to "
                                     f"ill-conditioned solver")
                reset_progress(userid, lynx="Done", class_init="Done",
                               compute_network="Done",
                               compute_statistics="In Progress",
                               groups=has_groups,
                               warning=user_warnings,
                               unconverted=unconverted, reading="Done",
                               version=version)
                try:
                    edge_colours, user_warnings = add_pcorrs(network, edge_colours, user_warnings,
                                                             groups, estimator="empirical",
                                                             significance_threshold=significance_threshold)
                except FloatingPointError:
                    user_warnings.append(
                        "Partial Correlations with empirical covariances failed. No partial correlations "
                        "could be computed!"
                    )
                    reset_progress(userid, lynx="Done", class_init="Done",
                                   compute_network="Done",
                                   compute_statistics="In Progress",
                                   groups=has_groups,
                                   warning=user_warnings,
                                   unconverted=unconverted, reading="Done",
                                   version=version)
                except PartialCorrelationError:
                    user_warnings.append(
                        "Invalid values in partial correlation calculation. Skipping..."
                    )
                    reset_progress(userid, lynx="Done", class_init="Done",
                                   compute_network="Done",
                                   compute_statistics="In Progress",
                                   groups=has_groups,
                                   warning=user_warnings,
                                   unconverted=unconverted, reading="Done",
                                   version=version)
            except PartialCorrelationError:
                user_warnings.append(
                    "Invalid values in partial correlation calculation. Skipping..."
                )
                reset_progress(userid, lynx="Done", class_init="Done",
                               compute_network="Done",
                               compute_statistics="In Progress",
                               groups=has_groups,
                               warning=user_warnings,
                               unconverted=unconverted, reading="Done",
                               version=version)
        except PartialCorrelationError:
            user_warnings.append(
                "Invalid values in partial correlation calculation. Skipping..."
            )
            reset_progress(userid, lynx="Done", class_init="Done",
                           compute_network="Done",
                           compute_statistics="In Progress",
                           groups=has_groups,
                           warning=user_warnings,
                           unconverted=unconverted, reading="Done",
                           version=version)


def compute_metrics(
    userid, network, groups, has_groups, nas,
    log_data, to_log, significance_threshold,
    test_method, compute_correlations,
    compute_partial_correlations, compute_pvalues,
    compute_fold_changes, user_warnings,
    unconverted, version
):
    node_colours = ["lipid_class", "desaturation", "chain_length",
                    "c_index", "db_index"]
    if version == 1:
        edge_colours = ["reaction_types"]
    elif version == 2:
        edge_colours = ["reaction_type"]
    else:
        raise ValueError(
            f"Illegal version number {version}. Please use only 1 and 2"
        )
    node_sizes = []
    # TODO: include edge size
    edge_sizes = []
    if compute_correlations:
        network.compute_correlations(significance=significance_threshold)
        if groups is not None:
            network.compute_correlation_changes()
            edge_colours.append("correlations")
            edge_colours.append("correlation_changes")
        else:
            edge_colours.append("correlations")
    if compute_partial_correlations:
        compute_pcors(
            userid, unconverted, network, groups, has_groups,
            nas, edge_colours, significance_threshold,
            user_warnings, version
        )
    # TODO: which errors can occur during pvalue computation?
    if compute_pvalues:
        network.compute_pvalues(method=test_method)
        node_colours.append("nlog_pvalues")
        node_sizes.append("nlog_pvalues")
    # TODO: which errors can occur during fold-change computation?
    if compute_fold_changes:
        if version == 1:
            network.compute_fold_changes(data_is_log=log_data,
                                         to_log=to_log)
        else:
            network.compute_fold_changes(data_is_log=log_data)
        node_colours.append("fold_changes")
        node_sizes.append("fold_changes")
    return node_colours, node_sizes, edge_colours, edge_sizes, network


def scores_to_plotly(scores: dict):
    # TODO: all components into one plot or one plot per component?
    return {
        component: {'y': list(c_scores), 'x': [i + 1 for i in range(len(c_scores))], 'type': 'scatter'}
        for component, c_scores in scores.items()
    }


def network_to_json_dict(network):
    return {
        'nodes': network['nodes'],
        'edges': network['edges']
    }
