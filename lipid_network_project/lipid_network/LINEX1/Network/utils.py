from ..exceptions import (
    CorrelationError, PartialCorrelationError,
    SignificanceTestNAError
)
from typing import Union, Callable, Tuple, Set, List
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from scipy import linalg, stats
from sklearn.covariance import (
    GraphicalLasso, LedoitWolf,
    ShrunkCovariance, empirical_covariance
)
from statsmodels.stats.moment_helpers import cov2corr
from statsmodels.stats.multitest import multipletests


def show_class_connections(class_transformations: dict,
                           figsize: tuple = (16, 9),
                           savepath: str = None,
                           return_graph: bool = False,
                           ax: plt.axis = None,
                           show: bool = True,
                           **kwargs):
    graph = nx.Graph()
    graph.add_nodes_from(class_transformations.keys())
    for key, vals in class_transformations.items():
        for val in vals:
            graph.add_edge(key, val)
    node_labels = dict(zip(class_transformations.keys(),
                           class_transformations.keys()))
    pos = nx.spring_layout(graph, k=.4)

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    plt.axis("off")

    nx.draw_networkx_nodes(
        graph, pos, ax=ax,
        node_color="white",
        edgecolors="k",
        node_size=1200
    )
    nx.draw_networkx_edges(graph, pos, ax=ax)
    nx.draw_networkx_labels(graph, pos, node_labels,
                            ax=ax)

    plt.tight_layout()

    if savepath is not None:
        if not kwargs:
            dpi = 400
        else:
            dpi = kwargs.pop("dpi", 400)
        plt.savefig(savepath, dpi=dpi, **kwargs)
        plt.close()
    elif show:
        plt.show()

    if return_graph:
        return graph
    else:
        return ax


def partial_correlations(data: pd.DataFrame,
                         estimator: str = "GraphLasso",
                         **kwargs) -> Tuple[pd.DataFrame, np.ndarray]:
    """
    Computing all feature partial correlations. Implementation based on
    Whittaker (1990) and the R package pppcor

    Parameters
    ----------
    :param data: pd.DataFrame, samples in rows
    :param estimator: str, 'empirical', 'GraphLasso', 'Shrinkage' or 'LedoitWolf', covariance estimator
    :return: tuple of pd.Series with partial correlations and pvalues, indices are data.columns
    """
    # computing covariance matrix
    if estimator == "empirical":
        try:
            cov = empirical_covariance(data, **kwargs)
        # TODO: what other errors can occur?
        except ValueError:
            raise PartialCorrelationError(
                f"Invalid values for empirical covariance calculation!",
                estimator
            )
    elif estimator == "GraphLasso":
        # TODO: should we use GraphicalLassoCV here?
        try:
            graph_lasso = GraphicalLasso()
            graph_lasso.fit(data)
            cov = graph_lasso.covariance_
        except ValueError:
            raise PartialCorrelationError(
                f"Invalid values for GraphLasso partial correlation calculations!",
                estimator
            )

    elif estimator == "LedoitWolf":
        try:
            lw = LedoitWolf()
            lw.fit(data)
            cov = lw.covariance_
        except ValueError:
            raise PartialCorrelationError(
                f"Invalid values for Ledoit-Wolf partial correlation calculations!",
                estimator
            )
    elif estimator == "GraphLasso":
        try:
            shrunk = ShrunkCovariance()
            shrunk.fit(data)
            cov = shrunk.covariance_
        except ValueError:
            raise PartialCorrelationError(
                f"Invalid values for GraphLasso partial correlation calculations!",
                estimator
            )
    else:
        raise ValueError(
            f"Estimator must be one of ['GraphLasso', 'Shrinkage', 'LedoitWolf', 'empirical'], not {estimator}"
        )
    # inverse covariance matrix
    try:
        inv_cov = linalg.inv(cov)
    except ValueError:
        raise PartialCorrelationError(
            f"Invalid values in covariance matrix inversion",
            estimator
        )
    # partial correlation matrix
    pcor = -cov2corr(inv_cov)
    np.fill_diagonal(pcor, 1)
    # p-values from fisher's z-transform
    n, m = data.shape
    fact = n - (m - 2) - 3
    pvals = -np.ones(pcor.shape)
    if fact > 0:
        stats_ = np.sqrt(fact) * .5 * np.log((1 + pcor)/(1 - pcor))
        # TODO: double check if normality actually is true
        pvals = 2 * stats.norm.cdf(1 - stats_)
    else:
        print("p-values could not be calculated due to too few samples")

    pcor_df = pd.DataFrame(pcor, index=data.columns,
                           columns=data.columns)
    return pcor_df, pvals


def correlations(data: pd.DataFrame,
                 method: str = "pearson",
                 **kwargs) -> Tuple[pd.DataFrame, np.ndarray]:
    methods = {
        "pearson": stats.pearsonr,
        "spearman": stats.spearmanr,
        "kendall": stats.kendalltau
    }

    if method not in methods.keys():
        raise ValueError(
            f"'method' must be one of {list(methods.keys())}"
        )

    if data.shape[1] < 2:
        cors = pd.DataFrame(columns=data.index,
                            index=data.index)
        pvals = pd.DataFrame(columns=data.index,
                             index=data.index)
        return cors, pvals.values

    cors = np.zeros(2*[data.shape[0]])
    pvals = np.zeros(2*[data.shape[0]])
    for i in range(data.shape[0] - 1):
        for j in range(i + 1, data.shape[0]):
            try:
                cor_, pval_ = methods[method](data.values[i, :],
                                              data.values[j, :],
                                              **kwargs)
            except ValueError:
                if data.iloc[i, :].isna().values.any() or data.iloc[j, :].isna().values.any():
                    cor_ = np.nan
                    pval_ = np.nan
                else:
                    lipids = (data.index[i], data.index[j])
                    raise CorrelationError(
                        f"Invalid values at lipids {lipids} during correlation calculations.\n",
                        lipids
                    )
            cors[i, j] = cor_
            cors[j, i] = cor_
            pvals[i, j] = pval_
            pvals[j, i] = pval_

    cor_df = pd.DataFrame(cors, columns=data.index,
                          index=data.index)
    return cor_df, pvals


def _matrix_pval_correction_(data: Union[pd.DataFrame, np.ndarray],
                             **kwargs) -> np.ndarray:
    """
    Multiple test correction on a symmetric(!) p-value matrix

    Parameters
    ----------
    data
    kwargs

    Returns
    -------

    """
    lidx = np.tril_indices_from(data, -1)

    if isinstance(data, pd.DataFrame):
        pval_corr = multipletests(data.values[lidx], **kwargs)
    else:
        pval_corr = multipletests(data[lidx], **kwargs)

    corrected = np.zeros(data.shape)
    # corrected[uidx] = pval_corr[1]
    corrected[lidx] = pval_corr[1]
    corrected += corrected.T

    return corrected


def generalised_log(data, c=1e-5,
                    log_fun=np.log2):
    return log_fun(data + np.sqrt(data ** 2 + c))


def fold_changes(data: pd.DataFrame,
                 groups: Union[np.ndarray, pd.Series],
                 compare_groups: Union[Tuple[str], List[str]],
                 data_is_log: bool = True,
                 log_func=generalised_log,
                 to_log: bool = False) -> pd.Series:
    """
    Helper function for computing (log-) fold changes of means of samples groups.

    Parameters
    ----------
    data : data from which to compute fold changes. Samples must be in columns
    groups : array of sample groups in the order of data.columns
    compare_groups : list of strings of the two groups to compare
    data_is_log : boolean, indicating whether input data is already log-transformed
                  if False, data will be in-line log-transformed
    log_func : function for log transformation
    to_log: boolean, whether to do log transformation if data_is_log is False

    Returns
    -------
    pd.Series of size data.shape[0] containing fold changes
    """
    group1 = data.loc[:, groups == compare_groups[0]]
    group2 = data.loc[:, groups == compare_groups[1]]
    if data_is_log:
        return group1.mean(axis=1) - group2.mean(axis=1)
    elif to_log:
        return log_func(group1 + 1).mean(axis=1) - log_func(group2 + 1).mean(axis=1)
    else:
        m1 = group1.mean(axis=1).values
        m2 = group2.mean(axis=1).values
        fcs = pd.Series(np.zeros(data.shape[0]),
                        index=data.index)
        for i in range(data.shape[0]):
            if m1[i] >= m2[i]:
                fcs[i] = m1[i]/m2[i]
            else:
                fcs[i] = -(m2[i] / m1[i])
        return fcs


def binary_test(data: pd.DataFrame,
                groups: Union[np.ndarray, pd.Series],
                compare_groups: Union[List[str], Set[str], Tuple[str]],
                method: Union[str, Callable] = "ttest",
                p_adjust_method: str = "fdr_bh",
                **kwargs) -> pd.Series:
    methods = {
        "wilcoxon": stats.wilcoxon,
        "ranksums": stats.ranksums,
        "ttest": stats.ttest_ind,
        "mannwhitneyu": stats.mannwhitneyu
    }
    if isinstance(method, str):
        test = methods[method]
    else:
        test = method
    # subsetting groups
    if np.unique(groups).size == 2:
        f_groups = groups
    else:
        f_groups = groups[np.isin(groups, compare_groups)]
    # subsetting data to only contain groups of interest
    sub_data = data.loc[:, f_groups.index]
    # computing feature-wise p-values
    try:
        pvals = sub_data.apply(lambda x: test(x.values[f_groups == compare_groups[0]],
                                              x.values[f_groups == compare_groups[1]],
                                              **kwargs)[1],
                               axis=1)
    except ValueError:
        if any(np.isnan(sub_data.values)):
            test_names = {
               "wilcoxon": "Wilcoxon signed-rank ",
               "ranksums": "Wilcoxon rank-sum ",
               "ttest": "t-",
               "mannwhitneyu": "Mann-Whitney U "
            }
            raise SignificanceTestNAError(
                "Input data contains empty cells, which is not allowed"
                f"with {test_names[method]}test"
            )
    # NOTE: nans cause multipletest to return nans only!
    pvals[pvals.isna()] = 1
    return pd.Series(multipletests(pvals.values, method=p_adjust_method)[1],
                     index=pvals.index)


def unique_elements(data: Union[pd.DataFrame, pd.Series]) -> np.ndarray:
    """
    Helper to extract all unique elements from a data frame or a series.<br>

    For a pd.Series the result is equivalent to pd.Series.unique().values

    Parameters
    ----------
    data : a general value

    Returns
    -------

    """
    nas = data.isna()
    na_free = data.values[np.invert(nas)]
    return np.unique(na_free)


def correlation_change(c_s_: pd.DataFrame,
                       c_t_: pd.DataFrame) -> pd.DataFrame:

    if not c_s_.shape == c_t_.shape:
        raise ValueError(
            "Both input matrices need to have the same input dimensions!"
        )
    c_s = c_s_.fillna(0, inplace=False)
    c_t = c_t_.fillna(0, inplace=False)
    x_ = np.tril_indices_from(c_s.values)
    x = c_s.values[x_]
    y = c_t.loc[c_s.index, c_s.columns].values[x_]

    inter = np.array(x.size * [""], dtype=object)
    # unsignificant
    inter[(x == 0) & (y == 0)] = "unsignificant"
    # significant vs. unsignificant
    inter[((x != 0) & (y == 0))] = "significant to unsignificant"
    inter[((x == 0) & (y != 0))] = "unsignificant to significant"
    # unchanged significant
    # NOTE: needs to come before the next two or else will overwrite
    inter[(x != 0) & (y != 0)] = "unchanged significant"
    # both significant
    inter[(x > 0) & (y < 0)] = "positive to negative"
    inter[(x < 0) & (y > 0)] = "negative to positive"

    ret = np.array(c_s.shape[1]*[c_s.shape[0]*[""]],
                   dtype=object)
    ret[x_] = inter
    ret.T[x_] = inter
    return pd.DataFrame(ret, columns=c_s.columns,
                        index=c_s.index)


def _get_sizes_report_nas_(attr: pd.Series, keys):
    sizes = {}
    na_keys = []
    for node in keys:
        val = abs(attr.loc[node])
        sizes[node] = val
        if np.isnan(val):
            na_keys.append(node)
    return sizes, na_keys
