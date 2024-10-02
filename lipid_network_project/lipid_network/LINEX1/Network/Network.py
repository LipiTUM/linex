from __future__ import annotations
from ..Lipids.LipidSpecies import (
    LipidSpecies,
    unify_name
)
from ..default_globals import (
    read_class_data,
    read_fa_data
)
from ..exceptions import (
    NotComputedError, LipidDataError,
    FaSettingError, LipidSettingError,
    NameConversionError, SignificanceTestError
)
from .utils import (
    show_class_connections,
    partial_correlations,
    correlations,
    _matrix_pval_correction_,
    fold_changes,
    binary_test,
    unique_elements,
    correlation_change,
    _get_sizes_report_nas_
)
from .colour_helper import _generate_colormap_
from .misc import (
    _aggregate_species_,
    _pandas_log_,
    _size_legend_,
    _infer_fatty_acids_,
    _range_scale_,
    _get_enzyme_,
    _lipid_to_object_,
    _tuple_to_string_
)
from .vis_utils import (
    VisParser, DynamicVisParser,
    STATIC_NODE_PROPERTIES, STATIC_EDGE_PROPERTIES
)
from warnings import warn
from typing import Union, Iterable, Callable, List, Tuple, Dict, Any
from itertools import combinations, permutations
import os
import re
import pandas as pd
import numpy as np
import networkx as nx
import networkx.drawing.layout as nx_layout
from pyvis.network import Network
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colorbar import Colorbar
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize, to_hex, to_rgba
import json

# NOTE: for some reason LipidLynxX changes the
# working directory in ALL of these three cases:
# 1. module import
# 2. Converter object initialisation
# 3. Converter function calls
wdir = os.getcwd()
from lynx import Converter

os.chdir(wdir)

# TODO: import lipid class colours here
LIPID_COLOURS = json.load(
    open(os.path.join(os.path.dirname(__file__), "templates/lipid_colours.json"),
         "rb")
)


class LipidNetwork:
    """
    Computation and visualisation of rule-based lipid networks.
    """
    _attribute_group_size_ = {
        # edge attributes
        "correlations": 1,
        "partial_correlations": 1,
        "reaction_types": 0,
        "reaction_enzymes": 0,
        "correlation_changes": 2,
        "partial_correlation_changes": 2,
        # node attributes
        "fold_changes": 2,
        "pvalues": 2,
        "nlog_pvalues": 2,
        "lipid_class": 0,
        "chain_length": 0,
        "desaturation": 0,
        "c_index": 0,
        "db_index": 0,
        "degree": 0,
        "closeness": 0,
        "betweenness": 0
    }
    _discrete_map_ = {
        # node features
        "lipid_class": True,
        "c_index": False,
        "db_index": False,
        "chain_length": False,
        "desaturation": False,
        "fold_changes": False,
        "pvalues": False,
        "nlog_pvalues": False,
        "betweenness": False,
        "closeness": False,
        "degree": False,
        # edge features
        "correlations": False,
        "partial_correlations": False,
        # FIXME: dynamic True/Fale
        "correlation_changes": True,
        "partial_correlation_changes": True,
        "reaction_types": True
    }
    # defining the slots yields no improvement, but is nice
    # for organisation IMO
    __slots__ = ("sample_axis", "lipid_axis",
                 "class_connections", "class_fas",
                 "class_enzymes", "lipid_level",
                 "fas", "fa_combs", "lipid_names",
                 "lipids", "n_lipids", "data",
                 "groups", "unique_groups",
                 "comparisons",
                 "network", "interaction_attributes",
                 "lipid_attributes", "layout",
                 "edge_colours", "node_colours",
                 "fa_restrictions", "reactions",
                 "directed", "unconverted",
                 "_name_map_")

    def __init__(self, class_file_path: str = None,
                 fa_file_path: str = None,
                 lipid_species: pd.Series = None,
                 lipid_data: pd.DataFrame = None,
                 sample_groups: pd.Series = None,
                 sample_axis: str = "columns",
                 lipid_level: str = "sum",
                 infer_fatty_acids: bool = False,
                 inference_level: str = "all",
                 unification_args: dict = None,
                 names_in_convention: bool = False,
                 sum_to_mol: bool = False,
                 aggregate_duplicate_species: bool = True,
                 return_unconverted: bool = False,
                 reference_group: str = None):
        """

        Parameters
        ----------
        class_file_path: str, optional
            Path to the file containing lipid class connections and
            number of fatty acids per lipid class.

        fa_file_path: str, optional
            Path to the file containing fatty acids, reaction rules
            and forbidden FA-reactions

        lipid_species: pd.Series, optional
            Lipids to be included in the network. At least one
            lipid_species and lipid_data must be provided. If both
            are given lipid_data will be extended to contain the
            union of both (wrt. to lipid_species) and values of
            species absent in lipid_data will be set to 0.

        lipid_data: pd.DataFrame, optional
            Lipid data to use for computing the network and subsequently
            measures. Required if lipid_species is not given.
            NOTE: no data-processing is implemented. If you intend to
            visualise p-values, fold-changes or similar please make sure
            your data is already processed properly.

        sample_groups: pd.Series, optional
            Vector specifying the group of each sample. Names of samples
            should be the series indices and match with the sample names
            in lipid_data.

        sample_axis: str, optional, default 'column'
            Specifies the orientation of samples in lipid_data.

        lipid_level: str, optional, default 'sum'
            'sum', 'molecular' or 'sn'.
            Specifies the 'resolution' of lipid species. Sum species
            will sum up the indices of all fatty acids in a lipid
            whereas MolecularSpecies retains individual fatty acid
            resolution without sn-positions.
            snSpecific will follow.

        infer_fatty_acids: bool, optional, default True
            If True the list of fatty acids specified in fa_file_path
            will be extended by all fatty acids that are present in
            lipid species with one fatty acid and, based on those
            fatty acids, must be present in lipids with two fatty
            acids (i.e. their sum-composition minus the inferred FAs).
            In the case of molecular species all FAs present in any
            lipid are considered.
            Currently NOT recommended as it cannot consider class-specifity
            and adds ALL inferred fatty acids to the list of un-specified FAs.

        inference_level: str, optional, default 'all'
            Which way to infer fatty acids (in addition to the user-defined
            fatty acids).
            Options are:
                * direct: use only fatty acids from species with one fatty acid
                * all: use all

        unification_args: dict, optional, default None
            Keyword arguments to pass to lynx.Converter.convert if
            lipid names are to be converted. Should not contain
            the lipid level, as this is fixed based on the class
            used and will throw an error.

        names_in_convention: bool, optional, default False
            Whether lipid names have to be converted to LipidLynxX
            format. If names are already in the correct format it
            is strongly recommended to set this parameter to True,
            as the name conversion takes up a lot more time than
            the actual network initialisation and computation.

        sum_to_mol: bool, optional, default False
            Whether to infer all molecular species possible from
            sum species composition based on possible fatty acids.
            Not yet supported => no mapping to data from sum species
            implemented so far

        aggregate_duplicate_species: bool, optional, default False
            Whether duplicate sum-species arising from different
            molecular species should be aggregated (i.e. summed).
            Currently this parameter has no effect!


        Raises
        ------
        ValueError
            If lipid_species and lipid_data are both not specified.
            If lipid_species or lipid_data are an incorrect types.

        NameConversionError
            If any of the provided lipid names is not in a format
            that can be handled by LipidLynxX


        Warns
        -----
        If some names in sample_groups are missing in lipid_data
        """
        # NOTE: supplying lipid_species and lipid_data is not adviced => slower!
        print("getting lipid classes utilities")
        #   Code used before custom reactions   #
        # class_connections, class_fas = define_global_classes(class_connect,
        #                                                      add_to_connections,
        #                                                      nFAs,
        #                                                      add_to_nFA)
        # self.class_connections = class_connections
        # self.class_fas = class_fas
        if class_file_path is None:
            class_file_path = os.path.join(os.path.dirname(__file__),
                                           "..", "lipid_classes.txt")
        try:
            class_data = read_class_data(class_file_path)
        except Exception as e:
            if isinstance(e, KeyboardInterrupt):
                raise e
            if isinstance(e, LipidSettingError):
                raise e
            raise LipidSettingError(
                "Lipid Setting file seems to be corrupted. failed at reading",
                "Unknown"
            )
        self.class_connections = class_data["class_connections"]
        self.class_fas = class_data["class_nFA"]
        self.class_enzymes = class_data["enzymes"]

        axis_invert = {"columns": "index",
                       "index": "columns"}
        if axis_invert.get(sample_axis) is None:
            raise ValueError("'sample_axis' must be 'columns' or 'index' "
                             "not {0}".format(sample_axis))
        self.sample_axis = sample_axis
        self.lipid_axis = axis_invert[sample_axis]

        self.network = None
        # converted lipids as instances of a subclass of LipidSpecies
        self.lipids = []
        # converted lipid names
        self.lipid_names = []
        self._name_map_ = {}
        # attributes of type pd.Series
        # * index/key => lipid
        # * value => attributes
        # e.g. logFCs, p-values etc. for each lipid
        self.lipid_attributes = {}

        # TODO: add checks for non-duplicate lipids
        print("loading data")
        if lipid_species is not None:
            self.unconverted = self._set_lipid_names_(lipid_species, lipid_level,
                                                      aggregate_duplicate_species,
                                                      names_in_convention,
                                                      return_unconverted)
            self.data = None
        if lipid_data is not None:
            # NOTE: lipid_data internally has samples in rows
            # we don't do any data pre-processing!
            if sample_axis == "columns":
                self.data = lipid_data.T
            else:
                self.data = lipid_data
            if aggregate_duplicate_species:
                try:
                    self.data = _aggregate_species_(self.data,
                                                    "columns")
                except Exception as e:
                    if isinstance(e, KeyboardInterrupt):
                        raise e
                    raise LipidDataError(
                        "Lipid data seems corrupted. Failed at aggregating sum species."
                    )
            if not self.lipids:
                self.unconverted = self._set_lipid_names_(lipid_data, lipid_level,
                                                          aggregate_duplicate_species,
                                                          names_in_convention,
                                                          return_unconverted)
                if self.unconverted:
                    self.data = self.data.loc[:, ~np.isin(self.data.columns, self.unconverted)]
                self.data.columns = self.lipid_names
                # NOTE: this step is necessary, because LipidLynxX can produce the same
                # output when e.g FAs are not recognised
                if self.data.columns.duplicated().any():
                    try:
                        self.data = _aggregate_species_(self.data,
                                                        "columns")
                        self.lipid_names = self.data.columns
                    except Exception as e:
                        if isinstance(e, KeyboardInterrupt):
                            raise e
                        raise LipidDataError(
                            "Lipid data seems corrupted. Failed at aggregating duplicates from conversion."
                        )
            else:
                if not names_in_convention:
                    # making sure names in data are unified the same way as in lipids
                    # TODO: think about whether this degree of flexibility is necessary
                    # => flexibility in terms of being able to specify lipids AND lipid_data
                    self._unify_lipid_data_indices(unification_args)

        if not self.lipids:
            raise ValueError(
                "either 'lipid_species' or 'lipid_data must be given"
            )
        self.lipids = [self.lipids[self._name_map_[lipid_name]]
                       for lipid_name in self.lipid_names]
        self.n_lipids = len(self.lipid_names)

        # defining the set of (possibly) present fatty acids if sum composition level
        if infer_fatty_acids:
            warn(
                """
                Fatty Acid inference is not recommended yet and should be
                use with caution! All fatty acids found will be added to the
                list of FAs used for classes which have no specific FA-pool.
                """
            )
            # sum level => need to infer FAs using Lysos, CEs etc.
            inf_lipids = _infer_fatty_acids_(self.lipids,
                                             inference_level)
            # fas, fa_combs = define_global_fas(inferred_fatty_acids=inf_lipids,
            #                                   pre_fatty_acids=fatty_acids,
            #                                   add_to_FAs=add_to_FAs)
        else:
            inf_lipids = []
            # fas, fa_combs = define_global_fas(pre_fatty_acids=fatty_acids,
            #                                   add_to_FAs=add_to_FAs)

        if fa_file_path is None:
            fa_file_path = os.path.join(os.path.dirname(__file__),
                                        "..", "fatty_acids.txt")
        try:
            fas = read_fa_data(fa_file_path)
        except Exception as e:
            if isinstance(e, KeyboardInterrupt):
                raise e
            if isinstance(e, FaSettingError):
                raise e
            raise FaSettingError(
                "Fatty Acid file seems to be corrupted. Failed at reading.",
                "Corrupted File"
            )
        self.fas = fas["fatty_acids"]
        # NOTE: inferred lipids are added to the "general" tab
        # FIXME: this HAS to come out in the function/parameter description
        self.fas.dict["general"] += inf_lipids
        self.fa_restrictions = fas["excluded_reactions"]
        self.reactions = fas["reactions"]
        # TODO: add fa combinations
        self.fa_combs = {}

        if isinstance(sample_groups, list) or isinstance(sample_groups, np.ndarray):
            self.groups = pd.Series(sample_groups)
        elif isinstance(sample_groups, pd.Series):
            self._set_data_groups_(groups=sample_groups)
            if self.data is not None:
                if not all(np.isin(sample_groups.index, self.data.index)):
                    warn(
                        "Not all sample names in 'sample_groups' were found in 'data'"
                    )
                self.groups = sample_groups[self.data.index]
        else:
            self.groups = None

        if sample_groups is not None:
            self.unique_groups = np.unique(sample_groups)
            if reference_group is not None:
                if reference_group not in self.unique_groups:
                    warn(f"{reference_group} is not one of the given groups ({list(self.unique_groups)})")
                    self.comparisons = list(combinations(self.unique_groups, 2))
                else:
                    self.comparisons = [(reference_group, group) for group in self.unique_groups
                                        if group != reference_group]
            else:
                self.comparisons = list(combinations(self.unique_groups, 2))
        else:
            self.unique_groups = np.array([])
            self.comparisons = None

        self.layout = None
        # TODO: find a good way of storing interactions
        # TODO: integrate attributes into networks?
        # maybe dict of shape
        # * key = property name
        # * value = pd.DataFrame of lipids x attributes
        # nest for different groups - first level == group name
        # e.g. (partial) correlations between lipids
        self.interaction_attributes = {
            "correlations": None,
            "partial_correlations": None,
            "correlation_changes": None,
            "partial_correlation_changes": None,
            "reaction_types": None
        }

        self.edge_colours = {}
        self.node_colours = {}
        # whether to draw a directed graph, default true
        self.directed = True

    # Alternative constructors from files instead of already
    # loaded pandas objects
    @classmethod
    def from_lipid_list(cls, lipid_species: str,
                        network_args: dict = None,
                        **kwargs) -> LipidNetwork:
        """
        Convenience function to generate a LipidNetwork instance from
        a list of lipids stored locally in a txt or csv file.

        Parameters
        ----------
        lipid_species: str
            Path to a file containing the lipids to be included in the network.

        network_args: dict, optional
            keyword arguments to pass to LipidNetwork.__init__

        kwargs:
            keyword arguments to pass to pd.read_csv

        Returns
        -------
        instance of LipidNetwork with the given lipids
        """
        species = pd.read_csv(lipid_species, **kwargs)
        if network_args is None:
            return cls(lipid_species=species)
        return cls(lipid_species=species, **network_args)

    @classmethod
    def from_lipid_data(cls, lipid_data: str,
                        sample_axis: str = "columns",
                        network_args: dict = None,
                        is_excel: bool = False,
                        **kwargs) -> LipidNetwork:
        """
        Convenience function to generate a LipidNetwork instance from
        a lipid dataset stored locally in a csv or txt file

        Parameters
        ----------
        lipid_data: str
            Path to a file containing the lipid data to be used in the network

        sample_axis: str, optional, default 'columns'
            Specifies the orientation of samples in lipid_data.

        is_excel: bool, optional, default True
            Whether to load data from a csv/txt file or an excel sheet

        network_args: dict, optional
            keyword arguments to pass to LipidNetwork.__init__

        kwargs:
            keyword arguments to pass to pd.read_csv/pd.read_excel

        Returns
        -------
        instance of LipidNetwork with the given lipid data
        """
        if is_excel:
            data = pd.read_excel(lipid_data, **kwargs)
        else:
            data = pd.read_csv(lipid_data, **kwargs)
        if network_args is None:
            return cls(lipid_data=data,
                       sample_axis=sample_axis)
        else:
            return cls(lipid_data=data,
                       sample_axis=sample_axis,
                       **network_args)

    def _set_lipid_names_(self, lipids: Union[pd.DataFrame, pd.Series, list],
                          lipid_level: str, remove_duplicates: bool,
                          names_in_convention: bool,
                          return_unconverted: bool):
        if isinstance(lipids, pd.Series) or isinstance(lipids, list):
            lipids_ = lipids
        elif isinstance(lipids, pd.DataFrame):
            lipids_ = getattr(lipids, self.lipid_axis)
        else:
            raise ValueError(
                "Lipid species must be provided as a list, pd.Series or in "
                f"a pd.DataFrame not {type(lipids).__name__}"
            )
        # converting lipids outside of LipidSpecies classes to avoid
        # mutliple initialisations of lynx.Converter
        # if not names_in_convention:
        #     lipids_ = unify_name(list(lipids_), self.class_fas)
        if remove_duplicates:
            lipids_ = _aggregate_species_(pd.Series(lipids_))
        # converting strings of lipid names to LipidSpecies objects
        # and saving some utilities for easier computation later on
        self.lipid_attributes["lipid_class"] = pd.Series()
        self.lipid_attributes["c_index"] = pd.Series()
        self.lipid_attributes["db_index"] = pd.Series()
        self.lipid_attributes["chain_length"] = pd.Series()
        self.lipid_attributes["desaturation"] = pd.Series()
        converter = Converter()
        if return_unconverted:
            unconverted = []
            for lipid in lipids_:
                try:
                    self.add_lipid(
                        lipid, lipid_level,
                        name_in_convention=names_in_convention,
                        converter=converter
                    )
                except NameConversionError:
                    unconverted.append(lipid)
                except ValueError as ve:
                    print(f"ValueError! {ve}")
                    raise NameConversionError(
                        f"Lipid name {lipid} seems corrupted",
                        lipid
                    )
            os.chdir(wdir)
            return unconverted
        else:
            for lipid in lipids_:
                self.add_lipid(
                    lipid, lipid_level,
                    name_in_convention=names_in_convention,
                    converter=converter
                )
            os.chdir(wdir)
            return None

    def _unify_lipid_data_indices(self, unification_args):
        # converting lipids outside of LipidSpecies classes to avoid
        # mutliple initialisations of lynx.Converter
        if unification_args is None:
            self.data.columns = unify_name(list(self.data.columns),
                                           self.class_fas)
        else:
            self.data.columns = unify_name(list(self.data.columns),
                                           self.class_fas,
                                           **unification_args)
        os.chdir(wdir)

    def _set_data_groups_(self, data: Union[pd.DataFrame, str] = None,
                          groups: Union[pd.Series, str] = None,
                          names_in_convention: bool = False,
                          unification_args: dict = None):
        if data is not None:
            if not hasattr(self, 'data'):
                self.data = None
            if self.data is not None:
                warn("'data' has already been set, skipping...")
            elif isinstance(data, pd.Series):
                self.data = data
            elif isinstance(data, str):
                self.data = pd.read_csv(data)
            else:
                raise ValueError(
                    f"'groups' must be str or pd.Series, not {type(data).__name__}"
                )
            # same unification as in __init__
            if not names_in_convention:
                self._unify_lipid_data_indices(unification_args)
            self.data.index = self.data.index.astype(str)
        if groups is not None:
            if not hasattr(self, 'groups'):
                self.groups = None
            if self.groups is not None:
                warn("'groups' has already been set, skipping...")
            if isinstance(groups, pd.Series):
                self.groups = groups
            elif isinstance(groups, str):
                self.groups = pd.read_csv(groups)
            else:
                raise ValueError(
                    f"'groups' must be str or pd.Series, not {type(groups).__name__}"
                )
            self.unique_groups = self.groups.unique()
            self.groups.index = self.groups.index.astype(str)
            # obsolete after changing the dict structure
            # self.interaction_attributes = {group: {} for group in self.unique_groups}

    def _check_attributes_(self,
                           edge_colour_attr: Union[str, None],
                           node_colour_attr: Union[str, None]):
        """
        Auxiliary Function to check if edge and node colouring attributes have been
        computed previously

        Parameters
        ----------
        edge_colour_attr : str
            indicating the attribute used for edge colouring
        node_colour_attr :
            str indicating the attribute used for node colouring

        Returns
        -------
        None
        """
        edge_colour_attr_ = {"correlations": "compute_correlations",
                             "partial_correlations": "compute_partial_correlations",
                             "correlation_changes": "compute_correlation_changes",
                             "partial_correlation_changes": "compute_correlation_changes",
                             "reaction_types": "compute_reaction_types",
                             "reaction_enzymes": None}
        # TODO: add stuff like neutral mass?
        node_colour_attr_ = {"fold_changes": "compute_fold_changes",
                             "pvalues": "compute_pvalues",
                             "nlog_pvalues": "compute_pvalues",
                             "lipid_class": None,
                             "c_index": None,
                             "db_index": None,
                             "chain_length": None,
                             "desaturation": None,
                             "betweenness": None,
                             "closeness": None,
                             "degree": None}

        not_in = "'{0}' must be one of {1}, not {2}"

        if edge_colour_attr is not None:
            if edge_colour_attr not in edge_colour_attr_.keys():
                raise ValueError(
                    not_in.format("edge_colour_attr", list(edge_colour_attr_.keys()),
                                  edge_colour_attr)
                )
            if self.interaction_attributes is not None:
                if edge_colour_attr not in self.interaction_attributes.keys():
                    raise NotComputedError(edge_colour_attr, edge_colour_attr_[edge_colour_attr])
            else:
                raise NotComputedError(edge_colour_attr, edge_colour_attr_[edge_colour_attr])
        if node_colour_attr is not None:
            if node_colour_attr not in node_colour_attr_.keys():
                raise ValueError(
                    not_in.format("node_colour_attr", list(node_colour_attr_.keys()),
                                  node_colour_attr)
                )
            if self.lipid_attributes is not None:
                if node_colour_attr not in self.lipid_attributes.keys():
                    raise NotComputedError(node_colour_attr, node_colour_attr_[node_colour_attr])
            else:
                raise NotComputedError(node_colour_attr, node_colour_attr_[node_colour_attr])

    def show_class_connections(self,
                               figsize: tuple = (16, 9),
                               savepath: str = None,
                               return_graph: bool = False,
                               ax: plt.axis = None,
                               **kwargs) -> Union[nx.Graph, plt.axis]:
        return show_class_connections(self.class_connections,
                                      figsize=figsize,
                                      savepath=savepath,
                                      return_graph=return_graph,
                                      ax=ax, **kwargs)

    def _add_to_network_(self, src: LipidSpecies, tgt: LipidSpecies):
        if src.is_transformable(tgt):
            self.network.add_edge(src.name, tgt.name)
            react_type = src.transform_type(tgt)
            self.interaction_attributes["reaction_types"].loc[src.name, tgt.name] = react_type
            enz = _get_enzyme_(react_type, self.class_enzymes,
                               src, tgt)
            if enz is not None:
                self.interaction_attributes["reaction_enzymes"].loc[src.name, tgt.name] = enz

    def add_lipid(self, species: str, lipid_level: str,
                  name_in_convention: bool = False,
                  **kwargs):
        """
        Adding a new lipid species

        Parameters
        ----------
        species: str
            Lipid name. Must be in the LipidLynxX naming convention
            or compatible with LipidLynxX conversion

        lipid_level: str
            'sum' or 'molecular' ('sn' not supported yet)
            Specifies the 'resolution' of lipid species. Sum species
            will sum up the indices of all fatty acids in a lipid
            whereas MolecularSpecies retains individual fatty acid
            resolution without sn-positions, while snSpecies includes
            them as well.

        name_in_convention: bool, optional, default False
            Whether lipid species name is already in the right format.

        Returns
        -------
        None
        """
        if not name_in_convention:
            if 'level' not in kwargs:
                species = unify_name(species, self.class_fas,
                                     level=lipid_level, **kwargs)
            else:
                species = unify_name(species, self.class_fas, **kwargs)
        lipid = _lipid_to_object_(species, lipid_level, self.class_fas)
        # TODO: add in the option to add lipid data as well
        if self.network is not None:
            for tgt in self.lipids:
                self._add_to_network_(lipid, tgt)
                self._add_to_network_(tgt, lipid)
        self.lipids.append(lipid)
        self.lipid_names.append(lipid.name)
        self._name_map_[lipid.name] = len(self.lipids) - 1
        self.lipid_attributes["lipid_class"][lipid.name] = lipid.lipid_class
        self.lipid_attributes["c_index"][lipid.name] = lipid.sum_composition.c_index / lipid.nFA
        self.lipid_attributes["db_index"][lipid.name] = lipid.sum_composition.db_index / lipid.nFA
        self.lipid_attributes["chain_length"][lipid.name] = lipid.sum_composition.c_index
        self.lipid_attributes["desaturation"][lipid.name] = lipid.sum_composition.db_index

    def compute_network(self, directed: bool = True):
        """
        Computing the lipid network from all lipids in self.lipids.
        Mainly relies on LipidSpecies.is_transformable and LipidSpecies.is_equal

        Parameters
        ----------
        directed: bool, optional, default True
            Whether to assume reactions to be directed. This means that
            lipid class connections are not symmetric, i.e. it is possible
            that i -> j but not j -> i.
        """
        self.directed = directed

        if self.network is not None:
            return
            # constructing edge list
        reactions = pd.DataFrame(columns=self.lipid_names,
                                 index=self.lipid_names)
        reaction_enzymes = pd.DataFrame(columns=self.lipid_names,
                                        index=self.lipid_names)
        edges = dict()
        # if directed we want to loop over all lipid pairs
        # else we only need to look at half the comparisons
        if directed:
            i_range = range(self.n_lipids)

            def j_range(obj, *args):
                return range(obj.n_lipids)
        else:
            i_range = range(self.n_lipids - 1)

            def j_range(obj, x):
                return range(x + 1, obj.n_lipids)

        for i in i_range:
            for j in j_range(self, i):
                if self.lipids[i].is_transformable(self.lipids[j],
                                                   self.class_connections,
                                                   self.fas):
                    li = self.lipids[i]
                    lj = self.lipids[j]
                    if edges.get(li.name) is None:
                        edges.setdefault(li.name,
                                         {lj.name: 1})
                    else:
                        edges[li.name][lj.name] = 1
                    rt = li.transform_type(lj)
                    enz = _get_enzyme_(rt, self.class_enzymes,
                                       li, lj)
                    if enz is not None:
                        reaction_enzymes.iloc[i, j] = enz
                    reactions.iloc[i, j] = rt
        # Graph from edge dictionary
        # TODO: test directionality
        if directed:
            self.network = nx.DiGraph(edges)
        else:
            self.network = nx.Graph(edges)
        # adding reaction types to edges
        nx.set_edge_attributes(
            self.network,
            {edge: reaction_enzymes.loc[edge[0], edge[1]] for edge in self.network.edges},
            name="reaction_enzymes"
        )
        # storing reaction types
        self.interaction_attributes["reaction_types"] = reactions
        self.interaction_attributes["reaction_enzymes"] = reaction_enzymes
        # adding unconnected lipids
        for lipid in self.lipids:
            if lipid.name not in self.network.nodes:
                self.network.add_node(lipid.name)
        # setting network-metrics
        self.lipid_attributes["degree"] = pd.Series(dict(self.network.degree()))
        self.lipid_attributes["betweenness"] = pd.Series(dict(nx.betweenness_centrality(self.network)))
        self.lipid_attributes["closeness"] = pd.Series(dict(nx.closeness_centrality(self.network)))
        # TODO: add edge "network" attributes

    def _add_network_attribute_(self, attr: str, nodes: bool = True,
                                group_subset: Union[tuple, list, set] = None,
                                overwrite: bool = True,
                                add_group: bool = False):
        if nodes:
            if self.lipid_attributes.get(attr) is None:
                raise NotComputedError(
                    f"{attr} has not been found in the list of computed "
                    "lipid properties. Please call the appropriate function first!"
                )
        else:
            if self.interaction_attributes.get(attr) is None:
                raise NotComputedError(
                    f"{attr} has not been found in the list of computed "
                    "interaction properties. Please call the appropriate function first!"
                )

        if self.network is None:
            self.compute_network()
        if not overwrite:
            if nodes:
                old = nx.get_node_attributes(self.network, attr)
            else:
                old = nx.get_edge_attributes(self.network, attr)
            if old is not None:
                if not isinstance(old, dict):
                    return
                # some attributes are dicts => need to check whether
                # they are empty
                if old:
                    return

        attr_size = self._attribute_group_size_[attr]
        if self.groups is not None:
            if attr_size > 0:
                # checking if correct number of groups are provided
                # if no groups are provided but the number of unique
                # groups in self are the same as the number of groups
                # required for the respective attribute, they are used
                if group_subset is None:
                    if self.unique_groups.size == attr_size:
                        group_subset = self.unique_groups
                    else:
                        raise ValueError(
                            f"'{attr}' requires exactly {attr_size} groups, "
                            f"but {self.unique_groups.size} are specified in the class instance. "
                            "Please use 'group_subset' to specify a subset of groups."
                        )
                elif len(group_subset) != attr_size:
                    raise ValueError(
                        f"'{attr}' requires exactly {attr_size} groups, "
                        f"but {len(group_subset)} are given."
                    )
                # actual computations with the correct number of groups
                if len(group_subset) == 1:
                    group = group_subset[0]
                    if add_group:
                        attr_ = f"{attr}_{group}"
                    else:
                        attr_ = attr
                    if nodes:
                        if not nx.get_node_attributes(self.network, attr_):
                            for node in self.network.nodes:
                                self.network.nodes[node][attr_] = {}
                        for node in self.network.nodes:
                            to_add = self.lipid_attributes[attr][group].loc[node]
                            if hasattr(to_add, "item"):
                                # json.serialize cannot handle numpy data types
                                to_add = to_add.item()
                            if add_group:
                                self.network.nodes[node][attr_] = to_add
                            else:
                                self.network.nodes[node][attr][group] = to_add
                    else:
                        if not nx.get_edge_attributes(self.network, attr_):
                            for edge in self.network.edges:
                                self.network.edges[edge][attr_] = {}
                        for edge in self.network.edges:
                            to_add = self.interaction_attributes[attr][group].loc[edge[0], edge[1]]
                            if hasattr(to_add, "item"):
                                # json.serialize cannot handle numpy data types
                                to_add = to_add.item()
                            if add_group:
                                self.network.edges[edge][attr_] = to_add
                            else:
                                self.network.edges[edge][attr][group] = to_add
                else:
                    # NOTE: only 1 or 2 possible (no more groups allowed)
                    # => would work for more than 2 groups as well!
                    attr_computed = False
                    for perm in permutations(group_subset):
                        if nodes:
                            if self.lipid_attributes[attr].get(perm) is not None:
                                attr_computed = True
                                break
                        else:
                            if self.interaction_attributes[attr].get(perm) is not None:
                                attr_computed = True
                                break
                    if not attr_computed:
                        raise NotComputedError(attr, group_subset, for_subset=True)
                    group = tuple(group_subset)
                    if group is None:
                        raise ValueError(
                            "Groups in 'group_subset' are not matching the group "
                            "group names provided in the class instance."
                        )
                    if add_group:
                        attr_ = f"{attr}_{_tuple_to_string_(group)}"
                    else:
                        attr_ = attr
                    if nodes:
                        nx.set_node_attributes(self.network,
                                               self.n_lipids * [""],
                                               name=attr_)
                        for node in self.network.nodes:
                            to_add = self.lipid_attributes[attr][group].loc[node]
                            if hasattr(to_add, "item"):
                                # json.serialize cannot handle numpy data types
                                to_add = to_add.item()
                            self.network.nodes[node][attr_] = to_add
                    else:
                        nx.set_edge_attributes(self.network,
                                               len(self.network.edges) * [""],
                                               name=attr_)
                        for edge in self.network.edges:
                            to_add = self.interaction_attributes[attr][group].loc[edge[0], edge[1]]
                            if hasattr(to_add, "item"):
                                # json.serialize cannot handle numpy data types
                                to_add = to_add.item()
                            self.network.edges[edge][attr_] = to_add
            else:
                if nodes:
                    for node in self.network.nodes:
                        to_add = self.lipid_attributes[attr].loc[node]
                        if hasattr(to_add, "item"):
                            # json.serialize cannot handle numpy data types
                            to_add = to_add.item()
                        self.network.nodes[node][attr] = to_add
                else:
                    for edge in self.network.edges:
                        to_add = self.interaction_attributes[attr].loc[edge[0], edge[1]]
                        if hasattr(to_add, "item"):
                            # json.serialize cannot handle numpy data types
                            to_add = to_add.item()
                        self.network.edges[edge][attr] = to_add
        else:
            if nodes:
                for node in self.network.nodes:
                    to_add = self.lipid_attributes[attr].loc[node]
                    if hasattr(to_add, "item"):
                        # json.serialize cannot handle numpy data types
                        to_add = to_add.item()
                    self.network.nodes[node][attr] = to_add
            else:
                for edge in self.network.edges:
                    to_add = self.interaction_attributes[attr].loc[edge[0], edge[1]]
                    if hasattr(to_add, "item"):
                        # json.serialize cannot handle numpy data types
                        to_add = to_add.item()
                    if attr == "reaction_enzymes" and not self.directed:
                        to_add += f"<br>{self.network.edges[(edge[1], edge[0])][attr]}"
                    self.network.edges[edge][attr] = to_add

    @staticmethod
    def _set_min_max_scale_(
        values: np.ndarray,
        vmax: float = None,
        min_to_zero: bool = True
    ) -> Tuple[float, float]:
        min_ = values.min()
        if min_ < 0:
            if vmax is None:
                max_ = max(abs(values.min()),
                           abs(values.max()))
            else:
                max_ = vmax
            min_ = -max_
        else:
            if min_to_zero:
                min_ = 0
            if vmax is None:
                max_ = values.max()
            else:
                max_ = vmax
        return min_, max_

    def _generate_edge_colours_(
        self, attr: str, discrete: bool = True,
        cmap: str = "tab10",
        group: Union[str, tuple] = None,
        vmax: float = None,
        colours_to_hex: bool = True,
        group_attr: str = None
    ) -> Tuple[dict, list]:
        # TODO: make this work with groups!
        #  if self.edge_colours.get(attr):
        #      return self.edge_colours.get(attr)
        if attr == "correlation_changes" or attr == "partial_correlation_changes":
            if colours_to_hex:
                edge_colour_map = {
                    "significant to unsignificant": to_hex("tab:cyan"),
                    "unsignificant to significant": to_hex("tab:green"),
                    "unchanged significant": to_hex("tab:red"),
                    "unsignificant": "#e1e5e7",
                    "positive to negative": to_hex("tab:blue"),
                    "negative to positive": to_hex("tab:orange")
                }
            else:
                edge_colour_map = {
                    "significant to unsignificant": "tab:cyan",
                    "unsignificant to significant": "tab:green",
                    "unchanged significant": "tab:red",
                    "unsignificant": "#e1e5e7",
                    "positive to negative": "tab:blue",
                    "negative to positive": "tab:orange"
                }
            edge_colours = [edge_colour_map[self.interaction_attributes[attr][group].loc[edge[0], edge[1]]]
                            for edge in self.network.edges]
            return edge_colour_map, edge_colours
        if self._attribute_group_size_[attr] == 0:
            group = None
        if discrete:
            if group is None:
                attr_set = self.interaction_attributes[attr]
                if isinstance(attr_set, dict):
                    keys = attr_set.keys()
                    if len(keys) == 1:
                        attr_set = attr_set[list(keys)[0]]
                    else:
                        raise ValueError(
                            f"Multiple group options found for {attr}."
                            " Please specify which should be used by using the 'groups' parameter"
                        )
                attr_val_set = unique_elements(attr_set)
            else:
                attr_val_set = unique_elements(self.interaction_attributes[attr][group])
            if attr_val_set.size > 10 and cmap == "tab10":
                cmap = "tab20"
            if attr_val_set.size > 20:
                custom_map = _generate_colormap_(attr_val_set.size)
                colours = [custom_map(i) for i in range(attr_val_set.size)]
            else:
                colours = [plt.get_cmap(cmap)(i)
                           for i in range(attr_val_set.size)]
            edge_colour_map = dict(zip(attr_val_set, colours))
            if group_attr is not None:
                attr = group_attr
            if colours_to_hex:
                edge_colours = [to_hex(edge_colour_map.get(self.network.edges[edge][attr], (0.84, 0.84, 0.84)))
                                for edge in self.network.edges]
            else:
                edge_colours = [edge_colour_map.get(self.network.edges[edge][attr], (0.84, 0.84, 0.84))
                                for edge in self.network.edges]
        else:
            if cmap == "tab10":
                # NOTE: seismic is bad here because 0 values will be non-visible
                cmap = "coolwarm"
            if group is None:
                vals = self.interaction_attributes[attr]
                if isinstance(vals, dict):
                    keys = vals.keys()
                    if len(keys) == 1:
                        vals = vals[list(keys)[0]]
                    else:
                        raise ValueError(
                            f"Multiple group options found for {attr}."
                            " Please specify which should be used by using the 'groups' parameter"
                        )
            else:
                vals = self.interaction_attributes[attr][group]
                if vals.isna().values.all():
                    vals = 5

            min_, max_ = self._set_min_max_scale_(vals.values, vmax)
            if np.isnan(min_):
                min_ = vals.min(skipna=True)
                if np.isnan(min_):
                    min_ = 0
            if np.isnan(max_) or min_ == max_:
                max_ = vals.max(skipna=True).max(skipna=True)
                if np.isnan(max_):
                    max_ = min_ + 1
            scm = ScalarMappable(
                norm=Normalize(
                    vmin=-max_,
                    vmax=max_
                ),
                cmap=cmap
            )
            if colours_to_hex:
                edge_colours = [to_hex(scm.to_rgba(vals.loc[edge[0], edge[1]]))
                                for edge in self.network.edges]
            else:
                edge_colours = [scm.to_rgba(vals.loc[edge[0], edge[1]])
                                for edge in self.network.edges]
            # edge_colour_map = {min_: plt.get_cmap("viridis")(0),
            #                    max_: plt.get_cmap("viridis")(1)}
            edge_colour_map = {"min": -max_,
                               "max": max_,
                               "cmap": cmap}

        self.edge_colours[attr] = (edge_colour_map, edge_colours)
        return edge_colour_map, edge_colours

    def _generate_node_colours(
        self, attr: str, discrete: bool = True,
        cmap: str = "tab10",
        group: Union[str, tuple] = None,
        vmin: float = None, vmax: float = None,
        colours_to_hex: bool = False,
        group_attr: str = None
    ) -> Tuple[dict, list]:
        if attr == "lipid_class":
            if colours_to_hex:
                # TODO: test mapping of unknown lipid colours
                node_colours = [
                    LIPID_COLOURS.get(self.network.nodes[node]["lipid_class"], "#7f7f7f")
                    for node in self.network.nodes
                ]
                # NOTE: this is necessary to also include lipid classes not covered by our default
                # colour scheme
                node_colour_map = {
                    lipid_class: LIPID_COLOURS.get(lipid_class, "#7f7f7f")
                    for lipid_class in self.lipid_attributes["lipid_class"].unique()
                }
            else:
                node_colours = [
                    to_rgba(LIPID_COLOURS.get(self.network.nodes[node]["lipid_class"], "tab:gray"))
                    for node in self.network.nodes
                ]
                # NOTE: this is necessary to also include lipid classes not covered by our default
                # colour scheme
                node_colour_map = {
                    lipid_class: to_rgba(LIPID_COLOURS.get(lipid_class, "tab:gray"))
                    for lipid_class in self.lipid_attributes["lipid_class"].unique()
                }
        else:
            if self._attribute_group_size_[attr] == 0:
                group = None
            if discrete:
                if group is None:
                    attr_val_set = unique_elements(self.lipid_attributes[attr])
                else:
                    attr_val_set = unique_elements(self.lipid_attributes[attr][group])
                if attr_val_set.size > 10 and cmap == "tab10":
                    cmap = "tab20"
                if attr_val_set.size > 20:
                    colours = _generate_colormap_(attr_val_set.size)
                else:
                    colours = [plt.get_cmap(cmap)(i)
                               for i in range(attr_val_set.size)]
                node_colour_map = dict(zip(attr_val_set, colours))
                if group_attr is not None:
                    attr = group_attr
                if colours_to_hex:
                    node_colours = [to_hex(node_colour_map[self.network.nodes[node][attr]])
                                    for node in self.network.nodes]
                else:
                    node_colours = [node_colour_map[self.network.nodes[node][attr]]
                                    for node in self.network.nodes]
            else:
                if cmap == "tab10":
                    cmap = "seismic"
                if group is None:
                    vals = self.lipid_attributes[attr]
                else:
                    vals = self.lipid_attributes[attr][group]

                index_attr = attr not in ["c_index", "db_index", "chain_length", "desaturation"]
                min_, max_ = self._set_min_max_scale_(vals.values, vmax,
                                                      min_to_zero=index_attr)
                if np.isnan(min_):
                    min_ = vals.min(skipna=True)
                    if np.isnan(min_):
                        min_ = 0
                if np.isnan(max_) or min_ == max_:
                    max_ = vals.max(skipna=True)
                    if np.isnan(max_):
                        max_ = min_ + 1
                if vmin is not None:
                    min_ = vmin
                scm = ScalarMappable(
                    norm=Normalize(
                        vmin=min_,
                        vmax=max_
                    ),
                    cmap=cmap
                )
                if colours_to_hex:
                    node_colours = [to_hex(scm.to_rgba(vals.loc[node]))
                                    for node in self.network.nodes]
                else:
                    node_colours = [scm.to_rgba(vals.loc[node])
                                    for node in self.network.nodes]
                # node_colour_map = {min_: plt.get_cmap("viridis")(0),
                #                    max_: plt.get_cmap("viridis")(1)}
                if np.isnan(vals.values).any():
                    node_colour_map = {"min": min_,
                                       "max": max_,
                                       "nan": np.nan,
                                       "cmap": cmap}
                else:
                    node_colour_map = {"min": min_,
                                       "max": max_,
                                       "cmap": cmap}

        self.node_colours[attr] = (node_colour_map, node_colours)
        return node_colour_map, node_colours

    def _discrete_legend_(self, attr: str, nodes: bool = True,
                          **kwargs) -> List[Line2D]:
        # TODO: adapt to group-specific
        if nodes:
            # node labels => scatters
            legend_map = self.node_colours[attr][0]
            markersize = kwargs.pop("markersize", 12)
            handles = [Line2D([0], [0], color="w", marker="o",
                              markerfacecolor=col, label=group,
                              markersize=markersize,
                              **kwargs)
                       for group, col in legend_map.items()]
        else:
            # edge labels => thick lines
            legend_map = self.edge_colours[attr][0]
            lw = kwargs.pop("lw", 4)
            handles = [Line2D([0], [0], color=col,
                              lw=lw, label=group,
                              **kwargs)
                       for group, col in legend_map.items()]

        return handles

    def _continuous_legend_(self, attr: str, ax: plt.axis = None,
                            nodes: bool = True,
                            **kwargs) -> Colorbar:
        if nodes:
            legend_map = self.node_colours[attr][0]
        else:
            legend_map = self.edge_colours[attr][0]
        scm = ScalarMappable(norm=Normalize(vmin=legend_map["min"],
                                            vmax=legend_map["max"]),
                             cmap=legend_map["cmap"])
        return plt.colorbar(scm, ax=ax, **kwargs)

    def save_legend(self, attr: str, file: str,
                    nodes: bool = True,
                    groups: Union[list, str] = None,
                    **kwargs):
        """
        Saving colour legend for a given attribute as an image

        Parameters
        ----------
        attr: str
            Name of the attribute for which the legend should
            be computed.

        file: str
            Filename to save the legend to.

        nodes: bool, optional, default True
            Whether attribute is considering nodes (True) or
            edge (False)

        groups: list or str, optional
            Currently not implemented!
            specifies which groups to use, e.g. for p-values,
            correlation changes etc.

        kwargs:
            keyword arguments to pass to plt.savefig
        """
        # TODO: adapt to group-specific
        discrete = self._discrete_map_[attr]
        if discrete:
            fig, ax = plt.subplots(figsize=(12, 9))
            handles = self._discrete_legend_(attr, nodes=nodes,
                                             markersize=100)
            ax.legend(handles=handles,
                      loc="center", fontsize=100)
            plt.axis("off")
            plt.savefig(file, **kwargs)
        else:
            fig, ax = plt.subplots(figsize=(7, 5))
            cbar = self._continuous_legend_(attr, ax=ax,
                                            nodes=nodes,
                                            fraction=1, pad=0)
            for tick in cbar.ax.get_yticklabels():
                tick.set_fontsize(23)
            # => actually fills the whole plot!
            plt.axis("off")
            plt.savefig(file, **kwargs)

    # TODO: use save_legend to put legends into interactive plot via html file modification?

    @staticmethod
    def save_size_legend(file: str, extrema: tuple,
                         scale: tuple = None,
                         nodes: bool = True, **kwargs):
        """
        Save legend for node/edge sizes to a file

        Parameters
        ----------
        file: str
            Filename to save the legend to.

        extrema: tuple
            Minimum and maximum value of the original values

        scale: tuple
            Minimum and maximum value of the scaled values

        nodes: bool, default True
            Whether attribute is considering nodes (True) or
            edge (False)

        kwargs:
            keyword arguments to pass to plt.savefig
        """
        fig, ax = plt.subplots(figsize=(12, 9))
        handles = _size_legend_(extrema, scale, nodes)
        ax.legend(handles=handles,
                  loc="center", fontsize=100)
        plt.axis("off")
        plt.savefig(file, **kwargs)

    @staticmethod
    def _scale_dict_(
        sizes: Dict[str, float],
        scale: Tuple[float, float],
        abs_: bool = False,
        map_as_ex: bool = False,
        vmax: float = None
    ) -> Tuple[Dict[str, float], Dict[float, float]]:
        scale_sizes = np.array(list(sizes.values()))
        if np.isnan(scale_sizes).all():
            # arbitraty default
            scale_sizes = np.zeros(scale_sizes.shape) + 5
        # if abs_:
        #     scale_sizes = abs(scale_sizes)  # NOTE: in this case min = 0
        # else:
        nans = np.isnan(scale_sizes)
        if vmax is None:
            max_ = scale_sizes[~nans].max()
        else:
            max_ = vmax
        if abs_:
            min_ = 0
            scale_sizes = abs(scale_sizes)
        else:
            min_ = scale_sizes[~nans].min()
            scale_sizes = scale_sizes - min_
        if max_ - min_ != 0:
            scale_sizes = scale_sizes / (max_ - min_)
        if abs_:
            scale_sizes = abs(scale_sizes * (scale[1] - scale[0]) + scale[0])
        else:
            scale_sizes = scale_sizes * (scale[1] - scale[0]) + scale[0]
        scaled_sizes = dict(zip(sizes.keys(),
                                list(scale_sizes)))
        final_min = scale_sizes[~nans].min()
        final_max = scale_sizes[~nans].max()
        if map_as_ex:
            if nans.any():
                scale_map = {
                    "min": (min_, final_min),
                    "max": (max_, final_max),
                    # arbitrary values to draw nans in legend
                    "nan": (1, 1)
                }
            else:
                scale_map = {
                    "min": (min_, final_min),
                    "max": (max_, final_max)
                }
        else:
            scale_sizes[nans] = final_min / 2
            scale_map = dict(zip(sizes.values(), scale_sizes))
        return scaled_sizes, scale_map

    # TODO: add parameters
    # * interactive (pyvis) vs. static (networkx)
    # use external functions?
    def plot_static_network(
        self,
        layout: Union[Callable, dict] = nx_layout.fruchterman_reingold_layout,
        overwrite_layout: bool = False,
        as_undirected: bool = True,
        edge_colour_attr: str = None,
        overwrite_edge_attr: bool = True,
        edge_map: str = "tab10",
        edge_legend: bool = True,
        edge_size_func: Union[Callable, int] = None,
        edge_size_scale: Tuple[int, int] = (3, 12),
        edge_size_legend: bool = True,
        edge_size_title: str = "Edge size",
        node_colour_attr: str = None,
        overwrite_node_attr: bool = True,
        node_map: str = "tab10",
        node_legend: bool = True,
        node_size_func: Union[Callable, int, dict] = None,
        node_size_scale: Tuple[int, int] = (3, 12),
        node_size_legend: bool = True,
        node_size_title: str = "Node size",
        edge_group_subset: Union[tuple, list, set] = None,
        node_group_subset: Union[tuple, list, set] = None,
        ax: plt.axis = None,
        figargs: dict = None,
        layout_args: dict = None,
        vmin: float = None,
        vmax: float = None,
        **kwargs
    ) -> plt.axis:
        """
        Plotting the inferred network in a static manner using networkx and matplotlib

        Parameters
        ----------
        layout: networkx.layout function or dict, optional, default nx_layout.fruchterman_reingold_layout
            If a dict specifying node position this is used to
            plot the network. If a callable it is used to
            compute the layout.
            If self.layout is None or overwrite_layout is True
            it will be stored in self.layout

        overwrite_layout: bool, optional, default False
            Whether to overwrite self.layout

        as_undirected: bool, optional, default True
            Whether to show the network as directed or undirected.
            If directed, lipid class reactions will not be symmetric,
            i.e. PC -> PE could be possible but not PE -> PC.

        edge_colour_attr: str, optional
            Attribute to use for edge colours.
            Available options: correlations, partial_correlations,
            reaction_types, correlation_changes, partial_correlation_changes.
            All metricx except reaction_types have to be computed beforehand.

        overwrite_edge_attr: bool, optional, default True
            Whether to overwrite attribute with the given group setting in
            self.interaction_attributes

        edge_map: str, optional, default 'tab10'
            Colourmap to use for edge colouring. Default for continuous
            attributes is 'seismic'

        edge_legend: bool, optional, default True
            Whether to show a legend for edge colours

        edge_size_func: Callable or int, optional
            Used to compute edge sizes if a function.
            Else set as constant edge size.

        edge_size_scale: tuple, optional, default (3, 12)
            Tuple of integers specifying the minimum and maximum
            of plotted edge sizes. The actual values will be scaled
            to fit this range.

        edge_size_legend: bool, optional, True
            Whether to show the legend for edge sizes

        edge_size_title: str, optional, default 'Edge size'
            Title of the edge size legend

        node_colour_attr: str, optional
            Specifies the attribute to use for determining node colours.
            Available options are: fold_changes, pvalues, nlog_pvalues,
            lipid_class, c_index, db_index

        overwrite_node_attr: bool, optional, default True
            Whether to overwrite attribute with the given group setting in
            self.lipid_attributes

        node_map: str, optional, default 'tab10'
            Colourmap to use for node colouring. Default for continuous
            attributes is 'seismic'

        node_legend: bool, optional, default True
            Whether to show node colour legend

        node_size_func: Callable or int, optional
            Used to compute node sizes if a function.
            Else set as constant node size.

        node_size_scale: tuple, optional, default (3, 12)
            Tuple of integers specifying the minimum and maximum
            of plotted node sizes. The actual values will be scaled
            to fit this range.

        node_size_legend: bool, optional, default True
            Whether to show the legend for node sizes

        node_size_title: str, optional, default 'Node size'
            Title of node size legend

        edge_group_subset: tuple or list of size 2, optional
            Group-comparison to show in edge colouring/sizing

        node_group_subset: tuple or list of size 2, optional
            Group-comparison to show in node colouring/sizing

        ax: plt.axis, optional
            matplotlib axis to plot onto

        figargs: dict, optional
            keyword arguments to pass to plt.figure.
            Only relevant if ax=None (default)

        layout_args: dict, optional
            keyword arguments to pass to layout.
            Only relevant if layout is not a function (default)

        vmin: float, optional
            Minimum value of the legend (for continuous scales).
            If not given it will be set to 0 (for only positive
            arrays) or to -max(abs(min), abs(max)) (for zero-centred
            colours).
            Setting vmin can lead to un-centred colour scales!

        vmax: float, optional
            Maximum value of the legend (for continuous scales).
            If vmin is not specified, the minimum value will be
            either 0 (if all values are positive) or -vmax

        kwargs:
            keyword arguments to pass to networkx.draw


        Returns
        -------
        plt.axis on which the network is plotted


        Raises
        ------
        NotComputedError
            If edge_colour_attr or node_colour_attr have not been
            computed before

        ValueError
            If edge_colour_attr or node_colour_attr are not in the
            list of available options
        """
        if self.network is None:
            self.compute_network()
        # TODO: change this to function properly
        # check if attributes have been computed
        edge_discrete = None
        node_discrete = None
        if edge_colour_attr is not None:
            self._check_attributes_(edge_colour_attr, node_colour_attr)
            # check if attributes have been added to network
            # edge attribute
            self._add_network_attribute_(edge_colour_attr, nodes=False,
                                         group_subset=edge_group_subset,
                                         overwrite=overwrite_edge_attr)
            edge_discrete = self._discrete_map_[edge_colour_attr]
            if edge_group_subset is None:
                edge_colours = self._generate_edge_colours_(edge_colour_attr,
                                                            edge_discrete,
                                                            edge_map)
            elif len(edge_group_subset) == 1:
                edge_colours = self._generate_edge_colours_(edge_colour_attr,
                                                            edge_discrete,
                                                            edge_map,
                                                            group=edge_group_subset[0])
            else:
                edge_colours = self._generate_edge_colours_(edge_colour_attr,
                                                            edge_discrete,
                                                            edge_map,
                                                            group=edge_group_subset)
        else:
            edge_legend = False
            edge_colours = ["tab:gray" for _ in range(len(self.network.edges))]
        # node attribute
        if node_colour_attr is not None:
            self._add_network_attribute_(node_colour_attr, nodes=True,
                                         group_subset=node_group_subset,
                                         overwrite=overwrite_node_attr)
            node_discrete = self._discrete_map_[node_colour_attr]
            if node_group_subset is None:
                node_colours = self._generate_node_colours(node_colour_attr,
                                                           node_discrete,
                                                           node_map,
                                                           vmax=vmax,
                                                           vmin=vmin)
            elif len(node_group_subset) == 1:
                node_colours = self._generate_node_colours(node_colour_attr,
                                                           node_discrete,
                                                           node_map,
                                                           group=node_group_subset[0],
                                                           vmax=vmax,
                                                           vmin=vmin)
            else:
                node_colours = self._generate_node_colours(node_colour_attr,
                                                           node_discrete,
                                                           node_map,
                                                           group=node_group_subset,
                                                           vmax=vmax,
                                                           vmin=vmin)
        else:
            node_legend = False
            node_colours = ["tab:blue" for _ in range(len(self.network.nodes))]

        # node sizes:
        # ===========
        node_size_legend_ = False
        node_sizes = None
        if kwargs.get("node_size") is not None:
            node_size = kwargs.pop("node_size")
            if isinstance(node_size, int) or isinstance(node_size, float):
                node_size = {node: node_size for node in self.network.nodes}
            elif isinstance(node_size, list) or isinstance(node_size, set):
                node_size = dict(zip(self.network.nodes, node_size))
        else:
            if node_size_func is None:
                sizes = nx.betweenness_centrality(self.network)
                scale_sizes = _range_scale_(x=np.array(list(sizes.values())),
                                            a=50, b=400)
                node_size = dict(zip(sizes.keys(),
                                     list(scale_sizes)))
                node_size_legend_ = True
                node_sizes = (min(list(sizes.values())),
                              max(list(sizes.values())))
            elif isinstance(node_size_func, int):
                node_size = dict(zip(self.network.nodes,
                                     len(self.network.nodes) * [node_size_func]))
            elif isinstance(node_size_func, Callable):
                node_size = node_size_func(self.network)
                node_size_legend_ = True
            else:
                node_size = dict(zip(self.network.nodes,
                                     len(self.network.nodes) * [23]))

        # edge sizes:
        # ===========
        edge_sizes = None
        if kwargs.get("width") is not None:
            edge_size = kwargs.pop("width")
            edge_size_legend_ = False
        else:
            edge_size_legend_ = False
            if edge_size_func is None:
                sizes = nx.edge_betweenness_centrality(self.network)
                scale_sizes = _range_scale_(x=np.array(list(sizes.values())),
                                            a=.2, b=3)
                edge_size = dict(zip(sizes.keys(),
                                     list(scale_sizes)))
                edge_size_legend_ = True
                edge_sizes = (min(list(sizes.values())),
                              max(list(sizes.values())))
            elif isinstance(edge_size_func, int):
                edge_size = dict(zip(self.network.edges,
                                     len(self.network.edges) * [edge_size_func]))
            elif isinstance(edge_size_func, Callable):
                edge_size = edge_size_func(self.network)
                if not isinstance(edge_size, dict):
                    edge_size = dict(zip(self.network.edges,
                                         edge_size))
                edge_size_legend_ = True
            else:
                edge_size = dict(zip(self.network.edges,
                                     len(self.network.edges) * [5]))

        if ax is None:
            if figargs is None:
                fig, ax = plt.subplots()
            else:
                fig, ax = plt.subplots(**figargs)
            ax.set_xticks([])
            ax.set_yticks([])

        if self.layout is None or overwrite_layout:
            if isinstance(layout, dict):
                self.layout = layout
            elif layout_args is None:
                self.layout = layout(self.network)
            else:
                self.layout = layout(self.network,
                                     **layout_args)
        # TODO: options
        # * node size
        # * edge width
        if as_undirected:
            if isinstance(self.network, nx.DiGraph):
                nx.draw(self.network.to_undirected(),
                        node_color=node_colours[1],
                        edge_color=edge_colours[1],
                        ax=ax,
                        node_size=[node_size[node] for node in self.network.nodes],
                        width=[edge_size[edge] for edge in self.network.edges],
                        pos=self.layout,
                        **kwargs)
            else:
                nx.draw(self.network,
                        node_color=node_colours[1],
                        edge_color=edge_colours[1],
                        ax=ax,
                        node_size=[node_size[node] for node in self.network.nodes],
                        width=[edge_size[edge] for edge in self.network.edges],
                        pos=self.layout,
                        **kwargs)
        else:
            nx.draw(self.network,
                    node_color=node_colours[1],
                    edge_color=edge_colours[1],
                    ax=ax,
                    node_size=[node_size[node] for node in self.network.nodes],
                    width=[edge_size[edge] for edge in self.network.edges],
                    pos=self.layout,
                    **kwargs)
        # colour legends
        # ==============
        if node_legend:
            if node_discrete:
                hands = self._discrete_legend_(node_colour_attr,
                                               nodes=True)
                if edge_discrete:
                    plt.gca().add_artist(
                        plt.legend(handles=hands,
                                   title=node_colour_attr,
                                   loc=1,
                                   bbox_to_anchor=(1.2, 1))
                    )
                else:
                    plt.gca().add_artist(
                        plt.legend(handles=hands,
                                   title=node_colour_attr,
                                   loc=1,
                                   bbox_to_anchor=(1.35, 1))
                    )
            else:
                nbar = self._continuous_legend_(node_colour_attr,
                                                ax, nodes=True)
                nbar.ax.set_title(node_colour_attr)
        if edge_legend:
            if edge_discrete:
                hands = self._discrete_legend_(edge_colour_attr,
                                               nodes=False)
                plt.gca().add_artist(
                    plt.legend(handles=hands,
                               title=edge_colour_attr,
                               loc=1,
                               bbox_to_anchor=(1.35, 1))
                )
            else:
                ebar = self._continuous_legend_(edge_colour_attr,
                                                ax, nodes=False)
                ebar.ax.set_title(edge_colour_attr)
        # size legends
        # ==============
        # TODO:
        # choose node/edge size in a way such that the legend
        # title can represent the actual property
        if node_size_legend and node_size_legend_:
            if node_sizes is None:
                node_sizes = (min(list(node_size.values())),
                              max(list(node_size.values())))
                node_size_scale = None
            plt.gca().add_artist(
                plt.legend(handles=_size_legend_(node_sizes,
                                                 node_size_scale,
                                                 nodes=True),
                           loc="upper left", bbox_to_anchor=(-.15, 1),
                           title=node_size_title)
            )
        if edge_size_legend and edge_size_legend_:
            if edge_sizes is None:
                edge_sizes = (min(list(edge_size.values())),
                              max(list(edge_size.values())))
                edge_size_scale = None
            plt.gca().add_artist(
                plt.legend(handles=_size_legend_(edge_sizes,
                                                 edge_size_scale,
                                                 nodes=False),
                           loc="upper left", bbox_to_anchor=(-.15, .7),
                           title=edge_size_title)
            )

        plt.tight_layout()
        return ax

    def network_to_pyvis(
        self,
        edge_colour_attr: str = None,
        overwrite_edge_attr: bool = True,
        edge_cmap: str = "tab10",
        edge_size_func: Union[Callable, int] = None,
        edge_size_scale: Tuple[int, int] = (.5, 2),
        node_colour_attr: str = None,
        overwrite_node_attr: bool = True,
        node_cmap: str = "tab10",
        node_size_func: Union[Callable, int, dict] = None,
        node_size_scale: Tuple[int, int] = (10, 40),
        edge_group_subset: Union[tuple, list, set] = None,
        node_group_subset: Union[tuple, list, set] = None,
        edge_size_legend: bool = True,
        node_size_legend: bool = True,
        net: Network = None,
        as_directed: bool = False,
        vmin: float = None, vmax: float = None,
        **kwargs
    ) -> VisParser:
        """
        Converting the generated network from networkx format to pyvis format

        Parameters
        ----------
        edge_colour_attr: str, optional
            Attribute to use for edge colours.
            Available options: correlations, partial_correlations,
            reaction_types, correlation_changes, partial_correlation_changes.
            All metricx except reaction_types have to be computed beforehand.

        overwrite_edge_attr: bool, optional, default True
            Whether to overwrite attribute with the given group setting in
            self.interaction_attributes

        edge_cmap: str, optional, default 'tab10'
            Colourmap to use for edge colouring. Default for continuous
            attributes is 'seismic'

        edge_size_func: Callable or int, optional
            Used to compute edge sizes if a function.
            Else set as constant edge size.

        edge_size_scale: tuple, optional, default (3, 12)
            Tuple of integers specifying the minimum and maximum
            of plotted edge sizes. The actual values will be scaled
            to fit this range.

        node_colour_attr: str, optional
            Specifies the attribute to use for determining node colours.
            Available options are: fold_changes, pvalues, nlog_pvalues,
            lipid_class, c_index, db_index

        overwrite_node_attr: bool, optional, default True
            Whether to overwrite attribute with the given group setting in
            self.lipid_attributes

        node_size_func: Callable or int, optional
            Used to compute node sizes if a function.
            Else set as constant node size.

        node_size_scale: tuple, optional, default (3, 12)
            Tuple of integers specifying the minimum and maximum
            of plotted node sizes. The actual values will be scaled
            to fit this range.

        node_cmap: str, optional, default 'tab10'
            Colourmap to use for edge colouring. Default for continuous
            attributes is 'seismic'

        edge_group_subset: tuple or list of size 2, optional
            Group-comparison to show in edge colouring/sizing

        node_group_subset: tuple or list of size 2, optional
            Group-comparison to show in node colouring/sizing

        node_size_legend: bool, optional, default True
            Whether to show node sizes in legend

        edge_size_legend: bool, optional, default True
            Whether to show edge sizes in legend

        net: pyvis.network.Network, optional
            If given this instance will be used instead of generating
            a new one within the function

        as_directed: bool, optional, default False
            Whether to show a direted graph

        vmin: float, optional
            Minimum value of the legend (for continuous scales).
            If not given it will be set to 0 (for only positive
            arrays) or to -max(abs(min), abs(max)) (for zero-centred
            colours).
            Setting vmin can lead to un-centred colour scales!

        vmax: float, optional
            Maximum value of the legend (for continuous scales).
            If vmin is not specified, the minimum value will be
            either 0 (if all values are positive) or -vmax

        kwargs:
            keyword arguments to pass to pyvis.network.Network. Only
            releveant if net=None (default)

        Returns
        -------
        vis_utils.VisParser
        """
        self._check_attributes_(edge_colour_attr, node_colour_attr)
        # check if attributes have been added to network
        # edge attribute
        self._add_network_attribute_(edge_colour_attr, nodes=False,
                                     group_subset=edge_group_subset,
                                     overwrite=overwrite_edge_attr)
        edge_discrete = self._discrete_map_[edge_colour_attr]
        # node attribute
        self._add_network_attribute_(node_colour_attr, nodes=True,
                                     group_subset=node_group_subset,
                                     overwrite=overwrite_node_attr)
        node_discrete = self._discrete_map_[node_colour_attr]

        # generate colours
        if edge_group_subset is None:
            edge_colours = self._generate_edge_colours_(edge_colour_attr,
                                                        edge_discrete,
                                                        edge_cmap)
        elif len(edge_group_subset) == 1:
            edge_colours = self._generate_edge_colours_(edge_colour_attr,
                                                        edge_discrete,
                                                        edge_cmap,
                                                        group=edge_group_subset[0])
        else:
            edge_colours = self._generate_edge_colours_(edge_colour_attr,
                                                        edge_discrete,
                                                        edge_cmap,
                                                        group=edge_group_subset)
        if node_group_subset is None:
            node_colours = self._generate_node_colours(node_colour_attr,
                                                       node_discrete,
                                                       node_cmap,
                                                       vmax=vmax,
                                                       vmin=vmin)
        elif len(node_group_subset) == 1:
            node_colours = self._generate_node_colours(node_colour_attr,
                                                       node_discrete,
                                                       node_cmap,
                                                       group=node_group_subset[0],
                                                       vmax=vmax,
                                                       vmin=vmin)
        else:
            node_colours = self._generate_node_colours(node_colour_attr,
                                                       node_discrete,
                                                       node_cmap,
                                                       group=node_group_subset,
                                                       vmax=vmax,
                                                       vmin=vmin)

        int_net = self.network.copy()
        node_size_dict = None
        if node_size_func is None:
            sizes = nx.betweenness_centrality(int_net)
            node_size, node_size_dict = self._scale_dict_(sizes, node_size_scale)
        elif isinstance(node_size_func, dict):
            node_size, node_size_dict = self._scale_dict_(node_size_func,
                                                          node_size_scale)
        elif isinstance(node_size_func, int):
            node_size = dict(zip(int_net.nodes,
                                 len(int_net.nodes) * [node_size_func]))
        elif isinstance(node_size_func, Callable):
            node_size = node_size_func(int_net)
        else:
            node_size = dict(zip(int_net.nodes,
                                 len(int_net.nodes) * [5]))
        for i, node in enumerate(int_net.nodes):
            int_net.nodes[node]["value"] = node_size[node]
            int_net.nodes[node]["color"] = to_hex(node_colours[1][i])
            node_title = f"<strong>{node}</strong>:"
            prime = "<br><br><strong>{0}</strong>:<br>{1}"
            for prop in ["fold_changes", "pvalues"]:
                if self.lipid_attributes.get(prop) is not None:
                    sup = re.sub("s$", "", prop)
                    group_vals = self.lipid_attributes[prop]
                    # if fold_changes and pvalues are computes groups
                    # must have been specified at some point
                    for group, vals in group_vals.items():
                        node_title += prime.format(sup, f"{group}: {vals[node]}")
            int_net.nodes[node]["title"] = node_title

        edge_size_dict = None
        if edge_size_func is None:
            sizes = nx.edge_betweenness_centrality(int_net)
            edge_size, edge_size_dict = self._scale_dict_(sizes, edge_size_scale)
        elif isinstance(edge_size_func, dict):
            edge_size, edge_size_dict = self._scale_dict_(edge_size_func,
                                                          edge_size_scale)
        elif isinstance(edge_size_func, int):
            edge_size = dict(zip(int_net.edges,
                                 len(int_net.edges) * [edge_size_func]))
        elif isinstance(edge_size_func, Callable):
            edge_size = edge_size_func(int_net)
        else:
            edge_size = dict(zip(int_net.edges,
                                 len(int_net.edges) * [5]))

        props_oi = [prop for prop, val in self.interaction_attributes.items()
                    if val is not None]
        if isinstance(edge_group_subset, list):
            if len(edge_group_subset) == 1:
                def get_prop(x):
                    return self.interaction_attributes[x][edge_group_subset[0]]
            else:
                def get_prop(x):
                    return self.interaction_attributes[x][tuple(edge_group_subset)]
        else:
            def get_prop(x):
                return self.interaction_attributes[x]
        for i, edge in enumerate(int_net.edges):
            int_net.edges[edge]["value"] = edge_size[edge]
            int_net.edges[edge]["color"] = to_hex(edge_colours[1][i])
            edge_title = ""
            prime = "<strong>{0}</strong>:<br>{1}<br><br>"
            for prop in props_oi:
                sup = re.sub("s$", "", prop)
                if self._attribute_group_size_[prop] == 0:
                    if prop != "reaction_enzymes":
                        val = self.interaction_attributes[prop].loc[edge[0], edge[1]]
                    else:
                        if self.directed:
                            renz = self.interaction_attributes[prop].loc[edge[0], edge[1]]
                            if isinstance(renz, str):
                                val = f"{renz.split(':')[1]}"
                            else:
                                val = str(renz)
                        else:
                            val_ij = self.interaction_attributes[prop].loc[edge[0], edge[1]]
                            val_ji = self.interaction_attributes[prop].loc[edge[1], edge[0]]
                            val = f"{val_ij}<br>{val_ji}"
                else:
                    try:
                        vals = get_prop(prop)
                        if isinstance(vals, dict):
                            if edge_group_subset is None:
                                raise ValueError(
                                    f"{prop} has been computed for groups. "
                                    "Please specify an edge group subset"
                                )
                            vals = vals[edge_group_subset]
                        val = vals.loc[edge[0], edge[1]]
                    # FIXME
                    # this is not a nice solution => make sure it's possible to
                    # correctly get edge values when supplying single groups
                    # for nodes
                    except KeyError:
                        val = None
                if val is not None:
                    if isinstance(val, float):
                        edge_title += prime.format(sup, round(val, ndigits=3))
                    else:
                        edge_title += prime.format(sup, val)
            int_net.edges[edge]["title"] = re.sub("<br>$", "", edge_title)
        # TODO: add edge width and node size as value
        # TODO: generate legend here (and return with network)
        if net is None:
            net = VisParser(directed=as_directed, **kwargs)
        else:
            net = VisParser.from_pyvis_network(net, **kwargs)
        net.inherit_edge_colors(False)
        if as_directed:
            net.from_nx(int_net.to_directed())
        else:
            net.from_nx(int_net)
        # TODO: test this properly
        if not node_size_legend:
            node_size_dict = None
        if not edge_size_legend:
            edge_size_dict = None
        net.generate_legend(
            node_colours=node_colours[0],
            node_sizes=node_size_dict,
            edge_colours=edge_colours[0],
            edge_sizes=edge_size_dict
        )
        return net

    def plot_interactive_network(self,
                                 file: str = "LipidNetwork.html",
                                 edge_colour_attr: str = None,
                                 overwrite_edge_attr: bool = True,
                                 edge_cmap: str = "tab10",
                                 edge_size_func: Union[Callable, int] = None,
                                 edge_size_scale: Tuple[int, int] = (.5, 2),
                                 node_colour_attr: str = None,
                                 overwrite_node_attr: bool = True,
                                 node_cmap: str = "tab10",
                                 node_size_func: Union[Callable, int, dict] = None,
                                 node_size_scale: Tuple[int, int] = (10, 40),
                                 edge_group_subset: Union[tuple, list, set] = None,
                                 node_group_subset: Union[tuple, list, set] = None,
                                 edge_size_legend: bool = True,
                                 node_size_legend: bool = True,
                                 save: bool = False,
                                 show: bool = False,
                                 do_return: bool = False,
                                 notebook: bool = False,
                                 net: Network = None,
                                 show_filters: Union[str, List[str]] = None,
                                 as_directed: bool = False,
                                 physics_options: dict = None,
                                 vmin: float = None, vmax: float = None,
                                 force_hierachical: bool = False,
                                 **kwargs) -> Union[Network, None]:
        """
        Plotting the inferred network interactively using pyvis/vis.js

        Parameters
        ----------
        file: str, optional, default 'LipidNetwork.html'
            Path to which the html file (self-contained) is saved to.

        edge_colour_attr: str, optional
            Attribute to use for edge colours.
            Available options: correlations, partial_correlations,
            reaction_types, correlation_changes, partial_correlation_changes.
            All metricx except reaction_types have to be computed beforehand.

        overwrite_edge_attr: bool, optional, default True
            Whether to overwrite attribute with the given group setting in
            self.interaction_attributes

        edge_cmap: str, optional, default 'tab10'
            Colourmap to use for edge colouring. Default for continuous
            attributes is 'seismic'

        edge_size_func: Callable or int, optional
            Used to compute edge sizes if a function.
            Else set as constant edge size.

        edge_size_scale: tuple, optional, default (3, 12)
            Tuple of integers specifying the minimum and maximum
            of plotted edge sizes. The actual values will be scaled
            to fit this range.

        node_colour_attr: str, optional
            Specifies the attribute to use for determining node colours.
            Available options are: fold_changes, pvalues, nlog_pvalues,
            lipid_class, c_index, db_index

        overwrite_node_attr: bool, optional, default True
            Whether to overwrite attribute with the given group setting in
            self.lipid_attributes

        node_size_func: Callable or int, optional
            Used to compute node sizes if a function.
            Else set as constant node size.

        node_size_scale: tuple, optional, default (3, 12)
            Tuple of integers specifying the minimum and maximum
            of plotted node sizes. The actual values will be scaled
            to fit this range.

        node_cmap: str, optional, default 'tab10'
            Colourmap to use for edge colouring. Default for continuous
            attributes is 'seismic'

        edge_group_subset: tuple or list of size 2, optional
            Group-comparison to show in edge colouring/sizing

        node_group_subset: tuple or list of size 2, optional
            Group-comparison to show in node colouring/sizing

        node_size_legend: bool, optional, default True
            Whether to show node sizes in legend

        edge_size_legend: bool, optional, default True
            Whether to show edge sizes in legend

        save: bool, optional, default False,
            Whether to (only) save the plot

        show: bool, optional, default False,
            Whether to open the saved plot in the standard browser

        do_return: bool, optional, default False
            Whether to return the pyvis.network.Network object

        notebook: bool, optional, default False
            Whether the function is called from within a jupyter
            notebook (required to show interactive plot inline)

        net: pyvis.network.Network, optional
            If given this instance will be used instead of generating
            a new one within the function

        show_filters: str or list, optional
            Specifies one or multiple filters to activate. See pyvis
            documentation for available filters.

        as_directed: bool, optional, default False
            Whether to show a direted graph

        vmin: float, optional
            Minimum value of the legend (for continuous scales).
            If not given it will be set to 0 (for only positive
            arrays) or to -max(abs(min), abs(max)) (for zero-centred
            colours).
            Setting vmin can lead to un-centred colour scales!

        vmax: float, optional
            Maximum value of the legend (for continuous scales).
            If vmin is not specified, the minimum value will be
            either 0 (if all values are positive) or -vmax

        physics_options: dict, optional
            Options for physics to use in visualisation. See
            visjs-network and pyvis documentation for available
            specifications.

        force_hierarchical : bool, optional, default False
            Whether to use hierarchical repulsion as a solver
            for network layout. Always True when the number
            of edges exceeds 1200 for performance reasons.

        kwargs:
            keyword arguments to pass to pyvis.network.Network. Only
            releveant if net=None (default)

        Returns
        -------
        pyvis.network.Network or None
        """
        net = self.network_to_pyvis(
            edge_colour_attr=edge_colour_attr,
            overwrite_edge_attr=overwrite_edge_attr,
            edge_cmap=edge_cmap,
            edge_size_func=edge_size_func,
            edge_size_scale=edge_size_scale,
            node_cmap=node_cmap,
            node_size_func=node_size_func,
            node_colour_attr=node_colour_attr,
            overwrite_node_attr=overwrite_node_attr,
            node_size_scale=node_size_scale,
            edge_group_subset=edge_group_subset,
            node_group_subset=node_group_subset,
            node_size_legend=node_size_legend,
            edge_size_legend=edge_size_legend,
            net=net,
            as_directed=as_directed,
            vmax=vmax, vmin=vmin,
            notebook=notebook,
            **kwargs
        )
        if as_directed:
            net.set_edge_smooth("straight_cross")
        if show_filters is not None:
            if isinstance(show_filters, str):
                if show_filters == "all":
                    net.show_buttons(filter_=True)
                else:
                    net.show_buttons(filter_=[show_filters])
            else:
                net.show_buttons(filter_=show_filters)
        if do_return:
            return net
        if show:
            net.show(
                file,
                physics_options=physics_options,
                force_hierachical=force_hierachical
            )
            return
        if save:
            net.write_vis(
                file, physics_options=physics_options,
                force_hierarchical=force_hierachical
            )

    def save_network(self, file: str, file_format: str,
                     weight: str = None, **kwargs):
        """
        Saving the lipid newtork to a file

        Parameters
        ----------
        file: str
            Path to the file in which the network should be saved

        file_format: str
            File format to use. Available options are:
            graphml, yaml, edgelist, adjacency

        weight: str, optional
            Attribute from self.interaction_attributes to use
            as edge weight

        kwargs:
            keyword arguments to pass to pd.DataFrame.to_csv (adjacency)
            or the respective networkx.write_* function
        """
        # use nx.write_graphml etc.
        formats = {"graphml": nx.write_graphml,
                   "yaml": nx.write_yaml,
                   "edgelist": nx.write_weighted_edgelist,
                   "adjacency": None}

        if file_format not in formats.keys():
            raise ValueError(
                f"{file_format} format is currently not supported "
                f"Please choose one of {list(formats.keys())}"
            )
        elif file_format == "adjacency":
            adj = nx.adjacency_matrix(self.network, weight=weight)
            adj_df = pd.DataFrame(adj.toarray(),
                                  columns=list(self.network.nodes),
                                  index=list(self.network.nodes))
            adj_df.to_csv(file, sep=kwargs.pop("sep", "\t"),
                          **kwargs)
        else:
            formats[file_format](self.network, file, **kwargs)

    def _correlation_computation_(self, attr: str,
                                  data: Union[pd.DataFrame, str] = None,
                                  groups: Union[pd.Series, str] = None,
                                  method: str = "pearson",
                                  estimator: str = "LedoitWolf",
                                  significance: float = 0.05,
                                  overwrite: bool = False,
                                  correct_args: dict = None,
                                  **kwargs):
        self._set_data_groups_(data=data, groups=groups)
        if self.groups is None:
            if not overwrite:
                if self.interaction_attributes.get(attr) is not None:
                    return
            if attr == "correlations":
                cors, pvals = correlations(self.data.T,
                                           method=method,
                                           **kwargs)
            else:
                cors, pvals = partial_correlations(self.data, estimator,
                                                   **kwargs)
            # multiple test correction
            if not (pvals < 0).all():
                if correct_args is None:
                    corr_pvals = _matrix_pval_correction_(pvals)
                else:
                    corr_pvals = _matrix_pval_correction_(pvals, **correct_args)
                sig_mask = corr_pvals >= significance
                cors.values[sig_mask] = 0
            else:
                corr_pvals = pvals
            self.interaction_attributes[attr] = cors
            self.interaction_attributes[f"{attr}_pvals"] = corr_pvals
        else:
            if self.interaction_attributes.get(attr) is None or overwrite:
                self.interaction_attributes[attr] = {}
                self.interaction_attributes[f"{attr}_pvals"] = {}
                group_correlations = self.unique_groups.size * [False]
            else:
                group_correlations = [self.interaction_attributes[attr].get(group)
                                      is not None for group in self.unique_groups]
                if all(group_correlations):
                    return

            for group in self.unique_groups[np.invert(group_correlations)]:
                group_data = self.data.loc[self.groups == group, :]
                if attr == "correlations":
                    cors, pvals = correlations(group_data.T,
                                               method=method,
                                               **kwargs)
                else:
                    cors, pvals = partial_correlations(group_data, estimator,
                                                       **kwargs)
                # multiple test correction
                if not (pvals < 0).all():
                    if correct_args is None:
                        corr_pvals = _matrix_pval_correction_(pvals)
                    else:
                        corr_pvals = _matrix_pval_correction_(pvals, **correct_args)
                    sig_mask = corr_pvals >= significance
                    cors.values[sig_mask] = 0
                    cors = cors.fillna(0)
                else:
                    corr_pvals = pvals
                self.interaction_attributes[attr][group] = cors
                self.interaction_attributes[f"{attr}_pvals"][group] = corr_pvals

    def compute_correlations(self, data: Union[pd.DataFrame, str] = None,
                             groups: Union[pd.Series, str] = None,
                             method: str = "pearson",
                             significance: float = 0.05,
                             overwrite: bool = False,
                             correct_args: dict = None,
                             **kwargs):
        """
        Computing correlations between all node pairs

        Parameters
        ----------
        data: pd.DataFrame or str, optional
            Used to overwrite self.data if given. Required
            if no lipid data was added before. If str it
            assumed to be a csv file with samples in columns.

        groups: pd.Series or str, optional
            Used to overwrite self.groups if given. Required
            if no groups were added before. If str it
            assumed to be a csv or txt file.

        method: str, optional, default 'pearson'
            Method to use for correlation calculation. Available
            options are: pearson, kendall and spearman

        significance: float, optional, default 0.05
            Significance threshold for corrected pvalues

        overwrite: bool, optional, default False
            If True previously computed correlations are overwritten

        correct_args: dict, optional
            Arguments to pass to statsmodels.stats.multitest.multipletests
            for multiple test correction

        kwargs:
            keyword arguments to pass to the respective scipy.stats method

        """
        self._correlation_computation_(attr="correlations",
                                       data=data, groups=groups,
                                       method=method, significance=significance,
                                       overwrite=overwrite, correct_args=correct_args,
                                       **kwargs)

    def compute_partial_correlations(self, data: Union[pd.DataFrame, str] = None,
                                     groups: Union[pd.Series, str] = None,
                                     estimator: str = "GraphLasso",
                                     significance: float = 0.05,
                                     overwrite: bool = False,
                                     correct_args: dict = None,
                                     **kwargs):
        """
        Computing correlations between all node pairs

        Parameters
        ----------
        data: pd.DataFrame or str, optional
            Used to overwrite self.data if given. Required
            if no lipid data was added before. If str it
            assumed to be a csv file with samples in columns.

        groups: pd.Series or str, optional
            Used to overwrite self.groups if given. Required
            if no groups were added before. If str it
            assumed to be a csv or txt file.

        estimator: str, optional, default 'GraphLasso'
            Method to use for partial correlation estimation. Available
            options are: GraphLasso, LedoitWolf and empirical

        significance: float, optional, default 0.05
            Significance threshold for corrected p-values. Only relevant
            if sample to feature ratio is large enough to compute
            Fisher's z-transform

        overwrite: bool, optional, default False
            If True previously computed correlations are overwritten

        correct_args: dict, optional
            Arguments to pass to statsmodels.stats.multitest.multipletests
            for multiple test correction

        kwargs:
            keyword arguments to pass to the respective scipy.stats method


        Raises
        ------
        ValueError
            If groups or data are not of the correct type


        Warns
        -----
        If 'data' or 'groups' are already set but provided again
        """
        self._correlation_computation_(attr="partial_correlations",
                                       data=data, groups=groups,
                                       estimator=estimator, significance=significance,
                                       overwrite=overwrite, correct_args=correct_args,
                                       **kwargs)

    def compute_correlation_changes(self, partial_corrs: bool = False,
                                    data: pd.DataFrame = None,
                                    groups: pd.Series = None,
                                    pval_change: bool = True,
                                    overwrite: bool = True):
        """
        Computing changes of partial correlations between different conditions

        Parameters
        ----------
        partial_corrs: bool, optional, default False
            Whether to use partial correlations

        data: pd.DataFrame or str, optional
            Used to overwrite self.data if given. Required
            if no lipid data was added before. If str it
            assumed to be a csv file with samples in columns.

        groups: pd.Series or str, optional
            Used to overwrite self.groups if given. Required
            if no groups were added before. If str it
            assumed to be a csv or txt file.

        pval_change: bool, optional, default True
            Whether to compute categories or absolute values

        overwrite: bool, optional, default False
            If True previously computed correlations are overwritten
        """
        if partial_corrs:
            attr = "partial_correlations"
            attr_change = "partial_correlation_changes"
        else:
            attr = "correlations"
            attr_change = "correlation_changes"
        if self.interaction_attributes[attr] is None:
            raise ValueError(
                f"{attr} have not been computed yet!"
            )
        if self.interaction_attributes[attr_change] is None:
            self.interaction_attributes[attr_change] = {}

        self._set_data_groups_(data=data, groups=groups)
        # compute how correlations/partial correlations change between groups
        # => should be mostly referring to from significant to unsignificant
        # or vice-versa
        if self.groups is None:
            raise ValueError("correlation changes cannot be computed when no groups are given")
        else:
            # TODO: which metrics to use
            for comb in self.comparisons:
                if not overwrite:
                    if self.interaction_attributes[attr_change].get(tuple(comb)) is not None:
                        continue
                    elif self.interaction_attributes[attr_change].get(tuple(comb[0], comb[1])) is not None:
                        continue
                if pval_change:
                    # classes:
                    # both significant
                    #     positive to negative
                    #     negative to positive
                    # significant to not significant
                    # both not significant
                    if (self.interaction_attributes[f"{attr}_pvals"][comb[0]] < 0).all() or \
                        (self.interaction_attributes[f"{attr}_pvals"][comb[1]] < 0).all():
                        raise SignificanceTestError(
                            "Too few samples for the given number of features to compute a fisher's z-transform."
                        )
                    corr_changes = correlation_change(
                        self.interaction_attributes[attr][comb[0]],
                        self.interaction_attributes[attr][comb[1]]
                    )
                else:
                    corr_changes = self.interaction_attributes[comb[0]] - \
                                   self.interaction_attributes[attr][comb[1]]
                self.interaction_attributes[attr_change][tuple(comb)] = corr_changes
            # raise NotImplementedError("correlation comparison currently not implemented")
            # group_correlations = [self.interaction_attributes[group].get("correlation_changes")
            #                       is not None for group in self.unique_groups]
            # if all(group_correlations):
            #     return
            # else:
            #     pass

    def compute_fold_changes(self, data: Union[pd.DataFrame, str] = None,
                             groups: Union[pd.Series, str] = None,
                             compare_groups: Iterable[str] = None,
                             data_is_log: bool = True,
                             round_vals: bool = False,
                             log_func=np.log2,
                             to_log: bool = False):
        """
        Computing fold changes between groups

        Parameters
        ----------
        data: pd.DataFrame or str, optional
            Used to overwrite self.data if given. Required
            if no lipid data was added before. If str it
            assumed to be a csv file with samples in columns.

        groups: pd.Series or str, optional
            Used to overwrite self.groups if given. Required
            if no groups were added before. If str it
            assumed to be a csv or txt file.

        compare_groups: Iterable of str, optional
            Which group comparison to compute. If more than two
            unique groups all pairwise comparisons will be
            computed.

        data_is_log: bool, optional, default True
            Whether data is already log transformed

        round_vals: bool, optional, default False
            If True fold-changes will be rounded to three digits

        log_func: function, optional, default np.log2
            Log function to use if data_is_log is False

        to_log: bool, optional, default False
            Whether to log-transform if data_is_log is false
        """
        self._set_data_groups_(data=data, groups=groups)
        if self.groups is None:
            raise ValueError("fold changes cannot be computed when no groups are given")
        if compare_groups is not None:
            # computing all pairwise combinations of compare_groups
            # (without replacement)
            comparisons = [comb for comb in self.comparisons
                           if comb[0] in compare_groups and comb[1] in compare_groups]
        else:
            comparisons = self.comparisons
        # computing actual fold changes
        if self.lipid_attributes.get("fold_changes") is None:
            fcs = {tuple(comparison): fold_changes(self.data.T, self.groups,
                                                   comparison, data_is_log,
                                                   log_func, to_log=to_log)
                   for comparison in comparisons}
            self.lipid_attributes["fold_changes"] = fcs
        else:
            for comparison in comparisons:
                # comp_st = "_".join(comparison)
                comp_st = tuple(comparison)
                if self.lipid_attributes["fold_changes"].get(comp_st) is None:
                    fc = fold_changes(self.data.T, self.groups,
                                      comparison, data_is_log,
                                      log_func, to_log=to_log)
                    self.lipid_attributes["fold_changes"][comp_st] = fc

    def compute_pvalues(self, data: Union[pd.DataFrame, str] = None,
                        groups: Union[pd.Series, str] = None,
                        compare_groups: Iterable[str] = None,
                        round_vals: bool = False,
                        method: str = "ttest"):
        """
        Computing p-values

        Parameters
        ----------
        data: pd.DataFrame or str, optional
            Used to overwrite self.data if given. Required
            if no lipid data was added before. If str it
            assumed to be a csv file with samples in columns.

        groups: pd.Series or str, optional
            Used to overwrite self.groups if given. Required
            if no groups were added before. If str it
            assumed to be a csv or txt file.

        compare_groups: Iterable of str, optional
            Which group comparison to compute. If more than two
            unique groups all pairwise comparisons will be
            computed.

        round_vals: bool, optional, default False
            If True pvalues will be rounded to four digits

        method: str, optional, default 'ttest'
            Which test to use. Available options are:
            ttest and wilcoxon
        """
        # NOTE: this is only meant for binary comparisons!
        # NOTE: this returns corrected p-values!
        self._set_data_groups_(data=data, groups=groups)
        if self.groups is None:
            raise ValueError("p-values cannot be computed when no groups are given")
        if compare_groups is not None:
            # computing all pairwise combinations of compare_groups
            # (without replacement)
            comparisons = [comb for comb in self.comparisons
                           if comb[0] in compare_groups and comb[1] in compare_groups]
        else:
            comparisons = self.comparisons
        # computing actual p-values
        if self.lipid_attributes.get("pvalues") is None:
            # pvals = {"_".join(comparison): binary_test(self.data.T, self.groups,
            #                                            comparison)
            #          for comparison in comparisons}
            pvals = {tuple(comparison): binary_test(self.data.T, self.groups,
                                                    comparison, method=method)
                     for comparison in comparisons}
            nlog_pvals = {(comp[0], comp[1]): -_pandas_log_(parr)
                          for comp, parr in pvals.items()}

            self.lipid_attributes["pvalues"] = pvals
            self.lipid_attributes["nlog_pvalues"] = nlog_pvals
        else:
            for comparison in comparisons:
                # comp_st = "_".join(comparison)
                comp_st = tuple(comparison)
                if self.lipid_attributes["pvalues"].get(comp_st) is None:
                    pval = binary_test(self.data.T, self.groups,
                                       comparison)
                    self.lipid_attributes["pvalues"][comp_st] = pval
                    self.lipid_attributes["nlog_pvalues"][comp_st] = -_pandas_log_(pval)

    def _get_attr_abs_max_(self, attr: str, nodes: bool) -> float:
        attributes = self.lipid_attributes if nodes else self.interaction_attributes
        if attr not in attributes.keys():
            raise KeyError(
                f'{attr} not found in computed attribute list!'
            )
        if not isinstance(attributes[attr], dict):
            raise ValueError(
                f'{attr} have not been computed in a group-specific manner'
            )
        max_ = 0.
        for group, values in attributes[attr].items():
            cmax = np.nanmax(abs(values.values),)
            if cmax > max_:
                max_ = cmax
        return max_

    def add_network_colours(
        self, attributes: List[str],
        nodes: bool = True,
        cmaps: Dict[str, str] = None
    ) -> Dict[str, Dict[str, Any]]:
        max_ = None
        if nodes:
            generator = self._generate_node_colours
            setter = nx.set_node_attributes
            keys = self.network.nodes
        else:
            generator = self._generate_edge_colours_
            setter = nx.set_edge_attributes
            keys = self.network.edges
        if cmaps is None:
            cmaps = {"lipid_class": "tab20"}
        to_string = {}
        legends_maps = {}
        for attr in attributes:
            if nodes:
                self._check_attributes_(None, attr)
            else:
                self._check_attributes_(attr, None)
            if self._attribute_group_size_[attr] == 0 or \
                self.unique_groups.size == 0:
                self._add_network_attribute_(
                    attr, nodes=nodes, group_subset=None,
                    overwrite=True
                )
                if attr in to_string:
                    pass
                attr_colours = generator(
                    attr, self._discrete_map_[attr],
                    cmaps.get(attr, "tab10"),
                    colours_to_hex=True
                )
                setter(
                    self.network,
                    dict(zip(keys, attr_colours[1])),
                    name=f"{attr}_colour"
                )
                legends_maps[f"{attr}_colour"] = attr_colours[0]
            elif self._attribute_group_size_[attr] == 1:
                if not self._discrete_map_[attr]:
                    max_ = self._get_attr_abs_max_(attr, nodes)
                for group in self.unique_groups:
                    self._add_network_attribute_(
                        attr, nodes=nodes, group_subset=[group],
                        overwrite=True, add_group=True
                    )
                    if not self._discrete_map_[attr]:
                        for attr_name, vmax in zip(
                                [f"{attr}_colour_{group}_individual", f"{attr}_colour_{group}_common"],
                                [None, max_]):
                            attr_colours = generator(
                                attr,
                                self._discrete_map_[attr],
                                cmaps.get(attr, "tab10"),
                                group=group,
                                colours_to_hex=True,
                                group_attr=f"{attr}_{group}",
                                vmax=vmax
                            )
                            setter(
                                self.network,
                                dict(zip(keys, attr_colours[1])),
                                name=attr_name
                            )
                            legends_maps[attr_name] = attr_colours[0]
                    else:
                        attr_colours = generator(
                            attr,
                            self._discrete_map_[attr],
                            cmaps.get(attr, "tab10"),
                            group=group,
                            colours_to_hex=True,
                            group_attr=f"{attr}_{group}"
                        )
                        attr_name = f"{attr}_colour_{group}"
                        setter(
                            self.network,
                            dict(zip(keys, attr_colours[1])),
                            name=attr_name
                        )
                        legends_maps[attr_name] = attr_colours[0]
            else:
                if not self._discrete_map_[attr]:
                    max_ = self._get_attr_abs_max_(attr, nodes)
                for comb in self.comparisons:
                    if not self._discrete_map_[attr]:
                        for attr_name, vmax in zip(
                                [f"{attr}_colour_{_tuple_to_string_(comb)}_individual",
                                 f"{attr}_colour_{_tuple_to_string_(comb)}_common"],
                                [None, max_]):
                            attr_colours = generator(
                                attr,
                                self._discrete_map_[attr],
                                cmaps.get(attr, "tab10"),
                                group=comb,
                                colours_to_hex=True,
                                group_attr=f"{attr}_{_tuple_to_string_(comb)}",
                                vmax=vmax
                            )
                            setter(
                                self.network,
                                dict(zip(keys, attr_colours[1])),
                                name=attr_name
                            )
                            legends_maps[attr_name] = attr_colours[0]
                    else:
                        self._add_network_attribute_(
                            attr, nodes=nodes, group_subset=comb,
                            overwrite=True, add_group=True,
                        )
                        attr_colours = generator(
                            attr,
                            self._discrete_map_[attr],
                            cmaps.get(attr, "tab10"),
                            group=comb,
                            colours_to_hex=True,
                            group_attr=f"{attr}_{_tuple_to_string_(comb)}"
                        )
                        attr_name = f"{attr}_colour_{_tuple_to_string_(comb)}"
                        setter(
                            self.network,
                            dict(zip(keys, attr_colours[1])),
                            name=attr_name
                        )
                        legends_maps[attr_name] = attr_colours[0]
        return legends_maps

    def add_node_sizes(
        self, attributes: List[str],
        scale: Tuple[float, float]
    ) -> Dict[str, Dict[str, Tuple[float, float]]]:
        legends_maps = {}
        abs_ = False
        max_ = None
        for attr in attributes:
            if attr == "fold_changes":
                abs_ = False
            self._check_attributes_(None, attr)
            if self._attribute_group_size_[attr] == 0 or \
                self.unique_groups.size == 0:
                sizes, na_keys = _get_sizes_report_nas_(self.lipid_attributes[attr],
                                                        self.network.nodes)
                attr_size = self._scale_dict_(sizes, scale,
                                              map_as_ex=True)
                if na_keys:
                    for node, val in attr_size[0].items():
                        if np.isnan(val):
                            attr_size[0][node] = scale[0] / 2
                nx.set_node_attributes(
                    self.network, attr_size[0],
                    name=f"{attr}_size"
                )
                legends_maps[f"{attr}_size"] = attr_size[1]
            elif self._attribute_group_size_[attr] == 1:
                if not self._discrete_map_[attr]:
                    max_ = self._get_attr_abs_max_(attr, nodes=True)
                # TODO:
                for group in self.unique_groups:
                    sizes, na_keys = _get_sizes_report_nas_(self.lipid_attributes[attr][group],
                                                            self.network.nodes)
                    if not self._discrete_map_[attr]:
                        for attr_name, vmax in zip(
                                [f"{attr}_size_{group}_individual", f"{attr}_size_{group}_common"],
                                [None, max_]):
                            attr_size = self._scale_dict_(sizes, scale,
                                                          map_as_ex=True,
                                                          vmax=vmax)
                            if na_keys:
                                for node, val in attr_size[0].items():
                                    if np.isnan(val):
                                        attr_size[0][node] = scale[0] / 2
                            nx.set_node_attributes(
                                self.network, attr_size[0],
                                name=attr_name
                            )
                            legends_maps[attr_name] = attr_size[1]
                    else:
                        attr_size = self._scale_dict_(sizes, scale,
                                                      map_as_ex=True)
                        if na_keys:
                            for node, val in attr_size[0].items():
                                if np.isnan(val):
                                    attr_size[0][node] = scale[0] / 2
                        nx.set_node_attributes(
                            self.network, attr_size[0],
                            name=f"{attr}_size_{group}"
                        )
                        legends_maps[f"{attr}_size_{group}"] = attr_size[1]
            else:
                if not self._discrete_map_[attr]:
                    max_ = self._get_attr_abs_max_(attr, nodes=True)
                # TODO: scales
                for comb in self.comparisons:
                    if self.lipid_attributes[attr].get(comb) is None:
                        comb = (comb[1], comb[0])
                    if attr == "fold_changes":
                        sizes, na_keys = _get_sizes_report_nas_(abs(self.lipid_attributes[attr][comb]),
                                                                self.network.nodes)
                    else:
                        sizes, na_keys = _get_sizes_report_nas_(self.lipid_attributes[attr][comb],
                                                                self.network.nodes)
                    if not self._discrete_map_[attr]:
                        for attr_name, vmax in zip(
                                [f"{attr}_size_{_tuple_to_string_(comb)}_individual",
                                 f"{attr}_size_{_tuple_to_string_(comb)}_common"],
                                [None, max_]):
                            attr_size = self._scale_dict_(sizes, scale,
                                                          map_as_ex=True,
                                                          vmax=vmax)
                            if na_keys:
                                for node, val in attr_size[0].items():
                                    if np.isnan(val):
                                        attr_size[0][node] = scale[0] / 2
                            nx.set_node_attributes(
                                self.network, attr_size[0],
                                name=attr_name
                            )
                            legends_maps[attr_name] = attr_size[1]
                    else:
                        attr_size = self._scale_dict_(sizes, scale,
                                                      map_as_ex=True)
                        if na_keys:
                            for node, val in attr_size[0].items():
                                if np.isnan(val):
                                    attr_size[0][node] = scale[0] / 2
                        nx.set_node_attributes(
                            self.network, attr_size[0],
                            name=f"{attr}_size_{_tuple_to_string_(comb)}"
                        )
                        legends_maps[f"{attr}_size_{_tuple_to_string_(comb)}"] = attr_size[1]
        return legends_maps

    # These are utilities for plotting into a dynamic vis file where colours
    # can be selected directly from within the html file -- also used for the
    # Web application
    def dynamic_network(
        self, node_colour_attributes: List[str],
        edge_colour_attributes: List[str],
        node_size_attributes: List[str],
        node_cmaps: Dict[str, str] = None,
        edge_cmaps: Dict[str, str] = None,
        node_size_scale: Tuple[int, int] = (10, 40),
        **kwargs
    ) -> DynamicVisParser:
        # adding static properties
        for node_prop in STATIC_NODE_PROPERTIES:
            if node_prop not in node_colour_attributes:
                node_colour_attributes.append(node_prop)
            if node_prop != "lipid_class":
                if node_prop not in node_size_attributes:
                    node_size_attributes.append(node_prop)
        for edge_prop in STATIC_EDGE_PROPERTIES:
            # reaction enzymes should not get any colours!
            if edge_prop != "reaction_enzymes":
                if edge_prop not in edge_colour_attributes:
                    edge_colour_attributes.append(edge_prop)
        # generating colours and adding to network attributes
        node_colour_legends = self.add_network_colours(node_colour_attributes,
                                                       nodes=True, cmaps=node_cmaps)
        edge_colour_legends = self.add_network_colours(edge_colour_attributes,
                                                       nodes=False, cmaps=edge_cmaps)
        node_size_legends = self.add_node_sizes(node_size_attributes, node_size_scale)

        dvp = DynamicVisParser(**kwargs)
        dvp.from_nx(self.network)
        if dvp.directed:
            dvp.set_edge_smooth("dynamic")
        dvp.generate_legend(
            node_colours=node_colour_legends,
            edge_colours=edge_colour_legends,
            node_sizes=node_size_legends,
            colours_to_hex=True
        )
        return dvp
