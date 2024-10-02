from ..Lipids import (
    LipidSpecies,
    SumSpecies,
    MolecularSpecies,
    snSpecies,
    FattyAcid
)
from typing import Union, List, Callable, Dict, Tuple
from itertools import combinations
import pandas as pd
import numpy as np
from matplotlib.pyplot import Line2D


def _collect_fatty_acids_(lipids: List[MolecularSpecies]) -> List[FattyAcid]:
    fas = set()
    for lipid in lipids:
        for fa in lipid.get_fatty_acids():
            fas.add(fa)
    return list(fas)


# TODO: test
def _infer_fatty_acids_(lipids: Union[List[SumSpecies], List[MolecularSpecies]],
                        level: str = "all") -> List[FattyAcid]:
    # TODO: this has to be changed once a flexible lipid level is possible
    # TODO: add the posibility to have class-specific fas
    if isinstance(lipids[0], MolecularSpecies):
        return _collect_fatty_acids_(lipids)
    # extracting the list of known fatty acids (i.e. lipids with one fatty acid)
    n1 = [species.sum_composition for species in lipids if species.nFA == 1]
    # TODO: add different levels
    # * only use directly observed FAs (implemented)
    # * infer based on same head group e.g. LPC and PC combined
    # * infer based on all classes (implemented)
    if level == "direct":
        return n1
    if level == "all":
        n2 = set()
        for species in lipids:
            if species.nFA != 2:
                continue
            for fa1 in n1:
                fa2 = species.sum_composition - fa1
                if fa2 is not None:
                    n2.add(fa2)
        return list(n2.union(n1))
    if level == "class":
        raise NotImplementedError("head group wise fatty acid inference"
                                  "is not yet implemented")


def _fa_combinations_(fas: List[FattyAcid], size: int) -> dict:
    # key: sum composition, value: list of molecular compositons
    # naive/brute force: go through all combinations and store
    combs = {}
    for comb in combinations(fas, size):
        combs.setdefault(sum(comb), []).append(comb)
    return combs


# TODO: should infer all possible molecular compositions
# from sum composition
def _sum_to_mol_(lipids: List[SumSpecies],
                 fas: dict) -> List[MolecularSpecies]:
    fa_combinations = {
        2: {classes: _fa_combinations_(c_fas, 2)
            for classes, c_fas in fas.items()},
        3: {classes: _fa_combinations_(c_fas, 3)
            for classes, c_fas in fas.items()}
    }
    mol_species = []
    for lipid in lipids:
        # => molecular and sum species are the same
        # for lipids with one fatty acid
        if lipid.nFA == 1:
            mol_species.append(
                MolecularSpecies._from_data_(
                    lipid_class=lipid.lipid_class,
                    fatty_acids=lipid.sum_composition,
                    class_nfa=1,
                    fa_restrictions=lipid.fa_restrictions,
                    reactions=lipid.reactions
                )
            )
        # getting the correct Fatty Acid list to build from
        class_key = False
        for ck in fas.keys():
            if lipid.lipid_class in ck:
                class_key = ck
        if not class_key:
            tmp_fas = fa_combinations[lipid.nFA]["general"]
        else:
            tmp_fas = fa_combinations[lipid.nFA][class_key]
        for fa_comb in tmp_fas:
            mol_species.append(
                MolecularSpecies._from_data_(
                    lipid_class=lipid.lipid_class,
                    fatty_acids=fa_comb,
                    class_nfa=lipid.nFA,
                    fa_restrictions=lipid.fa_restrictions,
                    reactions=lipid.reactions
                )
            )
    return mol_species


def _size_legend_(extrema: tuple, scale: tuple = None,
                  nodes: bool = True, **kwargs) -> List[Line2D]:
    def size_from_scale(step: int, scale: tuple, extrema: tuple):
        inter = (step - scale[0]) / (scale[1] - scale[0])
        return inter * (extrema[1] - extrema[0])

    if scale is None:
        scale_steps = [abs(extrema[1] - extrema[0]) / 3 * i
                       for i in range(1, 5)]
        sizes = scale_steps
    else:
        scale_steps = [abs(scale[1] - scale[0]) / 3 * i
                       for i in range(1, 5)]
        sizes = [size_from_scale(step, scale, extrema)
                 for step in scale_steps]
    if nodes:
        handles = [Line2D([0], [0],
                          markersize=int(ms),
                          marker="o",
                          markerfacecolor="r",
                          label=round(lab, ndigits=4),
                          **kwargs)
                   for ms, lab in zip(scale_steps, sizes)]
    else:
        handles = [Line2D([0], [0], lw=int(lw),
                          color="r",
                          label=round(lab, ndigits=4),
                          **kwargs)
                   for lw, lab in zip(scale_steps, sizes)]
    return handles


def _range_scale_(
        x: Union[np.ndarray, pd.Series],
        a: Union[int, float], b: Union[int, float],
) -> Union[np.ndarray, pd.Series]:
    min_ = x.min()
    max_ = x.max()
    return ((x - min_) / (max_ - min_) * (b - a)) + a


def _pandas_log_(data: Union[pd.DataFrame, pd.Series],
                 log_fun: Callable = np.log10) -> Union[pd.DataFrame, pd.Series]:
    logged = data.copy()
    log_vals = log_fun(logged.values)
    if isinstance(data, pd.DataFrame):
        return pd.DataFrame(log_vals, index=data.index,
                            columns=data.columns)
    elif isinstance(data, pd.Series):
        return pd.Series(log_vals, index=data.index,
                         name=data.name)
    else:
        raise ValueError(
            "_pandas_log_ is only implemented for pd.Series and pd.DataFrame "
            f"not {type(data).__name__}"
        )


def _reaction_enzyme_annotation(src: str, tgt: str, enzyme: str) -> str:
    return f"{src} -> {tgt}: {enzyme}"


def _aggregate_species_(
        data: Union[pd.DataFrame, pd.Series] = None,
        species_axis: str = "columns"
) -> Union[pd.DataFrame, pd.Series]:
    if isinstance(data, pd.DataFrame):
        species = getattr(data, species_axis)
        dupl_mask = species.duplicated()
        duplicates = species[dupl_mask]
        if species_axis == "columns":
            agg_data = data.copy(deep=True).loc[:, ~dupl_mask]
            for dup_idx in duplicates:
                agg_data.loc[:, dup_idx] = data.loc[:, dup_idx].sum(axis=1)
        elif species_axis == "index":
            agg_data = data.copy(deep=True).loc[~dupl_mask, :]
            for dup_idx in duplicates:
                agg_data.loc[dup_idx, :] = data.loc[dup_idx, :].sum(axis=1)
        else:
            raise ValueError(
                "species_axis must be 'columns' or 'index' "
                f"not {species_axis}"
            )
        return agg_data
    elif isinstance(data, pd.Series):
        return data[~data.duplicated()]
    else:
        raise ValueError(
            "Data must be a pandas DataFrame or Series "
            f"not a {type(data).__name__}"
        )


def _lipid_to_object_(
        lipid: str, lipid_level: str, class_fas: Dict[str, int],
) -> Union[SumSpecies, MolecularSpecies, snSpecies]:
    try:
        if not ("_" in lipid or "/" in lipid) or lipid_level == "sum":
            lipid_ = SumSpecies(lipid, class_fas,
                                names_in_convention=True)
        elif "/" not in lipid or lipid_level == "molecular":
            if lipid_level == "molecular":
                # "converting" sn to molecular species
                lipid = lipid.replace("/", "_")
            lipid_ = MolecularSpecies(lipid, class_fas,
                                      names_in_convention=True)
        else:
            lipid_ = snSpecies(lipid, class_fas,
                               names_in_convention=True)
    except ValueError as ve:
        raise ValueError(
            f"Lipid {lipid} caused the following exception:\n{ve}"
        )
    return lipid_


def _get_enzyme_(transform_type: str, enzymes: dict,
                 src: LipidSpecies, tgt: LipidSpecies):
    if transform_type == "Head Group Modification":
        enz = enzymes[src.lipid_class].get(tgt.lipid_class)
        if enz is not None:
            return _reaction_enzyme_annotation(src.lipid_class,
                                               tgt.lipid_class,
                                               enz)
        return None
    return None


def _tuple_to_string_(tup: Tuple[str, str]) -> str:
    return f"{tup[0]}_{tup[1]}"
