from __future__ import annotations
from .FattyAcid import (
    FattyAcid,
    sum_fatty_acids
)
from ..default_globals import _not_implemented_in_subclass, FaDict
from ..exceptions import (
    NameConversionError, MissingClassError,
    LipidDataError
)

from typing import Union, List, Dict, Iterable
from warnings import warn
import re
import numpy as np

import os

# NOTE: for some reason LipidLynxX changes the
# working directory in ALL of these three cases:
# 1. module import
# 2. Converter object initialisation
# 3. Converter function calls
wdir = os.getcwd()
from lynx import Converter

os.chdir(wdir)


def _ordered_string_(*args, sep: str = " - "):
    ordered_args = sorted(args)
    return sep.join(ordered_args)


def _sum_composition_in_fas_(lclass: str, sum_comp: FattyAcid,
                             fas: FaDict):
    # key_match = np.array([lclass in key for key in fas.keys()])
    # if any(key_match):
    #     # case: the larger lipid class has a specific FA-pool
    #     idx = np.where(key_match)[0][0]
    #     return sum_comp in fas[list(fas.keys())[idx]]
    # # case: no specific FA-pool given
    # return sum_comp in fas["general"]
    return sum_comp in fas.get(lclass)


# Decorator for instance checking specific function calls
def instance_check(n_args, *arg_types):
    def _check_instance(func):
        def wrapper(*args):
            if len(args) - 1 != n_args:
                raise ValueError(
                    #
                    f"Exactly {n_args} argument(s) allowed for {func.__name__}, "
                    f"but {len(args) - 1} provided"
                )
            if n_args > 1:
                for i in range(2, len(args)):
                    if not isinstance(args[i], arg_types[i - 2]):
                        raise ValueError(
                            # args[0] is self
                            f"argument {i} of {func.__name__} must be an instance "
                            f"of {arg_types[i - 2]}, not {type(args[i]).__name__}"
                        )
            return func(*args)

        return wrapper

    return _check_instance


# TODO: implement a nomenclature unifier
# => unify first - then get lipid class and fatty acids
# from a standard format
def unify_name(species_name: Union[str, list],
               n_fas: dict, converter: Converter = None,
               **kwargs) -> Union[str, list]:
    # warn("name unification has not been implemented properly yet")
    if converter is None:
        converter = Converter()

    if isinstance(species_name, Iterable) and not isinstance(species_name, str):
        return [unify_name(species, n_fas, converter=converter)
                for species in species_name]
    try:
        uni = converter.convert(species_name, **kwargs).dict()
    except ValueError:
        if species_name.startswith("P") or species_name.startswith("LP"):
            try:
                # NOTE: for some weir reason LipidLynxX fails for e.g. PC(16:0_16:0) (but converts to it)
                uni = converter.convert(species_name.replace("_", "/"), **kwargs).dict()
            except ValueError or IndexError:
                raise NameConversionError(
                    f"{species_name} is not in an incorrect format! "
                    "See LipidLynxX documentation for supported formats",
                    species_name
                )
        else:
            raise NameConversionError(
                f"{species_name} is not in any supported nomenclature format! "
                "See LipidLynxX documentation for supported formats",
                species_name
            )
    except IndexError:
        if species_name.startswith("P") or species_name.startswith("LP"):
            try:
                # NOTE: for some weir reason LipidLynxX fails for e.g. PC(16:0_16:0) (but converts to it)
                uni = converter.convert(species_name.replace("_", "/"), **kwargs).dict()
            except ValueError or IndexError:
                raise NameConversionError(
                    f"{species_name} is not in an incorrect format! "
                    "See LipidLynxX documentation for supported formats",
                    species_name
                )
        else:
            raise NameConversionError(
                f"{species_name} is not in an incorrect format! Maybe exchanging '/' and '_' resolves the issue. "
                "See LipidLynxX documentation for supported formats",
                species_name
            )
    except Exception as e:
        raise NameConversionError(
            f"An unknow error (see below) occured at {species_name}. "
            f"See LipidLynxX documentation for supported formats\n{e}",
            species_name
        )
    if uni["output"] == "":
        lc = re.findall("[a-zA-Z]*", uni["skipped"])[0]
        if n_fas.get(lc) == 1:
            kwargs.pop("level")
            uni = converter.convert(species_name,
                                    level="B0",
                                    **kwargs).dict()
        else:
            raise NameConversionError(
                f"Lipid {species_name} could not be converted into "
                "LipidLynxX format. Please make sure the nomenclature "
                "requirements are met.",
                species_name
            )
    if "ST(" in uni["output"]:
        return uni["output"].replace("ST", "CE")
    # NOTE: FAs are returned without bracekts => need to add
    # them in here!
    if re.match("^FA[0-9]", uni["output"]):
        return uni["output"].replace("FA", "FA(") + ")"
    return uni["output"]


class LipidSpecies:
    # TODO: revise once SumSpecies and MolecularSpecies are finished
    __slots__ = ("name", "lipid_class",
                 "sum_composition", "nFA",
                 "FA_combinations",
                 "possible_FAs", "fatty_acids",
                 "ether_lipid", "has_long_chain_base")

    def __init__(self, species_name: str,
                 n_fas: dict,
                 names_in_convention: bool = False,
                 **kwargs):
        # TODO: documentation
        # TODO: should name conversion be moved into the network class
        # (for efficiency reasons)?
        if not names_in_convention:
            self.name = unify_name(species_name, n_fas, **kwargs)
            os.chdir(wdir)
        else:
            self.name = species_name
        self.lipid_class = self.get_lipid_class().strip()

        if "O-" in self.name:
            self.ether_lipid = "O"
        elif "P-" in self.name:
            self.ether_lipid = "P"
        else:
            self.ether_lipid = ""

        if self.lipid_class == "SM" or self.lipid_class == "Cer":
            self.has_long_chain_base = True
        else:
            self.has_long_chain_base = False

        self.lipid_class = self.lipid_class + self.ether_lipid

    def __hash__(self):
        return hash(self.name)

    def is_equal(self, species: LipidSpecies) -> bool:
        if not isinstance(species, LipidSpecies):
            raise ValueError("'species' must be an instance of LipidSpecies"
                             " or a subclass")
        return self.class_equal(species) and self.fas_equal(species)

    # only defined for convenience - personally I feel like is_equal is more readable
    # most of the time, but sometimes not very handy
    def __eq__(self, species_name):
        return self.is_equal(species_name)

    def __repr__(self):
        return self.name

    def __str__(self):
        return self.name

    @instance_check(3, dict, FaDict)
    def is_transformable(self, species: LipidSpecies,
                         class_connections: dict,  # str: list[str]
                         fas: FaDict) -> bool:
        if not isinstance(species, LipidSpecies):
            raise ValueError(
                "'species' must be an instance of LipidSpecies"
                f" or a subclass, not {type(species).__name__}"
            )

        # TODO: double check!
        # NOTE: this only works in a correct way, if class_connections
        # is correctly specified in a way, such that transformable
        # classes with different numbers of fatty acids have the same
        # head group (e.g. PC and LPC, DG and TG)
        if self.nFA == species.nFA:
            #
            class_transform = self.class_transformable(species, class_connections) and \
                              self.fas_equal(species)
            fa_transform = self.class_equal(species) and \
                           self.fas_transformable(species, fas)
            # alternative: return class_transform is not fa_transform
            # xor statement
            return class_transform ^ fa_transform
        else:
            # we can skip the class_equal step, because we know it cannot
            # be True here
            # this is the step that implies that two classes with different
            # fatty acid numbers have the same head group
            # fas_transformable should only return True if the missing FA
            # is in fas and for molecular species if additionally other FAs
            # in the respective combination are the same
            return self.class_transformable(species, class_connections) and \
                   self.fas_transformable(species, fas)

    @instance_check(2, dict)
    def class_transformable(self, species: LipidSpecies,
                            class_connections: dict) -> bool:
        try:
            return species.lipid_class in class_connections[self.lipid_class]
        except KeyError:
            raise MissingClassError(
                f"{self.lipid_class} is missing in lipid class settings file. "
                "Please add all possible reactions to the settings file (see tutorial for format).",
                class_=self.lipid_class, lipid=self.name, error_location="Class Connections"
            )

    @instance_check(1)
    def class_equal(self, species: LipidSpecies) -> bool:
        return species.lipid_class == self.lipid_class

    @instance_check(1)
    def transform_type(self, species: LipidSpecies, level: str = "sum") -> Union[str, bool]:
        # NOTE: this function is meant to work on sum compositions only
        # TODO: should level options be renamed to have a less conflicting terminology?
        if level == "sum":
            # NOTE: transform_type does NOT double check whether
            # is actually true, but assumes it instead
            # => only use this function for known pairs!
            if self.lipid_class == species.lipid_class:
                # same class => check which FA transformation
                if level == "sum":
                    return self.sum_composition.transform_type(species.sum_composition)
                return _ordered_string_(self.sum_composition, species.sum_composition)
            if self.nFA == species.nFA:
                # different class, same number of fatty acids
                # => FAs must be equal for is_transformable to be true
                # => head group modification
                if level == "sum":
                    return "Head Group Modification"
                return _ordered_string_(self.lipid_class, species.lipid_class)
            # different number of fatty acids
            # => same head group, different class
            # => addition/removal of one FA
            return "FA addition"

    def get_lipid_class(self) -> str:
        if hasattr(self, "lipid_class"):
            return self.lipid_class
        if self.name.startswith("FA"):
            return re.split("\(|[0-9]", self.name)[0]
        return re.split("\(", self.name)[0]

    # Prototypes of functions to be defined in sub-classes
    def get_fatty_acids(self):
        _not_implemented_in_subclass("get_fatty_acids", "LipidSpecies")

    def fas_equal(self, species: LipidSpecies) -> bool:
        _not_implemented_in_subclass("fas_equal", "LipidSpecies")

    def fas_transformable(self, species: LipidSpecies, fas) -> bool:
        _not_implemented_in_subclass("fas_transformable", "LipidSpecies")

    def to_dict(self):
        _not_implemented_in_subclass("to_dict", "LipidSpecies")


# TODO add fa restrictions
class SumSpecies(LipidSpecies):
    def __init__(self, species_name: str,
                 class_nfa: Dict[str, int],
                 names_in_convention: bool = False,
                 reactions: dict = None,
                 fa_restrictions: dict = None):
        super(SumSpecies, self).__init__(
            species_name,
            n_fas=class_nfa,
            level="B0",
            names_in_convention=names_in_convention
        )

        try:
            self.sum_composition = self.get_fatty_acids()
        except IndexError:
            raise LipidDataError(
                f"{species_name} was not recognised as a lipid. "
                "Did you convert notation conventions?"
            )
        try:
            self.nFA = class_nfa[self.lipid_class]
        except KeyError:
            raise MissingClassError(
                f"{self.lipid_class} not found in lipid settings",
                self.lipid_class, self.name, "Class FA Numbers"
            )
        self.fa_restrictions = fa_restrictions
        self.reactions = reactions

        self.FA_combinations = None
        self.possible_FAs = None

    def compute_fa_combinations(self,
                                fa_combinations: Dict[int, List[FattyAcid]],
                                fa_restrictions: dict = None):
        # TODO: add FA restrictions
        # TODO: include hydroxy-groups => how?
        if self.FA_combinations is not None:
            return self.FA_combinations
        if self.nFA == 1:
            self.FA_combinations = [self.sum_composition]
            self.possible_FAs = self.FA_combinations
        else:
            self.FA_combinations = self._infer_mol_compositions_(fa_combinations,
                                                                 fa_restrictions)
            self.possible_FAs = {fa for fa_comb in self.FA_combinations
                                 for fa in fa_comb}

    def to_dict(self) -> dict:
        return {
            "LipidClass": self.lipid_class,
            "SumComposition": self.sum_composition,
            "FA_Combinations": self.FA_combinations
        }

    def _infer_mol_compositions_(self,
                                 fa_combinations: Dict[int, List[FattyAcid]],
                                 fa_restrictions: dict = None) -> list:
        if fa_restrictions is not None:
            warn("fa_restrictions have currently no effect!")
        matching_combs = [comb for comb in fa_combinations[self.nFA]
                          if sum_fatty_acids(comb) == self.sum_composition]
        return matching_combs

    def get_fatty_acids(self) -> FattyAcid:
        """

        Returns
        -------
        instance of FattyAcid with sum composition properties (c, db etc.)
        """
        if hasattr(self, "sum_composition"):
            return self.sum_composition
        fa_base = re.sub("[A-Za-z]+\(|\)", "", self.name)
        if "_" in fa_base or "/" in fa_base:
            return sum([FattyAcid(fa) for fa in re.split("[_/]", fa_base)])
        return FattyAcid(fa_base)

    @instance_check(1)
    def fas_equal(self, species) -> bool:
        return self.sum_composition == species.sum_composition

    @instance_check(2, FaDict)
    def fas_transformable(self, species,
                          fas: FaDict) -> bool:
        if self.nFA == 1 and species.nFA == 2:
            if isinstance(species, MolecularSpecies) or isinstance(species, snSpecies):
                diff = species.sum_composition - self.sum_composition
                if diff is None:
                    return False
                fa_match = diff in fas.get(species.lipid_class)
                sum_match = self.sum_composition in species.fatty_acids
                return sum_match and fa_match
        # transformable only when the difference is at max 1
        nfa_diff = self.nFA - species.nFA
        if abs(nfa_diff) > 1:
            return False
        # same number of fatty acids for both species
        if nfa_diff == 0:
            return self.sum_composition.is_transformable(species.sum_composition,
                                                         reactions=self.reactions,
                                                         excluded_reactions=self.fa_restrictions)
        # self has one FA more than species
        if nfa_diff > 0:
            return self._unequal_transformable_(self, species, fas)
        # self has one FA more than species
        return self._unequal_transformable_(species, self, fas)

    @staticmethod
    def _unequal_transformable_(l_species, s_species,
                                fas: FaDict) -> bool:
        # TODO: test
        # NOTE: use FA_combinations and essentially do the same as in MolecularSpecies
        sum_diff = l_species.sum_composition - s_species.sum_composition
        # this means there are no FAs fulfilling this criterion
        if sum_diff is None:
            return False
        return _sum_composition_in_fas_(l_species.lipid_class, sum_diff, fas)


class MolecularSpecies(LipidSpecies):
    def __init__(self, species_name: str,
                 class_nfa: Dict[str, int],
                 names_in_convention: bool = False,
                 reactions: dict = None,
                 fa_restrictions: dict = None):
        super(MolecularSpecies, self).__init__(
            species_name,
            n_fas=class_nfa,
            names_in_convention=names_in_convention,
            level="M0"
        )

        try:
            self.fatty_acids = self.get_fatty_acids()
        except IndexError:
            raise LipidDataError(
                f"{species_name} was not recognised as a lipid. "
                "Did you convert notation conventions?"
            )
        try:
            self.nFA = class_nfa[self.lipid_class]
        except KeyError:
            raise MissingClassError(
                f"{self.lipid_class} not found in lipid settings",
                self.lipid_class, self.name, "Class FA Numbers"
            )
        self.sum_composition = self.get_sum_composition()
        self.fa_restrictions = fa_restrictions
        self.reactions = reactions

    @classmethod
    def _from_data_(cls, lipid_class: str,
                    fatty_acids: Union[List[str], List[FattyAcid]],
                    class_nfa: int, fa_restrictions: dict,
                    reactions: dict, fas_as_str: bool = False):
        """
        Utility used when inferring molecular species from sum compositions.
        Not meant to be used outside the package unless the user is familiar,
        with the nitty gritty of this class initialisation and the required
        nomenclauture.

        Parameters
        ----------
        lipid_class
        fatty_acids
        class_nfa
        fa_restrictions
        reactions
        fas_as_str

        Returns
        -------
        MolecularSpecies

        """
        species_name = f"{lipid_class})"
        if fas_as_str:
            for n, fa in enumerate(fatty_acids):
                if n == len(fatty_acids):
                    species_name += f"{fa})"
                else:
                    species_name += f"{fa}_"
        else:
            for n, fa in enumerate(fatty_acids):
                if n == len(fatty_acids):
                    species_name += f"{fa.name})"
                else:
                    species_name += f"{fa.name}_"
        return MolecularSpecies(
            species_name, {lipid_class: class_nfa},
            names_in_convention=True,
            reactions=reactions,
            fa_restrictions=fa_restrictions
        )

    def to_dict(self) -> dict:
        return {
            "LipidClass": self.lipid_class,
            "FAs": self.fatty_acids,
            "SumComposition": self.sum_composition
        }

    def get_fatty_acids(self) -> list:
        if hasattr(self, "fatty_acids"):
            return self.fatty_acids
        # TODO: check if pattern works properly
        fa_base = re.sub("[A-Za-z]+\(|\)", "", self.name)
        fas = fa_base.split("_")
        return sorted([FattyAcid(fa) for fa in fas])

    def get_sum_composition(self) -> FattyAcid:
        if hasattr(self, "sum_composition"):
            return self.sum_composition
        return sum_fatty_acids(self.fatty_acids)

    # TODO: add position specificity for ether-lipids
    # TODO: add long-chain base specificity
    def fas_equal(self, species) -> bool:
        if isinstance(species, SumSpecies):
            return self.sum_composition == species.sum_composition
        # comparing each fatty acid
        # => MolecularSpecies does NOT consider sn-specificity
        #
        # this is a set rule!
        if self.nFA != species.nFA:
            return False
        if self.sum_composition != species.sum_composition:
            return False

        matches = _equal_pos_(self.fatty_acids,
                              species.fatty_acids)
        if matches is None:
            return False
        return len(matches) == self.nFA

    # TODO: add position specificity for ether-lipids?
    # TODO: add long-chain base specificity?
    def _fa_transform_sum_(self, species, fas: FaDict) -> bool:
        return species.fas_transformable(self, fas)

    def _fa_transform_mol_(self, species, fas: FaDict) -> bool:
        # NOTE: fas_transformable does NOT need to be implemented
        # for unequal FAs
        # this is a set rule!
        if abs(self.nFA - species.nFA) > 1:
            return False

        if self.nFA == species.nFA:
            if self.nFA == 1:
                return self.sum_composition.is_transformable(species.sum_composition,
                                                             reactions=self.reactions,
                                                             excluded_reactions=self.fa_restrictions)
            equal_matches = _equal_pos_(self.fatty_acids,
                                        species.fatty_acids)
            if equal_matches is None:
                return False
            target_n = self.nFA - 1
            if len(equal_matches) == target_n:
                # check if all but n-1 fatty acids match
                match_arr = np.array(equal_matches)
                fa1 = np.delete(self.fatty_acids, match_arr[:, 0])
                fa2 = np.delete(species.fatty_acids, match_arr[:, 1])
                if fa1.size != 1 or fa2.size != 1:
                    raise ValueError(
                        f"Something went wrong while matching {self.name} and {species.name}"
                    )
                # check if remaining FA pair is transformable
                return fa1[0].is_transformable(fa2[0],
                                               reactions=self.reactions,
                                               excluded_reactions=self.fa_restrictions)
            return False
        else:
            target_n = min(self.nFA, species.nFA)
            # NOTE: this is always positive (by FA subtraction implementation)
            sum_comp_diff = self.sum_composition - species.sum_composition
            # case: no subtraction possible (i.e. impossible FA result)
            if sum_comp_diff is None:
                return False
            # case: self is the "smaller" lipid
            if self.nFA < species.nFA:
                if not _sum_composition_in_fas_(species.lipid_class,
                                                sum_comp_diff, fas):
                    return False
            # case: self is the "larger" lipid
            elif self.nFA > species.nFA:
                if not _sum_composition_in_fas_(self.lipid_class,
                                                sum_comp_diff, fas):
                    return False
            # arriving here means theoretically the sum composition
            # difference is doable => now check whether individual
            # FAs are aligning
            equal_matches = _equal_pos_(self.fatty_acids,
                                        species.fatty_acids)
            if equal_matches is None:
                return False
            if len(equal_matches) == target_n:
                return True
            return False

    def fas_transformable(self, species: Union[SumSpecies, MolecularSpecies],
                          fas: FaDict) -> bool:
        if isinstance(species, SumSpecies):
            return self._fa_transform_sum_(species, fas)
        return self._fa_transform_mol_(species, fas)

    def transform_type(self, species: LipidSpecies, level="sum") -> Union[str, bool]:
        if level == "sum" or isinstance(species, SumSpecies):
            if self.lipid_class == species.lipid_class:
                return self.sum_composition.transform_type(species.sum_composition)
            if self.nFA == species.nFA:
                return "Head Group Modification"
            return "FA addition"

        if self.lipid_class == species.lipid_class:
            # find which FAs are the transformed pair
            if self.nFA == 1:
                fa_trans = [self.fatty_acids[0],
                            species.fatty_acids[0]]
            else:
                matching_pos = np.array(_equal_pos_(self.fatty_acids,
                                                    species.fatty_acids))
                fa_trans = [
                    np.delete(np.array(self.fatty_acids),
                              matching_pos[:, 0])[0].name,
                    np.delete(np.array(species.fatty_acids),
                              matching_pos[:, 1])[0].name
                ]
            return _ordered_string_(*fa_trans)
        if self.nFA == species.nFA:
            return _ordered_string_(self.lipid_class, species.lipid_class)
        add = self.sum_composition - species.sum_composition
        return f"{add} added"


class snSpecies(LipidSpecies):
    def __init__(self, species_name: str,
                 class_nfa: Dict[str, int],
                 names_in_convention: bool = False,
                 reactions: dict = None,
                 fa_restrictions: dict = None):
        super(snSpecies, self).__init__(
            species_name,
            n_fas=class_nfa,
            names_in_convention=names_in_convention,
            level="S0"
        )

        try:
            self.fatty_acids = self.get_fatty_acids()
        except IndexError:
            raise LipidDataError(
                f"{species_name} was not recognised as a lipid. "
                "Did you convert notation conventions?"
            )
        try:
            self.nFA = class_nfa[self.lipid_class]
        except KeyError:
            raise MissingClassError(
                f"{self.lipid_class} not found in lipid settings",
                self.lipid_class, self.name, "Class FA Numbers"
            )
        self.reactions = reactions
        self.fa_restrictions = fa_restrictions
        self.sum_composition = self.get_sum_composition()

    def to_dict(self) -> dict:
        return {
            "LipidClass": self.lipid_class,
            "FAs": self.fatty_acids,
            "SumComposition": self.sum_composition
        }

    def get_fatty_acids(self) -> list:
        if hasattr(self, "fatty_acids"):
            return self.fatty_acids
        # TODO: check if pattern works properly
        fa_base = re.sub("[A-Za-z]+\(|\)", "", self.name)
        fas = fa_base.split("/")
        return [FattyAcid(fa) for fa in fas]

    def get_sum_composition(self) -> FattyAcid:
        if hasattr(self, "sum_composition"):
            return self.sum_composition
        return sum_fatty_acids(self.fatty_acids)

    def fas_equal(self, species) -> bool:
        if isinstance(species, SumSpecies):
            return self.sum_composition == species.sum_composition
        elif isinstance(species, MolecularSpecies):
            return species.fas_equal(self)
        else:
            if len(self.fatty_acids) != len(species.fatty_acids):
                return False
            for s_fa, t_fa in zip(self.fatty_acids, species.fatty_acids):
                if s_fa != t_fa:
                    return False
            return True

    def fas_transformable(self, species,
                          fas: FaDict) -> bool:
        if not isinstance(species, snSpecies):
            return species.fas_transformable(self, fas)

        n_diff = self.nFA - species.nFA
        # set rule
        if abs(n_diff) > 1:
            return False
        if n_diff < 0:
            return species.fas_transformable(self, fas)
        fas_match = np.array([s_fa == t_fa for s_fa, t_fa in
                              zip(self.fatty_acids, species.fatty_acids)])
        if n_diff == 0:
            if self.nFA == 1:
                return self.fatty_acids[0] == species.fatty_acids[0]
            else:
                if fas_match.sum() != (self.nFA - 1):
                    return False
                # This is the case:
                # same number of fatty acids and all but one match
                # => check if this one is transformable
                trans_idx = np.where(np.invert(fas_match))[0][0]
                return self.fatty_acids[trans_idx].is_transformable(species.fatty_acids[trans_idx],
                                                                    reactions=self.reactions,
                                                                    excluded_reactions=self.fa_restrictions)
        # Case: on of the two species has one more fatty acid
        # => check if the "additional" fatty acid is in the list
        # of allowed fatty acids
        if n_diff < 0:
            # species has more FAs than self
            if fas_match.sum() != self.nFA:
                return False
            return species.fatty_acids[species.nFA] in fas.get(species.lipid_class)
        # self has more FAs than species
        if fas_match.sum() != species.nFA:
            return False
        return self.fatty_acids[self.nFA] in fas.get(self.lipid_class)

    def transform_type(self, species: "LipidSpecies", level: str = "sum") -> Union[str, bool]:
        if level == "sum" or isinstance(species, SumSpecies):
            if self.lipid_class == species.lipid_class:
                return self.sum_composition.transform_type(species.sum_composition)
            if self.nFA == species.nFA:
                return "Head Group Modification"
            return "FA addition"
        if isinstance(species, MolecularSpecies):
            return species.transform_type(self, level=level)
        # we assume here that the two lipid species are known to be transformable
        n_diff = self.nFA - species.nFA
        if n_diff == 0:
            for fa_pair in zip(self.fatty_acids, species.fatty_acids):
                if fa_pair[0] != fa_pair[1]:
                    return _ordered_string_(*fa_pair)
            return _ordered_string_(self.lipid_class, species.lipid_class)
        return f"{self.sum_composition - species.sum_composition} added"


def _equal_pos_(a, b, pairs: list = None,
                matched: list = None) -> Union[
    None, Union[List[List[int]]]
]:
    # TODO: revise
    if pairs is None:
        pairs = []
    if matched is None:
        matched = [[], []]
    for i, ai in enumerate(a):
        if i in matched[0]:
            continue
        for j, bj in enumerate(b):
            if j in matched[1]:
                continue
            if ai == bj:
                pairs.append([i, j])
                matched[0].append(i)
                matched[1].append(j)
                return _equal_pos_(a, b,
                                   pairs=pairs,
                                   matched=matched)
    if len(pairs) == 0:
        return None
    return pairs


if __name__ == "__main__":
    import os
    import pandas as pd

    os.chdir("/home/nik/PhD/LipidNetworkInference/ValidationStudies/Ansari2016NatureSD")

    data = pd.read_csv("lynx_format_data.tsv", index_col=0, sep='\t')
    lipid_names = list(data.columns[35:45])
    print(lipid_names)

    s_lipids = [SumSpecies(lipid, {"DG": 2}, names_in_convention=True) for lipid in lipid_names]
    print(f"{s_lipids[0].lipid_class}\t{s_lipids[0].sum_composition}")

    m_lipids = [MolecularSpecies(lipid, {"DG": 2}, names_in_convention=True) for lipid in lipid_names]
    print(f"{m_lipids[0].lipid_class}\t{[fa.name for fa in m_lipids[0].fatty_acids]}")

    lipid_names = [lipid.replace("_", "/") for lipid in lipid_names]
    sn_lipids = [snSpecies(lipid, {"DG": 2}, names_in_convention=True) for lipid in lipid_names]
    print(f"{sn_lipids[0].lipid_class}\t{[fa.name for fa in sn_lipids[0].fatty_acids]}")
