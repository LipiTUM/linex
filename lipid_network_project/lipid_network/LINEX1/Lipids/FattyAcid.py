from __future__ import annotations
from collections.abc import Iterable
from typing import Union, List
import re


def _pair_to_str_(idx_i: str, idx_j: str) -> str:
    try:
        id_i = int(idx_i)
        id_j = int(idx_j)
    except ValueError:
        raise ValueError("'idx_i' and 'idx_j' must be strings"
                         " of integers")
    primer = "{0} - {1}"
    if id_i < id_j:
        return primer.format(idx_i, idx_j)
    return primer.format(idx_j, idx_i)


class FattyAcid:
    # TODO: add double bond positions
    __slots__ = ("name", "c_index",
                 "db_index", "oh_index")

    def __init__(self, name: str):
        # TODO: documentation
        # if "O-" in name:
        #     self.ether_bond = name.count("O-")
        # elif "P-" in name:
        #     self.ether_bond = name.count("P-") * 2
        # else:
        #     self.ether_bond = 0

        name_ = re.sub("[A-Za-z]*| |-", "", name)
        indices = re.split("[:;<]", name_)

        self.c_index = int(indices[0])
        self.db_index = int(indices[1])
        if len(indices) == 3:
            if indices[2].startswith("{"):
                self.oh_index = int(re.sub("[>{}]", "", indices[2]))
            # this is the case when certain OH notation and 1 OH-group
            elif indices[2] == ">":
                self.oh_index = 1
            elif indices[2].endswith(">"):
                self.oh_index = int(indices[2].split(">")[0])
            else:
                try:
                    self.oh_index = int(indices[2])
                except ValueError:
                    try:
                        self.oh_index = int(indices[2].replace("O", ""))
                    except ValueError as ve:
                        if indices[2].replace("O", "") == "":
                            self.oh_index = 1
                        else:
                            raise ve
        else:
            self.oh_index = 0

        self.name = "{0}:{1};{2}".format(self.c_index, self.db_index,
                                         self.oh_index)

    def is_equal(self, fa: FattyAcid) -> Union[bool, list]:
        # TODO: documentation
        if isinstance(fa, FattyAcid):
            return self.c_index == fa.c_index and \
                   self.db_index == fa.db_index and \
                   self.oh_index == fa.oh_index
        elif isinstance(fa, Iterable):
            return [self.is_equal(fa_i) for fa_i in fa]
        else:
            raise ValueError(
                "fa must be an instance of Fatty Acid or "
                f"an iterable of FattyAcids, not {type(fa).__name__}"
            )

    # only for convenience
    def __str__(self):
        return self.name

    def __repr__(self):
        return f"FA {self.name}"

    def __eq__(self, fa: FattyAcid) -> bool:
        return self.is_equal(fa)

    def __add__(self, fa):
        c = self.c_index + fa.c_index
        db = self.db_index + fa.db_index
        oh = self.oh_index + fa.oh_index
        return FattyAcid("{0}:{1};{2}".format(c, db, oh))

    def __radd__(self, other: Union[FattyAcid, int]):
        if isinstance(other, int):
            if other != 0:
                raise ValueError(
                    "'FattyAcid' object cannnot be added to an integer"
                )
            return self.copy()
        c = self.c_index + other.c_index
        db = self.db_index + other.db_index
        oh = self.oh_index + other.oh_index
        return FattyAcid("{0}:{1};{2}".format(c, db, oh))

    def __sub__(self, fa):
        # subtract longer from shorter FA
        if self.c_index < fa.c_index:
            return fa.__sub__(self)
        # n_DB must be larger or equal
        if self.db_index < fa.db_index:
            # print("Number of double bonds of the smaller fatty acid",
            #       "must be smaller equal the number in the longer chain")
            return None
        # n_OH must be larger or equal
        if self.oh_index < fa.oh_index:
            # print("Number of hydroxy-groups of the smaller fatty acid",
            #       "must be smaller equal the number in the longer chain")
            return None

        c = self.c_index - fa.c_index
        db = self.db_index - fa.db_index
        oh = self.oh_index - fa.oh_index
        return FattyAcid("FA {0}:{1};{2}".format(c, db, oh))

    def __hash__(self):
        # TODO: should this implement a custom annotation?
        return hash((self.c_index, self.db_index, self.oh_index))

    def __lt__(self, fa) -> bool:
        if not isinstance(fa, FattyAcid):
            raise NotImplementedError("less/greater than FattyAcid only implemented"
                                      " for other instances of FattyAcid, not '{0}'".format(type(fa).__name__))
        if self == fa:
            return False
        # sorting by number of c, db, oh in that order
        if self.c_index < fa.c_index:
            return True
        elif self.c_index > fa.c_index:
            return False
        # DBs
        if self.db_index < fa.db_index:
            return True
        elif self.db_index > fa.db_index:
            return False
        # OHs
        if self.oh_index < fa.oh_index:
            return True
        elif self.oh_index > fa.oh_index:
            return False
        return False

    def __deepcopy__(self, memodict: dict = None) -> FattyAcid:
        return FattyAcid(self.name)

    def copy(self):
        return self.__deepcopy__()

    def is_transformable(
            self, fa: FattyAcid,
            reactions: dict = None,
            excluded_reactions: List[list] = None
    ) -> Union[bool, list]:
        # TODO: documentation
        """
        Parameters
        ----------
        fa
        reactions
        excluded_reactions

        Returns
        -------

        """
        if reactions is None:
            # These are the standard reactions
            reactions = {
                "Chain length": {"C": 2, "DB": 0, "OH": 0},
                "Desaturation": {"C": 0, "DB": 1, "OH": 0},
                "Hydroxylation": {"C": 0, "DB": 0, "OH": 1}
            }
        if isinstance(fa, FattyAcid):
            if excluded_reactions is not None:
                if sorted([self.name, fa.name]) in excluded_reactions:
                    return False
            c_transform = abs(self.c_index - fa.c_index)
            db_transform = abs(self.db_index - fa.db_index)
            oh_transform = abs(self.oh_index - fa.oh_index)
            for reaction in reactions.values():
                c = c_transform == reaction["C"]
                db = db_transform == reaction["DB"]
                oh = oh_transform == reaction["OH"]
                if c and db and oh:
                    return True
            return False
        elif isinstance(fa, Iterable):
            return [self.is_transformable(fa_i, reactions,
                                          excluded_reactions)
                    for fa_i in fa]
        else:
            raise ValueError("fa must be an instance of Fatty Acid or "
                             "an iterable of FattyAcids")

    def transform_type(
            self, fa: FattyAcid,
            specific_annotation: bool = False,
            reactions: dict = None
    ) -> Union[bool, str]:
        """
        Parameters
        ----------
        fa:
        specific_annotation:
        reactions

        Returns
        -------

        """
        # NOTE: this is assuming that input fas are transformable
        c_transform = abs(self.c_index - fa.c_index)
        db_transform = abs(self.db_index - fa.db_index)
        oh_transform = abs(self.oh_index - fa.oh_index)
        if reactions is None:
            # These are the standard reactions
            reactions = {
                "Chain length": {"C": 2, "DB": 0, "OH": 0},
                "Desaturation": {"C": 0, "DB": 1, "OH": 0},
                "Hydroxylation": {"C": 0, "DB": 0, "OH": 1}
            }
        for name, reaction in reactions.items():
            c = c_transform == reaction["C"]
            db = db_transform == reaction["DB"]
            oh = oh_transform == reaction["OH"]
            if c and db and oh:
                if specific_annotation:
                    react = ""
                    if c_transform == 0:
                        cs = [self.c_index, fa.c_index]
                        react += f"C {min(cs)} - {max(cs)}"
                    if db_transform == 0:
                        dbs = [self.db_index, fa.db_index]
                        if react == "":
                            react += f"DB {min(dbs)} - {max(dbs)}"
                        else:
                            react += f" | DB {min(dbs)} - {max(dbs)}"
                    if oh_transform == 0:
                        ohs = [self.oh_index, fa.oh_index]
                        if react == "":
                            react += f"OH {min(ohs)} - {max(ohs)}"
                        else:
                            react += f" | OH {min(ohs)} - {max(ohs)}"
                    return react
                return name
        # TODO: remove this case because should never happen?
        if c_transform == 0 and db_transform == 0 and oh_transform == 0:
            return "FAs equal"
        else:
            return False


def sum_fatty_acids(fa_iter: List[FattyAcid]) -> FattyAcid:
    if not isinstance(fa_iter, list):
        raise ValueError("'fa_iter' must be of type {0}, "
                         "not{1}".format(list,
                                         type(fa_iter).__name__))
    rsum = fa_iter[0]
    for i in range(1, len(fa_iter)):
        rsum += fa_iter[i]
    return rsum


if __name__ == "__main__":
    # implemented tests:
    # * is_equal (single + list)
    # * is_transformable (single + list)
    # * __hash__ => usage in sets
    # * __lt__ => usage in sort (especially consistency
    test_fa_strings = ["16:0:1", "18:1:0", "16:1:0", "16:0:0", "18:1:0"]
    test_fas = [FattyAcid(fa) for fa in test_fa_strings]

    print("Testing is_equal:")
    print(test_fas[0].is_equal(test_fas[1]) is False)
    print(test_fas[1].is_equal(test_fas[4]) is True)

    print("\nTesting is_transformable:")
    print(test_fas[0].is_transformable(test_fas[1]) is False)
    print(test_fas[1].is_transformable(test_fas[2]) is True)
    print(test_fas[1].is_transformable(test_fas[4]) is False)
    print(test_fas[3].is_transformable(test_fas[2]) is True)

    print("\ntesting __hash__:")
    print(test_fas[0].__hash__() == hash((16, 0, 1)))
    print([x.name for x in set(test_fas)])

    print("\ntesting __lt__ (sorting):")
    print([x.name for x in sorted(test_fas)])
