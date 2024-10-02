from typing import Iterable
from .exceptions import LipidSettingError, FaSettingError
from .Lipids.FattyAcid import FattyAcid


class FaDict:
    def __init__(self, dict_=None, **kwargs):
        if dict_ is None:
            self.dict = kwargs
        else:
            self.dict = dict_
        self.key_map = {}
        for key in self.dict.keys():
            if isinstance(key, tuple):
                for sub_key in key:
                    self.key_map[sub_key] = key

    def get(self, key: str):
        map_key = self.key_map.get(key)
        if map_key is None:
            return self.dict.get("general")
        else:
            return self.dict.get(map_key)


# aux
def _not_implemented_in_subclass(fname: str, cname: str):
    mess = "{0} is not implement in {1}, use a subclass of {1} instead"
    raise NotImplementedError(mess.format(fname, cname))


def _str_or_dict_error_(var: str, var_type: str):
    raise ValueError("'{0}' must a str or a dict, not a {1}".format(var, var_type))


def _update_dict_(to_update: dict, updater: dict) -> dict:
    # NOTE: this is only supposed to be used in cases where
    # ===== to_update is a dictionary containing lists as values.
    #       Otherwise more elegant options with python >= 3.5 are
    #       available.
    merged = {}
    for key, value in to_update.items():
        update_val = updater.get(key)
        if update_val is not None:
            if isinstance(update_val, Iterable):
                merged[key] = list(set(value + list(update_val)))
            else:
                merged[key] = value.append(update_val)
        else:
            merged[key] = value
    return merged


# TODO: add long-chain base
def read_class_data(class_file_path: str) -> dict:
    """

    Parameters
    ----------
    class_file_path

    Returns
    -------

    """
    key_map = {
        "ClassConnections": "class_connections",
        "ClassFANumbers": "class_nFA"
    }
    class_data = {
        "class_connections": {},
        "enzymes": {},
        "class_nFA": {}
    }
    curr_key = ""
    with open(class_file_path, "r") as file:
        for line_ in file:
            line = line_.strip().replace(" ", "")
            if line == "":
                continue
            if line.startswith(">"):
                key = line.replace(">", "")
                try:
                    curr_key = key_map[key]
                except KeyError:
                    raise LipidSettingError(
                        f"{key} is not an allowed keyword for lipid class input. "
                        f"Please use only {list(key_map.keys())}.",
                        error_type="class_key"
                    )
            elif curr_key == "class_connections":
                class_, connections = line.split(":")
                connected_classes = []
                enzymes = {}
                # checking if lipid class has no connections
                if connections.strip() != '':
                    for connected_class in connections.split(","):
                        if "|" in connected_class:
                            conn_class_, enzyme = connected_class.split("|")
                            enzymes[conn_class_] = enzyme
                            connected_classes.append(conn_class_)
                        else:
                            connected_classes.append(connected_class)
                class_data[curr_key][class_] = connected_classes
                class_data["enzymes"][class_] = enzymes
            elif curr_key == "class_nFA":
                class_, n = line.split(":")
                try:
                    class_data[curr_key][class_] = int(n)
                except ValueError:
                    raise LipidSettingError(
                        f"Class Fatty Acid numbers must be intergers. For {class_} "
                        f"{n} was provided, which could not be converted to an integer.",
                        error_type="fa_num"
                    )
    # TODO: make sure class connections are always symmetric
    return class_data


def read_fa_data(fa_file_path: str) -> dict:
    """

    Parameters
    ----------
    fa_file_path

    Returns
    -------

    """
    # TODO: add option for class-specific fatty acids
    key_map = {
        "FattyAcids": "fatty_acids",
        "ExcludedReactions": "excluded_reactions",
        "Reactions": "reactions"
    }
    fa_data = {
        "fatty_acids": FaDict(general=[]),
        "excluded_reactions": [],
        "reactions": {}
    }
    curr_key = ""
    fa_class_key = "general"
    with open(fa_file_path, "r") as file:
        for n, line_ in enumerate(file):
            line = line_.strip().replace(" ", "")
            if line == "" or line.startswith('#'):
                continue
            if line.startswith(">"):
                key = line.replace(">", "")
                try:
                    curr_key = key_map[key]
                except KeyError:
                    raise FaSettingError(
                        f"{key} is not an allowed keyword for fatty acid input.\n"
                        f"Please use only the following key words {list(key_map.keys())}.",
                        "invalid_keyword"
                    )
            elif curr_key == "fatty_acids":
                if line.startswith("Classes"):
                    lckey = line.split(":")[1]
                    if lckey == "":
                        fa_class_key = "general"
                    else:
                        fa_class_key = tuple(lckey.split(","))
                        fa_data[curr_key].dict[fa_class_key] = []
                else:
                    try:
                        fa_str = line.replace('l', '').replace('P-', '').replace('O-', '')
                        fa_data[curr_key].dict[fa_class_key].append(FattyAcid(fa_str))
                    except Exception as e:
                        if isinstance(e, KeyboardInterrupt):
                            raise e
                        else:
                            raise FaSettingError(
                                f"{line} in line {n} could not be read as a fatty acid internally.\n"
                                "Please make it is in the format 16:0;0 (where OH is optional). "
                                "For further details on file format check the tutorial page.",
                                "fa_format"
                            )
            elif curr_key == "excluded_reactions":
                reaction = [fa.replace('l', '').replace('P-', '').replace('O-', '')
                            for fa in sorted(line.split(","))]
                if len(reaction) != 2:
                    raise FaSettingError(
                        f"Excluded reaction {line} in line {n} is in the wrong format.\n"
                        "Please make sure it contains exactly two fatty acids split by a comma. "
                        "For further details on file format check the tutorial page.",
                        "excluded_reaction"
                    )
                fa_data[curr_key].append(reaction)
            elif curr_key == "reactions":
                reaction = {}
                name = ""
                line_elements = line.split(",")
                if len(line_elements) < 3:
                    raise FaSettingError(
                        f"Reaction '{line}' in line {n} is in the wrong format.\n"
                        "Please make sure it contains at least three fields with one integer value for C, DB and OH. "
                        "For further details on file format check the tutorial page.",
                        "reaction_rules"
                    )
                for i, attr in enumerate(line_elements):
                    if i > 2:
                        break
                    if i == 0:
                        try:
                            name, idx, val = attr.split(":")
                        except ValueError:
                            # optional naming
                            # TODO: does this cause downstream errors?
                            idx, val = attr.split(":")
                    else:
                        idx, val = attr.split(":")
                    try:
                        idx = idx.replace('l', '').replace('P-', '').replace('O-', '')
                        val = val.replace('l', '').replace('P-', '').replace('O-', '')
                        reaction[idx] = int(val)
                    except ValueError:
                        raise FaSettingError(
                            f"{val} in line {n} seems to be not convertable to an integer. "
                            "Please ensure all changes in reaction indices are integers.",
                            "reaction_rules"
                        )
                if len(reaction) != 3:
                    raise FaSettingError(
                        f"Reaction {line} in line {n} is in the wrong format. "
                        "Please make sure it contains one integer value for C, DB and OH. "
                        "For further details on file format check the tutorial page.",
                        "reaction_rules"
                    )
                fa_data[curr_key][name] = reaction
    return fa_data


if __name__ == "__main__":
    fa_data = read_fa_data("/home/nik/PhD/LipidNetworkInference/LipidNetworkInference/fatty_acids.txt")
    print("test finished")
