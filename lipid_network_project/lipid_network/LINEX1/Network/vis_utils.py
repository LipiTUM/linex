from __future__ import annotations
import os
from .misc import _tuple_to_string_
from typing import Union, Dict, Any, Tuple, List
from jinja2 import Template
from pyvis.network import Network
import networkx as nx
import webbrowser
import numpy as np
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize, to_hex
from warnings import warn
import json
from itertools import combinations


# Auxiliaries for computing annotations
STATIC_EDGE_PROPERTIES = {
    # edge attributes
    "reaction_types",
    "reaction_enzymes"
}
STATIC_NODE_PROPERTIES = {
    # node attributes
    "lipid_class",
    "chain_length",
    "desaturation",
    "c_index",
    "db_index",
    "degree",
    "betweenness",
    "closeness"
}
# For some reason integers are not JSON serialisable
# by jinja => casting to strings solves the issue
TO_STRING_ATTRIBUTES = {
    "chain_length",
    "desaturation",
    "c_index"
}
ATTR_NAMES = {
    # dynamic attributes
    "fold_changes": "Fold Changes",
    "nlog_pvalues": "-log10(FDR)",
    "correlation_changes": "Correlation Changes",
    "partial_correlation_changes": "Partial Correlation Changes",
    "correlations": "Correlations",
    "partial_correlations": "Partial Correlations",
    # static attributes
    "lipid_class": "Lipid Class",
    "reaction_types": "Reaction Types",
    "desaturation": "Desaturation",
    "db_index": "DB Index",
    "chain_length": "Chain Length",
    "c_index": "C Index",
    "degree": "Degree",
    "betweenness": "Betweenness Centrality",
    "closeness": "Closeness Centrality"
}


def _setup_scm_(values, n_steps):
    scm = ScalarMappable(
        norm=Normalize(vmin=values["min"],
                       vmax=values["max"]),
        cmap=values["cmap"]
    )
    vals = _continuous_steps_(values, n_steps)
    return scm, vals


def _continuous_steps_(
        colour_map: Dict[str, float],
        n_steps: int
) -> Dict[int, float]:
    legend = {0: colour_map["min"]}
    if colour_map["min"] == colour_map["max"]:
        if colour_map.get("nan") is not None:
            legend[1] = np.nan
        return legend
    step = (colour_map["max"] - colour_map["min"]) / (n_steps - 1)
    for i in range(1, n_steps - 1):
        legend[i] = legend[0] + (i * step)
    legend[n_steps - 1] = colour_map["max"]
    if colour_map.get("nan") is not None:
        legend[n_steps] = np.nan
    return legend


def _size_steps_(
        size_map: Dict[str, Tuple[float, float]],
        n_steps: int
) -> Dict[int, Tuple[float, float]]:
    legend = dict()
    legend[0] = size_map["min"]
    raw_step = (size_map["max"][0] - size_map["min"][0]) / (n_steps - 1)
    norm_step = (size_map["max"][1] - size_map["min"][1]) / (n_steps - 1)
    for i in range(1, n_steps - 1):
        legend[i] = (legend[0][0] + (i * raw_step), legend[0][1] + (i * norm_step))
    legend[n_steps - 1] = size_map["max"]
    if size_map.get("nan") is not None:
        legend[n_steps] = (np.nan, size_map["min"][1]/2)
    return legend


def _node_annotation_(node: dict):
    annotation = f"<strong>{node['label']}</strong><br><br>"
    for prop in STATIC_NODE_PROPERTIES:
        node_prop = node.get(prop)
        if node_prop is not None:
            if isinstance(node_prop, float):
                annotation += f"<b>{prop}</b>: {round(node_prop, ndigits=3)}<br>"
            else:
                annotation += f"<b>{prop}</b>: {node_prop}<br>"
    return annotation


def _edge_annotation_(edge: dict):
    annotation = ""
    for prop in STATIC_EDGE_PROPERTIES:
        edge_prop = edge.get(prop)
        if edge_prop is not None:
            if isinstance(edge_prop, float):
                annotation += f"<b>{prop}</b>: {round(edge_prop, ndigits=3)}<br>"
            else:
                annotation += f"<b>{prop}</b>: {edge_prop}<br>"
    return annotation


def _size_filter_(
        sizes: dict,
        n_steps: int
) -> Dict[float, float]:
    keys = np.array(list(sizes.keys()))
    key_dict = _continuous_steps_(
        {"min": np.min(keys),
         "max": np.max(keys)},
        n_steps=n_steps
    )
    vals = np.array(list(sizes.values()))
    val_dict = _continuous_steps_(
        {"min": np.min(vals),
         "max": np.max(vals)},
        n_steps=n_steps
    )
    return dict(zip(key_dict.values(), val_dict.values()))


class VisParser(Network):
    def __init__(self, template: Union[str, Template] = None,
                 **kwargs):
        super(VisParser, self).__init__(**kwargs)
        if template is not None:
            if isinstance(template, str):
                template = Template(open(template, "r").read())
            elif not isinstance(template, Template):
                raise ValueError(
                    "'template' must be a jinja2.Template object or string"
                )
        self.template = template
        self.legend = {"nodes": [], "edges": []}
        self._legend_size_ = 0
        # helper to generate multi column legend
        self._colour_legend_nodes_ = 5
        self._node_cols_ = 1
        self._node_rows_ = 8
        self._to_last_column_ = True

    @classmethod
    def from_pyvis_network(cls, network: Network,
                           template: Template = None,
                           notebook: bool = False) -> VisParser:
        vp = VisParser(template=template,
                       height=network.height,
                       width=network.width,
                       directed=network.directed,
                       bgcolor=network.bgcolor,
                       font_color=network.font_color)
        if notebook:
            vp.prep_notebook()
        return vp

    def _node_update_step_(
            self, attr, value,
            name, counter, colours_to_hex,
            scm, values, cont
    ):
        if attr == "size":
            self.legend["nodes"].append(
                {
                    attr: value,
                    "label": f"{round(name, ndigits=3)}",
                    "id": counter,
                    "x": 0,
                    "y": counter * 2,
                    "shape": "dot",
                    "hidden": False,
                    "fixed": True,
                    "physics": False,
                    "margin": {"top": 0, "right": 15,
                               "bottom": 0, "left": 0},
                }
            )
        elif attr == "color":
            if colours_to_hex:
                if cont:
                    colour = to_hex(scm.to_rgba(value))
                    name_ = f"{round(values[name], ndigits=3)}"
                    # name_ = "{:.3f}".format(values[name])
                else:
                    colour = to_hex(value)
                    name_ = str(name)
            else:
                colour = value
                name_ = str(name)
            self.legend["nodes"].append(
                {
                    "color": colour,
                    "label": name_,
                    "id": counter,
                    "x": 0,
                    "y": counter * 2,
                    "shape": "dot",
                    "hidden": False,
                    "fixed": True,
                    "physics": False,
                    "margin": {"top": 0, "right": 15,
                               "bottom": 0, "left": 0},
                }
            )
        else:
            self.legend["nodes"].append(
                {
                    attr: value,
                    "size": 40,
                    "label": name,
                    "id": counter,
                    "x": 0,
                    "y": counter * 2,
                    "shape": "dot",
                    "hidden": False,
                    "fixed": True,
                    "physics": False,
                    "margin": {"top": 0, "right": 15,
                               "bottom": 0, "left": 0},
                }
            )

    def _update_node_legend_(self, values, attr, i, colours_to_hex,
                             n_steps: int = 5):
        counter = i
        cont = False
        if "cmap" in values.keys():
            scm, values = _setup_scm_(values, n_steps)
            cont = True
        else:
            scm = None
        for name, value in values.items():
            self._node_update_step_(
                attr, value, name, counter,
                colours_to_hex, scm, values, cont
            )
            counter += 1
        return counter

    @staticmethod
    def _dummy_node_(xi, x, y):
        return {
            "id": xi,
            "label": f"dummy{xi}",
            "x": x,
            "y": y,
            "hidden": True,
            "fixed": True,
            "physics": False
        }

    def _update_edge_legend_(self, values, attr, i, colours_to_hex,
                             n_steps: int = 5):
        counter = i
        cont = False
        if "cmap" in values.keys():
            scm, values = _setup_scm_(values, n_steps)
            cont = True
        for name, value in values.items():
            xi = len(self.legend["nodes"])
            self.legend["nodes"].append(
                {
                    "id": xi,
                    "label": f"dummy{xi}",
                    "x": 1,
                    "y": counter * 2,
                    "hidden": True,
                    "fixed": True,
                    "physics": False,
                }
            )
            self.legend["nodes"].append(
                {
                    "id": xi + 1,
                    "label": f"dummy{xi + 1}",
                    "x": 2,
                    "y": counter * 2,
                    "hidden": True,
                    "fixed": True,
                    "physics": False
                }
            )
            if attr == "width":
                self.legend["edges"].append(
                    {
                        "value": value,
                        # "label": "{:.3f}".format(name),
                        "label": f"{round(name, ndigits=3)}",
                        "from": xi,
                        "to": xi + 1,
                        "hidden": False,
                        "fixed": True,
                        "physics": False,
                        "font": {
                            "align": "bottom",
                            "vadjust": 5
                        }
                    }
                )
            elif attr == "color":
                if colours_to_hex:
                    if cont:
                        colour = to_hex(scm.to_rgba(value))
                        # name_ = "{:.3f}".format(values[name])
                        name_ = f"{round(values[name], ndigits=3)}"
                    else:
                        colour = to_hex(value)
                        name_ = str(name)
                else:
                    colour = value
                    name_ = str(name)
                self.legend["edges"].append(
                    {
                        "color": colour,
                        "label": name_,
                        "from": xi,
                        "to": xi + 1,
                        "value": 1,
                        "hidden": False,
                        "fixed": True,
                        "physics": False,
                        "font": {
                            "align": "bottom",
                            "vadjust": 5
                        }
                    }
                )
            else:
                self.legend["edges"].append(
                    {
                        attr: value,
                        "label": name,
                        "from": xi,
                        "to": xi + 1,
                        "value": 1,
                        "hidden": False,
                        "fixed": True,
                        "physics": False,
                        "font": {
                            "align": "bottom",
                            "vadjust": 5
                        }
                    }
                )
            counter += 1
        return counter

    def generate_legend(self,
                        node_colours: dict = None,
                        node_sizes: dict = None,
                        node_shapes: dict = None,
                        edge_colours: dict = None,
                        edge_sizes: dict = None,
                        colours_to_hex: bool = True):
        # TODO: handling overwriting to emtpy
        self.legend = {"nodes": [], "edges": []}
        # nodes
        i = 0
        if node_colours is not None:
            i = self._update_node_legend_(node_colours, "color", i,
                                          colours_to_hex)
        if node_sizes is not None:
            # warn("node sizes are not yet supported in legends!")
            i = self._update_node_legend_(_size_filter_(node_sizes, 5),
                                          "size", i, colours_to_hex)
        if node_shapes is not None:
            i = self._update_node_legend_(node_shapes, "shape", i,
                                          colours_to_hex)
        n_nodes = i
        # edges
        i = 0
        if edge_colours is not None:
            i = self._update_edge_legend_(edge_colours, "color", i,
                                          colours_to_hex)
        if edge_sizes is not None:
            warn("edge sizes are not yet supported in legends!")
            i = self._update_edge_legend_(_size_filter_(edge_sizes, 5),
                                          "width", i,
                                          colours_to_hex)
        n_edges = i
        self._legend_size_ = n_nodes if n_nodes > n_edges else n_edges
        # TODO: check if it's possible to increase legend font size

    def _prepare_write_(
        self, template: Template = None,
        force_hierarchical: bool = False
    ):
        if self.legend is None:
            raise ValueError(
                "Please compute a legend before you call write_vis!"
            )
        if template is None:
            if self.template is None:
                template_path = os.path.join(
                    os.path.dirname(__file__),
                    "templates/vis_main.html"
                )
                template = Template(open(template_path, "r").read())
            else:
                if isinstance(self.template, str):
                    template = Template(open(self.template, "r").read())
                elif isinstance(self.template, Template):
                    template = self.template
                else:
                    raise ValueError(
                        "'template' must be a path or a jinja2 Template object"
                    )
        elif isinstance(template, str):
            template = Template(open(template, "r").read())
        elif not isinstance(template, Template):
            raise ValueError(
                "'template' must be a path or a jinja2 Template object"
            )
            # write network => same as pyvis - should do the trick
            # TODO: notebook compatibility
        use_link_template = False
        for n in self.nodes:
            title = n.get("title", None)
            if title:
                if "href" in title:
                    """
                    this tells the template to override default hover
                    mechanic, as the tooltip would move with the mouse
                    cursor which made interacting with hover data useless.
                    """
                    use_link_template = True
                    break
        nodes, edges, heading, height, width, options = self.get_network_data()

        for i in range(len(nodes)):
            node = nodes[i]
            for attr in ["degree", "betweenness", "closeness"]:
                node[attr] = f"{round(node[attr], ndigits=3)}"
        # check if physics is enabled
        if isinstance(self.options, dict):
            if 'physics' in self.options and 'enabled' in self.options['physics']:
                physics_enabled = self.options['physics']['enabled']
            else:
                physics_enabled = True
        else:
            physics_enabled = self.options.physics.enabled
        # TODO: add more flexibility to this?
        if len(self.edges) > 1200 or force_hierarchical:
            self.options["physics"].solver = "hierarchicalRepulsion"

        return {
            "nodes": nodes,
            "edges": edges,
            "heading": heading,
            "height": height,
            "width": width,
            "options": options,
            "physics_enabled": physics_enabled,
            "template": template,
            "use_link_template": use_link_template
        }

    def write_vis(
        self, name: str, template: Template = None,
        physics_options: dict = None,
        force_hierarchical: bool = False
    ):
        """
        Saving VisParser.network as a html file

        Parameters
        ----------
        name : str
            File name (full path) to write to

        template : jinja2.Template, optional
            Path to a custom html template to use.
            Must be compatible with pyvis and the template
            in '<>/templates/vis_main.html'.

        physics_options : optional
            physics options to pass to pyvis.Network

        force_hierarchical : bool, optional, default False
            Whether to use hierarchical repulsion as a solver
            for network layout. Always True when the number
            of edges exceeds 1200 for performance reasons.

        Returns
        -------
        None
        """
        network_params = self._prepare_write_(template, force_hierarchical)
        self.html = network_params["template"].render(
            nodes=network_params["nodes"],
            edges=network_params["edges"],
            options=network_params["options"],
            physics_enabled=network_params["physics_enabled"],
            use_DOT=self.use_DOT,
            dot_lang=self.dot_lang,
            widget=self.widget,
            bgcolor=self.bgcolor,
            conf=self.conf,
            tooltip_link=network_params["use_link_template"],
            legend_nodes=self.legend["nodes"],
            legend_edges=self.legend["edges"]
        )

        with open(name, "w+") as out:
            out.write(self.html)

    def show(self, name, **kwargs):
        """
        Showing the network in VisParser.network after saving as html

        Parameters
        ----------
        name : str
            File name (full path) to write to

        kwargs :
            Keyword arguments to pass to write_vis.
            Available options:
            template : jinja2.Template, optional
                Path to a custom html template to use.
                Must be compatible with pyvis and the template
                in '<>/templates/vis_main.html'.

            physics_options : optional
                physics options to pass to pyvis.Network

            force_hierarchical : bool, optional, default False
                Whether to use hierarchical repulsion as a solver
                for network layout. Always True when the number
                of edges exceeds 1200 for performance reasons.

        Returns
        -------
        None
        """
        self.write_vis(name, **kwargs)
        webbrowser.open(name)


class DynamicVisParser(VisParser):
    def __init__(
        self, template: Union[str, Template] = None,
        default_colour: str = "#0065bd",
        default_size: int = 3,
        random_seed: int = None, **kwargs
    ):
        if template is None:
            template = os.path.join(
                os.path.dirname(__file__),
                "templates", "vis_dynamic.html"
            )
        # TODO: change template to vis_dynamic.html
        super(DynamicVisParser, self).__init__(
            template=template, **kwargs
        )
        if random_seed is not None:
            layout = self.options.__dict__.get("layout", {})
            layout["randomSeed"] = random_seed
            self.options.__dict__["layout"] = layout
        self._default_colour_ = default_colour
        self._default_size_ = default_size

        self.legend = {
            "nodes": {"colour": [], "size": []},
            "edges": {"colour": [], "size": []},
        }
        # these are just helpers for ensuring visjs node id
        # consecutiveness requirement
        self._edges_by_node_ = {"colour": {}, "size": {}}
        self._legend_attributes_ = dict()
        self._visualisation_options_ = {
            "node_colours": [],
            "node_sizes": [],
            "edge_colours": []
        }

    @classmethod
    def from_pyvis_network(
        cls, network: Network,
        template: Template = None,
        notebook: bool = False
    ) -> DynamicVisParser:
        vp = DynamicVisParser(
            template=template,
            height=network.height,
            width=network.width,
            directed=network.directed,
            bgcolor=network.bgcolor,
            font_color=network.font_color
        )
        if notebook:
            vp.prep_notebook()
        return vp

    @classmethod
    def from_json(
        cls, network_file: str,
        legend_file: str, **kwargs
    ):
        vp = DynamicVisParser(**kwargs)
        # network
        network_data = json.load(open(network_file, "r"))
        vp.nodes = network_data["nodes"]
        vp.edges = network_data["edges"]
        # legend
        vp.legend = json.load(open(legend_file, "r"))
        return vp

    def _labels_to_ids_(self):
        """
        Auxiliary Function to modify edge and node label ids
        to be consecutive integers (otherwise this might lead to
        issues with vis). Also some node properties are explicitly
        set.
        """
        label_to_id = {}
        for i in range(len(self.nodes)):
            node = self.nodes[i]
            label_to_id[node["id"]] = i

            node["color"] = self._default_colour_
            node["physics"] = True
            node["id"] = i
            for attr in node.keys():
                if attr in TO_STRING_ATTRIBUTES:
                    node[attr] = str(node[attr])
            self.nodes[i] = node

        for j in range(len(self.edges)):
            edge = self.edges[j]
            edge["color"] = self._default_colour_
            edge["width"] = self._default_size_
            edge["from"] = label_to_id[edge["from"]]
            edge["to"] = label_to_id[edge["to"]]
            edge["id"] = j
            for attr in edge.keys():
                if attr in TO_STRING_ATTRIBUTES:
                    edge[attr] = str(edge[attr])
            self.edges[j] = edge

    def _update_node_colour_legend_(
        self, name: str, colours: Dict[Any, str],
        n_steps: int, colours_to_hex: bool
    ):
        n = len(self.legend["nodes"]["colour"])
        if name == "lipid_class_colour":
            nlc = len(colours)
            if nlc > 8:
                if nlc > 16:
                    self._node_cols_ = 3
                    if nlc > 24:
                        self._node_rows_ = nlc//3 + 1
                else:
                    self._node_cols_ = 2
                self._colour_legend_nodes_ = nlc
            elif nlc > 5:
                self._colour_legend_nodes_ = nlc
            self._to_last_column_ = 0 < self._colour_legend_nodes_ % self._node_rows_ >= 5

        cont = False
        if "cmap" in colours.keys():
            scm, colours = _setup_scm_(colours, n_steps)
            cont = True
        else:
            scm = None
        for i, (value, colour) in enumerate(colours.items()):
            if colours_to_hex:
                if cont:
                    colour = to_hex(scm.to_rgba(colour))
                else:
                    colour = to_hex(colour)
            if cont:
                # name_ = "{:.3f}".format(colours[i])
                name_ = f"{round(colours[i], ndigits=3)}"
            else:
                name_ = str(value)
            # NOTE: lipid_class comes first if present
            # therefore we do not need to check here!
            x = i // self._node_rows_
            y = i % self._node_rows_
            if i >= n:
                self.legend["nodes"]["colour"].append(
                    {
                        "id": i,
                        "label": "",
                        f"{name}_colour": to_hex(colour),
                        f"{name}_colour_label": name_,
                        "x": x,
                        "y": y * 2,
                        "shape": "dot",
                        "fixed": True,
                        "physics": False,
                        "color": self._default_colour_
                    }
                )
            else:
                elem = {
                    f"{name}_colour": to_hex(colour),
                    f"{name}_colour_label": name_
                }
                ui = (self._node_cols_ - 1) * self._node_rows_ + i
                if self._to_last_column_ or ui < self._colour_legend_nodes_:
                    self.legend["nodes"]["colour"][ui].update(elem)
                else:
                    self.legend["nodes"]["colour"].append(
                        {
                            "id": i,
                            "label": "",
                            f"{name}_colour": to_hex(colour),
                            f"{name}_colour_label": name_,
                            "x": self._node_cols_ - 1,
                            "y": i * 2,
                            "shape": "dot",
                            "fixed": True,
                            "physics": False,
                            "color": self._default_colour_
                        }
                    )

    def _update_node_size_legend_(
            self, name: str, sizes: Dict[str, Tuple[float, float]],
            n_steps: int = 5
    ):
        n = len(self.legend["nodes"]["size"])
        size_map = _size_steps_(sizes, n_steps)
        for i, (value, size) in size_map.items():
            if name == "fold_changes":
                size = abs(size)
            if i >= n:
                self.legend["nodes"]["size"].append(
                    {
                        "id": i,
                        "label": "",
                        f"{name}_size": size,
                        f"{name}_size_label": f"{round(value, ndigits=3)}",
                        "x": 1,
                        "y": i * 2,
                        "shape": "dot",
                        "fixed": True,
                        "physics": False,
                        "color": self._default_colour_
                    }
                )
            else:
                elem = {
                    f"{name}_size": size,
                    f"{name}_size_label": f"{round(value, ndigits=3)}"
                }
                self.legend["nodes"]["size"][i].update(elem)

    def _update_edge_colour_legend_(
        self, name: str, colours: Dict[Any, str],
        n_steps: int, colours_to_hex: bool
    ):
        n = len(self.legend["edges"]["colour"])
        cont = False
        if "cmap" in colours.keys():
            scm, colours = _setup_scm_(colours, n_steps)
            cont = True
        for i, (value, colour) in enumerate(colours.items()):
            xi = (2 * i) + 100
            if colours_to_hex:
                if cont:
                    colour = to_hex(scm.to_rgba(colour))
                else:
                    colour = to_hex(colour)
            if cont:
                name_ = f"{round(colours[i], ndigits=3)}"
            else:
                name_ = str(value)
            if i >= n:
                self.legend["nodes"]["colour"].append(
                    self._dummy_node_(xi, self._node_cols_, i * 2)
                )

                self.legend["nodes"]["colour"].append(
                    self._dummy_node_(xi + 1, self._node_cols_ + 1, i * 2)
                )
                self.legend["edges"]["colour"].append(
                    {
                        f"{name}_colour": colour,
                        f"{name}_colour_label": name_,
                        "from": xi,
                        "to": xi + 1,
                        "id": i,
                        "value": 1,
                        "hidden": False,
                        "fixed": True,
                        "physics": False,
                        "font": {
                            "align": "bottom",
                            "vadjust": 5
                        },
                        "color": self._default_colour_
                    }
                )
                self._edges_by_node_["colour"].setdefault(xi, []).append(i)
                self._edges_by_node_["colour"].setdefault(xi + 1, []).append(i)
            else:
                elem = {
                    f"{name}_colour": colour,
                    f"{name}_colour_label": name_
                }
                self.legend["edges"]["colour"][i].update(elem)

    def _update_edge_size_legend_(
        self, name: str, sizes: Dict[str, Tuple[float, float]],
        n_steps: int = 5
    ):
        n = len(self.legend["edges"]["size"])
        size_map = _size_steps_(sizes, n_steps)
        for i, (value, size) in enumerate(size_map.items()):
            xi = i + 100
            if i >= n:
                self.legend["nodes"]["size"].append(
                    self._dummy_node_(xi, 2, i * 2)
                )
                self.legend["nodes"]["size"].append(
                    self._dummy_node_(xi + 1, 3, i * 2)
                )
                self.legend["edges"]["size"].append(
                    {
                        f"{name}_size": size,
                        f"{name}_size_label": f"{round(value, ndigits=3)}",
                        "from": xi,
                        "to": xi + 1,
                        "value": 1,
                        "id": i,
                        "hidden": False,
                        "fixed": True,
                        "physics": False,
                        "font": {
                            "align": "bottom",
                            "vadjust": 5
                        },
                        "color": self._default_colour_
                    }
                )
                self._edges_by_node_["size"].setdefault(xi, []).append(i)
                self._edges_by_node_["size"].setdefault(xi + 1, []).append(i)
            else:
                elem = {
                    f"{name}_colour": size,
                    f"{name}_size_label": f"{round(value, ndigits=3)}"
                }
                self.legend["edges"]["size"][i].update(elem)

    def generate_legend(
            self,
            node_colours: Dict[str, Dict[Any, str]] = None,
            node_sizes: Dict[str, Dict[str, Tuple[float, float]]] = None,
            node_shapes: Dict[str, Dict[Any, str]] = None,
            edge_colours: Dict[str, Dict[Any, str]] = None,
            edge_sizes: Dict[str, Dict[str, Tuple[float, float]]] = None,
            colours_to_hex: bool = True,
            n_steps: int = 5
    ):
        # node legend
        if node_colours is not None:
            # lipid class is the only discrete node colour option
            # and can therefore be the only attribute with more than 5 legend nodes
            # to have proper multi-row legends it is easiest to do this one first
            class_colours = node_colours.pop("lipid_class_colour")
            if class_colours is not None:
                self._update_node_colour_legend_("lipid_class_colour", class_colours,
                                                 n_steps, colours_to_hex)
                self._legend_attributes_.update({"node_colour": "lipid_class_colour"})
                self._visualisation_options_["node_colours"].append("lipid_class_colour")
            for attr, colours_vals in node_colours.items():
                self._update_node_colour_legend_(attr, colours_vals,
                                                 n_steps, colours_to_hex)
                self._legend_attributes_.update({"node_colour": attr})
                self._visualisation_options_["node_colours"].append(attr)
        if node_sizes is not None:
            for attr, sizes_vals in node_sizes.items():
                self._update_node_size_legend_(attr, sizes_vals,
                                               n_steps)
                self._legend_attributes_.update({"node_size": attr})
                self._visualisation_options_["node_sizes"].append(attr)
        if node_shapes is not None:
            # TODO: implement node shape legend
            print("\nNode Shape legend currently not supported\n")
            # self._legend_attributes_.add("shape")
        # edge legend
        if edge_colours is not None:
            for attr, colours_vals in edge_colours.items():
                self._update_edge_colour_legend_(attr, colours_vals,
                                                 n_steps, colours_to_hex)
                self._legend_attributes_.update({"edge_colour": attr})
                self._visualisation_options_["edge_colours"].append(attr)
        if edge_sizes is not None:
            for attr, sizes_vals in edge_sizes.items():
                self._update_edge_size_legend_(attr, sizes_vals,
                                               n_steps)
                self._legend_attributes_.update({"edge_size": attr})
                self._visualisation_options_["edge_sizes"].append(attr)

    def _consecutive_nodes_(self, attr):
        for i in range(1, len(self.legend["nodes"][attr])):
            c_id = self.legend["nodes"][attr][i]["id"]
            if c_id - self.legend["nodes"][attr][i - 1]["id"] != 1:
                self.legend["nodes"][attr][i]["id"] = i
                for edge_id in self._edges_by_node_[attr].get(c_id, []):
                    c_edge: dict = self.legend["edges"][attr][edge_id]
                    if c_edge["from"] == c_id:
                        self.legend["edges"][attr][edge_id]["from"] = i
                    elif c_edge["to"] == c_id:
                        self.legend["edges"][attr][edge_id]["to"] = i

    def _add_annotations_(self):
        for i in range(1, len(self.nodes)):
            node = self.nodes[i]
            node["fixed_annotation"] = _node_annotation_(node)
            self.nodes[i] = node
        for i in range(1, len(self.edges)):
            edge = self.edges[i]
            edge["fixed_annotation"] = _edge_annotation_(edge)
            self.edges[i] = edge

    def write_vis(
        self, name: str,
        visualisation_options: Dict[str, List[str]],
        template: Union[Template, str] = None,
        physics_options: dict = None,
        force_hierarchical: bool = False,
        groups: Union[List[str], None] = None,
        as_dict: bool = False,
        group_combinations = None
    ):
        """
        Saving DynamicVisParser.network as a html file

        Parameters
        ----------
        name : str
            File name (full path) to write to

        visualisation_options: dict
            visualisation options to show in html
            TODO: describe properly

        template : jinja2.Template, optional
            Path to a custom html template to use.
            Must be compatible with pyvis and the template
            in '<>/templates/vis_main.html'.

        physics_options : optional
            physics options to pass to pyvis.Network

        force_hierarchical : bool, optional, default False
            Whether to use hierarchical repulsion as a solver
            for network layout. Always True when the number
            of edges exceeds 1200 for performance reasons.

        groups : list of str, optional
            groups for comparisons

        as_dict : bool, optional, default False
            do not use!!!

        Returns
        -------
        None
        """
        # NOTE: for some reason vis cannot handle non-consecutive ids
        # => Dummy nodes for edge legend are renamed at the time of saving.
        #    This enables the flexibility to later add other attributes and
        #    save again.
        self._labels_to_ids_()
        colour_legend = any(["colour" in attr for attr in self._legend_attributes_.keys()])
        size_legend = any(["size" in attr for attr in self._legend_attributes_.keys()])
        if colour_legend:
            self._consecutive_nodes_("colour")
        if size_legend:
            self._consecutive_nodes_("size")
        # adding fixed annotations to nodes and legends
        self._add_annotations_()
        # getting all relevant elements for saving
        network_params = self._prepare_write_(template, force_hierarchical)
        # ensuring compatibility when no group annotation is provided
        groups = [] if groups is None else groups
        if groups is None:
            group_combinations = []
        elif group_combinations is None:
            group_combinations = [_tuple_to_string_(comb) for comb in combinations(groups, 2)]
        else:
            group_combinations = [_tuple_to_string_(comb) for comb in group_combinations]

        if as_dict:
            return dict(
                nodes=network_params["nodes"],
                edges=network_params["edges"],
                options=network_params["options"],
                physics_enabled=network_params["physics_enabled"],
                dot_lang=self.dot_lang,
                widget=self.widget,
                bgcolor=self.bgcolor,
                conf=self.conf,
                tooltip_link=network_params["use_link_template"],
                legend_nodes=self.legend["nodes"]["colour"],
                legend_edges=self.legend["edges"]["colour"],
                legend_size_nodes=self.legend["nodes"]["size"],
                legend_size_edges=self.legend["edges"]["size"],
                legend_attributes=self._legend_attributes_,
                colour_legend=colour_legend,
                size_legend=size_legend,
                single_groups=groups,
                groups=group_combinations
            )

        self.html = network_params["template"].render(
            nodes=network_params["nodes"],
            edges=network_params["edges"],
            options=network_params["options"],
            physics_enabled=network_params["physics_enabled"],
            use_DOT=self.use_DOT,
            dot_lang=self.dot_lang,
            widget=self.widget,
            bgcolor=self.bgcolor,
            conf=self.conf,
            tooltip_link=network_params["use_link_template"],
            legend_colour_nodes=self.legend["nodes"]["colour"],
            legend_colour_edges=self.legend["edges"]["colour"],
            legend_size_nodes=self.legend["nodes"]["size"],
            legend_size_edges=self.legend["edges"]["size"],
            colour_legend=colour_legend,
            size_legend=size_legend,
            single_groups=list(groups),
            groups=group_combinations,
            node_colours={attr: ATTR_NAMES[attr] for attr in visualisation_options["node_colours"]},
            edge_colours={attr: ATTR_NAMES[attr] for attr in visualisation_options["edge_colours"]},
            node_sizes={attr: ATTR_NAMES[attr] for attr in visualisation_options["node_sizes"]}
        )

        with open(name, "w+") as out:
            out.write(self.html)

    def network_to_json(self, file: str, **kwargs):
        json.dump(
            {"nodes": self.nodes, "edges": self.edges},
            open(file, "w"),
            **kwargs
        )

    def to_networkx(self) -> Union[nx.Graph, nx.DiGraph]:
        excluded = {"physics", "size", "id", "width", "shape", "color"}
        graph = nx.DiGraph() if self.directed else nx.Graph()
        for node in self.nodes:
            graph.add_node(
                node["label"],
                **{name: node_attr for name, node_attr in node.items()
                    if name not in excluded}
            )
        for edge in self.edges:
            graph.add_edge(
                edge["from"], edge["to"],
                ** {name: edge_attr for name, edge_attr in edge.items()
                    if name not in excluded}
            )
        return graph

    def legend_to_json(self, file: str, **kwargs):
        json.dump(self.legend, open(file, "w"), **kwargs)

    def legend_to_networkx(self, legend: str = "colour") -> nx.Graph:
        if legend != "colour" and legend != "size":
            raise ValueError(
                f"'legend' must be 'colour' or 'size', not {legend}"
            )
        graph = nx.Graph()
        for node in self.legend["nodes"][legend]:
            graph.add_node(
                node["id"],
                **{name: node_attr for name, node_attr in node.items()
                   if name.endswith(legend) or name.endswith("label")}
            )
        for edge in self.legend["edges"][legend]:
            graph.add_edge(
                edge["from"], edge["to"],
                **{name: edge_attr for name, edge_attr in edge.items()
                   if name.endswith(legend) or name.endswith("label")}
            )
        return graph


if __name__ == "__main__":
    import networkx as nx

    g = nx.karate_club_graph()
    colors = {node: "red" if club == "Officer" else "blue"
              for node, club in nx.get_node_attributes(g, "club").items()}
    nx.set_node_attributes(g, colors, "color")

    vip = VisParser()
    vip.from_nx(g)
    vip.generate_legend(node_colours={"Officer": "red", "Mr. Hi": "blue", "Noop": "green"},
                        edge_colours={"freak": "lightgrey"})
    vip.show("test.html")

    vipc = VisParser(template="templates/vis_graph.html")
    vipc.from_nx(g)
    vipc.generate_legend(node_colours={"Officer": "red", "Mr. Hi": "blue", "Noop": "green"},
                         edge_colours={"freak": "lightgrey"})
    vipc.show("test_c.html")
    print()
