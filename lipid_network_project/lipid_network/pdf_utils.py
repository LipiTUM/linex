import math
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import to_rgba
from .LINEX1.Network.vis_utils import ATTR_NAMES as LINEX1_ATTR_NAMES
from linex2.vis_utils import ATTR_NAMES as LINEX2_ATTR_NAMES

MAX_LEGEND_ROWS = 10


def _scale_nodes_(node_sizes, scale=(250, 2100)):
    sizes = np.array(list(node_sizes.values()))
    min_ = sizes.min()
    max_ = sizes.max()
    if min_ == max_:
        sizes = np.ones(sizes.shape) * 400
    else:
        max_ -= min_
        sizes = (sizes - min_) / max_
        sizes = (sizes * (scale[1] - scale[0])) + scale[0]
    return (min_, max_), dict(zip(node_sizes.keys(), sizes))


def _check_node_(node):
    """
    Checks whether a given is to be shown in the node legend
    based on its label and hidden status.
    All dummy nodes (generated for legend completeness and
    to draw edge legend) and hidden nodes (not containing
    information on the currently displayed property) are
    hereby excluded.
    """
    return node["label"] != "" and \
           "dummy" not in node["label"] and \
           not node.get("hidden", False)


def node_legend(
        node_data,
        node_scale=None, attr="color",
        data_scale=(17, 45)):
    if attr == "color":
        return [
            Line2D(
                [0], [0], color="w", marker="o",
                markersize=25,
                markerfacecolor=to_rgba(node[attr]),
                markeredgecolor=to_rgba("#D3D3D3"),
                label=node["label"]
            )
            for node in node_data if _check_node_(node)
        ]
    elif attr == "size":
        def scale(x, src, tgt):
            return (((x - src[0]) / src[1]) * (tgt[1] - tgt[0])) + tgt[0]
        return [
            Line2D(
                [0], [0], color="w", marker="o",
                markersize=scale(node[attr], node_scale, data_scale),
                markerfacecolor=to_rgba(node["color"]),
                label=node["label"]
            )
            for node in node_data if _check_node_(node)
        ]
    else:
        raise ValueError(
            f"attr must be 'size' or 'color', not {attr}"
        )


def edge_legend(edge_data):
    return [
        Line2D(
            [0], [0], color=edge["color"],
            label=edge["label"], lw=7
        )
        for edge in edge_data if not edge["hidden"]
    ]


def plot_from_json(
    nodes, edges, network_type, ax, fs=7
):
    graph = nx.Graph()
    # extracting nodes
    positions = {}
    node_colors = {}
    vis_node_sizes = {}
    node_labels = {}
    node_shapes = {}
    for node in nodes:
        graph.add_node(node["id"], **node)
        positions[node["id"]] = np.array([node["x"], -node["y"]])
        node_colors[node["id"]] = node["color"]
        vis_node_sizes[node["id"]] = node.get("size", 400)
        node_labels[node["id"]] = node["label"]
        node_shapes[node["id"]] = node["shape"]
    # scaling node sizes from visjs to networkx
    node_scale, node_sizes = _scale_nodes_(vis_node_sizes)
    # extracting edges
    edge_colors = {}
    for edge in edges:
        if not edge.get("hidden", False):
            graph.add_edge(
                edge["from"], edge["to"],
                **edge
            )
            edge_id = tuple(sorted((edge["from"], edge["to"])))
            edge_colors[edge_id] = edge["color"]
    # generating label positions from node positions
    if network_type == 'bipartite_enrichment':
        label_positions = positions
    else:
        label_positions = {node: pos - np.array([0, vis_node_sizes[node] * 1.05])
                           for node, pos in positions.items()}
    if 'bipartite' in network_type:
        triangle_nodes = []
        default_nodes = []
        for node in graph.nodes:
            if node_shapes[node] == 'triangle':
                triangle_nodes.append(node)
            else:
                default_nodes.append(node)
        # for bipartite networks we need to draw different node shapes
        nx.draw_networkx_nodes(
            graph, positions, ax=ax,
            nodelist=triangle_nodes, node_shape='^',
            node_color=[node_colors[node] for node in triangle_nodes],
            node_size=[node_sizes[node] for node in triangle_nodes],
            edgecolors=to_rgba("#D3D3D3")
        )
        nx.draw_networkx_nodes(
            graph, positions, ax=ax,
            nodelist=default_nodes,
            node_color=[node_colors[node] for node in default_nodes],
            node_size=[node_sizes[node] for node in default_nodes],
            edgecolors=to_rgba("#D3D3D3")
        )
    else:
        nx.draw_networkx_nodes(
            graph, positions, ax=ax,
            node_color=[node_colors[node] for node in graph.nodes],
            node_size=[node_sizes[node] for node in graph.nodes],
            edgecolors=to_rgba("#D3D3D3")
        )
    nx.draw_networkx_edges(
        graph, positions, ax=ax,
        edge_color=[edge_colors.get(edge, 'tab:grey') for edge in graph.edges]
    )
    nx.draw_networkx_labels(
        graph, label_positions, ax=ax,
        labels=node_labels, font_size=fs
    )
    if 'enrichment' not in network_type:
        return node_scale


def plot_network_pdf(json_data, file, out_format="pdf", fs=7):

    fig, ax = plt.subplots(figsize=(23, 21))
    node_scale = plot_from_json(json_data['nodes'], json_data['edges'],
                                json_data['selectedType'], ax, fs=fs)
    # legends
    # NOTE: node colour is the only attribute with possibly a high number of
    # discrete classes in the legend (when many lipid classes were measured)
    node_colour_handles = node_legend(json_data["legendNodes"]["colour"])
    max_cols = math.ceil(len(node_colour_handles) / MAX_LEGEND_ROWS)
    max_rows = MAX_LEGEND_ROWS if max_cols > 1 else len(node_colour_handles)
    bbox_x = -.0875 * max(2, max_cols)

    attribute_map = LINEX1_ATTR_NAMES if int(json_data['version']) == 1 else LINEX2_ATTR_NAMES
    plt.gca().add_artist(
        plt.legend(
            handles=node_legend(json_data["legendNodes"]["colour"]),
            title=attribute_map[json_data["attributes"]["nodes"]["colour"]],
            loc="upper left",
            bbox_to_anchor=(bbox_x, 1),
            labelspacing=1, title_fontsize=21, fontsize=15,
            framealpha=.5, edgecolor="w", ncol=max_cols
        )
    )
    plt.gca().add_artist(
        plt.legend(
            handles=node_legend(json_data["legendNodes"]["size"],
                                node_scale=node_scale, attr="size"),
            title=attribute_map[json_data["attributes"]["nodes"]["size"]],
            loc="center left",
            bbox_to_anchor=(bbox_x, .605 + ((MAX_LEGEND_ROWS - max_rows) * .02)),
            labelspacing=1.8, title_fontsize=21, fontsize=15,
            framealpha=.5, edgecolor="w"
        )
    )
    plt.gca().add_artist(
        plt.legend(
            handles=edge_legend(json_data["legendEdges"]["colour"]),
            title=attribute_map[json_data["attributes"]["edges"]["colour"]],
            loc="lower left",
            bbox_to_anchor=(bbox_x, .315 + ((MAX_LEGEND_ROWS - max_rows) * .02)),
            labelspacing=1, title_fontsize=21, fontsize=15,
            framealpha=.5, edgecolor="w"
        )
    )
    # final plotting
    plt.axis("off")
    plt.savefig(file, format=out_format,
                bbox_inches="tight")
    plt.close(fig)


def plot_enrichment_scores(scores, ax):
    x = np.array(scores['x'], dtype=np.float64)
    y = np.array(scores['y'], dtype=np.float64)
    ax.plot(x, y, color='tab:blue')
    ax.scatter(x, y, color='tab:blue')
    ax.set_xlabel('Iteration')
    ax.set_ylabel('Score')
    # TODO: add p-value and max score (upper left?)
    ax.set_title('Score Progression')
    return np.nanmax(y)


def plot_enrichment_pdfs(json_data, file_base, out_format="pdf", fs=7):
    pdf_files = []
    for component in json_data['nodes'].keys():
        curr_file = f"{file_base}_{component.replace(' ', '_')}.pdf"
        nodes = json_data['nodes'][component]
        try:
            edges = json_data['edges'][component]
            scores = json_data['scores'][component][0]
            pval = json_data['p_values'][component]
        except KeyError:
            # TODO: how should we handle this case (should never happen)
            raise KeyError(
                f"component {component} not found in either edges or scores"
            )
        fig, axes = plt.subplots(figsize=(25, 9), ncols=2)
        # plotting network
        plot_from_json(nodes, edges, 'bipartite_enrichment', axes[0], fs=fs)
        axes[0].get_xaxis().set_visible(False)
        axes[0].get_yaxis().set_visible(False)
        axes[0].grid(False)
        axes[0].set_title('Enriched Subnetwork')
        # plotting scores
        max_score = plot_enrichment_scores(scores, axes[1])
        # formatting and saving figure
        fig.suptitle(
            f"{component.title()}\n"
            f"max. score: {round(max_score, ndigits=2)}, "
            f"p-value: {round(pval, ndigits=4)}"
        )
        plt.tight_layout()
        plt.savefig(curr_file, format=out_format)
        plt.close(fig)
        pdf_files.append(curr_file)
    return pdf_files
