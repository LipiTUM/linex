function extractElements(elements, properties) {
    let extraction = [];
    let i, j, elem, tmp;
    for (i = 0; i < elements.length; i++) {
        elem = elements.get(i);
        if (elem !== null) {
            tmp = {};
            for (j = 0; j < properties.length; j++) {
               tmp[properties[j]] = elem[properties[j]];
            }
        }
        extraction.push(tmp);
    }
    return extraction;
}

function addNodePositions(extractedNodes, network) {
    let i, tmp;
    let positions = network.getPositions();
    for (i = 0; i < extractedNodes.length; i++) {
        tmp = extractedNodes[i];
        if (tmp !== null) {
            tmp["x"] = positions[i].x;
            tmp["y"] = positions[i].y;
        }
        extractedNodes[i] = tmp;
    }
    return extractedNodes;
}

function extractNodes(nodes, network, addPositions=false) {
    let props = ["id", "label", "color", "size", "hidden", "shape"];
    let extraction = extractElements(nodes, props);
    if (addPositions) return addNodePositions(extraction, network);
    return extraction
}

function extractEdges(edges, enrichment=false) {
    // TODO: only use non-hidden edges in network
    let props;
    if (enrichment) props = ['from', 'to', 'color']
    else props = ["id", "label", "color", "size", "from", "to", "hidden"];
    return extractElements(edges, props);
}


function writeNetwork(
    nodes, edges, network, version, legend=false,
    legendColourNodes=null, legendColourEdges=null,
    legendSizeNodes=null,
    colourLegend=null, sizeLegend=null,
    selectedType=''
) {
    let extractedNodes = extractNodes(nodes, network, true);
    let extractedEdges = extractEdges(edges, false);

    if (legend) {
        let legendNodes = {
            "colour": extractNodes(legendColourNodes, colourLegend),
            "size": extractNodes(legendSizeNodes, sizeLegend)
        };
        let legendEdges = {
            "colour": extractEdges(legendColourEdges, false),
        };
        let attributes = {
            "nodes": {
                "colour": document.getElementById("node_colours").value,
                "size": document.getElementById("node_sizes").value
            },
            "edges": {
                "colour": document.getElementById("edge_colours").value
            }
        }
        return JSON.stringify({
            "nodes": extractedNodes,
            "edges": extractedEdges,
            "legendNodes": legendNodes,
            "legendEdges": legendEdges,
            "attributes": attributes,
            "selectedType": selectedType,
            "version": version
        });
    }
    return JSON.stringify(
        {"nodes": extractedNodes, "edges": extractedEdges,
               "selectedType": selectedType, "version": version}
    );
}


function writeEnrichment(enrichment, component) {
    // finding currently selected comparison
    if (enrichment === null || enrichment === undefined) {
        window.alert('No enrichment selected!');
    }
    let enrich_nodes = {};
    let enrich_edges = {};
    let curr_nodes, node_positions;
    node_positions = enrichment_vis_networks[enrichment][component].getPositions();
    curr_nodes = [];
    for (const [node, attrs] of Object.entries(enrichment_vis_nodes[enrichment][component]._data)) {
        curr_nodes.push(
            {
                'id': node,
                'label': attrs['label'],
                'color': attrs['color'],
                'shape': attrs['shape'],
                'x': node_positions[node].x,
                'y': node_positions[node].y
            }
        );
    }
    enrich_nodes[component] = curr_nodes;
    enrich_edges[component] = extractEdges(enrichment_vis_edges[enrichment][component]);
    return JSON.stringify(
        {
            // TODO: add nodes, edges and scores
            "nodes": enrich_nodes,
            "edges": enrich_edges,
            "scores": enrichment_vis_scores[enrichment],
            "p_values": enrichment_pvalues[enrichment],
            "selectedType": "enrichment",
            "version": 2
        }
    )
}
