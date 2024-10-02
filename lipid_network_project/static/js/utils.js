function fixNodes() {
    for (const id in nodes.getIds()) {
        let c_node = nodes.get(id);
        if (c_node != null) {
            c_node["physics"] = false;
            nodes.update(c_node);
        }
    }
}

function fixEdges() {
    for (const id in edges.getIds()) {
        let c_edge = edges.get(id);
        if (c_edge != null) {
            c_edge["physics"] = false;
            edges.update(c_edge);
        }
    }
}

function PhysicsToggle() {
    changeShownEdges()
    let checked = document.getElementById("physicsswitch");
    if (checked.checked){
        for (const id in nodes.getIds()) {
            let c_node = nodes.get(id);
            if (c_node != null) {
                c_node["physics"] = true;
                nodes.update(c_node);
            }
        }
        for (const id in edges.getIds()) {
            let c_edge = edges.get(id);
            if (c_edge != null) {
                if(c_edge["hidden"]){
                c_edge["physics"] = false;
                } else {
                c_edge["physics"] = true;
                }
                edges.update(c_edge);
            }
        }
    } else {
        fixNodes();
        fixEdges();
    }
}


// TODO: add properties and whether they need a group/comparison etc.
// TODO: add annotations
const propertyGroups = {
    "lipid_class": "none",
    "desaturation": "none",
    "chain_length": "none",
    "hydroxylation": "none",
    "c_index": "none",
    "db_index": "none",
    "oh_index": "none",
    "fold_changes": "comparison",
    "nlog_pvalues": "comparison",
    "degree": "none",
    "betweenness": "none",
    "closeness": "none",
    "correlation_changes": "comparison",
    "correlations": "group",
    "partial_correlation_changes": "comparison",
    "partial_correlations": "group",
    "reaction_types": "none",
    "reaction_type": "none"
};

const propertyContinuous = {
    "lipid_class": false,
    "desaturation": true,
    "chain_length": true,
    "hydroxylation": true,
    "c_index": true,
    "db_index": true,
    "oh_index": true,
    "fold_changes": true,
    "nlog_pvalues": true,
    "degree": true,
    "betweenness": true,
    "closeness": true,
    "correlation_changes": false,
    "correlations": true,
    "reaction_types": false,
    "reaction_type": false
};

function getSelectedProperty(selector) {
    const elem = document.getElementById(selector);
    if (elem !== null && elem !== undefined) {
        if (elem.options !== null) return elem.options[elem.selectedIndex].value;
    }
    return null;
}

function drawLegend(legend, legendNodes, legendEdges, options, containerId) {
    let legendData;
    let container = document.getElementById(containerId);

    // parsing and collecting nodes and edges from the python
    const step_x = 120;
    const step_y = 50;
    const x = -container.clientWidth / 2 + 50;
    const y = -container.clientHeight / 2 + 50;
    const ids = legendNodes.getIds();
    for (const id in ids) {
        const c_node = legendNodes.get(id)
        const new_x = x + c_node.x * step_x;
        const new_y = y + c_node.y * step_y;
        legendNodes.update(
            {
                id: c_node.id,
                x: new_x,
                y: new_y
            }
        );
    }

    legendData = {nodes: legendNodes, edges: legendEdges};
    legend = new vis.Network(container, legendData, options);

    return legend
}

function fitLegend(legend, legend_nodes) {
    legend.fit();
    let c_node;
    for (const id in legend_nodes.getIds()) {
        c_node = legend_nodes.get(id);
        if (c_node !== null) {
            // TODO: optimise this step
            c_node.x = c_node.x - 130;
            legend_nodes.update(c_node)
        }
    }
}

// annotations
function elementAnnotation(elements, attrAnnotation, hasDynamic) {
    let elem, val, opt;
    let annotation;
    for (const id in elements.getIds()) {
        elem = elements.get(id);
        if (elem != null) {
            annotation = elem["fixed_annotation"];
            if (hasDynamic) {
                for (const opt_i in Object.keys(attrAnnotation)) {
                    opt = Object.keys(attrAnnotation)[opt_i];
                    val = elem[attrAnnotation[opt]];
                    if (val !== undefined) {
                        if (typeof(val) === "number") {
                            // NOTE: Math.round rounds to integer
                            // => (x * 1000) / 1000 rounds to three decimal places
                            val = Math.round(val * 1000) / 1000;
                        }
                        annotation += "<br><b>" + opt + "</b>: " + val;
                    }
                }
            }
            elem["title"] = annotation;
            elements.update(elem);
        }
    }
    return elements;
}

function updateAnnotation() {
    // update node annotations based on a 'fixed'
    // annotation part and dynamic one dependent
    // on the groups
    // e.g fold-changes: fc_group1_group2
    let attrAnnotation = [{}];
    let pg, opt;
    let hasDynamic = false;
    if (document.getElementById("groups") !== null) {
        hasDynamic = true;
        const comp = getSelectedProperty("groups");
        // selecting single groups contained in comparison
        // only simple substring matching => not working for all group annotations!
        let compGroups = [];
        let singGroups = document.getElementById("single_groups");
        if (singGroups !== null) {
            let cGroup;
            for (let i = 0; i < singGroups.options.length; i++) {
                cGroup = singGroups.options[i].value;
                if (comp.includes(cGroup)) {
                    compGroups.push(cGroup);
                }
            }
        }
        let group = "";
        for (const opt_i in Object.keys(propertyGroups)) {
            opt = Object.keys(propertyGroups)[opt_i];
            pg = propertyGroups[opt];
            if (pg === "comparison") {
                attrAnnotation[opt] = opt + "_" + comp;
            } else if (pg === "group") {
                let optAnnot;
                for (let i = 0; i < compGroups.length; i++) {
                    group = compGroups[i];
                    optAnnot = opt + "_" + group;
                    attrAnnotation[optAnnot] = optAnnot;
                }
            }
        }
    } else {
        for (const opt_i in Object.keys(propertyGroups)) {
            opt = Object.keys(propertyGroups)[opt_i];
            pg = propertyGroups[opt];
            if (pg === "group") {
                hasDynamic = true;
                attrAnnotation[opt] = opt;
            }
        }
    }
    nodes = elementAnnotation(nodes, attrAnnotation, hasDynamic);
    edges = elementAnnotation(edges, attrAnnotation, hasDynamic);
}

function updateComparison() {
    if (propertyGroups[getSelectedProperty("node_colours")] === "comparison") {
        updateNodeColour();
    }
    if (propertyGroups[getSelectedProperty("node_sizes")] === "comparison") {
        updateNodeSize();
    }
    if (propertyGroups[getSelectedProperty("edge_colours")] === "comparison") {
        updateEdgeColour();
    }
    if (propertyGroups[getSelectedProperty("edge_sizes")] === "comparison") {
        updateEdgeSize();
    }
    updateAnnotation();
    fitLegend(colourLegend, legendNodeColours);
    // fitLegend(sizeLegend, legendNodeSizes);
    sizeLegend.fit();
}

function updateGroup() {
    if (propertyGroups[getSelectedProperty("node_colours")] === "group") {
        updateNodeColour();
    }
    if (propertyGroups[getSelectedProperty("node_sizes")] === "group") {
        updateNodeSize();
    }
    if (propertyGroups[getSelectedProperty("edge_colours")] === "group") {
        updateEdgeColour();
    }
    if (propertyGroups[getSelectedProperty("edge_sizes")] === "group") {
        updateEdgeSize();
    }
    updateAnnotation();
    fitLegend(colourLegend, legendNodeColours);
    // fitLegend(sizeLegend, legendNodeSizes);
    sizeLegend.fit()
}

function addGroupNames(property, addition = "") {
    if (propertyGroups[property] === "none") {
        if (addition === "")  return property;
        return property + "_" + addition;
    } else if (propertyGroups[property] === "group") {
        let select = getSelectedProperty("single_groups")
        // catching the case where no groups are given
        if (select === null) return property + "_" + addition;
        return property + "_" + addition + "_" + select;
    }
    return property + "_" + addition + "_" + getSelectedProperty("groups");
}

function addScale(property, selectedProp) {
    // checking if scale choice should be added
    // => needs to be continuous, group-dependent and groups have to be specified
    if (propertyGroups[selectedProp] === "none" ||
        !propertyContinuous[selectedProp] ||
        getSelectedProperty("single_groups") === null)
    {
        return property;
    }
    let commonScale = document.getElementById("scaleSwitch");
    if (commonScale.checked) return property + "_common";
    return property + "_individual";
}

function getDefault(dictObj, key, defaultValue) {
    let val = dictObj[key];
    return (typeof val !== "undefined") ? val : defaultValue;
}

function updateAttribute(
    attributeType, legendId, legend, elements, legendElements,
    defaultValue
) {
    // TODO: catch errors and return a pop-up message + stop message print
    document.getElementById('message_box_text').innerHTML =
        '<i class="fa fa-spinner slow-spin fa-fw" style="color: #0065bd"></i>&nbsp;Updating Network...'
    document.getElementById('foggy').style.display = 'inline';
    let legendProp = getSelectedProperty(legendId);
    let prop = addGroupNames(legendProp, attributeType);
    // adding whether the continuous, group dependent
    prop = addScale(prop, legendProp);
    let typeName;
    if (attributeType === "colour") {
        typeName = "color";
    } else {
        typeName = "size";
    }
    // network elements
    const ids = elements.getIds();
    let elem, borderColour, borderWidth;
    // this leaves the option to set a user-defined cutoff
    let thresh = -Math.log10(parseFloat(getCutoff("significanceThreshold")));
    let showSignificance = document.getElementById("significantSwitch").checked;
    let pvalColour = typeName === "color" && prop.startsWith("nlog_pvalues");
    for (const id in ids) {
        elem = elements.get(id);
        // FIXME: change prop with scaleProp
        if (elem != null) {
            if (pvalColour && showSignificance) {
                if (elem[prop.replace("_colour", "")] > thresh) {
                    borderColour = "#00cb2b";
                    borderWidth = 3;
                } else {
                    borderColour = "rgba(255,255,255,0.27)";
                    borderWidth = 0;
                }
                elem['color'] = {
                    'background': getDefault(elem, prop, defaultValue),
                    'border': borderColour
                }
                elem['borderWidth'] = borderWidth;
            } else {
                elem[typeName] = getDefault(elem, prop, defaultValue);
            }
            elements.update(elem);
        }
    }
    if (pvalColour) borderColourExplanation("significanceAnnotation", false);
    else if (legendId.startsWith("node_colours")) {
        borderColourExplanation("significanceAnnotation", true);
    }
    // legend elements
    // FIXME: add in scale_prop
    let propBak = prop + "_" + attributeType;
    const legendIds = legendElements.getIds();
    let legendElem;
    if (legendId === "node_sizes") {
        for (const id in legendIds) {
            legendElem = legendElements.get(id);
            if (legendElem != null) {
                const elem_prop =  legendElem[propBak];
                if (elem_prop != null) {
                    legendElem[typeName] = legendElem[propBak];
                    legendElem["label"] = legendElem[propBak + "_label"];
                    legendElem["hidden"] = false;
                    legendElements.update(legendElem);
                } else {
                    // These are the legends for enzyme/lipid species shapes
                    if (document.getElementById('network_type') !== null) {
                        legendElem["hidden"] =
                            document.getElementById('network_type').value.indexOf('bipartite') === -1;
                    } else {
                        legendElem["hidden"] = true;
                    }
                }
            }
            legendElements.update(legendElem);
        }
    } else {
        for (const id in legendIds) {
            legendElem = legendElements.get(id);
            if (legendElem != null) {
                const elem_prop =  legendElem[propBak];
                if (elem_prop != null) {
                    legendElem[typeName] = legendElem[propBak];
                    legendElem["label"] = legendElem[propBak + "_label"];
                    legendElem["hidden"] = false;
                    legendElements.update(legendElem);
                } else {
                    legendElem["hidden"] = true;
                }
            }
            legendElements.update(legendElem);
        }
    }
    if (legendProp !== "fold_changes") {
        if (attributeType === "colour") {
            document.getElementById("srcGroupFcColour").innerHTML = "";
            document.getElementById("tgtGroupFcColour").innerHTML = "";
        } else if (attributeType === "size") {
            document.getElementById("srcGroupFcSize").innerHTML = "";
            document.getElementById("tgtGroupFcSize").innerHTML = "";
        }
    } else {
        if (attributeType === "colour") {
            logFcExplanation(
                document.getElementById('srcGroupFcColour'),
                document.getElementById('tgtGroupFcColour'),
                attributeType
            );
        } else if (attributeType === "size") {
            logFcExplanation(
                document.getElementById('srcGroupFcSize'),
                document.getElementById('tgtGroupFcSize'),
                attributeType
            );
        }
	}
	document.getElementById('foggy').style.display = 'none';
}

function updateNodeColour() {
    updateAttribute(
        "colour", "node_colours", colourLegend,
        nodes, legendNodeColours, "#7f7f7f"
    );
    fitLegend(colourLegend, legendNodeColours);
}
function updateNodeSize() {
    updateAttribute(
        "size", "node_sizes", sizeLegend,
        nodes, legendNodeSizes, 10
    );
    // fitLegend(sizeLegend, legendNodeSizes);
    sizeLegend.fit()
}
function updateEdgeColour() {
    updateAttribute(
        "colour", "edge_colours", colourLegend,
        edges, legendEdgeColours, "#7f7f7f"
    );
    fitLegend(colourLegend, legendNodeColours);
}
function updateEdgeSize() {
    updateAttribute(
        "size", "edge_sizes", sizeLegend,
        edges, legendEdgeSizes, 1
    );
    // fitLegend(sizeLegend, legendNodeSizes);
    sizeLegend.fit();
}

function changeShownEdges(version) {
    let e = document.getElementById("showEdgeType").value;
    let checked = document.getElementById("physicsswitch").checked;
    if (e === "all"){
        for (const eid in edges.getIds()) {
            edges.update([{id: eid, hidden: false, physics: true}]);
        }
    } else if (e === "class") {
        for (const eid in edges.getIds()) {
            curr_edge = edges.get(eid);
            if (version === 1) {
                if (curr_edge.reaction_types === "Chain length" || curr_edge.reaction_types === "Desaturation") {
                    edges.update([{id: eid, hidden: true, physics: false}]);
                } else {
                    edges.update([{id: eid, hidden: false, physics: true}]);
                }
            } else {
                if (['L_FAmodify', 'L_FAether'].includes(curr_edge.reaction_type)) {
                    edges.update([{id: eid, hidden: true, physics: false}]);
                } else {
                    edges.update([{id: eid, hidden: false, physics: true}]);
                }
            }
        }
    } else if (e == "fa"){
        for (const eid in edges.getIds()) {
            curr_edge = edges.get(eid);
            if (version === 1) {
                if (curr_edge.reaction_types === "Chain length" || curr_edge.reaction_types === "Desaturation"){
                    edges.update([{id: eid, hidden: false, physics: true}]);
                } else {
                edges.update([{id: eid, hidden: true, physics: false}]);
                }
            } else {
                if (curr_edge.reaction_type.startsWith('L_FA')){
                    edges.update([{id: eid, hidden: false, physics: true}]);
                } else {
                edges.update([{id: eid, hidden: true, physics: false}]);
                }
            }
        }
    }
}

function findLipid() {
    let name = document.getElementById("lipid_search").value;
    let node, lipidId;
    for (const id in nodes.getIds()) {
        node = nodes.get(id);
        if (node !== null) {
            if (node["label"] === name)	{
                lipidId = node["id"];
                break;
            }
        }
    }
    network.focus(lipidId, {'scale': 2});
    let selection = [lipidId];
    network.selectNodes(selection);
}

function findSubLipid() {
    let name = document.getElementById("sub_search").value;
    let re = RegExp(name);
    let node;
    let selection = [];
    for (const id in nodes.getIds()) {
        node = nodes.get(id);
        if (node !== null) {
            // with regex
            if (re.test(node["label"])) {
                selection.push(node["id"]);
            }
        }
    }
    if (selection.length > 0) {
        network.selectNodes(selection);
    }
}

// edge selection grey-out
function fadeElements(
    legendSelection, legendElements,
    networkElements, propertySelector
) {
    let property = getSelectedProperty(propertySelector);
    property = addGroupNames(property, "colour");
    let legendProperty = property + "_colour";

    let selectedGroups = [];
    let elem;
    for (const i in legendSelection) {
        elem = legendElements.get(legendSelection[i]);
        if (elem !== null) {
            selectedGroups.push(elem[legendProperty]);
        }
    }

    for (const id in networkElements.getIds()) {
        elem = networkElements.get(id);
        if (elem !== null) {
            if (!selectedGroups.includes(elem[property])) {
                elem['color'] = "#e1e5e7";
            } else {
                elem['color'] = elem[property];
            }
            networkElements.update(elem);
        }
    }
}

function fadeUnselectedNodes() {
    let selection = colourLegend.getSelectedNodes();
    if (selection.length === 0) {
        updateNodeColour();
    } else {
        fadeElements(
            selection, legendNodeColours,
            nodes, "node_colours"
        );
    }
}
function fadeUnselectedEdges() {
    let selection = colourLegend.getSelectedEdges();
    if (selection.length === 0) {
        updateEdgeColour();
    } else {
        fadeElements(
            selection, legendEdgeColours,
            edges, "edge_colours"
        );
    }
}

function showLegendNavigation() {
    let navs = document.getElementsByClassName("vis-button");
    let i;
    if (document.getElementById("legend_nav").checked) {
        for (i = 0; i < navs.length; i++) {
            navs[i].style.display = 'none';
        }
    } else {
        for (i = 0; i < navs.length; i++) {
            navs[i].style.display = 'inline';
        }
    }
}

function getComparisonGroups() {
    let comp = getSelectedProperty("groups");
    let groupNames = document.getElementById("single_groups").options;
    let currGroup, srcGroup, tgtGroup;
    for (var i = 0; i < groupNames.length; i++) {
        currGroup = groupNames[i].value;
        if (comp.includes(currGroup + "_")) {
            srcGroup = currGroup;
        } else if (comp.includes("_" + currGroup)) {
            tgtGroup = currGroup;
        }
    }
    return {srcGroup, tgtGroup};
}

function logFcExplanation(srcElement, tgtElement, attributeType) {
    let compGroups = getComparisonGroups();
    let highStr, lowStr;
    if (attributeType === "colour") {
        highStr = "<b style=\"color: #00004c\">Blue Nodes</b>:<br>higher abundance in '" + compGroups["tgtGroup"] + "'";
        lowStr = "<b style=\"color: #800000\">Red Nodes</b>:<br>higher abundance in '" + compGroups["srcGroup"] + "'";
    } else {
        highStr = "<b>Large Nodes</b>:<br>higher abundance in '" + compGroups["tgtGroup"] + "'";
        lowStr = "<b>Small Nodes</b>:<br>higher abundance in '" + compGroups["srcGroup"] + "'";
    }
    srcElement.innerHTML = highStr;
    tgtElement.innerHTML = lowStr;
}

function getCutoff(cutoffElement) {
    return document.getElementById(cutoffElement).value;
}

function borderColourExplanation(elementId, remove=false) {
    if (remove) {
        document.getElementById(elementId).innerHTML = "";
    } else {
        let thresh = getCutoff("significanceThreshold");
        document.getElementById(elementId).innerHTML = "<b style='color: #00cb2b'>Node Border</b>: FDR < " + thresh;
    }
}

function updateCutoff() {
    updateAttribute(
        "colour", "node_colours", colourLegend,
        nodes, legendNodeColours
    )
}

function updateScale() {
    let nodeColourAttr = getSelectedProperty("node_colours");
    if (propertyGroups[nodeColourAttr] !== "none") {
        updateNodeColour();
    }
    let edgeColourAttr = getSelectedProperty("edge_colours");
    if (propertyGroups[edgeColourAttr] !== "none") {
        updateEdgeColour();
    }
    let nodeSizeAttr = getSelectedProperty("node_sizes");
    if (propertyGroups[nodeSizeAttr] !== "none") {
        updateNodeSize();
    }
    let edgeSizeAttr = getSelectedProperty("edge_sizes");
    if (propertyGroups[edgeSizeAttr] !== "none") {
        updateEdgeSize();
    }
}

/**
 * # ########################## #
 * # LINEX 2 specific additions #
 * # ########################## #
 */
function empty_network() {
    let i;
    for (i = 0; i < nodes.length(); i++) {
        nodes.delete(i);
    }
    for (i = 0; i < edges.length(); i++) {
        edges.delete(i);
    }
}


function updateNetworkType() {
    document.getElementById('message_box_text').innerHTML =
        '<i class="fa fa-spinner slow-spin fa-fw" style="color: #0065bd"></i>&nbsp;Loading network'
    document.getElementById('foggy').style.display = 'inline';
    let selected_type = document.getElementById('network_type').value;
    drawGraph(network_dict[selected_type]);
    updateNodeColour();
    updateNodeSize();
    updateEdgeColour();
    updateEdgeSize();
    updateAnnotation();
    current_network = selected_type;
    document.getElementById('foggy').style.display = 'none';
}
