function hide_chosen_group(src_id, tgt_id) {
    let selected_src = document.getElementById(src_id).value;
    let tgt = document.getElementById(tgt_id);
    let i;
    for (i = 0; i < tgt.options.length; i++) {
        tgt.options[i].disabled = tgt.options[i].label === selected_src;
    }
}


function getEnrichmentParameters() {
    let parameters = {
        "group1":  document.getElementById('group1').value,
        "group2":  document.getElementById('group2').value,
        "ratio_diff":  document.getElementById('ratio_diff').value,
        "force_overwrite":  document.getElementById('force_overwrite').value,
        "fa_penalty":  document.getElementById('fa_penalty').value,
        "min_size":  document.getElementById('min_size').value,
        "max_size":  document.getElementById('max_size').value,
        "max_iter":  document.getElementById('max_iter').value,
        "temp":  document.getElementById('temp').value,
        "nreac_penalty": document.getElementById('nreac_penalty').value
    }
    return parameters;
}


function plotEnrichment(enrichment) {
    let components = enrichment_subnetworks[enrichment]['subnetworks'];
    enrichment_vis_nodes[enrichment] = {};
    enrichment_vis_edges[enrichment] = {};
    enrichment_vis_networks[enrichment] = {};
    enrichment_pvalues[enrichment] = {};
    // plotting enrichment scores
    enrichment_vis_scores[enrichment] = {};
    enrichment_vis_scores_plots[enrichment] = {};
    let componentDiv;
    let componentDiv1;
    for (const [key, value] of Object.entries(components)) {

        componentDiv = enrichment + '_' + key + '_outer';
        // generating div to hold each component (title + subdiv with network + score plot)
        generateEnrichmentDiv(
            enrichment + '_ETab', componentDiv, 'outer'
        );
        // Title div
        generateEnrichmentTitleDiv(componentDiv, enrichment + key, value['max_score'], value['p_value'], key);
        enrichment_pvalues[enrichment][key] = value['p_value'];


        componentDiv1 = enrichment + '_' + key + '_inner' ;
        // generating div to hold network + score plot
        generateEnrichmentDiv(
            componentDiv, componentDiv1, 'component'
        );

        // generating network
        generateEnrichmentNetwork(
            componentDiv1, enrichment,
            key, value['nodes'], value['edges']
        );

        generateEnrichmentScorePlot(
            componentDiv1,
            enrichment,  key,
            enrichment_subnetworks[enrichment]['scores'][key]
        )

        let downloadButtonDiv = enrichment + '_' + key + '_download';
        generateEnrichmentDiv(componentDiv, downloadButtonDiv, '')

        let new_btn = false;
        let elem;
        for (elem of document.getElementById(downloadButtonDiv).children) {
            if (elem.tagName === "BUTTON") {
                new_btn = true;
                break;
            }
        }
        if (!new_btn) {
            let btn = document.createElement("button");
            btn.type = 'submit';
            btn.innerHTML = '<i class="fa fa-download fa-fw" aria-hidden="true"></i>Download View';
            document.getElementById(downloadButtonDiv).appendChild(btn);

            let infoText = document.createElement("p");
            document.getElementById(downloadButtonDiv).appendChild(infoText);
            if (to_string) {
                btn.onclick = function () {
                    getNetworkJSON(enrichment, key, true)
                };
                infoText.innerHTML = '' +
                    '<i class="fa fa-info-circle fa-fw" aria-hidden="true">' +
                    '</i>Downloads the network as a .json file, which can be uploaded to <a href="https://exbio.wzw.tum.de/linex/download"' +
						' target="_blank" rel="noopener noreferrer">https://exbio.wzw.tum.de/download/</a> to generate a pdf view';
            } else {
                btn.onclick = function () {
                    send_pdf_request(enrichment, key, true)
                };
                infoText.innerHTML = '' +
                    '<i class="fa fa-info-circle fa-fw" aria-hidden="true">' +
                    '</i>Downloads the network a pdf view</p>';
            }
        }
    }

    // showing in page
    document.getElementById(enrichment + '_ETabLink').hidden = false;
}


function parametersToURL(parameters, url) {
    let purl = url + '?';
    let key;
    for (key in parameters) {
        purl = purl + key + '=' + parameters[key] + '&';
    }
    return purl.substring(0, purl.length - 1);
}


function getEnrichmentUpdate(url) {
    let request = new XMLHttpRequest();
    request.open('GET', url, true);

    let csrftoken = getCookie('csrftoken');
    request.setRequestHeader("X-CSRFToken", csrftoken);
    request.setRequestHeader('Content-Type', 'text');
    request.onreadystatechange = function () {
        if (request.readyState === 4 && request.status === 200) {
            let net_data = JSON.parse(request.responseText);
            if (net_data['error']) {
                window.alert(net_data['message']);
                enrichment_progress = 1;
            } else {
                enrichment_progress = net_data['status'];
                if (net_data['status'] === 0) {
                    enrichment_subnetworks[net_data['enrichment']['groups']] = net_data['enrichment'];
                    latest_enrichment = net_data['enrichment']['groups'];
                    // plotting enrichment network
					plotEnrichment(latest_enrichment);
					// showing enrichment
					document.getElementById(net_data['enrichment']['groups'] + "_ETabLink").click();
                }
            }
        } else if (request.status !== 200) {
            window.alert(
                'Network enrichment failed with error ' + request.status
            );
            enrichment_progress = 1;
        }
    };
    request.send();
}


function sleep(ms) {
    return new Promise(resolve => setTimeout(resolve, ms));
}


async function updateEnrichmentProgress(url) {
	// document.getElementById('foggy').style.height = document.body.scrollHeight + 'px';
    document.getElementById('foggy').style.display = 'inline';
    while (enrichment_progress === 2) {
        await sleep(5000);
        getEnrichmentUpdate(url);
    }
    document.getElementById('foggy').style.display = 'none';
}

async function computeEnrichment(url) {
    let parameters = getEnrichmentParameters();
    let request_url = parametersToURL(parameters, url);

    let request = new XMLHttpRequest();
    request.open('POST', request_url, true);

    let csrftoken = getCookie('csrftoken');
    request.setRequestHeader("X-CSRFToken", csrftoken);
    request.setRequestHeader('Content-Type', 'text');

    document.getElementById('message_box_text').innerHTML =
        '<i class="fa fa-spinner slow-spin fa-fw" style="color: #0065bd"></i>&nbsp;Computing Enrichment<br>This might take a little while...'

    request.onreadystatechange = function () {
        if (request.readyState === 4 && request.status === 200) {
            let net_data = JSON.parse(request.responseText);
            if (net_data['error']) {
                window.alert(net_data['message']);
            } else {
                enrichment_progress = 2;
                updateEnrichmentProgress(request_url);
            }
        } else if (request.readyState === 4 && request.status !== 200) {
            window.alert(
                'Network enrichment failed with error ' + request.status
            );
        }
    };
    request.send();
}


/*
# ############################# #
# Utilities for enrichment tabs #
# ############################# #
 */
function formatGroupString(groupString) {
    let groupSplit = groupString;
    groupSplit = groupSplit.replace('[', '').replace(']', '').split(',');
    return groupSplit[0] + ' vs. ' + groupSplit[1];
}


function showEnrichmentSubnetwork(event=null, groupString=null) {
    showEnrichmentTabs(event, groupString + "_ETab");
}

function generateEnrichmentTitleDiv(parent, baseName, max_score=0, p_value=0, key=''){
    let div1;
        title_div_name = baseName + "_Title"
        if (document.getElementById(title_div_name) === null) {
            // creating new div if not existing
            div1 = document.createElement("div");
            div1.id = title_div_name;
        } else {
            div1 = document.getElementById(title_div_name);
        }
        div1.style.textAlign = 'center'
        div1.innerHTML ="<h5><b>" + key + "</b></h5><p> Max. score of subnetwork: " + max_score + " , p-value of subnetwork: " + p_value + "</p>";
        document.getElementById(parent).appendChild(div1);
}

function generateEnrichmentDiv(parent, baseName, plotType, make_title=false) {

    // Making plot divs
    let divName;
    if (plotType !== '') {
        if (plotType === 'outer') {
            divName = baseName;
        } else if (plotType !== 'component'){
            divName = baseName + '_' + plotType;
        } else {
            divName = baseName;
        }
    } else {
        divName = baseName;
    }
    let div;
    if (document.getElementById(divName) === null) {
        // creating new div if not existing
        div = document.createElement("div");
        div.id = divName;
    } else {
        div = document.getElementById(divName);
    }
    // setting div parent
    document.getElementById(parent).appendChild(div);

    if (plotType === 'network') {
        div.className = "enrichmentNetwork";
        return divName;
    } else if (plotType === 'score') {
        div.className = "enrichmentScore";
        return divName;
    } else if (plotType === 'tab-content') {
        div.className = "tabcontent";
    }  else if (plotType === 'tab-link') {
        div.className = "tablinks"
    } else if (plotType === 'component') {
        div.className = "component";
    } else if (plotType === 'outer') {
        div.className = "EnrichmentComponent"
    }//else {
      //  throw new Error('Unknown plot type ' + plotType);
    //}
}


function generateEnrichmentNetwork(parent, baseName, component, nodes, edges) {
    let containerId = generateEnrichmentDiv(parent, baseName + component, 'network', make_title=true);
    let container = document.getElementById(containerId);

    enrichment_vis_nodes[baseName][component] = new vis.DataSet(nodes);
    enrichment_vis_edges[baseName][component] = new vis.DataSet(edges);
    // adding nodes and edges to the graph
    let netData = {
        'nodes': enrichment_vis_nodes[baseName][component],
        'edges': enrichment_vis_edges[baseName][component]
    };
    enrichment_vis_networks[baseName][component] = new vis.Network(container, netData, enrichment_vis_options);
}


function generateEnrichmentScorePlot(parent, baseName, component, scores) {
    let containerId = generateEnrichmentDiv(parent, baseName + component, 'score', make_title=true);
    enrichment_vis_scores[baseName][component] = [scores];
    // adding nodes and edges to the graph
    // TODO: define options
    let layout = {
      title: 'Score Progression for ' + component,
      xaxis: {
        title: 'Iteration',
        showgrid: false,
        zeroline: false
      },
      yaxis: {
        title: 'Score',
        showline: false
      }
    };
    let config = {'responsive': true};
    enrichment_vis_scores_plots[baseName][component] = Plotly.newPlot(
        containerId, enrichment_vis_scores[baseName][component], layout,
        config
    );
}


function initComputedEnrichments() {
    let key;
    for (key of Object.keys(enrichment_subnetworks)) {
        plotEnrichment(key);
    }
}

