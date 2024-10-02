function showTab(contentName, linkName, event, tabId) {
    let i;
    let contents = document.getElementsByClassName(contentName);
    for (i = 0; i < contents.length; i++) {
        contents[i].style.display = "none"
    }
    let links = document.getElementsByClassName(linkName);
    for (i = 0; i < links.length; i++) {
        links[i].className = links[i].className.replace(" active", "");
    }
    document.getElementById(tabId).style.display = "block";
    event.currentTarget.className += " active";
}

function showMainTab(event, tabId) {
    showTab("main-tab", "main-tablink", event, tabId);
}

function showSummaryTab(event, tabId) {
    showTab('summary-tab', 'summary-link', event, tabId);
}

function showChainLengthTab(event, tabId) {
    showTab('chainlength-tab', 'chainlength-link', event, tabId);
}


function showEnrichmentTabs(event, tabId) {
    showTab('enrichment-tab', 'enrichment-link', event, tabId);
}
