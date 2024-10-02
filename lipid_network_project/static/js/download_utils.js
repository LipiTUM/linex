function getCookie(name) {
    let cookieValue = null;
    if (document.cookie && document.cookie !== '') {
        let cookies = document.cookie.split(';');
        for (let i = 0; i < cookies.length; i++) {
            let cookie = jQuery.trim(cookies[i]);
            if (cookie.substring(0, name.length + 1) === (name + '=')) {
                cookieValue = decodeURIComponent(cookie.substring(name.length + 1));
                break;
            }
        }
    }
    return cookieValue;
}


function complete_download(blob, contentTypeHeader, filename, newBlob=true) {
    let downloadLink = window.document.createElement('a');
    if (newBlob){
        downloadLink.href = window.URL.createObjectURL(
            new Blob([blob], {type: contentTypeHeader}));
    } else {
        downloadLink.href = window.URL.createObjectURL(blob);
    }
    downloadLink.download = filename;
    document.body.appendChild(downloadLink);
    downloadLink.click();
    document.body.removeChild(downloadLink);
    window.URL.revokeObjectURL((downloadLink.href));
}


function getFileName(xhr) {
    let filename = "";
    let disposition = xhr.getResponseHeader('Content-Disposition');
    // check if filename is given
    if (disposition && disposition.indexOf('attachment') !== -1) {
        let filenameRegex = /filename[^;=\n]*=((['"]).*?\2|[^;\n]*)/;
        let matches = filenameRegex.exec(disposition);
        if (matches != null && matches[1]) filename = matches[1].replace(/['"]/g, '');
    }
    return filename;
}


function processResponse(main, request) {
    if (main.status === 200) {
        let filename = getFileName(request);
        let blob = main.response;
        if (window.navigator.msSaveOrOpenBlob) {
            window.navigator.msSaveBlob(blob, filename);
        }
        else {
            complete_download(blob, request.getResponseHeader("Content-Type"), filename)
        }
    } else {
        alert('Download failed, please try again or contact the developers if the issue remains')
    }
}

function processFetchHeader(header) {
    let content = header.split(';');
    const token = `filename*=UTF-8''`;
    let fileName = 'downloaded.pdf';
    for (let val of content) {
        if (val.trim().indexOf(token) === 0) {
            fileName = decodeURIComponent(val.trim().replace(token, ''));
            break;
        }
    }
    return fileName;
}


function downloadData(url, errorMessage, method='GET', formData=null) {
    let foggy = document.getElementById("foggy");
    let fog = foggy !== null;
    if (fog) {
		// foggy.style.height = document.body.scrollHeight + 'px';
        foggy.style.display = 'inline';
    }
    let request = new XMLHttpRequest();
    request.open(method, url, true);
    let csrftoken = getCookie('csrftoken');
    request.setRequestHeader("X-CSRFToken", csrftoken);
    request.setRequestHeader('Content-Type', 'application/x-www-form-urlencoded');
    request.responseType = 'blob';
    request.onload = function (e) {
        processResponse(this, request);
        if (fog) document.getElementById('foggy').style.display = 'none';
        document.getElementById('message_box_text').innerHTML = '';
    };
    request.onreadystatechange = function () {
        // catching errors in the responses, that might cause the disabled status to freeze
        if (request.status !== 200 && request.readyState === 4) {
            window.alert(errorMessage);
        }
    }
    request.send();
}


function mainDownload(url, data, set_foggy=false, method='GET') {
    if (set_foggy) {
        document.getElementById('message_box_text').innerHTML =
            '<i class="fa fa-spinner slow-spin fa-fw" style="color: #0065bd"></i>&nbsp;Preparing Download...';
    }

    let i = 0;
    for (const [attr, choice] of Object.entries(data)) {
        if (i === 0) url += "?" + attr + "=" + choice;
        else url += "&" + attr + "=" + choice;
        i += 1;
    }

    downloadData(
        url, 'Unknown error while downloading substructure data',
        method
    );
}


function downloadSubstructure(url, coefficients) {
    document.getElementById('message_box_text').innerHTML =
        '<i class="fa fa-spinner slow-spin fa-fw" style="color: #0065bd"></i>&nbsp;Downloading substructure data...';
    downloadData(
        url + '?coefficients=' + coefficients,
        'Unknown error while downloading substructure data'
    );
}

function downloadChainLength(url, group) {
    document.getElementById('message_box_text').innerHTML =
        '<i class="fa fa-spinner slow-spin fa-fw" style="color: #0065bd"></i>&nbsp;Downloading chain length data...'
    if (group === null) {
        downloadData(
            url + '?group=' + '',
            'Unknown error while downloading substructure data'
        );
    } else {
        downloadData(
            url + '?group=' + group,
            'Unknown error while downloading substructure data'
        );
    }
}

function downloadClassPlot(url) {
    document.getElementById('message_box_text').innerHTML =
        '<i class="fa fa-spinner slow-spin fa-fw" style="color: #0065bd"></i>&nbsp;Downloading class data...';
    downloadData(
        url + '?dummy=',
        'Unknown error while downloading lipidome summary data'
    );
}

function getLionLists(url) {
    let comparison = '';
    let enrichment_tablinks = document.getElementsByClassName('enrichment-link')
    for (const tablink of enrichment_tablinks) {
        if (tablink.classList.value.indexOf('enrichment-link active') > -1) {
            comparison = tablink.id.split('_ETabLink')[0];
        }
    }
    document.getElementById('message_box_text').innerHTML =
        '<i class="fa fa-spinner slow-spin fa-fw" style="color: #0065bd"></i>&nbsp;Downloading LION target list data...'
    downloadData(
        url + '?comparison=' + comparison,
        'Unknown error while downloading substructure data'
    );
}
