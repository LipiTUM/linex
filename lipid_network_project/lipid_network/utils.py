from .models import (
    UserIds, UploadedData, UploadedGroups,
    UploadedClassSettings, UploadedFaSettings,
    ComputedNetwork, ComputationProgress,
    NetworkEnrichment, UploadedConfirmedSpecies,
    LipidomeSummary, ChainLengthAnalysis,
    SubstructureAnalysis
)
from django.db.models import Max
from django.db import connection
from django.http import (
    JsonResponse, FileResponse, Http404
)
from django.conf import settings
from django.utils.timezone import utc
from wsgiref.util import FileWrapper
from datetime import timedelta
from datetime import datetime
import pathlib
import os
import pandas as pd
import numpy as np
import json
from io import BytesIO


# TIMEOUT for deleting data from server in minutes
TIMEOUT_DELTA_MIN = 500

FILE_TYPES = ["tsv", "csv"]
PATH = settings.USER_DATA
EXAMPLE_PATH = settings.EXAMPLE_DATA
FORM2ATTR = {
    "corr_data": "correlations",
    "pcorr_data": "partial_correlations",
    "fc_data": "fold_changes",
    "pval_data": "pvalues",
}
ATTR_NODES = {
    "corr_data": False,
    "pcorr_data": False,
    "fc_data": True,
    "pval_data": True,
}
# TODO: add remaining from github issue
FORBIDDEN_CHARACTERS = {'&', "'", '"', '+', '/', '?'}


def assign_user_id():
    print("New user, assigning id")
    # session id generation
    tmp = UserIds.objects.aggregate(Max("usercount"))["usercount__max"]
    now = datetime.utcnow().replace(tzinfo=utc)
    if tmp is None:
        count = 2
    else:
        count = tmp + 1
    session_id = f"{count}_{now.strftime('%Y%m%d%H%M')}"
    # user model
    uid = UserIds()
    uid.usercount = count
    uid.timestamp = now
    uid.userid = session_id
    uid.save()
    return session_id


def check_user_id(request):
    session_id = request.session.get("session_id", 1)
    if session_id == 1 or UserIds.objects.filter(userid=session_id).count() == 0:
        request.session["session_id"] = assign_user_id()
        print(f"new user id assigned: {request.session['session_id']}")
    else:
        print(f"new request from user: {session_id} - {UserIds.objects.get(userid=session_id).timestamp}")


def check_running_task(session_id: str) -> bool:
    if UserIds.objects.filter(userid=session_id).count() == 0:
        return False
    else:
        userid_entry = UserIds.objects.get(userid=session_id)
        return userid_entry.running_task


def set_running_task(session_id: str, value: bool):
    if UserIds.objects.filter(userid=session_id).count() == 1:
        userid_entry = UserIds.objects.get(userid=session_id)
        userid_entry.running_task = value
        userid_entry.save()


def table_exists(name):
    return f"lipid_network_{name}" in connection.introspection.table_names()


def delete_db_by_time():
    print(f"time: {datetime.utcnow()}")
    print(f"timedelta: {datetime.utcnow() - timedelta(minutes=2)}")
    # Delete files
    if table_exists("uploadeddata"):
        UploadedData.objects.filter(
            timestamp__lt=datetime.utcnow().replace(tzinfo=utc) - timedelta(minutes=TIMEOUT_DELTA_MIN)).delete()
    if table_exists("uploadedgroups"):
        UploadedGroups.objects.filter(
            timestamp__lt=datetime.utcnow().replace(tzinfo=utc) - timedelta(minutes=TIMEOUT_DELTA_MIN)).delete()
    if table_exists("uploadedclasssettings"):
        UploadedClassSettings.objects.filter(
            timestamp__lt=datetime.utcnow().replace(tzinfo=utc) - timedelta(minutes=TIMEOUT_DELTA_MIN)).delete()
    if table_exists("uploadedfasettings"):
        UploadedFaSettings.objects.filter(
            timestamp__lt=datetime.utcnow().replace(tzinfo=utc) - timedelta(minutes=TIMEOUT_DELTA_MIN)).delete()
    if table_exists("uploadedconfirmedspecies"):
        UploadedConfirmedSpecies.objects.filter(
            timestamp__lt=datetime.utcnow().replace(tzinfo=utc) - timedelta(minutes=TIMEOUT_DELTA_MIN)).delete()
    # Delete networks
    if table_exists("computednetwork"):
        ComputedNetwork.objects.filter(
            timestamp__lt=datetime.utcnow().replace(tzinfo=utc) - timedelta(minutes=TIMEOUT_DELTA_MIN)).delete()
    # Progress
    if table_exists("computationprogress"):
        ComputationProgress.objects.filter(
            timestamp__lt=datetime.utcnow().replace(tzinfo=utc) - timedelta(minutes=TIMEOUT_DELTA_MIN)).delete()
    # Enrichment
    if table_exists("networkenrichment"):
        NetworkEnrichment.objects.filter(
            timestamp__lt=datetime.utcnow().replace(tzinfo=utc) - timedelta(minutes=TIMEOUT_DELTA_MIN)).delete()
        # substructure
    if table_exists("substructureanalysis"):
        SubstructureAnalysis.objects.filter(
            timestamp__lt=datetime.utcnow().replace(tzinfo=utc) - timedelta(minutes=TIMEOUT_DELTA_MIN)).delete()
    # summary
    if table_exists("lipidomesummary"):
        LipidomeSummary.objects.filter(
            timestamp__lt=datetime.utcnow().replace(tzinfo=utc) - timedelta(minutes=TIMEOUT_DELTA_MIN)).delete()
    # chainlength
    if table_exists("chainlengthanalysis"):
        ChainLengthAnalysis.objects.filter(
            timestamp__lt=datetime.utcnow().replace(tzinfo=utc) - timedelta(minutes=TIMEOUT_DELTA_MIN)).delete()


def delete_file_entries():
    for file in os.listdir(path=PATH):
        if file != "dummy.txt" and file != "example_data.zip" and \
                not file.endswith(".log"):
            fname = pathlib.Path(os.path.join(PATH, file))
            if fname.exists():
                mtime = datetime.fromtimestamp(fname.stat().st_mtime)
                if mtime < datetime.now() - timedelta(minutes=TIMEOUT_DELTA_MIN):
                    print(f"deleting file {os.path.join(PATH, file)}")
                    os.remove(os.path.join(PATH, file))


def delete_db_by_user(session_id):
    if table_exists("uploadeddata"):
        UploadedData.objects.filter(userid=session_id).delete()
    if table_exists("uploadedgroups"):
        UploadedGroups.objects.filter(userid=session_id).delete()
    if table_exists("uploadedclasssettings"):
        UploadedClassSettings.objects.filter(userid=session_id).delete()
    if table_exists("uploadedfasettings"):
        UploadedFaSettings.objects.filter(userid=session_id).delete()
    if table_exists("uploadedconfirmedspecies"):
        UploadedConfirmedSpecies.objects.filter(userid=session_id).delete()
    # Delete networks
    if table_exists("computednetwork"):
        ComputedNetwork.objects.filter(userid=session_id).delete()
    # Progress
    if table_exists("computationprogress"):
        ComputationProgress.objects.filter(userid=session_id).delete()
    # Enrichment
    if table_exists("networkenrichment"):
        NetworkEnrichment.objects.filter(userid=session_id).delete()
    # substructure
    if table_exists("substructureanalysis"):
        SubstructureAnalysis.objects.filter(userid=session_id).delete()
    # summary
    if table_exists("lipidomesummary"):
        LipidomeSummary.objects.filter(userid=session_id).delete()
    # chainlength
    if table_exists("chainlengthanalysis"):
        ChainLengthAnalysis.objects.filter(userid=session_id).delete()


def check_user_files(session_id):
    file_start = f"uploaded_file_{session_id}"
    group_file_start = f"uploaded_group_file_{session_id}"
    class_file_start = f"uploaded_class_file_{session_id}"
    fa_file_start = f"uploaded_fa_file_{session_id}"
    confirmed_file_start = f"uploaded_confirmed_species_file_{session_id}"
    for file in os.listdir(path=PATH):
        if file.startswith(file_start) or \
                file.startswith(group_file_start) or \
                file.startswith(class_file_start) or \
                file.startswith(fa_file_start) or \
                file.startswith(confirmed_file_start):
            return True
    return False


def delete_user_files(session_id):
    file_start = f"uploaded_file_{session_id}"
    group_file_start = f"uploaded_group_file_{session_id}"
    class_file_start = f"uploaded_class_file_{session_id}"
    fa_file_start = f"uploaded_fa_file_{session_id}"
    confirmed_file_start = f"uploaded_confirmed_species_file_{session_id}"
    vis_network_file = f"vis_network_{session_id}.pickle"
    network_file = f"lipid_network_{session_id}.pickle"
    for file in os.listdir(path=PATH):
        if file.startswith(file_start) or \
                file.startswith(group_file_start) or \
                file.startswith(class_file_start) or \
                file.startswith(fa_file_start) or \
                file.startswith(confirmed_file_start) or \
                file.startswith(vis_network_file) or \
                file.startswith(network_file):
            os.remove(os.path.join(PATH, file))


def check_file_type(filename):
    for ftype in FILE_TYPES:
        if filename.lower().endswith(ftype):
            return True
    return False


def init_model(request, model):
    session = request.session["session_id"]
    mod = model()
    mod.userid = session
    if isinstance(mod, UploadedGroups):
        # removing forbidden characters => characters interfering with
        # data transfer front end link generation etc.
        group_file_str = request.FILES["group_docfile"].file.read().decode('utf-8')
        contains_forbidden = False
        for char in FORBIDDEN_CHARACTERS:
            if not contains_forbidden:
                if char in group_file_str:
                    contains_forbidden = True
                    group_file_str = group_file_str.replace(char, '')
            else:
                group_file_str = group_file_str.replace(char, '')
        request.FILES["group_docfile"].file = BytesIO(
            bytes(group_file_str, 'utf-8')
        )
        mod.group_docfile = request.FILES["group_docfile"]
        mod.forbidden_characters = contains_forbidden
    elif isinstance(mod, UploadedClassSettings):
        mod.class_settings = request.FILES["class_file"]
    if isinstance(mod, UploadedFaSettings):
        mod.fa_settings = request.FILES["fa_file"]
    if isinstance(mod, UploadedConfirmedSpecies):
        mod.confirmed_species_file = request.FILES["confirmed_species_file"]
    else:
        mod.docfile = request.FILES["docfile"]
    mod.timestamp = datetime.utcnow().replace(tzinfo=utc)
    return mod


def checkbox_on(qdict, name, invert=False):
    if invert:
        return qdict.get(name) != "on"
    else:
        return qdict.get(name) == "on"


class NpEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, np.nan):
            return "nan"
        return super(NpEncoder, self).default(obj)


def to_json_series(data, safe=False, content_type="application/json", **kwargs):
    jr = JsonResponse(data, safe=safe,
                      content_type=content_type,
                      **kwargs)
    return jr.content.decode("utf-8")


def save_file_response(file, content_type, **kwargs):
    if os.path.exists(file):
        response = FileResponse(
            open(file, 'rb'),
            content_type=content_type,
            **kwargs
        )
        response['Content-Disposition'] = \
            f"attachment; filename={os.path.basename(file)}"

        os.remove(file)
        return response
    raise Http404


def delete_my_data(userid):
    """
    Delete all data from one user.
    """
    delete_db_by_user(userid)
    delete_user_files(userid)


def write_attribute_data(userid, network, form_attr):
    attr = FORM2ATTR[form_attr]
    if ATTR_NODES[form_attr]:
        data_ = network.lipid_attributes[attr]
        if isinstance(data_, dict):
            data_ = pd.DataFrame(
                {group: data for group, data in data_.items()}
            )
    else:
        data_ = network.interaction_attributes[attr]
        if isinstance(data_, dict):
            files = []
            for group, data in data_.items():
                file = os.path.join(PATH, f"{attr}_{group}_data_{userid}.csv")
                if data is not None:
                    data.to_csv(file)
                    files.append(file)
            return files
    if hasattr(data_, "to_csv"):
        file = os.path.join(PATH, f"{attr}_data_{userid}.csv")
        data_.to_csv(file)
        return file
    # this should mean 'None', but any 'unknown' data types
    # are not covered here
    return None


def save_node_metrics(network, file):
    pd.DataFrame(
        {
            "Degree": network.lipid_attributes["degree"],
            "Betweenness Centrality": network.lipid_attributes["betweenness"],
            "Closeness Centrality": network.lipid_attributes["closeness"]
        }
    ).to_csv(file)
