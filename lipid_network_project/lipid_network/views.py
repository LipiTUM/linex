import warnings
import pickle
from django.shortcuts import redirect
from .forms import (
    DataForm, NetworkDownload,
    DataProcessOptions, StatisticsOptions,
    NetworkOptions, StatisticsDownload,
    LynxDownload, Version2Specific,
    JSONDownload, ContributionForm
)
from django.shortcuts import render
from django.template.loader import render_to_string
from django.http import HttpResponse, JsonResponse, HttpResponseServerError
from json import JSONDecodeError, load
import networkx as nx
import datetime
from django.utils.timezone import utc
from django.conf import settings
from .models import (
    UploadedGroups, UploadedData,
    UploadedClassSettings, UploadedFaSettings,
    ComputedNetwork, ComputationProgress,
    UserIds, NetworkEnrichment, UploadedConfirmedSpecies,
    UserContribution, LipidomeSummary, SubstructureAnalysis,
    ChainLengthAnalysis
)
from .utils import (
    check_user_id, delete_db_by_time, delete_file_entries,
    check_user_files, delete_user_files, check_file_type,
    delete_db_by_user, init_model, to_json_series,
    save_file_response, checkbox_on, TIMEOUT_DELTA_MIN,
    delete_my_data, PATH, write_attribute_data,
    FORM2ATTR, save_node_metrics, NpEncoder, check_running_task,
    set_running_task
)
from .pdf_utils import plot_network_pdf, plot_enrichment_pdfs
from .view_utils import (
    get_network_view_as_dict,
    check_version2_network,
    build_function_call,
    enrichment_key,
    zip_files,
    bool_from_url_query,
    check_get_summary_model,
    save_summary_files
)
from .tasks import compute_network_1, compute_network_2, run_enrichment
import os
import json
from itertools import combinations
import base64

now = datetime.datetime.utcnow().replace(tzinfo=utc)

# ###
# Version number:
#   - (Please update with every commit)
#   - First number: Major version e.g. new publication
#   - Second number: Feature additions
#   - Third number: Bug fixes
LINEX_VERSION = "2.4.1"


# ###


def index(request):
    delete_db_by_time()
    delete_file_entries()
    check_user_id(request)
    return render(request, 'lipid_network/index.html')


def about(request):
    delete_db_by_time()
    delete_file_entries()
    check_user_id(request)
    return render(
        request, 'lipid_network/about.html',
        context={"timeout": TIMEOUT_DELTA_MIN,
                 "linex_version": LINEX_VERSION}
    )


def upload(request):
    delete_db_by_time()
    delete_file_entries()
    check_user_id(request)
    userid = UserIds.objects.get(userid=request.session['session_id'])
    # Handle file upload
    groups = False
    message = ""
    version = int(getattr(request, request.method).get('version', 0))
    print(f'Selected Version: {version}')
    if version > 2:
        message = f'Illegal version number {version}. Please choose between 1 and 2'
    if request.method == 'POST':
        data_form = DataForm(request.POST, request.FILES)
        userid.version = int(request.POST.get('version_number', 0))
        if userid.version == 0:
            # FIXME: proper error handling
            message = 'Illegal version number 0. Please choose between 1 and 2'
        elif userid.version > 2:
            # FIXME: proper error handling
            message = f'Illegal version number {userid.version}. Please choose between 1 and 2'
        else:
            userid.save()
        data_process = DataProcessOptions(request.POST)
        statistics = StatisticsOptions(request.POST)
        netopt = NetworkOptions(request.POST)
        if userid.version == 2:
            version2spec = Version2Specific(request.POST)
        else:
            version2spec = Version2Specific()
        # TODO: add a check for group_form in not empty
        class_settings = None
        fa_settings = None
        confirmed_species_filename = None
        if check_running_task(request.session['session_id']):
            message = 'You already have a computation running on the server ' \
                      '(Network generation or enrichment). ' \
                      'Please wait for the task to finish.'
        if message == "":
            if data_form.is_valid():
                delete_db_by_user(request.session["session_id"])
                if check_user_files(request.session["session_id"]):
                    delete_user_files(request.session["session_id"])
                newdoc = init_model(request, UploadedData)
                newdoc.save()
                if not os.path.exists(newdoc.docfile.path):
                    message = f"Data file '{newdoc.docfile}' could not be " \
                              "uploaded properly. Please try to remove " \
                              "special characters from the file and filename."
                    context = {
                        'data_form': data_form,
                        'statistics': statistics,
                        'data_process': data_process,
                        'netopt': netopt,
                        'groups': groups,
                        'version2specific': version2spec,
                        'message': message,
                        'upload_version': version
                    }
                    return render(
                        request, 'lipid_network/upload.html',
                        context=context)
                forbidden_chars = False
                if 'group_docfile' in request.FILES.keys():
                    groups = True
                    newgroup = init_model(request, UploadedGroups)
                    newgroup.save()
                    forbidden_chars = newgroup.forbidden_characters
                    check_type = check_file_type(str(newdoc.docfile)) and \
                                 check_file_type(str(newgroup.group_docfile))
                else:
                    check_type = check_file_type(str(newdoc.docfile))
                if 'fa_file' in request.FILES.keys():
                    fa_sets = init_model(request, UploadedFaSettings)
                    fa_sets.save()
                    fa_settings = f"uploaded_fa_file_{request.session['session_id']}.txt"
                if check_type:
                    if ComputationProgress.objects.filter(
                            userid=request.session["session_id"]).count() == 1:
                        ComputationProgress.objects.filter(
                            userid=request.session["session_id"]).delete()
                    progress = ComputationProgress()
                    progress.userid = request.session["session_id"]
                    progress.timestamp = userid.timestamp
                    progress.reading = "Waiting"
                    progress.groups = groups
                    progress.save()
                    if userid.version == 1:
                        if 'class_file' in request.FILES.keys():
                            class_sets = init_model(request,
                                                    UploadedClassSettings)
                            class_sets.save()
                            class_settings = f"uploaded_class_file_{request.session['session_id']}.txt"

                        # Set running task for user to True
                        set_running_task(request.session["session_id"],
                                         value=True)

                        compute_network_1(
                            request.session["session_id"],
                            has_groups=groups, repeat_until=None,
                            lipid_level=request.POST.get("lipid_level",
                                                         "molecular"),
                            # TODO: tests save setting files
                            lipid_settings=class_settings,
                            sample_axis="index" if checkbox_on(request.POST,
                                                               "lipids_as_columns") else "columns",
                            fa_settings=fa_settings,
                            is_in_lynx_format=checkbox_on(request.POST,
                                                          "lynx_format", True),
                            compute_correlations=checkbox_on(request.POST,
                                                             "corrs"),
                            compute_partial_correlations=checkbox_on(
                                request.POST, "pcorrs"),
                            compute_fold_changes=checkbox_on(request.POST,
                                                             "fcs"),
                            compute_pvalues=checkbox_on(request.POST, "pvals"),
                            test_method=request.POST.get("ptest", "wilcoxon"),
                            directed=checkbox_on(request.POST, "directed"),
                            log_data=checkbox_on(request.POST, "log_data"),
                            to_log=checkbox_on(request.POST, "log_ratios"),
                            reference_group=request.POST.get("ref_group"),
                            significance_threshold=float(
                                request.POST.get("corr_sig")),
                            forbidden_group_chars=forbidden_chars
                        )
                        # return upload_pending(request, version=1, redirected=True)
                        return redirect('upload-pending')
                    elif userid.version == 2:
                        data_form.fields['group_docfile'].required = True
                        data_form.fields['group_docfile'].label = \
                            "Select a sample group file"
                        if 'confirmed_species_file' in request.FILES.keys():
                            confirmed_species = init_model(request,
                                                           UploadedConfirmedSpecies)
                            confirmed_species.save()
                            confirmed_species_filename = f"uploaded_confirmed_species_file_{request.session['session_id']}.txt"
                        db_str = request.POST.get("database")
                        databases = db_str if db_str != "both" else (
                        'Rhea', 'Reactome')

                        if ComputationProgress.objects.filter(
                            userid=request.session["session_id"]).count() == 1:
                            ComputationProgress.objects.filter(
                                userid=request.session["session_id"]).delete()
                        progress = ComputationProgress()
                        progress.userid = request.session["session_id"]
                        progress.timestamp = userid.timestamp
                        progress.reading = "Waiting"
                        progress.groups = groups
                        progress.version = 2
                        progress.save()

                        # Set running task for user to True
                        set_running_task(request.session["session_id"],
                                         value=True)

                        compute_network_2(
                            request.session["session_id"],
                            has_groups=groups, repeat_until=None,
                            lipids_are_columns=checkbox_on(request.POST,
                                                           "lipids_as_columns"),
                            fa_settings=fa_settings,
                            is_in_lynx_format=checkbox_on(request.POST,
                                                          "lynx_format", True),
                            compute_correlations=checkbox_on(request.POST,
                                                             "corrs"),
                            compute_partial_correlations=checkbox_on(
                                request.POST, "pcorrs"),
                            compute_fold_changes=checkbox_on(request.POST,
                                                             "fcs"),
                            compute_pvalues=checkbox_on(request.POST, "pvals"),
                            test_method=request.POST.get("ptest", "wilcoxon"),
                            log_data=checkbox_on(request.POST, "log_data"),
                            to_log=checkbox_on(request.POST, "log_ratios"),
                            reference_group=request.POST.get("ref_group"),
                            significance_threshold=float(
                                request.POST.get("corr_sig")),
                            organism=request.POST.get("organism"),
                            database=databases,
                            ether_conversions=checkbox_on(request.POST,
                                                          "ether_conversions"),
                            forbidden_group_chars=forbidden_chars,
                            confirmed_species_file=confirmed_species_filename
                        )
                return redirect('upload-pending')
            else:
                if userid.version == 2:
                    data_form.group_docfile.required = True
                    data_form.group_docfile.label = \
                        "Select a sample group file"
                context = {
                    'data_form': data_form,
                    'statistics': statistics,
                    'data_process': data_process,
                    'netopt': netopt,
                    'groups': groups,
                    'version2specific': version2spec,
                    'message': "Upload Failed! Please make sure all uploaded data has the correct format and file type",
                    'upload_version': int(
                        request.POST.get('version_number', 0))
                }
                return render(request, 'lipid_network/upload.html', context)
        else:
            pass
    else:
        data_form = DataForm()
        data_process = DataProcessOptions()
        statistics = StatisticsOptions()
        netopt = NetworkOptions()
        version2spec = Version2Specific()

        # Example data execution ---------------------------------------------------------------------------------------
        exampledata_check = int(request.GET.get('exampledata', 0))
        if exampledata_check == 1:
            if not check_running_task(request.session['session_id']):
                # removing user entries from db
                delete_db_by_user(request.session["session_id"])
                # initiate progress
                progress = ComputationProgress()
                progress.userid = request.session["session_id"]
                progress.timestamp = userid.timestamp
                progress.reading = "Waiting"
                progress.groups = groups
                progress.version = 2
                progress.save()

                # Start task

                # Set running task for user to True
                set_running_task(request.session["session_id"], value=True)

                compute_network_2(
                    request.session["session_id"],
                    has_groups=True,
                    repeat_until=None,
                    lipids_are_columns=True,
                    fa_settings=None,
                    is_in_lynx_format=True,
                    compute_correlations=True,
                    compute_partial_correlations=False,
                    compute_fold_changes=True,
                    compute_pvalues=True,
                    test_method="ttest",
                    log_data=False,
                    to_log=True,
                    reference_group=None,
                    significance_threshold=0.05,
                    organism="HSA",
                    database=('Rhea', 'Reactome'),
                    ether_conversions=True,
                    forbidden_group_chars=False,
                    confirmed_species_file=None,
                    use_example_data=True
                )

                # Return to upload-pending
                return redirect('upload-pending')
                # return upload_pending(request, version=2, redirected=True)

            else:
                message = 'You already have a computation running on the server ' \
                          '(Network generation or enrichment). ' \
                          'Please wait for the task to finish.'
                version = 2

                data_form.group_docfile.required = True
                data_form.group_docfile.label = "Select a sample group file"

                context = {
                    'data_form': data_form,
                    'statistics': statistics,
                    'data_process': data_process,
                    'netopt': netopt,
                    'groups': groups,
                    'version2specific': version2spec,
                    'message': message,
                    'upload_version': version,
                    'pending_link': True
                }

                return render(request, 'lipid_network/upload.html',
                              context=context)

    if version == 2:
        data_form.fields['group_docfile'].required = True
        data_form.fields['group_docfile'].label = "Select a sample group file"
    if hasattr(data_form, "version_number"):
        version = data_form.version_number
    context = {
        'data_form': data_form,
        'statistics': statistics,
        'data_process': data_process,
        'netopt': netopt,
        'groups': groups,
        'version2specific': version2spec,
        'message': message,
        'upload_version': version
    }
    # Render list page with the documents and the data_form
    # TODO: split option form into different sections!
    if version == 0:
        return render(request, 'lipid_network/upload.html',
                      {'upload_version': version,
                       'message': message})
    return render(request, 'lipid_network/upload.html',
                  context=context)


def upload_pending(request, version=None, redirected=False):
    if request.method == 'POST' and not redirected:
        if ComputationProgress.objects.filter(
            userid=request.session["session_id"]).count() == 1:
            progress = ComputationProgress.objects.get(
                userid=request.session["session_id"])
            progress_dict = {
                'reading': progress.reading,
                'lynx': progress.lynx,
                'class_init': progress.class_init,
                'compute_network': progress.compute_network,
                'compute_statistics': progress.compute_statistics,
                'compute_views': progress.compute_views,
                'compute_summaries': progress.compute_summaries,
                'done': progress.done,
                'error': progress.error,
                'version': progress.version if version is None else version
            }
            if progress.unconverted:
                progress_dict['unconverted'] = progress.unconverted
            if progress.message:
                progress_dict['message'] = progress.message
            if progress.warning:
                progress_dict['warning'] = progress.warning
        else:
            progress_dict = {}
        return JsonResponse(progress_dict, safe=False)
    delete_db_by_time()
    delete_file_entries()
    check_user_id(request)
    if ComputationProgress.objects.filter(
        userid=request.session["session_id"]).count() == 1:
        progress = ComputationProgress.objects.get(
            userid=request.session["session_id"])
        context = {
            "pending": True, "timeout": TIMEOUT_DELTA_MIN,
            "progress": progress,
            "version": progress.version if version is None else version
        }
    else:
        context = {"pending": False, "timeout": TIMEOUT_DELTA_MIN,
                   "version": -1 if version is None else version}
    return render(request, 'lipid_network/upload-pending.html',
                  context=context)


def analysis(request, to_string: bool = False):
    # checking data availability
    delete_db_by_time()
    delete_file_entries()
    check_user_id(request)
    if ComputedNetwork.objects.filter(
        userid=request.session["session_id"]).count() != 1:
        return render(
            request, "lipid_network/analysis.html",
            context={"missing": True, "timeout": TIMEOUT_DELTA_MIN,
                     "host": request.get_host()}
        )
    network_entry = ComputedNetwork.objects.get(
        userid=request.session["session_id"])
    network = pickle.load(
        open(os.path.join(PATH, network_entry.network), "rb"))
    if network_entry.network == "":
        return render(
            request, "lipid_network/analysis.html",
            context={"missing": True, "timeout": TIMEOUT_DELTA_MIN,
                     "host": request.get_host()}
        )
    # actual analysis page rendering
    if network_entry.vis_network == "":
        if network_entry.version == 2:
            raise ValueError(
                'Something went wrong in the background task for version 2'
            )
        else:
            native_context, network = get_network_view_as_dict(
                network_entry, network, 1
            )
            vis_network_file = os.path.join(
                PATH, f"vis_newtork_{request.session['session_id']}.pickle"
            )
            pickle.dump(native_context, open(vis_network_file, "wb"))
            network_entry.vis_network = vis_network_file
            network_entry.save()
    vis_network = pickle.load(
        open(os.path.join(PATH, network_entry.vis_network), "rb")
    )
    if NetworkEnrichment.objects.filter(
        userid=request.session["session_id"]).count() == 1:
        enrichment = NetworkEnrichment.objects.get(
            userid=request.session["session_id"])
        if not vis_network.get('enrichments', {}):
            vis_network['enrichments'] = to_json_series(
                enrichment.computed_enrichments)
    else:
        vis_network['enrichments'] = {}

    # pre-computing groups in the format of network enrichment
    if not vis_network.get('enrichment_groups', {}):
        vis_network['enrichment_groups'] = {
            enrichment_key({'group1': group1, 'group2': group2})
            for group1, group2 in combinations(vis_network['single_groups'], 2)
        }

    vis_network['single_groups'] = list(vis_network['single_groups'])

    if to_string:
        # called when rendering to external .html file
        css = ''
        for file in os.listdir(os.path.join(settings.STATIC_DIR, 'css')):
            css += open(os.path.join(settings.STATIC_DIR, 'css', file),
                        'r').read()
        js = ''
        for file in os.listdir(os.path.join(settings.STATIC_DIR, 'js')):
            js += open(os.path.join(settings.STATIC_DIR, 'js', file),
                       'r').read()

        vis_network['css'] = css
        vis_network['js'] = js
        vis_network['to_string'] = True
        vis_network['host'] = request.get_host()
        encoded_string = base64.b64encode(
            open(os.path.join(settings.STATIC_DIR,
                              'images/logo_linex_small.png'), "rb").read()
        ).decode('ascii')
        vis_network['logo'] = encoded_string
        return render_to_string(
            os.path.join(settings.TEMPLATE_DIR, 'lipid_network',
                         'analysis.html'),
            context=vis_network
        )
    else:
        vis_network['to_string'] = False

    if network_entry.version == 2:
        # adding variables from substructure
        substruct = SubstructureAnalysis.objects.get(
            userid=request.session["session_id"])
        for key, vals in substruct.plot.items():
            vis_network[key] = vals

        # adding variables from lipidome summary
        summary = LipidomeSummary.objects.get(
            userid=request.session["session_id"])
        for key, vals in summary.data.items():
            vis_network[key] = vals

        # adding variables from chain-length analysis
        chain_length = ChainLengthAnalysis.objects.get(
            userid=request.session["session_id"])
        vis_network['chainlength_analysis'] = chain_length.plot

    return render(
        request, 'lipid_network/analysis.html',
        context=vis_network
    )


def pdf_download(request, network_data: dict = None, fs=7):
    # downloading pdf view
    if request.method == "POST":
        if network_data is None:
            try:
                network_data = json.loads(request.body)
            except JSONDecodeError:
                return HttpResponse(status=500)
        # ALL users can download, regardless of whether they are
        # currently using the web-tool. Necessary to avoid file
        # confusions
        check_user_id(request)
        if network_data['selectedType'] == "enrichment":
            file_base = f"LINEX_{request.session['session_id']}_enrichment_"
            pdf_files = plot_enrichment_pdfs(network_data,
                                             os.path.join(PATH, file_base),
                                             "pdf", fs=fs)
            # TODO return all files as a zip folder
            if len(pdf_files) == 1:
                return save_file_response(pdf_files[0], "application/pdf")
            else:
                zip_file = f"LINEX_enrichments_{request.session['session_id']}.zip"
                zip_files(zip_file, pdf_files)
                return save_file_response(zip_file,
                                          'application/force-download')
        else:
            file = f"LINEX_{request.session['session_id']}_network_view.pdf"

            pdf_file = os.path.join(PATH, file)
            plot_network_pdf(network_data, pdf_file, "pdf", fs=fs)

            with open(pdf_file, "rb") as pdf:
                response = HttpResponse(
                    pdf.read(), content_type="application/pdf")
                response['Content-Disposition'] = \
                    f"attachment; filename={os.path.basename(file)}"
                return response


def substructure_download(request):
    if request.GET:
        substruct = check_get_summary_model(
            request.session["session_id"], SubstructureAnalysis)
        if isinstance(substruct, HttpResponse):
            # this means error
            return substruct
        file = f"LINEX_{request.session['session_id']}_substructure_"
        if request.GET['coefficients'] == 'true':
            file = f'{file}coefficients.csv'
            substruct.data['coefficients'].to_csv(os.path.join(PATH, file))
        else:
            file = f'{file}data.csv'
            substruct.data['data'].to_csv(os.path.join(PATH, file))
        return save_file_response(os.path.join(PATH, file),
                                  "application/force-download")


def chainlength_download(request):
    if request.GET:
        cl_entry = check_get_summary_model(
            request.session["session_id"], ChainLengthAnalysis)
        if isinstance(cl_entry, HttpResponse):
            # this means error
            return cl_entry
        file = f"LINEX_{request.session['session_id']}_chainlength"
        files = []
        if request.GET['group'] != '':
            group = request.GET['group']
            for lclass, class_data in cl_entry.data[group].items():
                class_file = os.path.join(PATH, f'{file}_{group}_{lclass}.csv')
                class_data.to_csv(class_file)
                files.append(class_file)
        else:
            for lclass, class_data in cl_entry.data.items():
                class_file = os.path.join(PATH, f'{file}_{lclass}.csv')
                class_data.to_csv(class_file)
                files.append(class_file)
        zip_file = os.path.join(PATH, f"{file}.zip")
        zip_files(zip_file, files)
        return save_file_response(zip_file, "application/force-download")


def download_summary_data(request):
    if request.GET:
        userid = request.session["session_id"]
        network_entry = check_get_summary_model(
            userid, ComputedNetwork, "network")

        if isinstance(network_entry, HttpResponse):
            # this means error
            return network_entry

        network = pickle.load(
            open(os.path.join(PATH, network_entry.network), "rb"))
        with warnings.catch_warnings():
            warnings.simplefilter(action='ignore', category=FutureWarning)
            class_sum, length_sum, db_sum, c_db_sum = \
                network.lipidome_summary()
        # saving class summary
        files = [os.path.join(PATH, f"LINEX_{userid}_class_summary.csv")]
        class_sum.to_csv(files[0])
        # saving other summary files
        files += save_summary_files(
            length_sum, PATH, f"LINEX_{userid}_length_summary.csv")
        files += save_summary_files(
            db_sum, PATH, f"LINEX_{userid}_db_summary.zip")
        files += save_summary_files(
            c_db_sum, PATH, f"LINEX_{userid}_c_db_summary.csv")
        # saving to zip-files
        zip_file = os.path.join(PATH, f"LINEX_{userid}_lipidome_summary.zip")
        zip_files(zip_file, files)
        return save_file_response(zip_file, "application/force-download")


def change_network_type(request):
    # checking data availability
    delete_db_by_time()
    delete_file_entries()
    check_user_id(request)
    network_entry = check_version2_network(request)
    if not isinstance(network_entry, ComputedNetwork):
        return network_entry
    if int(network_entry.version) != 2:
        message = 'Networks computed with version 1 cannot change their network type.' \
                  'Please compute a version 2 network'
        return JsonResponse(
            {'error': True, 'message': message}, safe=False)
    try:
        network_type = request.GET.get('network_type')
        vis_network = pickle.load(open(network_entry.vis_network, "rb"))
        vis_context = vis_network.get(network_type)
        if vis_context is None:
            return HttpResponseServerError()
        print(vis_context)
        new_data = {
            'nodes': to_json_series(vis_context['nodes'], encoder=NpEncoder),
            'edges': to_json_series(vis_context['edges'], encoder=NpEncoder),
            'error': False
        }
        return JsonResponse(new_data, safe=True)
    except Exception as e:
        # TODO: error message as json
        print(e)
        message = f'Unknown internal error while computing network\n{e}'
        return JsonResponse(
            {'error': True, 'message': message}, safe=False
        )


def download(request):
    delete_db_by_time()
    delete_file_entries()
    check_user_id(request)
    has_lynx = False
    not_computed = True
    network_form = NetworkDownload()
    stats_form = StatisticsDownload()
    lynx_form = LynxDownload()
    json_form = JSONDownload()
    if request.GET:
        if ComputedNetwork.objects.filter(
                userid=request.session["session_id"]).count() != 1:
            return render(
                request, 'lipid_network/download.html',
                context={'has_lynx': has_lynx,
                         'not_computed': True,
                         'network_form': network_form,
                         'stats_form': stats_form,
                         'lynx_form': lynx_form}
            )
        # TODO: use data streams instead of writing to disk
        not_computed = False
        network_entry = ComputedNetwork.objects.get(
            userid=request.session["session_id"])
        network = pickle.load(open(network_entry.network, "rb"))
        multiple = sum(
            [bool_from_url_query(q) for q in request.GET.values()]) > 1
        multi = []
        # Saving the network
        if bool_from_url_query(request.GET.get("network_html")) or \
                bool_from_url_query(request.GET.get("colour_legend")) or \
                bool_from_url_query(request.GET.get("size_legend")):
            vis_net = network.dynamic_network(
                node_colour_attributes=network_entry.node_colours,
                edge_colour_attributes=network_entry.edge_colours,
                node_size_attributes=network_entry.node_sizes,
                directed=network_entry.directed, random_seed=42)
        else:
            vis_net = None
        if bool_from_url_query(request.GET.get("network_html")):
            # NOTE: this code saves as external file
            file = os.path.join(
                PATH, f"LINEX_network_{request.session['session_id']}.html")
            net_content = analysis(request, to_string=True)
            with open(file, 'w') as html_file:
                html_file.write(net_content)
            if not multiple:
                return save_file_response(file, "text/html")
            else:
                multi.append(file)
        if bool_from_url_query(request.GET.get("network_graphml")):
            file = os.path.join(
                PATH, f"LINEX_network_{request.session['session_id']}.graphml")
            try:
                nx.write_graphml(network.network, file)
            except nx.exception.NetworkXError or KeyError:
                # NOTE: this happens when we have a 'None' in the node/edge
                # data attributes
                for node, data in network.network.nodes(data=True):
                    for attr, value in data.items():
                        if value is None:
                            network.network.nodes[node][attr] = ''
                        elif isinstance(value, list) or isinstance(value,
                                                                   tuple):
                            network.network.nodes[node][attr] = ' ;'.join(
                                value)
                        elif isinstance(value, dict):
                            network.network.nodes[node][attr] = str(value)
                for src, tgt, data in network.network.edges(data=True):
                    for attr, value in data.items():
                        if value is None:
                            network.network.edges[(src, tgt)][attr] = ''
                        elif isinstance(value, list) or isinstance(value,
                                                                   tuple):
                            network.network.edges[(src, tgt)][
                                attr] = ' ;'.join(value)
                        elif isinstance(value, dict):
                            network.network.edges[(src, tgt)][attr] = str(
                                value)
                nx.write_graphml(network.network, file)
            if not multiple:
                return save_file_response(file, "text/plain")
            else:
                multi.append(file)
        # Saving network legends
        if bool_from_url_query(request.GET.get("colour_legend")):
            file = os.path.join(PATH,
                                f"LINEX_colour_legend_{request.session['session_id']}.graphml")
            nx.write_graphml(vis_net.legend_to_networkx(), file)
            if not multiple:
                return save_file_response(file, "text/plain")
            else:
                multi.append(file)
        if bool_from_url_query(request.GET.get("size_legend")):
            file = os.path.join(PATH,
                                f"LINEX_size_legend_{request.session['session_id']}.graphml")
            nx.write_graphml(vis_net.legend_to_networkx(legend="size"), file)
            if not multiple:
                return save_file_response(file, "text/plain")
            else:
                multi.append(file)
        # lipid class connection networks
        if bool_from_url_query(request.GET.get("connections_plot")):
            if network_entry.version == 1:
                file = os.path.join(PATH,
                                    f"LINEX_connections_plot_{request.session['session_id']}.png")
                network.show_class_connections(savepath=file, format="png")
                if not multiple:
                    return save_file_response(file, "image/png")
                else:
                    multi.append(file)
        if bool_from_url_query(request.GET.get("connections_graph")):
            if network_entry.version == 1:
                file = os.path.join(PATH,
                                    f"LINEX_connections_graph_{request.session['session_id']}.graphml")
                class_net = network.show_class_connections(savepath=None,
                                                           show=False,
                                                           return_graph=True)
                nx.write_graphml(class_net, file)
                if not multiple:
                    return save_file_response(file, "text/plain")
                else:
                    multi.append(file)
        # saving lipid lynx data
        if bool_from_url_query(request.GET.get("lynx_data")):
            file = os.path.join(
                PATH, f"LINEX_lynx_data_{request.session['session_id']}.csv")
            network.data.to_csv(file)
            if not multiple:
                return save_file_response(file, "text/csv")
            else:
                multi.append(file)
        if bool_from_url_query(request.GET.get("lynx_failed")):
            file = os.path.join(
                PATH, f"LINEX_lynx_failed_{request.session['session_id']}.txt")
            with open(file, "w") as f:
                if network_entry.version == 1:
                    for lipid in network.unconverted:
                        f.write(f"{lipid}\n")
                else:
                    for lipid in network._incompatible_lipids:
                        f.write(f"{lipid}\n")
            if not multiple:
                return save_file_response(file, "text/plain")
            else:
                multi.append(file)
        # saving failed molecular species conversion
        if bool_from_url_query(request.GET.get("molspec_failed")):
            if network_entry.version == 2:
                file = os.path.join(
                    PATH,
                    f"LINEX_molspec_failed_{request.session['session_id']}.txt"
                )
                with open(file, "w") as f:
                    for lipid in network.gln.failed_molspec():
                        f.write(f"{lipid}\n")
                if not multiple:
                    return save_file_response(file, "text/plain")
                else:
                    multi.append(file)
        # saving statistics
        # TODO: make this look nice => put all comparisons/groups into one
        #       table (at least for node attributes)
        for attr in FORM2ATTR.keys():
            if request.GET.get(attr):
                file = write_attribute_data(
                    request.session["session_id"], network, attr)
                if file:
                    # NOTE: this always results in a zip because we
                    # are not doing look-ahead here
                    if isinstance(file, str):
                        multi.append(file)
                        multiple = True
                    else:
                        if multiple:
                            multi += file
                        else:
                            multiple = True
                            multi = file
        if bool_from_url_query(request.GET.get("node_met")):
            file = f"LINEX_node_metrics_{request.session['session_id']}.csv"
            save_node_metrics(network, file)
            if multiple:
                multi.append(file)
            else:
                return save_file_response(file, "text/plain")
        # saving network statistics
        if multiple:
            zipped = os.path.join(
                PATH, f"LINEX_{request.session['session_id']}_download.zip")
            zip_files(zipped, multi)
            return save_file_response(zipped, "application/zip")
        return render(
            request, 'lipid_network/download.html',
            context={'has_lynx': has_lynx,
                     'not_computed': not_computed,
                     'network_form': network_form,
                     'stats_form': stats_form,
                     'lynx_form': lynx_form,
                     'json_form': json_form}
        )
    elif request.method == "POST":
        json_data = load(request.FILES['json_file'])
        form_data = request.POST.dict()
        fs = form_data.get('font_size')
        return pdf_download(request, network_data=json_data, fs=fs)
    else:
        if ComputedNetwork.objects.filter(
                userid=request.session["session_id"]).count() == 1:
            not_computed = False
            network_entry = ComputedNetwork.objects.get(
                userid=request.session["session_id"])
            has_lynx = network_entry.converted
            version = network_entry.version
        else:
            version = None
        # if has_lynx:
        #     form.lynx_data.disabled = False
        return render(
            request, 'lipid_network/download.html',
            context={'has_lynx': has_lynx,
                     'not_computed': not_computed,
                     'network_form': network_form,
                     'stats_form': stats_form,
                     'lynx_form': lynx_form,
                     'timeout': TIMEOUT_DELTA_MIN,
                     'version': version,
                     'json_form': json_form}
        )


def json_to_pdf(request, fs=7):
    if request.method == "POST":
        json_form = JSONDownload(request.POST)
        json_data = load(open(json_form.file, 'r'))
        return pdf_download(request, network_data=json_data, fs=fs)


def enrichment(request):
    delete_db_by_time()
    delete_file_entries()
    check_user_id(request)
    # enrichment request
    if request.method == "POST":
        print('received POST request')
        message = ''
        if check_running_task(request.session["session_id"]):
            message = 'You already have a computation running on the server ' \
                      '(Network generation or enrichment). ' \
                      'Please wait for the task to finish.'
            return JsonResponse(
                {'error': True, 'message': message}, safe=False)
        # parameter checks
        if int(request.GET.get('min_size', 5)) > int(
            request.GET.get('max_size', 30)):
            message = 'Min size cannot be larger than max size'
            return JsonResponse(
                {'error': True, 'message': message}, safe=False)
        if int(request.GET.get('max_iter', 100)) > 1000:
            message = 'Max iter cannot be larger than 1000'
            return JsonResponse(
                {'error': True, 'message': message}, safe=False)
        # network checks
        network_entry = check_version2_network(request)
        if not isinstance(network_entry, ComputedNetwork):
            return network_entry
        if int(network_entry.version) != 2:
            message = 'Networks computed with version 1 cannot be used for ' \
                      'network enrichment.' \
                      'Please compute a version 2 network.'
            return JsonResponse(
                {'error': True, 'message': message}, safe=False)

        if NetworkEnrichment.objects.filter(
            userid=request.session["session_id"]).count() == 0:
            enrichment_entry = NetworkEnrichment()
            enrichment_entry.userid = request.session['session_id']
            enrichment_entry.timestamp = UserIds.objects.get(
                userid=request.session["session_id"]).timestamp
            enrichment_entry.status = 2
            enrichment_entry.started = True
            enrichment_entry.computed_enrichments = {}
        else:
            enrichment_entry = NetworkEnrichment.objects.get(
                userid=request.session["session_id"])
            # check if user is already running a task
            if enrichment_entry.status == 2:
                return JsonResponse(
                    {'error': True, message: 'Enrichment already running!'})
            # Same groups
            if enrichment_key(
                request.GET) in enrichment_entry.computed_enrichments.keys():
                # Same parameters
                function_call = build_function_call(parameters=request.GET)
                if enrichment_entry.computed_enrichments[
                    enrichment_key(request.GET)][
                    'function_call'] == function_call:
                    if not checkbox_on(request.GET, 'force_overwrite', False):
                        return JsonResponse({'error': False,
                                             message: 'Enrichment already available, '
                                                      'proceed to GET query.'})

        enrichment_entry.save()

        enrichment_dict = dict(request.GET)
        enrichment_dict['session_id'] = request.session['session_id']

        # Set running task for user to True
        set_running_task(request.session["session_id"], value=True)

        run_enrichment(enrichment_dict)

        enrichment_entry.status = 2
        enrichment_entry.started = True
        enrichment_entry.save()

        print('running enrichment')
        return JsonResponse(
            {'error': False, message: 'Enrichment computation started'},
            safe=False)

    # progress update fetching
    if request.method == "GET":
        # report progress
        if NetworkEnrichment.objects.filter(
            userid=request.session["session_id"]).count() != 0:
            enrichment_entry = NetworkEnrichment.objects.get(
                userid=request.session["session_id"])

            if enrichment_entry.status == 0:
                if enrichment_key(
                    request.GET) in enrichment_entry.computed_enrichments.keys():
                    return JsonResponse({"error": False,
                                         "enrichment":
                                             enrichment_entry.computed_enrichments[
                                                 enrichment_key(request.GET)
                                             ],
                                         "status": 0})
                else:
                    return JsonResponse(
                        {'error': True, 'message': 'No enrichment available. '
                                                   'Please POST a request.',
                         'status': 1})

            elif enrichment_entry.status == 1:
                return JsonResponse({'error': True,
                                     'message': 'Error in enrichment analysis.',
                                     'status': 1})
            else:
                return JsonResponse({"error": False,
                                     'status': 2,
                                     'message': 'Computation in Progress'})
        else:
            return JsonResponse(
                {'error': True, 'message': 'No enrichment available. '
                                           'Please POST a request.',
                 'status': 1})


def lion_download(request):
    if request.GET:
        if NetworkEnrichment.objects.filter(
            userid=request.session["session_id"]).count() != 1 \
            or ComputedNetwork.objects.filter(
            userid=request.session["session_id"]).count() != 1:
            # TODO: error message
            return HttpResponse()
        enrichment_entry = NetworkEnrichment.objects.get(
            userid=request.session["session_id"])
        network_entry = ComputedNetwork.objects.get(
            userid=request.session["session_id"])
        vis_network = pickle.load(
            open(os.path.join(PATH, network_entry.vis_network), "rb")
        )
        key = request.GET['comparison']
        if key not in enrichment_entry.computed_enrichments.keys():
            return HttpResponse()
        subnet = enrichment_entry.computed_enrichments[key]['subnetworks']
        # filtering out lipid nodes
        # TODO: for each component separately
        enriched_nodes = set()
        for vis_graph in subnet.values():
            for node in vis_graph['nodes']:
                if node['shape'] == 'dot':
                    enriched_nodes.add(node['label'])
        nodes = set()
        for node in json.loads(vis_network['networks']['native']['nodes']):
            nodes.add(node['label'])
        background_nodes = nodes.difference(enriched_nodes)
        base_file = os.path.join(PATH,
                                 f"LINEX_{request.session['session_id']}_LION_enrichment")
        target_file = f'{base_file}_TargetList.csv'
        with open(target_file, 'w') as file:
            for node in enriched_nodes:
                file.write(f'{node}\n')
        background_file = f'{base_file}_BackgroundList.csv'
        with open(background_file, 'w') as file:
            for node in background_nodes:
                file.write(f'{node}\n')
        zip_file = os.path.join(PATH,
                                f'{base_file}_{key.replace(",", "_")}.zip')
        zip_files(zip_file, [target_file, background_file])
        return save_file_response(zip_file, "application/force-download")


def tutorial(request):
    delete_db_by_time()
    delete_file_entries()
    check_user_id(request)
    return render(
        request, 'lipid_network/tutorial.html',
        context={"timeout": TIMEOUT_DELTA_MIN}
    )


def request_data_delete(request):
    delete_db_by_time()
    delete_file_entries()
    check_user_id(request)
    return render(
        request, 'lipid_network/request-data-delete.html',
        context={"timeout": TIMEOUT_DELTA_MIN}
    )


def data_deleted(request):
    if not request.session.get("session_id"):
        return render(request, 'lipid_network/data-deleted.html')
    delete_db_by_time()
    delete_file_entries()
    delete_my_data(request.session["session_id"])
    return render(request, 'lipid_network/data-deleted.html')


def example_data(request):
    return save_file_response(
        os.path.join(settings.STATIC_DIR, "example_data.zip"),
        "application/force-download"
    )


def user_contribution(request):
    delete_db_by_time()
    delete_file_entries()
    check_user_id(request)
    contribform = ContributionForm()
    contrib_success = ""

    email = request.GET.get("email", "")
    text = request.GET.get("text", "")

    if email != "" and text != "":
        curr_contrib = UserContribution()
        curr_contrib.userid = request.session["session_id"]
        curr_contrib.email = email
        curr_contrib.text = text
        curr_contrib.save()
        contrib_success = "Contribution successfully submitted."

    return render(request, 'lipid_network/user-contribution.html',
                  context={'contribform': contribform,
                           "contrib_success": contrib_success})
