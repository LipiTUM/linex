from django import forms
from .models import contrib_max_length

DISABLE_PARTIAL_CORRELATIONS = False


class DataForm(forms.Form):
    lipids_as_columns = forms.BooleanField(label='Lipids are in columns',
                                           initial=True, required=False)
    docfile = forms.FileField(label='Select a lipid data file')
    group_docfile = forms.FileField(
        label='Select a sample group file (optional)',
        required=False
    )
    class_file = forms.FileField(
        label='Select Lipid Class Settings',
        required=False
    )
    fa_file = forms.FileField(
        label='Select Fatty Acid Settings (optional)',
        required=False
    )
    confirmed_species_file = forms.FileField(
        label='Additional confirmed molecular species (optional)',
        required=False
    )
    version_number = forms.ChoiceField(
        label='LINEX Version',
        required=True, choices=((1, 1), (2, 2))
    )


class DataProcessOptions(forms.Form):
    log_data = forms.BooleanField(label="Data is log-transformed",
                                  required=False)
    log_ratios = forms.BooleanField(label="Fold-changes as log-ratios",
                                    required=False, initial=True)
    lynx_format = forms.BooleanField(label="Convert to LipidLynxX Nomenclature",
                                     required=False,
                                     help_text="<b>NOTE:</b> converting names may take up to several minutes")
    lipid_level = forms.ChoiceField(
        label="Lipid Resolution",
        choices=(
            ("molecular", "Molecular Species"),
            ("sum", "Sum Species"),
            ("sn", "sn-specific")
        ),
        help_text="Setting the <b>highest</b> resolution to consider.<br>Please only change if all lipids are sum "
                  "species or some are confident sn-specific identifications."
    )


class Version2Specific(forms.Form):
    # https://reactome.org/content/schema/objects/Species
    organism = forms.CharField(
        label="Reactome Organism (three letter code)", max_length=3, min_length=3, initial='HSA',
        help_text='see available organisms <a href="https://reactome.org/content/schema/objects/Species">here</a>'
    )
    database = forms.ChoiceField(
        label="Database(s) to include",
        choices=(
            ('both', 'Rhea and Reatome'),
            ('Rhea', 'Rhea'),
            ('Reactome', 'Reactome')
        )
    )
    ether_conversions = forms.BooleanField(
        label="Allow Ether Conversions", initial=True,
        required=False)


class StatisticsOptions(forms.Form):
    corrs = forms.BooleanField(label="Correlations/Correlation Changes",
                               required=False, initial=True)
    # TODO: enable once errors are caught
    pcorrs = forms.BooleanField(label="Partial-Correlations/-Correlation Changes",
                                required=False, initial=False,
                                disabled=DISABLE_PARTIAL_CORRELATIONS)
    corr_sig = forms.FloatField(
        label="Significance threshold for (partial) correlations",
        min_value=0, max_value=1, initial=.05
    )

    fcs = forms.BooleanField(label="Fold Changes",
                             required=False, initial=True)
    ref_group = forms.CharField(
        label="Reference Group",
        required=False, empty_value=None,
        help_text="Reference Group to compute fold-changes against"
    )
    pvals = forms.BooleanField(label="p-values (binary statistical test)",
                               required=False, initial=True)
    ptest = forms.ChoiceField(
        label="Statistical Test",
        choices=(
            ("ttest", "t-test"),
            ("mannwhitneyu", "Mann-Whithney U rank test"),
            ("ranksums", "Wilcoxon rank-sum test"),
            ("wilcoxon", "Wilcoxon signed rank test")
        )
    )


class NetworkOptions(forms.Form):
    directed = forms.BooleanField(
        label="Directed Graph",
        required=False,
        help_text="<b>NOTE:</b> directed edges increase memory and processor uptake"
    )


class NetworkDownload(forms.Form):
    network_html = forms.BooleanField(label="Network as standalone .html",
                                      required=False)
    network_graphml = forms.BooleanField(
        label="Network as .graphml (currently disabled)",
        required=False,
        help_text='<b>NOTE:</b> To have all node and edge attributes included '
                  'you need to visit the <a '
                  'href="analysis">Analysis</a> site first. ',
        disabled=True
    )
    colour_legend = forms.BooleanField(label="Colour Legend as .graphml",
                                       required=False, disabled=False)
    size_legend = forms.BooleanField(label="Size Legend as .graphml",
                                     required=False, disabled=False)
    connections_plot = forms.BooleanField(
        label="Class Connection Network as png (Only Version 1)",
        required=False, disabled=False)
    connections_graph = forms.BooleanField(
        label="Class Connection Network as graphml file (Only Version 1)",
        required=False, disabled=False)


class StatisticsDownload(forms.Form):
    corr_data = forms.BooleanField(label="Lipid Correlation Matrix/Matrices",
                                   required=False, disabled=False)
    pcorr_data = forms.BooleanField(label="Lipid Partial-Correlation Matrix/Matrices",
                                    required=False, disabled=DISABLE_PARTIAL_CORRELATIONS)
    fc_data = forms.BooleanField(label="Lipid (log) Fold-Changes",
                                 required=False, disabled=False)
    pval_data = forms.BooleanField(label="Lipid FDR-values",
                                   required=False, disabled=False)
    node_met = forms.BooleanField(label="Lipid Node Graph Metrics",
                                  required=False, disabled=False)


class LynxDownload(forms.Form):
    lynx_data = forms.BooleanField(label="LipidLynxX Converted Data as .csv",
                                   required=False, disabled=False)
    lynx_failed = forms.BooleanField(label="Unconverted Lipid Species",
                                     required=False, disabled=False)

    molspec_failed = forms.BooleanField(
        label="Failed Molecular Species Inference",
        required=False, disabled=False)


class JSONDownload(forms.Form):
    json_file = forms.FileField(label="Select the .json data file")
    font_size = forms.IntegerField(label="Select the font size", initial=7)


class ContributionForm(forms.Form):
    email = forms.EmailField(label="Email", required=True)

    help_text = f"(Max {contrib_max_length} characters, for longer texts " \
                "send us an email directly)"
    text = forms.CharField(
        label="Contribution", widget=forms.Textarea,
        max_length=contrib_max_length,
        help_text=help_text)
