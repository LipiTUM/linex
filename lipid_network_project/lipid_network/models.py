from django.db import models
from django.core.validators import FileExtensionValidator
from picklefield.fields import PickledObjectField
from django.conf import settings
from django.core.files.storage import FileSystemStorage

# Max length for user contribution
contrib_max_length = 800


def update_filename(instance, filename, add=None, fe=None):
    if fe is None:
        fe = filename.split(".")[-1]
    if add is None:
        new_name = f"uploaded_file_{instance.userid}.{fe.lower()}"
    else:
        new_name = f"uploaded_{add}_file_{instance.userid}.{fe.lower()}"
    return new_name


def update_group_filename(instance, filename):
    return update_filename(instance, filename, add="group")


def update_class_filename(instance, filename):
    return update_filename(instance, filename,
                           add="class", fe="txt")


def update_fa_filename(instance, filename):
    return update_filename(instance, filename,
                           add="fa", fe="txt")


def update_confirmed_species_filename(instance, filename):
    return update_filename(instance, filename,
                           add="confirmed_species", fe="txt")


class UploadedData(models.Model):
    userid = models.TextField(default="")
    docfile = models.FileField(
        upload_to=update_filename,
        validators=[FileExtensionValidator(allowed_extensions=['csv'])],
        storage=FileSystemStorage(location=settings.USER_DATA)
    )
    timestamp = models.DateTimeField(auto_now_add=True)

    def __str__(self):
        return f"{self.docfile} - {str(self.timestamp)}"


class UploadedGroups(models.Model):
    userid = models.TextField(default="")
    group_docfile = models.FileField(
        upload_to=update_group_filename,
        validators=[FileExtensionValidator(allowed_extensions=['csv'])],
        storage=FileSystemStorage(location=settings.USER_DATA),
        blank=True
    )
    forbidden_characters = models.BooleanField(default=False)
    timestamp = models.DateTimeField(auto_now_add=True)

    def __str__(self):
        return f"{self.group_docfile} - {str(self.timestamp)}"


class UploadedClassSettings(models.Model):
    userid = models.TextField(default="")
    class_settings = models.FileField(
        upload_to=update_class_filename,
        validators=[FileExtensionValidator(allowed_extensions=['txt'])],
        storage=FileSystemStorage(location=settings.USER_DATA),
        blank=True
    )
    timestamp = models.DateTimeField(auto_now_add=True)

    def __str__(self):
        return f"{self.class_settings} - {str(self.timestamp)}"


class UploadedFaSettings(models.Model):
    userid = models.TextField(default="")
    fa_settings = models.FileField(
        upload_to=update_fa_filename,
        validators=[FileExtensionValidator(allowed_extensions=['txt'])],
        storage=FileSystemStorage(location=settings.USER_DATA),
        blank=True
    )
    timestamp = models.DateTimeField(auto_now_add=True)

    def __str__(self):
        return f"{self.fa_settings} - {str(self.timestamp)}"


class UploadedConfirmedSpecies(models.Model):
    userid = models.TextField(default="")
    confirmed_species_file = models.FileField(
        upload_to=update_confirmed_species_filename,
        validators=[FileExtensionValidator(allowed_extensions=['txt'])],
        storage=FileSystemStorage(location=settings.USER_DATA),
        blank=True
    )
    timestamp = models.DateTimeField(auto_now_add=True)

    def __str__(self):
        return f"{self.confirmed_species_file} - {str(self.timestamp)}"


class ComputedNetwork(models.Model):
    userid = models.TextField(default="")
    timestamp = models.DateTimeField(auto_now_add=True)
    message = models.TextField(default="")
    network = models.FilePathField(default="")
    vis_network = models.FilePathField(default="")
    node_colours = PickledObjectField(default=list)
    edge_colours = PickledObjectField(default=list)
    node_sizes = PickledObjectField(default=list)
    edge_sizes = PickledObjectField(default=list)
    directed = models.BooleanField(default=False)
    converted = models.BooleanField(default=False)
    version = models.IntegerField(default=0)
    reference_group = models.TextField(default="")
    is_log = models.BooleanField(default=False)

    # version 2 specific attributes
    native_mol_net = models.FilePathField(default="")
    bip_net = models.FilePathField(default="")
    bip_mol_net = models.FilePathField(default="")

    def __str__(self):
        return f"{self.userid} - {str(self.timestamp)} - {str(self.network)})"


class LipidomeSummary(models.Model):
    userid = models.TextField(default="")
    timestamp = models.DateTimeField(auto_now_add=True)
    data = PickledObjectField(default=dict)

    def __str__(self):
        return f"Lipidome Summary {self.userid} - {str(self.timestamp)})"


class ChainLengthAnalysis(models.Model):
    userid = models.TextField(default="")
    timestamp = models.DateTimeField(auto_now_add=True)
    reference_group = models.TextField(default="")
    plot = PickledObjectField(default=dict)
    data = PickledObjectField(default=dict)

    def __str__(self):
        return f"Chain-length analysis {self.userid} - {str(self.timestamp)})"


class SubstructureAnalysis(models.Model):
    userid = models.TextField(default="")
    timestamp = models.DateTimeField(auto_now_add=True)
    data = PickledObjectField(default=dict)
    plot = PickledObjectField(default=dict)

    def __str__(self):
        return f"Substructure analysis {self.userid} - {str(self.timestamp)})"


class UserIds(models.Model):
    usercount = models.IntegerField(default=1)
    userid = models.TextField(default="")
    timestamp = models.DateTimeField(auto_now_add=True)
    version = models.IntegerField(default=0)
    running_task = models.BooleanField(default=False)

    def __str__(self):
        return f"{self.userid} - {str(self.timestamp)}"


class ComputationProgress(models.Model):
    userid = models.TextField(default="")
    timestamp = models.DateTimeField(auto_now_add=True)

    version = models.IntegerField(default=1)

    reading = models.TextField(default="Waiting")
    lynx = models.TextField(default="Waiting")
    unconverted = PickledObjectField(default=list)
    failed_molspec = PickledObjectField(default=list)
    class_init = models.TextField(default="Waiting")
    compute_network = models.TextField(default="Waiting")
    compute_statistics = models.TextField(default="Waiting")
    compute_views = models.TextField(default="Waiting")
    compute_summaries = models.TextField(default="Waiting")
    groups = models.BooleanField(default=False)

    done = models.BooleanField(default=False)
    error = models.BooleanField(default=False)
    message = models.TextField(default="")
    warning = PickledObjectField(default=list)

    def __str__(self):
        return f"{self.userid} - {str(self.timestamp)}"


class NetworkEnrichment(models.Model):
    userid = models.TextField(default="")
    timestamp = models.DateTimeField(auto_now_add=True)

    started = models.BooleanField(default=False)
    # 0 => done
    # 1 => error
    # 2 => running
    status = models.IntegerField(default=0)
    # computed_enrichments = {'groups': {'subnetworks', 'nodes', 'final_scores', 'scores',
    # 'components', 'function_call'}}
    computed_enrichments = PickledObjectField(default=dict)

    def __str__(self):
        return f"{self.userid} - {str(self.timestamp)}"


class UserContribution(models.Model):
    userid = models.TextField()
    timestamp = models.DateTimeField(auto_now_add=True)
    email = models.EmailField()
    text = models.TextField(max_length=contrib_max_length)

    def __str__(self):
        return f"{self.userid} - {str(self.timestamp)} | \nEmail:{self.email} \nText: {self.text}"
