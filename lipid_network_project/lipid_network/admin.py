from django.contrib import admin

from .models import (
    ComputedNetwork, UploadedData, UserIds,
    ComputationProgress, UploadedGroups,
    UploadedFaSettings, UploadedClassSettings,
    UploadedConfirmedSpecies, NetworkEnrichment,
    UserContribution
)

admin.site.register(UploadedData)
admin.site.register(UploadedGroups)
admin.site.register(UploadedFaSettings)
admin.site.register(UploadedClassSettings)
admin.site.register(ComputedNetwork)
admin.site.register(UserIds)
admin.site.register(ComputationProgress)
admin.site.register(UploadedConfirmedSpecies)
admin.site.register(NetworkEnrichment)
admin.site.register(UserContribution)
