from django.urls import path
from django.conf.urls.static import static
from django.conf import settings
from . import views

#app_name = 'lipid_network'

urlpatterns = [
    path('', views.index, name='index'),
    path('about/', views.about, name='about'),
    path('about', views.about, name='about'),
    path('analysis', views.analysis, name='analysis'),
    path('download', views.download, name='download'),
    path('upload', views.upload, name='upload'),
    path('upload-pending', views.upload_pending, name='upload-pending'),
    path('tutorial', views.tutorial, name='tutorial'),
    path('request-data-delete', views.request_data_delete, name='request-data-delete'),
    path('data-deleted', views.data_deleted, name="data-deleted"),
    path('user-contribution', views.user_contribution, name="user-contribution"),
    path('change-network-type', views.change_network_type, name='change-network-type'),
    path('enrichment', views.enrichment, name='enrichment'),
    path('pdf-download', views.pdf_download, name='pdf-download'),
    path('substructure-download', views.substructure_download, name='substructure-download'),
    path('chainlength-download', views.chainlength_download, name='chainlength-download'),
    path('summary-download', views.download_summary_data, name='summary-download'),
    path('lion-download', views.lion_download, name='lion-download'),
    path('settings.STATIC_DIR/example_data.zip', views.example_data)
] + static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)
