"""lipid_network_project URL Configuration

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/3.1/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  path('', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  path('', Home.as_view(), name='home')
Including another URLconf
    1. Import the include() function: from django.urls import include, path
    2. Add a URL to urlpatterns:  path('blog/', include('blog.urls'))
"""
from django.contrib import admin
from django.urls import path
from django.urls import include
import os

root = os.environ.get('ROOT_DOMAIN', '').replace("/", "")
if root == '':
    urlpatterns = [
        path('', include('lipid_network.urls')),
        path('admin/', admin.site.urls)
    ]
else:
    urlpatterns = [
        path('', include('lipid_network.urls')),
        path(f'{root}/', include('lipid_network.urls')),
        path(f'{root}/admin/', admin.site.urls),
    ]

