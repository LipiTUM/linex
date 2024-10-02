from .shared_settings import *
import os


# Set in .env.sample in root
db_name = os.environ.get("POSTGRES_NAME", "postgres")
db_user = os.environ.get("POSTGRES_USER", "postgres")
db_password = os.environ.get("POSTGRES_PASSWORD", "postgres")

# Database
DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.postgresql',
        'NAME': db_name,
        'USER': db_user,
        'PASSWORD': db_password,
        'HOST': 'db',
        'PORT': 5432,
    }
}

DEBUG = False
ALLOWED_HOSTS = ['0.0.0.0', '127.0.0.1', 'exbio.wzw.tum.de']

# Static files (CSS, JavaScript, Images)
# https://docs.djangoproject.com/en/3.1/howto/static-files/
STATICFILES_DIRS = [STATIC_DIR, ]

STATIC_ROOT = os.path.join(BASE_DIR, 'run_static')
MEDIA_ROOT = MEDIA_DIR

MEDIA_URL = f'/linex/media/'
STATIC_URL = f'/linex/static/'


if os.environ.get("HTTPS", False):
    SESSION_COOKIE_PATH = f'{ROOT_DOMAIN}/'
    CSRF_COOKIE_PATH = f'{ROOT_DOMAIN}/'
CSRF_COOKIE_SESSION = True
SESSION_COOKIE_SECURE = True

# TODO
# ADMINS = (
#     ("", ""),
# )
# MANAGERS = ADMINS

# TODO
# server = {
#
# }

# TODO: logging

