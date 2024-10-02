from .shared_settings import *

# Database
# LOCAL DEPLOYMENT ----
DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.sqlite3',
        'NAME': os.path.join(BASE_DIR / 'db.sqlite3'),
    }
}

DEBUG = True
ALLOWED_HOSTS = ['0.0.0.0', '127.0.0.1']

# Static files (CSS, JavaScript, Images)
# https://docs.djangoproject.com/en/3.1/howto/static-files/
STATICFILES_DIRS = [STATIC_DIR, ]
STATIC_URL = '/static/'

MEDIA_ROOT = MEDIA_DIR
MEDIA_URL = '/media/'
