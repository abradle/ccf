from django.conf.urls import patterns, include, url
from django.conf import settings
from django.conf.urls.static import static
# Uncomment the next two lines to enable the admin:
#from django.contrib import admin
#admin.autodiscover()
import views
urlpatterns = patterns('',
    url(r'^$', views.index, name='index'),
    url(r'^Viewer/', include('Viewer.urls',namespace="Viewer")),
#    url(r'^SUPERSTAR/', include('SUPERSTAR.urls',namespace="SUPERSTAR")),
    #url(r'^GLOOP/', include('GLOOP.urls',namespace="GLOOP")),
    url(r'^OOMMPPAA/', include('OOMMPPAA.urls',namespace="OOMMPPAA")),
    #url(r'^LLOOMMPPAA/', include('LLOOMMPPAA.urls',namespace="LLOOMMPPAA")),
    url(r'^WONKA/', include('WONKA.urls',namespace="WONKA")),
#    url(r'^admin/', include(admin.site.urls)),
    url( r'^upload/', views.upload, name = 'jfu_upload' ),
    url(r'^run/$', views.run, name='run'),
)+static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)
