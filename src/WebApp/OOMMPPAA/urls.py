
from django.conf.urls import patterns, url

from OOMMPPAA import views

urlpatterns = patterns('',
    url(r'^(?P<target_id>\d+)/demoviewer/$', views.demoviewer, name='demoviewer'),
    url(r'^(?P<target_id>\d+)/viewer/$', views.viewer, name='viewer'),
    url(r'^$', views.index, name='index'),
    url(r'^fileupload/$', views.fileupload, name='fileupload'),
    url(r'^success/$', views.success, name='success'),
    url(r'^plugin/$', views.plugin, name='plugin'),
)
