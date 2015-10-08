from django.conf.urls import patterns, url

from Viewer import views


urlpatterns = patterns('',
    # ex: /polls/
    url(r'^$', views.index, name='index'),
    url(r'^loader/$', views.loader, name='loader'),
)
