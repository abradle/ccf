from django.conf.urls import patterns, url
from SUPERSTAR import views

urlpatterns = patterns('',
    url(r'^(?P<target_id>\d+)/ShowPlifs/$', views.ShowPlifs, name='ShowPlifs'),
    url(r'^GetPlifs/$', views.GetPlifs, name='GetPlifs'),
    url(r'^get_sim_map/$', views.get_sim_map, name='get_sim_map'),
    url(r'^loadmol/$', views.loadmol, name='loadmol'),
    url(r'^get_plif_json/$', views.get_plif_json, name='get_plif_json'),
#    url(r'^$', views.index, name='index'),
)