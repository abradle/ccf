from django.conf.urls import patterns, url
from GLOOP import views

urlpatterns = patterns('',
    url(r'^(?P<cluster_id>\d+)/get_cluster_info/$', views.get_cluster_info, name='get_cluster_info'),
    url(r'^(?P<cluster_id>\d+)/ShowClusterSummary/$', views.ShowClusterSummary, name='ShowClusterSummary'),
    url(r'^(?P<frag_id>\d+)/get_frag_sdf/$', views.get_frag_sdf, name='get_frag_sdf'),
    url(r'^(?P<mol_id>\d+)/ShowMolClusters/$', views.ShowMolClusters, name='ShowMolClusters'),
    url(r'^ShowMolClustersPDB/$', views.ShowMolClustersPDB, name='ShowMolClustersPDB'),
    )