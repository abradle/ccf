from django.conf.urls import patterns, url
from LLOOMMPPAA import views

urlpatterns = patterns('',
    url(r'^(?P<target_id>\d+)/div_viewer/$', views.div_viewer, name='div_viewer'),
    url(r'^(?P<target_id>\d+)/int_picker/$', views.int_picker, name='int_picker'),
    url(r'^(?P<target_id>\d+)/res_finder/$', views.res_finder, name='res_finder'),
    url(r'^(?P<target_id>\d+)/FollowUpMaker/$', views.FollowUpMaker, name='FollowUpMaker'),
    url(r'^(?P<frag_id>\d+)/get_sdf_info/$', views.get_sdf_info, name='get_sdf_info'),
    url(r'^(?P<frag_id>\d+)/(?P<choice>\d+)/FollowUpFrag/$', views.FollowUpFrag, name='FollowUpFrag'),
    url(r'^$', views.index, name='index'),
    url(r'^(?P<job_id>\d+)/CheckProgress/$', views.CheckProgress, name='CheckProgress'),
    url(r'^SubmitJob/$', views.SubmitJob, name='SubmitJob/'),
    url(r'^get_sdf_from_list/$', views.get_sdf_from_list, name='get_sdf_from_list/'),
    url(r'^(?P<job_id>\d+)/FollowUpReview/$', views.FollowUpReview, name='FollowUpReview'),
    url(r'^(?P<ra_id>\d+)/MakeDendrogram/$', views.MakeDendrogram, name='MakeDendrogram'),
    url(r'^(?P<ra_id>\d+)/GetMols/$', views.GetMols, name='GetMols'),
    url(r'^(?P<ra_id>\d+)/CompareMaps/$', views.CompareMaps, name='CompareMaps'),
    url(r'^(?P<ra_id>\d+)/CompareHists/$', views.CompareHists, name='CompareHists'),
    url(r'^(?P<job_id>\d+)/DownloadLibs/$', views.DownloadLibs, name='DownloadLibs'),
    url(r'^(?P<ra_id>\d+)/get_mol_sdf/$', views.get_mol_sdf, name='get_mol_sdf'),
    url(r'^(?P<ra_id>\d+)/get_mol_ids/$', views.get_mol_ids, name='get_mol_ids'),
    url(r'^(?P<target_id>\d+)/TargetSummary/$', views.TargetSummary, name='TargetSummary'),

)
