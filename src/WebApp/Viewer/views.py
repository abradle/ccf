from django.http import HttpResponse,HttpResponseRedirect
# Import the models
from IOhandle.models import Project,Protein
from django.shortcuts import render, get_object_or_404
from Viewer.models import ViewHandler
from djutils.decorators import async


def index(request):
    project_list = Project.objects.all()
    context = {'project_list': project_list}
    return render(request, 'Viewer/index.html', context)# Create your views here.


def loader(request):
    """Function to take requests for images and molecular files return them.
    Takes a request with GET attributes.
    Returns an HTPPResponse."""
    # The following are in all requests
    # The option in the function
    option = request.GET['choice']
    # Function choice
    function = request.GET['function']
    #Because sometimes it comes back like this (from Ubuntu) - the Saulo conundrum
    function = function.split(".pdb")[0]
    functions = ViewHandler()
    try:
        groupchoice = request.GET['groupchoice']
        HTML_out = functions.pdbroutine_list[function](option, groupchoice)
        return HttpResponse(HTML_out)
    except KeyError:
        pass
    # Options that are possible
    target_id = None
    mymap = None
    out_put = None
    extra = None
    try:
        mymap = request.GET['map']
    except KeyError:
        pass
    try:
        target_id = request.GET['target']
    except KeyError:
        pass
    try:
        out_put = request.GET['output']
    except KeyError:
        pass
    try:
        extra = request.GET['option']
    except KeyError:
        pass
    # Passing these options to the functions.
    HTML_out = functions.pdbroutine_list[function](option, maps=mymap, out_put=out_put, target_id=target_id, extra=extra)
    return HttpResponse(HTML_out)
