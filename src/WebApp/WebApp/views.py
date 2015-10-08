from django.http import HttpResponse,HttpResponseRedirect
from django.shortcuts import render, get_object_or_404
# Import the models
import os,sys
from django.conf import settings
from django.core.urlresolvers import reverse
from django.views.decorators.http import require_POST
from jfu.http import upload_receive, UploadResponse, JFUResponse
import tempfile

# Global variable to hold list of files for fileupload
arg_list = {}


def index(request):
    context = {"projects": None}
    return render(request, 'index.html', context)

@require_POST
def upload(request):
    """JFU function to handle uplaods
    The assumption here is that jQuery File Upload
    has been configured to send files one at a time.
    If multiple files can be uploaded simulatenously,
    'file' may be a list of files.
    """
    my_file = upload_receive(request)
    #instance = YourUploadModel( file_field = file )
    basename = my_file.name
    tfile = tempfile.NamedTemporaryFile("w", delete=False)
    tfile.write(my_file.read())
    arg_list[basename] = tfile.name
    #basename = os.path.basename( instance.file_field.file.name )
    file_dict = {
        'name': basename,
        'size': my_file.size,
        #instance.file_field.file.size,
        # The assumption is that file_field is a FileField that saves to
        # the 'media' directory.
        'url': settings.MEDIA_URL + basename,
        'thumbnail_url': settings.MEDIA_URL + basename,
        'delete_type': 'POST',
    }
    return UploadResponse(request, file_dict)


def run(request):
    """Function to handle running for the windows application.
    Takes a request. Handles the different lements of the requests
    Returns an HTTPResponse"""
    targ = request.GET.get('targ')
    mols = request.GET.get('mols')
    acts = request.GET.get('acts')
    prot = request.GET.get('prot')
    try:
        if targ is None or mols is None or acts is None:
            print "ERROR"
        # Protein addition is not essential
        elif prot is None:
            my_s = "start cmd /c " + os.path.join(os.path.split(sys.argv[0])[0], "winwrap.exe")
            my_s = my_s + " --targ " + targ + " --mols " + arg_list[mols] + " --acts "
            my_s = my_s + arg_list[acts] + "  --mmp True > out.log"
            os.system(my_s)
        else:
            my_s = "start cmd /c " + os.path.join(os.path.split(sys.argv[0])[0], "winwrap.exe")
            my_s = my_s + " --targ " + targ + " --mols " + arg_list[mols] + " --acts "
            my_s = my_s + arg_list[acts] + "  --mmp True  --prot " + arg_list[prot] + " > out.log"
            os.system(my_s)
    except KeyError:
        return HttpResponse("ERROR")
    return HttpResponse("LOADING")
