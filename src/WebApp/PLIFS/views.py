from django.http import HttpResponse,HttpResponseRedirect
from django.shortcuts import render, get_object_or_404

import WebApp.views as WV

from MMPMaker.functions import find_pharma_changes
from MMPMaker.models import MMPDiffMap
from IOhandle.models import Protein,Target,Project,Molecule


def index(request):
    targets = Target.objects.exclude(title="DUMMY")
    # No projects at the moment as this hasn't been set up yet
    context = {"targets": targets}
    return render(request, 'PLIFS/index.html', context)
