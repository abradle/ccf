from django.http import HttpResponse,HttpResponseRedirect
from django.shortcuts import render, get_object_or_404

import WebApp.views as WV

from MMPMaker.functions import find_pharma_changes
from MMPMaker.models import MMPDiffMap
from IOhandle.models import Protein,Target,Project,Molecule


def index(request):
    """View for the OOMMPPAA main page"""
    targets = Target.objects.exclude(title="DUMMY")
    # No projects at the moment as this hasn't been set up yet
    context = {"targets": targets}
    return render(request, 'OOMMPPAA/index.html', context)


def fileupload(request):
    """View for the OOMMPPAA file upload"""
    WV.arg_list = {}
    targets = Target.objects.exclude(title="DUMMY")
    # No projects at the moment as this hasn't been set up yet
    context = {"targets": targets}
    return render(request, 'OOMMPPAA/fileupload.html', context)


def success(request):
    """View for the OOMMPPAA success page"""
    targets = Target.objects.exclude(title="DUMMY")
    # No projects at the moment as this hasn't been set up yet
    context = {"targets": targets}
    return render(request, 'OOMMPPAA/success.html', context)


def demoviewer(request, target_id):
    """View for the OOMMPPAA demo page"""
    target = get_object_or_404(Target, pk=target_id)
    # Add in a filter for proteins with actual molecules
    targets = Target.objects.exclude(title="DUMMY")
    maps = MMPDiffMap.objects.filter(target_id=target_id).filter(activity_type=str("['IC50', 'Ki', 'TM']")).order_by("pk")
    context = {"target": target, "maps": maps, "targets": targets}
    # The class defining the list of ways of
    return render(request, 'OOMMPPAA/demoviewer.html', context)


def viewer(request, target_id):
    """View for the OOMMPPAA main viewer page"""
    target = get_object_or_404(Target, pk=target_id)
    target.mol = Molecule.objects.filter(prot_id__target_id=target_id).exclude(prot_id__code__startswith=target.title)[0].prot_id.code
    targets = Target.objects.exclude(title="DUMMY")
    mols = Molecule.objects.filter(prot_id__target_id=target_id).exclude(prot_id__code__startswith=target.title)
    # Add in a filter for proteins with actual molecules
    maps = MMPDiffMap.objects.filter(target_id=target_id).filter(activity_type=str("['IC50', 'Ki', 'TM']"))
    context = {"mols":  mols, "target": target, "maps": maps, "targets": targets}
    # The class defining the list of ways of
    return render(request, 'OOMMPPAA/viewer.html', context)


def plugin(request):
    """View for the PLUGIN failure page"""
    browser = request.GET["browser"]
    platform = request.GET["platform"]
    return render(request, 'OOMMPPAA/plugin.html', {"browser": browser, "platform": platform})