from WONKA.models import KeyCluster
from IOhandle.models import Molecule
from GLOOP.models import Transformation
import json
from rdkit import Chem
from django.shortcuts import render, redirect
from django.http import HttpResponse, HttpResponseRedirect
from Viewer.models import eucl_dist
from MMPMaker.models import MMPFrag
# "MOL_ONE": x.mol_one.pk, "MOL_TWO": x.mol_two.pk}


def get_cluster_info(request, cluster_id):
    """View to get the information for a given cluster"""
    dist = 2.0
    kc = KeyCluster.objects.get(pk=int(cluster_id))
    my_res = [{"frag_from": x.before_smiles, "smiles": x.after_smiles, "act_change": x.activity_change, "mol_from": x.mol_one.cmpd_id.pk ,"mol_to": x.mol_two.cmpd_id.pk} for x in Transformation.objects.filter(before_smiles=kc.function) if eucl_dist(x,kc) < dist]
    return HttpResponse(json.dumps(my_res), content_type="application/json")


def ShowClusterSummary(request, cluster_id):
    """View to show a summary of info for a cluster"""
    context = {"my_pk": cluster_id}
    return render(request, 'GLOOP/ShowClusterSummary.html', context)


def ShowMolClusters(request, mol_id):
    """"View to find all the clusters for a mol"""
    # First get all the clusters this mols is relate to
    my_clusts = KeyCluster.objects.filter(mmp_frag_id__mol_id__pk=mol_id)
    my_mol = Molecule.objects.get(pk=int(mol_id))
    target = my_mol.prot_id.target_id
    context = {"target": target,"my_pk": my_clusts[0].pk, "my_clusts": my_clusts, "mol_id": mol_id, "prot_code": my_mol.prot_id.code}
    return render(request, 'GLOOP/ShowMolClusters.html', context)


def ShowMolClustersPDB(request):
    """View to find a molecule from a PDB code"""
    pdb_code = request.GET["PDB_CODE"]
    mol_id = Molecule.objects.filter(prot_id__code=pdb_code)[0]
    return HttpResponse(mol_id.pk)


def get_frag_sdf(request, frag_id):
    """View to get the fragment information (in 3D) for a mol"""
    # Get the fragment
    mol_id = request.GET["MOL_ID"]
    my_frag = MMPFrag.objects.filter(keycluster__pk=int(frag_id), mol_id=int(mol_id))[0]
    in_mol = Chem.MolFromMolBlock(str(my_frag.sdf_info))
    # Remove the zero atoms
    out_mol = Chem.DeleteSubstructs(in_mol, Chem.MolFromSmarts('[#0]'))
    return HttpResponse(Chem.MolToMolBlock(out_mol))