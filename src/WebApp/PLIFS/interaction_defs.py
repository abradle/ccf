from PLIFS.models import InteractionScheme, Interaction
from IOhandle.models import Target
from PLIFS.functions import find_interactions

schemes = {"default":
            {
             "Arom5_Arom5":
             {"dist": 4.5,
              "angle_1": None,
              "angle_2": None},
             "Arom5_Arom6":
             {"dist": 4.5,
              "angle_1": None,
              "angle_2": None},
             "Arom6_Arom6":
             {"dist": 4.5,
              "angle_1": None,
              "angle_2": None},
            "AcidicGroup_BasicGroup":
             {"dist": 4.5,
              "angle_1": None,
              "angle_2": None},
            "Hydrophobe_Hydrophobe":
             {"dist": 4.0,
              "angle_1": None,
              "angle_2": None},
            "SingleAtomAcceptor_SingleAtomDonor":
             {"dist": 3.5,
              "angle_1": "120:360",
              "angle_2": None}
             },
           }


def init_schemes(mol_ids=None, refresh=False, react_id=None):
    """Function to initialise schemes.
    If they change - it then  needs to update all interactions for this scheme"""
    # Now go through the schemes
    targets = Target.objects.exclude(title="DUMMY")
    if react_id:
        targets = targets.filter(pk=react_id.mol_id.prot_id.target_id.pk)
    tot = len(targets)
    if react_id:
        react_id.proc_stage = "MAKE INT SCHEMES"
        react_id.save()
        react_id.tot_sch = len(schemes.keys())
    counter = 0
    for scheme in schemes:
        print "SCHEME:",scheme
        counter += 1
        if react_id:
            react_id.sch = counter
            react_id.save()
        i_schemes = InteractionScheme.objects.filter(scheme_name=scheme)
        # If it's a new scheme
        if len(i_schemes) == 0:
            new_scheme = InteractionScheme()
            new_scheme.scheme_name = scheme
            new_scheme.scheme_dict = str(schemes[scheme])
            new_scheme.save()
            # Now make it for all targets
            if mol_ids:
                # Find the target in the schemes
                target_ids = Target.objects.filter(pk__in=list(set(mol_ids.values_list("prot_id__target_id__pk", flat=True))))
                t_count = len(target_ids)
                for t_count, target in enumerate(target_ids):
                    find_interactions(target.pk, t_count, tot, scheme, mol_ids=mol_ids.filter(prot_id__target_id=target), react_id=react_id)
            else:
                for t_count, target in enumerate(targets):
                    find_interactions(target.pk, t_count, tot, scheme, react_id=react_id)
        # If it's not
        elif len(i_schemes) == 1:
            i_scheme = i_schemes[0]
            # If it's been updated
            if i_scheme.scheme_dict != str(schemes[scheme]) or refresh == True:
                i_scheme.scheme_dict = str(schemes[scheme])
                i_scheme.save()
                # Now delete all points left - only for the ones we're goint to remake
                if mol_ids:
                    i_objs = Interaction.objects.filter(scheme_id__scheme_name=scheme, mol_id__in=mol_ids)
                    [x.delete() for x in i_objs]
                else:
                    i_objs = Interaction.objects.filter(scheme_id__scheme_name=scheme)
                    [x.delete() for x in i_objs]
                # Now make it for all targets
                if mol_ids:
                    # Find the target in the schemes
                    target_ids = Target.objects.filter(pk__in=list(set(mol_ids.values_list("prot_id__target_id__pk", flat=True))))
                    t_count = len(target_ids)
                    for t_count, target in enumerate(target_ids):
                        find_interactions(target.pk, t_count, tot, scheme, mol_ids=mol_ids.filter(prot_id__target_id=target), react_id=react_id)
                else:
                    for t_count, target in enumerate(targets):
                        find_interactions(target.pk, t_count, tot, scheme, react_id=react_id)
            else:
                continue
        else:
            print "ERROR INVALID NUMBER OF SCHEMES - UNIQUE CONSTRAINT BREACHED"
