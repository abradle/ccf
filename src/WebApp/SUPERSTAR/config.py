from django.db import models


def find_allowed_ph4(item_dict):
    """Function to find allowed ph4"""
    # Create the out list to return
    out_list = []
    # Loop through the keys - except unused
    for item in item_dict:
        if item == "unused":
            continue
        out_list.extend(item_dict[item])
    return list(set(out_list))


class SuperStarConfig():
    """Model to store the required config information for SUPERSTAR"""
    def __init__(self):
        # First the path
        self.base_path = "W:\\Informatics\\Pharmacophore\\anthony\\DPhil\\CODE\\CHOC\\src\\WebApp\\SUPERSTAR\\PROTEINS"
        self.ins_path = "W:\\Informatics\\Pharmacophore\\anthony\\DPhil\\CODE\\CHOC\\src\\WebApp\\SUPERSTAR\\INS_FILES"
        self.out_path = "W:\\Informatics\\Pharmacophore\\anthony\\DPhil\\CODE\\CHOC\\src\\WebApp\\SUPERSTAR\\OUT_FILES"
        self.superstar_path = r"C:\Program Files (x86)\CCDC\SuperStar 2.1.2\bin\superstar_app.exe"
        self.probe_dict = {'ali': 'ALIPHATIC CH CARBON', 'aro':'AROMATIC CH CARBON', 'uncharged':'UNCHARGED NH NITROGEN', 'carbonyl_O':'CARBONYL OXYGEN', 'carboxylate': 'CARBOXYLATE OXYGEN', 'charged':'CHARGED NH NITROGEN', 'water':'WATER OXYGEN'}
        self.type_dict = {'unused': ['SmallHalogen','Cyano','Imidazole','Guanidine'],
 'water': [],
 'ali': ['RH6_6','iPropyl','Methyl','RH3_3', 'ThreeWayAttach'],
 'uncharged': ['SingleAtomDonor'],
 'carbonyl_O': ['SingleAtomAcceptor','Carbonyl'],
 'aro': ['Arom5','Arom6'],
 'carboxylate': ['AcidicGroup'],
 'charged': ['BasicGroup','PosN']
 }
        self.allowed_ph4 = find_allowed_ph4(self.type_dict)