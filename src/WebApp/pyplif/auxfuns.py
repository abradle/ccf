from LLOOMMPPAA.functions import PlifBitInstance
from django.core.exceptions import ValidationError
def my_fun(dj_mol,new_bit,new_plif,bit):
  new_pbi = PlifBitInstance()
  new_pbi.mol_id = dj_mol
  new_pbi.plif_bit_id = new_bit
  new_pbi.on= bit
  try:
    new_pbi.validate_unique()
    new_pbi.save()
  except ValidationError:
    new_pbi = PlifBitInstance.objects.get(mol_id= dj_mol,bit_id=new_bit)
    new_pbi.on = bit
    new_pbi.save()
    new_plif.bit_id.add(new_pbi)