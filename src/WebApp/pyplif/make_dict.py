def trans_res(res):
  """Function to take a RES of form blah./^V64 and output VAL64"""
  d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
  rev_d = dict((v, k) for (k, v) in d.iteritems())
  # Take the res
  out_res = rev_d[res.split("^")[1][0]]
  out_res += res.split("^")[1][1:]
  return out_res

# Now open the file
res_list = [x.rstrip() for x in open("myFirstFile.txt").readlines()]
# Make the output dictionary
out_d = {}
# write it
for res in res_list:
  if res[2:6] in out_d:
    out_d[res[2:6]].append(trans_res(res))
  else:
    out_d[res[2:6]] = [trans_res(res)]

# Now make a list of all in the same sequence space
# Make a union of these lists and a propensity dict
prot_list = []
tot_list = []
restrictive_list = []
tot_dict = {}
for prot in out_d:
  if out_d[prot][0] == "ILE10":
    prot_list.append(prot)
    tot_list.extend(out_d[prot])
    for res in  out_d[prot]:
      if res in tot_dict:
        tot_dict[res] +=1
      else:
        tot_dict[res] =1
list(set(tot_list))
import operator
sorted_x = sorted(tot_dict.iteritems(), key=operator.itemgetter(1),reverse=True)
core_res = [x[0] for x in sorted_x[:16]]
# Now make a restrictive list
for prot in out_d:
  if len([x for x in out_d[prot] if x in core_res]) == len(core_res):
    restrictive_list.append(prot)
restrictive_list[:30]