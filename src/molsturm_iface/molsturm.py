import molsturm_iface as __iface
import numpy as np

def __np_to_coords(arr):
  ret = __iface.CoordVector()
  for c in arr:
    if (len(c) != 3):
      raise ValueError("All items of the coordinates array need have exactly 3 items.")

    new = __iface.Coord()
    new[0] = c[0]
    new[1] = c[1]
    new[2] = c[2]
    ret.push_back(new)
  return ret

def __np_to_atnum(li):
  ret = __iface.AtomicNumberVector()
  for n in li:
    ret.push_back(n)
  return ret

def __np_to_nlm(arr):
  if len(arr.shape) != 2 or arr.shape[1] != 3:
    raise ValueError("Nlm array needs to be of shape n times 3")

  ret = __iface.NlmBasis()
  for b in arr:
    new = __iface.Nlm()
    new[0] = b[0]
    new[1] = b[1]
    new[2] = b[2]
    ret.push_back(new)
  return ret

# Special input parameters for which the above conversion functions need
# to be used before assignment
__params_transform_maps = {
  "coords": __np_to_coords,
  "atom_numbers": __np_to_atnum,
  "nlm_basis": __np_to_nlm,
}

"""The list of keys understood by the hartree_fock function"""
hartree_fock_keys = [ k for k in dir(__iface.Parameters) if k[0] != "_" ]

def hartree_fock(**kwargs):
  # The list of valid keys is the list of keys
  # with the special ones (starting with __) removed.
  params_keys = [ k for k in dir(__iface.Parameters) if k[0] != "_" ]
  res_keys = [ k for k in dir(__iface.HfResults) if k[0] != "_" ]

  # Build params and run:
  params = __iface.Parameters()
  for key in kwargs:
    if not key in params_keys:
      raise ValueError("Keyword " + key + " is unknown to hartree_fock")
    elif key in __params_transform_maps:
      setattr(params,key,__params_transform_maps[key](kwargs[key]))
    else:
      setattr(params,key,kwargs[key])
  res = __iface.hartree_fock(params)

  # Build output dictionary:
  out_arrays = [ "coeff_fb", "fock_ff", "orbital_energies_f",
                 "repulsion_integrals_ffff"]
  shape_lookup = { "f": res.n_orbs_alpha + res.n_orbs_beta,
                   "b": res.n_bas }

  out = { k :getattr(res,k) for k in res_keys if not k in out_arrays }
  for k in out_arrays:
    # Build the shape to cast the numpy arrays into from the
    # suffixes (e.g. _ffff, _bf) and the shape lookup object
    # we created above
    target_shape = tuple( shape_lookup[c] for c in k[k.rfind("_")+1:] )
    out[k] = np.array(getattr(res,k)).reshape(target_shape)

  return out

def print_mo_occupation(hfres, indention=""):
  """
  Print the MO occupation number and MO energies of the
  calculation results obtained with the function hartree_fock
  above.
  """
  fstr = indention
  if hfres["restricted"]:
    fstr+="{aocc:1s}  {ena:20}  {bocc:1s}"
  else:
    fstr+="{aocc:1s}  {ena:20}  |  {enb:20}  {bocc:1s}"
  print(fstr.format(aocc="a",bocc="b",enb="",ena=""))

  for i in range(max(hfres["n_orbs_alpha"],hfres["n_orbs_beta"])):
    kw = {
      "aocc": "*" if i < hfres["n_alpha"] else " ",
      "bocc": "*" if i < hfres["n_beta"]  else " ",
      "ena": "",
      "enb": "",
    }
    if i < hfres["n_orbs_alpha"]:
      kw["ena"] = hfres["orbital_energies_f"][i]
    if i < hfres["n_orbs_beta"]:
      kw["enb"] = hfres["orbital_energies_f"][i+hfres["n_orbs_alpha"]]
    print(fstr.format(**kw))

def print_energies(hfres, indention=""):
  """
  Print the energies of a calculation obtained with the
  hartree_fock function above.
  """
  energies = [ k for k in hfres if k.startswith("energy_") ]
  maxlen = max([ len(k) for k in energies ])
  fstr=indention + "{key:"+str(maxlen)+"s} = " + "{val:20}"

  for k in energies:
    print(fstr.format(key=k, val=hfres[k]))

def build_pyadc_input(hfres):
  """
  Take the results dictionary from a hf calculation and build
  the input dictionary for a pyadc run out of it.
  """
  params = { k:hfres[k] for k in hfres
             if not k in [ "n_bas", "energy_total" ] }
  params["energy_scf"] = hfres["energy_total"]
  return params
