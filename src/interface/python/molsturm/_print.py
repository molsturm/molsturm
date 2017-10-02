## vi: tabstop=4 shiftwidth=4 softtabstop=4 expandtab
## ---------------------------------------------------------------------
##
## Copyright (C) 2017 by the molsturm authors
##
## This file is part of molsturm.
##
## molsturm is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published
## by the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## molsturm is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with molsturm. If not, see <http://www.gnu.org/licenses/>.
##
## ---------------------------------------------------------------------

import sys
from ._scf import compute_derived_hartree_fock_energies
#
# Functions to help with printing nicely formatted output
# of the results obtained
#

def print_convergence_summary(hfres,out=sys.stdout, indention=6*" "):
  out.write("\nSCF converged after\n")
  out.write(indention+"{0:5d}  iterations".format(hfres["n_iter"])+"\n")
  out.write(indention+"{0:5d}  matrix applies".format(hfres["n_mtx_applies"])+"\n")

  label_key_map = {
    "Pulay error norm": "final_error_norm",
    "Final ΔE_1e":      "final_1e_energy_change",
    "Final ΔE_total":   "final_tot_energy_change",
  }
  maxlen = max([ len(k) for k in label_key_map ])
  fstr = indention+"{0:"+str(maxlen)+"s}  =  {1:10.8g}\n"

  out.write("\nFinal SCF errors:\n")
  for k in label_key_map:
    out.write(fstr.format(k,hfres[label_key_map[k]]))

def print_mo_occupation(hfres, out=sys.stdout, indention=6*" ", title="Orbital occupation:"):
  """
  Print the MO occupation number and MO energies of the
  calculation results obtained with the function hartree_fock
  above.
  """
  out.write("\n" + title+"\n")

  fstr = indention
  if hfres["restricted"]:
    fstr+="{aocc:1s}  {ena:20}  {bocc:1s}\n"
  else:
    fstr+="{aocc:1s}  {ena:20}  |  {enb:20}  {bocc:1s}\n"
  out.write(fstr.format(aocc="a",bocc="b",enb="",ena=""))

  for i in range(max(hfres["n_orbs_alpha"],hfres["n_orbs_beta"])):
    kw = {
      "aocc": "*" if i < hfres["n_alpha"] else " ",
      "bocc": "*" if i < hfres["n_beta"]  else " ",
      "ena": "",
      "enb": "",
    }
    numfmt="{0:12.8g}"
    if i < hfres["n_orbs_alpha"]:
      kw["ena"] = numfmt.format(hfres["orben_f"][i])
    if i < hfres["n_orbs_beta"]:
      kw["enb"] = numfmt.format(hfres["orben_f"][i+hfres["n_orbs_alpha"]])
    out.write(fstr.format(**kw))

def print_energies(hfres,out=sys.stdout, indention=6*" ", title="Final energies:"):
  """
  Print the energies of a calculation obtained with the
  hartree_fock function above.
  """

  energies = compute_derived_hartree_fock_energies(hfres)
  terms = energies["terms"]

  # Build print format
  maxlen = max([ len(k) for k in terms ])+1
  fstr=indention + "{key:"+str(maxlen)+"s} = {val:15.10g}\n"

  lines = [ "\n" + title + "\n" ]
  lines += [ fstr.format(key=k, val=terms[k]) for k in sorted(terms) ]
  lines.append("\n")
  lines.append(fstr.format(key="E_1e", val=energies["energy_1e"]))
  lines.append(fstr.format(key="E_2e", val=energies["energy_2e"]))
  lines.append(fstr.format(key="E electronic", val=energies["energy_electronic"]))
  lines.append("\n")
  lines.append(fstr.format(key="E_pot", val=energies["energy_potential"]))
  lines.append(fstr.format(key="E_kin", val=energies["energy_kinetic"]))
  lines.append(fstr.format(key="virial ratio", val=energies["virial_ratio"]))
  lines.append("\n")
  lines.append((indention+"{key:"+str(maxlen)+"s} = {val:20.15g}") \
               .format(key="E_total", val=hfres["energy_ground_state"])+"\n")
  out.writelines(lines)


def print_quote(hfres, out=sys.stdout):
  from numpy.random import randint

  # A list of dull Angus MacGyver quotes
  quotes = [ "Lord Cyril Cleeve: [rummaging through the scrolls] Where's the treasure?\n" \
             "Angus MacGyver:    I think you're looking at it.",
             "Atticus: [to MacGyver] You were always my brightest student!",
           ]
  quote = quotes[randint(0,len(quotes))]

  phrase="molsturm out   ...   and now for something completely different"
  width=max([len(phrase)+8] + [ len(line)+2 for line in quote.split("\n") ])
  out.write((width*"="+"\n"))
  out.write( ((width-len(phrase))//2)*" " + phrase + "\n")
  out.write((width*"="+"\n"))
  out.write(quote+"\n")
