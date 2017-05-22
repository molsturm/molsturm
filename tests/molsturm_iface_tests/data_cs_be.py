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
## vi: tabstop=2 shiftwidth=2 softtabstop=2 expandtab

import numpy as np

#
# Parameters for the run:
#
params = {
  "atom_numbers": [4],
  "coords":       [ [0,0,0] ],
  #
  "basis_type":     "atomic/cs_naive",
  "n_max":          3,
  "l_max":          2,
  "k_exp":          1.0,
  #
  "n_eigenpairs":   4,
  "error":          1e-9,
  "eigensolver":    "lapack",
  "guess_method":   "loewdin",
  "guess_esolver":  "lapack",
}

ref_convergence_result = {
  "n_bas":          14,
  "n_beta":         sum(params["atom_numbers"])//2,
  "n_alpha":        sum(params["atom_numbers"])//2,
  "n_orbs_beta":    params["n_eigenpairs"],
  "n_orbs_alpha":   params["n_eigenpairs"],
  "restricted":     True,
}

ref_n_iter = 13

ref_energies = {
  "energy_coulomb":              5.203944879648532,
  "energy_exchange":            -1.706882550460353,
  "energy_total":                -11.6237892481597,
  "energy_kinetic":              5.715093431464893,
  "energy_nuclear_attraction":  -20.83594500881277,
  "energy_nuclear_repulsion":    0.0,
}

ref_orbital_energies = np.array(2*[ -3.758063817830864, -0.3052996416547051,
                                     0.0704507645679821, 0.07045076456798217 ])
ref_coefficients = np.array([
  2*[-1.187354936022534,-0.272225607909177,-1.25614867939377e-17,-1.065780294534109e-18],
  2*[-1.068207213226528,0.9138515042491057,-2.512297358857309e-17,3.985396843869585e-17],
  2*[6.149715051798564e-18,-6.359308318245655e-17,-0.7670948435528844,0.4383497235019357],
  2*[0,0,-2.328726980435578e-07,9.233047072446519e-09],
  2*[-1.953585617658348e-19,-2.799409569385237e-17,0.4383497235019226,0.7670948435529125],
  2*[-0.5863470313711147,0.2315904281032498,1.428739539910476e-16,-7.442714201147541e-17],
  2*[5.574400463559955e-19,6.517221685441952e-18,0.2001337438555429,-0.114364699515069],
  2*[0,0,6.075609195179478e-08,-2.408886323982845e-09],
  2*[-3.55922209274283e-19,-8.687634987178679e-19,-0.1143646995150655,-0.2001337438555502],
  2*[0,0,3.1655448486224e-18,-2.361984507236245e-18],
  2*[0,0,3.873092284169681e-18,0],
  2*[-4.368646933475374e-14,1.079848195402297e-10,-9.087532021389557e-18,1.532680547196292e-18],
  2*[0,0,-5.228186199318049e-18,0],
  2*[0,0,7.966901144546446e-18,1.394176490338577e-17],
]).transpose()

def __make_ref_fock_aa():
  n_a = params["n_eigenpairs"]
  m = np.zeros(n_a*n_a)
  m[::n_a+1] = ref_orbital_energies[:n_a]
  return m.reshape(n_a,n_a)

ref_fock = { "aa": __make_ref_fock_aa(), "bb": __make_ref_fock_aa(), }

def __make_ref_repulsion_integrals_aaaa():
  n_a = params["n_eigenpairs"]
  return np.array([ 1.3201564829887236, -0.051927456174132254,
   1.0191882208517448e-17, -6.932670580196954e-18, -0.05192745617413223,
   0.4593375841694676, 1.698530801552094e-17, -1.0424022896950529e-17,
   1.0191882208517446e-17, 1.6985308015520938e-17, 0.39591672707582903,
   5.0306980803327406e-17, -6.9326705801969584e-18,
   -1.042402289695052e-17, 5.0306980803327406e-17, 0.39591672707582914,
   -0.05192745617413226, 0.011792639556088124, 7.593597510782904e-19,
   7.55534723624054e-19, 0.011792639556088122, 0.0001658252883217792,
   -6.206011871760161e-20, 8.758550173346611e-19,
   1.2039331081855872e-18, 5.328660528663701e-19, 0.0008587865083350711,
   3.7947076036992655e-19, -8.746295139447042e-19,
   -1.3056238171601812e-18, 2.168404344971009e-19, 0.0008587865083350716,
   1.0191882208517443e-17, 7.593597510782887e-19, 0.007637832871897834,
   0.004500667338809478, 1.2039331081855856e-18, 5.750816317314391e-18,
   0.007397011869553222, 0.004358761219897184, 0.008865240970987173,
   0.008585719770082817, 4.867074794059394e-18, 1.9147419675463928e-19,
   2.293746421337852e-18, 1.7574921538096497e-18,
   -1.0119584207858949e-19, 4.254666921509515e-18,
   -6.9326705801969515e-18, 7.555347236240549e-19, 0.004500667338809478,
   -0.007637832871897724, -8.746295139447029e-19, -2.699425808599076e-18,
   0.004358761219897184, -0.0073970118695531155, 2.279300163863037e-18,
   1.9905898152677444e-18, -2.2048170406582718e-18,
   -2.6950445618453333e-19, 0.008865240970987178, 0.008585719770082822,
   2.125088405303538e-19, -2.2333249555915615e-18, -0.051927456174132344,
   0.011792639556088093, 1.203933108185589e-18, -8.746295139447067e-19,
   0.011792639556088093, 0.00016582528832179155, 5.328660528663657e-19,
   -1.3056238171601756e-18, 7.593597510782918e-19,
   -6.206011871760581e-20, 0.0008587865083349958, 3.7947076036992655e-19,
   7.555347236240526e-19, 8.75855017334666e-19, 2.981555974335137e-19,
   0.0008587865083349964, 0.4593375841694676, 0.00016582528832178515,
   5.750816317314397e-18, -2.6994258085990804e-18,
   0.00016582528832178797, 0.36314078828371416, 1.0644636682373795e-17,
   -3.3230774782505913e-18, 5.750816317314387e-18,
   1.0644636682373797e-17, 0.32829995916392174, 2.42861286636753e-17,
   -2.6994258085990765e-18, -3.323077478250587e-18,
   2.7755575615628914e-17, 0.3282999591639219, 1.6985308015520935e-17,
   -6.20601187176043e-20, 0.007397011869553207, 0.004358761219897182,
   5.328660528663659e-19, 1.0644636682373797e-17, 0.0610161004769595,
   0.03595433091611151, 0.008585719770082788, 0.07082145458150146,
   -2.0000832530521282e-18, 6.6329225360475585e-19,
   -6.085408378962717e-18, 9.475441901805674e-18, 3.823770565503596e-18,
   5.8687378858838654e-18, -1.042402289695052e-17,
   8.758550173346625e-19, 0.004358761219897179, -0.007397011869553109,
   -1.3056238171601833e-18, -3.323077478250585e-18, 0.03595433091611152,
   -0.06101610047695864, -3.4565891887972663e-18, 7.433076717440049e-18,
   -4.675391685530693e-18, 2.9214451885946923e-18, 0.00858571977008282,
   0.0708214545815015, -5.648964257107835e-19, -4.451929434896007e-18,
   1.0191882208517437e-17, 1.2039331081855856e-18, 0.008865240970987188,
   -7.806255641895632e-18, 7.593597510782881e-19, 5.750816317314391e-18,
   0.008585719770082817, 5.692061405548898e-18, 0.0076378328718978445,
   0.00739701186955322, 4.867074794059396e-18, -1.0119584207858983e-19,
   0.0045006673388094925, 0.00435876121989718, 1.9147419675463896e-19,
   4.254666921509515e-18, 1.698530801552093e-17, 5.328660528663693e-19,
   0.008585719770082826, -6.071532165918825e-18, -6.206011871760373e-20,
   1.0644636682373791e-17, 0.07082145458150146, 2.6020852139652106e-18,
   0.007397011869553218, 0.0610161004769595, -2.0000832530521247e-18,
   3.823770565503596e-18, 0.004358761219897205, 0.035954330916111506,
   6.632922536047553e-19, 5.868737885883868e-18, 0.395916727075829,
   0.0008587865083350888, 4.8670747940593935e-18, -2.204817040658269e-18,
   0.0008587865083349722, 0.3282999591639216, -2.000083253052129e-18,
   -4.6753916855306996e-18, 4.867074794059394e-18,
   -2.0000832530521348e-18, 0.3205994401148098, 0.007826434364428984,
   -2.2048170406582718e-18, -4.675391685530698e-18, 0.007826434364429011,
   0.29403582491804736, -1.0704883119742235e-17, 4.502234862586063e-17,
   1.9147419675463456e-19, -2.695044561845294e-19,
   4.502234862586061e-17, 1.5886097591706314e-17, 6.632922536047577e-19,
   2.92144518859469e-18, -1.011958420785949e-19, 3.8237705655035985e-18,
   0.007826434364428989, -0.013281807598381132, 2.1250884053035677e-19,
   -5.648964257107848e-19, 0.022505412805277986, -0.007826434364428826,
   -6.932670580196954e-18, -8.746295139447021e-19,
   1.5178830414797062e-18, 0.008865240970987181, 7.555347236240555e-19,
   -2.6994258085990773e-18, 1.8973538018496328e-18, 0.008585719770082817,
   0.0045006673388094795, 0.004358761219897181, -2.2048170406582683e-18,
   2.1250884053035335e-19, -0.007637832871897725, -0.00739701186955311,
   -2.6950445618453323e-19, -2.2333249555915576e-18,
   -1.0424022896950527e-17, -1.3056238171601812e-18,
   6.938893903907228e-18, 0.008585719770082828, 8.758550173346604e-19,
   -3.323077478250585e-18, 1.734723475976807e-17, 0.07082145458150152,
   0.004358761219897198, 0.035954330916111506, -4.6753916855306965e-18,
   -5.648964257107841e-19, -0.007397011869553107, -0.06101610047695863,
   2.9214451885946916e-18, -4.4519294348960094e-18,
   7.843662647321178e-18, 3.796449982618446e-17,
   -1.0119584207859093e-19, 2.1250884053035485e-19,
   3.796449982618446e-17, 1.2770138435785416e-17,
   3.8237705655035985e-18, -5.6489642571078585e-19,
   1.9147419675463742e-19, 6.6329225360475845e-19, 0.007826434364429015,
   0.02250541280527798, -2.6950445618453212e-19, 2.9214451885946892e-18,
   -0.01328180759838113, -0.007826434364428798, 0.3959167270758293,
   0.0008587865083351027, 4.2546669215095164e-18,
   -2.2333249555915607e-18, 0.0008587865083350693, 0.3282999591639218,
   5.868737885883866e-18, -4.45192943489601e-18, 4.254666921509506e-18,
   5.868737885883872e-18, 0.29403582491804725, -0.007826434364428841,
   -2.233324955591556e-18, -4.4519294348960094e-18,
   -0.007826434364428841, 0.3205994401148097]).reshape((n_a,n_a,n_a,n_a))
ref_repulsion_integrals = {
  "aaaa": __make_ref_repulsion_integrals_aaaa(),
  "aabb": __make_ref_repulsion_integrals_aaaa(),
  "bbaa": __make_ref_repulsion_integrals_aaaa(),
  "bbbb": __make_ref_repulsion_integrals_aaaa(),
}
