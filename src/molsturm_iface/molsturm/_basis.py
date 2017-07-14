#!/usr/bin/env python3
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

import molsturm_iface as iface

available_basis_types = [ t for t in iface.available_basis_types() ]

def has_real_harmonics(**basis_args):
  """
  Does the basis have real spherical harmonics in the angular part(True)
  or complex ones(False)

  By default complex is assumed.
  """
  basis_type = basis_args["basis_type"]
  if basis_type == "gaussian/libint":
    return True
  elif basis_type.startswith("sturmian/atomic/"):
    m_max = basis_args.get("m_max", -1)
    # For m_max == 0 the harmonics are real because
    # the e^{I m \phi} term is identially 1 for all
    # basis functions.
    return m_max == 0
  else:
    return False

def is_sturmian(**basis_args):
  """
  Is the basis a Sturmian basis.
  """
  return basis_args["basis_type"].startswith("sturmian/")

def is_gaussian(**basis_args):
  """
  Is the basis a Gaussian basis.
  """
  return basis_args["basis_type"].startswith("gaussian/")
