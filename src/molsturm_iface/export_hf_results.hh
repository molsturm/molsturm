//
// Copyright (C) 2017 by the molsturm authors
//
// This file is part of molsturm.
//
// molsturm is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// molsturm is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with molsturm. If not, see <http://www.gnu.org/licenses/>.
//

#pragma once
#include "HfResults.hh"
#include "config.hh"
#include "gint/ERITensor_i.hh"

namespace molsturm {
namespace iface {

struct ExportParameters {
  double error = 0;
  bool export_repulsion_integrals = true;

  ExportParameters() = default;
  ExportParameters(double error, bool export_repulsion_integrals)
        : error(error), export_repulsion_integrals(export_repulsion_integrals) {}
};

template <typename State>
HfResults export_hf_results(const State& state, const gint::ERITensor_i<scalar_type>& eri,
                            const ExportParameters& params);

}  // namespace iface
}  // namespace molsturm
