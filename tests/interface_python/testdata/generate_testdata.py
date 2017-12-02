#!/usr/bin/env python3
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

import glob
import molsturm
import molsturm.posthf
import molsturm.yaml_utils
import os
import yaml


SKIPPED_SINCE_DONE = "Skipped since already done."
HASH_EXT = "sha256"  # The extension to check for if check_for_hash


def build_input_params():
    """Gather a dictionary of all input parameters"""
    dir_of_this_script = os.path.dirname(os.path.abspath(__file__))

    # Build the list of input data
    inputs = dict()
    for f in glob.iglob(os.path.join(dir_of_this_script, "*.in.yaml")):
        fn, _ = os.path.splitext(f)
        key = fn[:-3]

        with open(f, "r") as stream:
            inputs[key] = yaml.safe_load(stream)
    return inputs


def print_file_start(filename):
    print("#")
    print("# " + os.path.basename(filename))


def print_status(key, message):
    print("{0:15s} {1}".format(key + ":", message))


# Cache of calculations we already did
calculation_scf_cache = dict()


def run_scf_calculation(name, scfparams):
    if name not in calculation_scf_cache:
        print_status("run_hf", "Running HF calculation")
        calculation_scf_cache[name] = molsturm.self_consistent_field(scfparams)
    return calculation_scf_cache[name]


# --------------------------------------------------------------------


def job_dump_yaml(name, output, scfparams, dump_params):
    """Run a full calculation and dump the result as a yaml file"""
    res = run_scf_calculation(name, scfparams)

    # Remove keys which are given by the parameters
    for key in dump_params.get("remove_keys", []):
        if key in res:
            del res[key]
    molsturm.dump_state(res, output, type="yaml")


def job_dump_hdf5(name, output, scfparams, dump_params):
    """Run a full calculation and dump the result as a yaml file"""
    res = run_scf_calculation(name, scfparams)

    # Remove keys which are given by the parameters
    for key in dump_params.get("remove_keys", []):
        if key in res:
            del res[key]
    molsturm.dump_state(res, output, type="hdf5")


def job_posthf_mp2(name, output, scfparams, mp_params):
    hfres = run_scf_calculation(name, scfparams)
    print_status("posthf_mp2", "Running MP2")
    mp2 = molsturm.posthf.mp2(hfres, **mp_params)

    molsturm.yaml_utils.install_representers()
    with open(output, "w") as f:
        yaml.safe_dump(mp2, f)


def job_posthf_fci(name, output, scfparams, fci_params):
    hfres = run_scf_calculation(name, scfparams)
    print_status("posthf_fci", "Running Full-CI")
    fci = molsturm.posthf.fci(hfres, **fci_params)

    molsturm.yaml_utils.install_representers()
    with open(output, "w") as f:
        yaml.safe_dump(fci, f)


job_dump_yaml.extension = "hf.yaml"
job_dump_hdf5.extension = "hf.hdf5"
job_posthf_mp2.extension = "mp2.yaml"
job_posthf_fci.extension = "fci.yaml"

# --------------------------------------------------------------------


def run_job(name, jobfunction, scfparams, jobparams, check_for_hash):
    """
    Common stuff to do before and after a payload run for each job
    """
    output = name + "." + jobfunction.extension
    checkfile = output if not check_for_hash else output + "." + HASH_EXT

    if not os.path.exists(checkfile):
        jobfunction(name, output, scfparams, jobparams)

        # Create an empty checkfile if not yet done:
        if not os.path.exists(checkfile):
            open(checkfile, "w").close()
    else:
        # Name of the jobfunction with the "job_" prefix
        # removed is the jobname
        jobname = jobfunction.__name__[4:]
        print_status(jobname, message=SKIPPED_SINCE_DONE)


def work_on_case(name, scfparams, generator_params):
    print_file_start(name)
    store_on_server = generator_params.get("store_on_webserver", False)

    for job in generator_params["include"]:
        if isinstance(job, str):
            jobname = job
            jobparams = {}
        elif isinstance(job, dict):
            jobname, jobparams = job.popitem()
        else:
            raise SystemExit("Invalid entry in include list of type '" + type(job) +
                             "' found.")

        # Try to find the jobfunction and then execute it
        try:
            jobfunction = globals()["job_" + jobname]
        except KeyError:
            raise SystemExit("Unknown job name '" + jobname + "' in input file '" +
                             name + ".in.yaml'")

        # Run the job, use the file with the hash extension to see whether
        # running the payload jobfunction is needed if we will store
        # the file on the server (because in that case the actual file
        # will sit somewhere else in the filesystem)
        run_job(name, jobfunction, scfparams, jobparams, check_for_hash=store_on_server)


def main():
    inputs = build_input_params()
    for name in sorted(inputs):
        # Build the ScfParameters object to run the
        # SCF for producing the reference data.
        scfparams = molsturm.ScfParameters.from_dict(inputs[name]["input_parameters"])
        scfparams.normalise()

        # Work on the case and pass the generator arguments:
        work_on_case(name, scfparams, inputs[name]["generator"])


if __name__ == "__main__":
    main()
