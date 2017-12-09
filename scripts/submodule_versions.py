#!/usr/bin/env python3

import os
import re
import yaml
import collections
import argparse


class Config:
    def __init__(self, topdir=os.getcwd()):
        self.topdir = topdir
        self.config_path = os.path.join(self.topdir,
                                        "scripts/submodule_versions.cfg.yaml")
        with open(self.config_path, "r") as f:
            raw_config = yaml.safe_load(f)
        self.repos = raw_config["repos"]
        self.defaults = raw_config["defaults"]

        # Normalise repos:
        for key, repo in self.repos.items():
            repo.setdefault("name", key)

            path = repo.get("path", self.defaults["path"].format(**repo))
            if not os.path.isabs(path):
                path = os.path.realpath(os.path.join(self.topdir, path))
            repo["path"] = path

            repo.setdefault("pinfiles", self.defaults["pinfiles"])
            repo["pinfiles"] = [os.path.join(repo["path"], f) for f in repo["pinfiles"]]

            repo.setdefault("versionfile", self.defaults["versionfile"])
            repo["versionfile"] = os.path.join(repo["path"], repo["versionfile"])

    def find_dependents(self, repo):
        if not isinstance(repo, dict):
            repo = self.repos[repo]
        return [k for k, vdict in self.repos.items()
                if repo["name"] in vdict["depends_on"]]

    def find_dependencies(self, repo):
        if not isinstance(repo, dict):
            repo = self.repos[repo]
        return repo["depends_on"]


def search_repo_version(file, repo_name):
    """
    Search for the repository version in the specified file.
    If it is not found, returns None
    """
    version_re = re.compile(repo_name + ".*([0-9]+\.[0-9]+\.[0-9]+)")
    with open(file, "r") as f:
        for line in f:
            match = version_re.search(line)
            if match:
                return match.group(1)
    return None


def set_repo_version(file, repo_name, version):
    """
    Search for the repository version in the specified file
    and update it to the value supplied.
    """
    version_re = re.compile("(" + repo_name + ".*)([0-9]+\.[0-9]+\.[0-9]+)")
    newlines = []
    with open(file, "r") as f:
        for line in f:
            match = version_re.search(line)
            if match:
                newlines.append(line[:match.start()] + match.group(1) +
                                version + line[match.end():])
            else:
                newlines.append(line)
    with open(file, "w") as f:
        f.write("".join(newlines))


def extract_repo_version(repo):
    """
    Extract the actual version the repo is currently at
    """
    ver = search_repo_version(repo["versionfile"], repo["name"])
    if not ver:
        raise ValueError("Could not extract version for repo " + repo["name"] +
                         " from file " + repo["versionfile"] + ". Check versionfile")
    else:
        return ver


def extract_pinned_dependency_versions(repo):
    """
    Extract the version the repo pinnes its dependencies to.
    """
    pinned = {}
    for dep in repo["depends_on"]:
        for pin in repo["pinfiles"]:
            ver = search_repo_version(pin, dep)
            if not ver:
                continue

            if pinned.get(dep, ver) != ver:
                raise ValueError("Versions in pinfiles of repo " + repo["name"] +
                                 " for dependency " + dep + " do not agree.")
            else:
                pinned[dep] = ver

        if dep not in pinned:
            raise ValueError("Could not find any pinned version for dependency " + dep +
                             " in any of the pinfiles for repo " + repo["name"])
    return pinned


def set_pinned_dependency_version(repo, dependency, version):
    """
    Alter the version found in the pinfiles of repo, such that the version
    of dependency is consistent to version.
    """
    if dependency not in repo["depends_on"]:
        raise ValueError(dependency + " is not a dependency of " + repo["name"])

    any_found = False
    for pin in repo["pinfiles"]:
        ver = search_repo_version(pin, dependency)
        if not ver:
            continue
        any_found = True

        if ver != version:
            set_repo_version(pin, dependency, version)

    if not any_found:
        raise ValueError("Did not find dependency " + dependency +
                         " in any pinfile of repository " + repo["name"])


VersionDifferences = collections.namedtuple("VersionDifferences",
                                            ["repo_name", "dependency_name",
                                             "current_version", "pinned_version"])


def get_differing_versions(config):
    """
    Get the VersionDifferences for the repositories of the config
    """
    # Map from the name of the repository to its current version
    version_cache = {name: extract_repo_version(repo)
                     for name, repo in config.repos.items()}
    return [
        VersionDifferences(name, dep_name, version_cache[dep_name], pin_version)
        for name, repo in config.repos.items()
        for dep_name, pin_version in extract_pinned_dependency_versions(repo).items()
        if pin_version != version_cache[dep_name]
    ]


def check_versions(config):
    differences = get_differing_versions(config)
    if not differences:
        return

    maxlen = max(len(diff.repo_name) for diff in differences)
    for repo_name in config.repos:
        repolist = [
            diff.dependency_name + " (" + diff.pinned_version +
            " -> " + diff.current_version + ")"
            for diff in differences
            if diff.repo_name == repo_name
        ]

        # Skip repos without issues:
        if not repolist:
            continue

        string = ("{0:" + str(maxlen + 1) + "s} ").format(repo_name + ":")
        string += ", ".join(repolist)
        print(string)


def fix_versions(config, only_repos=[]):
    """
    Fix the version discrepancies between the repositories.
    If only_repos is not empty, only the listed repos are edited.
    """
    differences = get_differing_versions(config)
    if not differences:
        return

    if only_repos:
        differences = [d for d in differences if d.repo_name in only_repos]

    for repo_name, dep_name, dep_version, _ in differences:
        set_pinned_dependency_version(config.repos[repo_name], dep_name, dep_version)


def main():
    if not os.path.isfile("src/molsturm/CMakeLists.txt") or \
       not os.path.isfile("CMakeLists.txt"):
        raise SystemExit("Can only be run from the top directory "
                         "of the molsturm repository.")

    parser = argparse.ArgumentParser(
        description="Check and fix versions between the molsturm subprojects."
    )
    parser.add_argument("--fix", action="store_true", default=False,
                        help="Not only show differences but fix them as well")
    args = parser.parse_args()

    config = Config()
    check_versions(config)

    if args.fix:
        fix_versions(config)


if __name__ == "__main__":
    main()
