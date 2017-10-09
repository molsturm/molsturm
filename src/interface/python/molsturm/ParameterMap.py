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


class ParameterMap(dict):
    """
    Tree-like parameter dictionary mapping from a string to any value.

    If [], i.e. __getitem__ is used a ParameterMap is returned iff
    the string only points to a subtree instead of an actual value.
    """
    def __from_dict_inner(self, d, prefix):
        for k, v in d.items():
            fullkey = prefix + "/" + k
            if isinstance(v, dict):
                self.__from_dict_inner(v, fullkey)
            else:
                self[fullkey] = v

    @classmethod
    def from_dict(cls, d):
        params = ParameterMap()
        params.__from_dict_inner(d, "")
        return params

    def __copy__(self):
        """
        Return a shallow copy of a ParameterMap
        """
        return ParameterMap.from_dict(self)

    def copy(self):
        return self.__copy__()

    def __init__(self, init_values={}):
        """
        Construct a ParameterMap and optionally initialise some inner
        values using an initial dictionary init_values
        """
        if len(init_values) > 0:
            self.__from_dict_inner(init_values, "")

    def __recursive_setter(self, name, key, *args, **kwargs):
        if not isinstance(key, str):
            raise TypeError("key may only be a string")

        # Normalise and split into first part of the key
        # (before the /) and the second part. If only one part
        # is present, duplicate such that both key and rest are the same.
        splitted = key.strip("/").split("/", maxsplit=1)
        first, rest = (2 * splitted)[:2]

        namefct = getattr(super(), name)
        if len(splitted) == 2:
            subdict = super().setdefault(first, ParameterMap())
            namefct = getattr(subdict, name)
        return namefct(rest, *args, **kwargs)

    def __recursive_apply(self, name, key, *args, **kwargs):
        if not isinstance(key, str):
            raise TypeError("key may only be a string")

        splitted = key.strip("/").split("/", maxsplit=1)
        first, rest = (2 * splitted)[:2]

        namefct = getattr(super(), name)
        if len(splitted) == 2:
            try:
                subdict = super().__getitem__(first)
            except KeyError:
                raise KeyError(key)
            namefct = getattr(subdict, name)
        return namefct(rest, *args, **kwargs)

    def __recursive_delete(self, name, key, *args, **kwargs):
        if not isinstance(key, str):
            raise TypeError("key may only be a string")

        splitted = key.strip("/").split("/", maxsplit=1)
        first, rest = (2 * splitted)[:2]

        subdict = None
        namefct = getattr(super(), name)
        if len(splitted) == 2:
            try:
                subdict = super.__getitem__(first)
            except KeyError:
                raise KeyError(key)
            namefct = getattr(subdict, name)
        ret = namefct(rest, *args, **kwargs)

        if isinstance(subdict, ParameterMap) and len(subdict) == 0:
            super().__delitem__(key)
        return ret

    def __getitem__(self, key):
        return self.__recursive_apply("__getitem__", key)

    def __setitem__(self, key, value):
        return self.__recursive_setter("__setitem__", key, value)

    def __delitem__(self, key):
        return self.__recursive_delete("__delitem__", key)

    def __contains__(self, key):
        return self.__recursive_apply("__contains__", key)

    def __iter__(self):
        for k in super().__iter__():
            if isinstance(super().__getitem__(k), ParameterMap):
                for i in super().__getitem__(k):
                    yield k + "/" + i
            else:
                yield k

    def __len__(self):
        return super().__len__()
        pass

    def setdefault(self, key, value):
        return self.__recursive_setter("setdefault", key, value)

    def update(self, E=None, **F):
        if E is not None:
            if hasattr(E, "keys"):
                for k in E:
                    self.__setitem__(k, E[k])
            else:
                for k, v in E:
                    self.__setitem__(k, E[k])
        for k in F:
            self.__setitem__(k, F[k])

    def get(self, k, d=None):
        try:
            return self.__getitem__(k)
        except KeyError:
            return d

    def keys_recursive(self):
        return list(self.__iter__())

    def to_dict(self):
        """
        Return a dict of dicts which represents the same data
        """
        return {
            k: p.to_dict() if isinstance(p, ParameterMap) else p
            for k, p in super().items()
        }

    # TODO There is still stuff missing ...


if __name__ == "__main__":
    m = ParameterMap()
    m["/a/b/c"] = 5
    print(m["a"])
    try:
        print(m["b/c"])
    except KeyError as e:
        print(str(e))
    print(m)

    for k in m["a"]:
        print(k)

    m.setdefault("/a/b/c", 6)
    m.setdefault("/a/d/c", 6)

    print(m)

    m.update({"/a/b/c": 7, "d/e": 9})
    print(m)
