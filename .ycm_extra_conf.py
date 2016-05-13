## ---------------------------------------------------------------------
##
## Copyright (C) 2016 by the linalgwrap authors
##
## This file is part of linalgwrap.
##
## linalgwrap is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published
## by the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## linalgwrap is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with linalgwrap. If not, see <http://www.gnu.org/licenses/>.
##
## ---------------------------------------------------------------------

import os

# This file is loosely based upon the file
# cpp/ycm/.ycm_extra_conf.py from the youcompleteme daemon process
# available on github:
# https://github.com/Valloric/ycmd/blob/master/cpp/ycm/.ycm_extra_conf.py


# These are the compilation flags that will be used in case there's no
# compilation database set (by default, one is not set).
flags = [
    # Warnings: For a very detailed discussion about this
    # see the following stackexchange post:
    # https://programmers.stackexchange.com/questions/122608#124574
    '-Wall',
    '-Wextra',
    '-Werror',
    '-pedantic',
    # Generate unwind information
    '-fexceptions',
    # Compile debug code as well
    '-DDEBUG',
    # Compile extra code blocks:
    # '-DLINALGWRAP_HAVE_GLIBC_STACKTRACE',
    '-DLINALGWRAP_HAVE_ARMADILLO',
    # Compile as c++14
    '-std=c++14',
    # Treat .h header files as c++:
    '-x', 'c++',
    # Include other libraries and show errors and 
    # warnings within them
    # To suppress errors shown here, use "-isystem" 
    # instead of "-I"
    '-I', 'src',
    '-isystem', './../linalgwrap/src',
    '-isystem', './../linalgwrap/external/rapidcheck/include',
    '-isystem', './../linalgwrap/external/catch/include',
    '-isystem', './../gint/src',
    '-isystem', './../gscf/src',
]

def DirectoryOfThisScript():
  return os.path.dirname( os.path.abspath( __file__ ) )

SOURCE_EXTENSIONS = [ '.cpp', '.cxx', '.cc', '.c', '.C' ]

def MakeRelativePathsInFlagsAbsolute( flags, working_directory ):
  if not working_directory:
    return list( flags )
  new_flags = []
  make_next_absolute = False
  path_flags = [ '-isystem', '-I', '-iquote', '--sysroot=' ]
  for flag in flags:
    new_flag = flag

    if make_next_absolute:
      make_next_absolute = False
      if not flag.startswith( '/' ):
        new_flag = os.path.join( working_directory, flag )

    for path_flag in path_flags:
      if flag == path_flag:
        make_next_absolute = True
        break

      if flag.startswith( path_flag ):
        path = flag[ len( path_flag ): ]
        new_flag = path_flag + os.path.join( working_directory, path )
        break

    if new_flag:
      new_flags.append( new_flag )
  return new_flags


def IsHeaderFile( filename ):
  extension = os.path.splitext( filename )[ 1 ]
  return extension in [ '.h', '.hxx', '.hpp', '.hh' ]

def FlagsForFile( filename, **kwargs ):
  relative_to = DirectoryOfThisScript()
  final_flags = MakeRelativePathsInFlagsAbsolute( flags, relative_to )

  return {
    'flags': final_flags,
    'do_cache': True
  }
