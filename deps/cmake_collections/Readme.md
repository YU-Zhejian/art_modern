# Readme

## Copying

## Files Adapted From [GNU AutoConf](https://www.gnu.org/software/autoconf/)

The following files are adapted from GNU AutoConf Archive version [2.72](https://ftp.gnu.org/gnu/autoconf/autoconf-2.72.tar.xz).

The following files are adapted from `/autoconf-2.72/lib/autoconf/c.m4`:

- src/test_c11.c
- src/test_c11_globals.h
- src/test_c89.c
- src/test_c89_globals.h
- src/test_c99.c
- src/test_c99_globals.h

The copyright notice is as follows:

```text
# This file is part of Autoconf.			-*- Autoconf -*-
# Programming languages support.
# Copyright (C) 2001-2017, 2020-2021 Free Software Foundation, Inc.

# This file is part of Autoconf.  This program is free
# software; you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the
# Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# Under Section 7 of GPL version 3, you are granted additional
# permissions described in the Autoconf Configure Script Exception,
# version 3.0, as published by the Free Software Foundation.
#
# You should have received a copy of the GNU General Public License
# and a copy of the Autoconf Configure Script Exception along with
# this program; see the files COPYINGv3 and COPYING.EXCEPTION
# respectively.  If not, see <https://www.gnu.org/licenses/>.

# Written by David MacKenzie, with help from
# Akim Demaille, Paul Eggert,
# Franc,ois Pinard, Karl Berry, Richard Pixley, Ian Lance Taylor,
# Roland McGrath, Noah Friedman, david d zuhn, and many others.
```

## Files Adapted From [GNU AutoConf Archive](https://www.gnu.org/software/autoconf-archive/)

The following files are adapted from GNU AutoConf Archive version [2024.10.16](https://ftp.gnu.org/gnu/autoconf-archive/autoconf-archive-2024.10.16.tar.xz).

The following files are adapted from `/autoconf-archive-2024.10.16/m4/ax_cxx_compile_stdcxx.m4`:

- src/test_cxx11.cc
- src/test_cxx11_ax_globals.hh
- src/test_cxx11_globals.hh
- src/test_cxx14.cc
- src/test_cxx14_ax_globals.hh
- src/test_cxx17.cc
- src/test_cxx17_ax_globals.hh
- src/test_cxx20.cc
- src/test_cxx20_ax_globals.hh
- src/test_cxx98.cc
- src/test_cxx98_globals.hh

The copyringt notice is as follows:

```text
# ===========================================================================
#  https://www.gnu.org/software/autoconf-archive/ax_cxx_compile_stdcxx.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_CXX_COMPILE_STDCXX(VERSION, [ext|noext], [mandatory|optional])
#
# DESCRIPTION
#
#   Check for baseline language coverage in the compiler for the specified
#   version of the C++ standard.  If necessary, add switches to CXX and
#   CXXCPP to enable support.  VERSION may be '11', '14', '17', '20', or
#   '23' for the respective C++ standard version.
#
#   The second argument, if specified, indicates whether you insist on an
#   extended mode (e.g. -std=gnu++11) or a strict conformance mode (e.g.
#   -std=c++11).  If neither is specified, you get whatever works, with
#   preference for no added switch, and then for an extended mode.
#
#   The third argument, if specified 'mandatory' or if left unspecified,
#   indicates that baseline support for the specified C++ standard is
#   required and that the macro should error out if no mode with that
#   support is found.  If specified 'optional', then configuration proceeds
#   regardless, after defining HAVE_CXX${VERSION} if and only if a
#   supporting mode is found.
#
# LICENSE
#
#   Copyright (c) 2008 Benjamin Kosnik <bkoz@redhat.com>
#   Copyright (c) 2012 Zack Weinberg <zackw@panix.com>
#   Copyright (c) 2013 Roy Stogner <roystgnr@ices.utexas.edu>
#   Copyright (c) 2014, 2015 Google Inc.; contributed by Alexey Sokolov <sokolov@google.com>
#   Copyright (c) 2015 Paul Norman <penorman@mac.com>
#   Copyright (c) 2015 Moritz Klammler <moritz@klammler.eu>
#   Copyright (c) 2016, 2018 Krzesimir Nowak <qdlacz@gmail.com>
#   Copyright (c) 2019 Enji Cooper <yaneurabeya@gmail.com>
#   Copyright (c) 2020 Jason Merrill <jason@redhat.com>
#   Copyright (c) 2021 Jörn Heusipp <osmanx@problemloesungsmaschine.de>
#   Copyright (c) 2015, 2022, 2023, 2024 Olly Betts
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved.  This file is offered as-is, without any
#   warranty.
```
