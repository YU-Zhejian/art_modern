#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Format boost imports.

"""

import sys
from collections import defaultdict


if __name__ == "__main__":
    # The entries will be like 
    # #include <boost/algorithm/string/join.hpp> // src/libam_support/out/HeadlessBamReadOutput.cc
    d = defaultdict(list)
    for line in sys.stdin:
        line = line.strip()
        header = line.split("//")[0]
        file = line.split("//")[1].strip()
        d[header].append(file)
    for header in sorted(d.keys()):
        print(f"{header}")
        for file in sorted(d[header]): 
            print(f"{' ' * (len(header))}: {file}")

