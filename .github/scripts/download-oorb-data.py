# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
sbpy test environment helper
"""

import os
import shutil
import urllib.request

path = os.environ["OORB_DATA"]

base_url = "https://github.com/mkelley/pyoorb-experiment/raw/sbpy-testing/pyoorb/data"

files = [
    "ET-UT.dat",
    "TAI-UTC.dat",
    "de430.dat",
    "OBSCODE.dat",
]

# download files as needed
for file_name in files:
    fn = os.path.join(path, file_name)
    if os.path.exists(fn):
        continue

    url = "/".join((base_url, file_name))
    with urllib.request.urlopen(url) as response:
        with open(fn, "wb") as outf:
            shutil.copyfileobj(response, outf)
        print(url, "→", fn)
