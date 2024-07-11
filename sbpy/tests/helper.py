# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
sbpy test environment helper
"""

import os
import logging
import shutil
import urllib.request
from typing import List


def get_pyoorb_data() -> None:
    """Download pyoorb data for testing."""

    # only proceed if this is a tox environment
    path: str
    try:
        path = os.environ["TOX_ENV_DIR"]
    except KeyError:
        return

    # create necessary directories
    dir: str
    for dir in ["share", "share/oorb"]:
        if not os.path.exists(f"{path}/{dir}"):
            logging.info("sbpy.tests.helper creating directory: %s", f"{path}/{dir}")
            os.mkdir(f"{path}/{dir}")

    # download files as needed
    url: str = (
        "https://github.com/mkelley/pyoorb-experiment/raw/sbpy-testing/pyoorb/data"
    )
    files: List[str] = [
        "ET-UT.dat",
        "TAI-UTC.dat",
        "de430.dat",
        "OBSCODE.dat",
    ]
    for fn in files:
        if os.path.exists(f"{path}/share/oorb/{fn}"):
            continue

        logging.info("sbpy.tests.helper downloading pyoorb data file: %s", fn)
        with urllib.request.urlopen(f"{url}/{fn}") as response:
            with open(f"{path}/share/oorb/{fn}", "wb") as outf:
                shutil.copyfileobj(response, outf)
