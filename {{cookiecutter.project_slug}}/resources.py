# %%
import os
import sys

# %%
from findhere import *
from datatracker import *
tr = Tracker()
set_cloud(True)

import hail as hl

# %%
import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)
