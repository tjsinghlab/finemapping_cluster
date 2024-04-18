# %%
import pandas as pd
pd.set_option("display.max_colwidth", 80)
pd.set_option("display.precision", 3)

# %%
from {{ cookiecutter.project_slug }} import *
from findhere import *
from datatracker import *
from logzero import *
tr = Tracker()

# %%
tdf = tr.table

# %%