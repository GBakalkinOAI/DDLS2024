#!/usr/bin/env python
# coding: utf-8

# # Out-of-core (incremental) mean and variance calculation
# 
# This tutorial describes use of the cellxgene_census.experimental.pp API for calculating out-of-core mean and variance in the Census. The variance calculation is performed using [Welford's online algorithm](https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance).
# 
# **Contents**
# 
# 1. The mean and variance API.
# 2. Example: calculate mean and variance for a slice of the Census.
# 
# ## The mean and variance API
# 
# `mean_variance()` calculates the mean and the variance for an `ExperimentAxisQuery`. The following additional arguments are supported:
# 
# - `layer`: the X layer used for the calculation, defaults to `raw`
# - `axis`: the axis along which the calculation is performed. Specify 0 for the `var` axis and 1 for the `obs` axis
# - `calculate_mean`: if False, do not include the mean in the result
# - `calculate_variance`: if False, do not compute the variance.
# - `ddof`: _Delta Degrees of Freedom_: the divisor used in the calculation for variance is N - ddof, where N represents the number of elements.
# 
# 

# In[1]:


# Import packages
import cellxgene_census
import pandas as pd
import tiledbsoma as soma
from cellxgene_census.experimental.pp import mean_variance


# ## Example: calculate mean and variance for a slice of the Census

# As an example, we'll calculate the mean and variance along the `obs` axis for a subset of cells from the Mouse census.
# 
# The return value will be a Pandas dataframe indexed by `soma_joinid` (in this case, it will be relative to `obs`) and will contain the `mean` and `variance` columns.

# In[2]:


experiment_name = "mus_musculus"
obs_value_filter = 'is_primary_data == True and tissue_general == "skin of body"'

with cellxgene_census.open_soma(census_version="stable") as census:
    with census["census_data"][experiment_name].axis_query(
        measurement_name="RNA", obs_query=soma.AxisQuery(value_filter=obs_value_filter)
    ) as query:
        mv_df = mean_variance(
            query,
            axis=1,
            calculate_mean=True,
            calculate_variance=True,
        )

        obs_df = query.obs().concat().to_pandas()

mv_df


# We can now concatenate the resulting dataframe to `obs`:

# In[3]:


combined_df = pd.concat([obs_df.set_index("soma_joinid"), mv_df], axis=1)
combined_df

