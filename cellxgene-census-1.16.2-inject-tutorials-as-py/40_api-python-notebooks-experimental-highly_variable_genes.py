#!/usr/bin/env python
# coding: utf-8

# # Experimental Highly Variable Genes API
# 
# This tutorial describes use of the `cellxgene_census.experimental.pp` API for finding highly variable genes (HVGs) in the Census. The HVG algorithm implements the ranked normalized variance method `seurat_v3` described in [scanpy.pp.highly_variable_genes](https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.highly_variable_genes.html#scanpy.pp.highly_variable_genes).
# 
# There are two API available:
# 
# * `get_highly_variable_genes()` - high level function which accepts arguments similar to `cellxgene_census.get_anndata()`, and returns annotations for each `var` feature in a Pandas DataFrame.
# * `highly_variable_genes()` - lower level function which accepts a `tiledbsoma.ExperimentAxisQuery` and returns the same result.
# 
# Both functions accept common arguments to control ranking, with argument semantics matching the Scanpy API:
# 
# * `n_top_genes` - number of genes to rank.
# * `batch_key` - if specified, normalized ranking will be done in separate batches based upon the obs column value name specified, and then merged into the final result.
# * `span` - the fraction of the data (cells) used when estimating the variance in the [loess model fit](https://has2k1.github.io/scikit-misc/stable/generated/skmisc.loess.loess_model.html#skmisc.loess.loess_model).
# 
# In addition:
# 
# * `max_lowess_jitter` - maxmimum jitter (noise) to data if LOESS fails. Disable by setting to zero.
# 
# For more information, see the docstrings for both functions (e.g. `help(function)`)

# In[1]:


# Import packages
import cellxgene_census
import pandas as pd
import tiledbsoma as soma
from cellxgene_census.experimental.pp import (
    get_highly_variable_genes,
    highly_variable_genes,
)


# ## get_highly_variable_genes
# 
# This convenience function will meet most use cases, and is a wrapper around `highly_variable_genes`.  This demonstration requests the top 500 genes from the Mouse census where `tissue_general` is `heart`, and joins with the `var` dataframe.
# 
# The HVGs returned by get_highly_variable_genes are indexed by their `soma_joinid`.  Join with the `var` dataframe to have a merged view of var metadata.

# In[2]:


with cellxgene_census.open_soma(census_version="stable") as census:
    hvgs_df = get_highly_variable_genes(
        census,
        organism="mus_musculus",
        n_top_genes=500,
        obs_value_filter="""is_primary_data == True and tissue_general == 'heart'""",
    )

    # while the Census is open, also grab the var dataframe for the mouse
    var_df = cellxgene_census.get_var(census, "mus_musculus")

hvgs_df


# Concat the two dataframes for convenience:

# In[3]:


combined_df = pd.concat([var_df.set_index("soma_joinid"), hvgs_df], axis=1)
combined_df


# Select _only_ the highly_variable genes by using the `highly_variable` column value:

# In[4]:


combined_df[combined_df.highly_variable]


# ## highly_variable_genes
# 
# This API provides the same function as `get_highly_variable_genes`, but accepts any `tiledbsoma.ExperimentAxisQuery`.  It is intended for more advanced users who wish to use create and manage their own queries.

# In[5]:


with cellxgene_census.open_soma(census_version="stable") as census:
    experiment = census["census_data"]["mus_musculus"]
    with experiment.axis_query(
        measurement_name="RNA",
        obs_query=soma.AxisQuery(value_filter="""is_primary_data == True and tissue_general == 'heart'"""),
    ) as query:
        hvgs_df = highly_variable_genes(query, n_top_genes=500)

hvgs_df[hvgs_df.highly_variable]

