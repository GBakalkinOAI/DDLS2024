#!/usr/bin/env python
# coding: utf-8

# # Generating citations for Census slices
# 
# This notebook demonstrates how to generate a citation string for all datasets contained in a Census slice.
# 
# **Contents**
# 
# 1. Requirements
# 1. Generating citation strings
#    1. Via cell metadata query
#    1. Via an AnnData query 
# 
# ⚠️ Note that the Census RNA data includes duplicate cells present across multiple datasets. Duplicate cells can be filtered in or out using the cell metadata variable `is_primary_data` which is described in the [Census schema](https://github.com/chanzuckerberg/cellxgene-census/blob/main/docs/cellxgene_census_schema.md#repeated-data).
# 
# ## Requirements
# 
# This notebook requires:
# 
# - `cellxgene_census` Python package.
# - Census data release with [schema version](https://github.com/chanzuckerberg/cellxgene-census/blob/main/docs/cellxgene_census_schema.md) 1.3.0 or greater.
# 
# ## Generating citation strings
# 
# First we open a handle to the Census data. To ensure we open a data release with schema version 1.3.0 or greater, we use `census_version="latest"`

# In[1]:


import cellxgene_census

census = cellxgene_census.open_soma(census_version="latest")
census["census_info"]["summary"].read().concat().to_pandas()


# Then we load the dataset table which contains a column `"citation"` for each dataset included in Census. 

# In[2]:


datasets = census["census_info"]["datasets"].read().concat().to_pandas()
datasets["citation"]


# For cross-ref style citations you can look at the column `"collection_doi_label"`

# In[3]:


datasets["collection_doi_label"]


# And now we can use the column `"dataset_id"` present in both the dataset table and the Census cell metadata to create citation strings for any Census slice.
# 
# ### Via cell metadata query

# In[4]:


# Query cell metadata
cell_metadata = cellxgene_census.get_obs(
    census, "homo_sapiens", value_filter="tissue == 'cardiac atrium'", column_names=["dataset_id", "cell_type"]
)

# Get a citation string for the slice
slice_datasets = datasets[datasets["dataset_id"].isin(cell_metadata["dataset_id"])]
print(*set(slice_datasets["citation"]), sep="\n\n")


# In[5]:


print(*set(slice_datasets["collection_doi_label"]), sep="\n\n")


# ### Via AnnData query

# In[6]:


# Fetch an AnnData object
adata = cellxgene_census.get_anndata(
    census=census,
    organism="homo_sapiens",
    measurement_name="RNA",
    obs_value_filter="tissue == 'cardiac atrium'",
    var_value_filter="feature_name == 'MYBPC3'",
    obs_column_names=["dataset_id", "cell_type"],
)

# Get a citation string for the slice
slice_datasets = datasets[datasets["dataset_id"].isin(adata.obs["dataset_id"])]
print(*set(slice_datasets["citation"]), sep="\n\n")


# In[7]:


print(*set(slice_datasets["collection_doi_label"]), sep="\n\n")


# And don't forget to close the Census handle

# In[8]:


census.close()

