#!/usr/bin/env python
# coding: utf-8

# # Exploring the Census Datasets table
# 
# This tutorial demonstrates basic use of the `census_datasets` dataframe that contains metadata of the Census source datasets. This metadata can be joined to the cell metadata dataframe (`obs`) via the column `dataset_id`, 
# 
# **Contents**
# 
# 1. Fetching the datasets table.
# 2. Fetching the expression data from a single dataset.
# 3. Downloading the original source H5AD file of a dataset.
# 
# ⚠️ Note that the Census RNA data includes duplicate cells present across multiple datasets. Duplicate cells can be filtered in or out using the cell metadata variable `is_primary_data` which is described in the [Census schema](https://github.com/chanzuckerberg/cellxgene-census/blob/main/docs/cellxgene_census_schema.md#repeated-data).
# 
# ## Fetching the datasets table
# 
# 
# Each Census contains a top-level dataframe itemizing the datasets contained therein. You can read this into a `pandas.DataFrame`.

# In[1]:


import cellxgene_census

census = cellxgene_census.open_soma()
census_datasets = census["census_info"]["datasets"].read().concat().to_pandas()

# for convenience, indexing on the soma_joinid which links this to other census data.
census_datasets = census_datasets.set_index("soma_joinid")

census_datasets


# The sum cells across all datasets should match the number of cells across all SOMA experiments (human, mouse).

# In[2]:


# Count cells across all experiments
experiments_total_cells = 0
print("Count by experiment:")
for organism_name in census["census_data"].keys():
    num_cells = len(cellxgene_census.get_obs(census, organism_name, column_names=["soma_joinid"]))
    print(f"\t{num_cells} cells in {organism_name}")
    experiments_total_cells += num_cells

print(f"\nFound {experiments_total_cells} cells in all experiments.")

# Count cells across all datasets
print(f"Found {census_datasets.dataset_total_cell_count.sum()} cells in all datasets.")


# ## Fetching the expression data from a single dataset
# 
# Lets pick one dataset to slice out of the census, and turn into an [AnnData](https://anndata.readthedocs.io/en/latest/) in-memory object. This can be used with the [ScanPy](https://scanpy.readthedocs.io/en/stable/) toolchain. You can also save this AnnData locally using the AnnData [write](https://anndata.readthedocs.io/en/latest/api.html#writing) API.

# In[3]:


census_datasets[census_datasets.dataset_id == "0bd1a1de-3aee-40e0-b2ec-86c7a30c7149"]


# Create a query on the mouse experiment, "RNA" measurement, for the dataset_id.

# In[4]:


adata = cellxgene_census.get_anndata(
    census, organism="Mus musculus", obs_value_filter="dataset_id == '0bd1a1de-3aee-40e0-b2ec-86c7a30c7149'"
)

adata


# ## Downloading the original source H5AD file of a dataset.
# 
# You can download the original H5AD file for any given dataset. This is the same H5AD you can download from the [CZ CELLxGENE Discover](https://cellxgene.cziscience.com/), and may contain additional data-submitter provided information which was not included in the Census.
# 
# To do this you can fetch the location in the cloud or directly download to your system using the `cellxgene-census`

# In[5]:


# Option 1: Direct download
cellxgene_census.download_source_h5ad(
    "0bd1a1de-3aee-40e0-b2ec-86c7a30c7149", to_path="Tabula_Muris_Senis-bone_marrow.h5ad"
)


# In[6]:


# Option 2: Get location and download via preferred method
uri = cellxgene_census.get_source_h5ad_uri("0bd1a1de-3aee-40e0-b2ec-86c7a30c7149")
uri

# you can now download the H5AD in shell via AWS CLI e.g. `aws s3 cp uri ./`


# Close the census

# In[7]:


census.close()

