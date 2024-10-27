#!/usr/bin/env python
# coding: utf-8

# # Understanding and filtering out duplicate cells
# 
# This tutorial provides an explanation for the existence of duplicate cells in the Census, and it showcases different ways to handle these cells when performing queries on the Census using the `is_primary_data` cell metadata variable. 
# 
# **Contents**
# 
# 1. Why are there duplicate cells in the Census?
# 2. An example: duplicate cells in the Tabula Muris Senis data.
# 3. Filtering out duplicates cells.
#    1. Filtering out duplicate cells when reading the `obs` data frame.
#    2. Filtering out duplicate cells when creating an AnnData.
#    3. Filtering out duplicate cells for out-of-core operations.
#    
# ## Why are there duplicate cells in the Census?
# 
# Duplicate cells are labeled on the `is_primary_data` cell metadata variable as `False`. To learn more about this please take a look at the corresponding [section of the dataset schema](https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/3.0.0/schema.md#is_primary_data). 
# 
# The Census data is a concatenation of most RNA data from CZ CELLxGENE Discover and these data are ingested one dataset at a time. You can take a look at what data is included in the Census [here](https://chanzuckerberg.github.io/cellxgene-census/cellxgene_census_docsite_schema.html).
# 
# In some cases data from the same cell exists in different datasets, therefore cells can be duplicated throughout CELLxGENE Discover and by extension the Census. 
# 
# The following are a few examples where cells are duplicated in CELLxGENE Discover:
# 
# * There are datasets that combine data from other, pre-existing datasets.
# 
# > *For example [Tabula Sapiens](https://cellxgene.cziscience.com/collections/e5f58829-1a66-40b5-a624-9046778e74f5) has one dataset with all of its cells and separate datasets with cells divided by high-level lineage (i.e. immune, epithelial, stromal, endothelial)*
# 
# * A dataset may provide a meta-analysis of pre-existing datasets.
# 
# > *For example [Jin et al.](https://cellxgene.cziscience.com/collections/b9fc3d70-5a72-4479-a046-c2cc1ab19efc) performed a meta-analysis of COVID-19 data, and they included both the individual datasets as well as one concatenated dataset*
# 
# The Census has all of these data to allow for the execution of dataset-based queries, which would be otherwise be limited if only non-duplicate cells were included.
# 
# ## An example: duplicate cells in the Tabula Muris Senis data
# 
# Let's take a look at an example from the Census using the Tabula Muris Senis data. Some of its datasets contain duplicated cells.
# 
# We can obtain cell metadata for the **main** Tabula Muris Senis dataset: "All - A single-cell transcriptomic atlas characterizes ageing tissues in the mouse - 10x", which contains the original (non-duplicated) cells.
# 
# And remember we must include the `is_primary_data` column.

# In[1]:


import cellxgene_census

tabula_muris_dataset_id = "48b37086-25f7-4ecd-be66-f5bb378e3aea"
census = cellxgene_census.open_soma()

tabula_muris_obs = cellxgene_census.get_obs(
    census,
    "mus_musculus",
    value_filter=f"dataset_id == '{tabula_muris_dataset_id}'",
    column_names=["tissue", "is_primary_data"],
)


# Now let's take a look at counts for the unique combinations of values.
# 

# In[2]:


tabula_muris_obs.value_counts()


# You can see all cells across the tissues are labelled as `True` for `is_primary_data`.
# 
# But what if we select cells from the dataset that only contains cells from the liver: "Liver - A single-cell transcriptomic atlas characterizes ageing tissues in the mouse - 10x".
# 

# In[3]:


tabula_muris_liver_dataset_id = "6202a243-b713-4e12-9ced-c387f8483dea"

tabula_muris_liver_obs = cellxgene_census.get_obs(
    census,
    "mus_musculus",
    value_filter=f"dataset_id == '{tabula_muris_liver_dataset_id}'",
    column_names=["tissue", "is_primary_data"],
)


# And we take a look at counts for the unique combinations of values.

# In[4]:


tabula_muris_liver_obs.value_counts()


# You can see that:
# 
# 1. This dataset only contains cells from liver.
# 2. All cells are labelled as `False` for `is_primary_data`. **This is because the cells are marked as duplicate cells of the main Tabula Muris Senis dataset.**
# 
# ##  Filtering out duplicate cells
# 
# In some cases you may be interested in getting all cells for a specific biological context, for example *"all natural killer cells from blood of female cells with COVID-19"* but you need to be aware that there is a chance you end up with some duplicate cells.
# 
# We therefore recommend that you always look at `is_primary_data` and use that information based on your needs.
# 
# If you know *a priori* that you don't want duplicated cells this section shows you how to efficiently exclude them from your queries. 
# 
# ### Filtering out duplicate cells when reading the `obs` data frame.
# 
# Let's say you are interested in looking at the cell metadata of *"all natural killer cells from blood of female cells with COVID-19"* but you want to exclude duplicate cells, then you can use `value_filter` when reading the data frame to only include cells with `is_primary_data` as `True`.
# 
# Let's first read the cell metadata including **all** cells:

# In[5]:


nk_cells = cellxgene_census.get_obs(
    census,
    "mus_musculus",
    value_filter="cell_type == 'natural killer cell' "
    "and disease == 'COVID-19' "
    "and sex == 'female'"
    "and tissue_general == 'blood'",
)


# In[6]:


nk_cells.shape


# And now we repeat the query only using cells marked as `True` for `is_primary_data`.

# In[7]:


nk_cells_primary = cellxgene_census.get_obs(
    census,
    "mus_musculus",
    value_filter="cell_type == 'natural killer cell' "
    "and disease == 'COVID-19' "
    "and tissue_general == 'blood'"
    "and sex == 'female'"
    "and is_primary_data == True",
)


# In[8]:


nk_cells_primary.shape


# You can see a clear reduction in the number of cells.
# 
# ### Filtering out duplicate cells when creating an AnnData
# 
# You can also utilize `is_primary_data` on the `obs_value_filter` of `get_anndata`.
# 
# Let's repeat the process above. First querying by including **all** cells. To reduce the bandwidth and memory usage, let's just fetch data for one gene. 

# In[9]:


adata = cellxgene_census.get_anndata(
    census,
    organism="Homo sapiens",
    var_value_filter="feature_name == 'AQP5'",
    obs_value_filter="cell_type == 'natural killer cell' "
    "and disease == 'COVID-19' "
    "and sex == 'female'"
    "and tissue_general == 'blood'",
)


# In[10]:


len(adata.obs)


# And now we repeat the query only using cells marked as `True` for `is_primary_data`.

# In[11]:


adata_primary = cellxgene_census.get_anndata(
    census,
    organism="Homo sapiens",
    var_value_filter="feature_name == 'AQP5'",
    obs_value_filter="cell_type == 'natural killer cell' "
    "and disease == 'COVID-19' "
    "and sex == 'female' "
    "and tissue_general == 'blood'"
    "and is_primary_data == True",
)


# In[12]:


len(adata_primary.obs)


# In this case you can also observe a clear reduction in the number of cells.
# 
# #### Filtering out duplicate cells for out-of-core operations.
# 
# Finally we can utilize `is_primary_data` on the `value_filter` of `obs` of an "Axis Query" to perform out-of-core operations.
# 
# In this example we only include the version with duplicated cells removed.

# In[13]:


import tiledbsoma

human = census["census_data"]["homo_sapiens"]

# initialize lazy query
query = human.axis_query(
    measurement_name="RNA",
    obs_query=tiledbsoma.AxisQuery(
        value_filter="cell_type == 'natural killer cell' "
        "and disease == 'COVID-19' "
        "and tissue_general == 'blood' "
        "and sex == 'female' "
        "and is_primary_data == True"
    ),
)

# get iterator for X
iterator = query.X("raw").tables()

# iterate in chunks
for chunk in iterator:
    print(chunk)

    # since this is a demo we stop right away
    break

