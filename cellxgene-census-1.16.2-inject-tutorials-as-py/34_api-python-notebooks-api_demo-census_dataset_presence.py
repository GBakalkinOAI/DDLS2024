#!/usr/bin/env python
# coding: utf-8

# # Genes measured in each cell (dataset presence matrix)
# 
# The Census is a compilation of cells from multiple datasets that may differ by the sets of genes they measure. This notebook describes the way to identify the genes measured per dataset.
# 
# The presence matrix is a sparse boolean array, indicating which features (var) were present in each dataset.  The array has dimensions [n_datasets, n_var], and is stored in the SOMA Measurement `varp` collection. The first dimension is indexed by the `soma_joinid` in the `census_datasets` dataframe. The second is indexed by the `soma_joinid` in the `var` dataframe of the measurement.
# 
# As a reminder the `obs` data frame has a column `dataset_id` that can be used to link any cell in the Census to the presence matrix.
# 
# **Contents** 
# 
# 1. Opening the Census.
# 2. Fetching the IDs of the Census datasets.
# 3. Fetching the dataset presence matrix.
# 4. Identifying genes measured in a specific dataset.
# 5. Identifying datasets that measured specific genes.
# 6. Identifying all genes measured in a dataset.
# 
# ⚠️ Note that the Census RNA data includes duplicate cells present across multiple datasets. Duplicate cells can be filtered in or out using the cell metadata variable `is_primary_data` which is described in the [Census schema](https://github.com/chanzuckerberg/cellxgene-census/blob/main/docs/cellxgene_census_schema.md#repeated-data).
# 
# ## Opening the Census
# 
# The `cellxgene_census` python package contains a convenient API to open the latest version of the Census.

# In[1]:


import cellxgene_census

census = cellxgene_census.open_soma()


# ## Fetching the IDs of the Census datasets
# 
# Let's grab a table of all the datasets included in the Census and use this table in combination with the presence matrix below.

# In[2]:


# Grab the experiment containing human data, and the measurement therein with RNA
human = census["census_data"]["homo_sapiens"]
human_rna = human.ms["RNA"]

# The census-wide datasets
datasets_df = census["census_info"]["datasets"].read().concat().to_pandas()

datasets_df


# ## Fetching the dataset presence matrix
# 
# Now let's fetch the dataset presence matrix. 
# 
# For convenience, read the entire presence matrix (for Homo sapiens) into a SciPy array. There is a convenience API providing this capability, returning the matrix in a `scipy.sparse.array`.

# In[3]:


presence_matrix = cellxgene_census.get_presence_matrix(census, organism="Homo sapiens", measurement_name="RNA")

presence_matrix


# We also need the `var` dataframe, which is read into a Pandas DataFrame for convenient manipulation:

# In[4]:


var_df = human_rna.var.read().concat().to_pandas()

var_df


# ## Identifying genes measured in a specific dataset.
# 
# Now that we have the dataset table, the genes metadata table, and the dataset presence matrix, we can check if a gene or set of genes were measured in a specific dataset.
# 
# **Important:** the presence matrix is indexed by soma_joinid, and is *NOT* positionally indexed.  In other words:
# 
# * the first dimension of the presence matrix is the dataset's `soma_joinid`, as stored in the `census_datasets` dataframe.
# * the second dimension of the presence matrix is the feature's `soma_joinid`, as stored in the `var` dataframe.
# 
# Let's find out if the the gene `"ENSG00000286096"` was measured in the dataset with id `"97a17473-e2b1-4f31-a544-44a60773e2dd"`.
# 

# In[5]:


var_joinid = var_df.loc[var_df.feature_id == "ENSG00000286096"].soma_joinid
dataset_joinid = datasets_df.loc[datasets_df.dataset_id == "97a17473-e2b1-4f31-a544-44a60773e2dd"].soma_joinid
is_present = presence_matrix[dataset_joinid, var_joinid][0, 0]
print(f'Feature is {"present" if is_present else "not present"}.')


# ## Identifying datasets that measured specific genes
# 
# Similarly, we can determine the datasets that measured a specific gene or set of genes.

# In[6]:


# Grab the feature's soma_joinid from the var dataframe
var_joinid = var_df.loc[var_df.feature_id == "ENSG00000286096"].soma_joinid

# The presence matrix is indexed by the joinids of the dataset and var dataframes,
# so slice out the feature of interest by its joinid.
dataset_joinids = presence_matrix[:, var_joinid].tocoo().row

# From the datasets dataframe, slice out the datasets which have a joinid in the list
datasets_df.loc[datasets_df.soma_joinid.isin(dataset_joinids)]


# ## Identifying all genes measured in a dataset 
# 
# Finally, we can find the set of genes that were measured in the cells of a given dataset.

# In[7]:


# Slice the dataset(s) of interest, and get the joinid(s)
dataset_joinids = datasets_df.loc[datasets_df.collection_id == "17481d16-ee44-49e5-bcf0-28c0780d8c4a"].soma_joinid

# Slice the presence matrix by the first dimension, i.e., by dataset
var_joinids = presence_matrix[dataset_joinids, :].tocoo().col

# From the feature (var) dataframe, slice out features which have a joinid in the list.
var_df.loc[var_df.soma_joinid.isin(var_joinids)]

