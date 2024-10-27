#!/usr/bin/env python
# coding: utf-8

# # Querying and fetching the single-cell data and cell/gene metadata.
# 
# This tutorial showcases the easiest ways to query the expression data and cell/gene metadata from the Census, and load them into common in-memory Python objects, including `pandas.DataFrame` and `anndata.AnnData`.
# 
# **Contents**
# 
# 1. Opening the census.
# 2. Querying expression data.
# 3. Querying cell metadata (obs).
# 4. Querying gene metadata (var).
# 
# ⚠️ Note that the Census RNA data includes duplicate cells present across multiple datasets. Duplicate cells can be filtered in or out using the cell metadata variable `is_primary_data` which is described in the [Census schema](https://github.com/chanzuckerberg/cellxgene-census/blob/main/docs/cellxgene_census_schema.md#repeated-data).
# 
# ## Opening the census
# 
# The `cellxgene_census` python package contains a convenient API to open the latest version of the Census.

# In[1]:


import cellxgene_census

census = cellxgene_census.open_soma()


# You can learn more about the `cellxgene_census` methods by accessing their corresponding documentation via `help()`. For example `help(cellxgene_census.open_soma)`.

# ## Querying expression data
# 
# A convenient way to query and fetch expression data is to use the `get_anndata` method of the `cellxgene_census` API. This is a method that combines the column selection and value filtering we described above to obtain slices of the expression data based on metadata queries.
# 
# The method will return an `anndata.AnnData` object, it takes as an input a census object, the string for an organism, and for both cell and gene metadata we can specify filters and column selection as described above but with the following arguments:
# 
# - `obs_column_names` and `var_column_names` — a pair of arguments whose values are lists of strings indicating the columns to select for cell (`obs`) and gene (`var`) metadata respectively.
# - `obs_value_filter` —  python expression with selection conditions to fetch **cells** meeting a criteria. For full details see [tiledb.QueryCondition](https://tiledb-inc-tiledb.readthedocs-hosted.com/projects/tiledb-py/en/stable/python-api.html#query-condition).
# - `var_value_filter` —  python expression with selection conditions to fetch **genes** meeting a criteria. Details as above.  For full details see [tiledb.QueryCondition](https://tiledb-inc-tiledb.readthedocs-hosted.com/projects/tiledb-py/en/stable/python-api.html#query-condition).
# 
# 
# For example if we want to fetch the expression data for:
# 
# - Genes `"ENSG00000161798"` and `"ENSG00000188229"`.
# - All `"B cells"` of `"lung"` with `"COVID-19"` from non-duplicated cells.
# - With all gene metadata and adding `sex` cell metadata.

# In[2]:


adata = cellxgene_census.get_anndata(
    census=census,
    organism="Homo sapiens",
    var_value_filter="feature_id in ['ENSG00000161798', 'ENSG00000188229']",
    obs_value_filter="cell_type == 'B cell' and tissue_general == 'lung' and disease == 'COVID-19' and is_primary_data == True",
    obs_column_names=["sex"],
)


# And now we can take a look at the results.

# In[3]:


adata


# In[4]:


adata.obs


# In[5]:


adata.var


# For a full description of `get_anndata()` refer to `help(cellxgene_census.get_anndata)`
# 
# Don't forget to close the census!

# ## Querying cell metadata (obs)
# 
# The human gene metadata of the Census, for RNA assays, is located at `census["census_data"]["homo_sapiens"].obs`. This is a `SOMADataFrame` and as such it can be materialized as a `pandas.DataFrame` via the methods `read().concat().to_pandas()`. See also, the helper function `cellxgene_census.get_obs` which removes some boiler plate.
# 
# The mouse cell metadata is at `census["census_data"]["mus_musculus"].obs`.
# 
# For slicing the cell metadata there are two relevant arguments that can be passed through `read()`:
# 
# - `column_names` — list of strings indicating what metadata columns to fetch. 
# - `value_filter` — Python expression with selection conditions to fetch rows, it is similar to [pandas.DataFrame.query()](https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.query.html), for full details see [tiledb.QueryCondition](https://tiledb-inc-tiledb.readthedocs-hosted.com/projects/tiledb-py/en/stable/python-api.html#query-condition) shortly:
#    - Expressions are one or more comparisons
#    - Comparisons are one of `<column> <op> <value>` or `<column> <op> <column>`
#    - Expressions can combine comparisons using and, or, & or |
#    - op is one of < | > | <= | >= | == | != or in
# 
# To learn what metadata columns are available for fetching and filtering we can directly look at the keys of the cell metadata.

# In[6]:


keys = list(census["census_data"]["homo_sapiens"].obs.keys())

keys


# `soma_joinid` is a special `SOMADataFrame` column that is used for join operations. The definition for all other columns can be found at the [Census schema](https://github.com/chanzuckerberg/cellxgene-census/blob/main/docs/cell_census_schema.md#cell-metadata--census_objcensus_dataorganismobs--somadataframe).
# 
# All of these can be used to fetch specific columns or specific rows matching a condition. For the latter we need to know the values we are looking for _a priori_.
# 
# For example let's see what are the possible values available for `sex`. To this we can load all cell metadata but fetching only for the column `sex`. 

# In[7]:


sex_cell_metadata = cellxgene_census.get_obs(census, "homo_sapiens", column_names=["sex"])

sex_cell_metadata.drop_duplicates()


# As you can see there are only three different values for `sex`, that is `"male"`, `"female"` and `"unknown"`. 
# 
# With this information we can fetch all cell metatadata for a specific `sex` value, for example `"unknown"`.

# In[8]:


cell_metadata_all_unknown_sex = cellxgene_census.get_obs(census, "homo_sapiens", value_filter="sex == 'unknown'")

cell_metadata_all_unknown_sex


# You can use both `column_names` and `value_filter` to perform specific queries. For example let's fetch the `disease` columns for the `cell_type` `"B cell"` in the `tissue_general` `"lung"` and from non-duplicated cells. 

# In[9]:


cell_metadata_b_cell = cellxgene_census.get_obs(
    census,
    "homo_sapiens",
    value_filter="cell_type == 'B cell' and tissue_general == 'lung' and is_primary_data==True",
    column_names=["disease"],
)

cell_metadata_b_cell.value_counts()


# ## Querying gene metadata (var)
# 
# The human gene metadata of the Census is located at `census["census_data"]["homo_sapiens"].ms["RNA"].var`. Similarly to the cell metadata, it is a `SOMADataFrame` and thus we can also use its method `read()`.
# 
# The mouse gene metadata is at `census["census_data"]["mus_musculus"].ms["RNA"].var`.
# 
# Let's take a look at the metadata available for column selection and row filtering.

# In[10]:


keys = list(census["census_data"]["homo_sapiens"].ms["RNA"].var.keys())

keys


# With the exception of `soma_joinid` these columns are defined in the [Census schema](https://github.com/chanzuckerberg/cellxgene-census/blob/main/docs/cell_census_schema_0.1.0.md). Similarly to the cell metadata, we can use the same operations to learn and fetch gene metadata.
# 
# For example, to get the `feature_name` and `feature_length` of the genes `"ENSG00000161798"` and `"ENSG00000188229"` we can do the following.

# In[11]:


gene_metadata = cellxgene_census.get_var(
    census,
    "homo_sapiens",
    value_filter="feature_id in ['ENSG00000161798', 'ENSG00000188229']",
    column_names=["feature_name", "feature_length"],
)

gene_metadata


# In[12]:


census.close()

