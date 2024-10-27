#!/usr/bin/env python
# coding: utf-8

# # Summarizing cell and gene metadata
# 
# This notebook provides examples for basic axis metadata handling using Pandas. The Census stores `obs` (cell) and `var` (gene) metadata in `SOMADataFrame` objects via the [TileDB-SOMA API](https://github.com/single-cell-data/TileDB-SOMA) ([documentation](https://tiledbsoma.readthedocs.io/en/latest/)), which can be queried and read as a Pandas `DataFrame` using `TileDB-SOMA`. 
# 
# Note that Pandas `DataFrame` is an in-memory object, therefore queries should be small enough for results to fit in memory.
# 
# **Contents**
# 
# 1. Opening the Census
# 1. Summarizing cell metadata
#    1. Example: Summarize all cell types
#    1. Example: Summarize a subset of cell types, selected with a `value_filter`
# 1. Full Census metadata stats
# 
# ⚠️ Note that the Census RNA data includes duplicate cells present across multiple datasets. Duplicate cells can be filtered in or out using the cell metadata variable `is_primary_data` which is described in the [Census schema](https://github.com/chanzuckerberg/cellxgene-census/blob/main/docs/cellxgene_census_schema.md#repeated-data).
# 
# ## Opening the Census
# 
# The `cellxgene_census` python package contains a convenient API to open the latest version of the Census. If you open the Census, you should close it. `open_soma()` returns a context, so you can open/close it in several ways, like a Python file handle. The context manager is preferred, as it will automatically close upon an error raise.
# 
# You can learn more about the `cellxgene_census` methods by accessing their corresponding documentation via `help()`. For example `help(cellxgene_census.open_soma)`.

# In[1]:


import cellxgene_census

# Preferred: use a Python context manager
with cellxgene_census.open_soma() as census:
    ...

# or, directly open the census (don't forget to close it!)
census = cellxgene_census.open_soma()


# ## Summarizing cell metadata
# 
# Once the Census is open you can use its `TileDB-SOMA` methods as it is itself a `SOMACollection`. You can thus access the metadata `SOMADataFrame` objects encoding cell and gene metadata.
# 
# Tips:
# 
# - You can read an _entire_ `SOMADataFrame` into a Pandas `DataFrame` using `soma_df.read().concat().to_pandas()`, allowing the use of the standard Pandas API.
# - Queries will be much faster if you request only the DataFrame columns required for your analysis (e.g., `column_names=["cell_type_ontology_term_id"]`).
# - You can also further refine query results by using a `value_filter`, which will filter the census for matching records.
# 
# ### Example: Summarize all cell types
# 
# This example reads the cell metadata (`obs`) into a Pandas DataFrame, and summarizes in a variety of ways using Pandas API.

# In[2]:


# Read entire _obs_ into a pandas dataframe.
obs_df = cellxgene_census.get_obs(census, "homo_sapiens", column_names=["cell_type_ontology_term_id"])

# Use Pandas API to find all unique values in the `cell_type_ontology_term_id` column.
unique_cell_type_ontology_term_id = obs_df.cell_type_ontology_term_id.unique()

# Display only the first 10, as there are a LOT!
print(
    f"There are {len(unique_cell_type_ontology_term_id)} cell types in the Census! The first 10 are:",
    unique_cell_type_ontology_term_id[0:10].tolist(),
)

# Using Pandas API, count the instances of each cell type term and return the top 10.
top_10 = obs_df.cell_type_ontology_term_id.value_counts()[0:10]
print("\nThe top 10 cell types and their counts are:")
print(top_10)


# ### Example: Summarize a subset of cell types, selected with a `value_filter`
# 
# This example utilizes a SOMA "value filter" to read the subset of cells with `tissue_ontology_term_id` equal to `UBERON:0002048` (lung tissue), and summarizes the query result using Pandas.

# In[3]:


# Count cell_type occurrences for cells with tissue == 'lung'

# Read cell_type terms for cells which have a specific tissue term
LUNG_TISSUE = "UBERON:0002048"

obs_df = cellxgene_census.get_obs(
    census,
    "homo_sapiens",
    column_names=["cell_type_ontology_term_id"],
    value_filter=f"tissue_ontology_term_id == '{LUNG_TISSUE}'",
)

# Use Pandas API to find all unique values in the `cell_type_ontology_term_id` column.
unique_cell_type_ontology_term_id = obs_df.cell_type_ontology_term_id.unique()

print(
    f"There are {len(unique_cell_type_ontology_term_id)} cell types in the Census where tissue_ontology_term_id == {LUNG_TISSUE}! The first 10 are:",
    unique_cell_type_ontology_term_id[0:10].tolist(),
)

# Use Pandas API to count, and grab 10 most common
top_10 = obs_df.cell_type_ontology_term_id.value_counts()[0:10]
print(f"\nTop 10 cell types where tissue_ontology_term_id == {LUNG_TISSUE}")
print(top_10)


# You can also define much more complex value filters. For example:
# 
# * combine terms with `and` and `or`
# * use the `in` operator to query on multiple values

# In[4]:


# You can also do more complex queries, such as testing for inclusion in a list of values and "and" operations
VENTRICLES = ["UBERON:0002082", "UBERON:OOO2084", "UBERON:0002080"]

obs_df = cellxgene_census.get_obs(
    census,
    "homo_sapiens",
    column_names=["cell_type_ontology_term_id"],
    value_filter=f"tissue_ontology_term_id in {VENTRICLES} and is_primary_data == True",
)

# Use Pandas API to summarize
top_10 = obs_df.cell_type_ontology_term_id.value_counts()[0:10]
display(top_10)


# ## Full Census metadata stats
# 
# This example queries all organisms in the Census, and summarizes the diversity of various metadata lables.

# In[5]:


COLS_TO_QUERY = [
    "cell_type_ontology_term_id",
    "assay_ontology_term_id",
    "tissue_ontology_term_id",
]

obs_df = {
    name: cellxgene_census.get_obs(census, name, column_names=COLS_TO_QUERY) for name in census["census_data"].keys()
}

# Use Pandas API to summarize each organism
print(f"Complete census contains {sum(len(df) for df in obs_df.values())} cells.")
for organism, df in obs_df.items():
    print(organism)
    for col in COLS_TO_QUERY:
        print(f"\tUnique {col} values: {len(df[col].unique())}")


# Close the census

# In[6]:


census.close()

