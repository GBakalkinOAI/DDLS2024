#!/usr/bin/env python
# coding: utf-8

# # Exploring pre-calculated summary cell counts
# 
# This tutorial describes how to access pre-calculated summary cell counts. Each Census contains a top-level dataframe summarizing counts of various cell labels, this is the `census_summary_cell_counts` dataframe . You can read this into a Pandas DataFrame
# 
# **Contents**
# 
# 1. Fetching the `census_summary_cell_counts` dataframe.
# 2. Creating summary counts beyond pre-calculated values.
# 
# ⚠️ Note that the Census RNA data includes duplicate cells present across multiple datasets. Duplicate cells can be filtered in or out using the cell metadata variable `is_primary_data` which is described in the [Census schema](https://github.com/chanzuckerberg/cellxgene-census/blob/main/docs/cellxgene_census_schema.md#repeated-data).
# 
# ## Fetching the `census_summary_cell_counts` dataframe

# In[1]:


import cellxgene_census

census = cellxgene_census.open_soma()
census_summary_cell_counts = census["census_info"]["summary_cell_counts"].read().concat().to_pandas()

# Dropping the soma_joinid column as it isn't useful in this demo
census_summary_cell_counts = census_summary_cell_counts.drop(columns=["soma_joinid"])

census_summary_cell_counts


# ## Creating summary counts beyond pre-calculated values.
# 
# The dataframe above is precomputed from the experiments in the Census, providing a quick overview of the Census contents.
# 
# You can do similar group statistics using Pandas `groupby` functions. 
# 
# The code below reproduces the above counts using full `obs` dataframe in the `Homo_sapiens` experiment.
# 
# Keep in mind that the Census is very large, and any queries will return significant amount of data. You can manage that by narrowing the query request using `column_names` and `value_filter` in your query.

# In[2]:


human = census["census_data"]["homo_sapiens"]
obs_df = human.obs.read(column_names=["cell_type_ontology_term_id", "cell_type"]).concat().to_pandas()
obs_df.groupby(by=["cell_type_ontology_term_id", "cell_type"], as_index=False, observed=True).size()


# Close the census when complete. 

# In[3]:


census.close()

