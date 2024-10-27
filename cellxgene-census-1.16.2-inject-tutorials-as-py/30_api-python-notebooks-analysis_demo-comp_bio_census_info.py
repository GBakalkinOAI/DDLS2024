#!/usr/bin/env python
# coding: utf-8

# # Learning about the CZ CELLxGENE Census
# 
# This notebook showcases the Census contents and how to obtain high-level information about it. It covers the organization of data within the Census, what cell and gene metadata are available, and it provides simple demonstrations to summarize cell counts across cell metadata. 
# 
# **Contents**
# 
# - Opening the census
# - Census organization
# - Cell metadata
# - Gene metadata
# - Census summary content tables
# - Understanding Census contents beyond the summary tables
# 
# ⚠️ Note that the Census RNA data includes duplicate cells present across multiple datasets. Duplicate cells can be filtered in or out using the cell metadata variable `is_primary_data` which is described in the [Census schema](https://github.com/chanzuckerberg/cellxgene-census/blob/main/docs/cellxgene_census_schema.md#repeated-data).
# 
# ## Opening the Census
# 
# The `cellxgene_census` python package contains a convenient API to open the latest version of the Census. If you open the census, you should close it. `open_soma()` returns a context, so you can open/close it in several ways, like a Python file handle. The context manager is preferred, as it will automatically close upon an error raise.
# 

# In[1]:


import cellxgene_census

# Preferred: use a Python context manager
with cellxgene_census.open_soma() as census:
    ...

# or
census = cellxgene_census.open_soma()
...
census.close()


# You can learn more about the `cellxgene_census` methods by accessing their corresponding documentation via `help()`. For example `help(cellxgene_census.open_soma)`.
# 

# In[2]:


census = cellxgene_census.open_soma()


# ## Census organization
# 
# The [Census schema](https://chanzuckerberg.github.io/cellxgene-census/cellxgene_census_docsite_schema.html) defines the structure of the Census. In short, you can think of the Census as a structured collection of items that stores different pieces of information. All of these items and the parent collection are SOMA objects of various types and can all be accessed with the [TileDB-SOMA API](https://github.com/single-cell-data/TileDB-SOMA) ([documentation](https://tiledbsoma.readthedocs.io/en/latest/)).
# 
# 
# The `cellxgene_census` package contains some convenient wrappers of the `TileDB-SOMA` API. An example of this is the function we used to open the Census: `cellxgene_census.open_soma()`
# 
# ### Main Census components
# 
# With the command above you created `census`, which is a `SOMACollection`. It is analogous to a Python dictionary, and it has two items: `census_info` and `census_data`.
# 
# #### Census summary info
# 
# - `census["census_info"]` A collection of tables providing information of the census as a whole.
#   - `census["census_info"]["summary"]`: A data frame with high-level information of this Census, e.g. build date, total cell count, etc.
#   - `census["census_info"]["datasets"]`: A data frame with all datasets from [CELLxGENE Discover](https://cellxgene.cziscience.com/) used to create the Census.
#   - `census["census_info"]["summary_cell_counts"]`: A data frame with cell counts stratified by **relevant** cell metadata
# 
# #### Census data
# 
# Data for each organism is stored in independent `SOMAExperiment` objects which are a specialized form of a `SOMACollection`. Each of these store a data matrix (cell by genes), cell metadata, gene metadata, and some other useful components not covered in this notebook.
# 
# This is how the data is organized for one organism -- _Homo sapiens_:
# 
# - `census_obj["census_data"]["homo_sapiens"].obs`: Cell metadata
# - `census_obj["census_data"]["homo_sapiens"].ms["RNA"].X:` Data matrices, currently only raw counts exist `X["raw"]`
# - `census_obj["census_data"]["homo_sapiens"].ms["RNA"].var:` Gene Metadata

# ## Cell metadata
# 
# You can obtain all cell metadata variables by directly querying the columns of the corresponding `SOMADataFrame`.
# 
# All of these variables can be used for querying the Census in case you want to work with specific cells.
# 

# In[3]:


keys = list(census["census_data"]["homo_sapiens"].obs.keys())

keys


# All of these variables are defined in the [CELLxGENE dataset schema](https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/3.0.0/schema.md#obs-cell-metadata) except for the following:
# 
# - `soma_joinid`: a SOMA-defined value use for join operations.
# - `dataset_id`: the dataset id as encoded in `census["census-info"]["datasets"]`.
# - `tissue_general` and `tissue_general_ontology_term_id`: the high-level tissue mapping.
# 

# ## Gene metadata
# 
# Similarly, we can obtain all gene metadata variables by directly querying the columns of the corresponding `SOMADataFrame`.
# 
# These are the variables you can use for querying the Census in case there are specific genes you are interested in.
# 

# In[4]:


keys = list(census["census_data"]["homo_sapiens"].ms["RNA"].var.keys())

keys


# All of these variables are defined in the [CELLxGENE dataset schema](https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/3.0.0/schema.md#var-and-rawvar-gene-metadata) except for the following:
# 
# - `soma_joinid`: a SOMA-defined value use for join operations.
# - `feature_length`: the length in base pairs of the gene.
# 

# In[5]:


census_info = census["census_info"]["summary"].read().concat().to_pandas()

census_info


# ## Census summary content tables
# 
# You can take a quick look at the high-level Census information by looking at `census["census_info"]["summary"]`
# 
# Of special interest are the `label`-`value` combinations for :
# 
# - `total_cell_count` is the total number of cells in the Census.
# - `unique_cell_count` is the number of unique cells, as some cells may be present twice due to meta-analysis or consortia-like data.
# - `number_donors_homo_sapiens` and `number_donors_mus_musculus` are the number of individuals for human and mouse. These are not guaranteed to be unique as one individual ID may be present or identical in different datasets.
# 
# ### Cell counts by cell metadata
# 
# By looking at `census["summary_cell_counts"]` you can get a general idea of cell counts stratified by **some relevant** cell metadata. Not all cell metadata is included in this table, you can take a look at all cell and gene metadata available in the sections below "Cell metadata" and "Gene metadata".
# 
# The line below retrieves this table and casts it into a `pandas.DataFrame`.
# 

# In[6]:


census_counts = census["census_info"]["summary_cell_counts"].read().concat().to_pandas()

census_counts


# For each combination of `organism` and values for each `category` of cell metadata you can take a look at `total_cell_count` and `unique_cell_count` for the cell counts of that combination.
# 
# The values for each `category` are specified in `ontology_term_id` and `label`, which are the value's IDs and labels, respectively.
# 
# #### Example: cell metadata included in the summary counts table
# 
# To get all the available cell metadata in the summary counts table you can do the following. Remember this is not all the cell metadata available, as some variables were omitted in the creation of this table.
# 

# In[7]:


census_counts[["organism", "category"]].value_counts(sort=False)


# #### Example: cell counts for each sequencing assay in human data
# 
# To get the cell counts for each sequencing assay type in human data, you can perform the following `pandas.DataFrame` operations:
# 

# In[8]:


census_human_assays = census_counts.query("organism == 'Homo sapiens' & category == 'assay'")
census_human_assays.sort_values("total_cell_count", ascending=False)


# #### Example: number of microglial cells in the Census
# 
# If you have a specific term from any of the categories shown above you can directly find out the number of cells for that term.
# 

# In[9]:


census_counts.query("label == 'microglial cell'")


# ## Understanding Census contents beyond the summary tables
# 
# While using the pre-computed tables in `census["census_info"]` is an easy and quick way to understand the contents of the Census, it falls short if you want to learn more about certain slices of the Census.
# 
# For example, you may want to learn more about:
# 
# - What are the cell types available for human liver?
# - What are the total number of cells in all lung datasets stratified by sequencing technology?
# - What is the sex distribution of all cells from brain in mouse?
# - What are the diseases available for T cells?
# 
# All of these questions can be answered by directly querying the cell metadata as shown in the examples below.
# 
# ### Example: all cell types available in human
# 
# To exemplify the process of accessing and slicing cell metadata for summary stats, let's start with a trivial example and take a look at all human cell types available in the Census:
# 

# In[10]:


human_cell_types = (
    census["census_data"]["homo_sapiens"].obs.read(column_names=["cell_type", "is_primary_data"]).concat().to_pandas()
)
human_cell_types


# The number of rows is the total number of cells for humans. Now, if you wish to get the cell counts per cell type we can perform some `pandas` operations on this object.
# 
# In addition, we will only focus on cells that are marked with `is_primary_data=True` as this ensures we de-duplicate cells that appear more than once in CELLxGENE Discover.
# 

# In[11]:


human_cell_types = (
    census["census_data"]["homo_sapiens"]
    .obs.read(column_names=["cell_type"], value_filter="is_primary_data == True")
    .concat()
    .to_pandas()
)

human_cell_types = human_cell_types[["cell_type"]]
human_cell_types.shape


# This is the number of unique cells. Now let's look at the counts per cell type:
# 

# In[12]:


human_cell_type_counts = human_cell_types.value_counts()
human_cell_type_counts


# This shows you that the most abundant cell types are "glutamatergic neuron", "CD8-positive, alpha-beta T cell", and "CD4-positive, alpha-beta T cell".
# 
# Now let's take a look at the number of unique cell types:
# 

# In[13]:


human_cell_type_counts.shape


# That is the total number of different cell types for human.
# 
# All the information in this example can be quickly obtained from the summary table at `census["census-info"]["summary_cell_counts"]`.
# 
# The examples below are more complex and can only be achieved by accessing the cell metadata.
# 
# ### Example: cell types available in human liver
# 
# Similar to the example above, we can learn what cell types are available for a specific tissue, e.g. liver.
# 
# To achieve this goal we just need to limit our cell metadata to that tissue. We will use the information in the cell metadata variable `tissue_general`. This variable contains the high-level tissue label for all cells in the Census:
# 

# In[14]:


human_liver_cell_types = (
    census["census_data"]["homo_sapiens"]
    .obs.read(column_names=["cell_type"], value_filter="is_primary_data == True and tissue_general == 'liver'")
    .concat()
    .to_pandas()
)

human_liver_cell_types["cell_type"].value_counts()


# These are the cell types and their cell counts in the human liver.
# 
# ### Example: diseased T cells in human tissues
# 
# In this example we are going to get the counts for all diseased cells annotated as T cells. For the sake of the example we will focus on "CD8-positive, alpha-beta T cell" and "CD4-positive, alpha-beta T cell":
# 

# In[15]:


t_cells_list = ["CD8-positive, alpha-beta T cell", "CD4-positive, alpha-beta T cell"]

t_cells_diseased = (
    census["census_data"]["homo_sapiens"]
    .obs.read(
        column_names=["disease", "tissue_general"],
        value_filter=f"is_primary_data == True and cell_type in {t_cells_list} and disease != 'normal'",
    )
    .concat()
    .to_pandas()
)

t_cells_diseased = t_cells_diseased[["disease", "tissue_general"]].value_counts(sort=False)
t_cells_diseased


# These are the cell counts annotated with the indicated disease across human tissues for "CD8-positive, alpha-beta T cell" or "CD4-positive, alpha-beta T cell".
# 
# And, don't forget to close the census!
# 

# In[16]:


census.close()
del census

