#!/usr/bin/env python
# coding: utf-8

# # Querying data using the gget cellxgene module
# 
# *By Laura Luebbert, lauralubbert@gmail.com.*
# 
# [gget](https://github.com/pachterlab/gget) is a free, open-source command-line tool and Python package that enables efficient querying of genomic databases. gget consists of a collection of separate but interoperable modules, each designed to facilitate one type of database querying in a single line of code.
# 
# The [gget cellxgene](https://pachterlab.github.io/gget/en/cellxgene.html) module builds on the [CZ CELLxGENE Discover Census](https://chanzuckerberg.github.io/cellxgene-census/) to query data from [CZ CELLxGENE Discover](https://cellxgene.cziscience.com/). This notebook briefly introduces the [gget cellxgene](https://pachterlab.github.io/gget/en/cellxgene.html) module by providing one simple example for each supported query type.
# 
# If you use gget cellxgene in a publication, please [cite gget](https://pachterlab.github.io/gget/en/cite.html) in addition to [citing CZ CELLxGENE](https://cellxgene.cziscience.com/docs/08__Cite%20cellxgene%20in%20your%20publications).
# 
# You can also [open this notebook in Google Colab](https://colab.research.google.com/github/chanzuckerberg/cellxgene-census/blob/main/api/python/notebooks/api_demo/census_gget_demo.ipynb).
# 
# **Contents**  
# 
# 1. Install gget.
# 2. Fetch an [AnnData](https://anndata.readthedocs.io/en/latest/) object by selecting gene(s), tissue(s) and cell type(s).
# 3. Plot a dot plot similar to those shown on the  CZ CELLxGENE Discover [Gene Expression](https://cellxgene.cziscience.com/gene-expression).
# 4. Fetch only cell metadata (corresponds to AnnData.obs).
# 5. Use [gget cellxgene](https://pachterlab.github.io/gget/en/cellxgene.html) from the command line.

# ## Install gget and set up cellxgene module

# In[1]:


# The cellxgene module was added to gget in version 0.25.7
get_ipython().system('pip install -q gget >=0.25.7')

import gget

gget.setup("cellxgene")


# In[2]:


# Display all options of the cellxgene gget module
help(gget.cellxgene)


# ## Fetch an [AnnData](https://anndata.readthedocs.io/en/latest/) object by selecting gene(s), tissue(s) and cell type(s)
# You can use all of the options listed above to filter for data of interest. Here, we will demonstrate the module by fetching a small dataset containing only three genes and two lung cell types:

# In[3]:


# Fetch AnnData object based on specified genes, tissue and cell types
adata = gget.cellxgene(
    gene=["ACE2", "ABCA1", "SLC5A1"], tissue="lung", cell_type=["mucus secreting cell", "neuroendocrine cell"]
)


# Let's look at some of the features of the AnnData object we just fetched:

# In[4]:


adata


# A few thousand cells from CZ CELLxGENE Discover matched the filters specified above and their ACE2, ABCA1, and SLC5A1 expression matrix in lung mucus secreting and neuroendocrine cells was fetched. The `.var` and `.obs` layers contain additional information about each gene and cell, respectively:

# In[5]:


adata.var


# In[6]:


adata.obs


# ## Plot a dot plot similar to those shown on the  CZ CELLxGENE Discover [Gene Expression](https://cellxgene.cziscience.com/gene-expression)
# Using the data we just fetched, we can plot a dot plot using [scanpy](https://scanpy.readthedocs.io/en/stable/):

# In[7]:


import scanpy as sc

# retina increases the resolution of plots displayed in notebooks
get_ipython().run_line_magic('config', 'InlineBackend.figure_format="retina"')


# In[8]:


sc.pl.dotplot(adata, adata.var["feature_name"].values, groupby="cell_type", gene_symbols="feature_name")


# ## Fetch only cell metadata (corresponds to AnnData.obs)
# By setting `meta_only=True` and again filtering by the cell metadata attributes listed above, you can also fetch only the cell metadata:

# In[9]:


df = gget.cellxgene(
    meta_only=True,
    census_version="2023-05-15",  # Specify Census version for reproducibility over time
    gene="ENSMUSG00000015405",
    ensembl=True,  # Setting 'ensembl=True' here since the gene is passed as an Ensembl ID
    tissue="lung",
    species="mus_musculus",  # Let's switch up the species
)

df


# ## Use [gget cellxgene](https://pachterlab.github.io/gget/en/cellxgene.html) from the command line
# All gget modules support use from the command line. Note that the command line interface requires the `-o/--out` argument to specify a path to save the fetched data. Here are the command line versions of the queries demonstrated above:

# In[10]:


# # Fetch AnnData object based on specified genes, tissue and cell types
# !gget cellxgene --gene ACE2 ABCA1 SLC5A1 --tissue lung --cell_type 'mucus secreting cell' 'neuroendocrine cell' -o example_adata.h5ad


# In[11]:


# # Fetch only metadata
# !gget cellxgene --meta_only --gene ENSMUSG00000015405 --ensembl --tissue lung --species mus_musculus -o example_meta.csv

