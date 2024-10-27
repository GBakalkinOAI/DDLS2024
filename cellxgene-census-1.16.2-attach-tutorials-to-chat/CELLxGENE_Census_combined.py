
# Section: api_demo-census_datasets

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




# Section: api_demo-census_dataset_presence

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




# Section: analysis_demo-comp_bio_geneformer_prediction

#!/usr/bin/env python
# coding: utf-8

# # Geneformer for cell class prediction and data projection
# 
# This notebook provides examples to utilize the CELLxGENE collaboration fine-tuned Geneformer model with user data. For more information on the model please refer to the [Census model page](https://cellxgene.cziscience.com/census-models).
# 
# **IMPORTANT:** This tutorial requires cellxgene-census package version 1.9.1 or later.
# 
# **Contents**
# 
# 1. Requirements.
# 1. Preparing data and model.
# 1. Using the Geneformer fine-tuned model for **cell subclass inference**.
# 1. Using the Geneformer fine-tuned model for **data projection**.
# 
# > ⚠️ Note "cell subclass" is a high-level grouping of cell types as annotated in CELLxGENE Discover via the CL ontology see [https://cellxgene.cziscience.com/collections](https://cellxgene.cziscience.com/collections
# 
# > ⚠️ Note that the Census RNA data includes duplicate cells present across multiple datasets. Duplicate cells can be filtered in or out using the cell metadata variable `is_primary_data` which is described in the [Census schema](https://github.com/chanzuckerberg/cellxgene-census/blob/main/docs/cellxgene_census_schema.md#repeated-data).
# 
# ## Requirements
# 
# ### System requirements
# 
# To run this notebook the following are required:
# 
# - Unix system.
# - A system with one or more GPUs is highly recommended.
# - [AWS command-line interface](https://docs.aws.amazon.com/cli/latest/userguide/getting-started-install.html).
# - [Geneformer Python package](https://huggingface.co/ctheodoris/Geneformer) and its dependencies.
# - CELLxGENE Census package.
# 
# ### Downloading example data
# 
# Throughout the notebook the 10X PBMC 3K dataset will be used, you can download it via the following shell commands.
# 

# In[1]:


get_ipython().system('mkdir -p data')
get_ipython().system('wget -nv -O data/pbmc3k_filtered_gene_bc_matrices.tar.gz http://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz')
get_ipython().system('tar -xzf data/pbmc3k_filtered_gene_bc_matrices.tar.gz -C data/')


# ### Downloading the fine-tuned Geneformer model

# The model is currently hosted in S3, you can find out more deatails in the [Census model page](https://cellxgene.cziscience.com/census-models).
# 
# Additional information, including its S3 URI, is also included in the metadata of the corresponding embeddings inside Census. These metadata can be obtained as follows.

# In[2]:


import cellxgene_census
import cellxgene_census.experimental

census_version = "2023-12-15"
organism = "homo_sapiens"

census = cellxgene_census.open_soma(census_version=census_version)

geneformer_info = cellxgene_census.experimental.get_embedding_metadata_by_name(
    embedding_name="geneformer",
    organism=organism,
    census_version=census_version,
)


# In[3]:


geneformer_info["model_link"]


# And we can download it via the AWS CLI.

# In[4]:


get_ipython().system('aws s3 sync --no-sign-request  --no-progress --only-show-errors s3://cellxgene-contrib-public/models/geneformer/2023-12-15/homo_sapiens/fined-tuned-model/ ./fine_tuned_geneformer')


# ### Importing required packages

# Finally all the required packages are loaded.

# In[5]:


import warnings

warnings.filterwarnings("ignore")

import json
import os

import cellxgene_census
import datasets
import numpy as np
import scanpy as sc
from geneformer import (
    DataCollatorForCellClassification,
    EmbExtractor,
    TranscriptomeTokenizer,
)
from transformers import BertForSequenceClassification, Trainer


# ## Preparing data and model
# 
# ### Preparing single-cell data
# 
# Let's load the test data. In preparation to use with Geneformer we do the following: 
# 
# - Set the index as the ENSEMBL gene ID and stores it in the `obs` column `"ensembl_id"`
#   - e.g. `ENSG00000139618` (*without* a version number suffix)
# - Add read counts to the `obs` column `"n_counts"`
# - Add an ID column to be used for joining later in the  `obs` column `"joinid"`
# 
# Then we write the resulting H5AD file to disk.

# In[6]:


adata = sc.read_10x_mtx("data/filtered_gene_bc_matrices/hg19/", var_names="gene_ids")
adata.var["ensembl_id"] = adata.var.index
adata.obs["n_counts"] = adata.X.sum(axis=1)
adata.obs["joinid"] = list(range(adata.n_obs))

h5ad_dir = "./data/h5ad/"

if not os.path.exists(h5ad_dir):
    os.makedirs(h5ad_dir)

adata.write(h5ad_dir + "pbmcs.h5ad")


# Now we can tokenize the test data using Geneformer's tokenizer, while keeping track of `"joinid"` for future joining.

# In[7]:


token_dir = "data/tokenized_data/"

if not os.path.exists(token_dir):
    os.makedirs(token_dir)

tokenizer = TranscriptomeTokenizer(custom_attr_name_dict={"joinid": "joinid"})
tokenizer.tokenize_data(
    data_directory=h5ad_dir,
    output_directory=token_dir,
    output_prefix="pbmc",
    file_format="h5ad",
)


# ### Preparing data from model
# 
# Then let's fetch the mapping dictionary between Geneformer IDs and the associated cell subclass labels. This information is stored along the fine-tuned model.

# In[8]:


model_dir = "./fine_tuned_geneformer/"
label_mapping_dict_file = os.path.join(model_dir, "label_to_cell_subclass.json")

with open(label_mapping_dict_file) as fp:
    label_mapping_dict = json.load(fp)


# This dictionary contains all the possible cell labels available for the model, and the predictions on the section below will use these labels.

# In[9]:


label_mapping_dict


# ## Using the Geneformer fine-tuned model for cell subclass inference
# 
# ### Loading tokenized data

# Let's load the tokenized test data.

# In[10]:


dataset = datasets.load_from_disk(token_dir + "pbmc.dataset")
dataset


# We add a dummy cell metadata column `"label"` needed for Geneformer to make predictions.

# In[11]:


dataset
dataset = dataset.add_column("label", [0] * len(dataset))


# ### Performing inference of cell subclass

# Now we can load the model and run the inference workflow.
# 
# > ⚠️ Note, this step will be slow with CPUs, a machine with one GPU is recommended

# In[12]:


# reload pretrained model
model = BertForSequenceClassification.from_pretrained(model_dir)
# create the trainer
trainer = Trainer(model=model, data_collator=DataCollatorForCellClassification())
# use trainer
predictions = trainer.predict(dataset)


# And finally we select the most likely cell class based on the probability vector from the predictions of each cell in our test data.

# In[13]:


predicted_label_ids = np.argmax(predictions.predictions, axis=1)
predicted_logits = [predictions.predictions[i][predicted_label_ids[i]] for i in range(len(predicted_label_ids))]
predicted_labels = [label_mapping_dict[str(i)] for i in predicted_label_ids]


# ### Inspecting inference results

# Then we add the prediction back to our loaded AnnData test dataset.

# In[14]:


adata.obs["predicted_cell_subclass"] = predicted_labels
adata.obs["predicted_cell_subclass_probability"] = np.exp(predicted_logits) / (1 + np.exp(predicted_logits))


# And it's ready for inspecting the predictions. Let's visualize the predictions on the UMAP space, the following is a basic processing workflow to derive a UMAP representation, of the data.

# In[15]:


sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata = adata[:, adata.var.highly_variable]
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver="arpack")
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)


# Let's also add the original cell type annotations as obtained in [Scapy's annotation tutorial](https://scanpy-tutorials.readthedocs.io/en/latest/pbmc3k.html) of the same data. 

# In[16]:


sc.tl.leiden(adata)
original_cell_types = [
    "CD4-positive, alpha-beta T cell (1)",
    "CD4-positive, alpha-beta T cell (2)",
    "CD14-positive, monocyte",
    "B cell (1)",
    "CD8-positive, alpha-beta T cell",
    "FCGR3A-positive, monocyte",
    "natural killer cell",
    "dendritic cell",
    "megakaryocyte",
    "B cell (2)",
]
adata.rename_categories("leiden", original_cell_types)


# These are the original annotations.

# In[17]:


sc.pl.umap(adata, color="leiden", title="Original Annotations")


# And these are the predicted annotations.

# In[18]:


sc.pl.umap(
    adata,
    color=["predicted_cell_subclass_probability", "predicted_cell_subclass"],
    title="Predicted Geneformer Annotations",
)


# ## Using the Geneformer fine-tuned model for data projection
# 
# ### Generating Geneformer embeddings for 10X PBMC 3K data
# 
# To project new data, for example the 10X PBMC 3K data, into the Census embedding space from Geneformer's fine-tune model, we can use `EmbExtractor` from the [Geneformer](https://huggingface.co/ctheodoris/Geneformer) package as follows.
# 
# We first need to get the number of categories (cell subclasses) present in the model.

# In[19]:


n_classes = len(label_mapping_dict)


# Then we can run the `EmbExtractor`, which randomize the cells during the process and thus we keep track of `"joinid"`.
# 
# > ⚠️ Note, this step will be slow with CPUs, a machine with one GPU is recommended

# In[20]:


output_dir = "data/geneformer_embeddings"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

embex = EmbExtractor(
    model_type="CellClassifier",
    num_classes=n_classes,
    max_ncells=None,
    emb_label=["joinid"],
    emb_layer=0,
    forward_batch_size=30,
    nproc=8,
)

embs = embex.extract_embs(
    model_directory=model_dir,
    input_data_file=token_dir + "pbmc.dataset",
    output_directory=output_dir,
    output_prefix="emb",
)


# Then we simply re-order the embeddings based on `"joinid"` and then merge them to the original AnnData

# In[21]:


embs = embs.sort_values("joinid")
adata.obsm["geneformer"] = embs.drop(columns="joinid").to_numpy()


# Let's take a look at these Geneformer embeddings in a UMAP representation

# In[22]:


sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40, use_rep="geneformer")
sc.tl.umap(adata)


# In[23]:


sc.pl.umap(adata, color="predicted_cell_subclass", title="10X PBMC 3K in Geneformer")


# ### Joining Geneformer embeddings from 10X PBMC 3K data with other Census datasets
# 
# There are multiple datasets in Census from PBMCs, and all human Census data has pre-calculated Geneformer embeddings, so now we can join the embeddings we generated above from the 10X PBMC 3K dataset with Census data.
# 
# Let's grab a few PBMC datasets from Census and request the Geneformer embeddings.

# In[24]:


census = cellxgene_census.open_soma(census_version="2023-12-15")

# Some PBMC data from these collections
# 1. https://cellxgene.cziscience.com/collections/c697eaaf-a3be-4251-b036-5f9052179e70
# 2. https://cellxgene.cziscience.com/collections/f2a488bf-782f-4c20-a8e5-cb34d48c1f7e

dataset_ids = [
    "fa8605cf-f27e-44af-ac2a-476bee4410d3",
    "3c75a463-6a87-4132-83a8-c3002624394d",
]

adata_census = cellxgene_census.get_anndata(
    census=census,
    measurement_name="RNA",
    organism="Homo sapiens",
    obs_value_filter=f"dataset_id in {dataset_ids}",
    obs_embeddings=["geneformer"],
)


# To simplify let's select the genes that are also present in the 10X PBMC 3K dataset.

# In[25]:


adata_census.var_names = adata_census.var["feature_id"]
shared_genes = list(set(adata.var_names) & set(adata_census.var_names))
adata_census = adata_census[:, shared_genes]


# And take a subset of these cells, let's take 3K cells to match the size of the test data. 

# In[26]:


index_subset = np.random.choice(adata_census.n_obs, size=3000, replace=False)
adata_census = adata_census[index_subset, :]


# Now we can join these Census data to the 10X PBMC 3K data

# In[27]:


adata_census.obs["dataset"] = "Census - " + adata_census.obs["dataset_id"].astype(str)
adata.obs["dataset"] = "10X PBMC 3K"
adata.obs["cell_type"] = "Predicted - " + adata.obs["predicted_cell_subclass"].astype(str)

adata_joined = sc.concat([adata, adata_census], join="outer", label="batch")


# Let's now inspect all of the cells in the UMAP space.

# In[28]:


sc.pp.neighbors(adata_joined, n_neighbors=10, n_pcs=40, use_rep="geneformer")
sc.tl.umap(adata_joined)


# In[29]:


sc.pl.umap(adata_joined, color="dataset")


# In[30]:


sc.pl.umap(adata_joined, color="cell_type")




# Section: analysis_demo-comp_bio_summarize_axis_query

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




# Section: experimental-highly_variable_genes

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




# Section: api_demo-census_and_cell_guide_example

#!/usr/bin/env python
# coding: utf-8

# # Query Census utilizing cell metadata ontologies 
# 
# This notebook demonstrates how to utilize the [CELLxGENE Ontology Guide](https://github.com/chanzuckerberg/cellxgene-ontology-guide/) API to leverage the ontological structure of Census cell metadata to perform flexible and succinct queries. For example, how to get all Census single cell for all T cells and their descendants.
# 
# The notebook focuses specifically on cell types, but the same principles can be used to any of the other metadata fields that utilize ontologies (e.g. tissue and developmental stage).
# 
# **IMPORTANT:** This tutorial requires [cellxgene-ontology-guide package](https://chanzuckerberg.github.io/cellxgene-ontology-guide/cellxgene_ontology_guide.html) version 1.0.0 or later.
# 
# **Contents**
# 
# 1. Obtain all descendant terms of a cell type
# 1. Query Census for all descendant terms of a cell type
# 
# ⚠️ Note that the Census RNA data includes duplicate cells present across multiple datasets. Duplicate cells can be filtered in or out using the cell metadata variable `is_primary_data` which is described in the [Census schema](https://github.com/chanzuckerberg/cellxgene-census/blob/main/docs/cellxgene_census_schema.md#repeated-data).
# 
# 
# ## Obtain all descendant terms of a cell type
# 
# The [CELLxGENE Ontology Guide](https://github.com/chanzuckerberg/cellxgene-ontology-guide/) provides a programatic interface to access and and parse the ontologies used by CELLxGENE for its cell and gene metadata. You can perform basic ontological operations, for example obtaining all ascendants or descendants of a cell type according to the Cell Ontology (CL).
# 
# To learn more about all the capabilities of the CELLxGENE Ontology Guide, please visit its [API documentation](https://chanzuckerberg.github.io/cellxgene-ontology-guide/cellxgene_ontology_guide.html).
# 
# In this notebook we will showcase an example to obtain all descendants of muscle cells. Let's first load the API, open the Census and fetch its **dataset** schema version.

# In[1]:


import warnings

import cellxgene_census
import scanpy as sc
from cellxgene_ontology_guide.ontology_parser import OntologyParser

warnings.filterwarnings("ignore")


with cellxgene_census.open_soma(census_version="latest") as census:
    census = cellxgene_census.open_soma(census_version="latest")
    census_info = census["census_info"]["summary"].read().concat().to_pandas()

census_info


# In[2]:


schema_version = census_info.iloc[2, 2]


# Now we can instantiate an `OntologyParser` to access the appropriate ontology versions according the dataset schema version.

# In[3]:


ontology_parser = OntologyParser(schema_version=schema_version)


# The `OntologyParser` class has methods to parse any available ontology in CELLxGENE. 
# 
# Let's identify all the descendants of muscle cells.

# In[4]:


muscle_cell = "CL:0000187"
mappings = ontology_parser.map_term_descendants([muscle_cell])

# top 5
mappings[muscle_cell][:5]


# ## Query Census for all descendant terms of a cell type
# 
# Now that we have all descendant terms of muscle cells, we can utilize those directly in Census to fetch all single-cell data for the all-encompassing set of muscle cells and its descendants.
# 
# To keep it small we will limit the query to only muscle cells of eye and esophagus.

# In[5]:


tissues = ["esophagus", "eye"]

with cellxgene_census.open_soma(census_version="2023-12-15") as census:
    census = cellxgene_census.open_soma(census_version="2023-12-15")
    value_filter = f"cell_type_ontology_term_id in {mappings[muscle_cell]} and tissue_general in {tissues} and is_primary_data == True"

    adata = cellxgene_census.get_anndata(
        census,
        organism="homo_sapiens",
        measurement_name="RNA",
        obs_value_filter=value_filter,
        obs_embeddings=["scvi"],
    )


# In[6]:


adata


# In[7]:


sc.pp.neighbors(adata, use_rep="scvi")
sc.tl.umap(adata)
sc.pl.umap(adata, color=["tissue_general", "cell_type"])




# Section: analysis_demo-comp_bio_scvi_model_use

#!/usr/bin/env python
# coding: utf-8

# # scVI for cell type prediction and data projection
# 
# This notebook provides examples to utilize the pretrained scVI model with user data. For more information on the model please refer to the [Census model page](https://cellxgene.cziscience.com/census-models).
# 
# **IMPORTANT:** This tutorial requires cellxgene-census package version 1.9.1 or later.
# 
# **Contents**
# 
# 1. Requirements.
# 1. Preparing data and model.
# 1. Using the scVI pretrained model for **data projection**.
# 1. Using the scVI pretrained model for **cell type inference**.
# 
# 
# ⚠️ Note that the Census RNA data includes duplicate cells present across multiple datasets. Duplicate cells can be filtered in or out using the cell metadata variable `is_primary_data` which is described in the [Census schema](https://github.com/chanzuckerberg/cellxgene-census/blob/main/docs/cellxgene_census_schema.md#repeated-data).
# 
# ## Requirements
# 
# ### System requirements
# 
# To run this notebook the following are required:
# 
# - Unix system.
# - A system with one or more GPUs is highly recommended.
# - [AWS command-line interface](https://docs.aws.amazon.com/cli/latest/userguide/getting-started-install.html).
# - [scvi-tools](https://github.com/scverse/scvi-tools) and its dependencies.
# - CELLxGENE Census package
# 
# ### Downloading example data
# 
# Throughout the notebook the 10X PBMC 3K dataset will be used, you can download it via the following shell commands.
# 

# In[1]:


get_ipython().system('mkdir -p data')
get_ipython().system('wget -nv -O data/pbmc3k_filtered_gene_bc_matrices.tar.gz http://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz')
get_ipython().system('tar -xzf data/pbmc3k_filtered_gene_bc_matrices.tar.gz -C data/')


# ### Downloading the trained scVI model
# 
# The model is currently hosted in S3, you can find out more deatails in the [Census model page](https://cellxgene.cziscience.com/census-models).
# 
# Additional information, including its S3 URI, is also included in the metadata of the corresponding embeddings inside Census. These metadata can be obtained as follows.

# In[2]:


import cellxgene_census
import cellxgene_census.experimental

census_version = "2023-12-15"
organism = "homo_sapiens"

census = cellxgene_census.open_soma(census_version=census_version)

scvi_info = cellxgene_census.experimental.get_embedding_metadata_by_name(
    embedding_name="scvi",
    organism=organism,
    census_version=census_version,
)


# In[3]:


scvi_info["model_link"]


# In[4]:


get_ipython().system('aws s3 cp --no-sign-request --no-progress --only-show-errors s3://cellxgene-contrib-public/models/scvi/2024-02-12/homo_sapiens/model.pt 2024-02-12-scvi-homo-sapiens/scvi.model/')


# ## Using the scVI pretrained model for **data projection**
# Import all the required packages for this demonstration

# In[5]:


import warnings

warnings.filterwarnings("ignore")

import anndata
import cellxgene_census
import numpy as np
import scanpy as sc
import scvi
from sklearn.ensemble import RandomForestClassifier


# Load the example query dataset (the 10X pbmc3k data).

# In[6]:


adata = sc.read_10x_mtx("data/filtered_gene_bc_matrices/hg19/", var_names="gene_ids")
adata.var["ensembl_id"] = adata.var.index
adata.obs["n_counts"] = adata.X.sum(axis=1)
adata.obs["joinid"] = list(range(adata.n_obs))
# initialize the batch to be unassigned. This could be any dummy value.
adata.obs["batch"] = "unassigned"


# Load the scVI model and prepare the query data

# In[7]:


folder = "2024-02-12-scvi-homo-sapiens"

model_filename = f"{folder}/scvi.model"
scvi.model.SCVI.prepare_query_anndata(adata, model_filename)


# Load the query data into the model, set "is_trained" to True to trick the model into thinking it was already trained, and do a forward pass through the model to get the latent reprsentation of the query data.

# In[8]:


vae_q = scvi.model.SCVI.load_query_data(
    adata,
    model_filename,
)

# This allows for a simple forward pass
vae_q.is_trained = True
latent = vae_q.get_latent_representation()
adata.obsm["scvi"] = latent

# filter out missing features
adata = adata[:, adata.var["gene_symbols"].notnull().values].copy()
adata.var.set_index("gene_symbols", inplace=True)


# Run UMAP

# In[9]:


sc.pp.neighbors(adata, n_neighbors=15, use_rep="scvi")
sc.tl.umap(adata)


# Run leiden clustering

# In[10]:


sc.tl.leiden(adata)


# Normalize and log-transform the expression data

# In[11]:


sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)


# Using the marker genes from the [Scanpy pbmc3k vignette](https://scanpy-tutorials.readthedocs.io/en/latest/pbmc3k.html), we can map our leiden clusters to the corresponding cell type labels used in the tutorial. Our Leiden clustering does not match up perfectly so we need to visualize the marker genes to appropriately map the clusters to the original cell type annotation.
# ![image.png](attachment:58d99df3-de02-4bea-8670-7b2d95a6f8c0.png)

# In[12]:


markers_row1 = ["IL7R", "CD14", "LYZ", "MS4A1", "CD8A", "GNLY"]
markers_row2 = ["NKG7", "FCGR3A", "MS4A7", "FCER1A", "CST3", "PPBP"]

with warnings.catch_warnings():
    warnings.filterwarnings("ignore")

    sc.pl.violin(adata, markers_row1, groupby="leiden")
    sc.pl.violin(adata, markers_row2, groupby="leiden")


# Based on the expression of the provided marker genes, we can map the following Leiden clusters to these cell type labels:
# 
#  - 2,3,4,5,10 = CD4 T cells
#  - 0 = CD14+ monocytes
#  - 1 = B cells
#  - 6,9 = CD8 T cells
#  - 8 = NK cells
#  - 7 = FCGR3A+ Monocytes
#  - 11 = dendritic cells
#  - 12 = megakaryocytes

# In[13]:


original_cell_types = [
    "CD14+ monocytes",
    "B cells",
    "CD4 T cells",
    "CD4 T cells",
    "CD4 T cells",
    "CD4 T cells",
    "CD8 T cells",
    "FCGR3A+ Monocytes",
    "NK cells",
    "CD8 T cells",
    "CD4 T cells",
    "dendritic cells",
    "megakaryocytes",
]
label_mapping = dict(zip(range(len(original_cell_types)), original_cell_types))
adata.obs["original_cell_type"] = adata.obs["leiden"].apply(lambda x: label_mapping[int(x)])


# In[14]:


sc.pl.umap(adata, color=["original_cell_type"])


# Display the scatter plot

# ## Using the scVI pretrained model for **cell cell type inference**.

# Fetch the reference scVI embeddings corresponding to some example PBMC data from Census

# In[15]:


census = cellxgene_census.open_soma(census_version="2023-12-15")

# Some PBMC data from these collections
# 1. https://cellxgene.cziscience.com/collections/c697eaaf-a3be-4251-b036-5f9052179e70
# 2. https://cellxgene.cziscience.com/collections/f2a488bf-782f-4c20-a8e5-cb34d48c1f7e
dataset_ids = [
    "fa8605cf-f27e-44af-ac2a-476bee4410d3",
    "3c75a463-6a87-4132-83a8-c3002624394d",
]

adata_census = cellxgene_census.get_anndata(
    census=census,
    measurement_name="RNA",
    organism="Homo sapiens",
    obs_value_filter=f"dataset_id in {dataset_ids}",
    obs_embeddings=["scvi"],
)
adata_census.var.set_index("feature_id", inplace=True)


# Let's run UMAP on a subset of the reference combined with the query dataset and plot the UMAP, coloring by dataset ID.

# In[16]:


adata.obs["dataset_id"] = "QUERY"
# Subset the reference dataset to have a similar number of cells to the query dataset
index_subset = np.random.choice(adata_census.n_obs, size=3000, replace=False)
adata_census_subset = adata_census[index_subset, :]

adata_combined = anndata.concat([adata_census_subset, adata])
sc.pp.neighbors(adata_combined, n_neighbors=15, use_rep="scvi", metric="correlation")
sc.tl.umap(adata_combined)
sc.pl.umap(adata_combined, color=["dataset_id"])


# Fit a Random Forest Classifier on the reference scVI embedding fetched from Census and use it to predict cell type labels on the projected scVI embedding for the query dataset.

# In[17]:


rfc = RandomForestClassifier()
rfc.fit(adata_census.obsm["scvi"], adata_census.obs["cell_type"].values)
adata.obs["predicted_cell_type"] = rfc.predict(adata.obsm["scvi"])

# let's get confidence scores
probabilities = rfc.predict_proba(adata.obsm["scvi"])
confidence = np.zeros(adata.n_obs)
for i in range(adata.n_obs):
    confidence[i] = probabilities[i][rfc.classes_ == adata.obs["predicted_cell_type"][i]]


# In[18]:


probabilities[i]


# In[19]:


# let's get confidence scores
probabilities = rfc.predict_proba(adata.obsm["scvi"])
confidence = np.zeros(adata.n_obs)
for i in range(adata.n_obs):
    confidence[i] = probabilities[i][rfc.classes_ == adata.obs["predicted_cell_type"][i]]

adata.obs["predicted_cell_type_probability"] = confidence


# Plot the results and compare the annotations

# In[20]:


sc.pl.umap(adata, color="original_cell_type")


# In[21]:


sc.pl.umap(adata, color=["predicted_cell_type_probability", "predicted_cell_type"])


# Let's look at the predicted cell type annotations on the combined query and reference datasets

# In[22]:


adata_combined.obs["cell_type"] = (
    adata_census_subset.obs["cell_type"].tolist() + adata.obs["predicted_cell_type"].tolist()
)
sc.pl.umap(adata_combined, color=["dataset_id", "cell_type"])




# Section: api_demo-census_duplicated_cells

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




# Section: experimental-mean_variance

#!/usr/bin/env python
# coding: utf-8

# # Out-of-core (incremental) mean and variance calculation
# 
# This tutorial describes use of the cellxgene_census.experimental.pp API for calculating out-of-core mean and variance in the Census. The variance calculation is performed using [Welford's online algorithm](https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance).
# 
# **Contents**
# 
# 1. The mean and variance API.
# 2. Example: calculate mean and variance for a slice of the Census.
# 
# ## The mean and variance API
# 
# `mean_variance()` calculates the mean and the variance for an `ExperimentAxisQuery`. The following additional arguments are supported:
# 
# - `layer`: the X layer used for the calculation, defaults to `raw`
# - `axis`: the axis along which the calculation is performed. Specify 0 for the `var` axis and 1 for the `obs` axis
# - `calculate_mean`: if False, do not include the mean in the result
# - `calculate_variance`: if False, do not compute the variance.
# - `ddof`: _Delta Degrees of Freedom_: the divisor used in the calculation for variance is N - ddof, where N represents the number of elements.
# 
# 

# In[1]:


# Import packages
import cellxgene_census
import pandas as pd
import tiledbsoma as soma
from cellxgene_census.experimental.pp import mean_variance


# ## Example: calculate mean and variance for a slice of the Census

# As an example, we'll calculate the mean and variance along the `obs` axis for a subset of cells from the Mouse census.
# 
# The return value will be a Pandas dataframe indexed by `soma_joinid` (in this case, it will be relative to `obs`) and will contain the `mean` and `variance` columns.

# In[2]:


experiment_name = "mus_musculus"
obs_value_filter = 'is_primary_data == True and tissue_general == "skin of body"'

with cellxgene_census.open_soma(census_version="stable") as census:
    with census["census_data"][experiment_name].axis_query(
        measurement_name="RNA", obs_query=soma.AxisQuery(value_filter=obs_value_filter)
    ) as query:
        mv_df = mean_variance(
            query,
            axis=1,
            calculate_mean=True,
            calculate_variance=True,
        )

        obs_df = query.obs().concat().to_pandas()

mv_df


# We can now concatenate the resulting dataframe to `obs`:

# In[3]:


combined_df = pd.concat([obs_df.set_index("soma_joinid"), mv_df], axis=1)
combined_df




# Section: analysis_demo-comp_bio_data_integration_scvi

#!/usr/bin/env python
# coding: utf-8

# # Integrating multi-dataset slices of data
# 
# The Census contains data from multiple studies providing an opportunity to perform inter-dataset analysis. To this end integration of data has to be performed first to account for batch effects.
# 
# This notebook provides a demonstration for integrating two Census datasets using [scvi-tools](https://docs.scvi-tools.org/en/stable/index.html). **The goal is not to provide an exhaustive guide on proper integration, but to showcase what information in the Census can inform data integration.**
# 
# **Contents**
# 
# 1. Finding and fetching data from mouse liver (10X Genomics and Smart-Seq2).
# 1. Gene-length normalization of Smart-Seq2 data.
# 1. Integration with `scvi-tools`.
#    1. Inspecting data prior to integration.
#    1. Integration with batch defined as `dataset_id`.
#    1. Integration with batch defined as `dataset_id` + `donor_id`.
#    1. Integration with batch defined as `dataset_id` + `donor_id` + `assay_ontology_term_id` + `suspension_type`.
#    
# ⚠️ Note that the Census RNA data includes duplicate cells present across multiple datasets. Duplicate cells can be filtered in or out using the cell metadata variable `is_primary_data` which is described in the [Census schema](https://github.com/chanzuckerberg/cellxgene-census/blob/main/docs/cellxgene_census_schema.md#repeated-data).
#  For this notebook we will focus on individual datasets, therefore we can ignore this variable.
# 
# ## Finding and fetching data from mouse liver (10X Genomics and Smart-Seq2)
# 
# Let's load all modules needed for this notebook.

# In[1]:


import cellxgene_census
import numpy as np
import scanpy as sc
import scvi
from scipy.sparse import csr_matrix


# Now we can open the Census, if you are not familiar with the basics of the Census API you should take a look at the notebook [Learning about the CELLxGENE Census](https://cellxgene-census.readthedocs.io/en/latest/notebooks/analysis_demo/comp_bio_census_info.html).

# In[2]:


census = cellxgene_census.open_soma(census_version="latest")


# In this notebook we will use Tabula Muris Senis data from the liver as it contains cells from both 10X Genomics and Smart-Seq2 technologies.
# 
# Let's query the `datasets` table of the Census by filtering on `collection_name` for "Tabula Muris Senis" and `dataset_title` for "liver". 

# In[3]:


census_datasets = (
    census["census_info"]["datasets"].read(value_filter="collection_name == 'Tabula Muris Senis'").concat().to_pandas()
)
tabula_liver = census_datasets["dataset_title"].str.contains("liver", case=False)
census_datasets.loc[tabula_liver,]


# Now we can use the values from `dataset_id` to query and load an AnnData object with all the cells from those datasets. 

# In[4]:


tabula_muris_liver_ids = [
    "4546e757-34d0-4d17-be06-538318925fcd",
    "6202a243-b713-4e12-9ced-c387f8483dea",
]

adata = cellxgene_census.get_anndata(
    census,
    organism="Mus musculus",
    obs_value_filter=f"dataset_id in {tabula_muris_liver_ids}",
)


# Close the census

# In[5]:


census.close()
del census


# We can check the cell counts for both 10X Genomics and Smart-Seq2 data by looking at `obs["assay"]`.

# In[6]:


adata.obs.assay.value_counts()


# ##  Gene-length normalization of Smart-Seq2 data.
# 
# Smart-seq2 read counts have to be normalized by gene length. For full details on gene-length normalization take a look at the notebook [Normalizing full-length gene sequencing data from the Census](https://cellxgene-census.readthedocs.io/en/latest/notebooks/analysis_demo/comp_bio_normalizing_full_gene_sequencing.html).
# 
# Let's first get the gene lengths from `var.feature_length`.

# In[7]:


smart_seq_gene_lengths = adata.var[["feature_length"]].to_numpy()


# Now we create a copy of a slice of the expression matrix only containing Smart-Seq2 data.

# In[8]:


smart_seq_index = np.where(adata.obs.assay == "Smart-seq2")[0]
smart_seq_index
smart_seq_X = adata.X[smart_seq_index, :].copy()


# We proceed to normalize it using the gene lengths.

# In[9]:


smart_seq_X = csr_matrix((smart_seq_X.T / smart_seq_gene_lengths).T)
smart_seq_X = smart_seq_X.ceil()


# And now we put it back into the AnnData object.

# In[10]:


adata.X[smart_seq_index, :] = smart_seq_X


# ## Integration with scvi-tools
# 
# From its documentation `scvi-tools` is described as a package for end-to-end analysis of single-cell omics data primarily developed and maintained by the Yosef Lab at UC Berkeley.
# 
# Here we will use the "single-cell Variational Inference" model or scVI which uses a deep generative model for the integration of spatial transcriptomic data and scRNA-seq data.
# 
# **For comprehensive usage and best practices of scVI please refer to the [doc site](https://docs.scvi-tools.org/en/stable/index.html) of scvi-tools.**
# 
# ### Inspecting data prior to integration
# 
# Let's take a look at the strength of batch effects in our data. For that we will perform bread-and-butter normalization, neighbor graph calculation, and embedding visualization via UMAP.
# 
# But first let's save the read counts in a different layer as we will need them for integration

# In[11]:


adata.layers["counts"] = adata.X.copy()


# Let's do basic data normalization.

# In[12]:


sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.scale(adata, max_value=10)


# Then selection of highly variable genes.

# In[13]:


sc.pp.highly_variable_genes(
    adata,
    n_top_genes=1000,
    flavor="seurat_v3",
    layer="counts",
    batch_key="dataset_id",
    subset=True,
)


# And finally neighbor graph calculation as well as embedding visualization via UMAP.

# In[14]:


sc.tl.pca(adata)
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)


# In[15]:


sc.pl.umap(adata, color="assay")


# In[16]:


sc.pl.umap(adata, color="cell_type")


# You can see that **batch effects are strong** as cells cluster primarily by `assay` and then by `cell_type`. Properly integrated embedding would in principle cluster primarily by `cell_type`, `assay` should at best randomly distributed. 

# ### Data integration with scVI
# 
# Whenever you query and fetch Census data from multiple datasets then integration needs to be performed as evidenced by the batch effects we observed.
# 
# **The paramaters for SCVI used in this notebook were selected to the model run quickly.** For best practices on integration of single-cell data using [scvi-tools](https://docs.scvi-tools.org/en/stable/index.html) please refer to their documentation page. 
# 
# Additionally we recommend reading the article [An integrated cell atlas of the human lung in health and disease](https://www.biorxiv.org/content/10.1101/2022.03.10.483747v1) by Sikkema et al. whom perfomed integration of 43 datasets from Lung.
# 
# Here we focus on the metadata from the Census that can be as batch information for integration.
# 
# #### Integration with batch defined as `dataset_id`
# 
# All cells in the Census are annotated with the dataset they come from in `obs["dataset_id"]`. This is a great place to start for integration.
# 
# So let's run an scVI model and obtain the latent embeddings. First we define our model with batch set as `dataset_id`.

# In[17]:


scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key="dataset_id")
vae = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb", n_hidden=50)


# Now let's train the model with default parameters.

# In[18]:


vae.train(max_epochs=100)


# And finally let's get the latent representations as cell embeddings and use those for UMAP visualization.

# In[19]:


adata.obsm["X_scVI"] = vae.get_latent_representation()


# In[20]:


sc.pp.neighbors(adata, use_rep="X_scVI")
sc.tl.umap(adata)
sc.pl.umap(adata, color="assay")


# In[21]:


sc.pl.umap(adata, color="cell_type")


# Great! You can see that the clustering is no longer mainly driven by assay, albeit still contributing to it.
# 
# 
# #### Integration with batch defined as `dataset_id` + `donor_id`
# 
# Similar to `dataset_id`, all cells in Census are annotated with `donor_id`. The definition of `donor_id` depends on the dataset and it is left to the discretion of data curators. However it is still rich in information and can be used as a batch variable during integration.
# 
# Because `donor_id` is not guaranteed to be unique across all cells of the Census, we strongly recommend concatenating `dataset_id` and `donor_id` and use that as the batch key for scVI.

# In[22]:


adata.obs["dataset_id_donor_id"] = adata.obs["dataset_id"].astype("str") + adata.obs["donor_id"].astype("str")


# In[23]:


scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key="dataset_id_donor_id")
vae = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb", n_hidden=50)


# Now we can train the model with the new batch definition.

# In[24]:


vae.train(max_epochs=100)


# And we get the latent variables as embeddings.

# In[25]:


adata.obsm["X_scVI"] = vae.get_latent_representation()


# In[26]:


sc.pp.neighbors(adata, use_rep="X_scVI")
sc.tl.umap(adata)
sc.pl.umap(adata, color="assay")


# In[27]:


sc.pl.umap(adata, color="cell_type")


# As you can see using `dataset_id` and `donor_id` as batch the cells now mostly cluster by cell type. 
# 
# #### Integration with batch defined as `dataset_id` + `donor_id` + `assay_ontology_term_id` + `suspension_type`
# 
# In some cases one dataset may contain multiple assay types and/or multiple suspension types (cell vs nucleus), and for those it is important to consider these metadata as batches. 
# 
# Therefore, the most comprehensive definition of batch in the Census can be accomplished by combining the cell metadata of `dataset_id`, `donor_id`, `assay_ontology_term_id` and `suspension_type`, the latter will encode the `EFO` ids for assay types.
# 
# In our example, the two datasets that we used only contain cells from one assay each, and one suspension type for all of them. Thus it would not make a difference to include these metadata as part of batch.



# Section: experimental-pytorch

#!/usr/bin/env python
# coding: utf-8

# # Training a PyTorch Model
# 
# This tutorial shows how to train a Logistic Regression model in PyTorch using the Census API's `experimental.ml.ExperimentDataPipe` class. This is intended only to demonstrate the use of the `ExperimentDataPipe`, and not as an example of how to train a biologically useful model.
# 
# This tutorial assumes a basic familiarity with PyTorch and the Census API. See the [Querying and fetching the single-cell data and cell/gene metadata](https://chanzuckerberg.github.io/cellxgene-census/notebooks/api_demo/census_query_extract.html) notebook tutorial for a quick primer on Census API usage.
# 
# **Contents**
# 
# * [Open the Census](#Open-the-Census)
# * [Create a DataLoader](#Create-a-DataLoader)
# * [Define the model](#Define-the-model)
# * [Train the model](#Train-the-model)
# * [Make predictions with the model](#Make-predictions-with-the-model)
# 

# ## Open the Census
# 
# First, obtain a handle to the Census data, in the usual manner:

# In[1]:


import cellxgene_census

census = cellxgene_census.open_soma()


# ## Create an ExperimentDataPipe
# 
# To train a model in PyTorch using this `census` data object, first instantiate an `ExperimentDataPipe` as follows:

# In[2]:


import cellxgene_census.experimental.ml as census_ml
import tiledbsoma as soma

experiment = census["census_data"]["homo_sapiens"]

experiment_datapipe = census_ml.ExperimentDataPipe(
    experiment,
    measurement_name="RNA",
    X_name="raw",
    obs_query=soma.AxisQuery(value_filter="tissue_general == 'tongue' and is_primary_data == True"),
    obs_column_names=["cell_type"],
    batch_size=128,
    shuffle=True,
    soma_chunk_size=10_000,
)


# ### `ExperimentDataPipe` class explained
# 
# This class provides an implementation of PyTorch's [DataPipe interface](https://pytorch.org/data/main/torchdata.datapipes.iter.html), which defines a common mechanism for wrapping and accessing training data from any underlying source. The `ExperimentDataPipe` class encapsulates the details of querying and retrieving Census data from a single SOMA `Experiment` and returning it to the caller as PyTorch Tensors. Most importantly, it retrieves the data lazily from the Census in batches, avoiding having to load the entire training dataset into memory at once. (Note: PyTorch also provides `DataSet` as a legacy interface for wrapping and accessing training data sources, but a `DataPipe` can be used interchangeably.)
# 

# ### `ExperimentDataPipe` parameters explained
# 
# The constructor only requires a single parameter, `experiment`, which is a `soma.Experiment` containing the data of the organism to be used for training.
# 
# To retrieve a subset of the Experiment's data, along either the `obs` or `var` axes, you may specify query filters via the `obs_query` and `var_query` parameters, which are both `soma.AxisQuery` objects.
# 
# The values for the prediction label(s) that you intend to use for training are specified via the `obs_column_names` array.
# 
# The `batch_size` allows you to specify the number of obs rows (cells) to be returned by each return PyTorch tensor. You may exclude this parameter if you want single rows (`batch_size=1`).
# 
# The `shuffle` flag allows you to randomize the ordering of the training data for each training epoch. Note:
# * You should use this flag instead of the `DataLoader` `shuffle` flag, as `DataLoader` does not support shuffling when used with an `IterDataPipe` dataset.
# * PyTorch's TorchData library provides a [Shuffler](https://pytorch.org/data/main/generated/torchdata.datapipes.iter.Shuffler.html) `DataPipe`, which is alternate mechanism one can use to perform shuffling of an `IterDataPipe`. However, the `Shuffler` will not "globally" randomize the training data, as it only "locally" randomizes the ordering of the training data within fixed-size "windows". Due to the layout of Census data, a given "window" of Census data may be highly homogeneous in terms of its `obs` axis attribute values, and so this shuffling strategy may not provide sufficient randomization for certain types of models.
# 
# The `soma_chunk_size` sets the number of rows of data that are retrieved from the Census and held in memory at a given time. This controls
#  the maximum memory usage of the `ExperimentDataPipe`. Smaller values will require less memory but will also result in lower read performance. If you are running out of memory when training a model, try reducing this value. The default is set to retrieve ~1GB of data per chunk, which takes into account how many `var` (gene) columns are being requested. This parameter also affects the granularity of the "global" shuffling step when `shuffle=True` (see ``shuffle`` parameter API docs for details).

# You can inspect the shape of the full dataset, without causing the full dataset to be loaded:

# In[3]:


experiment_datapipe.shape


# ## Split the dataset
# 
# You may split the overall dataset into the typical training, validation, and test sets by using the PyTorch [RandomSplitter](https://pytorch.org/data/main/generated/torchdata.datapipes.iter.RandomSplitter.html#torchdata.datapipes.iter.RandomSplitter) `DataPipe`. Using PyTorch's functional form for chaining `DataPipe`s, this is done as follows:

# In[4]:


train_datapipe, test_datapipe = experiment_datapipe.random_split(weights={"train": 0.8, "test": 0.2}, seed=1)


# ## Create the DataLoader
# 
# With the full set of DataPipe operations chained together, we can now instantiate a PyTorch [DataLoader](https://pytorch.org/docs/stable/data.html#torch.utils.data.DataLoader) on the training data. 

# In[5]:


experiment_dataloader = census_ml.experiment_dataloader(train_datapipe)


# Alternately, you can instantiate a `DataLoader` object directly via its constructor. However, many of the parameters are not usable with iterable-style DataPipes, which is the case for `ExperimentDataPipe`. In particular, the `shuffle`, `batch_size`, `sampler`, `batch_sampler`, `collate_fn` parameters should not be specified. Using `experiment_dataloader` helps enforce correct usage.

# ## Define the model
# 
# With the training data retrieval code now in place, we can move on to defining a simple logistic regression model, using PyTorch's `torch.nn.Linear` class:

# In[6]:


import torch


class LogisticRegression(torch.nn.Module):
    def __init__(self, input_dim, output_dim):
        super(LogisticRegression, self).__init__()  # noqa: UP008
        self.linear = torch.nn.Linear(input_dim, output_dim)

    def forward(self, x):
        outputs = torch.sigmoid(self.linear(x))
        return outputs


# Next, we define a function to train the model for a single epoch:

# In[7]:


def train_epoch(model, train_dataloader, loss_fn, optimizer, device):
    model.train()
    train_loss = 0
    train_correct = 0
    train_total = 0

    for batch in train_dataloader:
        optimizer.zero_grad()
        X_batch, y_batch = batch

        X_batch = X_batch.float().to(device)

        # Perform prediction
        outputs = model(X_batch)

        # Determine the predicted label
        probabilities = torch.nn.functional.softmax(outputs, 1)
        predictions = torch.argmax(probabilities, axis=1)

        # Compute the loss and perform back propagation

        y_batch = y_batch.flatten()
        y_batch = y_batch.to(device)

        train_correct += (predictions == y_batch).sum().item()
        train_total += len(predictions)

        loss = loss_fn(outputs, y_batch.long())
        train_loss += loss.item()
        loss.backward()
        optimizer.step()

    train_loss /= train_total
    train_accuracy = train_correct / train_total
    return train_loss, train_accuracy


# Note the line, `X_batch, y_batch = batch`. Since the `train_dataloader` was configured with `batch_size=16`, these variables will hold tensors of rank 2. The `X_batch` tensor will appear, for example, as:
# 
# ```
# tensor([[0., 0., 0.,  ..., 1., 0., 0.],
#         [0., 0., 2.,  ..., 0., 3., 0.],
#         [0., 0., 0.,  ..., 0., 0., 0.],
#         ...,
#         [0., 0., 0.,  ..., 0., 0., 0.],
#         [0., 1., 0.,  ..., 0., 0., 0.],
#         [0., 0., 0.,  ..., 0., 0., 8.]])
#       
# ```
# 
# For `batch_size=1`, the tensors will be of rank 1. The `X_batch` tensor will appear, for example, as:
# 
# ```
# tensor([0., 0., 0.,  ..., 1., 0., 0.])
# ```
#     
# For `y_batch`, this will contain the user-specified `obs` `cell_type` training labels. By default, these are encoded using a LabelEncoder and it will be a matrix where each column represents the encoded values of each column specified in `obs_column_names` when creating the datapipe (in this case, only the cell type). It will look like this:
# 
# ```
# tensor([1, 1, 3, ..., 2, 1, 4])
# 
# ```
# Note that cell type values are integer-encoded values, which can be decoded using `experiment_datapipe.obs_encoders` (more on this below).
# 

# ## Train the model
# 
# Finally, we are ready to train the model. Here we instantiate the model, a loss function, and an optimization method and then iterate through the desired number of training epochs. Note how the `train_dataloader` is passed into `train_epoch`, where for each epoch it will provide a new iterator through the training dataset.

# In[8]:


device = torch.device("cuda") if torch.cuda.is_available() else torch.device("cpu")

# The size of the input dimension is the number of genes
input_dim = experiment_datapipe.shape[1]

# The size of the output dimension is the number of distinct cell_type values
cell_type_encoder = experiment_datapipe.obs_encoders["cell_type"]
output_dim = len(cell_type_encoder.classes_)

model = LogisticRegression(input_dim, output_dim).to(device)
loss_fn = torch.nn.CrossEntropyLoss()
optimizer = torch.optim.Adam(model.parameters(), lr=1e-05)

for epoch in range(10):
    train_loss, train_accuracy = train_epoch(model, experiment_dataloader, loss_fn, optimizer, device)
    print(f"Epoch {epoch + 1}: Train Loss: {train_loss:.7f} Accuracy {train_accuracy:.4f}")


# ## Make predictions with the model
# 
# To make predictions with the model, we first create a new `DataLoader` using the `test_datapipe`, which provides the "test" split of the original dataset. For this example, we will only make predictions on a single batch of data from the test split.

# In[9]:


experiment_dataloader = census_ml.experiment_dataloader(test_datapipe)
X_batch, y_batch = next(iter(experiment_dataloader))


# Next, we invoke the model on the `X_batch` input data and extract the predictions:

# In[10]:


model.eval()

model.to(device)
outputs = model(X_batch.to(device))

probabilities = torch.nn.functional.softmax(outputs, 1)
predictions = torch.argmax(probabilities, axis=1)

display(predictions)


# The predictions are returned as the encoded values of `cell_type` label. To recover the original cell type labels as strings, we decode using the encoders from `experiment_datapipe.obs_encoders`.
# 
# At inference time, if the model inputs are not obtained via an `ExperimentDataPipe`, one could pickle the encoder at training time and save it along with the model. Then, at inference time it can be unpickled and used as shown below.

# In[11]:


cell_type_encoder = experiment_datapipe.obs_encoders["cell_type"]

predicted_cell_types = cell_type_encoder.inverse_transform(predictions.cpu())

display(predicted_cell_types)


# Finally, we create a Pandas DataFrame to examine the predictions:

# In[12]:


import pandas as pd

display(
    pd.DataFrame(
        {
            "actual cell type": cell_type_encoder.inverse_transform(y_batch.ravel().numpy()),
            "predicted cell type": predicted_cell_types,
        }
    )
)


# In[ ]:







# Section: analysis_demo-comp_bio_explore_and_load_lung_data

#!/usr/bin/env python
# coding: utf-8

# # Exploring all data from a tissue
# 
# This tutorial provides a series of examples for how to explore and query the Census in the context of a single tissue, lung. We will summarize cell and gene metadata, then fetch the single-cell expression counts and perform some basic data explorations via [Scanpy](https://scanpy.readthedocs.io/en/stable/) 
# 
# 
# **Contents**
# 
# 1. Learning about the human lung data.
#    1. Learning about cells of the lung.
#    2. Learning about genes of the lung .
# 2. Fetching all single-cell human lung data from the Census.
# 3. Calculating QC metrics of the lung data.
# 4. Creating a normalized expression layer and embeddings.
# 
# ⚠️ Note that the Census RNA data includes duplicate cells present across multiple datasets. Duplicate cells can be filtered in or out using the cell metadata variable `is_primary_data` which is described in the [Census schema](https://github.com/chanzuckerberg/cellxgene-census/blob/main/docs/cellxgene_census_schema.md#repeated-data).
# 
# ## Learning about the lung data in the Census
# 
# First we will open the Census. If you are not familiar with the basics of the Census API you should take a look at notebook [Learning about the CZ CELLxGENE Census](https://cellxgene-census.readthedocs.io/en/latest/notebooks/analysis_demo/comp_bio_census_info.html)
# 

# In[1]:


import cellxgene_census
import numpy as np
import pandas as pd
import scanpy as sc

census = cellxgene_census.open_soma()


# Let's first take a look at the number of cells from human lung:
# 

# In[2]:


summary_table = census["census_info"]["summary_cell_counts"].read().concat().to_pandas()

summary_table.query("organism == 'Homo sapiens' & category == 'tissue_general' & label =='lung'")


# There you can see the total of cells of under `total_cell_count` and the unique number cells under `unique_cell_count` (i.e. after removing cells that were included in multiple datasets).
# 
# Let's now take a look at the cell and gene information of this slice of the Census.
# 
# ### Learning about cells of lung data
# 
# Let's load the cell metadata for all lung cells and select only the unique cells using `is_primary_data`.
# 

# In[3]:


lung_obs = cellxgene_census.get_obs(
    census, "homo_sapiens", value_filter="tissue_general == 'lung' and is_primary_data == True"
)
lung_obs


# You can see that the number or rows represents the total number of unique lung cells in the Census. Now let's take a deeper dive into the characteristics of these cells.
# 
# #### Datasets
# 
# First let's start by looking at what are the datasets and collections from [CELLxGENE Discover](https://cellxgene.cziscience.com/collections) contributing to lung. For this we will use the dataset table at `census["census-info"]["datasets"]` that contains metadata of all datasets used to build this Census.
# 

# In[4]:


census_datasets = (
    census["census_info"]["datasets"]
    .read(column_names=["collection_name", "dataset_title", "dataset_id", "soma_joinid"])
    .concat()
    .to_pandas()
)
census_datasets = census_datasets.set_index("dataset_id")
census_datasets


# The `obs` cell metadata `pandas.DataFrame` contains a column `dataset_id` that can be used for joining to the `census_dataset` `pandas.DataFrame` we just created.
# 
# So let's take a look at the cell counts per `dataset_id` of the lung slice and then join to the dataset table to append the human-readable labels.
# 

# In[5]:


dataset_cell_counts = pd.DataFrame(lung_obs[["dataset_id"]].value_counts())
dataset_cell_counts = dataset_cell_counts.rename(columns={0: "cell_counts"})
dataset_cell_counts = dataset_cell_counts.merge(census_datasets, on="dataset_id")

dataset_cell_counts


# These are all the datasets lung cells whose counts are reprensented in the column `cell_counts`. The top collections with lung data are:
# 
# 1. [The integrated Human Lung Cell Atlas](https://cellxgene.cziscience.com/collections/6f6d381a-7701-4781-935c-db10d30de293).
# 2. [A human cell atlas of fetal gene expression](https://cellxgene.cziscience.com/collections/c114c20f-1ef4-49a5-9c2e-d965787fb90c).
# 3. [High-resolution single-cell atlas reveals diversity and plasticity of tumor-associated neutrophils in non-small cell lung cancer](https://cellxgene.cziscience.com/collections/edb893ee-4066-4128-9aec-5eb2b03f8287).
# 4. [HTAN MSK - Single cell profiling reveals novel tumor and myeloid subpopulations in small cell lung cancer](https://cellxgene.cziscience.com/collections/62e8f058-9c37-48bc-9200-e767f318a8ec).
# 5. [A human fetal lung cell atlas uncovers proximal-distal gradients of differentiation and key regulators of epithelial fates.](https://cellxgene.cziscience.com/collections/2d2e2acd-dade-489f-a2da-6c11aa654028).
# 
# #### Assays
# 
# Let's use similar logic to take a look at all the assays available for human lung data. This tells us that most assays are from 10x technologies and sci-RNA-seq.
# 

# In[6]:


lung_obs[["assay"]].value_counts()


# #### Disease
# 
# And now let's take a look at diseased cell counts, with `normal` indicating non-diseased cells.
# 

# In[7]:


lung_obs[["disease"]].value_counts()


# #### Sex
# 
# There doesn't seem to be strong biases for sex.
# 

# In[8]:


lung_obs[["sex"]].value_counts()


# #### Cell vs nucleus
# 
# The majority of data are from cells and not nucleus.
# 

# In[9]:


lung_obs[["suspension_type"]].value_counts()


# #### Cell types
# 
# Let's take a look at the counts of the top 20 cell types.
# 

# In[10]:


lung_obs[["cell_type"]].value_counts().head(20)


# #### Sub-tissues
# 
# We can look at the original tissue annotations that were mapped to "lung".
# 

# In[11]:


lung_obs[["tissue"]].value_counts()


# ### Learning about genes of lung data
# 
# Let's load the gene metadata of the Census.
# 

# In[12]:


lung_var = cellxgene_census.get_var(census, "homo_sapiens")
lung_var


# You can see the total number of genes represented by the number of rows. This number is actually misleading because it is the join of all genes in the Census. However we know that the lung data comes from a subset of datasets.
# 
# So let's take a look at the number of genes that were measured in each of those datasets.
# 
# To accomplish this we can use the "dataset presence matrix" at `census["census_data"]["homo_sapiens"].ms["RNA"]["feature_dataset_presence_matrix"]`. This is a boolean matrix `N x M` where `N` is the number of datasets and `M` is the number of genes in the Census.
# 
# So we can select the rows corresponding to the lung datasets and perform a row-wise sum.
# 

# In[13]:


presence_matrix = cellxgene_census.get_presence_matrix(census, "Homo sapiens", "RNA")
presence_matrix = presence_matrix[dataset_cell_counts.soma_joinid, :]


# In[14]:


presence_matrix.sum(axis=1).A1


# In[15]:


genes_measured = presence_matrix.sum(axis=1).A1
dataset_cell_counts["genes_measured"] = genes_measured
dataset_cell_counts


# You can see the genes measured in each dataset represented in `genes_measured`. Now lets get the **genes that were measured in all datasets**.
# 

# In[16]:


var_somaid = np.nonzero(presence_matrix.sum(axis=0).A1 == presence_matrix.shape[0])[0].tolist()


# In[17]:


lung_var = lung_var.query(f"soma_joinid in {var_somaid}")
lung_var


# The number of rows represents the genes that were measured in all lung datasets.
# 
# ### Summary of lung metadata
# 
# In the previous sections, using the Census we learned the following information:
# 
# - The total number of unique lung cells and their composition for:
#   - Number of datasets.
#   - Number sequencing technologies, most of which are 10x
#   - Mostly human data, but some diseases exist, primarily "lung adenocarcinoma" and "COVID-19 infected"
#   - No sex biases.
#   - Mostly data from cells (\~80%) rather than nucleus (\~20%)
# - A total of **~12k** genes were measured across all cells.
# 
# ##  Fetching all single-cell human lung data from the Census
# 
# Since loading the entire lung data is resource-intensive, for the sake of this exercise let's load a subset of the lung data into an `anndata.AnnData` object and perform some exploratory analysis. 
# 
# We will subset to 100,000 random unique cells using the `lung_obs` `pandas.DataFrame` we previously created.

# In[18]:


lung_cell_subsampled_n = 100000
lung_cell_subsampled_ids = lung_obs["soma_joinid"].sample(lung_cell_subsampled_n, random_state=1).tolist()


# Now we can directly use the values of `soma_joinid` for querying the Census data and obtaining an `AnnData` object.

# In[19]:


lung_gene_ids = lung_var["soma_joinid"].to_numpy()
lung_adata = cellxgene_census.get_anndata(
    census,
    organism="Homo sapiens",
    obs_coords=lung_cell_subsampled_ids,
    var_coords=lung_gene_ids,
)

lung_adata.var_names = lung_adata.var["feature_name"]


# In[20]:


lung_adata


# We are done with the census, so close it

# In[21]:


census.close()
del census


# ## Calculating QC metrics of the lung data
# 
# Now let's take a look at some QC metrics
# 
# **Top genes per cell**
# 

# In[22]:


sc.pl.highest_expr_genes(lung_adata, n_top=20)


# **Number of sequenced genes by assay**
# 

# In[23]:


sc.pp.calculate_qc_metrics(lung_adata, percent_top=None, log1p=False, inplace=True)
sc.pl.violin(lung_adata, "n_genes_by_counts", groupby="assay", jitter=0.4, rotation=90)


# **Total counts by assay**
# 

# In[24]:


sc.pl.violin(lung_adata, "total_counts", groupby="assay", jitter=0.4, rotation=90)


# You can see that Smart-Seq2 is an outlier for the total counts per cell, so let's exlcude it to see how the rest of the assays look like
# 

# In[25]:


sc.pl.violin(
    lung_adata[lung_adata.obs["assay"] != "Smart-seq2",],
    "total_counts",
    groupby="assay",
    jitter=0.4,
    rotation=90,
)


# ## Creating a normalized expression layer and embeddings
# 
# Let's perform a bread and butter normalization and take a look at UMAP embeddings, but for all the data below we'll exclude Smart-seq2 as this requires an extra step to normalize based on gene lengths
# 

# In[26]:


lung_adata = lung_adata[lung_adata.obs["assay"] != "Smart-seq2",].copy()
lung_adata.layers["counts"] = lung_adata.X


# Now let's do some basic normalization:
# 
# - Normalize by sequencing depth
# - Transform to log-scale
# - Select 500 highly variable genes
# - Scale values across the gene axis
# 

# In[27]:


sc.pp.normalize_total(lung_adata, target_sum=1e4)
sc.pp.log1p(lung_adata)
sc.pp.highly_variable_genes(lung_adata, n_top_genes=500, flavor="seurat_v3", layer="counts")
lung_adata = lung_adata[:, lung_adata.var.highly_variable]
sc.pp.scale(lung_adata, max_value=10)


# And reduce dimensionality by obtaining UMAP embeddings.
# 

# In[28]:


sc.tl.pca(lung_adata)
sc.pp.neighbors(lung_adata)
sc.tl.umap(lung_adata)


# And plot these embeddings.
# 

# In[29]:


n_cell_types = len(lung_adata.obs["cell_type"].drop_duplicates())

from random import randint

colors = []

for i in range(len(lung_adata.obs["cell_type"].drop_duplicates())):
    colors.append("#%06X" % randint(0, 0xFFFFFF))


# In[30]:


sc.pl.umap(lung_adata, color="cell_type", palette=colors, legend_loc=None)


# Let's color by assay.
# 

# In[31]:


sc.pl.umap(lung_adata, color="assay")


# Given the high number of cell types it makes it hard to visualize, so let's look at the top 20 most abundant cell types.
# 

# In[32]:


top_cell_types = lung_adata.obs["cell_type"].value_counts()
top_cell_types = list(top_cell_types.reset_index().head(20)["cell_type"])


# In[33]:


lung_adata_top_cell_types = lung_adata[[i in top_cell_types for i in lung_adata.obs["cell_type"]], :]
sc.pl.umap(lung_adata_top_cell_types, color="cell_type")


# Let's color by assay of this subset of the data.
# 

# In[34]:


sc.pl.umap(lung_adata_top_cell_types, color="assay")




# Section: api_demo-census_query_extract

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




# Section: analysis_demo-comp_bio_normalizing_full_gene_sequencing

#!/usr/bin/env python
# coding: utf-8

# # Normalizing full-length gene sequencing data
# 
# This tutorial shows you how to fetch full-length gene sequencing data from the Census and normalize it to account for gene length.
# 
# **Contents**
# 
# 1. Opening the census
# 2. Fetching example full-length sequencing data (Smart-Seq2)
# 3. Normalizing expression to account for gene length
# 4. Validation through clustering exploration
# 
# ⚠️ Note that the Census RNA data includes duplicate cells present across multiple datasets. Duplicate cells can be filtered in or out using the cell metadata variable `is_primary_data` which is described in the [Census schema](https://github.com/chanzuckerberg/cellxgene-census/blob/main/docs/cellxgene_census_schema.md#repeated-data).
# For this notebook we will focus on individual datasets, therefore we can ignore this variable.
# 
# ## Opening the census
# 
# First we open the Census, if you are not familiar with the basics of the Census API you should take a look at notebook [Learning about the CZ CELLxGENE Census](https://chanzuckerberg.github.io/cellxgene-census/notebooks/analysis_demo/comp_bio_census_info.html)

# In[1]:


import cellxgene_census
import scanpy as sc
from scipy.sparse import csr_matrix

census = cellxgene_census.open_soma()


# You can learn more about the all of the `cellxgene_census` methods by accessing their corresponding documention via `help()`. For example `help(cellxgene_census.open_soma)`. 

# ## Fetching full-length example sequencing data (Smart-Seq)
# 
# Let's get some example data, in this case we'll fetch all cells from a relatively small dataset derived from the Smart-Seq2 technology which performs full-length gene sequencing:
# 
# - Collection: [Tabula Muris Senis](https://cellxgene.cziscience.com/collections/0b9d8a04-bb9d-44da-aa27-705bb65b54eb)
# - Dataset: [Liver - A single-cell transcriptomic atlas characterizes ageing tissues in the mouse - Smart-seq2](https://cellxgene.cziscience.com/e/524179b0-b406-4723-9c46-293ffa77ca81.cxg/)

# Let's first find this dataset's id by using the dataset table of the Census

# In[2]:


liver_dataset = (
    census["census_info"]["datasets"]
    .read(
        value_filter="dataset_title == 'Liver - A single-cell transcriptomic atlas characterizes ageing tissues in the mouse - Smart-seq2'"
    )
    .concat()
    .to_pandas()
)
liver_dataset


# Now we can use this id to fetch the data.

# In[3]:


liver_dataset_id = "4546e757-34d0-4d17-be06-538318925fcd"
liver_adata = cellxgene_census.get_anndata(
    census,
    organism="Mus musculus",
    obs_value_filter=f"dataset_id=='{liver_dataset_id}'",
)


# Let's make sure this data only contains Smart-Seq2 cells

# In[4]:


liver_adata.obs["assay"].value_counts()


# Great! As you can see this a small dataset only containing **2,859** cells. Now let's proceed to normalize by gene lengths.
# 
# For good practice let's add the number of genes expressed by cell, and the total counts per cell.

# Don't forget to close the census

# In[5]:


census.close()


# ## Normalizing expression to account for gene length
# 
# By default `cellxgene_census.get_anndata()` fetches all genes in the Census. So let's first identify the genes that were measured in this dataset and subset the `AnnData` to only include those. 
# 
# To this goal we can use the "Dataset Presence Matrix" in `census["census_data"]["mus_musculus"].ms["RNA"]["feature_dataset_presence_matrix"]`. This is a boolean matrix `N x M` where `N` is the number of datasets and `M` is the number of genes in the Census, `True` indicates that a gene was measured in a dataset.

# In[6]:


liver_adata.n_vars


# Let's get the genes measured in this dataset.

# In[7]:


liver_dataset_id = liver_dataset["soma_joinid"][0]

presence_matrix = cellxgene_census.get_presence_matrix(census, "Mus musculus", "RNA")
presence_matrix = presence_matrix[liver_dataset_id, :]
gene_presence = presence_matrix.nonzero()[1]

liver_adata = liver_adata[:, gene_presence].copy()


# In[8]:


liver_adata.n_vars


# We can see that out of the all genes in the Census **17,992** were measured in this dataset. 
# 
# Now let's normalize these genes by gene length. We can easily do this because the Census has gene lengths included in the gene metadata under `feature_length`.

# In[9]:


liver_adata.X[:5, :5].toarray()


# In[10]:


gene_lengths = liver_adata.var[["feature_length"]].to_numpy()
liver_adata.X = csr_matrix((liver_adata.X.T / gene_lengths).T)


# In[11]:


liver_adata.X[:5, :5].toarray()


# All done! You can see that we now have real numbers instead of integers.
# 
# ## Validation through clustering exploration
# 
# Let's perform some basic clustering analysis to see if cell types cluster as expected using the normalized counts.
# 
# First we do some basic filtering of cells and genes.

# In[12]:


sc.pp.filter_cells(liver_adata, min_genes=500)
sc.pp.filter_genes(liver_adata, min_cells=5)


# Then we normalize to account for sequencing depth and transform data to log scale.

# In[13]:


sc.pp.normalize_total(liver_adata, target_sum=1e4)
sc.pp.log1p(liver_adata)


# Then we subset to highly variable genes.

# In[14]:


sc.pp.highly_variable_genes(liver_adata, n_top_genes=1000)
liver_adata = liver_adata[:, liver_adata.var.highly_variable].copy()


# And finally we scale values across the gene axis.

# In[15]:


sc.pp.scale(liver_adata, max_value=10)


# Now we can proceed to do clustering analysis

# In[16]:


sc.tl.pca(liver_adata)
sc.pp.neighbors(liver_adata)
sc.tl.umap(liver_adata)
sc.pl.umap(liver_adata, color="cell_type")


# With a few exceptions we can see that all cells from the same cell type cluster near each other which serves as a sanity check for the gene-length normalization that we applied.



# Section: api_demo-census_summary_cell_counts

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




# Section: api_demo-census_embedding

#!/usr/bin/env python
# coding: utf-8

# # Access CELLxGENE-hosted embeddings
# 
# This notebook demonstrates basic access to CELLxGENE-hosted embeddings of the Census. **CELLxGENE-hosted embeddings have been contributed by the community**, CELLxGENE Discover does not actively maintain or update them. Find out more about these in the [Census model page](https://cellxgene.cziscience.com/census-models). 
# 
# **IMPORTANT:** This tutorial requires cellxgene-census package version 1.9.1 or later.
# 
# ## Contents
# 
# 1. Background
# 2. Quick start
# 3. Query cells and load associated embeddings
# 4. Load embeddings and fetch associated Census data
# 5. Embedding metadata
# 
# > ⚠️ Note that the Census RNA data includes duplicate cells present across multiple datasets. Duplicate cells can be filtered in or out using the cell metadata variable `is_primary_data` which is described in the [Census schema](https://github.com/chanzuckerberg/cellxgene-census/blob/main/docs/cellxgene_census_schema.md#repeated-data).
# 
# 
# ## Background
# 
# This notebook demonstrates access to CELLxGENE-hosted embeddings of the Census. The Census has multiple releases, named by a `census_version`, which normally looks like an ISO date, e.g., `2023-02-01`. A CELLxGENE-hosted embedding is a 2D sparse matrix of cell embeddings _for a given census version_, encoded as a SOMA SparseNDArray.
# 
# ⚠️ Note that embeddings may be available for one or both organisms, see the [Census model page](https://cellxgene.cziscience.com/census-models) for the latest availability.
# 
# ⚠️ **IMPORTANT:** embeddings are only meaningful in the context of the Census from which they were created. Each embedding contains a metadata field indicating the source Census, suitable for confirming embedding lineage.
# 
# ## Quick start
# 
# The easiest way to access Census CELLxGENE-hosted embeddings is by calling the `get_anndata` function with an `obs_embeddings` or `var_embeddings` parameter.
# 
# Let's start by exploring what embeddings are available:

# In[1]:


from cellxgene_census.experimental import get_all_available_embeddings

CENSUS_VERSION = "2023-12-15"

for e in get_all_available_embeddings(CENSUS_VERSION):
    print(f"{e['embedding_name']:15} {e['experiment_name']:15} {e['data_type']:15}")


# These can also be viewed on the [CELLxGENE Census Models page](https://cellxgene.cziscience.com/census-models). 
# 
# For this example, we'll use `scgpt`. We can call `get_anndata` with the `obs_embeddings` parameter:

# In[2]:


import cellxgene_census
from cellxgene_census.experimental import get_embedding

CENSUS_VERSION = "2023-12-15"

with cellxgene_census.open_soma(census_version=CENSUS_VERSION) as census:
    adata = cellxgene_census.get_anndata(
        census,
        organism="homo_sapiens",
        measurement_name="RNA",
        obs_value_filter="tissue == 'tongue'",
        obs_embeddings=["scgpt"],
    )


# And now you can use these data for downstream analysis.

# In[3]:


adata


# In[4]:


adata.obsm


# ⚠️ **IMPORTANT**: `get_embedding` will fill in missing values in the embedding matrix with `NaN`, therefore missing cells are represented with rows where all values are `NaN`. To learn more use `help(get_embedding)`.

# ## Storage format
# 
# Each embedding is encoded as a SOMA SparseNDArray, where:
# 
# * dimension 0 (`soma_dim_0`) encodes the cell (obs) `soma_joinid` value
# * dimension 1 (`soma_dim_1`) encodes the embedding feature, and is in the range [0, N) where N is the number of featues in the embedding
# * data (`soma_data`) is float32
# 
# ⚠️ **IMPORTANT:** CELLxGENE-hosted embeddings may embed a subset of the cells in any given Census version. If a cell has an embedding, it will be explicitly stored in the sparse array, _even if the embedding value is zero_. In other words, missing array values values imply that the cell was not embedded, whereas zero valued embeddings are explicitly stored. Put another way, the `nnz` of the embedding array indicate the number of embedded cells, not the number of non-zero values.
# 
# The first axis of the embedding array will have the same shape as the corresponding `obs` DataFrame for the Census build and experiment. The second axis of the embedding will have a shape (0, N) where N is the number of features in the embedding.
# 
# Embedding values, while stored as a float32, are precision reduced. Currently they are equivalent to a bfloat16, i.e., have 8 bits of exponent and 7 bits of mantissa.
# 

# ## Query cells and load associated embeddings
# 
# This section demonstrates two methods to query cells from the Census by `obs` metadata, and then fetch CELLxGENE-hosted embeddings associated with each cell.
# 
# 1. Load an embedding into an AnnData `obsm` slot
# 2. Load an embedding into a dense NumPy array
# 
# Let's first do a few imports and utility functions used throughout this notebook.

# In[5]:


# A few imports and utility functions used throughout this notebook


import warnings

import cellxgene_census
import numpy as np
import scanpy
import tiledbsoma as soma
from cellxgene_census.experimental import get_embedding_metadata

warnings.filterwarnings("ignore")

# The Census version to utilize
CENSUS_VERSION = "2023-12-15"
EXPERIMENT_NAME = "homo_sapiens"
MEASUREMENT_NAME = "RNA"

# The location of the embedding.  See <https://cellxgene.cziscience.com/census-models>
# for available CELLxGENE-hosted embeddings.
EMBEDDING_URI = "s3://cellxgene-contrib-public/contrib/cell-census/soma/2023-12-15/CxG-contrib-1/"


# ### Load an embedding into an AnnData `obsm` slot
# 
# There are two main ways to load hosted embeddings into an AnnData.
# 
# 1. Via `cellxgene_census.get_anndata()`, followed by merging embeddings.
# 2. With a lazy query via `ExperimentAxisQuery`, followed by merging embeddings.
# 
# 
# #### AnnData embeddings via cellxgene_census.get_anndata()
# 
# This is the simplest way of getting the embeddings. In this example we create an AnnData for all "central nervous system" cells, and use the `obs_embeddings` parameter to add scGPT embeddings to the `obsm` slot.
# 

# In[6]:


with cellxgene_census.open_soma(census_version=CENSUS_VERSION) as census:
    adata = cellxgene_census.get_anndata(
        census,
        organism=EXPERIMENT_NAME,
        measurement_name=MEASUREMENT_NAME,
        obs_value_filter="tissue_general == 'central nervous system'",
        obs_column_names=["cell_type", "soma_joinid"],
        obs_embeddings=["scgpt"],
    )


# ⚠️ **IMPORTANT**: `get_embedding` will fill in missing values in the embedding matrix with `NaN`, therefore missing cells are represented with rows where all values are `NaN`. To learn more use `help(get_embedding)`.

# Then we can take a quick look at the embeddings in a 2D scatter plot via UMAP.

# In[7]:


scanpy.pp.neighbors(adata, use_rep="scgpt")
scanpy.tl.umap(adata)
scanpy.pl.umap(adata, color="cell_type", title="scGPT")


# #### AnnData embeddings via `ExperimentAxisQuery`
# 
# Using an `ExperimentAxisQuery` to get embeddings into an AnnData has the main advantage of inspecting the query in a lazy manner before loading all data into AnnData.
# 
# As a reminder this class offers a lazy interface to query Census based on cell and gene metadata, and provides access to the correspondong expression data, and cell/gene metadata.
# 
# Let's initiate a lazy query with the same filters as the previous example.

# In[8]:


census = cellxgene_census.open_soma(census_version=CENSUS_VERSION)

experiment = census["census_data"][EXPERIMENT_NAME]
query = experiment.axis_query(
    measurement_name="RNA",
    obs_query=soma.AxisQuery(value_filter="tissue_general == 'central nervous system'"),
)


# Now, before downloading all the data we can take a look at different attributes, for example the number of cells in our query.

# In[9]:


query.n_obs


# And for example grab all the `soma_joinid` values for the cells in this query.

# In[10]:


soma_joinids = query.obs_joinids().to_numpy()


# Then create an AnnData.

# In[11]:


adata = query.to_anndata(X_name="raw", column_names={"obs": ["cell_type"]})


# And finally add the embedding matrix via `get_embedding` for the cells of the query using the `soma_joinid` values.

# In[12]:


adata.obsm["scgpt"] = get_embedding(
    census_version=CENSUS_VERSION,
    embedding_uri=EMBEDDING_URI,
    obs_soma_joinids=soma_joinids,
)


# ⚠️ **IMPORTANT**: `get_embedding` will fill in missing values in the embedding matrix with `NaN`, therefore missing cells are represented with rows where all values are `NaN`. To learn more use `help(get_embedding)`.

# In[13]:


adata


# In[14]:


query.close()
census.close()


# ### Load an embedding into a dense NumPy array
# 
# To load a  embeddinng into a stand-alone numpy array you can select cells from the Census based on obs metadata, then given the resulting cells, use the `soma_joinid` values to download an embedding, and finally save as a dense NDArray.
# 
# Let's first select cells based on cell metadata

# In[15]:


census = cellxgene_census.open_soma(census_version=CENSUS_VERSION)

obs_df = cellxgene_census.get_obs(
    census,
    EXPERIMENT_NAME,
    value_filter="tissue_general == 'exocrine gland'",
    column_names=["soma_joinid", "cell_type"],
)


# Now you can use the `soma_joinid` values to download the corresponding rows of the embedding matrix via `get_embedding`.

# In[16]:


embeddings = get_embedding(CENSUS_VERSION, EMBEDDING_URI, obs_df.soma_joinid.to_numpy())
embeddings[0:3, 0:4]


# ⚠️ **IMPORTANT**: `get_embedding` will fill in missing values in the embedding matrix with `NaN`, therefore missing cells are represented with rows where all values are `NaN`.  

# In[17]:


embeddings.shape


# ## Load embeddings and fetch associated Census data
# 
# This section describes a more advanced use case. Here we showcase how to load large slices of an embeding matrix, and then append cell metadata to them.
# 
# The method starts with the loaded embedding, and for each embedded cell loads metadata or X data.

# In[18]:


# Load a portion of the embedding (caution: embeddings can be quite large)

import cellxgene_census
import tiledbsoma as soma

# Fetch first 500_000 joinids from the embedding.
# NOTE: will fail if the there are no cells embedded within this obs joinid range
embedding_slice = (slice(500_000),)

emb_data = []
emb_joinids = []
ctx = {"vfs.s3.region": "us-west-2", "vfs.s3.no_sign_request": True}
with soma.open(EMBEDDING_URI, context=soma.SOMATileDBContext(tiledb_config=ctx)) as E:
    # read embedding and obs joinids for each embedded cell
    for d, (obs_joinids, _) in (
        E.read(coords=embedding_slice).blockwise(axis=0, size=2**20, reindex_disable_on_axis=1).scipy()
    ):
        embedding_presence_mask = d.getnnz(axis=1) != 0
        emb_joinids.append(obs_joinids[embedding_presence_mask])
        emb_data.append(d[embedding_presence_mask, :].toarray())

    # concat
    embedding_data = np.vstack(emb_data)
    embedding_joinids = np.concatenate(emb_joinids)

# Load the associated metadata - in this case, obs.suspension_type
with cellxgene_census.open_soma(census_version=CENSUS_VERSION) as census:
    experiment = census["census_data"][EXPERIMENT_NAME]
    with experiment.axis_query(
        measurement_name=MEASUREMENT_NAME,
        obs_query=soma.AxisQuery(coords=(embedding_joinids,)),
    ) as query:
        obs_df = query.obs(column_names=["soma_joinid", "suspension_type"]).concat().to_pandas()

# Display the first few cells, and the first few columns of their embeddings
display(obs_df[0:3])
display(embedding_data[0:3, 0:4])


# ## Embedding Metadata
# 
# Each embedding contains descriptive information stored in the SOMA `metadata` slot, encoded as a JSON string. This metadata includes:
# 
# * census_version - the Census which is embedded. It is critical to confirm this matches the Census in use, or the embeddings will be meaningless.
# * experiment_name - the Census experiment embedded, e.g., `homo_sapiens` or `mus_musculus`.
# * measurement_name - the Census measurement embedded, e.g., `RNA`
# 
# There are a variety of other metadata values [documented here](https://github.com/chanzuckerberg/cellxgene-census/blob/main/tools/census_contrib/embedding_metadata.md).
# 
# The following example demonstrates how to access and decode the metadata into a Python dictionary.

# In[19]:


embedding_metadata = get_embedding_metadata(EMBEDDING_URI)

embedding_metadata


# A common use case for accessing embedding metadata is the prevention of nonsense results, which will occur when the embeddings were not derived from the Census or experiment in use.
# 
# This demonstrates an easy way to grab the embedding metadata and confirm that it matches your expectations.

# In[20]:


embedding_metadata = get_embedding_metadata(EMBEDDING_URI)

assert embedding_metadata["census_version"] == CENSUS_VERSION
assert embedding_metadata["experiment_name"] == EXPERIMENT_NAME
assert embedding_metadata["measurement_name"] == MEASUREMENT_NAME

display("all good!")




# Section: api_demo-census_access_maintained_embeddings

#!/usr/bin/env python
# coding: utf-8

# # Access CELLxGENE collaboration embeddings (scVI, Geneformer)
# 
# This notebook demonstrates basic access to CELLxGENE collaboration embeddings of CELLxGENE Discover Census. Currently, embeddings from scVI and a fine-tuned Geneformer model are maintained by CELLxGENE Discover. There are other CELLxGENE-hosted embeddings contributed by the community to CELLxGENE Discover, find out more about these in the [Census model page](https://cellxgene.cziscience.com/census-models).
# 
# **IMPORTANT:** This tutorial requires cellxgene-census package version 1.9.1 or later.
# 
# **Contents**
# 
# 1. Quick start
# 1. Storage format
# 1. Query cells and load associated embeddings
# 
# ⚠️ Note that the Census RNA data includes duplicate cells present across multiple datasets. Duplicate cells can be filtered in or out using the cell metadata variable `is_primary_data` which is described in the [Census schema](https://github.com/chanzuckerberg/cellxgene-census/blob/main/docs/cellxgene_census_schema.md#repeated-data).
# 
# 
# ## Quick start
# 
# CELLxGENE collaboration embeddings can easily be exported into an AnnData as shown below for any slice of Census. This example queries all cells from tongue tissue. 
# 
# ⚠️ Note that Geneformer embeddings are only available for human data

# In[1]:


import cellxgene_census

emb_names = ["scvi", "geneformer"]
census_version = "2023-12-15"

with cellxgene_census.open_soma(census_version=census_version) as census:
    adata = cellxgene_census.get_anndata(
        census,
        organism="homo_sapiens",
        measurement_name="RNA",
        obs_value_filter="tissue == 'tongue'",
        obs_embeddings=emb_names,
    )


# In[2]:


adata


# In[3]:


adata.obsm


# ## Storage format
# 
# Each embedding is encoded as a SOMA SparseNDArray, where:
# 
# * dimension 0 (`soma_dim_0`) encodes the cell (obs) `soma_joinid` value
# * dimension 1 (`soma_dim_1`) encodes the embedding feature, and is in the range [0, N) where N is the number of features in the embedding
# * data (`soma_data`) is float32
# 
# The first axis of the embedding array will have the same shape as the corresponding `obs` DataFrame for the Census build and experiment. The second axis of the embedding will have a shape (0, N) where N is the number of features in the embedding.
# 
# Embedding values, while stored as a float32, are precision reduced. Currently they are equivalent to a bfloat16, i.e., have 8 bits of exponent and 7 bits of mantissa.

# ## Query cells and load associated embeddings
# 
# This section demonstrates several methods to query cells from the Census by `obs` metadata, and then fetch embeddings associated with each cell.
# 

# ### Loading embeddings into an AnnData `obsm` slot

# There are two main ways to load CELLxGENE collaboration embeddings into an AnnData.
# 
# 1. Via `cellxgene_census.get_anndata()`.
# 2. With a lazy query via `ExperimentAxisQuery`.

# #### AnnData embeddings via `cellxgene_census.get_anndata()`

# This is the simplest way of getting the embeddings. In this example we create an AnnData for all "central nervous system" cells.

# In[4]:


import warnings

import cellxgene_census
import scanpy

warnings.filterwarnings("ignore")

census_version = "2023-12-15"
census = cellxgene_census.open_soma(census_version=census_version)

emb_names = ["scvi", "geneformer"]

adata = cellxgene_census.get_anndata(
    census,
    organism="homo_sapiens",
    measurement_name="RNA",
    obs_value_filter="tissue_general == 'central nervous system'",
    obs_column_names=["cell_type"],
    obs_embeddings=emb_names,
)

census.close()


# Then we can take a quick look at the embeddings in a 2D scatter plot via UMAP.

# In[5]:


scanpy.pp.neighbors(adata, use_rep="scvi")
scanpy.tl.umap(adata)
scanpy.pl.umap(adata, color="cell_type", title="scVI")


# In[6]:


scanpy.pp.neighbors(adata, use_rep="geneformer")
scanpy.tl.umap(adata)
scanpy.pl.umap(adata, color="cell_type", title="Geneformer")


# #### AnnData embeddings via `ExperimentAxisQuery`

# Using an `ExperimentAxisQuery` to get embeddings into an AnnData has the main advantage of inspecting the query in a lazy manner before loading all data into AnnData.
# 
# As a reminder this class offers a lazy interface to query Census based on cell and gene metadata, and provides access to the correspondong expression data, cell/gene metadata, and the embeddings.
# 
# Let's initiate a lazy query with the same filters as the previous example.

# In[7]:


import cellxgene_census
import scanpy
import tiledbsoma as soma

census_version = "2023-12-15"
organism = "homo_sapiens"

census = cellxgene_census.open_soma(census_version=census_version)

experiment = census["census_data"][organism]
query = experiment.axis_query(
    measurement_name="RNA",
    obs_query=soma.AxisQuery(value_filter="tissue_general == 'central nervous system'"),
)


# Now, before downloading all the data we can take a look at different attributes, for example the number of cells in our query.

# In[8]:


query.n_obs


# Then we create an AnnData, retrieve the embeddings and add them to the AnnData's obsm slot:

# In[9]:


from cellxgene_census.experimental import get_embedding, get_embedding_metadata_by_name

emb_names = ["scvi", "geneformer"]

adata = query.to_anndata(X_name="raw", column_names={"obs": ["cell_type"]})

for embedding_name in ["scvi", "geneformer"]:
    metadata = get_embedding_metadata_by_name(embedding_name, "homo_sapiens", census_version=census_version)
    embedding_uri = (
        f"s3://cellxgene-contrib-public/contrib/cell-census/soma/{metadata['census_version']}/{metadata['id']}"
    )
    embedding = get_embedding(metadata["census_version"], embedding_uri, query.obs_joinids().to_numpy())
    adata.obsm[embedding_name] = embedding

adata


# In[10]:


adata


# In[11]:


query.close()
census.close()


# ### Load an embedding into a dense NumPy array
# 
# To load a  embeddinng into a stand-alone numpy array you can select cells from the Census based on obs metadata, then given the resulting cells, use the `soma_joinid` values to download an embedding, and finally save as a dense NDArray.
# 
# Let's first select cells based on cell metadata.

# In[12]:


import cellxgene_census
import tiledbsoma as soma

census_version = "2023-12-15"
experiment_name = "homo_sapiens"

census = cellxgene_census.open_soma(census_version=census_version)

obs_df = cellxgene_census.get_obs(
    census,
    experiment_name,
    value_filter="tissue_general == 'central nervous system'",
    column_names=["soma_joinid", "cell_type"],
)


# Now you can use the `soma_joinid` values to download the corresponding rows of the embedding matrix via `get_embedding`.

# In[13]:


metadata = get_embedding_metadata_by_name("scvi", experiment_name, census_version=census_version)
embedding_uri = f"s3://cellxgene-contrib-public/contrib/cell-census/soma/{metadata['census_version']}/{metadata['id']}"
embedding = get_embedding(metadata["census_version"], embedding_uri, obs_df.soma_joinid.to_numpy())


# And now we have a dense matrix with the embedding data.

# In[14]:


embedding


# In[16]:


embedding.shape


# In[17]:


query.close()
census.close()




# Section: experimental-pca

#!/usr/bin/env python
# coding: utf-8

# # Incremental PCA
# 
# The notebook demonstrates the use of [scikit-learn IncrementalPCA](https://scikit-learn.org/stable/modules/generated/sklearn.decomposition.IncrementalPCA.html#sklearn.decomposition.IncrementalPCA) to perform PCA on Census data.
# 
# Approach:
# 
# * Use a SOMA query to define the cells to be embedded,
# * From these cells, select N top genes using the `experimental.pp.highly_variable_genes` method,
# * Incrementally train over the selected cells and the N top genes,
# * Compute components, and annotate the `obs` dataframe.
# 
# Depending on the number of cells and genes selected, this can be a resource intensive computation. It is known to complete succesfully when trained on the top 5000 genes for all cells in the human and mouse Census data, but requires a large host. For example, the full human PCA has been succesfully demonstrated on an AWS EC2 c6id.32xlarge instance.

# In[3]:


import cellxgene_census
import numpy as np
import tiledbsoma as soma
from cellxgene_census.experimental.pp import highly_variable_genes
from sklearn.decomposition import IncrementalPCA

"""
Configuration - the dataset and computational parameters.
"""
census_version = "latest"  # which Census version is used
experiment_name = "mus_musculus"  # which organism: mus_musculus or homo_sapiens
obs_value_filter = "tissue_general == 'heart'"  # the subset of cells (both train and embed). Set to None if all cells.
n_components = 30  # number of components to keep in the final result
n_top_genes = 3000  # number of genes to use as analysis input


# In[5]:


with cellxgene_census.open_soma(census_version=census_version) as census:
    exp = census["census_data"][experiment_name]

    with exp.axis_query(
        measurement_name="RNA",
        obs_query=soma.AxisQuery(value_filter=obs_value_filter),
    ) as query:
        print(f"{query.n_obs} cells selected")
        print("Beginning HVG calculation")
        hvgs = highly_variable_genes(query, n_top_genes=n_top_genes)
        var_soma_joinids = hvgs[hvgs.highly_variable].index.to_numpy()
        del hvgs
        print("Finished HVG calculation")

    with exp.axis_query(
        measurement_name="RNA",
        obs_query=soma.AxisQuery(value_filter=obs_value_filter),
        var_query=soma.AxisQuery(coords=(var_soma_joinids,)),
    ) as query:
        print("Start training")
        pca = IncrementalPCA(n_components=n_components)
        training_chunk_size = 2000
        for n, (chunk, _) in enumerate(query.X("raw").blockwise(axis=0).scipy()):
            for i in range(0, chunk.shape[0], training_chunk_size):
                training_chunk = chunk[i : i + training_chunk_size, :].toarray()
                pca.partial_fit(training_chunk)
        print("End training")

        obs = query.obs(column_names=["soma_joinid"]).concat().to_pandas().set_index("soma_joinid")
        for colname in (f"X_pca_{n}" for n in range(0, n_components)):
            obs[colname] = np.zeros((len(obs),), dtype=np.float64)

        print("Start transform")
        for n, (chunk, (obs_join_ids, _)) in enumerate(query.X("raw").blockwise(axis=0).scipy()):
            chunk_trnsfm = pca.transform(chunk.toarray())
            for c in range(n_components):
                obs.loc[obs_join_ids, f"X_pca_{c}"] = chunk_trnsfm[:, c]
        print("Complete")

obs


# In[7]:


chunk_trnsfm.dtype




# Section: api_demo-census_citation_generation

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




# Section: analysis_demo-comp_bio_census_info

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




# Section: api_demo-census_compute_over_X

#!/usr/bin/env python
# coding: utf-8

# # Computing on X using online (incremental) algorithms
# 
# This tutorial showcases computing a variety of per-gene and per-cell statistics for a user-defined query using out-of-core operations.
# 
# *NOTE*: when query results are small enough to fit in memory, it may be easier to use the `SOMAExperiment` Query class to extract an AnnData, and then just compute over that. This tutorial shows means of incrementally processing larger-than-core (RAM) data, where incremental (online) algorithms are used.
# 
# **Contents**
# 
# 1. Incremental count and mean calculation.
# 2. Incremental variance calculation.
# 3. Counting cells per gene, grouped by `dataset_id`.
# 
# ⚠️ Note that the Census RNA data includes duplicate cells present across multiple datasets. Duplicate cells can be filtered in or out using the cell metadata variable `is_primary_data` which is described in the [Census schema](https://github.com/chanzuckerberg/cellxgene-census/blob/main/docs/cellxgene_census_schema.md#repeated-data).

# In[1]:


import cellxgene_census
import numpy as np
import pandas as pd
import tiledbsoma as soma
from tiledbsoma.experiment_query import X_as_series


# ## Incremental count and mean calculation.
# 
# Many statistics, such as `mean`, are easy to calculate incrementally.  This cell demonstrates a query on the `X['raw']` sparse nD array, which will return results in batches. Accumulate the sum and count incrementally, into `raw_sum` and `raw_n`, and then compute mean.
# 
# First define a query - in this case a slice over the obs axis for cells with a specific tissue & sex value, and all genes on the var axis.  The `query.X()` method returns an iterator of results, each as a PyArrow Table.  Each table will contain the sparse X data and obs/var coordinates, using standard SOMA names:
# 
# * `soma_data` - the X value (float32)
# * `soma_dim_0` - the obs coordinate (int64)
# * `soma_dim_1` - the var coordinate (int64)
# 
# **Important**: the X matrices are joined to var/obs axis DataFrames by an integer join "id" (aka `soma_joinid`). They are *NOT* positionally indexed, and any given cell or gene may have a `soma_joinid` of any value (e.g., a large integer). In other words, for any given `X` value, the `soma_dim_0` corresponds to the `soma_joinid` in the `obs` dataframe, and the `soma_dim_1` coordinate corresponds to the `soma_joinid` in the `var` dataframe.
# 
# For convenience, the query package contains a utility function to simplify operations on query slices.  `query.indexer` returns an indexer that can be used to wrap the output of `query.X()`, converting from `soma_joinids` to positional indexing. Positions are `[0, N)`, where `N` are the number of results on the query for any given axis (equivalent to the Pandas `.iloc` of the axis dataframe).
# 
# Key points:
# 
# * it is expensive to query and read the results - so rather than make multiple passes over the data, read it once and perform multiple computations.
# * by default, data in the census is indexed by `soma_joinid` and not positionally. Use `query.indexer` if you want positions.

# In[2]:


with cellxgene_census.open_soma() as census:
    mouse = census["census_data"]["mus_musculus"]
    with mouse.axis_query(
        measurement_name="RNA",
        obs_query=soma.AxisQuery(value_filter="tissue=='brain' and sex=='male' and is_primary_data==True"),
    ) as query:
        var_df = query.var().concat().to_pandas().set_index("soma_joinid")
        n_vars = len(var_df)

        raw_n = np.zeros((n_vars,), dtype=np.int64)  # accumulate number of non-zero X values
        raw_sum = np.zeros((n_vars,), dtype=np.float64)  # accumulate the sum of expression

        # query.X() returns an iterator of pyarrow.Table, with X data in COO format.
        # You can request an indexer from the query that will map it to positional indices
        indexer = query.indexer
        for arrow_tbl in query.X("raw").tables():
            var_dim = indexer.by_var(arrow_tbl["soma_dim_1"])
            data = arrow_tbl["soma_data"]
            np.add.at(raw_n, var_dim, 1)
            np.add.at(raw_sum, var_dim, data)

    with np.errstate(divide="ignore", invalid="ignore"):
        raw_mean = raw_sum / query.n_obs
    raw_mean[np.isnan(raw_mean)] = 0

    var_df = var_df.assign(raw_n=pd.Series(data=raw_n, index=var_df.index))
    var_df = var_df.assign(raw_mean=pd.Series(data=raw_mean, index=var_df.index))

    display(var_df)


# ## Incremental variance calculation
# 
# Other statistics are not as simple when implemented as an online algorithm. This cell demonstrates an implementation of an online computation of `variance`, using [Welford's online calculation of mean and variance](https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Welford's_online_algorithm).

# In[3]:


import numba
import numpy.typing as npt


class OnlineMatrixMeanVariance:
    n_samples: int
    n_variables: int

    def __init__(self, n_samples: int, n_variables: int):
        """Compute mean and variance for n_variables over n_samples, encoded
        in a COO format. Equivalent to:
            numpy.mean(data, axis=0)
            numpy.var(data, axix=0)
        where the input `data` is of shape (n_samples, n_variables)
        """
        self.n_samples = n_samples
        self.n_variables = n_variables

        self.n_a = np.zeros((n_variables,), dtype=np.int32)
        self.u_a = np.zeros((n_variables,), dtype=np.float64)
        self.M2_a = np.zeros((n_variables,), dtype=np.float64)

    def update(self, coord_vec: npt.NDArray[np.int64], value_vec: npt.NDArray[np.float32]) -> None:
        _mean_variance_update(coord_vec, value_vec, self.n_a, self.u_a, self.M2_a)

    def finalize(self) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
        """Returns tuple containing mean and variance"""
        u, M2 = _mean_variance_finalize(self.n_samples, self.n_a, self.u_a, self.M2_a)

        # compute sample variance
        var = M2 / max(1, (self.n_samples - 1))

        return u, var


@numba.jit(nopython=True)
def _mean_variance_update(
    col_arr: npt.NDArray[np.int64],
    val_arr: npt.NDArray[np.float32],
    n: npt.NDArray[np.int32],
    u: npt.NDArray[np.float64],
    M2: npt.NDArray[np.float64],
):
    """Incrementally accumulate mean and sum of square of distance from mean using
    Welford's online method.
    """
    for col, val in zip(col_arr, val_arr):
        u_prev = u[col]
        M2_prev = M2[col]
        n[col] += 1
        u[col] = u_prev + (val - u_prev) / n[col]
        M2[col] = M2_prev + (val - u_prev) * (val - u[col])


@numba.jit(nopython=True)
def _mean_variance_finalize(
    n_samples: int,
    n_a: npt.NDArray[np.int32],
    u_a: npt.NDArray[np.float64],
    M2_a: npt.NDArray[np.float64],
):
    """Finalize incremental values, acconting for missing elements (due to sparse input).
    Non-sparse and sparse combined using Chan's parallel adaptation of Welford's.
    The code assumes the sparse elements are all zero and ignores those terms.
    """
    n_b = n_samples - n_a
    delta = -u_a  # assumes u_b == 0
    u = (n_a * u_a) / n_samples
    M2 = M2_a + delta**2 * n_a * n_b / n_samples  # assumes M2_b == 0
    return u, M2


with cellxgene_census.open_soma() as census:
    mouse = census["census_data"]["mus_musculus"]
    with mouse.axis_query(
        measurement_name="RNA",
        obs_query=soma.AxisQuery(value_filter="tissue=='brain' and sex=='male' and is_primary_data==True"),
    ) as query:
        var_df = query.var().concat().to_pandas().set_index("soma_joinid")
        n_vars = len(var_df)

        indexer = query.indexer
        mvn = OnlineMatrixMeanVariance(query.n_obs, n_vars)
        for arrow_tbl in query.X("raw").tables():
            var_dim = indexer.by_var(arrow_tbl["soma_dim_1"])
            data = arrow_tbl["soma_data"].to_numpy()
            mvn.update(var_dim, data)

        u, v = mvn.finalize()

    var_df = var_df.assign(raw_mean=pd.Series(data=u, index=var_df.index))
    var_df = var_df.assign(raw_variance=pd.Series(data=v, index=var_df.index))

    display(var_df)


# ## Counting cells per gene, grouped by `dataset_id`
# 
# This example demonstrates a more complex example where the goal is to count the number of cells per gene, grouped by cell dataset_id.  The result is a Pandas DataFrame indexed by `obs.dataset_id` and `var.feature_id`, containing the number of cells per pair.
# 
# This example does not use positional indexing, but rather demonstrates the use of Pandas DataFrame `join` to join on the `soma_joinid`. For the sake of this example we will query only 4 genes, but this can be expanded to all genes.

# In[4]:


with cellxgene_census.open_soma() as census:
    mouse = census["census_data"]["mus_musculus"]

    with mouse.axis_query(
        measurement_name="RNA",
        obs_query=soma.AxisQuery(value_filter="tissue=='brain'"),
        var_query=soma.AxisQuery(value_filter="feature_name in ['Malat1', 'Ptprd', 'Dlg2', 'Pcdh9']"),
    ) as query:
        obs_df = query.obs(column_names=["soma_joinid", "dataset_id"]).concat().to_pandas().set_index("soma_joinid")
        var_df = query.var().concat().to_pandas().set_index("soma_joinid")
        n_cells_by_dataset = pd.Series(
            0,
            index=pd.MultiIndex.from_product(
                (var_df.index, obs_df.dataset_id.unique()),
                names=["soma_joinid", "dataset_id"],
            ),
            dtype=np.int64,
            name="n_cells",
        )

        for X_tbl in query.X("raw").tables():
            # Group by dataset_id and count unique (genes, dataset_id)
            value_counts = (
                X_as_series(X_tbl)
                .to_frame()
                .join(obs_df[["dataset_id"]], on="soma_dim_0")
                .reset_index(level=1)
                .drop(columns=["soma_data"])
                .value_counts()
            )
            np.add.at(
                n_cells_by_dataset,
                n_cells_by_dataset.index.get_indexer(value_counts.index),
                value_counts.to_numpy(),
            )

    # drop any combinations that are not observed
    n_cells_by_dataset = n_cells_by_dataset[n_cells_by_dataset > 0]

    # and join with var_df to pick up feature_id and feature_name
    n_cells_by_dataset = (
        n_cells_by_dataset.to_frame()
        .reset_index(level=1)
        .join(var_df[["feature_id", "feature_name"]])
        .set_index(["dataset_id", "feature_id"])
    )

    display(n_cells_by_dataset)




# Section: api_demo-census_embedding_search

#!/usr/bin/env python
# coding: utf-8

# # Find the most similar Census cells in any other data with embeddings vector search

# This tutorial demonstrates the experimental API in the Python `cellxgene_census` package to search Census embeddings using [TileDB-Vector-Search](https://github.com/TileDB-Inc/TileDB-Vector-Search) indexes. We will generate [scVI embeddings](https://docs.scvi-tools.org/en/1.0.0/tutorials/notebooks/api_overview.html) for some test cells, search the Census scVI embeddings for nearest neighbors, and use them to predict cell type and tissue of the test cells.
# 
# To reproduce this notebook, `pip install 'cellxgene_census[experimental]' scvi-tools` and obtain the [pbmc3k blood cells data](https://scanpy-tutorials.readthedocs.io/en/latest/pbmc3k.html) and scVI model:

# **Contents**
# 
# 1. Downloading data and Census scVI model.
# 2. Loading and embedding pbmc3k cells.
# 3. Search for similar Census cells.
# 4. Predicting cell metadata.
# 
# ⚠️ Note that the Census RNA data includes duplicate cells present across multiple datasets. Duplicate cells can be filtered in or out using the cell metadata variable `is_primary_data` which is described in the [Census schema](https://github.com/chanzuckerberg/cellxgene-census/blob/main/docs/cellxgene_census_schema.md#repeated-data).

# ## Downloading data and Census scVI model
# 
# Download the data first. This the 10X PBMC 3K dataset.

# In[1]:


import warnings

warnings.filterwarnings("ignore")

import anndata
import cellxgene_census
import cellxgene_census.experimental
import pandas as pd
import scanpy as sc
import scvi

CENSUS_VERSION = "2024-07-01"


# Download the test data which will be used to look for their Census most similar cells. This is the test 10X PMBC dataset.

# In[2]:


get_ipython().system('mkdir -p data')
get_ipython().system('wget --no-check-certificate -q -O data/pbmc3k_filtered_gene_bc_matrices.tar.gz http://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz')
get_ipython().system('tar -xzf data/pbmc3k_filtered_gene_bc_matrices.tar.gz -C data/')


# Now download the model corresponding to the census version used in the notebook. 
# 
# First find S3 location of the model.

# In[3]:


with cellxgene_census.open_soma(census_version=CENSUS_VERSION) as census:
    census = cellxgene_census.open_soma(census_version=CENSUS_VERSION)

    scvi_info = cellxgene_census.experimental.get_embedding_metadata_by_name(
        embedding_name="scvi",
        organism="homo_sapiens",
        census_version=CENSUS_VERSION,
    )

scvi_info["model_link"]


# Now use that path to download via HTTPs, make sure to update the URL if you change the Census version.

# In[4]:


get_ipython().system('mkdir -p scvi-human-2024-07-01')
get_ipython().system('wget --no-check-certificate -q -O scvi-human-2024-07-01/model.pt https://cellxgene-contrib-public.s3.us-west-2.amazonaws.com/models/scvi/2024-07-01/homo_sapiens/model.pt')


# ## Loading and embedding pbmc3k cells

# Load the pmbc3k cell data into an AnnData:

# In[5]:


adata = sc.read_10x_mtx("data/filtered_gene_bc_matrices/hg19/", var_names="gene_ids")
adata.var["ensembl_id"] = adata.var.index
adata.obs["n_counts"] = adata.X.sum(axis=1)
adata.obs["joinid"] = list(range(adata.n_obs))
adata.obs["batch"] = "unassigned"


# Run them through the scVI forward pass and extract their latent representation (embedding):

# In[6]:


scvi.model.SCVI.prepare_query_anndata(adata, "scvi-human-2024-07-01")
vae_q = scvi.model.SCVI.load_query_data(
    adata,
    "scvi-human-2024-07-01",
)

# This allows for a simple forward pass
vae_q.is_trained = True
latent = vae_q.get_latent_representation()
adata.obsm["scvi"] = latent


# The scVI embedding vectors for each cell are now stored in the `scvi` obsm layer. Lastly, clean up the AnnData a little:

# In[7]:


# filter out missing features
adata = adata[:, adata.var["gene_symbols"].notnull().values].copy()
adata.var.set_index("gene_symbols", inplace=True)
# assign placeholder cell_type and tissue_general labels
adata.var_names = adata.var["ensembl_id"]
adata.obs["cell_type"] = "Query - PBMC 10X"
adata.obs["tissue_general"] = "Query - PBMC 10X"


# And for visualization, compute leiden clusters on the scVI embeddings. These will come at handy later on.

# In[8]:


sc.pp.neighbors(adata, n_neighbors=15, use_rep="scvi")
sc.tl.umap(adata)
sc.tl.leiden(adata)
sc.pl.umap(adata, color="leiden")


# ## Search for similar Census cells
# 
# Use the CELLxGENE Census experimental API to search the vector index of scVI embeddings.

# In[9]:


get_ipython().run_cell_magic('time', '', 'neighbors = cellxgene_census.experimental.find_nearest_obs(\n    "scvi", "homo_sapiens", CENSUS_VERSION, query=adata, k=30, memory_GiB=8, nprobe=20\n)\n')


# This accessed the cell embeddings in the `scvi` obsm layer and searched the latent space for the *k* nearest neighbors (by Euclidean distance) among the Census cell embeddings, returning the distances and obs `soma_joinid`s (*k* for each query cell).

# In[10]:


neighbors


# To explore these results, fetch an AnnData with each query cell's (single) nearest neighbor in scVI's latent space, including their embedding vectors:

# In[11]:


with cellxgene_census.open_soma(census_version=CENSUS_VERSION) as census:
    neighbors_adata = cellxgene_census.get_anndata(
        census,
        "homo_sapiens",
        "RNA",
        obs_coords=sorted(neighbors.neighbor_ids[:, 0].tolist()),
        obs_embeddings=["scvi"],
        X_name="normalized",
        column_names={"obs": ["soma_joinid", "tissue", "tissue_general", "cell_type"]},
    )
neighbors_adata.var_names = neighbors_adata.var["feature_id"]


# Make a UMAP visualization of these nearest neighbors:

# In[12]:


sc.pp.neighbors(neighbors_adata, n_neighbors=15, use_rep="scvi")
sc.tl.umap(neighbors_adata)
sc.pl.umap(neighbors_adata, color="tissue_general")


# As expected, the nearest neighbors are largely Census blood cells, with some distinct cell type clusters.
# 
# Now display the pbmc3k query cells together with their Census nearest neighbors:

# In[13]:


adata_concat = anndata.concat([adata, neighbors_adata])
sc.pp.neighbors(adata_concat, n_neighbors=15, use_rep="scvi")
sc.tl.umap(adata_concat)
sc.pl.umap(adata_concat, color=["tissue_general"])


# The neighbor clusters appear to correspond to cell type heterogeneity also present in the query.
# 
# ## Predicting cell metadata
# 
# The experimental API also has a method to predict metadata attributes of the query cells, like `tissue_general` and `cell_type`, based on the Census nearest neighbors.

# In[14]:


predictions = cellxgene_census.experimental.predict_obs_metadata(
    "homo_sapiens", CENSUS_VERSION, neighbors, ["tissue_general", "cell_type"]
)
predictions


# Unlike the above visualizations of *single* nearest neighbors, these predictions are informed by each query cell's *k*=30 neighbors, with an attached confidence score.
# 
# These predictions can be added to the original AnnData object for visualization.

# In[15]:


predictions.index = adata.obs.index
predictions = predictions.rename(columns={"cell_type": "predicted_cell_type"})
adata.obs = pd.concat([adata.obs, predictions], axis=1)


# In[16]:


sc.pl.umap(adata, color="predicted_cell_type")


# Now the leiden clusters calculated at the beginning can be annotated by popular vote, whereby each cluster gets assigned the most common predicted cell type from the previous step.

# In[17]:


adata.obs["predicted_consolidated_cell_type"] = ""
for leiden_cluster in adata.obs["leiden"].drop_duplicates():
    most_popular_type = (
        adata.obs.loc[adata.obs["leiden"] == leiden_cluster,].value_counts("predicted_cell_type").index[0]
    )
    adata.obs.loc[adata.obs["leiden"] == leiden_cluster, "predicted_consolidated_cell_type"] = most_popular_type


# In[18]:


sc.pl.umap(adata, color="predicted_consolidated_cell_type")




# Section: analysis_demo-comp_bio_embedding_exploration

#!/usr/bin/env python
# coding: utf-8

# # Exploring biologically relevant clusters in Census embeddings
# 
# In this notebook, we explore biologically relevant clusters in Census embeddings using UMAP as a visualization tool. This demonstration assumes knowledge of how to access both, collaboration and hosted (community) Census embeddings. To learn the basics on accessing these data please visit the [Census model page](https://cellxgene.cziscience.com/census-models).
# 
# **IMPORTANT:** This tutorial requires cellxgene-census package version 1.9.1 or later.
# 
# **Contents**
# 
# 1. Background
# 1. Requirements
# 1. Imports and function definitions
# 1. Melanocytes in eye
# 1. Retinal bipolar neurons in eye
# 1. Dopaminergic neurons in brain
# 1. Pulmonary ionocytes in lung
# 
# >⚠️ Note that the Census RNA data includes duplicate cells present across multiple datasets. Duplicate cells can be filtered in or out using the cell metadata variable is_primary_data which is described in the [Census schema](https://github.com/chanzuckerberg/cellxgene-census/blob/main/docs/cellxgene_census_schema.md#repeated-data).
# 
# ## Background
# 
# The journey from a gene expression matrix to a 2D scatterplot involves numerous highly nonlinear transformations. Such transformations can introduce artifacts that affect both the global and local structures in the visualized manifold. 
# 
# Common issues like overclustering and clustering by batch are typical artifacts resulting from these dimensionality reduction methods. With that in mind, these embeddings and their UMAP visualizations are best used as tools for generating hypotheses. They should not be the final word in analysis. Instead, we recommend focusing on the full representation of the embedding matrices and ultimately returning to the underlying gene expressions to investigate the reasons behind the observed clustering patterns.
# 
# One of the key objectives of foundation models in single-cell RNA sequencing is to embed cells within a universal coordinate system that minimizes the impact of technical variations, such as batch effects. However, as we will see in the examples presented in this notebook, cells often cluster by batch. This clustering could be biologically driven, as certain cell types or states might be unique to specific batches (e.g. datasets). In other cases, the separation might be purely due to systematic technical biases or a combination of both biological and technical factors. Complicating the matter, the techniques used for nearest neighbor graph construction and 2D projection can themselves amplify batch effects. Rigorous benchmarking is necessary to fully assess each model's capability in integrating data within their respective latent spaces. This complexity highlights that data integration in single-cell RNA sequencing remains a challenging and unsolved problem. 
# 
# In this tutorial, we briefly highlight a few simple case studies that illustrate the capacity of these embeddings to capture intriguing biological phenomena.
# 
# **Disclaimers** 
# 
# 1. These embeddings were explored in-depth in a [cellxgene](https://github.com/chanzuckerberg/cellxgene) instance and not all the insights gleaned there will be expanded on here.
# 2. Most of the following examples utilize UMAP to visualize embeddings in a 2D scatter plot, however as shown [here](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1011288) and [here](https://www.cell.com/cell-systems/abstract/S2405-4712%2823%2900209-0), biological interpretations from these visualizations may be inaccurate.
# 
# 
# 
# ## Requirements
#  - cellxgene-census
#  - scanpy
#  - numpy
#  - scipy
#  - leidenalg
#  - hdbscan
#  - pandas
#  - scikit-learn

# ## Imports and function definitions

# In[1]:


import warnings
from typing import List

import anndata
import cellxgene_census
import numpy as np
import scanpy as sc

warnings.filterwarnings("ignore")


def remove_missing_embedding_cells(adata: anndata.AnnData, emb_names: List[str]):
    """Embeddings with missing data contain all NaN,
    so we must find the intersection of non-NaN rows in the fetched embeddings
    and subset the AnnData accordingly.
    """
    filt = np.ones(adata.shape[0], dtype="bool")
    for key in emb_names:
        nan_row_sums = np.sum(np.isnan(adata.obsm[key]), axis=1)
        total_columns = adata.obsm[key].shape[1]
        filt = filt & (nan_row_sums != total_columns)
    adata = adata[filt].copy()

    return adata


def generate_umaps_from_embeddings(adata: anndata.AnnData, emb_names: list, metric="euclidean"):
    """Generate UMAPs from embeddings stored in `adata.obsm`.
    `emb_names` is a list that contains keys present in `adata.obsm`.
    """
    adata = adata.copy()
    for emb_name in emb_names:
        print(f"Generating UMAP for {emb_name}")
        sc.pp.neighbors(
            adata,
            n_neighbors=15,
            use_rep=emb_name,
            method="umap",
            key_added=emb_name,
            metric=metric,
        )
        sc.tl.umap(adata, neighbors_key=emb_name)
        X_emb_name = emb_name if emb_name[:2] == "X_" else f"X_{emb_name}"
        if metric != "euclidean":
            X_emb_name += f"_{metric}"
        adata.obsm[f"{X_emb_name}_umap"] = adata.obsm["X_umap"]
        del adata.obsm["X_umap"]

    adata.var_names = adata.var["feature_name"]
    adata.raw = adata.copy()
    sc.pp.normalize_total(adata, target_sum=10000)
    sc.pp.log1p(adata)

    return adata


# In[2]:


# human embeddings
CENSUS_VERSION = "2023-12-15"
EXPERIMENT_NAME = "homo_sapiens"

# These are embeddings available to this Census version
embedding_names = ["geneformer", "scvi", "scgpt", "uce"]


# ## Melanocytes in eye
# 
# ### Sample and fetch 150k cells from eye tissue

# In[3]:


census = cellxgene_census.open_soma(census_version=CENSUS_VERSION)

# Let's find our cells of interest
obs_value_filter = "tissue_general=='eye' and is_primary_data == True"

obs_df = cellxgene_census.get_obs(census, EXPERIMENT_NAME, value_filter=obs_value_filter, column_names=["soma_joinid"])

print(obs_df.shape[0], "cells in", obs_value_filter)

# Let's subset to 150K
n_subset_cells = 150000

print("Selecting", n_subset_cells, "random cells")
idx_rand = np.random.choice(obs_df.shape[0], size=n_subset_cells, replace=False)
soma_joinids_subset = obs_df["soma_joinid"].values[idx_rand].tolist()


# In[4]:


# Let's get the AnnData
adata = cellxgene_census.get_anndata(
    census=census,
    organism=EXPERIMENT_NAME,
    obs_coords=soma_joinids_subset,
    obs_embeddings=embedding_names,
)

adata = remove_missing_embedding_cells(adata, embedding_names)
adata = generate_umaps_from_embeddings(adata, embedding_names)


# ### Observations
# 
# In the study of melanocytes within the eye, the following observations are made across various embeddings:
# 
#  - Melanocytes are distinctly clustered in all embeddings, with OCA2 as a noted marker.
#  - KIT, [identified as a marker for mature melanocytes](https://iovs.arvojournals.org/article.aspx?articleid=2743921), shows varying degrees of separation. In SCVI and UCE embeddings, mature and immature melanocytes are clearly separable based on KIT expression. The scGPT embedding shows a slight extension from the main melanocyte cluster, indicative of some separation, where KIT expression is concentrated. In Geneformer, cells expressing KIT are primarily found at one end of the larger melanocyte cluster.
#  - The UCE embedding demonstrates potential signs of overclustering. An example is seen in retinal bipolar neurons, which separate into numerous small satellite clusters without clear gene expression signatures. The high degree of local structure in the UCE manifold is probably due to the presence of many disconnected components in the graph constructed by UMAP. This could also indicate that the embedding may capture less global structure.
#  - Across all embeddings, assays tend to cluster separately. The extent to which this reflects biological variability versus technical variation is unclear.
#  - Qualitatively, SCVI appears to offer the best integration across different datasets, with other embeddings showing more pronounced clustering by dataset.
#  

# In[5]:


sc.pl.scatter(
    adata,
    basis="geneformer_umap",
    color=["OCA2", "KIT", "cell_type"],
    size=10,
    use_raw=False,
)
sc.pl.scatter(
    adata,
    basis="scgpt_umap",
    color=["OCA2", "KIT", "cell_type"],
    size=10,
    use_raw=False,
)
sc.pl.scatter(adata, basis="uce_umap", color=["OCA2", "KIT", "cell_type"], size=10, use_raw=False)
sc.pl.scatter(adata, basis="scvi_umap", color=["OCA2", "KIT", "cell_type"], size=10, use_raw=False)


# In[6]:


sc.pl.scatter(
    adata,
    basis="geneformer_umap",
    color=["dataset_id", "assay"],
    size=10,
    use_raw=False,
)
sc.pl.scatter(adata, basis="scgpt_umap", color=["dataset_id", "assay"], size=10, use_raw=False)
sc.pl.scatter(adata, basis="uce_umap", color=["dataset_id", "assay"], size=10, use_raw=False)
sc.pl.scatter(adata, basis="scvi_umap", color=["dataset_id", "assay"], size=10, use_raw=False)


# ## Retinal bipolar neurons in eye
# 
# In a more detailed analysis of retinal bipolar neurons in the eye, we focus on subclustering within this cell type across various embeddings. This involves rerunning UMAP specifically for retinal bipolar neurons and applying Leiden clustering to each embedding. Additionally, we employ [HDBSCAN](https://hdbscan.readthedocs.io/en/latest/how_hdbscan_works.html), a density-based clustering algorithm, on a full pairwise Euclidean distance matrix calculated from each embedding to compare the clustering results.
# 
# Key findings from this analysis include:
# 
#  - In the scGPT embedding, the nearest neighbor graph construction reveals 25 distinct clusters, but the density-based HDBSCAN approach identifies only 3 clusters, indicating a significant difference in clustering patterns.
#  - The UCE and SCVI embeddings show good agreement between graph-based and density-based clustering methods.
#  - For Geneformer, the subclusters observed in UMAP appear to be an artifact of the graph construction, as HDBSCAN results in only one primary cluster.
# 
# When assessing the Normalized Mutual Information (NMI) score between Leiden and HDBSCAN cluster assignments, it's found that:
# 
#  - All embeddings yield Leiden clusters with generally good agreement across methods (NMI > 0.65), indicating a consistent clustering pattern.
#  - HDBSCAN clusterings are more method-specific, reflecting inherent differences in how each embedding interprets distances and densities. Geneformer, in particular, shows minimal agreement with other methods as HDBSCAN only identified one main cluster. This is expected as Geneformer was finetuned on a cell subclass prediction task, which will homogenize cells belonging to the same label.
#  - Additionally, all methods, except SCVI, distinctly separate retinal bipolar neurons by batch (dataset ID), underscoring the presence of batch effects.
# 
# From this analysis, we can draw a couple conclusions:
# 
#  - The construction of k-nearest neighbor graphs in embeddings like UMAP can lead to the identification of subclusters that are not evident when examining pairwise distances directly in the original embedding spaces. This suggests that the reliance of UMAP (and many other methods spanning a wide variety of tasks) on k-nearest neighbor graphs may introduce biologically unjustifiable subclusters. This is a known phenomenon and further supports the recommendation to cross-reference findings using fuller representations of the data.
#  - In this example, UCE clusters the data much more than other methods. This observation holds true in both graph- and density-based clustering. The small clusters identified in UCE often lack unique or biologically relevant gene expression signatures, necessitating further investigation to understand any biological relevance or lack thereof for each of the identified clusters.

# In[7]:


import hdbscan
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from sklearn.metrics import normalized_mutual_info_score

# subset anndata
adata_rbn = adata[adata.obs["cell_type"] == "retinal bipolar neuron"].copy()

# generate UMAPs
adata_rbn = generate_umaps_from_embeddings(adata_rbn, embedding_names)

# run clustering methods
for embedding in embedding_names:
    sc.tl.leiden(adata_rbn, obsp=f"{embedding}_connectivities", key_added=f"{embedding}_leiden")

    points = adata_rbn.obsm[embedding]

    # calculate full pairwise distance matrix
    pairwise_dist = squareform(pdist(points, "euclidean"))
    # run HDBSCAN
    adata_rbn.obs[f"{embedding}_hdbscan"] = (
        hdbscan.HDBSCAN(min_cluster_size=5, min_samples=5, metric="precomputed")
        .fit_predict(pairwise_dist)
        .astype("int")
        .astype("str")
    )

# display UMAPs and report normalized mutual information scores between leiden and hdbscan cluster assignments
for embedding in embedding_names:
    print("Normalized mutual information between Leiden and HDBSCAN clusters:")
    print(normalized_mutual_info_score(adata_rbn.obs[f"{embedding}_leiden"], adata_rbn.obs[f"{embedding}_hdbscan"]))
    sc.pl.scatter(
        adata_rbn,
        basis=f"{embedding}_umap",
        color=[f"{embedding}_leiden", f"{embedding}_hdbscan", "dataset_id"],
        size=10,
        use_raw=False,
    )


# compare leiden and hdbscan cluster assignments across methods and display the similarity tables
embedding_keys = embedding_names
sim_scores_leiden = np.zeros((len(embedding_keys), len(embedding_keys)))
sim_scores_hdbscan = np.zeros((len(embedding_keys), len(embedding_keys)))
for i, embedding_i in enumerate(embedding_keys):
    for j, embedding_j in enumerate(embedding_keys):
        sim_scores_leiden[i, j] = normalized_mutual_info_score(
            adata_rbn.obs[f"{embedding_i}_leiden"],
            adata_rbn.obs[f"{embedding_j}_leiden"],
        )
        sim_scores_hdbscan[i, j] = normalized_mutual_info_score(
            adata_rbn.obs[f"{embedding_i}_hdbscan"],
            adata_rbn.obs[f"{embedding_j}_hdbscan"],
        )

sim_scores_leiden_table = pd.DataFrame(data=sim_scores_leiden, index=embedding_keys, columns=embedding_keys)
sim_scores_hdbscan_table = pd.DataFrame(data=sim_scores_hdbscan, index=embedding_keys, columns=embedding_keys)

print("Leiden:")
print(sim_scores_leiden_table)
print("")
print("HDBSCAN:")
print(sim_scores_hdbscan_table)


# ## Dopaminergic neurons in brain
# ### Sample and fetch 150k cells from brain tissue

# In[8]:


# Let's find our cells of interest
obs_value_filter = "tissue_general=='brain' and is_primary_data == True"

obs_df = cellxgene_census.get_obs(census, EXPERIMENT_NAME, value_filter=obs_value_filter, column_names=["soma_joinid"])

print(obs_df.shape[0], "cells in", obs_value_filter)

# Let's subset to 150K
n_subset_cells = 150000
print("Selecting ", n_subset_cells, " random cells")
idx_rand = np.random.choice(obs_df.shape[0], size=n_subset_cells, replace=False)
soma_joinids_subset = obs_df["soma_joinid"].values[idx_rand].tolist()


# In[9]:


# Let's get the AnnData
adata = cellxgene_census.get_anndata(
    census=census,
    organism=EXPERIMENT_NAME,
    obs_coords=soma_joinids_subset,
    obs_embeddings=embedding_names,
)

adata = remove_missing_embedding_cells(adata, embedding_names)
adata = generate_umaps_from_embeddings(adata, embedding_names)


# ### Observations
# Here, we visualize a randomly selected subset of cells in the brain from CELLxGENE Census. We can observe that dopaminergic neurons, marked by TH expression, separate into distinct clusters in Geneformer and SCVI latent spaces, whereas in UCE and scGPT embeddings, they are grouped at one end of a larger neuron cluster.
# 
# All embeddings show a tendency to cluster by assay, indicating a consistent pattern across different models. Conditions like glioblastoma are clearly separated in all embeddings, while pilocytic astrocytoma is distinctly clustered in Geneformer and UCE and more mixed in others.
# 
# In the UCE embedding, we observe many small satellite glioblastoma clusters outside the main cluster that do not have distinct gene expression signatures. This is similar to what we observed previously in the eye (e.g. for the retinal bipolar neurons).

# In[10]:


sc.pl.scatter(
    adata,
    basis="geneformer_umap",
    color=["TH", "assay", "disease"],
    size=10,
    use_raw=False,
)
sc.pl.scatter(adata, basis="scgpt_umap", color=["TH", "assay", "disease"], size=10, use_raw=False)
sc.pl.scatter(adata, basis="uce_umap", color=["TH", "assay", "disease"], size=10, use_raw=False)
sc.pl.scatter(adata, basis="scvi_umap", color=["TH", "assay", "disease"], size=10, use_raw=False)


# ## Pulmonary ionocytes in lung (Tabula Sapiens)
# ### Fetch lung cells from Tabula Sapiens
# 

# In[11]:


obs_value_filter = "tissue_general=='lung' and dataset_id=='53d208b0-2cfd-4366-9866-c3c6114081bc'"

adata = cellxgene_census.get_anndata(
    census=census,
    organism=EXPERIMENT_NAME,
    obs_value_filter=obs_value_filter,
    obs_embeddings=embedding_names,
)

adata = remove_missing_embedding_cells(adata, embedding_names)
adata = generate_umaps_from_embeddings(adata, embedding_names)


# ### Observations
# For the case study focusing on pulmonary ionocytes in lung tissue, as part of the Tabula Sapiens project, the following observations are noted:
# 
# - In all embeddings, except for SCVI, a clear separation is seen between SmartSeq data and 10x data. The distinction is most pronounced in the scGPT embedding.
# - CFTR, a marker for pulmonary ionocytes, identifies a rare cell type in the lung. This cell type is distinctly recognizable in all the embeddings.

# In[12]:


sc.pl.scatter(
    adata,
    basis="geneformer_umap",
    color=["CFTR", "assay", "cell_type"],
    size=10,
    use_raw=False,
)
sc.pl.scatter(
    adata,
    basis="scgpt_umap",
    color=["CFTR", "assay", "cell_type"],
    size=10,
    use_raw=False,
)
sc.pl.scatter(
    adata,
    basis="uce_umap",
    color=["CFTR", "assay", "cell_type"],
    size=10,
    use_raw=False,
)
sc.pl.scatter(
    adata,
    basis="scvi_umap",
    color=["CFTR", "assay", "cell_type"],
    size=10,
    use_raw=False,
)


# In[13]:


census.close()




# Section: api_demo-census_gget_demo

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



