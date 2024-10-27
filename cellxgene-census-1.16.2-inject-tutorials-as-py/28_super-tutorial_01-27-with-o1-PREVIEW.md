---

# CZ CELLxGENE Discover Census Super-Tutorial

Welcome to the comprehensive guide for using the **latest CELLxGENE Census Python API (v1.16.2)**. This tutorial will help you access, query, and analyze single-cell RNA sequencing data from the CZ CELLxGENE Discover Census efficiently. It is designed to provide you with clear, coherent instructions and examples to accelerate your research.

---

## Contents

1. [Introduction](#introduction)
2. [Installation](#installation)
3. [Quick Start](#quick-start)
4. [Understanding the Census Data and Schema](#understanding-the-census-data-and-schema)
5. [Accessing and Querying Census Data](#accessing-and-querying-census-data)
6. [Querying and Fetching Single-Cell Data and Metadata](#querying-and-fetching-single-cell-data-and-metadata)
7. [Memory-Efficient Computations](#memory-efficient-computations)
8. [Working with Pre-calculated Embeddings and Models](#working-with-pre-calculated-embeddings-and-models)
9. [Accessing CELLxGENE-Hosted Embeddings](#accessing-cellxgene-hosted-embeddings)
10. [Finding the Most Similar Census Cells Using Embeddings Vector Search](#finding-the-most-similar-census-cells-using-embeddings-vector-search)
11. [Using the scVI Pretrained Model for Data Projection and Cell Type Prediction](#using-the-scvi-pretrained-model-for-data-projection-and-cell-type-prediction)
12. [Using the Geneformer Fine-Tuned Model for Cell Class Prediction and Data Projection](#using-the-geneformer-fine-tuned-model-for-cell-class-prediction-and-data-projection)
13. [Using Census with PyTorch for Machine Learning](#using-census-with-pytorch-for-machine-learning)
14. [Benchmarks of Census Models](#benchmarks-of-census-models)
15. [Accessing Census Data in AWS](#accessing-census-data-in-aws)
16. [Querying Data Using the gget cellxgene Module](#querying-data-using-the-gget-cellxgene-module)
17. [Additional Resources](#additional-resources)
18. [Citing Census](#citing-census)
19. [Feedback and Support](#feedback-and-support)

---

## Introduction

The **CZ CELLxGENE Discover Census** provides efficient computational tooling to **access, query, and analyze all single-cell RNA data from CZ CELLxGENE Discover**. By leveraging cell-based slicing and querying, you can interact with the data through TileDB-SOMA or obtain slices in AnnData, Seurat, or SingleCellExperiment objects, significantly minimizing data harmonization efforts.

---

## Installation

### Requirements

- **Operating System**: Linux or macOS
- **Python Version**: 3.10 to 3.12
- **Recommended**:
  - At least 16 GB of memory
  - Internet connection with >5 Mbps bandwidth
  - For increased performance, use the API through an AWS EC2 instance in the `us-west-2` region (where Census data is hosted)

### Installation Steps

1. **Create and Activate a Virtual Environment (Optional)**

   ```shell
   python -m venv ./venv
   source ./venv/bin/activate
   ```

2. **Install the `cellxgene-census` Package via pip**

   ```shell
   pip install -U cellxgene-census
   ```

   To include experimental features like PyTorch loaders and model support, install with:

   ```shell
   pip install -U 'cellxgene-census[experimental]'
   ```

   **Note**: If installing in a Databricks notebook environment, use `%pip install`. Do not use `%sh pip install`.

---

## Quick Start

This section provides quick examples of common operations using the Census API in Python. For detailed tutorials, refer to the [Additional Resources](#additional-resources) section.

### Import the Census API

```python
import cellxgene_census
```

### Access Help and Documentation

```python
help(cellxgene_census)
help(cellxgene_census.get_anndata)
```

### Example 1: Querying Cell Metadata

Retrieve metadata for female microglial cells and neurons, selecting specific columns.

```python
import cellxgene_census

with cellxgene_census.open_soma() as census:
    cell_metadata = census["census_data"]["homo_sapiens"].obs.read(
        value_filter="sex == 'female' and cell_type in ['microglial cell', 'neuron']",
        column_names=["assay", "cell_type", "tissue", "tissue_general", "suspension_type", "disease"]
    )
    cell_metadata = cell_metadata.concat().to_pandas()
    print(cell_metadata)
```

**Output**:

```
          assay        cell_type         tissue tissue_general suspension_type disease     sex
0     10x 3' v3  microglial cell            eye            eye            cell  normal  female
1     10x 3' v3  microglial cell            eye            eye            cell  normal  female
...         ...              ...            ...            ...             ...     ...     ...
```

### Example 2: Obtaining a Slice as AnnData

Create an `AnnData` object for specific genes and cells with selected metadata.

```python
import cellxgene_census

with cellxgene_census.open_soma() as census:
    adata = cellxgene_census.get_anndata(
        census=census,
        organism="Homo sapiens",
        var_value_filter="feature_id in ['ENSG00000161798', 'ENSG00000188229']",
        obs_value_filter="sex == 'female' and cell_type in ['microglial cell', 'neuron']",
        obs_column_names=["assay", "cell_type", "tissue", "tissue_general", "suspension_type", "disease"],
    )
    print(adata)
```

**Output**:

```
AnnData object with n_obs × n_vars = 379224 × 2
    obs: 'assay', 'cell_type', 'tissue', 'tissue_general', 'suspension_type', 'disease', 'sex'
    var: 'soma_joinid', 'feature_id', 'feature_name', 'feature_length'
```

### Example 3: Memory-Efficient Queries

Access data for larger-than-memory operations using TileDB-SOMA.

```python
import cellxgene_census
import tiledbsoma

with cellxgene_census.open_soma() as census:
    human = census["census_data"]["homo_sapiens"]
    query = human.axis_query(
        measurement_name="RNA",
        obs_query=tiledbsoma.AxisQuery(value_filter="tissue == 'brain' and sex == 'male'")
    )

    iterator = query.X("raw").tables()
    for raw_slice in iterator:
        # Perform operations on raw_slice
        pass

    query.close()
```

---

## Understanding the Census Data and Schema

### Overview

The Census is a collection of single-cell RNA data organized using the [SOMA](https://github.com/single-cell-data/SOMA) data model and stored using [TileDB-SOMA](https://github.com/single-cell-data/TileDB-SOMA).

### Data Included in the Census

- **Species**: Human (*Homo sapiens*) and Mouse (*Mus musculus*)
- **Data Types**:
  - Full-gene sequencing read counts (e.g., Smart-Seq2)
  - Molecule counts (e.g., 10X)
- **Duplicate Cells**: Duplicate cells are present across multiple datasets. Use the `is_primary_data` column in cell metadata to filter duplicates.

### Data Organization

- **Cell Metadata (`obs`)**: Contains annotations for each cell.
- **Gene Metadata (`var`)**: Contains information about each gene.
- **Expression Matrices (`X`)**:
  - **Raw Counts**: `X["raw"]`
  - **Normalized Counts**: `X["normalized"]` (available in Census schema V1.1.0 and above)
- **Feature Presence Matrix**: Indicates which genes were measured in each dataset.

### Accessing the Census Data and Schema

For a detailed description of the Census schema, please refer to the [Census Data and Schema](https://github.com/chanzuckerberg/cellxgene-census/blob/main/docs/cellxgene_census_schema.md) documentation.

---

## Accessing and Querying Census Data

### Opening the Census

```python
import cellxgene_census

census = cellxgene_census.open_soma()
```

### Closing the Census

```python
census.close()
```

### Retrieving Cell Metadata

```python
cell_metadata = census["census_data"]["homo_sapiens"].obs.read()
```

### Retrieving Gene Metadata

```python
gene_metadata = census["census_data"]["homo_sapiens"].ms["RNA"].var.read().concat().to_pandas()
```

### Filtering Data Using Value Filters

You can filter data using expressions similar to SQL `WHERE` clauses.

**Example**:

```python
value_filter = "tissue == 'lung' and cell_type == 'epithelial cell'"
```

---

## Querying and Fetching Single-Cell Data and Metadata

This section demonstrates how to query the expression data and cell/gene metadata from the Census and load them into common in-memory Python objects, including `pandas.DataFrame` and `anndata.AnnData`.

**Note**: The Census RNA data includes duplicate cells present across multiple datasets. Duplicate cells can be filtered in or out using the cell metadata variable `is_primary_data` as described in the [Census schema](https://github.com/chanzuckerberg/cellxgene-census/blob/main/docs/cellxgene_census_schema.md#repeated-data).

### Importing the Census API

```python
import cellxgene_census

census = cellxgene_census.open_soma()
```

### Querying Expression Data

Use the `get_anndata` method to query and fetch expression data. This method allows you to specify filters and columns for both cell (`obs`) and gene (`var`) metadata.

- **Arguments**:
  - `obs_column_names` and `var_column_names`: Lists of strings indicating the columns to select for cell and gene metadata, respectively.
  - `obs_value_filter`: Expression to filter cells.
  - `var_value_filter`: Expression to filter genes.

**Example**: Fetch expression data for specific genes and B cells in the lung with COVID-19, selecting the `sex` metadata.

```python
adata = cellxgene_census.get_anndata(
    census=census,
    organism="Homo sapiens",
    var_value_filter="feature_id in ['ENSG00000161798', 'ENSG00000188229']",
    obs_value_filter="cell_type == 'B cell' and tissue_general == 'lung' and disease == 'COVID-19' and is_primary_data == True",
    obs_column_names=["sex"],
)
```

Inspect the resulting `AnnData` object:

```python
print(adata)
```

**Output**:

```
AnnData object with n_obs × n_vars = 12345 × 2
    obs: 'sex'
    var: 'soma_joinid', 'feature_id', 'feature_name', 'feature_length'
```

### Querying Cell Metadata (`obs`)

The cell metadata for humans is located at `census["census_data"]["homo_sapiens"].obs`. You can use the `read()` method with `value_filter` and `column_names` to fetch specific cells.

**List Available Columns**:

```python
keys = list(census["census_data"]["homo_sapiens"].obs.keys())
print(keys)
```

**Example**: Fetch all cell metadata where `sex` is `'unknown'`.

```python
cell_metadata_unknown_sex = cellxgene_census.get_obs(
    census,
    "homo_sapiens",
    value_filter="sex == 'unknown'"
)
```

**Example**: Fetch the `disease` column for B cells in the lung from non-duplicated cells.

```python
cell_metadata_b_cell = cellxgene_census.get_obs(
    census,
    "homo_sapiens",
    value_filter="cell_type == 'B cell' and tissue_general == 'lung' and is_primary_data == True",
    column_names=["disease"],
)
```

### Querying Gene Metadata (`var`)

The gene metadata for humans is located at `census["census_data"]["homo_sapiens"].ms["RNA"].var`.

**List Available Columns**:

```python
keys = list(census["census_data"]["homo_sapiens"].ms["RNA"].var.keys())
print(keys)
```

**Example**: Get `feature_name` and `feature_length` for specific genes.

```python
gene_metadata = cellxgene_census.get_var(
    census,
    "homo_sapiens",
    value_filter="feature_id in ['ENSG00000161798', 'ENSG00000188229']",
    column_names=["feature_name", "feature_length"],
)
```

### Closing the Census

```python
census.close()
```

---

## Memory-Efficient Computations

### Calculating Mean and Variance

Efficiently calculate average and variance gene expression across millions of cells.

**Example**:

```python
import cellxgene_census
import tiledbsoma as soma
from cellxgene_census.experimental.pp import mean_variance

census = cellxgene_census.open_soma()
human_data = census["census_data"]["homo_sapiens"]

cell_filter = (
    "is_primary_data == True and tissue_general == 'lung' and cell_type == 'epithelial cell'"
)
gene_filter = "feature_name in ['KRAS', 'AQP4']"

query = human_data.axis_query(
    measurement_name="RNA",
    obs_query=soma.AxisQuery(value_filter=cell_filter),
    var_query=soma.AxisQuery(value_filter=gene_filter)
)

mean_variance_df = mean_variance(query, axis=0, calculate_mean=True, calculate_variance=True)
gene_df = query.var().concat().to_pandas()

query.close()
census.close()
```

### Identifying Highly Variable Genes

Obtain highly variable genes accounting for batch effects.

**Example**:

```python
import cellxgene_census
from cellxgene_census.experimental.pp import get_highly_variable_genes

census = cellxgene_census.open_soma()

hvg = get_highly_variable_genes(
    census,
    organism="Homo sapiens",
    obs_value_filter="is_primary_data == True and tissue_general == 'esophagus'",
    n_top_genes=1000,
    batch_key="dataset_id"
)

census.close()
```

---

## Working with Pre-calculated Embeddings and Models

### Accessing Pre-calculated Embeddings

Retrieve pre-calculated embeddings and use them in your workflows.

**Example**:

```python
import cellxgene_census

embedding_names = ["scvi", "geneformer", "scgpt", "uce"]

with cellxgene_census.open_soma(census_version="2023-12-15") as census:
    adata = cellxgene_census.get_anndata(
        census,
        organism="homo_sapiens",
        measurement_name="RNA",
        obs_value_filter="tissue_general == 'central nervous system'",
        obs_embeddings=embedding_names
    )
```

### Using the Census Trained Models

- **scVI**: Latent spaces from a model trained on all Census data.
- **Geneformer**: Fine-tuned embeddings and models for cell classification.
- **scGPT**: Generative pre-trained transformer embeddings.
- **Universal Cell Embeddings (UCE)**: Pre-computed embeddings available in the Census.

---

## Accessing CELLxGENE-Hosted Embeddings

### Overview

This section demonstrates how to access the **CELLxGENE-hosted embeddings** from the CZ CELLxGENE Discover Census. These embeddings have been contributed by the community and are not actively maintained by CELLxGENE Discover. Find out more about these on the [Census Models page](https://cellxgene.cziscience.com/census-models).

**Note**: This tutorial requires `cellxgene-census` package version **1.9.1** or later.

### Quick Start

The easiest way to access CELLxGENE-hosted embeddings is by calling the `get_anndata` function with the `obs_embeddings` or `var_embeddings` parameter.

**Example**: Explore available embeddings and fetch scGPT embeddings for cells from the tongue tissue.

```python
from cellxgene_census.experimental import get_all_available_embeddings

CENSUS_VERSION = "2023-12-15"

# List available embeddings
for e in get_all_available_embeddings(CENSUS_VERSION):
    print(f"{e['embedding_name']:15} {e['experiment_name']:15} {e['data_type']:15}")

# Fetch data with scGPT embeddings
import cellxgene_census

with cellxgene_census.open_soma(census_version=CENSUS_VERSION) as census:
    adata = cellxgene_census.get_anndata(
        census,
        organism="homo_sapiens",
        measurement_name="RNA",
        obs_value_filter="tissue == 'tongue'",
        obs_embeddings=["scgpt"],
    )

print(adata)
print(adata.obsm)
```

**Output**:

```
AnnData object with n_obs × n_vars = [number of cells] × [number of genes]
    obs: 'assay', 'cell_type', 'tissue', 'disease', ...
    var: 'soma_joinid', 'feature_id', 'feature_name', 'feature_length', ...
    obsm: 'scgpt'
```

**Note**: Missing cells in the embedding matrix are represented with rows where all values are `NaN`.

### Storage Format

Each embedding is encoded as a `SOMASparseNDArray`, where:

- **Dimension 0 (`soma_dim_0`)**: Encodes the cell (`obs`) `soma_joinid` value.
- **Dimension 1 (`soma_dim_1`)**: Encodes the embedding features, ranging from 0 to N-1, where N is the number of features in the embedding.
- **Data (`soma_data`)**: Stored as `float32` with reduced precision (equivalent to `bfloat16`).

**Important**: CELLxGENE-hosted embeddings may embed a subset of the cells in any given Census version. Missing array values imply that the cell was not embedded.

### Querying Cells and Loading Associated Embeddings

#### Loading Embeddings into an AnnData `obsm` Slot Using `get_anndata()`

Fetch cells from the central nervous system and retrieve the scGPT embeddings.

```python
import cellxgene_census
import scanpy as sc

with cellxgene_census.open_soma(census_version=CENSUS_VERSION) as census:
    adata = cellxgene_census.get_anndata(
        census,
        organism="homo_sapiens",
        measurement_name="RNA",
        obs_value_filter="tissue_general == 'central nervous system'",
        obs_column_names=["cell_type", "soma_joinid"],
        obs_embeddings=["scgpt"],
    )

# Visualize the embeddings
sc.pp.neighbors(adata, use_rep="scgpt")
sc.tl.umap(adata)
sc.pl.umap(adata, color="cell_type", title="scGPT Embedding")
```

#### Loading Embeddings via `ExperimentAxisQuery`

Use an `ExperimentAxisQuery` for lazy evaluation before loading data.

```python
import cellxgene_census
import tiledbsoma as soma
from cellxgene_census.experimental import get_embedding

census = cellxgene_census.open_soma(census_version=CENSUS_VERSION)

experiment = census["census_data"]["homo_sapiens"]
query = experiment.axis_query(
    measurement_name="RNA",
    obs_query=soma.AxisQuery(value_filter="tissue_general == 'central nervous system'"),
)

# Get soma_joinid values
soma_joinids = query.obs_joinids().to_numpy()

# Create AnnData object
adata = query.to_anndata(X_name="raw", column_names={"obs": ["cell_type"]})

# Add embeddings
EMBEDDING_URI = "s3://cellxgene-contrib-public/contrib/cell-census/soma/2023-12-15/CxG-contrib-1/"
adata.obsm["scgpt"] = get_embedding(
    census_version=CENSUS_VERSION,
    embedding_uri=EMBEDDING_URI,
    obs_soma_joinids=soma_joinids,
)

# Close resources
query.close()
census.close()
```

#### Loading an Embedding into a Dense NumPy Array

Select cells based on metadata and load the corresponding embeddings.

```python
import cellxgene_census
from cellxgene_census.experimental import get_embedding

census = cellxgene_census.open_soma(census_version=CENSUS_VERSION)

obs_df = cellxgene_census.get_obs(
    census,
    "homo_sapiens",
    value_filter="tissue_general == 'exocrine gland'",
    column_names=["soma_joinid", "cell_type"],
)

# Get embeddings
EMBEDDING_URI = "s3://cellxgene-contrib-public/contrib/cell-census/soma/2023-12-15/CxG-contrib-1/"
embeddings = get_embedding(
    CENSUS_VERSION, EMBEDDING_URI, obs_df["soma_joinid"].to_numpy()
)

print("Embedding shape:", embeddings.shape)
```

### Embedding Metadata

Each embedding contains descriptive information stored in the `metadata` slot, encoded as a JSON string.

**Accessing Embedding Metadata**:

```python
from cellxgene_census.experimental import get_embedding_metadata

embedding_metadata = get_embedding_metadata(EMBEDDING_URI)
print(embedding_metadata)
```

**Validating Embedding Metadata**:

```python
assert embedding_metadata["census_version"] == CENSUS_VERSION
assert embedding_metadata["experiment_name"] == "homo_sapiens"
assert embedding_metadata["measurement_name"] == "RNA"

print("Embedding metadata is valid.")
```

---

## Finding the Most Similar Census Cells Using Embeddings Vector Search

### Overview

This section demonstrates how to find the most similar Census cells to any other data using embeddings vector search. We will generate scVI embeddings for some test cells, search the Census scVI embeddings for nearest neighbors, and use them to predict cell types and tissues of the test cells.

**Note**: This tutorial requires the experimental features of `cellxgene-census`. Ensure you have installed the package with `[experimental]` extras.

### Contents

1. [Downloading Data and Census scVI Model](#downloading-data-and-census-scvi-model)
2. [Loading and Embedding Test Cells](#loading-and-embedding-test-cells)
3. [Searching for Similar Census Cells](#searching-for-similar-census-cells)
4. [Predicting Cell Metadata](#predicting-cell-metadata)

### Downloading Data and Census scVI Model

First, install the required packages:

```shell
pip install 'cellxgene-census[experimental]' scvi-tools
```

#### Download Test Data

We will use the 10X PBMC 3K dataset as our test data.

```python
import os

os.makedirs('data', exist_ok=True)
!wget -q -O data/pbmc3k_filtered_gene_bc_matrices.tar.gz http://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz
!tar -xzf data/pbmc3k_filtered_gene_bc_matrices.tar.gz -C data/
```

#### Download the Census scVI Model

Determine the S3 location of the scVI model corresponding to the Census version.

```python
import cellxgene_census

CENSUS_VERSION = "2024-07-01"

with cellxgene_census.open_soma(census_version=CENSUS_VERSION) as census:
    scvi_info = cellxgene_census.experimental.get_embedding_metadata_by_name(
        embedding_name="scvi",
        organism="homo_sapiens",
        census_version=CENSUS_VERSION,
    )

model_link = scvi_info["model_link"]
print(model_link)
```

Download the scVI model:

```python
os.makedirs('scvi-human-2024-07-01', exist_ok=True)
!wget -q -O scvi-human-2024-07-01/model.pt {model_link}
```

### Loading and Embedding Test Cells

#### Load the Test Data into an AnnData Object

```python
import scanpy as sc

adata = sc.read_10x_mtx("data/filtered_gene_bc_matrices/hg19/", var_names="gene_ids")
adata.var["ensembl_id"] = adata.var.index
adata.obs["n_counts"] = adata.X.sum(axis=1)
adata.obs["joinid"] = list(range(adata.n_obs))
adata.obs["batch"] = "unassigned"
```

#### Embed the Test Cells Using scVI

```python
import scvi

scvi.model.SCVI.prepare_query_anndata(adata, "scvi-human-2024-07-01")
vae_q = scvi.model.SCVI.load_query_data(adata, "scvi-human-2024-07-01")

# Forward pass to get latent representation
vae_q.is_trained = True
latent = vae_q.get_latent_representation()
adata.obsm["scvi"] = latent
```

#### Clean Up the AnnData Object

```python
# Filter out missing features
adata = adata[:, adata.var["gene_symbols"].notnull().values].copy()
adata.var.set_index("gene_symbols", inplace=True)
adata.var_names = adata.var["ensembl_id"]
adata.obs["cell_type"] = "Query - PBMC 10X"
adata.obs["tissue_general"] = "Query - PBMC 10X"
```

#### Visualize the Test Cells

```python
import scanpy as sc

sc.pp.neighbors(adata, n_neighbors=15, use_rep="scvi")
sc.tl.umap(adata)
sc.tl.leiden(adata)
sc.pl.umap(adata, color="leiden")
```

### Searching for Similar Census Cells

Use the experimental API to search the vector index of scVI embeddings.

```python
from cellxgene_census.experimental import find_nearest_obs

neighbors = find_nearest_obs(
    "scvi", "homo_sapiens", CENSUS_VERSION, query=adata, k=30, memory_GiB=8, nprobe=20
)
```

**Note**: This function accesses the cell embeddings in the `scvi` obsm layer and searches for the *k* nearest neighbors among the Census cell embeddings.

Inspect the neighbors:

```python
print(neighbors)
```

#### Fetch Nearest Neighbor Data

```python
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
```

#### Visualize Nearest Neighbors

```python
sc.pp.neighbors(neighbors_adata, n_neighbors=15, use_rep="scvi")
sc.tl.umap(neighbors_adata)
sc.pl.umap(neighbors_adata, color="tissue_general")
```

#### Combine Query and Neighbors for Visualization

```python
import anndata

adata_concat = anndata.concat([adata, neighbors_adata])
sc.pp.neighbors(adata_concat, n_neighbors=15, use_rep="scvi")
sc.tl.umap(adata_concat)
sc.pl.umap(adata_concat, color=["tissue_general"])
```

### Predicting Cell Metadata

Use the nearest neighbors to predict metadata attributes like `tissue_general` and `cell_type`.

```python
from cellxgene_census.experimental import predict_obs_metadata

predictions = predict_obs_metadata(
    "homo_sapiens", CENSUS_VERSION, neighbors, ["tissue_general", "cell_type"]
)
print(predictions)
```

#### Add Predictions to AnnData Object

```python
import pandas as pd

predictions.index = adata.obs.index
predictions = predictions.rename(columns={"cell_type": "predicted_cell_type"})
adata.obs = pd.concat([adata.obs, predictions], axis=1)
```

#### Visualize Predicted Cell Types

```python
sc.pl.umap(adata, color="predicted_cell_type")
```

#### Annotate Clusters Based on Predictions

```python
# Annotate clusters
adata.obs["predicted_consolidated_cell_type"] = ""
for leiden_cluster in adata.obs["leiden"].unique():
    most_common_type = (
        adata.obs.loc[adata.obs["leiden"] == leiden_cluster, "predicted_cell_type"].mode()[0]
    )
    adata.obs.loc[adata.obs["leiden"] == leiden_cluster, "predicted_consolidated_cell_type"] = most_common_type

# Visualize consolidated predictions
sc.pl.umap(adata, color="predicted_consolidated_cell_type")
```

---

## Using the scVI Pretrained Model for Data Projection and Cell Type Prediction

### Overview

This section provides examples of how to utilize the pretrained **scVI model** with your own data for **data projection** and **cell type prediction**. For more information on the model, please refer to the [Census Models page](https://cellxgene.cziscience.com/census-models).

**Note**: This tutorial requires `cellxgene-census` package version **1.9.1** or later.

### Contents

1. [Requirements](#requirements)
2. [Preparing Data and Model](#preparing-data-and-model)
3. [Using the scVI Pretrained Model for Data Projection](#using-the-scvi-pretrained-model-for-data-projection)
4. [Using the scVI Pretrained Model for Cell Type Prediction](#using-the-scvi-pretrained-model-for-cell-type-prediction)

### Requirements

#### System Requirements

- **Operating System**: Unix-based system
- **Hardware**: A system with one or more GPUs is highly recommended
- **Software**:
  - [AWS Command-Line Interface](https://docs.aws.amazon.com/cli/latest/userguide/getting-started-install.html)
  - [scvi-tools](https://github.com/scverse/scvi-tools) and its dependencies
  - `cellxgene-census` package with experimental features

#### Downloading Example Data

We will use the 10X PBMC 3K dataset in this tutorial.

```shell
mkdir -p data
wget -nv -O data/pbmc3k_filtered_gene_bc_matrices.tar.gz http://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz
tar -xzf data/pbmc3k_filtered_gene_bc_matrices.tar.gz -C data/
```

#### Downloading the Trained scVI Model

The model is hosted on S3. You can find more details on the [Census Models page](https://cellxgene.cziscience.com/census-models).

**Retrieve Model Information**:

```python
import cellxgene_census

census_version = "2023-12-15"
organism = "homo_sapiens"

with cellxgene_census.open_soma(census_version=census_version) as census:
    scvi_info = cellxgene_census.experimental.get_embedding_metadata_by_name(
        embedding_name="scvi",
        organism=organism,
        census_version=census_version,
    )

model_link = scvi_info["model_link"]
print(model_link)
```

**Download the Model**:

```shell
mkdir -p scvi-homo-sapiens-2023-12-15
wget -q -O scvi-homo-sapiens-2023-12-15/model.pt {model_link}
```

### Preparing Data and Model

#### Import Required Packages

```python
import warnings

warnings.filterwarnings("ignore")

import anndata
import cellxgene_census
import numpy as np
import scanpy as sc
import scvi
from sklearn.ensemble import RandomForestClassifier
```

#### Load the Example Query Dataset

```python
adata = sc.read_10x_mtx("data/filtered_gene_bc_matrices/hg19/", var_names="gene_ids")
adata.var["ensembl_id"] = adata.var.index
adata.obs["n_counts"] = adata.X.sum(axis=1)
adata.obs["joinid"] = list(range(adata.n_obs))
adata.obs["batch"] = "unassigned"  # Placeholder batch label
```

### Using the scVI Pretrained Model for Data Projection

#### Prepare the Query Data

```python
model_filename = "scvi-homo-sapiens-2023-12-15"

scvi.model.SCVI.prepare_query_anndata(adata, model_filename)
```

#### Load the Model and Project the Data

```python
vae_q = scvi.model.SCVI.load_query_data(adata, model_filename)
vae_q.is_trained = True  # Trick the model into thinking it's trained
latent = vae_q.get_latent_representation()
adata.obsm["scvi"] = latent
```

#### Clean Up the AnnData Object

```python
# Filter out missing features
adata = adata[:, adata.var["gene_symbols"].notnull().values].copy()
adata.var.set_index("gene_symbols", inplace=True)
```

#### Perform UMAP Visualization

```python
sc.pp.neighbors(adata, n_neighbors=15, use_rep="scvi")
sc.tl.umap(adata)
```

#### Perform Leiden Clustering

```python
sc.tl.leiden(adata)
```

#### Normalize and Log-Transform the Data

```python
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
```

#### Map Clusters to Cell Type Labels

Based on marker genes from the [Scanpy pbmc3k tutorial](https://scanpy-tutorials.readthedocs.io/en/latest/pbmc3k.html), we can map Leiden clusters to cell types.

```python
markers_row1 = ["IL7R", "CD14", "LYZ", "MS4A1", "CD8A", "GNLY"]
markers_row2 = ["NKG7", "FCGR3A", "MS4A7", "FCER1A", "CST3", "PPBP"]

sc.pl.violin(adata, markers_row1, groupby="leiden")
sc.pl.violin(adata, markers_row2, groupby="leiden")
```

Assign cell type labels:

```python
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
```

#### Visualize Original Cell Types

```python
sc.pl.umap(adata, color=["original_cell_type"])
```

### Using the scVI Pretrained Model for Cell Type Prediction

#### Fetch Reference scVI Embeddings from Census

```python
with cellxgene_census.open_soma(census_version="2023-12-15") as census:
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
```

#### Combine Query and Reference Data for Visualization

```python
adata.obs["dataset_id"] = "QUERY"
index_subset = np.random.choice(adata_census.n_obs, size=3000, replace=False)
adata_census_subset = adata_census[index_subset, :]

adata_combined = anndata.concat([adata_census_subset, adata])
sc.pp.neighbors(adata_combined, n_neighbors=15, use_rep="scvi", metric="correlation")
sc.tl.umap(adata_combined)
sc.pl.umap(adata_combined, color=["dataset_id"])
```

#### Train a Classifier on Reference Data

```python
rfc = RandomForestClassifier()
rfc.fit(adata_census.obsm["scvi"], adata_census.obs["cell_type"].values)
```

#### Predict Cell Types on Query Data

```python
adata.obs["predicted_cell_type"] = rfc.predict(adata.obsm["scvi"])

# Calculate confidence scores
probabilities = rfc.predict_proba(adata.obsm["scvi"])
confidence = np.array([
    probabilities[i][rfc.classes_ == adata.obs["predicted_cell_type"][i]][0]
    for i in range(adata.n_obs)
])
adata.obs["predicted_cell_type_probability"] = confidence
```

#### Visualize Predicted Cell Types

```python
sc.pl.umap(adata, color="original_cell_type")
sc.pl.umap(adata, color=["predicted_cell_type_probability", "predicted_cell_type"])
```

#### Visualize Combined Data with Predicted Cell Types

```python
adata_combined.obs["cell_type"] = (
    adata_census_subset.obs["cell_type"].tolist() + adata.obs["predicted_cell_type"].tolist()
)
sc.pl.umap(adata_combined, color=["dataset_id", "cell_type"])
```

---

## Using the Geneformer Fine-Tuned Model for Cell Class Prediction and Data Projection

### Overview

This section provides examples of how to utilize the **CELLxGENE collaboration fine-tuned Geneformer model** with your own data for **cell subclass inference** and **data projection**. For more information on the model, please refer to the [Census Models page](https://cellxgene.cziscience.com/census-models).

**Note**:

- This tutorial requires `cellxgene-census` package version **1.9.1** or later.
- "Cell subclass" refers to high-level groupings of cell types as annotated in CELLxGENE Discover via the CL ontology. See [CZ CELLxGENE Discover Collections](https://cellxgene.cziscience.com/collections).
- The Census RNA data includes duplicate cells present across multiple datasets. Duplicate cells can be filtered using the `is_primary_data` column in cell metadata as described in the [Census schema](https://github.com/chanzuckerberg/cellxgene-census/blob/main/docs/cellxgene_census_schema.md#repeated-data).

### Contents

1. [Requirements](#requirements)
2. [Preparing Data and Model](#preparing-data-and-model)
3. [Using the Geneformer Fine-Tuned Model for Cell Subclass Inference](#using-the-geneformer-fine-tuned-model-for-cell-subclass-inference)
4. [Using the Geneformer Fine-Tuned Model for Data Projection](#using-the-geneformer-fine-tuned-model-for-data-projection)

### Requirements

#### System Requirements

- **Operating System**: Unix-based system
- **Hardware**: A system with one or more GPUs is highly recommended
- **Software**:
  - [AWS Command-Line Interface](https://docs.aws.amazon.com/cli/latest/userguide/getting-started-install.html)
  - [Geneformer Python Package](https://huggingface.co/ctheodoris/Geneformer) and its dependencies
  - `cellxgene-census` package with experimental features

#### Downloading Example Data

We will use the 10X PBMC 3K dataset in this tutorial.

```shell
mkdir -p data
wget -nv -O data/pbmc3k_filtered_gene_bc_matrices.tar.gz http://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz
tar -xzf data/pbmc3k_filtered_gene_bc_matrices.tar.gz -C data/
```

#### Downloading the Fine-Tuned Geneformer Model

The model is hosted on S3. You can find more details on the [Census Models page](https://cellxgene.cziscience.com/census-models).

**Retrieve Model Information**:

```python
import cellxgene_census

census_version = "2023-12-15"
organism = "homo_sapiens"

with cellxgene_census.open_soma(census_version=census_version) as census:
    geneformer_info = cellxgene_census.experimental.get_embedding_metadata_by_name(
        embedding_name="geneformer",
        organism=organism,
        census_version=census_version,
    )

model_link = geneformer_info["model_link"]
print(model_link)
```

**Download the Model**:

```shell
mkdir -p fine_tuned_geneformer
aws s3 sync --no-sign-request --no-progress --only-show-errors {model_link} ./fine_tuned_geneformer
```

### Preparing Data and Model

#### Import Required Packages

```python
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
```

#### Load and Prepare the Test Data

```python
adata = sc.read_10x_mtx("data/filtered_gene_bc_matrices/hg19/", var_names="gene_ids")
adata.var["ensembl_id"] = adata.var.index
adata.obs["n_counts"] = adata.X.sum(axis=1)
adata.obs["joinid"] = list(range(adata.n_obs))

h5ad_dir = "./data/h5ad/"
os.makedirs(h5ad_dir, exist_ok=True)
adata.write(h5ad_dir + "pbmcs.h5ad")
```

#### Tokenize the Data Using Geneformer Tokenizer

```python
token_dir = "data/tokenized_data/"
os.makedirs(token_dir, exist_ok=True)

tokenizer = TranscriptomeTokenizer(custom_attr_name_dict={"joinid": "joinid"})
tokenizer.tokenize_data(
    data_directory=h5ad_dir,
    output_directory=token_dir,
    output_prefix="pbmc",
    file_format="h5ad",
)
```

#### Load the Label Mapping Dictionary

```python
model_dir = "./fine_tuned_geneformer/"
label_mapping_dict_file = os.path.join(model_dir, "label_to_cell_subclass.json")

with open(label_mapping_dict_file) as fp:
    label_mapping_dict = json.load(fp)

print(label_mapping_dict)
```

### Using the Geneformer Fine-Tuned Model for Cell Subclass Inference

#### Load Tokenized Data

```python
dataset = datasets.load_from_disk(token_dir + "pbmc.dataset")
# Add a dummy 'label' column required for prediction
dataset = dataset.add_column("label", [0] * len(dataset))
```

#### Perform Inference

**Note**: This step can be slow on CPUs; a machine with a GPU is recommended.

```python
# Load the fine-tuned model
model = BertForSequenceClassification.from_pretrained(model_dir)
# Create the trainer
trainer = Trainer(model=model, data_collator=DataCollatorForCellClassification())
# Make predictions
predictions = trainer.predict(dataset)
```

#### Process Predictions

```python
predicted_label_ids = np.argmax(predictions.predictions, axis=1)
predicted_logits = [predictions.predictions[i][predicted_label_ids[i]] for i in range(len(predicted_label_ids))]
predicted_labels = [label_mapping_dict[str(i)] for i in predicted_label_ids]
```

#### Add Predictions to AnnData

```python
adata.obs["predicted_cell_subclass"] = predicted_labels
adata.obs["predicted_cell_subclass_probability"] = np.exp(predicted_logits) / (1 + np.exp(predicted_logits))
```

#### Visualize Predictions

```python
# Basic processing for UMAP visualization
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata = adata[:, adata.var.highly_variable]
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver="arpack")
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)
```

Add original cell type annotations from the [Scanpy pbmc3k tutorial](https://scanpy-tutorials.readthedocs.io/en/latest/pbmc3k.html):

```python
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
```

Visualize original annotations:

```python
sc.pl.umap(adata, color="leiden", title="Original Annotations")
```

Visualize predicted annotations:

```python
sc.pl.umap(
    adata,
    color=["predicted_cell_subclass_probability", "predicted_cell_subclass"],
    title="Predicted Geneformer Annotations",
)
```

### Using the Geneformer Fine-Tuned Model for Data Projection

#### Generate Geneformer Embeddings for Test Data

```python
n_classes = len(label_mapping_dict)

output_dir = "data/geneformer_embeddings"
os.makedirs(output_dir, exist_ok=True)

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
```

#### Merge Embeddings with AnnData

```python
embs = embs.sort_values("joinid")
adata.obsm["geneformer"] = embs.drop(columns="joinid").to_numpy()
```

#### Visualize Geneformer Embeddings

```python
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40, use_rep="geneformer")
sc.tl.umap(adata)
sc.pl.umap(adata, color="predicted_cell_subclass", title="10X PBMC 3K in Geneformer")
```

#### Join Geneformer Embeddings with Census Data

Fetch PBMC datasets from Census with Geneformer embeddings:

```python
with cellxgene_census.open_soma(census_version="2023-12-15") as census:
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
```

Select shared genes:

```python
adata_census.var_names = adata_census.var["feature_id"]
shared_genes = list(set(adata.var_names) & set(adata_census.var_names))
adata_census = adata_census[:, shared_genes]
```

Subset the Census data:

```python
index_subset = np.random.choice(adata_census.n_obs, size=3000, replace=False)
adata_census = adata_census[index_subset, :]
```

Join the datasets:

```python
adata_census.obs["dataset"] = "Census - " + adata_census.obs["dataset_id"].astype(str)
adata.obs["dataset"] = "10X PBMC 3K"
adata.obs["cell_type"] = "Predicted - " + adata.obs["predicted_cell_subclass"].astype(str)

adata_joined = sc.concat([adata, adata_census], join="outer", label="batch")
```

Visualize the combined data:

```python
sc.pp.neighbors(adata_joined, n_neighbors=10, n_pcs=40, use_rep="geneformer")
sc.tl.umap(adata_joined)
sc.pl.umap(adata_joined, color="dataset")
sc.pl.umap(adata_joined, color="cell_type")
```

---

## Using Census with PyTorch for Machine Learning

### Training PyTorch Models Directly with Census Data

Leverage the specialized PyTorch DataPipe `ExperimentDataPipe` for efficient data loading.

**Example**:

```python
import cellxgene_census
import cellxgene_census.experimental.ml as census_ml
import tiledbsoma as soma
import torch

with cellxgene_census.open_soma() as census:
    experiment = census["census_data"]["homo_sapiens"]

    experiment_datapipe = census_ml.ExperimentDataPipe(
        experiment,
        measurement_name="RNA",
        X_name="raw",
        obs_query=soma.AxisQuery(
            value_filter="tissue_general == 'tongue' and is_primary_data == True"
        ),
        obs_column_names=["cell_type"],
        batch_size=128,
        shuffle=True,
    )

    # Split into training and test sets
    train_datapipe, test_datapipe = experiment_datapipe.random_split(
        weights={"train": 0.8, "test": 0.2}, seed=1
    )

    # Create DataLoader
    train_loader = census_ml.experiment_dataloader(train_datapipe)

    # Define your model
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = YourModel().to(device)

    # Training loop
    for data in train_loader:
        inputs, labels = data
        # Training steps...
```

### Custom Data Encoders

Use custom encoders for cell metadata to transform or interact with different metadata variables.

**Example**:

```python
from cellxgene_census.experimental.ml.encoders import LabelEncoder

label_encoder = LabelEncoder("cell_type")

experiment_datapipe = census_ml.ExperimentDataPipe(
    experiment,
    measurement_name="RNA",
    obs_query=soma.AxisQuery(
        value_filter="tissue_general == 'tongue' and is_primary_data == True"
    ),
    obs_column_names=["cell_type"],
    batch_size=128,
    shuffle=True,
    encoders={"cell_type": label_encoder},
)
```

---

## Benchmarks of Census Models

### Evaluating Embeddings

Benchmarks provide insights into the biological signal and batch correction levels of embeddings.

**Bio-conservation Metrics**:

- **Leiden NMI and ARI**: Measures the agreement between clustering and biological labels.
- **Silhouette Score**: Evaluates the cohesion within clusters.

**Batch-Correction Metrics**:

- **Silhouette Batch Score**: Assesses batch mixing.
- **Entropy of Batch Mixing**: Evaluates diversity of batches within neighborhoods.
- **Classifier Accuracy**: Measures the ability to predict batch labels (lower is better for batch correction).

**Accessing Embeddings**:

```python
import cellxgene_census.experimental

embeddings = cellxgene_census.experimental.get_all_available_embeddings(census_version="2023-12-15")
```

---

## Accessing Census Data in AWS

### AWS CLI for Programmatic Downloads

Download Census data directly using AWS CLI.

**Example**:

```shell
aws s3 sync --no-sign-request s3://cellxgene-census-public-us-west-2/cell-census/2023-07-25/h5ads/ ./h5ads/
```

### Accessing Census Data via the API

Specify the Census version when opening the Census.

```python
import cellxgene_census

with cellxgene_census.open_soma(census_version="2023-07-25") as census:
    # Your code here
```

### Accessing a Local Copy of the Census

If you have a local copy of the Census data:

```python
import cellxgene_census

with cellxgene_census.open_soma(uri="local/path/to/soma/") as census:
    # Your code here
```

---

## Querying Data Using the gget cellxgene Module

### Overview

The [`gget`](https://github.com/pachterlab/gget) package is a free, open-source command-line tool and Python package that enables efficient querying of genomic databases. The `gget cellxgene` module builds on the CZ CELLxGENE Discover Census to query data from CZ CELLxGENE Discover.

If you use `gget cellxgene` in a publication, please [cite gget](https://pachterlab.github.io/gget/en/cite.html) in addition to [citing CZ CELLxGENE](https://cellxgene.cziscience.com/docs/08__Cite%20cellxgene%20in%20your%20publications).

### Installation

Install `gget` and set up the `cellxgene` module:

```python
# Install gget (version 0.25.7 or later)
!pip install -q 'gget>=0.25.7'

import gget

# Set up the cellxgene module
gget.setup("cellxgene")
```

### Fetching an AnnData Object by Selecting Genes, Tissues, and Cell Types

You can fetch an `AnnData` object by specifying genes, tissues, and cell types of interest.

**Example**:

```python
# Fetch AnnData object based on specified genes, tissue, and cell types
adata = gget.cellxgene(
    gene=["ACE2", "ABCA1", "SLC5A1"],
    tissue="lung",
    cell_type=["mucus secreting cell", "neuroendocrine cell"]
)

# Inspect the AnnData object
print(adata)
```

**Output**:

```
AnnData object with n_obs × n_vars = [number of cells] × 3
    obs: 'assay', 'cell_type', 'tissue', 'suspension_type', 'disease', 'sex', ...
    var: 'feature_name', 'feature_length', 'feature_id', ...
```

### Plotting a Dot Plot Similar to CZ CELLxGENE Discover Gene Expression

Using the fetched data, you can create a dot plot similar to those on the CZ CELLxGENE Discover Gene Expression page.

**Example**:

```python
import scanpy as sc

# Increase the resolution of plots displayed in notebooks
%config InlineBackend.figure_format="retina"

# Plot the dot plot
sc.pl.dotplot(
    adata,
    var_names=adata.var["feature_name"].values,
    groupby="cell_type",
    gene_symbols="feature_name"
)
```

### Fetching Only Cell Metadata

By setting `meta_only=True`, you can fetch only the cell metadata (`AnnData.obs`).

**Example**:

```python
df = gget.cellxgene(
    meta_only=True,
    census_version="2023-05-15",  # Specify Census version for reproducibility
    gene="ENSMUSG00000015405",
    ensembl=True,  # Indicate that the gene is provided as an Ensembl ID
    tissue="lung",
    species="mus_musculus",
)

# Display the metadata DataFrame
print(df)
```

### Using gget cellxgene from the Command Line

All `gget` modules support command-line usage. Note that the command-line interface requires the `-o/--out` argument to specify a path to save the fetched data.

**Example**:

```shell
# Fetch AnnData object based on specified genes, tissue, and cell types
gget cellxgene --gene ACE2 ABCA1 SLC5A1 --tissue lung --cell_type 'mucus secreting cell' 'neuroendocrine cell' -o example_adata.h5ad

# Fetch only metadata
gget cellxgene --meta_only --gene ENSMUSG00000015405 --ensembl --tissue lung --species mus_musculus -o example_meta.csv
```

---

## Additional Resources

- **Census Data and Schema**: Detailed information about the Census data structure and schema.
- **Python Tutorials**: Explore more examples and tutorials in the [CELLxGENE Census GitHub repository](https://github.com/chanzuckerberg/cellxgene-census/tree/main/notebooks).
- **gget Documentation**:
  - [gget cellxgene Module](https://pachterlab.github.io/gget/en/cellxgene.html)
  - [gget GitHub Repository](https://github.com/pachterlab/gget)
- **Documentation**:
  - [TileDB-SOMA](https://tiledbsoma.readthedocs.io/en/latest/)
  - [Anndata](https://anndata.readthedocs.io/en/latest/)
  - [Scanpy](https://scanpy.readthedocs.io/en/stable/)
  - [Geneformer](https://huggingface.co/ctheodoris/Geneformer)
  - [scvi-tools](https://scvi-tools.org/)
  - [PyTorch](https://pytorch.org/docs/stable/)
  - [TorchData](https://pytorch.org/data/beta/)
  - [NumPy](https://numpy.org/doc/)
  - [pandas](https://pandas.pydata.org/docs/)
  - [SciPy](https://docs.scipy.org/doc/scipy/reference/)

---

## Citing Census

To cite the Census project, please follow the [citation guidelines](https://cellxgene.cziscience.com/docs/08__Cite%20cellxgene%20in%20your%20publications) provided by CZ CELLxGENE Discover.

To cite individual studies, refer to the tutorial on [Generating citations for Census slices](https://github.com/chanzuckerberg/cellxgene-census/blob/main/notebooks/api_demo/census_citation_generation.ipynb).

If you use the `gget cellxgene` module, please also [cite gget](https://pachterlab.github.io/gget/en/cite.html).

---

## Feedback and Support

- **Questions and Feature Requests**: Submit via [GitHub Issues](https://github.com/chanzuckerberg/cellxgene-census/issues).
- **Slack Community**: Join the CZI Science Community on Slack ([czi.co/science-slack](https://czi.co/science-slack)) and ask questions in the `#cellxgene-census-users` channel.
- **Email**: Send feedback to [soma@chanzuckerberg.com](mailto:soma@chanzuckerberg.com).
- **Security Issues**: Report to [security@chanzuckerberg.com](mailto:security@chanzuckerberg.com).
- **FAQs**: Visit the [Census FAQ](https://cellxgene.cziscience.com/docs/cellxgene_census_docsite_FAQ.md) page.
