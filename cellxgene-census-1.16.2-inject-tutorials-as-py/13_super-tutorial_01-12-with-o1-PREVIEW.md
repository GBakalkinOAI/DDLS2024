# CZ CELLxGENE Discover Census Super-Tutorial

Welcome to the comprehensive guide for using the **latest CELLxGENE Census Python API (v1.16.2)**. This tutorial will help you access, query, and analyze single-cell RNA sequencing data from the CZ CELLxGENE Discover Census efficiently. It is designed to provide you with clear, coherent instructions and examples to accelerate your research.

---

## Contents

1. [Introduction](#introduction)
2. [Installation](#installation)
3. [Quick Start](#quick-start)
4. [Understanding the Census Data and Schema](#understanding-the-census-data-and-schema)
5. [Accessing and Querying Census Data](#accessing-and-querying-census-data)
6. [Memory-Efficient Computations](#memory-efficient-computations)
7. [Working with Pre-calculated Embeddings and Models](#working-with-pre-calculated-embeddings-and-models)
8. [Using Census with PyTorch for Machine Learning](#using-census-with-pytorch-for-machine-learning)
9. [Benchmarks of Census Models](#benchmarks-of-census-models)
10. [Accessing Census Data in AWS](#accessing-census-data-in-aws)
11. [Additional Resources](#additional-resources)
12. [Citing Census](#citing-census)
13. [Feedback and Support](#feedback-and-support)

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

   To include experimental features like PyTorch loaders, install with:

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
        column_names={"obs": ["assay", "cell_type", "tissue", "tissue_general", "suspension_type", "disease"]},
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
- **Geneformer**: Fine-tuned embeddings.
- **scGPT**: Generative pre-trained transformer embeddings.
- **Universal Cell Embeddings (UCE)**: Pre-computed embeddings available in the Census.

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

## Additional Resources

- **Census Data and Schema**: Detailed information about the Census data structure and schema.
- **Python Tutorials**: Explore more examples and tutorials in the [CELLxGENE Census GitHub repository](https://github.com/chanzuckerberg/cellxgene-census/tree/main/notebooks).
- **Documentation**:
  - [TileDB-SOMA](https://tiledbsoma.readthedocs.io/en/latest/)
  - [Anndata](https://anndata.readthedocs.io/en/latest/)
  - [Scanpy](https://scanpy.readthedocs.io/en/stable/)
  - [PyTorch](https://pytorch.org/docs/stable/)
  - [TorchData](https://pytorch.org/data/beta/)
  - [NumPy](https://numpy.org/doc/)
  - [pandas](https://pandas.pydata.org/docs/)
  - [SciPy](https://docs.scipy.org/doc/scipy/reference/)

---

## Citing Census

To cite the Census project, please follow the [citation guidelines](https://cellxgene.cziscience.com/docs/08__Cite%20cellxgene%20in%20your%20publications) provided by CZ CELLxGENE Discover.

To cite individual studies, refer to the tutorial on [Generating citations for Census slices](https://github.com/chanzuckerberg/cellxgene-census/blob/main/notebooks/api_demo/census_citation_generation.ipynb).

---

## Feedback and Support

- **Questions and Feature Requests**: Submit via [GitHub Issues](https://github.com/chanzuckerberg/cellxgene-census/issues).
- **Slack Community**: Join the CZI Science Community on Slack ([czi.co/science-slack](https://czi.co/science-slack)) and ask questions in the `#cellxgene-census-users` channel.
- **Email**: Send feedback to [soma@chanzuckerberg.com](mailto:soma@chanzuckerberg.com).
- **Security Issues**: Report to [security@chanzuckerberg.com](mailto:security@chanzuckerberg.com).
- **FAQs**: Visit the [Census FAQ](https://cellxgene.cziscience.com/docs/cellxgene_census_docsite_FAQ.md) page.

---

**Happy Analyzing!**
