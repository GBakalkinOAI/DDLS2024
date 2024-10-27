
# Section: cellxgene_census_aws_open_data

# CZ CELLxGENE Discover Census in AWS

The single-cell data from [CZ CELLxGENE Discover Census](cellxgene_census_docsite_landing.md) are available for public access via Amazon Web Services (AWS).

This page describes what Census data are available in AWS and how to access them.

Contents

- [Census data available in AWS](#census-data-available-in-aws)
- [How to access AWS Census data](#how-to-access-aws-census-data)

## Census data available in AWS

The single-cell data from CZ CELLxGENE Discover included in Census (see [inclusion criteria](cellxgene_census_docsite_schema.md#data-included-in-the-census)) are available either as Census-wide TileDB files or individual H5AD files of the source datasets.

### Data specifications

<table class="custom-table">
<thead>
    <th>Data</th>
    <th>Format</th>
    <th>Access API</th>
    <th>Data schema</th>
    <th>Root S3 bucket</th>
    <th>Regions</th>
  </tr>
</thead>
<tbody>
  <tr>
    <td rowspan="2">Census-wide</td>
    <td rowspan="2"><a href="https://github.com/TileDB-Inc/TileDB">TileDB</a></td>
    <td><a href="https://github.com/chanzuckerberg/cellxgene-census/tree/main">CELLxGENE-Census</a></td>
    <td rowspan="2"><a href="https://github.com/chanzuckerberg/cellxgene-census/blob/main/docs/cellxgene_census_schema.md">CZ CELLxGENE Discover <b>Census</b> Schema</a></td>
    <td rowspan="2">s3://cellxgene-census-public-us-west-2/cell-census/<code>[tag]</code>/soma/</td>
    <td rowspan="3">us-west-2</td>
  </tr>
  <tr>
    <td><a href="https://github.com/single-cell-data/TileDB-SOMA">TileDB-SOMA</a></td>
  </tr>
  <tr>
    <td rowspan>Source datasets</td>
    <td rowspan"><a href="https://anndata.readthedocs.io/en/latest/fileformat-prose.html#elements">H5AD</a></td>
    <td><a href="https://anndata.readthedocs.io/en/latest/index.html">AnnData</a></td>
    <td rowspan><a href="https://github.com/chanzuckerberg/single-cell-curation/tree/main/schema">CZ CELLxGENE Discover <b>Dataset</b> Schema</a></td>
    <td rowspan>s3://cellxgene-census-public-us-west-2/cell-census/<code>[tag]</code>/h5ads/</td>
  </tr>
</tbody>
</table>

See the next section for a definition of `[tag]`.

### Data release versioning

A data release is a Census build that is publicly hosted in AWS. A Census build is a TileDB-SOMA collection and its corresponding source H5AD files with the Census data from CZ CELLxGENE Discover.

Any given Census build is named with a unique `[tag]`, normally the date of build, e.g. "2023-05-15".

The are two types of data releases:

- Long-Term Supported (LTS).
- Weekly.

For more information and for a list of all LTS Census data releases available please refer to [Census data releases](cellxgene_census_docsite_data_release_info.md).

## How to access AWS Census data

### AWS CLI for programatic downloads

Users can bulk-download Census data via the [AWS CLI](https://aws.amazon.com/cli/).

For example, to download the H5ADs files of the Census LTS release `2023-07-25`, users can execute the following from a shell session:

```bash
aws s3 sync --no-sign-request s3://cellxgene-census-public-us-west-2/cell-census/2023-07-25/h5ads/ ./h5ads/
```

And to download the TileDB files:

```bash
aws s3 sync --no-sign-request s3://cellxgene-census-public-us-west-2/cell-census/2023-07-25/soma/ ./soma/
```

### CELLxGENE Census API (Python and R)

This is the recommend method for accessing Census data. Please follow the [Census API quick start guide](cellxgene_census_docsite_quick_start.md) for a full guide.

For example, in Python users can create an iterator for the cell metadata Data Frame as follows:

``` python
import cellxgene_census

with cellxgene_census.open_soma() as census:
    cell_metadata = census["census_data"]["homo_sapiens"].obs.read(
        value_filter = "sex == 'female' and cell_type in ['microglial cell', 'neuron']",
        column_names = ["assay", "cell_type", "tissue", "tissue_general", "suspension_type", "disease"]
    )
```

If a **local copy** of the Census data exists, users can access it by providing the path to the `soma/` folder.

``` python
import cellxgene_census

with cellxgene_census.open_soma(uri="local/path/to/soma/") as census:
   ...
```

If a copy of the Census data exists in a **private S3 bucket**, users can access it by providing the URI `soma/`
folder in the S3 bucket. This will also require customizing TileDB configuration options to specify the
bucket's AWS region and that signed requests should be used for S3 API calls. This can be done as follows:

``` python
import cellxgene_census

uri = "s3://my-private-data-bucket/cell-census/2023-07-25/soma/"

tiledb_config={"vfs.s3.no_sign_request": "false",
               "vfs.s3.region": "us-east-1"}

with cellxgene_census.open_soma(uri=uri, tiledb_config=tiledb_config) as census:
   ...
```

### TileDB-SOMA API (Python and R)

The Census API provides convenience wrappers for TileDB-SOMA to access the Census Data hosted at AWS. Users can interact directly with the Census TileDB data directly via the TileDB-SOMA APIs. Please refer to the [TileDb-SOMA documentation](https://tiledbsoma.readthedocs.io/en/latest/) for full details on usage.

For example, in Python users can create an iterator for the cell metadata Data Frame as follows:

``` python
import cellxgene_census
import tiledbsoma

uri = "s3://cellxgene-census-public-us-west-2/cell-census/2023-07-25/soma/"
ctx = cellxgene_census.get_default_soma_context()

with tiledbsoma.open(uri, context=ctx) as census:
    cell_metadata = census["census_data"]["homo_sapiens"].obs.read(
        value_filter = "sex == 'female' and cell_type in ['microglial cell', 'neuron']",
        column_names = ["assay", "cell_type", "tissue", "tissue_general", "suspension_type", "disease"]
    )
```



# Section: cell_census_pypi_redirect-README

# cell-census is now cellxgene-census

This package has been renamed. Use `pip install cellxgene-census` instead.

New package: [https://pypi.org/project/cellxgene-census/](https://pypi.org/project/cellxgene-census/)



# Section: cellxgene_census_docsite_landing

<span style="background-color: #f3bfcb; color; font-size: 18px"> 🚀 New to the Census: Train PyTorch models directly with Census data with our efficient and easy-to-use PyTorch loaders. [Learn more](https://chanzuckerberg.github.io/cellxgene-census/articles/2024/20240709-pytorch.html)!

<span style="background-color: #f3bfcb; color; font-size: 18px"> 💻 Explore benchmarks of Census models and embeddings. [See the report](https://chanzuckerberg.github.io/cellxgene-census/articles/2024/20240710_embedding_metrics_dec_2023_lts.html)!

 </span>

# CZ CELLxGENE Discover Census

The Census provides efficient computational tooling to **access, query, and analyze all single-cell RNA data from CZ CELLxGENE Discover**. Using a new access paradigm of cell-based slicing and querying, you can interact with the data through TileDB-SOMA, or get slices in AnnData, Seurat, or SingleCellExperiment objects, thus accelerating your research by significantly minimizing data harmonization.

Get started:

- [Quick start (Python and R)](cellxgene_census_docsite_quick_start.md)
- [Census data and schema](cellxgene_census_docsite_schema.md)
- [Python tutorials](examples.rst)
- [R tutorials](https://chanzuckerberg.github.io/cellxgene-census/r/articles/)
- [Github repository](https://github.com/chanzuckerberg/cellxgene-census)

![image](cellxgene_census_docsite_workflow.svg)

## Citing Census

To cite the project please follow the [citation guidelines](https://cellxgene.cziscience.com/docs/08__Cite%20cellxgene%20in%20your%20publications) offered by CZ CELLxGENE Discover.

To cite individual studies please refer to the tutorial [Generating citations for Census slices](notebooks/api_demo/census_citation_generation.ipynb).

## Census Capabilities

The Census is a data object publicly hosted online and an API to open it. The object is built using the [SOMA](https://github.com/single-cell-data/SOMA) API specification and data model, and it is implemented via [TileDB-SOMA](https://github.com/single-cell-data/TileDB-SOMA). As such, the Census has all the data capabilities offered by TileDB-SOMA including:

**Data access at scale:**

- Cloud-based data access.
- Efficient access for larger-than-memory slices of data.
- Query and access data based on cell or gene metadata at low latency.

**Interoperability with existing single-cell toolkits:**

- Load and create [AnnData](https://anndata.readthedocs.io/en/latest/) objects.
- Load and create [Seurat](https://satijalab.org/seurat/) objects.
- Load and create [SingleCellExperiment](https://bioconductor.org/packages/release/bioc/html/SingleCellExperiment.html) objects.

**Interoperability with existing Python or R data structures:**

- From Python create [PyArrow](https://arrow.apache.org/docs/python/index.html) objects, SciPy sparse matrices, NumPy arrays, and pandas data frames.
- From R create [R Arrow](https://arrow.apache.org/docs/r/index.html) objects, sparse matrices (via the [Matrix](https://cran.r-project.org/package=Matrix) package), and standard data frames and (dense) matrices.

## Census Data and Schema

A description of the Census data and its schema is detailed [here](cellxgene_census_docsite_schema.md).

⚠️ Note that the data includes:

- **Full-gene sequencing read counts** (e.g. Smart-Seq2) and **molecule counts** (e.g. 10X).
- **Duplicate cells** present across multiple datasets, these can be filtered in or out using the cell metadata variable `is_primary_data`.

## Census Data Releases

The Census data release plans are detailed [here](cellxgene_census_docsite_data_release_info.md).

Starting May 15th, 2023, Census data releases with long-term support will be published every six months. These releases will be publicly accessible for at least five years. In addition, weekly releases may be published without any guarantee of permanence.

## Questions, Feedback and Issues

- Users are encouraged to submit questions and feature requests about the Census via [github issues](https://github.com/chanzuckerberg/cellxgene-census/issues).
- For quick support, you can join the CZI Science Community on Slack ([czi.co/science-slack](https://czi.co/science-slack)) and ask questions in the `#cellxgene-census-users` channel.
- Users are encouraged to share their feedback by emailing <soma@chanzuckerberg.com>.
- Bugs can be submitted via [github issues](https://github.com/chanzuckerberg/cellxgene-census/issues).
- If you believe you have found a security issue, please disclose it by contacting <security@chanzuckerberg.com>.
- Additional FAQs can be found [here](cellxgene_census_docsite_FAQ.md).

## Coming Soon!

- We are currently working on creating the tooling necessary to perform data modeling at scale with seamless integration of the Census and [PyTorch](https://pytorch.org/).
- To increase the usability of the Census for research, in 2023 and 2024 we are planning to explore the following areas:
  - Include organism-wide normalized layers.
  - Include organism-wide embeddings.
  - On-demand information-rich subsampling.

## Projects and Tools Using Census

If you are interested in listing a project here, please reach out to us at <soma@chanzuckerberg.com>



# Section: notebooks-README

# ReadMe

Demonstration notebooks for the CZ CELLxGENE Discover Census. There are two kinds of notebooks:

1. **API mechanics** — Located under `api_demo` these notebooks provide technical demonstrations of the Census API capabilities.
2. **Computational biology analysis** — Located under `analysis_demo` these notebooks provide an overview of the data in the Census, how to access it and how to use the it in an analytical framework.

## Dependencies

You must be on a Linux or MacOS system, with the following installed:

* Python 3.10 to 3.12
* Jupyter or some other means of running notebooks (e.g., vscode)

For now, it is recommended that you do all this on a host with sufficient memory,
and a high bandwidth connection to AWS S3 in the us-west-2 region, e.g., an m6i.8xlarge.
If you utilize AWS, Ubuntu 20 or 22 AMI are recommended (AWS AMI should work fine, but has
not been tested).

I also recommend you use a `d` instance type, and mount all of the NVME drives as swap,
as it will keep you from running out of RAM.

## Set up Python environment

1. (optional, but highly recommended) In your working directory, make and activate a virtual environment. For example:

    ```shell
      python -m venv ./venv
      source ./venv/bin/activate
    ```

2. Install the required dependencies:

    ```shell
      pip install -U -r cellxgene-census/api/python/notebooks/requirements.txt
    ```

## Verify your installation

Check that your installation works - this make take a few seconds, as it loads metadata from S3:

```shell
$ python -c 'import cellxgene_census; print(cellxgene_census.open_soma().soma_type)'
SOMACollection
```

## Run notebooks

Run notebooks, which you can find in the `cellxgene-census/api/python/notebooks` directory.

## For more help

If you have difficulties or questions, feel to reach out to us using Github issues, or any of the other means described in the [README](../../../README.md).



# Section: articles-2023-20231012-normalized_layer_precalc_stats

# Introducing a normalized layer and pre-calculated cell and gene statistics in Census

<!-- markdownlint-disable MD036 -->
*Published: October 12, 2023*

*By: [Maximilian Lombardo](mailto:mlombardo@chanzuckerberg.com) and [Pablo Garcia-Nieto](mailto:pgarcia-nieto@chanzuckerberg.com)*
<!-- markdownlint-enable MD036 -->

The Census team is happy to announce the introduction of two new data features, tailored to empower your single-cell research: a library-size, normalized expression layer and pre-calculated cell and gene statistics. This work is reflected in changes introduced in the [Census schema V1.1.0](https://github.com/chanzuckerberg/cellxgene-census/blob/main/docs/cellxgene_census_schema.md).

With these additions users can easily:

- Expand their Census query filters to select genes or cells based on the new metadata. For example, selecting only cells with N number of genes expressed.
- Export the expanded cell and gene metadata for downstream analysis.
- Export a normalized expression data matrix for downstream analysis.

These features are currently exclusive to the "latest" versions of the Census data release and they will be available in the next LTS data release. We invite your feedback as you explore these novel functionalities.

Keep on reading to find out more about these features!

## Description of new data added to Census

All of the following changes were introduced in the Census schema V1.1.0.

### Added a new library-size normalized layer

We have introduced a library-size normalized X layer for the RNA measurements of both the human and mouse experiments available as `X["normalized"]`.  The normalized layer is built by dividing each value in the raw count matrix by its corresponding row sum (i.e. size normalization).

To reduce data size and improve performance, normalized values are stored with a reduced floating point precision. In addition, to ensure that small count values do not round to zero, a small sigma has been added. You will see the effect of these artifacts in row (per-cell) values not summing to precisely 1.0.

### Enhanced gene metadata

The `ms["RNA"].var` DataFrame for both the human and mouse experiments has been enriched with two new metadata fields:

- `nnz` — the number of explicitly stored values, effectively the number of cells expressing this gene.
- `n_measured_obs` — the "measured" cells for this gene, effectively the number of cells for which this gene was measured in their respective dataset.

### Enhanced cell metadata

The `obs` DataFrame for both the human and mouse experiments is now augmented with the following new metadata, allowing users to forego common calculations used in early data pre-processing. For each cell:

- `raw_sum` — the sum of the raw counts, derived from `X["raw"]`.
- `nnz` — the number of explicitly stored values, effectively the number of genes expressed on this cell.
- `raw_mean_nnz` — the average counts from explicitly stored values.
- `raw_variance_nnz` — the variance of the counts from explicitly stored values.
- `n_measured_vars` — the "measured" genes, effectively the number of genes measured in the dataset from which the cell originated, thus all cells from the same dataset have the same value for this variable.

## How to use the new features

### Exporting the normalized data to existing single-cell toolkits

In Python, the normalized data can be exported into AnnData specifying the `X_name = "normalized"` argument of the `cellxgene.get_anndata()` method.

```python
import cellxgene_census

with cellxgene_census.open_soma(census_version = "latest") as census
    adata = cellxgene_census.get_anndata
        census = census,
        organism = "Homo sapiens",
        var_value_filter = "feature_id in ['ENSG00000161798', 'ENSG00000188229']",
        obs_value_filter = "cell_type == 'sympathetic neuron'",
        column_names = {"obs": ["tissue", "cell_type"]},
        X_name = "normalized" # Specify the normalized layer for this query
    )
```

Similarly, in R we can export the data to Seurat or SingleCellExperiment objects with the argument `X_layers` of the functions `get_seurat()` and `get_single_cell_experiment()`.

```r
library("cellxgene.census")
library("Seurat")

census <- open_soma(census_version = "latest")

organism <- "Homo sapiens"
gene_filter <- "feature_id %in% c('ENSG00000107317', 'ENSG00000106034')"
cell_filter <-  "cell_type == 'sympathetic neuron'"
cell_columns <- c("tissue", "cell_type")
layers <- c(data = "normalized")

seurat_obj <- get_seurat(
   census = census,
   organism = organism,
   var_value_filter = gene_filter,
   obs_value_filter = cell_filter,
   obs_column_names = cell_columns,
   X_layers = layers
)

#Single Cell Experiment

library("SingleCellExperiment")

sce_obj <- get_single_cell_experiment(
   census = census,
   organism = organism,
   var_value_filter = gene_filter,
   obs_value_filter = cell_filter,
   obs_column_names = cell_columns,
   X_layers = layers
)

```

### Accessing library-size normalized data layer via TileDB-SOMA

For memory-efficient data retrieval, you can use TileDB-SOMA as outlined below. In Python this looks like the following.

```python

import cellxgene_census
import tiledbsoma

# Open context manager
with cellxgene_census.open_soma(census_version = "latest") as census:

    # Access human SOMA object
    human = census["census_data"]["homo_sapiens"]

    query = human.axis_query(
       measurement_name = "RNA",
       obs_query = tiledbsoma.AxisQuery(
           value_filter = "tissue == 'brain' and sex == 'male'"
       )
    )

    # Set iterable for normalized matrix
    iterator = query.X("normalized").tables()
    
    # Iterate over the normalized matrix.
    # Get an iterative slice as pyarrow.Table
    raw_slice = next(iterator)
    
    # Perform analysis
    
    # close the query
    query.close()
```

And the equivalent code in R.

```r
library("cellxgene.census")
library("tiledbsoma")

human <-  census$get("census_data")$get("homo_sapiens")
query <-  human$axis_query(
  measurement_name = "RNA",
  obs_query = SOMAAxisQuery$new(
    value_filter = "tissue == 'brain' & sex == 'male'"
  )
)

# Set iterable for normalized matrix
iterator <-  query$X("normalized")$tables()

# Iterate over the normalized matrix.
# Get an iterative slice as an Arrow Table
raw_slice <-  iterator$read_next()

# Perform analysis
```

### Utilizing pre-calculated stats for querying `obs` and `var`

To filter cells or genes based on pre-calculated statistics and export to AnnData, you can use the new metadata variables as value filters.

For example, you can add a filter to query cells with more than 500 genes expressed, along with other filters. In Python this looks like the following.

```python
import cellxgene_census

with cellxgene_census.open_soma(census_version = "latest") as census:
    adata = cellxgene_census.get_anndata(
        census = census,
        organism = "Homo sapiens",
        obs_value_filter = "nnz > 500 and cell_type == 'sympathetic neuron'",
        column_names = {"obs": ["tissue", "cell_type"]},
        var_value_filter = "feature_id in ['ENSG00000161798', 'ENSG00000188229']",
    )

    print(adata)
```

In R, the equivalent code looks as follows.

```r
#Seurat
library("cellxgene.census")
library("Seurat")

census <- open_soma(census_version = "latest")

organism <- "Homo sapiens"
gene_filter <- "feature_id %in% c('ENSG00000107317', 'ENSG00000106034')"
cell_filter <-  "nnz > 500 & cell_type == 'sympathetic neuron'"
cell_columns <- c("tissue", "cell_type")

seurat_obj <- get_seurat(
   census = census,
   organism = organism,
   var_value_filter = gene_filter,
   obs_value_filter = cell_filter,
   obs_column_names = cell_columns
)

#Single Cell Experiment

library("SingleCellExperiment")

sce_obj <- get_single_cell_experiment(
   census = census,
   organism = organism,
   var_value_filter = gene_filter,
   obs_value_filter = cell_filter,
   obs_column_names = cell_columns
)
```

## Help us improve these data additions

We encourage you to engage with these new features in the Census API and share your feedback. This input is invaluable for the ongoing enhancement of the Census project.

For further information on any new feature, please reach out to us at [soma@chanzuckerberg.com](soma@chanzuckerberg.com). To report issues or for additional feedback, refer to our [Census GitHub repository](https://github.com/chanzuckerberg/cellxgene-census/issues).



# Section: cellxgene_census_docsite_schema

# Census data and schema

This page provides a user-friendly overview of the Census contents and its schema, in case you are interested you can find the full schema specification [here](https://github.com/chanzuckerberg/cellxgene-census/blob/main/docs/cellxgene_census_schema.md).

**Contents:**

1. [Schema](#schema)
2. [Data included in the Census](#data-included-in-the-census)
3. [SOMA objects](#soma-objects)

## Schema

The Census is a collection of a variety of **[SOMA objects](#soma-objects)** organized with the following hierarchy.

![image](cellxgene_census_docsite_model.svg)

As you can see the Census data is a `SOMACollection` with two high-level items:

1. `"census_info"` for the census summary info.
2. `"census_data"` for the single-cell data and metadata.

### Census summary info `"census_info"`

A `SOMAcollection` with tables providing information of the census as a whole, it has the following items:

- `"summary"`: high-level information of this Census, e.g. build date, total cell count, etc.
- `"datasets"`: A table with all datasets from CELLxGENE Discover used to create the Census.
- `"summary_cell_counts"`: Cell counts stratified by relevant cell metadata.

### Census single-cell data `"census_data"`

Data for each organism is stored in independent `SOMAExperiment` objects which are a specialized form of a `SOMACollection`. Each of these store a data matrix (cell by genes), cell metadata, gene metadata, and feature presence matrix:

This is how the data is organized for one organism – *Homo sapiens*:

- `["homo_sapiens"].obs`: Cell metadata.
- `["homo_sapiens"].ms["RNA"].X`: Data matrices: raw counts in `X["raw"]`, and library-size normalized counts in `X["normalized"]` (only avialble in Census schema V1.1.0 and above).
- `["homo_sapiens"].ms["RNA"].var`: Gene Metadata.
- `["homo_sapiens"].ms["RNA"]["feature_dataset_presence_matrix"]`: a sparse boolean array indicating which genes were measured in each dataset.

## Data included in the Census

All data from [CZ CELLxGENE Discover](https://cellxgene.cziscience.com/) that adheres to the following criteria is included in the Census:

- Cells from human or mouse.
- Non-spatial RNA data, see full list of sequencing technologies included [here](https://github.com/chanzuckerberg/cellxgene-census/blob/main/docs/cellxgene_census_schema.md#assays).
- Raw counts.
- Only standardized cell and gene metadata as described in the CELLxGENE Discover dataset [schema](https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/3.0.0/schema.md).

⚠️ Note that the data includes:

- **Full-gene sequencing read counts** (e.g. Smart-Seq2) and **molecule counts** (e.g. 10X).
- **Duplicate cells** present across multiple datasets, these can be filtered in or out using the cell metadata variable `is_primary_data`.

## SOMA objects

You can find the full SOMA specification [here](https://github.com/single-cell-data/SOMA/blob/main/abstract_specification.md#foundational-types).

The following is short description of the main SOMA objects used by the Census:

- `DenseNDArray` is a dense, N-dimensional array, with offset (zero-based) integer indexing on each dimension.
- `SparseNDArray` is the same as `DenseNDArray` but sparse, and supports point indexing (disjoint index access).
- `DataFrame` is a multi-column table with a user-defined columns names and value types, with support for point indexing.
- `Collection` is a persistent container of named SOMA objects.
- `Experiment` is a class that represents a single-cell experiment. It always contains two objects:
  - `obs`: a  `DataFrame` with primary annotations on the observation axis.
  - `ms`: a  `Collection` of measurements, each composed of `X` matrices and axis annotation matrices or data frames (e.g. `var`, `varm`, `obsm`, etc).



# Section: articles-2023-20230919-out_of_core_methods

# Memory-efficient implementations of commonly used single-cell methods

*Published:* *September 18, 2023*

*By:* *[Pablo Garcia-Nieto](mailto:pgarcia-nieto@chanzuckerberg.com)*

The Census team is thrilled to officially announce memory-efficient implementations of some of the most widely used single-cell algorithms.

With just a few lines of code, using the Census Python API, users can now perform the following computational tasks across tens of millions of cells using a conventional laptop with 8GB of memory:

* Calculating average and variance gene expression for cells or genes. See example below or the full [tutorial here](../../notebooks/experimental/mean_variance.ipynb).
* Obtaining batch-corrected highly variable genes.  See example below or the full [tutorial here](../../notebooks/experimental/highly_variable_genes).

These implementations are interwoven with the way users query slices of Census data, which means that these tasks can be seamlessly applied to any slice of the 33M+ cells available in Census.

Continue reading for more implementation details and usage examples.

## Efficient calculation of average and variance gene expression across millions of cells

With `cellxgene_census.experimental.pp.mean_variance` users can now get gene expression average and variance for all genes or cells in a given Census query.

### How it works

Calculations are done in an accumulative incremental fashion, meaning that only a small fraction of the total data is processed at any given time.

The Census data is downloaded in increments and average and variance accumulators are updated at each incremental step. The implementation also takes advantage of CPU-based multiprocessing to speed up the process.

Currently, the mean and variance are calculated using the full population of cells/genes, including those with a zero valued measurement. In the future, we will enable calculation of mean including only the population of non-zero cells/genes.

### Example: *KRAS* and *AQP4* average and variance expression in lung epithelial cells

The following calculates the average and variance values for the genes *KRAS* and *AQP4* in all epithelial cells of the human lung.

Users can easily switch the calculation, and obtain average and variance for each cell across the genes in the query. This is controlled by the `axis` argument of `mean_variance`.

```python
import cellxgene_census
import tiledbsoma as soma
from cellxgene_census.experimental.pp import mean_variance

# Open the Census
census = cellxgene_census.open_soma()
human_data = census["census_data"]["homo_sapiens"]

# Set filters
cell_filter = (
  "is_primary_data == True "
  "and tissue_general == 'lung' "
  "and cell_type == 'epithelial cell'"
 )
gene_filter = "feature_name in ['KRAS', 'AQP4']"

# Perform query
query = human_data.axis_query(
  measurement_name="RNA",
  obs_query=soma.AxisQuery(value_filter= cell_filter),
  var_query=soma.AxisQuery(value_filter= gene_filter)
)

# Calculate mean and average per gene
mean_variance_df = mean_variance(query, axis=0, calculate_mean=True, calculate_variance=True)

# Get gene metadata of query
gene_df = query.var().concat().to_pandas()

query.close()
census.close()
```

Which results in:

```python
mean_variance_df
# soma_joinid      mean     variance
# 8624         3.071926  5741.242485
# 16437        8.233282   452.119153

gene_df
#    soma_joinid       feature_id feature_name  feature_length
# 0         8624  ENSG00000171885         AQP4            5943
# 1        16437  ENSG00000133703         KRAS            6845
```

## Efficient calculation of highly variable genes across millions of cells

With `cellxgene_census.experimental.pp.get_highly_variable_genes` users can get the most highly variable genes of a Census query while accounting for batch effects.

This is usually the first pre-processing step necessary for other downstream tasks, for example data integration.

### How it works

The Census algorithm is based on the scanpy method `scanpy.pp.highly_variable_genes`, and in particular the Seurat V3 method, which is designed for raw counts and can account for batch effects.

The Census implementation utilizes the same incremental paradigm used in `cellxgene_census.experimental.pp.mean_variance` (see above), calculating incremental-based mean and variance accumulators with some tweaks to comply to the Seurat V3 method.

### Example: Finding highly variable genes for all cells of the human esophagus

The following example identifies the top 1000 highly variable genes for all human esophagus cells. As a general rule of thumb it is good to use `dataset_id` as the batch variable.

```python
import cellxgene_census
from cellxgene_census.experimental.pp import get_highly_variable_genes

census = cellxgene_census.open_soma()

hvg = get_highly_variable_genes(
  census,
  organism="Homo sapiens",
  obs_value_filter="is_primary_data == True and tissue_general == 'esophagus'",
  n_top_genes = 1000,
  batch_key = "dataset_id"
)

census.close()
```

Which results in:

```python
hvg
# soma_joinid    means  variances  ...  variances_norm  highly_variable
# 0            0.003692   0.004627  ...        0.748221            False
# 1            0.003084   0.003203  ...        0.898657            False
# 2            0.014962   0.037395  ...        0.513473            False
# 3            0.218865   1.547648  ...        4.786928             True
# 4            0.002142   0.002242  ...        0.894955            False
# ...               ...        ...  ...             ...              ...
# 60659        0.000000   0.000000  ...        0.000000            False
# 60660        0.000000   0.000000  ...        0.000000            False
# 60661        0.000000   0.000000  ...        0.000000            False
# 60662        0.000000   0.000000  ...        0.000000            False
# 60663        0.000000   0.000000  ...        0.000000            False
```



# Section: cellxgene_census-tests-README

# Test README

This directory contains tests of the cellxgene-census package API, _and_ the use of the API on the
live "corpus", i.e., data in the public Census S3 bucket. The tests use Pytest, and have
Pytest marks to control which tests are run.

Tests can be run in the usual manner. First, ensure you have cellxgene-census installed, e.g., from the top-level repo directory:

> pip install -e ./api/python/cellxgene_census/

Then run the tests:

> pytest ./api/python/cellxgene_census/

## Pytest Marks

There are various Pytest marks you can use from the command line:

- live_corpus: tests that directly access the `latest` version of the Census. Enabled by default.
- expensive: tests that are expensive (ie., cpu, memory, time). Disabled by default - enable with `--expensive`. Some of these tests are _very_ expensive, ie., require a very large memory host to succeed.
- experimental: tests that are for code in the `experimental` package. Disabled by default - enable with `--experimental`. These tests require installation the optional Python packages installed via pip `pip install -e ./api/python/cellxgene_census/[experimental]`

By default, only relatively cheap & fast tests are run. To enable `expensive` tests:

> pytest --expensive ...

To enable `experimental` tests:

> pytest --experimental ...

To disable `live_corpus` tests:

> pytest -m 'not live_corpus'

You can also combine them, e.g.,

> pytest -m 'not live_corpus' --expensive --experimental

## Acceptance (expensive) tests

These tests are periodically run, and are not part of CI due to their overhead.

When run, please record the results in this file (below) and commit the change to git. Please include the following information:

- date
- config:
  - EC2 instance type and any system config (i.e., swap)
  - host and OS as reported by `uname -a`. **Please remove IP address**
  - Python & package versions and OS - suggest capturing the output of `tiledbsoma.show_package_versions()`
  - The Census version used for the test (i.e., the version aliased as `latest`). This can be easily captured using `cellxgene_census.get_census_version_description('latest')`
  - the cellxgene_census package version (ie., `cellxgene_census.__version__`)
- any run notes
- full output of: `pytest -v --durations=0 --experimental --expensive ./api/python/cellxgene_census/tests/`

### 2023-10-23

- Host: EC2 instance type: `r6id.32xlarge`, all nvme mounted as swap.
- Uname: Linux 6.2.0-1015-aws #15~22.04.1-Ubuntu SMP Fri Oct  6 21:37:24 UTC 2023 x86_64 x86_64 x86_64 GNU/Linux
- Python & census versions:

```python
>>> import cellxgene_census, tiledbsoma
>>> tiledbsoma.show_package_versions()
tiledbsoma.__version__        1.5.1
TileDB-Py tiledb.version()    (0, 23, 4)
TileDB core version           2.17.4
libtiledbsoma version()       libtiledb=2.17.4
python version                3.10.12.final.0
OS version                    Linux 6.2.0-1015-aws
```

**Pytest output:**

```text
============================= test session starts ==============================
platform linux -- Python 3.10.12, pytest-7.4.3, pluggy-1.3.0 -- /usr/bin/python3
cachedir: .pytest_cache
rootdir: /home/ssm-user/cellxgene-census/api/python/cellxgene_census
configfile: pyproject.toml
plugins: requests-mock-1.11.0
collecting ... collected 411 items

api/python/cellxgene_census/tests/test_acceptance.py::test_load_axes[homo_sapiens] PASSED [  0%]
api/python/cellxgene_census/tests/test_acceptance.py::test_load_axes[mus_musculus] PASSED [  0%]
api/python/cellxgene_census/tests/test_acceptance.py::test_incremental_read_obs[2-None-homo_sapiens] PASSED [  0%]
api/python/cellxgene_census/tests/test_acceptance.py::test_incremental_read_obs[2-None-mus_musculus] PASSED [  0%]
api/python/cellxgene_census/tests/test_acceptance.py::test_incremental_read_obs[None-ctx_config1-homo_sapiens] PASSED [  1%]
api/python/cellxgene_census/tests/test_acceptance.py::test_incremental_read_obs[None-ctx_config1-mus_musculus] PASSED [  1%]
api/python/cellxgene_census/tests/test_acceptance.py::test_incremental_read_var[2-None-homo_sapiens] PASSED [  1%]
api/python/cellxgene_census/tests/test_acceptance.py::test_incremental_read_var[2-None-mus_musculus] PASSED [  1%]
api/python/cellxgene_census/tests/test_acceptance.py::test_incremental_read_var[None-ctx_config1-homo_sapiens] PASSED [  2%]
api/python/cellxgene_census/tests/test_acceptance.py::test_incremental_read_var[None-ctx_config1-mus_musculus] PASSED [  2%]
api/python/cellxgene_census/tests/test_acceptance.py::test_incremental_read_X[2-None-coords0-homo_sapiens] PASSED [  2%]
api/python/cellxgene_census/tests/test_acceptance.py::test_incremental_read_X[2-None-coords0-mus_musculus] PASSED [  2%]
api/python/cellxgene_census/tests/test_acceptance.py::test_incremental_read_X[None-ctx_config1-coords1-homo_sapiens] PASSED [  3%]
api/python/cellxgene_census/tests/test_acceptance.py::test_incremental_read_X[None-ctx_config1-coords1-mus_musculus] PASSED [  3%]
api/python/cellxgene_census/tests/test_acceptance.py::test_incremental_read_X[None-ctx_config2-coords2-homo_sapiens] PASSED [  3%]
api/python/cellxgene_census/tests/test_acceptance.py::test_incremental_read_X[None-ctx_config2-coords2-mus_musculus] PASSED [  3%]
api/python/cellxgene_census/tests/test_acceptance.py::test_incremental_read_X[None-ctx_config3-coords3-homo_sapiens] PASSED [  4%]
api/python/cellxgene_census/tests/test_acceptance.py::test_incremental_read_X[None-ctx_config3-coords3-mus_musculus] PASSED [  4%]
api/python/cellxgene_census/tests/test_acceptance.py::test_incremental_query[2-tissue=='aorta'-homo_sapiens] PASSED [  4%]
api/python/cellxgene_census/tests/test_acceptance.py::test_incremental_query[2-tissue=='aorta'-mus_musculus] PASSED [  4%]
api/python/cellxgene_census/tests/test_acceptance.py::test_incremental_query[2-tissue=='brain'-homo_sapiens] PASSED [  5%]
api/python/cellxgene_census/tests/test_acceptance.py::test_incremental_query[2-tissue=='brain'-mus_musculus] PASSED [  5%]
api/python/cellxgene_census/tests/test_acceptance.py::test_incremental_query[None-tissue=='aorta'-homo_sapiens] PASSED [  5%]
api/python/cellxgene_census/tests/test_acceptance.py::test_incremental_query[None-tissue=='aorta'-mus_musculus] PASSED [  5%]
api/python/cellxgene_census/tests/test_acceptance.py::test_incremental_query[None-tissue=='brain'-homo_sapiens] PASSED [  6%]
api/python/cellxgene_census/tests/test_acceptance.py::test_incremental_query[None-tissue=='brain'-mus_musculus] PASSED [  6%]
api/python/cellxgene_census/tests/test_acceptance.py::test_get_anndata[tissue=='aorta'-None-ctx_config0-homo_sapiens] PASSED [  6%]
api/python/cellxgene_census/tests/test_acceptance.py::test_get_anndata[tissue=='aorta'-None-ctx_config0-mus_musculus] PASSED [  6%]
api/python/cellxgene_census/tests/test_acceptance.py::test_get_anndata[First 10K cells-homo_sapiens] PASSED [  7%]
api/python/cellxgene_census/tests/test_acceptance.py::test_get_anndata[First 10K cells-mus_musculus] PASSED [  7%]
api/python/cellxgene_census/tests/test_acceptance.py::test_get_anndata[First 100K cells-homo_sapiens] PASSED [  7%]
api/python/cellxgene_census/tests/test_acceptance.py::test_get_anndata[First 100K cells-mus_musculus] PASSED [  7%]
api/python/cellxgene_census/tests/test_acceptance.py::test_get_anndata[First 250K cells-homo_sapiens] PASSED [  8%]
api/python/cellxgene_census/tests/test_acceptance.py::test_get_anndata[First 250K cells-mus_musculus] PASSED [  8%]
api/python/cellxgene_census/tests/test_acceptance.py::test_get_anndata[First 500K cells-homo_sapiens] PASSED [  8%]
api/python/cellxgene_census/tests/test_acceptance.py::test_get_anndata[First 500K cells-mus_musculus] PASSED [  8%]
api/python/cellxgene_census/tests/test_acceptance.py::test_get_anndata[First 750K cells-homo_sapiens] PASSED [  9%]
api/python/cellxgene_census/tests/test_acceptance.py::test_get_anndata[First 750K cells-mus_musculus] PASSED [  9%]
api/python/cellxgene_census/tests/test_acceptance.py::test_get_anndata[First 1M cells-homo_sapiens] PASSED [  9%]
api/python/cellxgene_census/tests/test_acceptance.py::test_get_anndata[First 1M cells-mus_musculus] PASSED [  9%]
api/python/cellxgene_census/tests/test_acceptance.py::test_get_anndata[tissue_general=='brain'-obs_coords7-ctx_config7-homo_sapiens] PASSED [  9%]
api/python/cellxgene_census/tests/test_acceptance.py::test_get_anndata[tissue_general=='brain'-obs_coords7-ctx_config7-mus_musculus] PASSED [ 10%]
api/python/cellxgene_census/tests/test_acceptance.py::test_get_anndata[cell_type=='neuron'-obs_coords8-ctx_config8-homo_sapiens] PASSED [ 10%]
api/python/cellxgene_census/tests/test_acceptance.py::test_get_anndata[cell_type=='neuron'-obs_coords8-ctx_config8-mus_musculus] PASSED [ 10%]
api/python/cellxgene_census/tests/test_acceptance.py::test_get_anndata[is_primary_data==True-obs_coords9-ctx_config9-homo_sapiens] PASSED [ 10%]
api/python/cellxgene_census/tests/test_acceptance.py::test_get_anndata[is_primary_data==True-obs_coords9-ctx_config9-mus_musculus] PASSED [ 11%]
api/python/cellxgene_census/tests/test_acceptance.py::test_get_anndata[None-obs_coords10-ctx_config10-homo_sapiens] PASSED [ 11%]
api/python/cellxgene_census/tests/test_acceptance.py::test_get_anndata[None-obs_coords10-ctx_config10-mus_musculus] PASSED [ 11%]
api/python/cellxgene_census/tests/test_directory.py::test_get_census_version_directory PASSED [ 11%]
api/python/cellxgene_census/tests/test_directory.py::test_get_census_version_directory__lts_only PASSED [ 12%]
api/python/cellxgene_census/tests/test_directory.py::test_get_census_version_directory__exclude_lts PASSED [ 12%]
api/python/cellxgene_census/tests/test_directory.py::test_get_census_version_directory__include_retracted PASSED [ 12%]
api/python/cellxgene_census/tests/test_directory.py::test_get_census_version_directory__retraction_info PASSED [ 12%]
api/python/cellxgene_census/tests/test_directory.py::test_get_census_version_description_errors PASSED [ 13%]
api/python/cellxgene_census/tests/test_directory.py::test_get_census_mirrors_directory PASSED [ 13%]
api/python/cellxgene_census/tests/test_directory.py::test_live_directory_contents PASSED [ 13%]
api/python/cellxgene_census/tests/test_get_anndata.py::test_get_anndata_value_filter PASSED [ 13%]
api/python/cellxgene_census/tests/test_get_anndata.py::test_get_anndata_coords PASSED [ 14%]
api/python/cellxgene_census/tests/test_get_anndata.py::test_get_anndata_allows_missing_obs_or_var_filter PASSED [ 14%]
api/python/cellxgene_census/tests/test_get_anndata.py::test_get_anndata_x_layer[raw] PASSED [ 14%]
api/python/cellxgene_census/tests/test_get_anndata.py::test_get_anndata_x_layer[normalized] PASSED [ 14%]
api/python/cellxgene_census/tests/test_get_anndata.py::test_get_anndata_two_layers[layers0] PASSED [ 15%]
api/python/cellxgene_census/tests/test_get_anndata.py::test_get_anndata_two_layers[layers1] PASSED [ 15%]
api/python/cellxgene_census/tests/test_get_anndata.py::test_get_anndata_wrong_layer_names PASSED [ 15%]
api/python/cellxgene_census/tests/test_get_helpers.py::test_get_experiment PASSED [ 15%]
api/python/cellxgene_census/tests/test_get_helpers.py::test_get_presence_matrix[homo_sapiens] PASSED [ 16%]
api/python/cellxgene_census/tests/test_get_helpers.py::test_get_presence_matrix[mus_musculus] PASSED [ 16%]
api/python/cellxgene_census/tests/test_lts_compat.py::test_open[stable] PASSED [ 16%]
api/python/cellxgene_census/tests/test_lts_compat.py::test_read_dataframe[stable] PASSED [ 16%]
api/python/cellxgene_census/tests/test_lts_compat.py::test_read_arrays[stable] PASSED [ 17%]
api/python/cellxgene_census/tests/test_open.py::test_open_soma_stable PASSED [ 17%]
api/python/cellxgene_census/tests/test_open.py::test_open_soma_latest PASSED [ 17%]
api/python/cellxgene_census/tests/test_open.py::test_open_soma_with_context FAILED [ 17%]
api/python/cellxgene_census/tests/test_open.py::test_open_soma_invalid_args PASSED [ 18%]
api/python/cellxgene_census/tests/test_open.py::test_open_soma_errors PASSED [ 18%]
api/python/cellxgene_census/tests/test_open.py::test_open_soma_uses_correct_mirror PASSED [ 18%]
api/python/cellxgene_census/tests/test_open.py::test_open_soma_rejects_non_s3_mirror PASSED [ 18%]
api/python/cellxgene_census/tests/test_open.py::test_open_soma_works_if_no_relative_uri_specified PASSED [ 18%]
api/python/cellxgene_census/tests/test_open.py::test_open_soma_defaults_to_stable PASSED [ 19%]
api/python/cellxgene_census/tests/test_open.py::test_get_source_h5ad_uri PASSED [ 19%]
api/python/cellxgene_census/tests/test_open.py::test_get_source_h5ad_uri_errors PASSED [ 19%]
api/python/cellxgene_census/tests/test_open.py::test_download_source_h5ad PASSED [ 19%]
api/python/cellxgene_census/tests/test_open.py::test_download_source_h5ad_errors PASSED [ 20%]
api/python/cellxgene_census/tests/test_open.py::test_opening_census_without_anon_access_fails_with_bogus_creds PASSED [ 20%]
api/python/cellxgene_census/tests/test_open.py::test_can_open_with_anonymous_access PASSED [ 20%]
api/python/cellxgene_census/tests/test_util.py::test_uri_join PASSED     [ 20%]
api/python/cellxgene_census/tests/experimental/ml/test_pytorch.py::test_non_batched[6-3-pytorch_x_value_gen-True] PASSED [ 21%]
api/python/cellxgene_census/tests/experimental/ml/test_pytorch.py::test_non_batched[6-3-pytorch_x_value_gen-False] PASSED [ 21%]
api/python/cellxgene_census/tests/experimental/ml/test_pytorch.py::test_batching__all_batches_full_size[6-3-pytorch_x_value_gen-True] PASSED [ 21%]
api/python/cellxgene_census/tests/experimental/ml/test_pytorch.py::test_batching__all_batches_full_size[6-3-pytorch_x_value_gen-False] PASSED [ 21%]
api/python/cellxgene_census/tests/experimental/ml/test_pytorch.py::test_unique_soma_joinids[obs_range0-3-pytorch_x_value_gen-True] PASSED [ 22%]
api/python/cellxgene_census/tests/experimental/ml/test_pytorch.py::test_unique_soma_joinids[obs_range1-3-pytorch_x_value_gen-False] PASSED [ 22%]
api/python/cellxgene_census/tests/experimental/ml/test_pytorch.py::test_batching__partial_final_batch_size[5-3-pytorch_x_value_gen-True] PASSED [ 22%]
api/python/cellxgene_census/tests/experimental/ml/test_pytorch.py::test_batching__partial_final_batch_size[5-3-pytorch_x_value_gen-False] PASSED [ 22%]
api/python/cellxgene_census/tests/experimental/ml/test_pytorch.py::test_batching__exactly_one_batch[3-3-pytorch_x_value_gen-True] PASSED [ 23%]
api/python/cellxgene_census/tests/experimental/ml/test_pytorch.py::test_batching__exactly_one_batch[3-3-pytorch_x_value_gen-False] PASSED [ 23%]
api/python/cellxgene_census/tests/experimental/ml/test_pytorch.py::test_batching__empty_query_result[6-3-pytorch_x_value_gen-True] PASSED [ 23%]
api/python/cellxgene_census/tests/experimental/ml/test_pytorch.py::test_batching__empty_query_result[6-3-pytorch_x_value_gen-False] PASSED [ 23%]
api/python/cellxgene_census/tests/experimental/ml/test_pytorch.py::test_sparse_output__non_batched[6-3-pytorch_x_value_gen-True] PASSED [ 24%]
api/python/cellxgene_census/tests/experimental/ml/test_pytorch.py::test_sparse_output__non_batched[6-3-pytorch_x_value_gen-False] PASSED [ 24%]
api/python/cellxgene_census/tests/experimental/ml/test_pytorch.py::test_sparse_output__batched[6-3-pytorch_x_value_gen-True] PASSED [ 24%]
api/python/cellxgene_census/tests/experimental/ml/test_pytorch.py::test_sparse_output__batched[6-3-pytorch_x_value_gen-False] PASSED [ 24%]
api/python/cellxgene_census/tests/experimental/ml/test_pytorch.py::test_batching__partial_soma_batches_are_concatenated[10-1-pytorch_x_value_gen-True] PASSED [ 25%]
api/python/cellxgene_census/tests/experimental/ml/test_pytorch.py::test_batching__partial_soma_batches_are_concatenated[10-1-pytorch_x_value_gen-False] PASSED [ 25%]
api/python/cellxgene_census/tests/experimental/ml/test_pytorch.py::test_encoders[3-3-pytorch_x_value_gen] PASSED [ 25%]
api/python/cellxgene_census/tests/experimental/ml/test_pytorch.py::test_multiprocessing__returns_full_result[6-3-pytorch_x_value_gen] PASSED [ 25%]
api/python/cellxgene_census/tests/experimental/ml/test_pytorch.py::test_distributed__returns_data_partition_for_rank[6-3-pytorch_x_value_gen] PASSED [ 26%]
api/python/cellxgene_census/tests/experimental/ml/test_pytorch.py::test_distributed_and_multiprocessing__returns_data_partition_for_rank[12-3-pytorch_x_value_gen] PASSED [ 26%]
api/python/cellxgene_census/tests/experimental/ml/test_pytorch.py::test_experiment_dataloader__non_batched[3-3-pytorch_x_value_gen-True] PASSED [ 26%]
api/python/cellxgene_census/tests/experimental/ml/test_pytorch.py::test_experiment_dataloader__non_batched[3-3-pytorch_x_value_gen-False] PASSED [ 26%]
api/python/cellxgene_census/tests/experimental/ml/test_pytorch.py::test_experiment_dataloader__batched[6-3-pytorch_x_value_gen-True] PASSED [ 27%]
api/python/cellxgene_census/tests/experimental/ml/test_pytorch.py::test_experiment_dataloader__batched[6-3-pytorch_x_value_gen-False] PASSED [ 27%]
api/python/cellxgene_census/tests/experimental/ml/test_pytorch.py::test__X_tensor_dtype_matches_X_matrix[6-3-pytorch_x_value_gen-True] PASSED [ 27%]
api/python/cellxgene_census/tests/experimental/ml/test_pytorch.py::test__X_tensor_dtype_matches_X_matrix[6-3-pytorch_x_value_gen-False] PASSED [ 27%]
api/python/cellxgene_census/tests/experimental/ml/test_pytorch.py::test__pytorch_splitting[10-1-pytorch_x_value_gen] PASSED [ 27%]
api/python/cellxgene_census/tests/experimental/ml/test_pytorch.py::test__shuffle[16-1-pytorch_seq_x_value_gen] PASSED [ 28%]
api/python/cellxgene_census/tests/experimental/ml/test_pytorch.py::test_experiment_dataloader__multiprocess_sparse_matrix__fails SKIPPED [ 28%]
api/python/cellxgene_census/tests/experimental/ml/test_pytorch.py::test_experiment_dataloader__multiprocess_dense_matrix__ok SKIPPED [ 28%]
api/python/cellxgene_census/tests/experimental/ml/test_pytorch.py::test_experiment_dataloader__unsupported_params__fails PASSED [ 28%]
api/python/cellxgene_census/tests/experimental/ml/huggingface/test_geneformer.py::test_GeneformerTokenizer[4] SKIPPED [ 29%]
api/python/cellxgene_census/tests/experimental/ml/huggingface/test_geneformer.py::test_GeneformerTokenizer[100000] SKIPPED [ 29%]
api/python/cellxgene_census/tests/experimental/ml/huggingface/test_geneformer.py::test_GeneformerTokenizer_docstring_example SKIPPED [ 29%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-None-None-50-mus_musculus-is_primary_data == True and tissue_general == "liver"] PASSED [ 29%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-None-None-50-mus_musculus-is_primary_data == True and tissue_general == "skin of body"] PASSED [ 30%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-None-None-50-mus_musculus-is_primary_data == True and tissue_general in ["heart", "lung"]] PASSED [ 30%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-None-None-50-mus_musculus-is_primary_data == True and assay == "Smart-seq2"] PASSED [ 30%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-None-None-500-mus_musculus-is_primary_data == True and tissue_general == "liver"] PASSED [ 30%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-None-None-500-mus_musculus-is_primary_data == True and tissue_general == "skin of body"] PASSED [ 31%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-None-None-500-mus_musculus-is_primary_data == True and tissue_general in ["heart", "lung"]] PASSED [ 31%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-None-None-500-mus_musculus-is_primary_data == True and assay == "Smart-seq2"] PASSED [ 31%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-None-dataset_id-50-mus_musculus-is_primary_data == True and tissue_general == "liver"] PASSED [ 31%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-None-dataset_id-50-mus_musculus-is_primary_data == True and tissue_general == "skin of body"] PASSED [ 32%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-None-dataset_id-50-mus_musculus-is_primary_data == True and tissue_general in ["heart", "lung"]] PASSED [ 32%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-None-dataset_id-50-mus_musculus-is_primary_data == True and assay == "Smart-seq2"] PASSED [ 32%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-None-dataset_id-500-mus_musculus-is_primary_data == True and tissue_general == "liver"] PASSED [ 32%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-None-dataset_id-500-mus_musculus-is_primary_data == True and tissue_general == "skin of body"] PASSED [ 33%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-None-dataset_id-500-mus_musculus-is_primary_data == True and tissue_general in ["heart", "lung"]] PASSED [ 33%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-None-dataset_id-500-mus_musculus-is_primary_data == True and assay == "Smart-seq2"] PASSED [ 33%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-None-batch_key2-50-mus_musculus-is_primary_data == True and tissue_general == "liver"] PASSED [ 33%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-None-batch_key2-50-mus_musculus-is_primary_data == True and tissue_general == "skin of body"] PASSED [ 34%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-None-batch_key2-50-mus_musculus-is_primary_data == True and tissue_general in ["heart", "lung"]] PASSED [ 34%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-None-batch_key2-50-mus_musculus-is_primary_data == True and assay == "Smart-seq2"] PASSED [ 34%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-None-batch_key2-500-mus_musculus-is_primary_data == True and tissue_general == "liver"] PASSED [ 34%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-None-batch_key2-500-mus_musculus-is_primary_data == True and tissue_general == "skin of body"] PASSED [ 35%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-None-batch_key2-500-mus_musculus-is_primary_data == True and tissue_general in ["heart", "lung"]] PASSED [ 35%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-None-batch_key2-500-mus_musculus-is_primary_data == True and assay == "Smart-seq2"] PASSED [ 35%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-None-batch_key3-50-mus_musculus-is_primary_data == True and tissue_general == "liver"] PASSED [ 35%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-None-batch_key3-50-mus_musculus-is_primary_data == True and tissue_general == "skin of body"] PASSED [ 36%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-None-batch_key3-50-mus_musculus-is_primary_data == True and tissue_general in ["heart", "lung"]] PASSED [ 36%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-None-batch_key3-50-mus_musculus-is_primary_data == True and assay == "Smart-seq2"] PASSED [ 36%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-None-batch_key3-500-mus_musculus-is_primary_data == True and tissue_general == "liver"] PASSED [ 36%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-None-batch_key3-500-mus_musculus-is_primary_data == True and tissue_general == "skin of body"] PASSED [ 36%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-None-batch_key3-500-mus_musculus-is_primary_data == True and tissue_general in ["heart", "lung"]] PASSED [ 37%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-None-batch_key3-500-mus_musculus-is_primary_data == True and assay == "Smart-seq2"] PASSED [ 37%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-0.5-None-50-mus_musculus-is_primary_data == True and tissue_general == "liver"] PASSED [ 37%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-0.5-None-50-mus_musculus-is_primary_data == True and tissue_general == "skin of body"] PASSED [ 37%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-0.5-None-50-mus_musculus-is_primary_data == True and tissue_general in ["heart", "lung"]] PASSED [ 38%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-0.5-None-50-mus_musculus-is_primary_data == True and assay == "Smart-seq2"] PASSED [ 38%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-0.5-None-500-mus_musculus-is_primary_data == True and tissue_general == "liver"] PASSED [ 38%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-0.5-None-500-mus_musculus-is_primary_data == True and tissue_general == "skin of body"] PASSED [ 38%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-0.5-None-500-mus_musculus-is_primary_data == True and tissue_general in ["heart", "lung"]] PASSED [ 39%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-0.5-None-500-mus_musculus-is_primary_data == True and assay == "Smart-seq2"] PASSED [ 39%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-0.5-dataset_id-50-mus_musculus-is_primary_data == True and tissue_general == "liver"] PASSED [ 39%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-0.5-dataset_id-50-mus_musculus-is_primary_data == True and tissue_general == "skin of body"] PASSED [ 39%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-0.5-dataset_id-50-mus_musculus-is_primary_data == True and tissue_general in ["heart", "lung"]] PASSED [ 40%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-0.5-dataset_id-50-mus_musculus-is_primary_data == True and assay == "Smart-seq2"] PASSED [ 40%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-0.5-dataset_id-500-mus_musculus-is_primary_data == True and tissue_general == "liver"] PASSED [ 40%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-0.5-dataset_id-500-mus_musculus-is_primary_data == True and tissue_general == "skin of body"] PASSED [ 40%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-0.5-dataset_id-500-mus_musculus-is_primary_data == True and tissue_general in ["heart", "lung"]] PASSED [ 41%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-0.5-dataset_id-500-mus_musculus-is_primary_data == True and assay == "Smart-seq2"] PASSED [ 41%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-0.5-batch_key2-50-mus_musculus-is_primary_data == True and tissue_general == "liver"] PASSED [ 41%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-0.5-batch_key2-50-mus_musculus-is_primary_data == True and tissue_general == "skin of body"] PASSED [ 41%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-0.5-batch_key2-50-mus_musculus-is_primary_data == True and tissue_general in ["heart", "lung"]] PASSED [ 42%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-0.5-batch_key2-50-mus_musculus-is_primary_data == True and assay == "Smart-seq2"] PASSED [ 42%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-0.5-batch_key2-500-mus_musculus-is_primary_data == True and tissue_general == "liver"] PASSED [ 42%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-0.5-batch_key2-500-mus_musculus-is_primary_data == True and tissue_general == "skin of body"] PASSED [ 42%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-0.5-batch_key2-500-mus_musculus-is_primary_data == True and tissue_general in ["heart", "lung"]] PASSED [ 43%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-0.5-batch_key2-500-mus_musculus-is_primary_data == True and assay == "Smart-seq2"] PASSED [ 43%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-0.5-batch_key3-50-mus_musculus-is_primary_data == True and tissue_general == "liver"] PASSED [ 43%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-0.5-batch_key3-50-mus_musculus-is_primary_data == True and tissue_general == "skin of body"] PASSED [ 43%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-0.5-batch_key3-50-mus_musculus-is_primary_data == True and tissue_general in ["heart", "lung"]] PASSED [ 44%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-0.5-batch_key3-50-mus_musculus-is_primary_data == True and assay == "Smart-seq2"] PASSED [ 44%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-0.5-batch_key3-500-mus_musculus-is_primary_data == True and tissue_general == "liver"] PASSED [ 44%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-0.5-batch_key3-500-mus_musculus-is_primary_data == True and tissue_general == "skin of body"] PASSED [ 44%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-0.5-batch_key3-500-mus_musculus-is_primary_data == True and tissue_general in ["heart", "lung"]] PASSED [ 45%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-0.5-batch_key3-500-mus_musculus-is_primary_data == True and assay == "Smart-seq2"] PASSED [ 45%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-None-None-50-mus_musculus-is_primary_data == True and tissue_general == "liver"] PASSED [ 45%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-None-None-50-mus_musculus-is_primary_data == True and tissue_general == "skin of body"] PASSED [ 45%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-None-None-50-mus_musculus-is_primary_data == True and tissue_general in ["heart", "lung"]] PASSED [ 45%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-None-None-50-mus_musculus-is_primary_data == True and assay == "Smart-seq2"] PASSED [ 46%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-None-None-500-mus_musculus-is_primary_data == True and tissue_general == "liver"] PASSED [ 46%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-None-None-500-mus_musculus-is_primary_data == True and tissue_general == "skin of body"] PASSED [ 46%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-None-None-500-mus_musculus-is_primary_data == True and tissue_general in ["heart", "lung"]] PASSED [ 46%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-None-None-500-mus_musculus-is_primary_data == True and assay == "Smart-seq2"] PASSED [ 47%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-None-dataset_id-50-mus_musculus-is_primary_data == True and tissue_general == "liver"] PASSED [ 47%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-None-dataset_id-50-mus_musculus-is_primary_data == True and tissue_general == "skin of body"] PASSED [ 47%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-None-dataset_id-50-mus_musculus-is_primary_data == True and tissue_general in ["heart", "lung"]] PASSED [ 47%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-None-dataset_id-50-mus_musculus-is_primary_data == True and assay == "Smart-seq2"] PASSED [ 48%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-None-dataset_id-500-mus_musculus-is_primary_data == True and tissue_general == "liver"] PASSED [ 48%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-None-dataset_id-500-mus_musculus-is_primary_data == True and tissue_general == "skin of body"] PASSED [ 48%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-None-dataset_id-500-mus_musculus-is_primary_data == True and tissue_general in ["heart", "lung"]] PASSED [ 48%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-None-dataset_id-500-mus_musculus-is_primary_data == True and assay == "Smart-seq2"] PASSED [ 49%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-None-batch_key2-50-mus_musculus-is_primary_data == True and tissue_general == "liver"] PASSED [ 49%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-None-batch_key2-50-mus_musculus-is_primary_data == True and tissue_general == "skin of body"] PASSED [ 49%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-None-batch_key2-50-mus_musculus-is_primary_data == True and tissue_general in ["heart", "lung"]] PASSED [ 49%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-None-batch_key2-50-mus_musculus-is_primary_data == True and assay == "Smart-seq2"] PASSED [ 50%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-None-batch_key2-500-mus_musculus-is_primary_data == True and tissue_general == "liver"] PASSED [ 50%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-None-batch_key2-500-mus_musculus-is_primary_data == True and tissue_general == "skin of body"] PASSED [ 50%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-None-batch_key2-500-mus_musculus-is_primary_data == True and tissue_general in ["heart", "lung"]] PASSED [ 50%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-None-batch_key2-500-mus_musculus-is_primary_data == True and assay == "Smart-seq2"] PASSED [ 51%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-None-batch_key3-50-mus_musculus-is_primary_data == True and tissue_general == "liver"] PASSED [ 51%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-None-batch_key3-50-mus_musculus-is_primary_data == True and tissue_general == "skin of body"] PASSED [ 51%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-None-batch_key3-50-mus_musculus-is_primary_data == True and tissue_general in ["heart", "lung"]] PASSED [ 51%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-None-batch_key3-50-mus_musculus-is_primary_data == True and assay == "Smart-seq2"] PASSED [ 52%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-None-batch_key3-500-mus_musculus-is_primary_data == True and tissue_general == "liver"] PASSED [ 52%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-None-batch_key3-500-mus_musculus-is_primary_data == True and tissue_general == "skin of body"] PASSED [ 52%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-None-batch_key3-500-mus_musculus-is_primary_data == True and tissue_general in ["heart", "lung"]] PASSED [ 52%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-None-batch_key3-500-mus_musculus-is_primary_data == True and assay == "Smart-seq2"] PASSED [ 53%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-0.5-None-50-mus_musculus-is_primary_data == True and tissue_general == "liver"] PASSED [ 53%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-0.5-None-50-mus_musculus-is_primary_data == True and tissue_general == "skin of body"] PASSED [ 53%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-0.5-None-50-mus_musculus-is_primary_data == True and tissue_general in ["heart", "lung"]] PASSED [ 53%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-0.5-None-50-mus_musculus-is_primary_data == True and assay == "Smart-seq2"] PASSED [ 54%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-0.5-None-500-mus_musculus-is_primary_data == True and tissue_general == "liver"] PASSED [ 54%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-0.5-None-500-mus_musculus-is_primary_data == True and tissue_general == "skin of body"] PASSED [ 54%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-0.5-None-500-mus_musculus-is_primary_data == True and tissue_general in ["heart", "lung"]] PASSED [ 54%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-0.5-None-500-mus_musculus-is_primary_data == True and assay == "Smart-seq2"] PASSED [ 54%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-0.5-dataset_id-50-mus_musculus-is_primary_data == True and tissue_general == "liver"] PASSED [ 55%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-0.5-dataset_id-50-mus_musculus-is_primary_data == True and tissue_general == "skin of body"] PASSED [ 55%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-0.5-dataset_id-50-mus_musculus-is_primary_data == True and tissue_general in ["heart", "lung"]] PASSED [ 55%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-0.5-dataset_id-50-mus_musculus-is_primary_data == True and assay == "Smart-seq2"] PASSED [ 55%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-0.5-dataset_id-500-mus_musculus-is_primary_data == True and tissue_general == "liver"] PASSED [ 56%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-0.5-dataset_id-500-mus_musculus-is_primary_data == True and tissue_general == "skin of body"] PASSED [ 56%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-0.5-dataset_id-500-mus_musculus-is_primary_data == True and tissue_general in ["heart", "lung"]] PASSED [ 56%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-0.5-dataset_id-500-mus_musculus-is_primary_data == True and assay == "Smart-seq2"] PASSED [ 56%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-0.5-batch_key2-50-mus_musculus-is_primary_data == True and tissue_general == "liver"] PASSED [ 57%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-0.5-batch_key2-50-mus_musculus-is_primary_data == True and tissue_general == "skin of body"] PASSED [ 57%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-0.5-batch_key2-50-mus_musculus-is_primary_data == True and tissue_general in ["heart", "lung"]] PASSED [ 57%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-0.5-batch_key2-50-mus_musculus-is_primary_data == True and assay == "Smart-seq2"] PASSED [ 57%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-0.5-batch_key2-500-mus_musculus-is_primary_data == True and tissue_general == "liver"] PASSED [ 58%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-0.5-batch_key2-500-mus_musculus-is_primary_data == True and tissue_general == "skin of body"] PASSED [ 58%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-0.5-batch_key2-500-mus_musculus-is_primary_data == True and tissue_general in ["heart", "lung"]] PASSED [ 58%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-0.5-batch_key2-500-mus_musculus-is_primary_data == True and assay == "Smart-seq2"] PASSED [ 58%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-0.5-batch_key3-50-mus_musculus-is_primary_data == True and tissue_general == "liver"] PASSED [ 59%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-0.5-batch_key3-50-mus_musculus-is_primary_data == True and tissue_general == "skin of body"] PASSED [ 59%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-0.5-batch_key3-50-mus_musculus-is_primary_data == True and tissue_general in ["heart", "lung"]] PASSED [ 59%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-0.5-batch_key3-50-mus_musculus-is_primary_data == True and assay == "Smart-seq2"] PASSED [ 59%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-0.5-batch_key3-500-mus_musculus-is_primary_data == True and tissue_general == "liver"] PASSED [ 60%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-0.5-batch_key3-500-mus_musculus-is_primary_data == True and tissue_general == "skin of body"] PASSED [ 60%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-0.5-batch_key3-500-mus_musculus-is_primary_data == True and tissue_general in ["heart", "lung"]] PASSED [ 60%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-0.5-batch_key3-500-mus_musculus-is_primary_data == True and assay == "Smart-seq2"] PASSED [ 60%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_get_highly_variable_genes[mus_musculus-Mus musculus-tissue_general == "liver" and is_primary_data == True-None-None] PASSED [ 61%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_get_highly_variable_genes[mus_musculus-Mus musculus-is_primary_data == True and tissue_general == "heart"-dataset_id-None] PASSED [ 61%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_get_highly_variable_genes[mus_musculus-Mus musculus-is_primary_data == True and tissue_general == "heart"-batch_key2-None] PASSED [ 61%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_get_highly_variable_genes[mus_musculus-Mus musculus-is_primary_data == True-dataset_id-obs_coords3] PASSED [ 61%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_get_highly_variable_genes[homo_sapiens-Homo sapiens-is_primary_data == True-dataset_id-obs_coords4] FAILED [ 62%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_error_cases PASSED [ 62%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_max_loess_jitter_error PASSED [ 62%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_user_defined_batch_key_func[None] PASSED [ 62%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_user_defined_batch_key_func[suspension_type] PASSED [ 63%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_user_defined_batch_key_func[batch_key2] PASSED [ 63%]
api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_hvg_user_defined_batch_key_func[batch_key3] PASSED [ 63%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_meanvar[True-0-1200-511-1-101] PASSED [ 63%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_meanvar[True-0-1200-511-1-53] PASSED [ 63%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_meanvar[True-0-1200-511-3-101] PASSED [ 64%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_meanvar[True-0-1200-511-3-53] PASSED [ 64%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_meanvar[True-0-1200-511-11-101] PASSED [ 64%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_meanvar[True-0-1200-511-11-53] PASSED [ 64%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_meanvar[True-0-1200-511-101-101] PASSED [ 65%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_meanvar[True-0-1200-511-101-53] PASSED [ 65%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_meanvar[True-0-100001-57-1-101] PASSED [ 65%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_meanvar[True-0-100001-57-1-53] PASSED [ 65%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_meanvar[True-0-100001-57-3-101] PASSED [ 66%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_meanvar[True-0-100001-57-3-53] PASSED [ 66%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_meanvar[True-0-100001-57-11-101] PASSED [ 66%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_meanvar[True-0-100001-57-11-53] PASSED [ 66%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_meanvar[True-0-100001-57-101-101] PASSED [ 67%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_meanvar[True-0-100001-57-101-53] PASSED [ 67%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_meanvar[True-1-1200-511-1-101] PASSED [ 67%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_meanvar[True-1-1200-511-1-53] PASSED [ 67%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_meanvar[True-1-1200-511-3-101] PASSED [ 68%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_meanvar[True-1-1200-511-3-53] PASSED [ 68%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_meanvar[True-1-1200-511-11-101] PASSED [ 68%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_meanvar[True-1-1200-511-11-53] PASSED [ 68%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_meanvar[True-1-1200-511-101-101] PASSED [ 69%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_meanvar[True-1-1200-511-101-53] PASSED [ 69%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_meanvar[True-1-100001-57-1-101] PASSED [ 69%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_meanvar[True-1-100001-57-1-53] PASSED [ 69%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_meanvar[True-1-100001-57-3-101] PASSED [ 70%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_meanvar[True-1-100001-57-3-53] PASSED [ 70%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_meanvar[True-1-100001-57-11-101] PASSED [ 70%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_meanvar[True-1-100001-57-11-53] PASSED [ 70%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_meanvar[True-1-100001-57-101-101] PASSED [ 71%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_meanvar[True-1-100001-57-101-53] PASSED [ 71%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_meanvar[True-100-1200-511-1-101] PASSED [ 71%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_meanvar[True-100-1200-511-1-53] PASSED [ 71%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_meanvar[True-100-1200-511-3-101] PASSED [ 72%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_meanvar[True-100-1200-511-3-53] PASSED [ 72%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_meanvar[True-100-1200-511-11-101] PASSED [ 72%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_meanvar[True-100-1200-511-11-53] PASSED [ 72%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_meanvar[True-100-1200-511-101-101] PASSED [ 72%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_meanvar[True-100-1200-511-101-53] PASSED [ 73%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_meanvar[True-100-100001-57-1-101] PASSED [ 73%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_meanvar[True-100-100001-57-1-53] PASSED [ 73%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_meanvar[True-100-100001-57-3-101] PASSED [ 73%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_meanvar[True-100-100001-57-3-53] PASSED [ 74%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_meanvar[True-100-100001-57-11-101] PASSED [ 74%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_meanvar[True-100-100001-57-11-53] PASSED [ 74%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_meanvar[True-100-100001-57-101-101] PASSED [ 74%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_meanvar[True-100-100001-57-101-53] PASSED [ 75%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_meanvar[False-0-1200-511-1-101] PASSED [ 75%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_meanvar[False-0-1200-511-1-53] PASSED [ 75%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_meanvar[False-0-1200-511-3-101] PASSED [ 75%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_meanvar[False-0-1200-511-3-53] PASSED [ 76%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_meanvar[False-0-1200-511-11-101] PASSED [ 76%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_meanvar[False-0-1200-511-11-53] PASSED [ 76%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_meanvar[False-0-1200-511-101-101] PASSED [ 76%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_meanvar[False-0-1200-511-101-53] PASSED [ 77%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_meanvar[False-0-100001-57-1-101] PASSED [ 77%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_meanvar[False-0-100001-57-1-53] PASSED [ 77%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_meanvar[False-0-100001-57-3-101] PASSED [ 77%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_meanvar[False-0-100001-57-3-53] PASSED [ 78%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_meanvar[False-0-100001-57-11-101] PASSED [ 78%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_meanvar[False-0-100001-57-11-53] PASSED [ 78%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_meanvar[False-0-100001-57-101-101] PASSED [ 78%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_meanvar[False-0-100001-57-101-53] PASSED [ 79%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_meanvar[False-1-1200-511-1-101] PASSED [ 79%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_meanvar[False-1-1200-511-1-53] PASSED [ 79%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_meanvar[False-1-1200-511-3-101] PASSED [ 79%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_meanvar[False-1-1200-511-3-53] PASSED [ 80%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_meanvar[False-1-1200-511-11-101] PASSED [ 80%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_meanvar[False-1-1200-511-11-53] PASSED [ 80%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_meanvar[False-1-1200-511-101-101] PASSED [ 80%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_meanvar[False-1-1200-511-101-53] PASSED [ 81%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_meanvar[False-1-100001-57-1-101] PASSED [ 81%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_meanvar[False-1-100001-57-1-53] PASSED [ 81%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_meanvar[False-1-100001-57-3-101] PASSED [ 81%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_meanvar[False-1-100001-57-3-53] PASSED [ 81%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_meanvar[False-1-100001-57-11-101] PASSED [ 82%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_meanvar[False-1-100001-57-11-53] PASSED [ 82%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_meanvar[False-1-100001-57-101-101] PASSED [ 82%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_meanvar[False-1-100001-57-101-53] PASSED [ 82%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_meanvar[False-100-1200-511-1-101] PASSED [ 83%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_meanvar[False-100-1200-511-1-53] PASSED [ 83%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_meanvar[False-100-1200-511-3-101] PASSED [ 83%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_meanvar[False-100-1200-511-3-53] PASSED [ 83%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_meanvar[False-100-1200-511-11-101] PASSED [ 84%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_meanvar[False-100-1200-511-11-53] PASSED [ 84%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_meanvar[False-100-1200-511-101-101] PASSED [ 84%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_meanvar[False-100-1200-511-101-53] PASSED [ 84%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_meanvar[False-100-100001-57-1-101] PASSED [ 85%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_meanvar[False-100-100001-57-1-53] PASSED [ 85%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_meanvar[False-100-100001-57-3-101] PASSED [ 85%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_meanvar[False-100-100001-57-3-53] PASSED [ 85%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_meanvar[False-100-100001-57-11-101] PASSED [ 86%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_meanvar[False-100-100001-57-11-53] PASSED [ 86%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_meanvar[False-100-100001-57-101-101] PASSED [ 86%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_meanvar[False-100-100001-57-101-53] PASSED [ 86%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_meanvar_nnz_only_batches_fails[1200-511] PASSED [ 87%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_meanvar_nnz_only_batches_fails[100001-57] PASSED [ 87%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_mean[1200-511-101] PASSED [ 87%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_mean[1200-511-53] PASSED [ 87%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_mean[100001-57-101] PASSED [ 88%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_mean[100001-57-53] PASSED [ 88%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_counts[1200-511-1-101] PASSED [ 88%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_counts[1200-511-1-53] PASSED [ 88%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_counts[1200-511-3-101] PASSED [ 89%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_counts[1200-511-3-53] PASSED [ 89%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_counts[1200-511-11-101] PASSED [ 89%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_counts[1200-511-11-53] PASSED [ 89%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_counts[1200-511-101-101] PASSED [ 90%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_counts[1200-511-101-53] PASSED [ 90%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_counts[100001-57-1-101] PASSED [ 90%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_counts[100001-57-1-53] PASSED [ 90%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_counts[100001-57-3-101] PASSED [ 90%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_counts[100001-57-3-53] PASSED [ 91%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_counts[100001-57-11-101] PASSED [ 91%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_counts[100001-57-11-53] PASSED [ 91%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_counts[100001-57-101-101] PASSED [ 91%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_counts[100001-57-101-53] PASSED [ 92%]
api/python/cellxgene_census/tests/experimental/pp/test_online.py::test_mean_fails_no_variables_or_samples PASSED [ 92%]
api/python/cellxgene_census/tests/experimental/pp/test_stats.py::test_mean_variance[mus_musculus-tissue_general == "liver" and is_primary_data == True-obs_coords0-True-True-0] PASSED [ 92%]
api/python/cellxgene_census/tests/experimental/pp/test_stats.py::test_mean_variance[mus_musculus-tissue_general == "liver" and is_primary_data == True-obs_coords0-True-True-1] PASSED [ 92%]
api/python/cellxgene_census/tests/experimental/pp/test_stats.py::test_mean_variance[mus_musculus-tissue_general == "liver" and is_primary_data == True-obs_coords0-True-False-0] PASSED [ 93%]
api/python/cellxgene_census/tests/experimental/pp/test_stats.py::test_mean_variance[mus_musculus-tissue_general == "liver" and is_primary_data == True-obs_coords0-True-False-1] PASSED [ 93%]
api/python/cellxgene_census/tests/experimental/pp/test_stats.py::test_mean_variance[mus_musculus-tissue_general == "liver" and is_primary_data == True-obs_coords0-False-True-0] PASSED [ 93%]
api/python/cellxgene_census/tests/experimental/pp/test_stats.py::test_mean_variance[mus_musculus-tissue_general == "liver" and is_primary_data == True-obs_coords0-False-True-1] PASSED [ 93%]
api/python/cellxgene_census/tests/experimental/pp/test_stats.py::test_mean_variance[mus_musculus-is_primary_data == True and tissue_general == "heart"-obs_coords1-True-True-0] PASSED [ 94%]
api/python/cellxgene_census/tests/experimental/pp/test_stats.py::test_mean_variance[mus_musculus-is_primary_data == True and tissue_general == "heart"-obs_coords1-True-True-1] PASSED [ 94%]
api/python/cellxgene_census/tests/experimental/pp/test_stats.py::test_mean_variance[mus_musculus-is_primary_data == True and tissue_general == "heart"-obs_coords1-True-False-0] PASSED [ 94%]
api/python/cellxgene_census/tests/experimental/pp/test_stats.py::test_mean_variance[mus_musculus-is_primary_data == True and tissue_general == "heart"-obs_coords1-True-False-1] PASSED [ 94%]
api/python/cellxgene_census/tests/experimental/pp/test_stats.py::test_mean_variance[mus_musculus-is_primary_data == True and tissue_general == "heart"-obs_coords1-False-True-0] PASSED [ 95%]
api/python/cellxgene_census/tests/experimental/pp/test_stats.py::test_mean_variance[mus_musculus-is_primary_data == True and tissue_general == "heart"-obs_coords1-False-True-1] PASSED [ 95%]
api/python/cellxgene_census/tests/experimental/pp/test_stats.py::test_mean_variance[mus_musculus-is_primary_data == True-obs_coords2-True-True-0] PASSED [ 95%]
api/python/cellxgene_census/tests/experimental/pp/test_stats.py::test_mean_variance[mus_musculus-is_primary_data == True-obs_coords2-True-True-1] PASSED [ 95%]
api/python/cellxgene_census/tests/experimental/pp/test_stats.py::test_mean_variance[mus_musculus-is_primary_data == True-obs_coords2-True-False-0] PASSED [ 96%]
api/python/cellxgene_census/tests/experimental/pp/test_stats.py::test_mean_variance[mus_musculus-is_primary_data == True-obs_coords2-True-False-1] PASSED [ 96%]
api/python/cellxgene_census/tests/experimental/pp/test_stats.py::test_mean_variance[mus_musculus-is_primary_data == True-obs_coords2-False-True-0] PASSED [ 96%]
api/python/cellxgene_census/tests/experimental/pp/test_stats.py::test_mean_variance[mus_musculus-is_primary_data == True-obs_coords2-False-True-1] PASSED [ 96%]
api/python/cellxgene_census/tests/experimental/pp/test_stats.py::test_mean_variance[homo_sapiens-is_primary_data == True-obs_coords3-True-True-0] FAILED [ 97%]
api/python/cellxgene_census/tests/experimental/pp/test_stats.py::test_mean_variance[homo_sapiens-is_primary_data == True-obs_coords3-True-True-1] FAILED [ 97%]
api/python/cellxgene_census/tests/experimental/pp/test_stats.py::test_mean_variance[homo_sapiens-is_primary_data == True-obs_coords3-True-False-0] FAILED [ 97%]
api/python/cellxgene_census/tests/experimental/pp/test_stats.py::test_mean_variance[homo_sapiens-is_primary_data == True-obs_coords3-True-False-1] FAILED [ 97%]
api/python/cellxgene_census/tests/experimental/pp/test_stats.py::test_mean_variance[homo_sapiens-is_primary_data == True-obs_coords3-False-True-0] FAILED [ 98%]
api/python/cellxgene_census/tests/experimental/pp/test_stats.py::test_mean_variance[homo_sapiens-is_primary_data == True-obs_coords3-False-True-1] FAILED [ 98%]
api/python/cellxgene_census/tests/experimental/pp/test_stats.py::test_mean_variance_nnz_only[mus_musculus-obs_coords0-True-True-0] PASSED [ 98%]
api/python/cellxgene_census/tests/experimental/pp/test_stats.py::test_mean_variance_nnz_only[mus_musculus-obs_coords0-True-True-1] PASSED [ 98%]
api/python/cellxgene_census/tests/experimental/pp/test_stats.py::test_mean_variance_no_flags PASSED [ 99%]
api/python/cellxgene_census/tests/experimental/pp/test_stats.py::test_mean_variance_empty_query[mus_musculus] PASSED [ 99%]
api/python/cellxgene_census/tests/experimental/pp/test_stats.py::test_mean_variance_wrong_axis PASSED [ 99%]
api/python/cellxgene_census/tests/experimental/util/test_csr_iter.py::test_X_sparse_iter PASSED [ 99%]
api/python/cellxgene_census/tests/experimental/util/test_csr_iter.py::test_X_sparse_iter_unsupported PASSED [100%]

============================== slowest durations ===============================
336.88s call     tests/test_acceptance.py::test_get_anndata[None-obs_coords10-ctx_config10-homo_sapiens]
316.66s call     tests/test_acceptance.py::test_get_anndata[tissue_general=='brain'-obs_coords7-ctx_config7-homo_sapiens]
276.94s call     tests/test_acceptance.py::test_get_anndata[None-obs_coords10-ctx_config10-mus_musculus]
235.17s call     tests/test_acceptance.py::test_get_anndata[tissue_general=='brain'-obs_coords7-ctx_config7-mus_musculus]
235.12s call     tests/test_acceptance.py::test_get_anndata[is_primary_data==True-obs_coords9-ctx_config9-mus_musculus]
162.02s call     tests/experimental/pp/test_hvg.py::test_get_highly_variable_genes[mus_musculus-Mus musculus-is_primary_data == True-dataset_id-obs_coords3]
151.34s call     tests/test_acceptance.py::test_get_anndata[First 1M cells-homo_sapiens]
120.43s call     tests/test_acceptance.py::test_get_anndata[First 750K cells-homo_sapiens]
114.92s call     tests/test_acceptance.py::test_get_anndata[First 1M cells-mus_musculus]
83.01s call     tests/test_acceptance.py::test_incremental_query[None-tissue=='brain'-homo_sapiens]
82.79s call     tests/test_acceptance.py::test_get_anndata[First 750K cells-mus_musculus]
82.26s call     tests/test_acceptance.py::test_incremental_read_X[None-ctx_config1-coords1-homo_sapiens]
76.06s call     tests/test_acceptance.py::test_get_anndata[First 500K cells-homo_sapiens]
73.13s call     tests/test_acceptance.py::test_incremental_read_X[None-ctx_config2-coords2-homo_sapiens]
71.89s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-None-None-50-mus_musculus-is_primary_data == True and tissue_general in ["heart", "lung"]]
63.83s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-None-None-50-mus_musculus-is_primary_data == True and tissue_general in ["heart", "lung"]]
52.57s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-None-None-50-mus_musculus-is_primary_data == True and assay == "Smart-seq2"]
51.95s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-None-None-500-mus_musculus-is_primary_data == True and tissue_general in ["heart", "lung"]]
51.87s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-0.5-None-50-mus_musculus-is_primary_data == True and assay == "Smart-seq2"]
51.51s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-0.5-batch_key3-50-mus_musculus-is_primary_data == True and assay == "Smart-seq2"]
51.41s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-None-dataset_id-50-mus_musculus-is_primary_data == True and assay == "Smart-seq2"]
51.37s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-None-None-500-mus_musculus-is_primary_data == True and assay == "Smart-seq2"]
51.14s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-None-dataset_id-500-mus_musculus-is_primary_data == True and assay == "Smart-seq2"]
50.82s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-None-batch_key2-500-mus_musculus-is_primary_data == True and assay == "Smart-seq2"]
50.71s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-None-dataset_id-50-mus_musculus-is_primary_data == True and assay == "Smart-seq2"]
50.48s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-0.5-batch_key2-50-mus_musculus-is_primary_data == True and assay == "Smart-seq2"]
50.41s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-None-batch_key2-500-mus_musculus-is_primary_data == True and assay == "Smart-seq2"]
50.28s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-None-None-500-mus_musculus-is_primary_data == True and tissue_general in ["heart", "lung"]]
50.21s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-None-batch_key3-500-mus_musculus-is_primary_data == True and assay == "Smart-seq2"]
50.14s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-0.5-batch_key2-500-mus_musculus-is_primary_data == True and assay == "Smart-seq2"]
50.05s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-None-None-50-mus_musculus-is_primary_data == True and assay == "Smart-seq2"]
49.70s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-None-dataset_id-500-mus_musculus-is_primary_data == True and assay == "Smart-seq2"]
49.68s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-0.5-dataset_id-500-mus_musculus-is_primary_data == True and assay == "Smart-seq2"]
49.62s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-0.5-dataset_id-50-mus_musculus-is_primary_data == True and assay == "Smart-seq2"]
49.59s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-0.5-batch_key3-500-mus_musculus-is_primary_data == True and assay == "Smart-seq2"]
49.58s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-None-None-500-mus_musculus-is_primary_data == True and assay == "Smart-seq2"]
49.51s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-None-batch_key3-500-mus_musculus-is_primary_data == True and assay == "Smart-seq2"]
49.25s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-None-batch_key2-50-mus_musculus-is_primary_data == True and assay == "Smart-seq2"]
49.20s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-0.5-dataset_id-50-mus_musculus-is_primary_data == True and assay == "Smart-seq2"]
49.20s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-0.5-None-500-mus_musculus-is_primary_data == True and assay == "Smart-seq2"]
49.19s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-None-batch_key3-50-mus_musculus-is_primary_data == True and assay == "Smart-seq2"]
48.98s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-0.5-dataset_id-500-mus_musculus-is_primary_data == True and assay == "Smart-seq2"]
48.91s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-0.5-None-50-mus_musculus-is_primary_data == True and assay == "Smart-seq2"]
48.80s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-0.5-batch_key2-500-mus_musculus-is_primary_data == True and assay == "Smart-seq2"]
48.70s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-0.5-batch_key3-500-mus_musculus-is_primary_data == True and tissue_general in ["heart", "lung"]]
48.64s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-0.5-batch_key3-50-mus_musculus-is_primary_data == True and tissue_general in ["heart", "lung"]]
48.60s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-None-batch_key3-50-mus_musculus-is_primary_data == True and assay == "Smart-seq2"]
48.59s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-0.5-batch_key2-500-mus_musculus-is_primary_data == True and tissue_general in ["heart", "lung"]]
48.43s call     tests/test_acceptance.py::test_incremental_read_X[None-ctx_config1-coords1-mus_musculus]
48.31s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-0.5-batch_key3-50-mus_musculus-is_primary_data == True and assay == "Smart-seq2"]
48.15s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-None-dataset_id-50-mus_musculus-is_primary_data == True and tissue_general in ["heart", "lung"]]
48.10s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-0.5-batch_key3-500-mus_musculus-is_primary_data == True and assay == "Smart-seq2"]
47.99s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-None-dataset_id-500-mus_musculus-is_primary_data == True and tissue_general in ["heart", "lung"]]
47.92s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-None-batch_key2-50-mus_musculus-is_primary_data == True and assay == "Smart-seq2"]
47.81s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-0.5-batch_key2-50-mus_musculus-is_primary_data == True and assay == "Smart-seq2"]
47.75s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-0.5-batch_key2-50-mus_musculus-is_primary_data == True and tissue_general in ["heart", "lung"]]
47.70s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-None-batch_key2-500-mus_musculus-is_primary_data == True and tissue_general in ["heart", "lung"]]
47.65s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-0.5-None-500-mus_musculus-is_primary_data == True and assay == "Smart-seq2"]
47.56s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-None-batch_key2-50-mus_musculus-is_primary_data == True and tissue_general in ["heart", "lung"]]
47.50s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-None-batch_key2-500-mus_musculus-is_primary_data == True and tissue_general in ["heart", "lung"]]
47.37s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-None-dataset_id-500-mus_musculus-is_primary_data == True and tissue_general in ["heart", "lung"]]
47.35s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-None-batch_key3-500-mus_musculus-is_primary_data == True and tissue_general in ["heart", "lung"]]
47.15s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-None-batch_key2-50-mus_musculus-is_primary_data == True and tissue_general in ["heart", "lung"]]
47.10s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-None-batch_key3-50-mus_musculus-is_primary_data == True and tissue_general in ["heart", "lung"]]
46.89s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-0.5-dataset_id-500-mus_musculus-is_primary_data == True and tissue_general in ["heart", "lung"]]
46.85s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-0.5-batch_key3-50-mus_musculus-is_primary_data == True and tissue_general in ["heart", "lung"]]
46.62s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-None-batch_key3-500-mus_musculus-is_primary_data == True and tissue_general in ["heart", "lung"]]
46.32s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-None-dataset_id-50-mus_musculus-is_primary_data == True and tissue_general in ["heart", "lung"]]
45.57s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-0.5-None-500-mus_musculus-is_primary_data == True and tissue_general in ["heart", "lung"]]
45.47s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-0.5-batch_key3-500-mus_musculus-is_primary_data == True and tissue_general in ["heart", "lung"]]
45.37s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-0.5-None-500-mus_musculus-is_primary_data == True and tissue_general in ["heart", "lung"]]
45.36s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-0.5-dataset_id-50-mus_musculus-is_primary_data == True and tissue_general in ["heart", "lung"]]
45.02s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-0.5-batch_key2-50-mus_musculus-is_primary_data == True and tissue_general in ["heart", "lung"]]
44.96s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-None-batch_key3-50-mus_musculus-is_primary_data == True and tissue_general in ["heart", "lung"]]
44.81s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-0.5-dataset_id-50-mus_musculus-is_primary_data == True and tissue_general in ["heart", "lung"]]
44.50s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-0.5-dataset_id-500-mus_musculus-is_primary_data == True and tissue_general in ["heart", "lung"]]
44.27s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-0.5-None-50-mus_musculus-is_primary_data == True and tissue_general in ["heart", "lung"]]
43.98s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-0.5-batch_key2-500-mus_musculus-is_primary_data == True and tissue_general in ["heart", "lung"]]
43.54s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-0.5-None-50-mus_musculus-is_primary_data == True and tissue_general in ["heart", "lung"]]
43.32s call     tests/test_acceptance.py::test_get_anndata[First 500K cells-mus_musculus]
41.88s call     tests/test_acceptance.py::test_get_anndata[First 250K cells-homo_sapiens]
40.84s call     tests/test_acceptance.py::test_incremental_query[2-tissue=='brain'-homo_sapiens]
40.68s call     tests/test_acceptance.py::test_incremental_read_X[None-ctx_config2-coords2-mus_musculus]
33.42s call     tests/experimental/pp/test_stats.py::test_mean_variance[mus_musculus-is_primary_data == True-obs_coords2-True-True-0]
31.24s call     tests/test_acceptance.py::test_incremental_query[2-tissue=='brain'-mus_musculus]
26.90s call     tests/test_acceptance.py::test_get_anndata[is_primary_data==True-obs_coords9-ctx_config9-homo_sapiens]
26.41s call     tests/test_acceptance.py::test_incremental_query[None-tissue=='brain'-mus_musculus]
23.43s call     tests/experimental/pp/test_hvg.py::test_get_highly_variable_genes[mus_musculus-Mus musculus-tissue_general == "liver" and is_primary_data == True-None-None]
23.29s call     tests/experimental/pp/test_hvg.py::test_hvg_user_defined_batch_key_func[None]
22.18s call     tests/test_acceptance.py::test_get_anndata[First 100K cells-homo_sapiens]
21.76s call     tests/experimental/util/test_csr_iter.py::test_X_sparse_iter
21.70s call     tests/experimental/pp/test_hvg.py::test_get_highly_variable_genes[mus_musculus-Mus musculus-is_primary_data == True and tissue_general == "heart"-dataset_id-None]
19.91s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-None-None-50-mus_musculus-is_primary_data == True and tissue_general == "liver"]
19.86s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-None-None-50-mus_musculus-is_primary_data == True and tissue_general == "skin of body"]
19.11s call     tests/experimental/pp/test_stats.py::test_mean_variance[mus_musculus-is_primary_data == True-obs_coords2-True-True-1]
18.78s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-None-batch_key2-500-mus_musculus-is_primary_data == True and tissue_general == "skin of body"]
18.71s call     tests/test_acceptance.py::test_get_anndata[First 250K cells-mus_musculus]
18.26s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-None-batch_key2-50-mus_musculus-is_primary_data == True and tissue_general == "liver"]
17.95s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-None-batch_key3-500-mus_musculus-is_primary_data == True and tissue_general == "skin of body"]
17.91s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-None-dataset_id-500-mus_musculus-is_primary_data == True and tissue_general == "skin of body"]
17.82s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-0.5-batch_key2-500-mus_musculus-is_primary_data == True and tissue_general == "liver"]
17.68s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-None-batch_key2-50-mus_musculus-is_primary_data == True and tissue_general == "skin of body"]
17.68s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-0.5-dataset_id-50-mus_musculus-is_primary_data == True and tissue_general == "skin of body"]
17.64s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-None-dataset_id-50-mus_musculus-is_primary_data == True and tissue_general == "skin of body"]
17.56s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-0.5-batch_key2-50-mus_musculus-is_primary_data == True and tissue_general == "skin of body"]
17.54s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-0.5-batch_key3-500-mus_musculus-is_primary_data == True and tissue_general == "skin of body"]
17.47s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-None-dataset_id-50-mus_musculus-is_primary_data == True and tissue_general == "skin of body"]
17.46s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-None-dataset_id-500-mus_musculus-is_primary_data == True and tissue_general == "liver"]
17.43s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-None-dataset_id-50-mus_musculus-is_primary_data == True and tissue_general == "liver"]
17.42s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-None-batch_key2-500-mus_musculus-is_primary_data == True and tissue_general == "skin of body"]
17.41s call     tests/experimental/pp/test_stats.py::test_mean_variance[mus_musculus-is_primary_data == True-obs_coords2-True-False-1]
17.41s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-0.5-batch_key3-500-mus_musculus-is_primary_data == True and tissue_general == "liver"]
17.40s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-None-None-500-mus_musculus-is_primary_data == True and tissue_general == "liver"]
17.39s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-None-batch_key3-50-mus_musculus-is_primary_data == True and tissue_general == "skin of body"]
17.36s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-None-None-50-mus_musculus-is_primary_data == True and tissue_general == "skin of body"]
17.31s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-0.5-dataset_id-500-mus_musculus-is_primary_data == True and tissue_general == "skin of body"]
17.28s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-None-batch_key3-50-mus_musculus-is_primary_data == True and tissue_general == "liver"]
17.17s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-None-batch_key3-500-mus_musculus-is_primary_data == True and tissue_general == "skin of body"]
17.17s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-0.5-batch_key2-50-mus_musculus-is_primary_data == True and tissue_general == "liver"]
17.16s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-None-None-500-mus_musculus-is_primary_data == True and tissue_general == "skin of body"]
17.13s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-None-dataset_id-500-mus_musculus-is_primary_data == True and tissue_general == "skin of body"]
17.11s call     tests/experimental/pp/test_stats.py::test_mean_variance[mus_musculus-is_primary_data == True-obs_coords2-False-True-1]
17.06s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-None-batch_key2-500-mus_musculus-is_primary_data == True and tissue_general == "liver"]
17.06s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-0.5-batch_key3-500-mus_musculus-is_primary_data == True and tissue_general == "skin of body"]
17.04s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-0.5-dataset_id-50-mus_musculus-is_primary_data == True and tissue_general == "skin of body"]
17.03s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-0.5-batch_key2-50-mus_musculus-is_primary_data == True and tissue_general == "skin of body"]
17.01s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-0.5-dataset_id-500-mus_musculus-is_primary_data == True and tissue_general == "liver"]
17.00s call     tests/experimental/pp/test_stats.py::test_mean_variance[mus_musculus-is_primary_data == True-obs_coords2-True-False-0]
17.00s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-None-batch_key2-50-mus_musculus-is_primary_data == True and tissue_general == "liver"]
16.98s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-0.5-batch_key2-500-mus_musculus-is_primary_data == True and tissue_general == "skin of body"]
16.93s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-None-dataset_id-50-mus_musculus-is_primary_data == True and tissue_general == "liver"]
16.91s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-0.5-dataset_id-500-mus_musculus-is_primary_data == True and tissue_general == "skin of body"]
16.89s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-None-batch_key2-50-mus_musculus-is_primary_data == True and tissue_general == "skin of body"]
16.89s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-0.5-batch_key3-50-mus_musculus-is_primary_data == True and tissue_general == "liver"]
16.85s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-0.5-batch_key2-50-mus_musculus-is_primary_data == True and tissue_general == "liver"]
16.70s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-0.5-batch_key2-500-mus_musculus-is_primary_data == True and tissue_general == "skin of body"]
16.69s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-0.5-dataset_id-50-mus_musculus-is_primary_data == True and tissue_general == "liver"]
16.67s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-None-None-50-mus_musculus-is_primary_data == True and tissue_general == "liver"]
16.62s call     tests/test_get_anndata.py::test_get_anndata_two_layers[layers0]
16.61s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-None-batch_key3-50-mus_musculus-is_primary_data == True and tissue_general == "liver"]
16.61s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-None-None-500-mus_musculus-is_primary_data == True and tissue_general == "skin of body"]
16.60s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-None-batch_key2-500-mus_musculus-is_primary_data == True and tissue_general == "liver"]
16.58s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-None-batch_key3-500-mus_musculus-is_primary_data == True and tissue_general == "liver"]
16.56s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-None-batch_key3-50-mus_musculus-is_primary_data == True and tissue_general == "skin of body"]
16.54s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-0.5-batch_key3-50-mus_musculus-is_primary_data == True and tissue_general == "skin of body"]
16.52s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-0.5-dataset_id-50-mus_musculus-is_primary_data == True and tissue_general == "liver"]
16.51s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-0.5-batch_key3-50-mus_musculus-is_primary_data == True and tissue_general == "liver"]
16.49s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-None-dataset_id-500-mus_musculus-is_primary_data == True and tissue_general == "liver"]
16.37s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-None-batch_key3-500-mus_musculus-is_primary_data == True and tissue_general == "liver"]
16.33s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-0.5-dataset_id-500-mus_musculus-is_primary_data == True and tissue_general == "liver"]
16.33s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-0.5-None-50-mus_musculus-is_primary_data == True and tissue_general == "liver"]
16.24s call     tests/experimental/pp/test_stats.py::test_mean_variance[mus_musculus-is_primary_data == True-obs_coords2-False-True-0]
16.15s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-0.5-batch_key3-50-mus_musculus-is_primary_data == True and tissue_general == "skin of body"]
16.06s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-0.5-None-500-mus_musculus-is_primary_data == True and tissue_general == "skin of body"]
15.99s call     tests/experimental/pp/test_hvg.py::test_get_highly_variable_genes[mus_musculus-Mus musculus-is_primary_data == True and tissue_general == "heart"-batch_key2-None]
15.97s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-0.5-batch_key3-500-mus_musculus-is_primary_data == True and tissue_general == "liver"]
15.94s call     tests/experimental/pp/test_stats.py::test_mean_variance[mus_musculus-is_primary_data == True and tissue_general == "heart"-obs_coords1-True-True-0]
15.76s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-0.5-batch_key2-500-mus_musculus-is_primary_data == True and tissue_general == "liver"]
15.60s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-0.5-None-500-mus_musculus-is_primary_data == True and tissue_general == "liver"]
15.52s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-None-None-500-mus_musculus-is_primary_data == True and tissue_general == "liver"]
15.51s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-0.5-None-50-mus_musculus-is_primary_data == True and tissue_general == "liver"]
15.49s call     tests/experimental/pp/test_stats.py::test_mean_variance[mus_musculus-tissue_general == "liver" and is_primary_data == True-obs_coords0-True-True-0]
15.40s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[stable-0.5-None-50-mus_musculus-is_primary_data == True and tissue_general == "skin of body"]
15.31s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-0.5-None-500-mus_musculus-is_primary_data == True and tissue_general == "skin of body"]
15.13s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-0.5-None-500-mus_musculus-is_primary_data == True and tissue_general == "liver"]
14.94s call     tests/experimental/pp/test_hvg.py::test_hvg_vs_scanpy[latest-0.5-None-50-mus_musculus-is_primary_data == True and tissue_general == "skin of body"]
14.55s call     tests/test_get_anndata.py::test_get_anndata_two_layers[layers1]
13.72s call     tests/experimental/pp/test_stats.py::test_mean_variance[mus_musculus-is_primary_data == True and tissue_general == "heart"-obs_coords1-True-False-0]
13.30s call     tests/test_get_anndata.py::test_get_anndata_allows_missing_obs_or_var_filter
13.27s call     tests/experimental/pp/test_stats.py::test_mean_variance[mus_musculus-is_primary_data == True and tissue_general == "heart"-obs_coords1-True-True-1]
13.22s call     tests/experimental/pp/test_stats.py::test_mean_variance[mus_musculus-is_primary_data == True and tissue_general == "heart"-obs_coords1-True-False-1]
13.20s call     tests/experimental/pp/test_stats.py::test_mean_variance[mus_musculus-is_primary_data == True and tissue_general == "heart"-obs_coords1-False-True-1]
12.98s call     tests/experimental/pp/test_stats.py::test_mean_variance[mus_musculus-is_primary_data == True and tissue_general == "heart"-obs_coords1-False-True-0]
12.61s call     tests/test_acceptance.py::test_incremental_read_X[None-ctx_config3-coords3-homo_sapiens]
12.55s call     tests/experimental/pp/test_hvg.py::test_max_loess_jitter_error
12.31s call     tests/experimental/pp/test_hvg.py::test_hvg_user_defined_batch_key_func[suspension_type]
11.96s call     tests/experimental/pp/test_hvg.py::test_hvg_user_defined_batch_key_func[batch_key2]
11.41s call     tests/test_acceptance.py::test_incremental_query[2-tissue=='aorta'-mus_musculus]
10.73s call     tests/experimental/pp/test_hvg.py::test_hvg_user_defined_batch_key_func[batch_key3]
10.66s call     tests/experimental/pp/test_stats.py::test_mean_variance[mus_musculus-tissue_general == "liver" and is_primary_data == True-obs_coords0-False-True-1]
10.59s call     tests/test_acceptance.py::test_get_anndata[First 100K cells-mus_musculus]
10.59s call     tests/experimental/pp/test_stats.py::test_mean_variance[mus_musculus-tissue_general == "liver" and is_primary_data == True-obs_coords0-True-True-1]
10.44s call     tests/experimental/pp/test_stats.py::test_mean_variance[mus_musculus-tissue_general == "liver" and is_primary_data == True-obs_coords0-True-False-0]
10.41s call     tests/test_acceptance.py::test_get_anndata[tissue=='aorta'-None-ctx_config0-homo_sapiens]
10.18s call     tests/test_acceptance.py::test_get_anndata[First 10K cells-homo_sapiens]
10.07s call     tests/experimental/pp/test_stats.py::test_mean_variance[mus_musculus-tissue_general == "liver" and is_primary_data == True-obs_coords0-True-False-1]
10.02s call     tests/experimental/pp/test_stats.py::test_mean_variance[mus_musculus-tissue_general == "liver" and is_primary_data == True-obs_coords0-False-True-0]
9.61s call     tests/test_acceptance.py::test_incremental_query[2-tissue=='aorta'-homo_sapiens]
9.19s call     tests/test_acceptance.py::test_incremental_query[None-tissue=='aorta'-mus_musculus]
8.98s call     tests/test_acceptance.py::test_get_anndata[cell_type=='neuron'-obs_coords8-ctx_config8-mus_musculus]
8.74s call     tests/test_acceptance.py::test_load_axes[homo_sapiens]
8.62s call     tests/test_lts_compat.py::test_read_arrays[stable]
8.23s call     tests/test_get_anndata.py::test_get_anndata_x_layer[raw]
8.18s call     tests/test_acceptance.py::test_get_anndata[tissue=='aorta'-None-ctx_config0-mus_musculus]
8.10s call     tests/test_acceptance.py::test_incremental_query[None-tissue=='aorta'-homo_sapiens]
7.92s call     tests/test_acceptance.py::test_incremental_read_X[None-ctx_config3-coords3-mus_musculus]
7.90s call     tests/test_acceptance.py::test_incremental_read_X[2-None-coords0-homo_sapiens]
7.88s call     tests/test_get_anndata.py::test_get_anndata_x_layer[normalized]
7.37s call     tests/test_lts_compat.py::test_read_dataframe[stable]
6.98s call     tests/test_acceptance.py::test_get_anndata[cell_type=='neuron'-obs_coords8-ctx_config8-homo_sapiens]
6.05s call     tests/experimental/pp/test_stats.py::test_mean_variance_nnz_only[mus_musculus-obs_coords0-True-True-0]
5.92s call     tests/test_directory.py::test_live_directory_contents
5.90s call     tests/test_get_anndata.py::test_get_anndata_value_filter
5.72s call     tests/experimental/pp/test_stats.py::test_mean_variance_nnz_only[mus_musculus-obs_coords0-True-True-1]
5.20s call     tests/test_acceptance.py::test_get_anndata[First 10K cells-mus_musculus]
4.72s call     tests/test_get_helpers.py::test_get_presence_matrix[homo_sapiens]
4.55s call     tests/test_open.py::test_get_source_h5ad_uri
4.54s call     tests/experimental/ml/test_pytorch.py::test_multiprocessing__returns_full_result[6-3-pytorch_x_value_gen]
4.19s call     tests/test_acceptance.py::test_incremental_read_X[2-None-coords0-mus_musculus]
3.46s call     tests/test_acceptance.py::test_load_axes[mus_musculus]
3.27s call     tests/experimental/pp/test_hvg.py::test_get_highly_variable_genes[homo_sapiens-Homo sapiens-is_primary_data == True-dataset_id-obs_coords4]
3.25s call     tests/test_acceptance.py::test_incremental_read_obs[None-ctx_config1-homo_sapiens]
2.88s call     tests/test_get_anndata.py::test_get_anndata_coords
2.72s call     tests/test_get_helpers.py::test_get_presence_matrix[mus_musculus]
2.50s call     tests/test_acceptance.py::test_incremental_read_obs[None-ctx_config1-mus_musculus]
2.44s call     tests/test_lts_compat.py::test_open[stable]
2.39s call     tests/experimental/ml/test_pytorch.py::test_batching__partial_soma_batches_are_concatenated[10-1-pytorch_x_value_gen-False]
2.04s call     tests/experimental/ml/test_pytorch.py::test_batching__partial_soma_batches_are_concatenated[10-1-pytorch_x_value_gen-True]
1.90s call     tests/test_acceptance.py::test_incremental_read_obs[2-None-homo_sapiens]
1.88s call     tests/experimental/pp/test_online.py::test_meanvar[True-0-100001-57-1-101]
1.85s call     tests/test_acceptance.py::test_incremental_read_obs[2-None-mus_musculus]
1.81s call     tests/experimental/pp/test_stats.py::test_mean_variance[homo_sapiens-is_primary_data == True-obs_coords3-True-True-0]
1.80s call     tests/test_directory.py::test_get_census_version_directory__include_retracted
1.80s call     tests/test_acceptance.py::test_incremental_read_var[2-None-homo_sapiens]
1.61s call     tests/test_open.py::test_download_source_h5ad
1.51s call     tests/test_acceptance.py::test_incremental_read_var[None-ctx_config1-homo_sapiens]
1.51s call     tests/test_acceptance.py::test_incremental_read_var[None-ctx_config1-mus_musculus]
1.42s call     tests/test_acceptance.py::test_incremental_read_var[2-None-mus_musculus]
1.40s call     tests/experimental/pp/test_stats.py::test_mean_variance[homo_sapiens-is_primary_data == True-obs_coords3-True-True-1]
1.40s call     tests/experimental/ml/test_pytorch.py::test_batching__partial_final_batch_size[5-3-pytorch_x_value_gen-False]
1.39s call     tests/experimental/pp/test_stats.py::test_mean_variance[homo_sapiens-is_primary_data == True-obs_coords3-True-False-0]
1.33s call     tests/experimental/pp/test_stats.py::test_mean_variance[homo_sapiens-is_primary_data == True-obs_coords3-False-True-1]
1.32s call     tests/experimental/pp/test_stats.py::test_mean_variance[homo_sapiens-is_primary_data == True-obs_coords3-True-False-1]
1.29s call     tests/experimental/pp/test_stats.py::test_mean_variance_empty_query[mus_musculus]
1.27s call     tests/experimental/ml/test_pytorch.py::test_batching__partial_final_batch_size[5-3-pytorch_x_value_gen-True]
1.21s call     tests/experimental/pp/test_stats.py::test_mean_variance[homo_sapiens-is_primary_data == True-obs_coords3-False-True-0]
1.16s setup    tests/test_open.py::test_download_source_h5ad
1.10s call     tests/experimental/ml/test_pytorch.py::test_experiment_dataloader__batched[6-3-pytorch_x_value_gen-False]
1.07s call     tests/experimental/ml/test_pytorch.py::test_batching__exactly_one_batch[3-3-pytorch_x_value_gen-False]
1.07s setup    tests/test_open.py::test_download_source_h5ad_errors
1.06s call     tests/experimental/ml/test_pytorch.py::test_experiment_dataloader__non_batched[3-3-pytorch_x_value_gen-False]
1.05s call     tests/experimental/ml/test_pytorch.py::test_batching__all_batches_full_size[6-3-pytorch_x_value_gen-False]
1.04s call     tests/experimental/ml/test_pytorch.py::test_unique_soma_joinids[obs_range1-3-pytorch_x_value_gen-False]
0.98s call     tests/test_open.py::test_get_source_h5ad_uri_errors
0.95s call     tests/test_open.py::test_open_soma_stable
0.95s call     tests/experimental/ml/test_pytorch.py::test_batching__exactly_one_batch[3-3-pytorch_x_value_gen-True]
0.95s call     tests/experimental/ml/test_pytorch.py::test_distributed__returns_data_partition_for_rank[6-3-pytorch_x_value_gen]
0.94s call     tests/experimental/ml/test_pytorch.py::test__pytorch_splitting[10-1-pytorch_x_value_gen]
0.93s call     tests/experimental/ml/test_pytorch.py::test_batching__all_batches_full_size[6-3-pytorch_x_value_gen-True]
0.93s call     tests/experimental/ml/test_pytorch.py::test_distributed_and_multiprocessing__returns_data_partition_for_rank[12-3-pytorch_x_value_gen]
0.93s call     tests/experimental/ml/test_pytorch.py::test_experiment_dataloader__batched[6-3-pytorch_x_value_gen-True]
0.93s call     tests/experimental/ml/test_pytorch.py::test_unique_soma_joinids[obs_range0-3-pytorch_x_value_gen-True]
0.92s call     tests/experimental/ml/test_pytorch.py::test_experiment_dataloader__non_batched[3-3-pytorch_x_value_gen-True]
0.90s call     tests/experimental/ml/test_pytorch.py::test__shuffle[16-1-pytorch_seq_x_value_gen]
0.87s call     tests/experimental/util/test_csr_iter.py::test_X_sparse_iter_unsupported
0.77s call     tests/experimental/pp/test_hvg.py::test_hvg_error_cases
0.72s call     tests/experimental/ml/test_pytorch.py::test_sparse_output__batched[6-3-pytorch_x_value_gen-False]
0.72s call     tests/experimental/ml/test_pytorch.py::test_non_batched[6-3-pytorch_x_value_gen-True]
0.71s call     tests/experimental/ml/test_pytorch.py::test_sparse_output__non_batched[6-3-pytorch_x_value_gen-False]
0.71s call     tests/experimental/ml/test_pytorch.py::test_non_batched[6-3-pytorch_x_value_gen-False]
0.68s call     tests/experimental/ml/test_pytorch.py::test__X_tensor_dtype_matches_X_matrix[6-3-pytorch_x_value_gen-False]
0.68s call     tests/experimental/pp/test_online.py::test_meanvar[True-0-100001-57-1-53]
0.67s call     tests/test_get_helpers.py::test_get_experiment
0.66s call     tests/experimental/pp/test_online.py::test_meanvar[True-100-100001-57-1-53]
0.62s call     tests/experimental/ml/test_pytorch.py::test_sparse_output__non_batched[6-3-pytorch_x_value_gen-True]
0.62s setup    tests/test_get_anndata.py::test_get_anndata_two_layers[layers0]
0.59s call     tests/experimental/ml/test_pytorch.py::test_encoders[3-3-pytorch_x_value_gen]
0.58s call     tests/experimental/ml/test_pytorch.py::test__X_tensor_dtype_matches_X_matrix[6-3-pytorch_x_value_gen-True]
0.58s call     tests/experimental/ml/test_pytorch.py::test_sparse_output__batched[6-3-pytorch_x_value_gen-True]
0.56s call     tests/experimental/ml/test_pytorch.py::test_batching__empty_query_result[6-3-pytorch_x_value_gen-False]
0.55s call     tests/experimental/ml/test_pytorch.py::test_batching__empty_query_result[6-3-pytorch_x_value_gen-True]
0.54s call     tests/experimental/pp/test_online.py::test_meanvar[True-1-100001-57-1-53]
0.54s call     tests/test_open.py::test_open_soma_with_context
0.54s call     tests/experimental/pp/test_online.py::test_meanvar[True-1-100001-57-1-101]
0.53s call     tests/experimental/pp/test_online.py::test_meanvar[True-100-100001-57-1-101]
0.53s call     tests/test_get_anndata.py::test_get_anndata_wrong_layer_names
0.50s setup    tests/test_get_anndata.py::test_get_anndata_x_layer[raw]
0.49s setup    tests/test_get_anndata.py::test_get_anndata_allows_missing_obs_or_var_filter
0.48s setup    tests/test_get_anndata.py::test_get_anndata_value_filter
0.48s setup    tests/test_get_anndata.py::test_get_anndata_two_layers[layers1]
0.48s setup    tests/test_get_anndata.py::test_get_anndata_x_layer[normalized]
0.46s call     tests/test_open.py::test_open_soma_latest
0.46s setup    tests/test_get_anndata.py::test_get_anndata_coords
0.45s setup    tests/test_get_anndata.py::test_get_anndata_wrong_layer_names
0.43s call     tests/test_open.py::test_can_open_with_anonymous_access
0.35s call     tests/test_open.py::test_opening_census_without_anon_access_fails_with_bogus_creds
0.24s setup    tests/experimental/ml/test_pytorch.py::test_non_batched[6-3-pytorch_x_value_gen-True]
0.22s call     tests/experimental/pp/test_online.py::test_meanvar[True-0-1200-511-1-101]
0.19s setup    tests/experimental/ml/test_pytorch.py::test_experiment_dataloader__non_batched[3-3-pytorch_x_value_gen-False]
0.18s setup    tests/experimental/ml/test_pytorch.py::test_experiment_dataloader__batched[6-3-pytorch_x_value_gen-False]
0.18s setup    tests/experimental/ml/test_pytorch.py::test__X_tensor_dtype_matches_X_matrix[6-3-pytorch_x_value_gen-True]
0.18s setup    tests/experimental/ml/test_pytorch.py::test_sparse_output__non_batched[6-3-pytorch_x_value_gen-True]
0.18s setup    tests/experimental/ml/test_pytorch.py::test_experiment_dataloader__batched[6-3-pytorch_x_value_gen-True]
0.18s setup    tests/experimental/ml/test_pytorch.py::test_experiment_dataloader__non_batched[3-3-pytorch_x_value_gen-True]
0.18s setup    tests/experimental/ml/test_pytorch.py::test_batching__exactly_one_batch[3-3-pytorch_x_value_gen-True]
0.18s setup    tests/experimental/ml/test_pytorch.py::test_batching__exactly_one_batch[3-3-pytorch_x_value_gen-False]
0.18s setup    tests/experimental/ml/test_pytorch.py::test__pytorch_splitting[10-1-pytorch_x_value_gen]
0.18s setup    tests/experimental/ml/test_pytorch.py::test_distributed_and_multiprocessing__returns_data_partition_for_rank[12-3-pytorch_x_value_gen]
0.18s call     tests/test_directory.py::test_get_census_version_description_errors
0.18s setup    tests/experimental/ml/test_pytorch.py::test_batching__empty_query_result[6-3-pytorch_x_value_gen-True]
0.18s setup    tests/experimental/ml/test_pytorch.py::test_encoders[3-3-pytorch_x_value_gen]
0.18s setup    tests/experimental/ml/test_pytorch.py::test_batching__partial_soma_batches_are_concatenated[10-1-pytorch_x_value_gen-False]
0.17s setup    tests/experimental/ml/test_pytorch.py::test_sparse_output__batched[6-3-pytorch_x_value_gen-True]
0.17s setup    tests/experimental/ml/test_pytorch.py::test_batching__empty_query_result[6-3-pytorch_x_value_gen-False]
0.17s setup    tests/experimental/ml/test_pytorch.py::test_batching__all_batches_full_size[6-3-pytorch_x_value_gen-True]
0.17s setup    tests/experimental/ml/test_pytorch.py::test__shuffle[16-1-pytorch_seq_x_value_gen]
0.17s setup    tests/experimental/ml/test_pytorch.py::test_batching__partial_final_batch_size[5-3-pytorch_x_value_gen-False]
0.17s setup    tests/experimental/ml/test_pytorch.py::test_sparse_output__batched[6-3-pytorch_x_value_gen-False]
0.17s setup    tests/experimental/ml/test_pytorch.py::test_non_batched[6-3-pytorch_x_value_gen-False]
0.17s setup    tests/experimental/ml/test_pytorch.py::test_unique_soma_joinids[obs_range0-3-pytorch_x_value_gen-True]
0.17s setup    tests/experimental/ml/test_pytorch.py::test_sparse_output__non_batched[6-3-pytorch_x_value_gen-False]
0.17s setup    tests/experimental/ml/test_pytorch.py::test_unique_soma_joinids[obs_range1-3-pytorch_x_value_gen-False]
0.17s setup    tests/experimental/ml/test_pytorch.py::test_batching__partial_final_batch_size[5-3-pytorch_x_value_gen-True]
0.17s setup    tests/experimental/ml/test_pytorch.py::test_distributed__returns_data_partition_for_rank[6-3-pytorch_x_value_gen]
0.17s setup    tests/experimental/ml/test_pytorch.py::test_batching__all_batches_full_size[6-3-pytorch_x_value_gen-False]
0.17s setup    tests/experimental/ml/test_pytorch.py::test_batching__partial_soma_batches_are_concatenated[10-1-pytorch_x_value_gen-True]
0.17s setup    tests/experimental/ml/test_pytorch.py::test__X_tensor_dtype_matches_X_matrix[6-3-pytorch_x_value_gen-False]
0.17s setup    tests/experimental/ml/test_pytorch.py::test_multiprocessing__returns_full_result[6-3-pytorch_x_value_gen]
0.16s call     tests/test_open.py::test_open_soma_errors
0.15s call     tests/test_open.py::test_open_soma_works_if_no_relative_uri_specified
0.15s call     tests/test_open.py::test_open_soma_defaults_to_stable
0.15s call     tests/test_open.py::test_open_soma_rejects_non_s3_mirror
0.14s call     tests/experimental/pp/test_online.py::test_meanvar[True-0-100001-57-101-53]
0.14s call     tests/experimental/pp/test_online.py::test_meanvar[True-1-100001-57-101-53]
0.14s call     tests/experimental/pp/test_online.py::test_meanvar[True-100-100001-57-101-53]
0.14s call     tests/experimental/pp/test_online.py::test_meanvar[False-100-100001-57-101-53]
0.14s call     tests/experimental/pp/test_online.py::test_meanvar[False-1-100001-57-101-53]
0.14s call     tests/experimental/pp/test_online.py::test_meanvar[False-0-100001-57-101-53]
0.13s call     tests/experimental/pp/test_online.py::test_meanvar[True-0-100001-57-101-101]
0.13s call     tests/experimental/pp/test_online.py::test_meanvar[True-1-100001-57-101-101]
0.13s call     tests/experimental/pp/test_online.py::test_meanvar[True-100-100001-57-101-101]
0.13s call     tests/experimental/pp/test_online.py::test_meanvar[False-100-100001-57-101-101]
0.13s call     tests/experimental/pp/test_online.py::test_meanvar[False-0-100001-57-101-101]
0.13s call     tests/experimental/pp/test_online.py::test_meanvar[False-1-100001-57-101-101]
0.11s call     tests/experimental/pp/test_online.py::test_meanvar[True-100-100001-57-11-53]
0.11s call     tests/experimental/pp/test_online.py::test_counts[100001-57-101-53]
0.11s call     tests/experimental/pp/test_online.py::test_meanvar[True-0-100001-57-11-53]
0.11s call     tests/experimental/pp/test_online.py::test_meanvar[True-1-100001-57-11-53]
0.11s call     tests/experimental/pp/test_online.py::test_meanvar[False-0-100001-57-3-53]
0.11s call     tests/experimental/pp/test_online.py::test_meanvar[False-100-100001-57-11-53]
0.11s call     tests/experimental/pp/test_online.py::test_meanvar[False-1-100001-57-11-53]
0.11s call     tests/experimental/pp/test_online.py::test_meanvar[False-0-100001-57-11-53]
0.11s call     tests/experimental/pp/test_online.py::test_meanvar[False-1-100001-57-3-53]
0.11s call     tests/experimental/pp/test_online.py::test_meanvar[True-0-100001-57-3-101]
0.11s call     tests/experimental/pp/test_online.py::test_meanvar[True-100-100001-57-3-53]
0.11s call     tests/experimental/pp/test_online.py::test_meanvar[True-0-100001-57-3-53]
0.11s call     tests/experimental/pp/test_online.py::test_meanvar[True-1-100001-57-3-53]
0.10s call     tests/experimental/pp/test_online.py::test_meanvar[True-0-100001-57-11-101]
0.10s call     tests/experimental/pp/test_online.py::test_meanvar[False-100-100001-57-3-53]
0.10s call     tests/experimental/pp/test_online.py::test_counts[100001-57-101-101]
0.10s call     tests/experimental/pp/test_online.py::test_meanvar[True-1-100001-57-11-101]
0.10s call     tests/experimental/pp/test_online.py::test_meanvar[False-100-100001-57-11-101]
0.10s call     tests/experimental/pp/test_online.py::test_meanvar[False-0-100001-57-11-101]
0.10s call     tests/experimental/pp/test_online.py::test_meanvar[True-100-100001-57-11-101]
0.10s call     tests/experimental/pp/test_online.py::test_meanvar[False-1-100001-57-11-101]
0.10s call     tests/experimental/pp/test_online.py::test_meanvar[True-1-100001-57-3-101]
0.09s call     tests/experimental/pp/test_online.py::test_meanvar[True-100-100001-57-3-101]
0.09s call     tests/experimental/pp/test_online.py::test_meanvar[False-0-100001-57-3-101]
0.09s call     tests/experimental/pp/test_online.py::test_meanvar[False-100-100001-57-3-101]
0.09s call     tests/experimental/pp/test_online.py::test_meanvar[False-1-100001-57-3-101]
0.09s call     tests/experimental/pp/test_online.py::test_meanvar[False-100-100001-57-1-53]
0.09s call     tests/experimental/pp/test_online.py::test_meanvar[False-1-100001-57-1-53]
0.09s call     tests/experimental/pp/test_online.py::test_meanvar[False-0-100001-57-1-53]
0.09s call     tests/experimental/pp/test_online.py::test_counts[100001-57-3-53]
0.09s call     tests/experimental/pp/test_online.py::test_meanvar[False-100-100001-57-1-101]
0.09s call     tests/experimental/pp/test_online.py::test_counts[100001-57-11-53]
0.08s call     tests/experimental/pp/test_online.py::test_meanvar[False-0-100001-57-1-101]
0.08s call     tests/experimental/pp/test_online.py::test_meanvar[False-1-100001-57-1-101]
0.08s call     tests/experimental/pp/test_online.py::test_counts[100001-57-11-101]
0.07s call     tests/experimental/pp/test_online.py::test_counts[100001-57-1-53]
0.07s call     tests/experimental/pp/test_online.py::test_counts[100001-57-3-101]
0.06s call     tests/experimental/pp/test_online.py::test_counts[100001-57-1-101]
0.03s setup    tests/experimental/pp/test_online.py::test_meanvar[True-0-100001-57-1-101]
0.03s call     tests/experimental/pp/test_online.py::test_meanvar[True-100-1200-511-1-101]
0.03s call     tests/experimental/pp/test_online.py::test_meanvar[True-100-1200-511-1-53]
0.03s call     tests/experimental/pp/test_online.py::test_meanvar[True-1-1200-511-1-101]
0.03s call     tests/experimental/pp/test_online.py::test_meanvar[True-1-1200-511-1-53]
0.03s call     tests/experimental/pp/test_online.py::test_meanvar[True-0-1200-511-1-53]
0.03s call     tests/experimental/pp/test_online.py::test_meanvar[True-100-1200-511-101-53]
0.03s call     tests/experimental/pp/test_online.py::test_meanvar[False-100-1200-511-101-53]
0.03s call     tests/experimental/pp/test_online.py::test_meanvar[True-100-1200-511-101-101]
0.03s call     tests/experimental/pp/test_online.py::test_meanvar[False-100-1200-511-101-101]
0.03s setup    tests/experimental/pp/test_online.py::test_counts[100001-57-1-101]
0.03s setup    tests/experimental/pp/test_online.py::test_meanvar[True-0-100001-57-1-53]
0.03s call     tests/experimental/pp/test_online.py::test_meanvar[True-1-1200-511-101-53]
0.03s call     tests/experimental/pp/test_online.py::test_meanvar[False-1-1200-511-101-53]
0.03s call     tests/experimental/pp/test_online.py::test_meanvar[False-0-1200-511-101-53]
0.03s setup    tests/experimental/pp/test_online.py::test_counts[100001-57-3-53]
0.03s setup    tests/experimental/pp/test_online.py::test_counts[100001-57-11-53]
0.03s call     tests/experimental/pp/test_online.py::test_meanvar[True-1-1200-511-101-101]
0.03s setup    tests/experimental/pp/test_online.py::test_counts[100001-57-101-53]
0.03s call     tests/experimental/pp/test_online.py::test_meanvar[False-1-1200-511-101-101]
0.03s call     tests/experimental/pp/test_online.py::test_meanvar[True-0-1200-511-101-53]
0.03s setup    tests/experimental/pp/test_online.py::test_meanvar[True-0-100001-57-3-53]
0.03s call     tests/experimental/pp/test_online.py::test_meanvar[False-0-1200-511-101-101]
0.03s setup    tests/experimental/pp/test_online.py::test_meanvar[False-100-100001-57-1-101]
0.03s setup    tests/experimental/pp/test_online.py::test_meanvar[True-1-100001-57-3-53]
0.03s setup    tests/experimental/pp/test_online.py::test_meanvar[False-100-100001-57-101-101]
0.03s setup    tests/experimental/pp/test_online.py::test_meanvar[False-1-100001-57-1-101]
0.03s setup    tests/experimental/pp/test_online.py::test_meanvar[False-100-100001-57-11-101]
0.03s setup    tests/experimental/pp/test_online.py::test_meanvar[True-100-100001-57-3-101]
0.03s setup    tests/experimental/pp/test_online.py::test_meanvar[False-0-100001-57-11-101]
0.03s setup    tests/experimental/pp/test_online.py::test_meanvar[False-1-100001-57-101-101]
0.03s setup    tests/experimental/pp/test_online.py::test_meanvar[False-1-100001-57-3-53]
0.03s setup    tests/experimental/pp/test_online.py::test_meanvar[True-1-100001-57-11-101]
0.03s call     tests/experimental/pp/test_online.py::test_meanvar[True-0-1200-511-101-101]
0.03s setup    tests/experimental/pp/test_online.py::test_meanvar[False-0-100001-57-1-101]
0.03s setup    tests/experimental/pp/test_online.py::test_meanvar[False-0-100001-57-101-101]
0.03s setup    tests/experimental/pp/test_online.py::test_meanvar[True-100-100001-57-3-53]
0.03s setup    tests/experimental/pp/test_online.py::test_meanvar[True-1-100001-57-3-101]
0.03s setup    tests/experimental/pp/test_online.py::test_meanvar[True-100-100001-57-11-101]
0.03s setup    tests/experimental/pp/test_online.py::test_meanvar[True-0-100001-57-11-101]
0.03s setup    tests/experimental/pp/test_online.py::test_meanvar[True-100-100001-57-11-53]
0.03s setup    tests/experimental/pp/test_online.py::test_meanvar[True-100-100001-57-1-101]
0.03s setup    tests/experimental/pp/test_online.py::test_meanvar[False-1-100001-57-11-101]
0.03s setup    tests/experimental/pp/test_online.py::test_meanvar[True-100-100001-57-101-101]
0.03s setup    tests/experimental/pp/test_online.py::test_meanvar[True-0-100001-57-3-101]
0.03s setup    tests/experimental/pp/test_online.py::test_meanvar[True-1-100001-57-1-101]
0.03s setup    tests/experimental/pp/test_online.py::test_meanvar[True-1-100001-57-1-53]
0.03s setup    tests/experimental/pp/test_online.py::test_counts[100001-57-1-53]
0.03s setup    tests/experimental/pp/test_online.py::test_meanvar[False-0-100001-57-3-53]
0.03s setup    tests/experimental/pp/test_online.py::test_meanvar[True-100-100001-57-1-53]
0.03s setup    tests/experimental/pp/test_online.py::test_meanvar[False-0-100001-57-11-53]
0.03s setup    tests/experimental/pp/test_online.py::test_counts[100001-57-101-101]
0.03s setup    tests/experimental/pp/test_online.py::test_meanvar[False-100-100001-57-11-53]
0.03s setup    tests/experimental/pp/test_online.py::test_meanvar[False-1-100001-57-11-53]
0.03s setup    tests/experimental/pp/test_online.py::test_counts[100001-57-11-101]
0.03s setup    tests/experimental/pp/test_online.py::test_meanvar[True-0-100001-57-11-53]
0.03s setup    tests/experimental/pp/test_online.py::test_meanvar[True-1-100001-57-11-53]
0.03s setup    tests/experimental/pp/test_online.py::test_meanvar[True-1-100001-57-101-101]
0.03s setup    tests/experimental/pp/test_online.py::test_meanvar[True-0-100001-57-101-101]
0.03s setup    tests/experimental/pp/test_online.py::test_mean[100001-57-53]
0.02s setup    tests/experimental/pp/test_online.py::test_meanvar[True-0-100001-57-101-53]
0.02s setup    tests/experimental/pp/test_online.py::test_meanvar[False-0-100001-57-101-53]
0.02s setup    tests/experimental/pp/test_online.py::test_meanvar[False-100-100001-57-101-53]
0.02s setup    tests/experimental/pp/test_online.py::test_meanvar[True-1-100001-57-101-53]
0.02s setup    tests/experimental/pp/test_online.py::test_meanvar[True-100-100001-57-101-53]
0.02s setup    tests/experimental/pp/test_online.py::test_meanvar[False-1-100001-57-101-53]
0.02s setup    tests/experimental/pp/test_online.py::test_meanvar_nnz_only_batches_fails[100001-57]
0.02s setup    tests/experimental/pp/test_online.py::test_meanvar[False-0-100001-57-3-101]
0.02s setup    tests/experimental/pp/test_online.py::test_mean[100001-57-101]
0.02s setup    tests/experimental/pp/test_online.py::test_counts[100001-57-3-101]
0.02s setup    tests/experimental/pp/test_online.py::test_meanvar[False-1-100001-57-3-101]
0.02s setup    tests/experimental/pp/test_online.py::test_meanvar[False-1-100001-57-1-53]
0.02s call     tests/experimental/pp/test_online.py::test_counts[1200-511-101-53]
0.02s setup    tests/experimental/pp/test_online.py::test_meanvar[False-0-100001-57-1-53]
0.02s setup    tests/experimental/pp/test_online.py::test_meanvar[False-100-100001-57-3-101]
0.02s setup    tests/experimental/pp/test_online.py::test_meanvar[False-100-100001-57-1-53]
0.02s call     tests/experimental/pp/test_online.py::test_counts[1200-511-101-101]
0.02s setup    tests/experimental/pp/test_online.py::test_meanvar[False-100-100001-57-3-53]
0.02s call     tests/experimental/pp/test_online.py::test_mean[100001-57-53]
0.01s call     tests/experimental/pp/test_online.py::test_mean[100001-57-101]
0.01s call     tests/experimental/pp/test_online.py::test_meanvar[False-0-1200-511-3-53]
0.01s call     tests/experimental/pp/test_online.py::test_meanvar[True-100-1200-511-3-53]
0.01s call     tests/experimental/pp/test_online.py::test_meanvar[False-100-1200-511-11-53]
0.01s call     tests/experimental/pp/test_online.py::test_meanvar[True-100-1200-511-11-53]
0.01s call     tests/experimental/pp/test_online.py::test_meanvar[True-0-1200-511-11-53]
0.01s call     tests/experimental/pp/test_online.py::test_meanvar[True-1-1200-511-11-53]
0.01s call     tests/experimental/pp/test_online.py::test_meanvar[False-0-1200-511-11-53]
0.01s call     tests/experimental/pp/test_online.py::test_meanvar[False-1-1200-511-11-53]
0.01s call     tests/experimental/pp/test_online.py::test_meanvar[False-100-1200-511-3-53]
0.01s call     tests/experimental/pp/test_online.py::test_meanvar[True-0-1200-511-3-101]
0.01s call     tests/experimental/pp/test_online.py::test_counts[1200-511-11-53]
0.01s call     tests/experimental/pp/test_online.py::test_meanvar[True-100-1200-511-11-101]
0.01s call     tests/experimental/pp/test_online.py::test_meanvar[False-100-1200-511-11-101]
0.01s call     tests/experimental/pp/test_online.py::test_meanvar[True-1-1200-511-11-101]
0.01s call     tests/experimental/pp/test_online.py::test_meanvar[True-0-1200-511-11-101]
0.01s call     tests/experimental/pp/test_online.py::test_meanvar[True-1-1200-511-3-53]
0.01s call     tests/experimental/pp/test_online.py::test_meanvar[False-0-1200-511-11-101]
0.01s call     tests/experimental/pp/test_online.py::test_meanvar[False-1-1200-511-3-53]
0.01s call     tests/experimental/pp/test_online.py::test_meanvar[False-1-1200-511-11-101]
0.01s call     tests/experimental/pp/test_online.py::test_meanvar[True-0-1200-511-3-53]
0.01s call     tests/experimental/pp/test_online.py::test_meanvar[True-100-1200-511-3-101]
0.01s call     tests/experimental/pp/test_online.py::test_meanvar[False-100-1200-511-3-101]
0.01s call     tests/experimental/pp/test_online.py::test_meanvar[True-1-1200-511-3-101]
0.01s call     tests/experimental/pp/test_online.py::test_counts[1200-511-11-101]
0.01s call     tests/experimental/pp/test_online.py::test_counts[1200-511-3-53]
0.01s call     tests/experimental/pp/test_online.py::test_meanvar[False-0-1200-511-3-101]
0.01s call     tests/experimental/pp/test_online.py::test_meanvar[False-1-1200-511-3-101]
0.01s call     tests/experimental/pp/test_online.py::test_counts[1200-511-3-101]
0.01s call     tests/experimental/pp/test_online.py::test_counts[1200-511-1-53]
0.01s call     tests/experimental/pp/test_online.py::test_meanvar[False-1-1200-511-1-53]
0.01s call     tests/experimental/pp/test_online.py::test_meanvar[False-0-1200-511-1-53]
0.01s call     tests/experimental/pp/test_online.py::test_meanvar[False-1-1200-511-1-101]
0.01s call     tests/experimental/pp/test_online.py::test_meanvar[False-0-1200-511-1-101]
0.01s call     tests/experimental/pp/test_online.py::test_meanvar[False-100-1200-511-1-53]
0.01s call     tests/experimental/pp/test_online.py::test_meanvar[False-100-1200-511-1-101]
0.01s call     tests/experimental/pp/test_online.py::test_counts[1200-511-1-101]
0.01s setup    tests/experimental/pp/test_online.py::test_meanvar[True-0-1200-511-1-101]

(739 durations < 0.005s hidden.  Use -vv to show these durations.)
=========================== short test summary info ============================
FAILED api/python/cellxgene_census/tests/test_open.py::test_open_soma_with_context
FAILED api/python/cellxgene_census/tests/experimental/pp/test_hvg.py::test_get_highly_variable_genes[homo_sapiens-Homo sapiens-is_primary_data == True-dataset_id-obs_coords4]
FAILED api/python/cellxgene_census/tests/experimental/pp/test_stats.py::test_mean_variance[homo_sapiens-is_primary_data == True-obs_coords3-True-True-0]
FAILED api/python/cellxgene_census/tests/experimental/pp/test_stats.py::test_mean_variance[homo_sapiens-is_primary_data == True-obs_coords3-True-True-1]
FAILED api/python/cellxgene_census/tests/experimental/pp/test_stats.py::test_mean_variance[homo_sapiens-is_primary_data == True-obs_coords3-True-False-0]
FAILED api/python/cellxgene_census/tests/experimental/pp/test_stats.py::test_mean_variance[homo_sapiens-is_primary_data == True-obs_coords3-True-False-1]
FAILED api/python/cellxgene_census/tests/experimental/pp/test_stats.py::test_mean_variance[homo_sapiens-is_primary_data == True-obs_coords3-False-True-0]
FAILED api/python/cellxgene_census/tests/experimental/pp/test_stats.py::test_mean_variance[homo_sapiens-is_primary_data == True-obs_coords3-False-True-1]
====== 8 failed, 398 passed, 5 skipped, 95 warnings in 7698.98s (2:08:18) ======
```

Note: failures have been identified to be related to the tests themselves, so we won't block the LTS release.

### 2023-07-26

- Host: EC2 instance type: `r6id.32xlarge`, all nvme mounted as swap.
- Uname: Linux 5.19.0-1028-aws #29~22.04.1-Ubuntu SMP Tue Jun 20 19:12:11 UTC 2023 x86_64 x86_64 x86_64 GNU/Linux
- Python & census versions:

```python
>>> import cellxgene_census, tiledbsoma
>>> tiledbsoma.show_package_versions()
tiledbsoma.__version__        1.2.7
TileDB-Py tiledb.version()    (0, 21, 3)
TileDB core version           2.15.2
libtiledbsoma version()       libtiledb=2.15.2
python version                3.10.6.final.0
OS version                    Linux 5.19.0-1028-aws
```

**Pytest output:**

```text
============================= test session starts ==============================
platform linux -- Python 3.10.6, pytest-7.1.3, pluggy-1.0.0 -- /home/ubuntu/venv/bin/python
cachedir: .pytest_cache
rootdir: /home/ubuntu/repos/cellxgene-census/api/python/cellxgene_census, configfile: pyproject.toml
plugins: anyio-3.6.2, requests-mock-1.11.0
collecting ... collected 274 items / 202 deselected / 72 selected

api/python/cellxgene_census/tests/test_acceptance.py::test_load_axes[homo_sapiens] PASSED [  1%]
api/python/cellxgene_census/tests/test_acceptance.py::test_load_axes[mus_musculus] PASSED [  2%]
api/python/cellxgene_census/tests/test_acceptance.py::test_incremental_read_obs[2-None-homo_sapiens] PASSED [  4%]
api/python/cellxgene_census/tests/test_acceptance.py::test_incremental_read_obs[2-None-mus_musculus] PASSED [  5%]
api/python/cellxgene_census/tests/test_acceptance.py::test_incremental_read_obs[None-ctx_config1-homo_sapiens] PASSED [  6%]
api/python/cellxgene_census/tests/test_acceptance.py::test_incremental_read_obs[None-ctx_config1-mus_musculus] PASSED [  8%]
api/python/cellxgene_census/tests/test_acceptance.py::test_incremental_read_var[2-None-homo_sapiens] PASSED [  9%]
api/python/cellxgene_census/tests/test_acceptance.py::test_incremental_read_var[2-None-mus_musculus] PASSED [ 11%]
api/python/cellxgene_census/tests/test_acceptance.py::test_incremental_read_var[None-ctx_config1-homo_sapiens] PASSED [ 12%]
api/python/cellxgene_census/tests/test_acceptance.py::test_incremental_read_var[None-ctx_config1-mus_musculus] PASSED [ 13%]
api/python/cellxgene_census/tests/test_acceptance.py::test_incremental_read_X[2-None-homo_sapiens] PASSED [ 15%]
api/python/cellxgene_census/tests/test_acceptance.py::test_incremental_read_X[2-None-mus_musculus] PASSED [ 16%]
api/python/cellxgene_census/tests/test_acceptance.py::test_incremental_read_X[None-ctx_config1-homo_sapiens] PASSED [ 18%]
api/python/cellxgene_census/tests/test_acceptance.py::test_incremental_read_X[None-ctx_config1-mus_musculus] PASSED [ 19%]
api/python/cellxgene_census/tests/test_acceptance.py::test_incremental_read_X[None-ctx_config2-homo_sapiens] PASSED [ 20%]
api/python/cellxgene_census/tests/test_acceptance.py::test_incremental_read_X[None-ctx_config2-mus_musculus] PASSED [ 22%]
api/python/cellxgene_census/tests/test_acceptance.py::test_incremental_query[2-tissue=='aorta'-homo_sapiens] PASSED [ 23%]
api/python/cellxgene_census/tests/test_acceptance.py::test_incremental_query[2-tissue=='aorta'-mus_musculus] PASSED [ 25%]
api/python/cellxgene_census/tests/test_acceptance.py::test_incremental_query[2-tissue=='brain'-homo_sapiens] PASSED [ 26%]
api/python/cellxgene_census/tests/test_acceptance.py::test_incremental_query[2-tissue=='brain'-mus_musculus] PASSED [ 27%]
api/python/cellxgene_census/tests/test_acceptance.py::test_incremental_query[None-tissue=='aorta'-homo_sapiens] PASSED [ 29%]
api/python/cellxgene_census/tests/test_acceptance.py::test_incremental_query[None-tissue=='aorta'-mus_musculus] PASSED [ 30%]
api/python/cellxgene_census/tests/test_acceptance.py::test_incremental_query[None-tissue=='brain'-homo_sapiens] PASSED [ 31%]
api/python/cellxgene_census/tests/test_acceptance.py::test_incremental_query[None-tissue=='brain'-mus_musculus] PASSED [ 33%]
api/python/cellxgene_census/tests/test_acceptance.py::test_get_anndata[tissue=='aorta'-None-ctx_config0-homo_sapiens] PASSED [ 34%]
api/python/cellxgene_census/tests/test_acceptance.py::test_get_anndata[tissue=='aorta'-None-ctx_config0-mus_musculus] PASSED [ 36%]
api/python/cellxgene_census/tests/test_acceptance.py::test_get_anndata[First 10K cells-homo_sapiens] PASSED [ 37%]
api/python/cellxgene_census/tests/test_acceptance.py::test_get_anndata[First 10K cells-mus_musculus] PASSED [ 38%]
api/python/cellxgene_census/tests/test_acceptance.py::test_get_anndata[First 100K cells-homo_sapiens] PASSED [ 40%]
api/python/cellxgene_census/tests/test_acceptance.py::test_get_anndata[First 100K cells-mus_musculus] PASSED [ 41%]
api/python/cellxgene_census/tests/test_acceptance.py::test_get_anndata[First 250K cells-homo_sapiens] PASSED [ 43%]
api/python/cellxgene_census/tests/test_acceptance.py::test_get_anndata[First 250K cells-mus_musculus] PASSED [ 44%]
api/python/cellxgene_census/tests/test_acceptance.py::test_get_anndata[First 500K cells-homo_sapiens] PASSED [ 45%]
api/python/cellxgene_census/tests/test_acceptance.py::test_get_anndata[First 500K cells-mus_musculus] PASSED [ 47%]
api/python/cellxgene_census/tests/test_acceptance.py::test_get_anndata[First 750K cells-homo_sapiens] PASSED [ 48%]
api/python/cellxgene_census/tests/test_acceptance.py::test_get_anndata[First 750K cells-mus_musculus] PASSED [ 50%]
api/python/cellxgene_census/tests/test_acceptance.py::test_get_anndata[First 1M cells-homo_sapiens] PASSED [ 51%]
api/python/cellxgene_census/tests/test_acceptance.py::test_get_anndata[First 1M cells-mus_musculus] PASSED [ 52%]
api/python/cellxgene_census/tests/test_acceptance.py::test_get_anndata[tissue=='brain'-None-ctx_config7-homo_sapiens] PASSED [ 54%]
api/python/cellxgene_census/tests/test_acceptance.py::test_get_anndata[tissue=='brain'-None-ctx_config7-mus_musculus] PASSED [ 55%]
api/python/cellxgene_census/tests/test_acceptance.py::test_get_anndata[cell_type=='neuron'-None-ctx_config8-homo_sapiens] PASSED [ 56%]
api/python/cellxgene_census/tests/test_acceptance.py::test_get_anndata[cell_type=='neuron'-None-ctx_config8-mus_musculus] PASSED [ 58%]
api/python/cellxgene_census/tests/test_acceptance.py::test_get_anndata[is_primary_data==True-None-ctx_config9-homo_sapiens] PASSED [ 59%]
api/python/cellxgene_census/tests/test_acceptance.py::test_get_anndata[is_primary_data==True-None-ctx_config9-mus_musculus] PASSED [ 61%]
api/python/cellxgene_census/tests/test_acceptance.py::test_get_anndata[None-None-ctx_config10-homo_sapiens] PASSED [ 62%]
api/python/cellxgene_census/tests/test_acceptance.py::test_get_anndata[None-None-ctx_config10-mus_musculus] PASSED [ 63%]
api/python/cellxgene_census/tests/test_directory.py::test_get_census_version_directory PASSED [ 65%]
api/python/cellxgene_census/tests/test_directory.py::test_get_census_version_description_errors PASSED [ 66%]
api/python/cellxgene_census/tests/test_directory.py::test_live_directory_contents PASSED [ 68%]
api/python/cellxgene_census/tests/test_get_anndata.py::test_get_anndata_value_filter PASSED [ 69%]
api/python/cellxgene_census/tests/test_get_anndata.py::test_get_anndata_coords PASSED [ 70%]
api/python/cellxgene_census/tests/test_get_anndata.py::test_get_anndata_allows_missing_obs_or_var_filter PASSED [ 72%]
api/python/cellxgene_census/tests/test_get_helpers.py::test_get_experiment PASSED [ 73%]
api/python/cellxgene_census/tests/test_get_helpers.py::test_get_presence_matrix[homo_sapiens] PASSED [ 75%]
api/python/cellxgene_census/tests/test_get_helpers.py::test_get_presence_matrix[mus_musculus] PASSED [ 76%]
api/python/cellxgene_census/tests/test_open.py::test_open_soma_stable PASSED [ 77%]
api/python/cellxgene_census/tests/test_open.py::test_open_soma_latest PASSED [ 79%]
api/python/cellxgene_census/tests/test_open.py::test_open_soma_with_context PASSED [ 80%]
api/python/cellxgene_census/tests/test_open.py::test_open_soma_invalid_args PASSED [ 81%]
api/python/cellxgene_census/tests/test_open.py::test_open_soma_errors PASSED [ 83%]
api/python/cellxgene_census/tests/test_open.py::test_open_soma_defaults_to_latest_if_missing_stable PASSED [ 84%]
api/python/cellxgene_census/tests/test_open.py::test_open_soma_defaults_to_stable PASSED [ 86%]
api/python/cellxgene_census/tests/test_open.py::test_get_source_h5ad_uri PASSED [ 87%]
api/python/cellxgene_census/tests/test_open.py::test_get_source_h5ad_uri_errors PASSED [ 88%]
api/python/cellxgene_census/tests/test_open.py::test_download_source_h5ad PASSED [ 90%]
api/python/cellxgene_census/tests/test_open.py::test_download_source_h5ad_errors PASSED [ 91%]
api/python/cellxgene_census/tests/test_open.py::test_opening_census_without_anon_access_fails_with_bogus_creds PASSED [ 93%]
api/python/cellxgene_census/tests/test_open.py::test_can_open_with_anonymous_access PASSED [ 94%]
api/python/cellxgene_census/tests/test_util.py::test_uri_join PASSED     [ 95%]
api/python/cellxgene_census/tests/experimental/pp/test_stats.py::test_mean_variance_no_flags PASSED [ 97%]
api/python/cellxgene_census/tests/experimental/pp/test_stats.py::test_mean_variance_empty_query[mus_musculus] PASSED [ 98%]
api/python/cellxgene_census/tests/experimental/pp/test_stats.py::test_mean_variance_wrong_axis PASSED [100%]

============================== slowest durations ===============================
8283.70s call     tests/test_acceptance.py::test_get_anndata[None-None-ctx_config10-homo_sapiens]
1767.98s call     tests/test_acceptance.py::test_get_anndata[is_primary_data==True-None-ctx_config9-homo_sapiens]
1304.89s call     tests/test_acceptance.py::test_incremental_read_X[None-ctx_config1-homo_sapiens]
953.35s call     tests/test_acceptance.py::test_incremental_read_X[None-ctx_config2-homo_sapiens]
903.42s call     tests/test_acceptance.py::test_get_anndata[cell_type=='neuron'-None-ctx_config8-homo_sapiens]
309.04s call     tests/test_acceptance.py::test_get_anndata[None-None-ctx_config10-mus_musculus]
195.72s call     tests/test_acceptance.py::test_get_anndata[is_primary_data==True-None-ctx_config9-mus_musculus]
135.65s call     tests/test_acceptance.py::test_incremental_read_X[None-ctx_config1-mus_musculus]
115.40s call     tests/test_acceptance.py::test_incremental_read_X[None-ctx_config2-mus_musculus]
65.46s call     tests/test_acceptance.py::test_get_anndata[First 1M cells-mus_musculus]
53.50s call     tests/test_acceptance.py::test_get_anndata[First 750K cells-mus_musculus]
44.52s call     tests/test_acceptance.py::test_get_anndata[cell_type=='neuron'-None-ctx_config8-mus_musculus]
38.62s call     tests/test_acceptance.py::test_get_anndata[tissue=='brain'-None-ctx_config7-homo_sapiens]
38.56s call     tests/test_acceptance.py::test_get_anndata[First 500K cells-mus_musculus]
34.90s call     tests/test_acceptance.py::test_get_anndata[First 1M cells-homo_sapiens]
31.58s call     tests/test_acceptance.py::test_incremental_query[None-tissue=='brain'-homo_sapiens]
26.10s call     tests/test_acceptance.py::test_get_anndata[First 750K cells-homo_sapiens]
23.89s call     tests/test_acceptance.py::test_get_anndata[First 500K cells-homo_sapiens]
23.75s call     tests/test_acceptance.py::test_get_anndata[First 250K cells-mus_musculus]
19.00s call     tests/test_acceptance.py::test_incremental_query[2-tissue=='brain'-homo_sapiens]
16.80s call     tests/test_acceptance.py::test_get_anndata[First 250K cells-homo_sapiens]
12.47s call     tests/test_acceptance.py::test_get_anndata[tissue=='brain'-None-ctx_config7-mus_musculus]
12.25s call     tests/test_acceptance.py::test_get_anndata[First 100K cells-mus_musculus]
10.88s call     tests/test_acceptance.py::test_get_anndata[First 100K cells-homo_sapiens]
10.52s call     tests/test_get_anndata.py::test_get_anndata_allows_missing_obs_or_var_filter
9.96s call     tests/test_acceptance.py::test_incremental_query[2-tissue=='aorta'-homo_sapiens]
9.67s call     tests/test_acceptance.py::test_incremental_read_X[2-None-homo_sapiens]
9.45s call     tests/test_acceptance.py::test_incremental_query[2-tissue=='brain'-mus_musculus]
9.11s call     tests/test_acceptance.py::test_load_axes[homo_sapiens]
8.81s call     tests/test_acceptance.py::test_incremental_query[None-tissue=='aorta'-homo_sapiens]
7.96s call     tests/test_acceptance.py::test_incremental_query[None-tissue=='brain'-mus_musculus]
7.77s call     tests/test_acceptance.py::test_get_anndata[tissue=='aorta'-None-ctx_config0-homo_sapiens]
7.14s call     tests/test_directory.py::test_live_directory_contents
5.47s call     tests/test_get_anndata.py::test_get_anndata_value_filter
5.30s call     tests/test_acceptance.py::test_get_anndata[First 10K cells-homo_sapiens]
4.82s call     tests/test_acceptance.py::test_get_anndata[First 10K cells-mus_musculus]
4.69s call     tests/test_open.py::test_get_source_h5ad_uri
4.47s call     tests/test_get_helpers.py::test_get_presence_matrix[homo_sapiens]
4.46s call     tests/test_acceptance.py::test_incremental_query[2-tissue=='aorta'-mus_musculus]
4.26s call     tests/test_acceptance.py::test_incremental_query[None-tissue=='aorta'-mus_musculus]
3.89s call     tests/test_acceptance.py::test_incremental_read_X[2-None-mus_musculus]
3.80s call     tests/test_acceptance.py::test_get_anndata[tissue=='aorta'-None-ctx_config0-mus_musculus]
3.39s call     tests/test_acceptance.py::test_incremental_read_obs[None-ctx_config1-homo_sapiens]
3.24s call     tests/test_get_helpers.py::test_get_presence_matrix[mus_musculus]
3.02s call     tests/test_get_anndata.py::test_get_anndata_coords
2.94s call     tests/test_open.py::test_download_source_h5ad
2.54s call     tests/test_acceptance.py::test_load_axes[mus_musculus]
1.80s call     tests/test_acceptance.py::test_incremental_read_obs[2-None-homo_sapiens]
1.62s call     tests/test_acceptance.py::test_incremental_read_obs[2-None-mus_musculus]
1.54s call     tests/test_acceptance.py::test_incremental_read_obs[None-ctx_config1-mus_musculus]
1.40s call     tests/experimental/pp/test_stats.py::test_mean_variance_empty_query[mus_musculus]
1.37s call     tests/test_acceptance.py::test_incremental_read_var[2-None-mus_musculus]
1.34s call     tests/test_acceptance.py::test_incremental_read_var[None-ctx_config1-homo_sapiens]
1.34s call     tests/test_acceptance.py::test_incremental_read_var[2-None-homo_sapiens]
1.33s setup    tests/test_open.py::test_download_source_h5ad_errors
1.28s call     tests/test_acceptance.py::test_incremental_read_var[None-ctx_config1-mus_musculus]
1.20s setup    tests/test_open.py::test_download_source_h5ad
1.13s call     tests/test_open.py::test_get_source_h5ad_uri_errors
0.76s call     tests/test_open.py::test_open_soma_with_context
0.71s call     tests/test_open.py::test_open_soma_stable
0.56s call     tests/test_directory.py::test_get_census_version_description_errors
0.48s call     tests/test_get_helpers.py::test_get_experiment
0.47s setup    tests/test_get_anndata.py::test_get_anndata_allows_missing_obs_or_var_filter
0.43s setup    tests/test_get_anndata.py::test_get_anndata_coords
0.39s call     tests/test_open.py::test_open_soma_defaults_to_latest_if_missing_stable
0.34s setup    tests/test_get_anndata.py::test_get_anndata_value_filter
0.34s call     tests/test_open.py::test_open_soma_latest
0.33s call     tests/test_open.py::test_can_open_with_anonymous_access
0.26s call     tests/test_open.py::test_opening_census_without_anon_access_fails_with_bogus_creds
0.04s setup    tests/test_directory.py::test_get_census_version_directory
0.03s teardown tests/test_acceptance.py::test_get_anndata[None-None-ctx_config10-homo_sapiens]
0.03s setup    tests/test_acceptance.py::test_get_anndata[None-None-ctx_config10-mus_musculus]
0.02s teardown tests/test_acceptance.py::test_get_anndata[is_primary_data==True-None-ctx_config9-homo_sapiens]
0.02s call     tests/test_directory.py::test_get_census_version_directory
0.02s setup    tests/test_acceptance.py::test_get_anndata[is_primary_data==True-None-ctx_config9-mus_musculus]
0.01s setup    tests/experimental/pp/test_stats.py::test_mean_variance_empty_query[mus_musculus]
0.01s teardown tests/test_acceptance.py::test_get_anndata[cell_type=='neuron'-None-ctx_config8-homo_sapiens]
0.01s setup    tests/test_acceptance.py::test_get_anndata[cell_type=='neuron'-None-ctx_config8-mus_musculus]
0.01s call     tests/experimental/pp/test_stats.py::test_mean_variance_no_flags

(137 durations < 0.005s hidden.  Use -vv to show these durations.)
=============== 72 passed, 202 deselected in 14584.03s (4:03:04) ===============
```

### 2023-06-23

- Host: EC2 instance type: `r6id.32xlarge`, all nvme mounted as swap.
- Uname: Linux 5.19.0-1025-aws #26~22.04.1-Ubuntu SMP Mon Apr 24 01:58:15 UTC 2023 x86_64 x86_64 x86_64 GNU/Linux
- Python & census versions:

```python
>>> import cellxgene_census, tiledbsoma
>>> tiledbsoma.show_package_versions()
tiledbsoma.__version__        1.2.5
TileDB-Py tiledb.version()    (0, 21, 5)
TileDB core version           2.15.4
libtiledbsoma version()       libtiledb=2.15.2
python version                3.10.6.final.0
OS version
>>> cellxgene_census.__version__
'1.2.1'
>>> cellxgene_census.get_census_version_description('latest')
{'release_date': None, 'release_build': '2023-06-20', 'soma': {'uri': 's3://cellxgene-data-public/cell-census/2023-06-20/soma/', 's3_region': 'us-west-2'}, 'h5ads': {'uri': 's3://cellxgene-data-public/cell-census/2023-06-20/h5ads/', 's3_region': 'us-west-2'}, 'alias': 'latest'}
```

**Pytest output:**

```text
============================= test session starts ==============================
platform linux -- Python 3.10.6, pytest-7.4.0, pluggy-1.2.0 -- /home/ubuntu/venv/bin/python3
cachedir: .pytest_cache
rootdir: /home/ubuntu/repos/cellxgene-census/api/python/cellxgene_census
configfile: pyproject.toml
plugins: requests-mock-1.11.0
collecting ... collected 180 items / 111 deselected / 69 selected

test_acceptance.py::test_load_axes[homo_sapiens] PASSED                  [  1%]
test_acceptance.py::test_load_axes[mus_musculus] PASSED                  [  2%]
test_acceptance.py::test_incremental_read_obs[2-None-homo_sapiens] PASSED [  4%]
test_acceptance.py::test_incremental_read_obs[2-None-mus_musculus] PASSED [  5%]
test_acceptance.py::test_incremental_read_obs[None-ctx_config1-homo_sapiens] PASSED [  7%]
test_acceptance.py::test_incremental_read_obs[None-ctx_config1-mus_musculus] PASSED [  8%]
test_acceptance.py::test_incremental_read_var[2-None-homo_sapiens] PASSED [ 10%]
test_acceptance.py::test_incremental_read_var[2-None-mus_musculus] PASSED [ 11%]
test_acceptance.py::test_incremental_read_var[None-ctx_config1-homo_sapiens] PASSED [ 13%]
test_acceptance.py::test_incremental_read_var[None-ctx_config1-mus_musculus] PASSED [ 14%]
test_acceptance.py::test_incremental_read_X[2-None-homo_sapiens] PASSED  [ 15%]
test_acceptance.py::test_incremental_read_X[2-None-mus_musculus] PASSED  [ 17%]
test_acceptance.py::test_incremental_read_X[None-ctx_config1-homo_sapiens] PASSED [ 18%]
test_acceptance.py::test_incremental_read_X[None-ctx_config1-mus_musculus] PASSED [ 20%]
test_acceptance.py::test_incremental_read_X[None-ctx_config2-homo_sapiens] PASSED [ 21%]
test_acceptance.py::test_incremental_read_X[None-ctx_config2-mus_musculus] PASSED [ 23%]
test_acceptance.py::test_incremental_query[2-tissue=='aorta'-homo_sapiens] PASSED [ 24%]
test_acceptance.py::test_incremental_query[2-tissue=='aorta'-mus_musculus] PASSED [ 26%]
test_acceptance.py::test_incremental_query[2-tissue=='brain'-homo_sapiens] PASSED [ 27%]
test_acceptance.py::test_incremental_query[2-tissue=='brain'-mus_musculus] PASSED [ 28%]
test_acceptance.py::test_incremental_query[None-tissue=='aorta'-homo_sapiens] PASSED [ 30%]
test_acceptance.py::test_incremental_query[None-tissue=='aorta'-mus_musculus] PASSED [ 31%]
test_acceptance.py::test_incremental_query[None-tissue=='brain'-homo_sapiens] PASSED [ 33%]
test_acceptance.py::test_incremental_query[None-tissue=='brain'-mus_musculus] PASSED [ 34%]
test_acceptance.py::test_get_anndata[tissue=='aorta'-None-ctx_config0-homo_sapiens] PASSED [ 36%]
test_acceptance.py::test_get_anndata[tissue=='aorta'-None-ctx_config0-mus_musculus] PASSED [ 37%]
test_acceptance.py::test_get_anndata[First 10K cells-homo_sapiens] PASSED [ 39%]
test_acceptance.py::test_get_anndata[First 10K cells-mus_musculus] PASSED [ 40%]
test_acceptance.py::test_get_anndata[First 100K cells-homo_sapiens] PASSED [ 42%]
test_acceptance.py::test_get_anndata[First 100K cells-mus_musculus] PASSED [ 43%]
test_acceptance.py::test_get_anndata[First 250K cells-homo_sapiens] PASSED [ 44%]
test_acceptance.py::test_get_anndata[First 250K cells-mus_musculus] PASSED [ 46%]
test_acceptance.py::test_get_anndata[First 500K cells-homo_sapiens] PASSED [ 47%]
test_acceptance.py::test_get_anndata[First 500K cells-mus_musculus] PASSED [ 49%]
test_acceptance.py::test_get_anndata[First 750K cells-homo_sapiens] PASSED [ 50%]
test_acceptance.py::test_get_anndata[First 750K cells-mus_musculus] PASSED [ 52%]
test_acceptance.py::test_get_anndata[First 1M cells-homo_sapiens] PASSED [ 53%]
test_acceptance.py::test_get_anndata[First 1M cells-mus_musculus] PASSED [ 55%]
test_acceptance.py::test_get_anndata[tissue=='brain'-None-ctx_config7-homo_sapiens] PASSED [ 56%]
test_acceptance.py::test_get_anndata[tissue=='brain'-None-ctx_config7-mus_musculus] PASSED [ 57%]
test_acceptance.py::test_get_anndata[cell_type=='neuron'-None-ctx_config8-homo_sapiens] PASSED [ 59%]
test_acceptance.py::test_get_anndata[cell_type=='neuron'-None-ctx_config8-mus_musculus] PASSED [ 60%]
test_acceptance.py::test_get_anndata[is_primary_data==True-None-ctx_config9-homo_sapiens] PASSED [ 62%]
test_acceptance.py::test_get_anndata[is_primary_data==True-None-ctx_config9-mus_musculus] PASSED [ 63%]
test_acceptance.py::test_get_anndata[None-None-ctx_config10-homo_sapiens] PASSED [ 65%]
test_acceptance.py::test_get_anndata[None-None-ctx_config10-mus_musculus] PASSED [ 66%]
test_directory.py::test_get_census_version_directory PASSED              [ 68%]
test_directory.py::test_get_census_version_description_errors PASSED     [ 69%]
test_directory.py::test_live_directory_contents PASSED                   [ 71%]
test_get_anndata.py::test_get_anndata_value_filter PASSED                [ 72%]
test_get_anndata.py::test_get_anndata_coords PASSED                      [ 73%]
test_get_anndata.py::test_get_anndata_allows_missing_obs_or_var_filter PASSED [ 75%]
test_get_helpers.py::test_get_experiment PASSED                          [ 76%]
test_get_helpers.py::test_get_presence_matrix[homo_sapiens] PASSED       [ 78%]
test_get_helpers.py::test_get_presence_matrix[mus_musculus] PASSED       [ 79%]
test_open.py::test_open_soma_stable PASSED                               [ 81%]
test_open.py::test_open_soma_latest PASSED                               [ 82%]
test_open.py::test_open_soma_with_context PASSED                         [ 84%]
test_open.py::test_open_soma_invalid_args PASSED                         [ 85%]
test_open.py::test_open_soma_errors PASSED                               [ 86%]
test_open.py::test_open_soma_defaults_to_latest_if_missing_stable PASSED [ 88%]
test_open.py::test_open_soma_defaults_to_stable PASSED                   [ 89%]
test_open.py::test_get_source_h5ad_uri PASSED                            [ 91%]
test_open.py::test_get_source_h5ad_uri_errors PASSED                     [ 92%]
test_open.py::test_download_source_h5ad PASSED                           [ 94%]
test_open.py::test_download_source_h5ad_errors PASSED                    [ 95%]
test_open.py::test_opening_census_without_anon_access_fails_with_bogus_creds PASSED [ 97%]
test_open.py::test_can_open_with_anonymous_access PASSED                 [ 98%]
test_util.py::test_uri_join PASSED                                       [100%]

============================== slowest durations ===============================
8331.89s call     tests/test_acceptance.py::test_get_anndata[None-None-ctx_config10-homo_sapiens]
2755.15s call     tests/test_acceptance.py::test_get_anndata[is_primary_data==True-None-ctx_config9-homo_sapiens]
1146.66s call     tests/test_acceptance.py::test_incremental_read_X[None-ctx_config1-homo_sapiens]
891.72s call     tests/test_acceptance.py::test_incremental_read_X[None-ctx_config2-homo_sapiens]
650.78s call     tests/test_acceptance.py::test_get_anndata[cell_type=='neuron'-None-ctx_config8-homo_sapiens]
287.95s call     tests/test_acceptance.py::test_get_anndata[None-None-ctx_config10-mus_musculus]
190.27s call     tests/test_acceptance.py::test_get_anndata[is_primary_data==True-None-ctx_config9-mus_musculus]
118.72s call     tests/test_acceptance.py::test_incremental_read_X[None-ctx_config1-mus_musculus]
102.51s call     tests/test_acceptance.py::test_incremental_read_X[None-ctx_config2-mus_musculus]
55.81s call     tests/test_acceptance.py::test_get_anndata[First 1M cells-mus_musculus]
48.63s call     tests/test_acceptance.py::test_get_anndata[First 750K cells-mus_musculus]
43.77s call     tests/test_acceptance.py::test_get_anndata[cell_type=='neuron'-None-ctx_config8-mus_musculus]
36.64s call     tests/test_acceptance.py::test_get_anndata[tissue=='brain'-None-ctx_config7-homo_sapiens]
35.32s call     tests/test_acceptance.py::test_get_anndata[First 1M cells-homo_sapiens]
34.39s call     tests/test_acceptance.py::test_get_anndata[First 500K cells-mus_musculus]
29.83s call     tests/test_acceptance.py::test_get_anndata[First 750K cells-homo_sapiens]
27.46s call     tests/test_acceptance.py::test_incremental_query[None-tissue=='brain'-homo_sapiens]
23.21s call     tests/test_acceptance.py::test_get_anndata[First 500K cells-homo_sapiens]
20.32s call     tests/test_acceptance.py::test_get_anndata[First 250K cells-mus_musculus]
18.22s call     tests/test_acceptance.py::test_incremental_query[2-tissue=='brain'-homo_sapiens]
16.58s call     tests/test_acceptance.py::test_get_anndata[First 250K cells-homo_sapiens]
11.48s call     tests/test_acceptance.py::test_get_anndata[tissue=='brain'-None-ctx_config7-mus_musculus]
10.18s call     tests/test_acceptance.py::test_incremental_query[2-tissue=='aorta'-homo_sapiens]
9.88s call     tests/test_acceptance.py::test_get_anndata[First 100K cells-homo_sapiens]
9.47s call     tests/test_acceptance.py::test_incremental_query[2-tissue=='brain'-mus_musculus]
9.44s call     tests/test_acceptance.py::test_get_anndata[First 100K cells-mus_musculus]
8.86s call     tests/test_get_anndata.py::test_get_anndata_allows_missing_obs_or_var_filter
8.70s call     tests/test_acceptance.py::test_incremental_read_X[2-None-homo_sapiens]
7.92s call     tests/test_acceptance.py::test_load_axes[homo_sapiens]
7.32s call     tests/test_acceptance.py::test_incremental_query[None-tissue=='brain'-mus_musculus]
7.22s call     tests/test_acceptance.py::test_incremental_query[None-tissue=='aorta'-homo_sapiens]
6.99s call     tests/test_acceptance.py::test_get_anndata[tissue=='aorta'-None-ctx_config0-homo_sapiens]
6.20s call     tests/test_directory.py::test_live_directory_contents
5.57s call     tests/test_get_anndata.py::test_get_anndata_value_filter
4.58s call     tests/test_acceptance.py::test_incremental_query[2-tissue=='aorta'-mus_musculus]
4.53s call     tests/test_acceptance.py::test_get_anndata[First 10K cells-homo_sapiens]
4.19s call     tests/test_open.py::test_get_source_h5ad_uri
4.03s call     tests/test_get_helpers.py::test_get_presence_matrix[homo_sapiens]
4.01s call     tests/test_acceptance.py::test_incremental_query[None-tissue=='aorta'-mus_musculus]
3.98s call     tests/test_acceptance.py::test_get_anndata[First 10K cells-mus_musculus]
3.48s call     tests/test_acceptance.py::test_get_anndata[tissue=='aorta'-None-ctx_config0-mus_musculus]
3.41s call     tests/test_acceptance.py::test_incremental_read_obs[None-ctx_config1-homo_sapiens]
2.24s call     tests/test_get_anndata.py::test_get_anndata_coords
2.19s call     tests/test_acceptance.py::test_incremental_read_X[2-None-mus_musculus]
2.17s call     tests/test_get_helpers.py::test_get_presence_matrix[mus_musculus]
2.04s call     tests/test_acceptance.py::test_load_axes[mus_musculus]
2.01s call     tests/test_open.py::test_download_source_h5ad
1.38s call     tests/test_acceptance.py::test_incremental_read_obs[2-None-homo_sapiens]
1.25s call     tests/test_acceptance.py::test_incremental_read_obs[None-ctx_config1-mus_musculus]
1.24s call     tests/test_acceptance.py::test_incremental_read_obs[2-None-mus_musculus]
1.16s call     tests/test_acceptance.py::test_incremental_read_var[2-None-mus_musculus]
1.15s call     tests/test_acceptance.py::test_incremental_read_var[None-ctx_config1-mus_musculus]
1.12s call     tests/test_acceptance.py::test_incremental_read_var[2-None-homo_sapiens]
1.08s call     tests/test_acceptance.py::test_incremental_read_var[None-ctx_config1-homo_sapiens]
1.01s setup    tests/test_open.py::test_download_source_h5ad_errors
1.01s setup    tests/test_open.py::test_download_source_h5ad
0.90s call     tests/test_open.py::test_get_source_h5ad_uri_errors
0.68s call     tests/test_open.py::test_open_soma_stable
0.64s call     tests/test_open.py::test_open_soma_with_context
0.56s call     tests/test_get_helpers.py::test_get_experiment
0.36s call     tests/test_open.py::test_open_soma_defaults_to_latest_if_missing_stable
0.36s setup    tests/test_get_anndata.py::test_get_anndata_allows_missing_obs_or_var_filter
0.36s setup    tests/test_get_anndata.py::test_get_anndata_coords
0.35s setup    tests/test_get_anndata.py::test_get_anndata_value_filter
0.34s call     tests/test_open.py::test_open_soma_latest
0.32s call     tests/test_directory.py::test_get_census_version_description_errors
0.32s call     tests/test_open.py::test_can_open_with_anonymous_access
0.25s call     tests/test_open.py::test_opening_census_without_anon_access_fails_with_bogus_creds
0.03s setup    tests/test_directory.py::test_get_census_version_directory
0.03s call     tests/test_directory.py::test_get_census_version_directory
0.02s teardown tests/test_acceptance.py::test_get_anndata[cell_type=='neuron'-None-ctx_config8-homo_sapiens]
0.02s teardown tests/test_acceptance.py::test_get_anndata[None-None-ctx_config10-homo_sapiens]
0.01s teardown tests/test_acceptance.py::test_get_anndata[is_primary_data==True-None-ctx_config9-homo_sapiens]
0.01s setup    tests/test_acceptance.py::test_get_anndata[is_primary_data==True-None-ctx_config9-mus_musculus]
0.01s setup    tests/test_acceptance.py::test_get_anndata[None-None-ctx_config10-mus_musculus]
0.01s setup    tests/test_acceptance.py::test_get_anndata[cell_type=='neuron'-None-ctx_config8-mus_musculus]

(131 durations < 0.005s hidden.  Use -vv to show these durations.)
=============== 69 passed, 111 deselected in 15040.59s (4:10:40) ===============
```

### 2023-05-16

- Host: EC2 instance type: `r6id.32xlarge`, all nvme mounted as swap.
- Uname: Linux 5.19.0-1022-aws #23~22.04.1-Ubuntu SMP Fri Mar 17 15:38:24 UTC 2023 x86_64 x86_64 x86_64 GNU/Linux
- Python & census versions:

```python
>>> import cellxgene_census, tiledbsoma
>>> tiledbsoma.show_package_versions()
tiledbsoma.__version__        1.2.3
TileDB-Py tiledb.version()    (0, 21, 3)
TileDB core version           2.15.2
libtiledbsoma version()       libtiledb=2.15.2
python version                3.10.6.final.0
OS version                    Linux 5.19.0-1022-aws
>>> cellxgene_census.__version__
  '1.0.2.dev2+g1598cfd'
>>> cellxgene_census.get_census_version_description('latest')
{'release_date': None, 'release_build': '2023-05-15', 'soma': {'uri': 's3://cellxgene-data-public/cell-census/2023-05-15/soma/', 's3_region': 'us-west-2'}, 'h5ads': {'uri': 's3://cellxgene-data-public/cell-census/2023-05-15/h5ads/', 's3_region': 'us-west-2'}, 'alias': 'latest'}
```

**Pytest output:**

```text
============================= test session starts ==============================
platform linux -- Python 3.10.6, pytest-7.3.1, pluggy-1.0.0 -- /home/ubuntu/venv-cellxgene-census/bin/python3
cachedir: .pytest_cache
rootdir: /home/ubuntu/cellxgene-census/api/python/cellxgene_census
configfile: pyproject.toml
plugins: requests-mock-1.10.0
collecting ... collected 51 items

tests/test_acceptance.py::test_load_axes[homo_sapiens] PASSED            [  1%]
tests/test_acceptance.py::test_load_axes[mus_musculus] PASSED            [  3%]
tests/test_acceptance.py::test_incremental_read[homo_sapiens] PASSED     [  5%]
tests/test_acceptance.py::test_incremental_read[mus_musculus] PASSED     [  7%]
tests/test_acceptance.py::test_incremental_query[2-tissue=='aorta'-homo_sapiens] PASSED [  9%]
tests/test_acceptance.py::test_incremental_query[2-tissue=='aorta'-mus_musculus] PASSED [ 11%]
tests/test_acceptance.py::test_incremental_query[2-tissue=='brain'-homo_sapiens] PASSED [ 13%]
tests/test_acceptance.py::test_incremental_query[2-tissue=='brain'-mus_musculus] PASSED [ 15%]
tests/test_acceptance.py::test_incremental_query[None-tissue=='aorta'-homo_sapiens] PASSED [ 17%]
tests/test_acceptance.py::test_incremental_query[None-tissue=='aorta'-mus_musculus] PASSED [ 19%]
tests/test_acceptance.py::test_incremental_query[None-tissue=='brain'-homo_sapiens] PASSED [ 21%]
tests/test_acceptance.py::test_incremental_query[None-tissue=='brain'-mus_musculus] PASSED [ 23%]
tests/test_acceptance.py::test_get_anndata[tissue=='aorta'-None-ctx_config0-homo_sapiens] PASSED [ 25%]
tests/test_acceptance.py::test_get_anndata[tissue=='aorta'-None-ctx_config0-mus_musculus] PASSED [ 27%]
tests/test_acceptance.py::test_get_anndata[First 10K cells-homo_sapiens] PASSED [ 29%]
tests/test_acceptance.py::test_get_anndata[First 10K cells-mus_musculus] PASSED [ 31%]
tests/test_acceptance.py::test_get_anndata[First 100K cells-homo_sapiens] PASSED [ 33%]
tests/test_acceptance.py::test_get_anndata[First 100K cells-mus_musculus] PASSED [ 35%]
tests/test_acceptance.py::test_get_anndata[First 1M cells-homo_sapiens] PASSED [ 37%]
tests/test_acceptance.py::test_get_anndata[First 1M cells-mus_musculus] PASSED [ 39%]
tests/test_acceptance.py::test_get_anndata[cell_type=='neuron'-None-ctx_config4-homo_sapiens] PASSED [ 41%]
tests/test_acceptance.py::test_get_anndata[cell_type=='neuron'-None-ctx_config4-mus_musculus] PASSED [ 43%]
tests/test_acceptance.py::test_get_anndata[tissue=='brain'-None-ctx_config5-homo_sapiens] PASSED [ 45%]
tests/test_acceptance.py::test_get_anndata[tissue=='brain'-None-ctx_config5-mus_musculus] PASSED [ 47%]
tests/test_acceptance.py::test_get_anndata[is_primary_data==True-None-ctx_config6-homo_sapiens] PASSED [ 49%]
tests/test_acceptance.py::test_get_anndata[is_primary_data==True-None-ctx_config6-mus_musculus] PASSED [ 50%]
tests/test_acceptance.py::test_get_anndata[None-None-ctx_config7-homo_sapiens] PASSED [ 52%]
tests/test_acceptance.py::test_get_anndata[None-None-ctx_config7-mus_musculus] PASSED [ 54%]
tests/test_directory.py::test_get_census_version_directory PASSED        [ 56%]
tests/test_directory.py::test_get_census_version_description_errors PASSED [ 58%]
tests/test_directory.py::test_live_directory_contents PASSED             [ 60%]
tests/test_get_anndata.py::test_get_anndata_value_filter PASSED          [ 62%]
tests/test_get_anndata.py::test_get_anndata_coords PASSED                [ 64%]
tests/test_get_anndata.py::test_get_anndata_allows_missing_obs_or_var_filter PASSED [ 66%]
tests/test_get_helpers.py::test_get_experiment PASSED                    [ 68%]
tests/test_get_helpers.py::test_get_presence_matrix[homo_sapiens] PASSED [ 70%]
tests/test_get_helpers.py::test_get_presence_matrix[mus_musculus] PASSED [ 72%]
tests/test_open.py::test_open_soma_stable PASSED                         [ 74%]
tests/test_open.py::test_open_soma_latest PASSED                         [ 76%]
tests/test_open.py::test_open_soma_with_context PASSED                   [ 78%]
tests/test_open.py::test_open_soma_invalid_args PASSED                   [ 80%]
tests/test_open.py::test_open_soma_errors PASSED                         [ 82%]
tests/test_open.py::test_open_soma_defaults_to_latest_if_missing_stable PASSED [ 84%]
tests/test_open.py::test_open_soma_defaults_to_stable PASSED             [ 86%]
tests/test_open.py::test_get_source_h5ad_uri PASSED                      [ 88%]
tests/test_open.py::test_get_source_h5ad_uri_errors PASSED               [ 90%]
tests/test_open.py::test_download_source_h5ad PASSED                     [ 92%]
tests/test_open.py::test_download_source_h5ad_errors PASSED              [ 94%]
tests/test_open.py::test_opening_census_without_anon_access_fails_with_bogus_creds PASSED [ 96%]
tests/test_open.py::test_can_open_with_anonymous_access PASSED           [ 98%]
tests/test_util.py::test_uri_join PASSED                                 [100%]

============================== slowest durations ===============================
6905.70s call     tests/test_acceptance.py::test_get_anndata[None-None-ctx_config7-homo_sapiens]
2222.81s call     tests/test_acceptance.py::test_get_anndata[is_primary_data==True-None-ctx_config6-homo_sapiens]
743.26s call     tests/test_acceptance.py::test_get_anndata[cell_type=='neuron'-None-ctx_config4-homo_sapiens]
223.36s call     tests/test_acceptance.py::test_get_anndata[None-None-ctx_config7-mus_musculus]
174.85s call     tests/test_acceptance.py::test_get_anndata[is_primary_data==True-None-ctx_config6-mus_musculus]
51.53s call     tests/test_acceptance.py::test_get_anndata[First 1M cells-mus_musculus]
50.58s call     tests/test_acceptance.py::test_get_anndata[tissue=='brain'-None-ctx_config5-homo_sapiens]
39.31s call     tests/test_acceptance.py::test_get_anndata[cell_type=='neuron'-None-ctx_config4-mus_musculus]
37.15s call     tests/test_acceptance.py::test_get_anndata[First 1M cells-homo_sapiens]
34.09s call     tests/test_acceptance.py::test_get_anndata[tissue=='brain'-None-ctx_config5-mus_musculus]
30.23s call     tests/test_acceptance.py::test_incremental_query[None-tissue=='brain'-homo_sapiens]
19.14s call     tests/test_acceptance.py::test_incremental_query[2-tissue=='brain'-homo_sapiens]
14.76s call     tests/test_directory.py::test_live_directory_contents
13.29s call     tests/test_acceptance.py::test_get_anndata[First 100K cells-mus_musculus]
9.64s call     tests/test_acceptance.py::test_incremental_query[2-tissue=='brain'-mus_musculus]
9.48s call     tests/test_acceptance.py::test_get_anndata[First 100K cells-homo_sapiens]
9.48s call     tests/test_acceptance.py::test_incremental_read[homo_sapiens]
9.47s call     tests/test_get_anndata.py::test_get_anndata_allows_missing_obs_or_var_filter
8.13s call     tests/test_acceptance.py::test_get_anndata[tissue=='aorta'-None-ctx_config0-homo_sapiens]
7.94s call     tests/test_acceptance.py::test_incremental_query[None-tissue=='brain'-mus_musculus]
7.82s call     tests/test_acceptance.py::test_load_axes[homo_sapiens]
7.52s call     tests/test_acceptance.py::test_incremental_query[None-tissue=='aorta'-homo_sapiens]
7.32s call     tests/test_acceptance.py::test_incremental_query[2-tissue=='aorta'-homo_sapiens]
5.51s call     tests/test_acceptance.py::test_get_anndata[First 10K cells-homo_sapiens]
5.44s call     tests/test_get_anndata.py::test_get_anndata_value_filter
4.89s call     tests/test_acceptance.py::test_get_anndata[First 10K cells-mus_musculus]
4.84s call     tests/test_acceptance.py::test_get_anndata[tissue=='aorta'-None-ctx_config0-mus_musculus]
4.73s call     tests/test_acceptance.py::test_incremental_read[mus_musculus]
4.17s call     tests/test_open.py::test_get_source_h5ad_uri
3.92s call     tests/test_acceptance.py::test_incremental_query[2-tissue=='aorta'-mus_musculus]
3.59s call     tests/test_acceptance.py::test_incremental_query[None-tissue=='aorta'-mus_musculus]
3.29s call     tests/test_get_helpers.py::test_get_presence_matrix[homo_sapiens]
3.11s call     tests/test_get_anndata.py::test_get_anndata_coords
2.58s call     tests/test_get_helpers.py::test_get_presence_matrix[mus_musculus]
1.99s call     tests/test_acceptance.py::test_load_axes[mus_musculus]
1.89s call     tests/test_open.py::test_download_source_h5ad
1.06s call     tests/test_open.py::test_open_soma_with_context
1.05s setup    tests/test_open.py::test_download_source_h5ad_errors
1.00s call     tests/test_open.py::test_open_soma_stable
0.96s call     tests/test_open.py::test_get_source_h5ad_uri_errors
0.89s setup    tests/test_open.py::test_download_source_h5ad
0.75s call     tests/test_get_helpers.py::test_get_experiment
0.42s setup    tests/test_get_anndata.py::test_get_anndata_value_filter
0.39s call     tests/test_open.py::test_open_soma_defaults_to_latest_if_missing_stable
0.34s setup    tests/test_get_anndata.py::test_get_anndata_allows_missing_obs_or_var_filter
0.34s call     tests/test_open.py::test_can_open_with_anonymous_access
0.33s call     tests/test_open.py::test_open_soma_latest
0.32s call     tests/test_directory.py::test_get_census_version_description_errors
0.31s setup    tests/test_get_anndata.py::test_get_anndata_coords
0.25s call     tests/test_open.py::test_opening_census_without_anon_access_fails_with_bogus_creds
0.04s setup    tests/test_directory.py::test_get_census_version_directory
0.03s call     tests/test_directory.py::test_get_census_version_directory
0.01s teardown tests/test_acceptance.py::test_get_anndata[None-None-ctx_config7-homo_sapiens]
0.01s setup    tests/test_acceptance.py::test_get_anndata[None-None-ctx_config7-mus_musculus]
0.01s teardown tests/test_acceptance.py::test_get_anndata[is_primary_data==True-None-ctx_config6-homo_sapiens]
0.01s setup    tests/test_acceptance.py::test_get_anndata[is_primary_data==True-None-ctx_config6-mus_musculus]

(97 durations < 0.005s hidden.  Use -vv to show these durations.)
======================= 51 passed in 10696.20s (2:58:16) =======================
```

### 2023-03-29

**Config:**

- Host: EC2 instance type: `r6id.32xlarge`, all nvme mounted as swap.
- Uname: Linux bruce.aegea 5.15.0-1033-aws #37~20.04.1-Ubuntu SMP Fri Mar 17 11:39:30 UTC 2023 x86_64 x86_64 x86_64 GNU/Linux
- Python & census versions:

```python
In [1]: import cell_census, tiledbsoma

In [2]: tiledbsoma.show_package_versions()
tiledbsoma.__version__        1.2.1
TileDB-Py tiledb.version()    (0, 21, 1)
TileDB core version           2.15.0
libtiledbsoma version()       libtiledbsoma=;libtiledb=2.15.0
python version                3.9.16.final.0
OS version                    Linux 5.15.0-1033-aws

In [3]: cell_census.get_census_version_description('latest')
Out[3]:
{'release_date': None,
 'release_build': '2023-03-16',
 'soma': {'uri': 's3://cellxgene-data-public/cell-census/2023-03-16/soma/',
  's3_region': 'us-west-2'},
 'h5ads': {'uri': 's3://cellxgene-data-public/cell-census/2023-03-16/h5ads/',
  's3_region': 'us-west-2'}}

In [4]: cell_census.__version__
Out[4]: '0.12.0'
```

**Run notes:**

The test `test_acceptance.py::test_get_anndata[None-homo_sapiens]` manifest a large amount of paging activity.

**Pytest output:**

```text
$ pytest -v --durations=0 --expensive ./api/python/cell_census/tests/
==================================================== test session starts =====================================================
platform linux -- Python 3.9.16, pytest-7.2.2, pluggy-1.0.0 -- /home/bruce/cell-census/venv/bin/python
cachedir: .pytest_cache
rootdir: /home/bruce/cell-census/api/python/cell_census, configfile: pyproject.toml
plugins: requests-mock-1.10.0, anyio-3.6.2
collected 45 items

api/python/cell_census/tests/test_acceptance.py::test_load_axes[homo_sapiens] PASSED                                   [  2%]
api/python/cell_census/tests/test_acceptance.py::test_load_axes[mus_musculus] PASSED                                   [  4%]
api/python/cell_census/tests/test_acceptance.py::test_incremental_read[homo_sapiens] PASSED                            [  6%]
api/python/cell_census/tests/test_acceptance.py::test_incremental_read[mus_musculus] PASSED                            [  8%]
api/python/cell_census/tests/test_acceptance.py::test_incremental_query[2-tissue=='aorta'-homo_sapiens] PASSED         [ 11%]
api/python/cell_census/tests/test_acceptance.py::test_incremental_query[2-tissue=='aorta'-mus_musculus] PASSED         [ 13%]
api/python/cell_census/tests/test_acceptance.py::test_incremental_query[2-tissue=='brain'-homo_sapiens] PASSED         [ 15%]
api/python/cell_census/tests/test_acceptance.py::test_incremental_query[2-tissue=='brain'-mus_musculus] PASSED         [ 17%]
api/python/cell_census/tests/test_acceptance.py::test_incremental_query[None-tissue=='aorta'-homo_sapiens] PASSED      [ 20%]
api/python/cell_census/tests/test_acceptance.py::test_incremental_query[None-tissue=='aorta'-mus_musculus] PASSED      [ 22%]
api/python/cell_census/tests/test_acceptance.py::test_incremental_query[None-tissue=='brain'-homo_sapiens] PASSED      [ 24%]
api/python/cell_census/tests/test_acceptance.py::test_incremental_query[None-tissue=='brain'-mus_musculus] PASSED      [ 26%]
api/python/cell_census/tests/test_acceptance.py::test_get_anndata[tissue=='aorta'-None-ctx_config0-homo_sapiens] PASSED [ 28%]
api/python/cell_census/tests/test_acceptance.py::test_get_anndata[tissue=='aorta'-None-ctx_config0-mus_musculus] PASSED [ 31%]
api/python/cell_census/tests/test_acceptance.py::test_get_anndata[First 10K cells-homo_sapiens] PASSED                 [ 33%]
api/python/cell_census/tests/test_acceptance.py::test_get_anndata[First 10K cells-mus_musculus] PASSED                 [ 35%]
api/python/cell_census/tests/test_acceptance.py::test_get_anndata[First 100K cells-homo_sapiens] PASSED                [ 37%]
api/python/cell_census/tests/test_acceptance.py::test_get_anndata[First 100K cells-mus_musculus] PASSED                [ 40%]
api/python/cell_census/tests/test_acceptance.py::test_get_anndata[First 1M cells-homo_sapiens] PASSED                  [ 42%]
api/python/cell_census/tests/test_acceptance.py::test_get_anndata[First 1M cells-mus_musculus] PASSED                  [ 44%]
api/python/cell_census/tests/test_acceptance.py::test_get_anndata[cell_type=='neuron'-None-ctx_config4-homo_sapiens] PASSED [ 46%]
api/python/cell_census/tests/test_acceptance.py::test_get_anndata[cell_type=='neuron'-None-ctx_config4-mus_musculus] PASSED [ 48%]
api/python/cell_census/tests/test_acceptance.py::test_get_anndata[tissue=='brain'-None-ctx_config5-homo_sapiens] PASSED [ 51%]
api/python/cell_census/tests/test_acceptance.py::test_get_anndata[tissue=='brain'-None-ctx_config5-mus_musculus] PASSED [ 53%]
api/python/cell_census/tests/test_acceptance.py::test_get_anndata[is_primary_data==True-None-ctx_config6-homo_sapiens] PASSED [ 55%]
api/python/cell_census/tests/test_acceptance.py::test_get_anndata[is_primary_data==True-None-ctx_config6-mus_musculus] PASSED [ 57%]
api/python/cell_census/tests/test_acceptance.py::test_get_anndata[None-None-ctx_config7-homo_sapiens] PASSED           [ 60%]
api/python/cell_census/tests/test_acceptance.py::test_get_anndata[None-None-ctx_config7-mus_musculus] PASSED           [ 62%]
api/python/cell_census/tests/test_directory.py::test_get_census_version_directory PASSED                               [ 64%]
api/python/cell_census/tests/test_directory.py::test_get_census_version_description_errors PASSED                      [ 66%]
api/python/cell_census/tests/test_directory.py::test_live_directory_contents PASSED                                    [ 68%]
api/python/cell_census/tests/test_get_anndata.py::test_get_anndata_value_filter PASSED                                 [ 71%]
api/python/cell_census/tests/test_get_anndata.py::test_get_anndata_coords PASSED                                       [ 73%]
api/python/cell_census/tests/test_get_anndata.py::test_get_anndata_allows_missing_obs_or_var_filter PASSED             [ 75%]
api/python/cell_census/tests/test_get_helpers.py::test_get_experiment PASSED                                           [ 77%]
api/python/cell_census/tests/test_get_helpers.py::test_get_presence_matrix[homo_sapiens] PASSED                        [ 80%]
api/python/cell_census/tests/test_get_helpers.py::test_get_presence_matrix[mus_musculus] PASSED                        [ 82%]
api/python/cell_census/tests/test_open.py::test_open_soma_latest PASSED                                                [ 84%]
api/python/cell_census/tests/test_open.py::test_open_soma_with_context PASSED                                          [ 86%]
api/python/cell_census/tests/test_open.py::test_open_soma_errors PASSED                                                [ 88%]
api/python/cell_census/tests/test_open.py::test_get_source_h5ad_uri PASSED                                             [ 91%]
api/python/cell_census/tests/test_open.py::test_get_source_h5ad_uri_errors PASSED                                      [ 93%]
api/python/cell_census/tests/test_open.py::test_download_source_h5ad PASSED                                            [ 95%]
api/python/cell_census/tests/test_open.py::test_download_source_h5ad_errors PASSED                                     [ 97%]
api/python/cell_census/tests/test_util.py::test_uri_join PASSED                                                        [100%]

===================================================== slowest durations ======================================================
5455.14s call     tests/test_acceptance.py::test_get_anndata[None-None-ctx_config7-homo_sapiens]
1388.18s call     tests/test_acceptance.py::test_get_anndata[is_primary_data==True-None-ctx_config6-homo_sapiens]
400.45s call     tests/test_acceptance.py::test_get_anndata[cell_type=='neuron'-None-ctx_config4-homo_sapiens]
183.85s call     tests/test_acceptance.py::test_get_anndata[None-None-ctx_config7-mus_musculus]
110.33s call     tests/test_acceptance.py::test_get_anndata[is_primary_data==True-None-ctx_config6-mus_musculus]
63.52s call     tests/test_acceptance.py::test_get_anndata[First 1M cells-mus_musculus]
44.27s call     tests/test_acceptance.py::test_get_anndata[First 1M cells-homo_sapiens]
35.95s call     tests/test_acceptance.py::test_get_anndata[tissue=='brain'-None-ctx_config5-homo_sapiens]
25.85s call     tests/test_acceptance.py::test_incremental_query[None-tissue=='brain'-homo_sapiens]
24.19s call     tests/test_acceptance.py::test_get_anndata[cell_type=='neuron'-None-ctx_config4-mus_musculus]
22.38s call     tests/test_acceptance.py::test_incremental_query[2-tissue=='brain'-homo_sapiens]
13.23s call     tests/test_acceptance.py::test_get_anndata[tissue=='brain'-None-ctx_config5-mus_musculus]
11.56s call     tests/test_get_anndata.py::test_get_anndata_allows_missing_obs_or_var_filter
9.32s call     tests/test_acceptance.py::test_incremental_query[2-tissue=='brain'-mus_musculus]
9.31s call     tests/test_acceptance.py::test_get_anndata[First 100K cells-homo_sapiens]
8.39s call     tests/test_acceptance.py::test_incremental_read[homo_sapiens]
8.14s call     tests/test_acceptance.py::test_incremental_query[2-tissue=='aorta'-homo_sapiens]
7.60s call     tests/test_acceptance.py::test_get_anndata[tissue=='aorta'-None-ctx_config0-homo_sapiens]
7.25s call     tests/test_acceptance.py::test_incremental_query[None-tissue=='brain'-mus_musculus]
7.25s call     tests/test_acceptance.py::test_incremental_query[None-tissue=='aorta'-homo_sapiens]
7.23s call     tests/test_acceptance.py::test_load_axes[homo_sapiens]
6.91s call     tests/test_acceptance.py::test_get_anndata[First 100K cells-mus_musculus]
6.25s setup    tests/test_open.py::test_download_source_h5ad
5.88s call     tests/test_acceptance.py::test_incremental_query[None-tissue=='aorta'-mus_musculus]
5.58s call     tests/test_acceptance.py::test_incremental_query[2-tissue=='aorta'-mus_musculus]
5.14s call     tests/test_directory.py::test_live_directory_contents
5.13s call     tests/test_acceptance.py::test_get_anndata[First 10K cells-homo_sapiens]
4.89s call     tests/test_open.py::test_get_source_h5ad_uri
4.59s call     tests/test_open.py::test_open_soma_latest
4.35s call     tests/test_acceptance.py::test_incremental_read[mus_musculus]
4.23s call     tests/test_get_anndata.py::test_get_anndata_value_filter
3.96s call     tests/test_acceptance.py::test_get_anndata[tissue=='aorta'-None-ctx_config0-mus_musculus]
3.66s call     tests/test_get_helpers.py::test_get_presence_matrix[homo_sapiens]
3.37s call     tests/test_acceptance.py::test_get_anndata[First 10K cells-mus_musculus]
2.97s call     tests/test_get_helpers.py::test_get_presence_matrix[mus_musculus]
2.62s call     tests/test_get_anndata.py::test_get_anndata_coords
2.35s call     tests/test_open.py::test_download_source_h5ad
2.04s call     tests/test_acceptance.py::test_load_axes[mus_musculus]
1.94s setup    tests/test_get_anndata.py::test_get_anndata_coords
1.21s call     tests/test_open.py::test_get_source_h5ad_uri_errors
0.99s setup    tests/test_open.py::test_download_source_h5ad_errors
0.55s call     tests/test_get_helpers.py::test_get_experiment
0.51s call     tests/test_open.py::test_open_soma_with_context
0.25s setup    tests/test_get_anndata.py::test_get_anndata_value_filter
0.23s setup    tests/test_get_anndata.py::test_get_anndata_allows_missing_obs_or_var_filter
0.06s call     tests/test_directory.py::test_get_census_version_description_errors
0.04s setup    tests/test_directory.py::test_get_census_version_directory
0.02s call     tests/test_directory.py::test_get_census_version_directory
0.01s teardown tests/test_acceptance.py::test_get_anndata[is_primary_data==True-None-ctx_config6-homo_sapiens]
0.01s setup    tests/test_acceptance.py::test_get_anndata[is_primary_data==True-None-ctx_config6-mus_musculus]
0.01s teardown tests/test_acceptance.py::test_get_anndata[None-None-ctx_config7-homo_sapiens]

(84 durations < 0.005s hidden.  Use -vv to show these durations.)
============================================== 45 passed in 7924.13s (2:12:04) ===============================================
```



# Section: cellxgene_census_docsite_data_release_info

# Census data releases

**Last edited**: July 8th, 2024.

**Contents:**

1. [What is a Census data release?](#what-is-a-census-data-release)
1. [List of LTS Census data releases](#list-of-lts-census-data-releases)
1. [Compatibility with package versions](#compatibility-with-package-versions)

## What is a Census data release?

It is a Census build that is publicly hosted online. A Census build is
a [TileDB-SOMA](https://github.com/single-cell-data/TileDB-SOMA) collection with the Census data from [CZ CELLxGENE Discover](https://cellxgene.cziscience.com/) as specified in the [Census schema](cellxgene_census_docsite_schema.md).

Any given Census build is named with a unique tag, normally the date of build, e.g., `"2023-05-15"`.

### Long-term supported (LTS) Census releases

To enable data stability and scientific reproducibility, [CZ CELLxGENE Discover](https://cellxgene.cziscience.com/) plans to perform regular LTS Census data releases:

* Published online every six months for public access, starting on May 15, 2023.
* Available for public access for at least 5 years upon publication.

The most recent LTS Census data release is the default opened by the APIs and recognized as `census_version = "stable"`. To open previous LTS Census data releases, you can directly specify the version via its build date `census_version = "[YYYY]-[MM]-[DD]"`.

Python

```python
import cellxgene_census
census = cellxgene_census.open_soma(census_version = "stable")
```

R

```r
library("cellxgene.census")
census <- open_soma(census_version = "stable")
```

### Weekly Census releases (latest)

[CZ CELLxGENE Discover](https://cellxgene.cziscience.com/) ingests a handful of new datasets every week. To quickly enable access to these new data via the Census, CZ CELLxGENE Discover plans to perform weekly Census data releases:

* Available for public access for 1 month.

The most recent weekly release can be opened by the APIs by specifying `census_version = "latest"`.

Python

```python
import cellxgene_census
census = cellxgene_census.open_soma(census_version = "latest")
```

R

```r
library("cellxgene.census")
census <- open_soma(census_version = "latest")
```

## List of LTS Census data releases

### LTS 2024-07-01

Open this data release by specifying `census_version = "2024-07-01"` in future calls to `open_soma()`.

#### Version information

| Information                       | Value      |
|-----------------------------------|------------|
| Census schema version             | [2.0.1](https://github.com/chanzuckerberg/cellxgene-census/blob/fad674674e5070b735a29bc069d1d3dc21d2e5e8/docs/cellxgene_census_schema.md) |
| Census build date                 | 2024-05-20 |
| Dataset schema version            | [5.0.0](https://github.com/chanzuckerberg/cellxgene-census/blob/fad674674e5070b735a29bc069d1d3dc21d2e5e8/docs/cellxgene_census_schema.md)      |
| Number of datasets                | 812        |

#### Cell and donor counts

| Type              | _Homo sapiens_ | _Mus musculus_ |
|-------------------|----------------|----------------|
| Total cells       | 74,322,510     | 41,233,630     |
| Unique cells      | 44,265,932     | 16,332,034     |
| Number of donors  | 17,651         | 4,216          |

#### Cell metadata

| Category                | _Homo sapiens_ | _Mus musculus_ |
|-------------------------|----------------|----------------|
| Assay                   | 24             | 11             |
| Cell type               | 698            | 364            |
| Development stage       | 176            | 48             |
| Disease                 | 109            | 7              |
| Self-reported ethnicity | 31             | _NA_           |
| Sex                     | 3              | 3              |
| Suspension type         | 2              | 2              |
| Tissue                  | 267            | 84             |
| Tissue general          | 55             | 29             |

#### Embbedings

Find out more in the [Census model page](https://cellxgene.cziscience.com/census-models).

Available embeddings can be accessed via [`cellxgene_census.experimental.get_embedding()`](https://chanzuckerberg.github.io/cellxgene-census/_autosummary/cellxgene_census.experimental.get_embedding.html#cellxgene_census.experimental.get_embedding), or by specifying the `obs_embeddings`/`var_embeddings` field in [`cellxgene_census.get_anndata()`](https://chanzuckerberg.github.io/cellxgene-census/_autosummary/cellxgene_census.get_anndata.html#cellxgene_census.get_anndata).

##### Cells

| Method                    | _Homo sapiens_ | _Mus musculus_ |
|---------------------------|----------------|----------------|
| scVI                      | `scvi`         | `scvi`         |
| Geneformer                | `geneformer`   | _NA_           |

### LTS 2023-12-15

Open this data release by specifying `census_version = "2023-12-15"` in future calls to `open_soma()`.

#### Version information

| Information                       | Value      |
|-----------------------------------|------------|
| Census schema version             | [1.2.0](https://github.com/chanzuckerberg/cellxgene-census/blob/3ff1033135b3a9365c239a9442798d88aae94d03/docs/cellxgene_census_schema.md) |
| Census build date                 | 2023-12-15 |
| Dataset schema version            | [3.1.0](https://github.com/chanzuckerberg/single-cell-curation/blob/8ae36ef3fb5a826511dc657d1b8c6d4a772d32e8/schema/3.1.0/schema.md)      |
| Number of datasets                | 651        |

#### Cell and donor counts

| Type              | _Homo sapiens_ | _Mus musculus_ |
|-------------------|----------------|----------------|
| Total cells       | 62,998,417     | 5,684,805      |
| Unique cells      | 36,227,903     | 4,128,230     |
| Number of donors  | 15,588         | 1,990          |

#### Cell metadata

| Category                | _Homo sapiens_ | _Mus musculus_ |
|-------------------------|----------------|----------------|
| Assay                   | 20             | 10              |
| Cell type               | 631            | 248            |
| Development stage       | 173            | 36             |
| Disease                 | 72             | 5              |
| Self-reported ethnicity | 30             | _NA_           |
| Sex                     | 3              | 3              |
| Suspension type         | 2              | 2              |
| Tissue                  | 230            | 74             |
| Tissue general          | 53             | 27             |

#### Embbedings

Find out more in the [Census model page](https://cellxgene.cziscience.com/census-models).

Available embeddings can be accessed via [`cellxgene_census.experimental.get_embedding()`](https://chanzuckerberg.github.io/cellxgene-census/_autosummary/cellxgene_census.experimental.get_embedding.html#cellxgene_census.experimental.get_embedding), or by specifying the `obs_embeddings`/`var_embeddings` field in [`cellxgene_census.get_anndata()`](https://chanzuckerberg.github.io/cellxgene-census/_autosummary/cellxgene_census.get_anndata.html#cellxgene_census.get_anndata).

##### Cells

| Method                    | _Homo sapiens_ | _Mus musculus_ |
|---------------------------|----------------|----------------|
| scVI                      | `scvi`         | `scvi`         |
| Fine-tuned Geneformer     | `geneformer`   | _NA_           |
| scGPT                     | `scgpt`        | _NA_           |
| Universal Cell Embeddings | `uce`          | _NA_           |
| NMF                       | `nmf`          | _NA_           |

##### Features

| Method                    | _Homo sapiens_ | _Mus musculus_ |
|---------------------------|----------------|----------------|
| NMF                       | `nmf`          | _NA_           |

### LTS 2023-07-25

Open this data release by specifying `census_version = "2023-07-25"` in future calls to `open_soma()`.

#### Version information

| Information                       | Value      |
|-----------------------------------|------------|
| Census schema version             | [1.0.0](https://github.com/chanzuckerberg/cellxgene-census/blob/f06bcebb6471735681fd84734d2d581c44e049e7/docs/cellxgene_census_schema.md) |
| Census build date                 | 2023-07-25 |
| Dataset schema version            | [3.0.0](https://github.com/chanzuckerberg/single-cell-curation/blob/a64ac9eb70e3e777ee34098ae82120c2d21692b0/schema/3.0.0/schema.md)      |
| Number of datasets                | 593        |

#### Cell and donor counts

| Type              | _Homo sapiens_ | _Mus musculus_ |
|-------------------|----------------|----------------|
| Total cells       | 56,400,873     | 5,255,245      |
| Unique cells      | 33,364,242     | 4,083,531     |
| Number of donors  | 13,035         | 1,417          |

#### Cell metadata

| Category                | _Homo sapiens_ | _Mus musculus_ |
|-------------------------|----------------|----------------|
| Assay                   | 19             | 9              |
| Cell type               | 613            | 248            |
| Development stage       | 164            | 33             |
| Disease                 | 64             | 5              |
| Self-reported ethnicity | 26             | _NA_           |
| Sex                     | 3              | 3              |
| Suspension type         | 2              | 2              |
| Tissue                  | 220            | 66             |
| Tissue general          | 54             | 27             |

### LTS 2023-05-15

Open this data release by specifying `census_version = "2023-05-15"` in future calls to `open_soma()`.

#### 🔴 Errata 🔴

##### Duplicate observations with  `is_primary_data = True`

In order to prevent duplicate data in analyses, each observation (cell) should be marked `is_primary data = True` exactly once in the Census. Since this LTS release, 243,569 observations have been identified that are represented at least twice with `is_primary_data = True`.

This issue will be corrected in the following LTS data release, by identifying and marking only one cell out of the duplicates as  `is_primary_data = True`.

If you wish to use this data release, you can consider filtering out all of these 243,569 cells by using the `soma_joinids` provided in this file [duplicate_cells_census_LTS_2023-05-15.csv.zip](https://github.com/chanzuckerberg/cellxgene-census/raw/773edab79bbdc78eccb26ec4f8211a9b4c98a71a/tools/cell_dup_check/duplicate_cells_census_LTS_2023-05-15.csv.zip). You can filter specific cells by using the `value_filter` or `obs_value_filter` of the querying API functions, for more information follow this [tutorial](https://chanzuckerberg.github.io/cellxgene-census/notebooks/api_demo/census_query_extract.html).

#### Version information

| Information                       | Value      |
|-----------------------------------|------------|
| Census schema version             | [1.0.0](https://github.com/chanzuckerberg/cellxgene-census/blob/f06bcebb6471735681fd84734d2d581c44e049e7/docs/cellxgene_census_schema.md) |
| Census build date                 | 2023-05-15 |
| Dataset schema version            | [3.0.0](https://github.com/chanzuckerberg/single-cell-curation/blob/a64ac9eb70e3e777ee34098ae82120c2d21692b0/schema/3.0.0/schema.md)      |
| Number of datasets                | 562        |

#### Cell and donor counts

| Type              | _Homo sapiens_ | _Mus musculus_ |
|-------------------|----------------|----------------|
| Total cells       | 53,794,728     | 4,086,032      |
| Unique cells      | 33,758,887     | 2,914,318      |
| Number of donors  | 12,493         | 1,362          |

#### Cell metadata

| Category                | _Homo sapiens_ | _Mus musculus_ |
|-------------------------|----------------|----------------|
| Assay                   | 20             | 9              |
| Cell type               | 604            | 226            |
| Development stage       | 164            | 30             |
| Disease                 | 68             | 5              |
| Self-reported ethnicity | 26             | _NA_           |
| Sex                     | 3              | 3              |
| Suspension type         | 2              | 2              |
| Tissue                  | 227            | 51             |
| Tissue general          | 61             | 27             |

## Compatibility with package versions

Due to the nature of the Census storage backend, the format version will change from time to time. Format upgrades are always backwards compatible, but they're not always forwards compatible, which means that reading a recent Census data version using an older version of the package might result in an error.
We aim to guarantee the following policy:

* Every Census package version released after an LTS will be able to read _every_ Census data release until the next LTS.

The current LTS release (2023-12-15) is compatible with the following package versions:

* 1.10.x
* 1.11.x
* 1.12.x
* 1.13.x



# Section: cellxgene_census_docsite_FAQ

# FAQ

Last updated: Sept, 2024.

## Why should I use the Census?

The Census provides efficient low-latency access via Python and R APIs to most single-cell RNA data from [CZ CELLxGENE Discover](https://cellxgene.cziscience.com/). To accelerate computational research, the Census enables researchers to:

- Access slices of data from more than 500 single-cell datasets spanning about 33M unique cells (50M total) from >60K genes from human or mice.
- Access to data with standardized cell and gene metadata with harmonized labels.
- Easily load multi-dataset slices into Scanpy or Seurat.
- Implement out-of-core (a.k.a online) operations for larger-than-memory processes.

For example, a user can easily get “*all T-cells from Lung with COVID-19*” into [AnnData](https://anndata.readthedocs.io/en/latest/), [Seurat](https://satijalab.org/seurat/), or into memory-sufficient data chunks via [PyArrow](https://arrow.apache.org/docs/python/index.html) or [R Arrow](https://arrow.apache.org/docs/r/).

The Census is not suited for:

- Access to non-standardized cell metadata and gene metadata available in the original [datasets](https://cellxgene.cziscience.com/datasets).
- Access to the author-contributed normalized expression values or embeddings.
- Access to all data from just one dataset.
- Access to non-RNA or spatial data present in CZ CELLxGENE Discover as it is not yet supported in the Census.

If you’d like to perform any of the above tasks, you can access web downloads directly from the [CZ CELLxGENE Discover Datasets](https://cellxgene.cziscience.com/datasets) feature. [Click here](https://cellxgene.cziscience.com/docs/03__Download%20Published%20Data) for more information about downloading published data on CELLxGENE Discover.

## What data is contained in the Census?

Most RNA non-spatial data from [CZ CELLxGENE Discover](https://cellxgene.cziscience.com/) is included. You can see a general description of these data and their organization in the [schema description](cellxgene_census_docsite_schema.md) or you can use the APIs to explore the data as indicated in this [tutorial](notebooks/analysis_demo/comp_bio_census_info.ipynb).

## How do I cite the use of the Census for a publication?

Please follow the [citation guidelines](https://cellxgene.cziscience.com/docs/08__Cite%20cellxgene%20in%20your%20publications) offered by CZ CELLxGENE Discover.

## Why does the Census not have a normalized layer or embeddings?

The Census does not have normalized counts or embeddings because:

- The original normalized values and embeddings are not harmonized or integrated across datasets and are therefore numerically incompatible.
- We have not implemented a general-purpose normalization or embedding generation method to be used across all Census data.

If you have any suggestions for methods that our team should explore, please share them with us via a [feature request in the github repository](https://github.com/chanzuckerberg/cellxgene-census/issues/new?assignees=&labels=user+request&template=feature-request.md&title=).

## How does the Census differentiate from other tools?

The Census differentiates from existing single-cell tools by providing fast, efficient access to the largest corpus of standardized single-cell data from CZ CELLxGENE Discover via [TileDB-SOMA](https://github.com/single-cell-data/TileDB-SOMA/issues/new/choose).  Thus, single-cell data from about 33M unique cells (50M total) across >60 K genes, with 11 standardized cell metadata variables and harmonized GENCODE annotations are ready for:

- Opening and reading data at low latency from the cloud.
- Querying and accessing data using metadata filters.
- Loading and creating AnnData objects.
- Loading and creating Seurat objects.
- From Python, creating PyArrow objects, SciPy sparse matrices, NumPy arrays, and Pandas data frames.
- From R, creating R Arrow objects, sparse matrices (via the Matrix package), and standard data frames and (dense) matrices.

## Can I query human and mouse data in a single query?

It is not possible to query both mouse and human data in a single query. This is due to the data from these organisms using different [organism-specific gene annotations](https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/3.0.0/schema.md#required-gene-annotations).

## Where are the Census data hosted?

The Census data is publicly hosted free-of-cost in an Amazon Web Services (AWS) S3 bucket in the [`us-west-2` region](https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/using-regions-availability-zones.html#concepts-available-regions).

## Can I retrieve the original H5AD datasets from which the Census was built?

Yes, you can use the API function `download_source_h5ad` to do so. For usage, please see the reference documentation at the [doc-site](https://chanzuckerberg.github.io/cellxgene-census/) or directly from Python or R:

Python

```python
import cellxgene_census
help(cellxgene_census.download_source_h5ad)
```

R

```r
library(cellxgene.census)
?download_source_h5ad
```

## How can I increase the performance of my queries?

Since the access patterns are via the internet, usually the main limiting step for data queries is bandwidth and client location. We recommend the following tactics to increase query efficiency:

- Utilize a computer connected to high-speed internet.
- Utilize an ethernet connection and not a wifi connection.
- If possible utilize online computing located in the west coast of the US.
- Highly recommended: [EC2 AWS instances](https://aws.amazon.com/ec2/) in the `us-west-2` region.

## Can I use conda to install the Census Python API?

There is not a conda package available for `cellxgene-census`. However you can use conda in combination with `pip` to install the package in a conda environment:

```bash
conda create -n census_env python=3.10
conda activate census_env
pip install cellxgene-census
```

## How can I ask for support?

You can either submit a [github issue](https://github.com/chanzuckerberg/cellxgene-census/issues/new/choose), or for quick support, you can join the CZI Science Community on Slack ([czi.co/science-slack](https://czi.co/science-slack)) and ask questions in the `#cellxgene-census-users` channel.

## How can I ask for new features?

You can submit a [feature request in the github repository](https://github.com/chanzuckerberg/cellxgene-census/issues/new?assignees=&labels=user+request&template=feature-request.md&title=).

## How can I contribute my data to the Census?

To inquire about submitting your data to CZ CELLxGENE Discover, [click here](https://cellxgene.cziscience.com/docs/032__Contribute%20and%20Publish%20Data). If your data request is accepted, the data will automatically be included in the Census if it meets the [biological criteria defined in the Census schema](https://github.com/chanzuckerberg/cellxgene-census/blob/main/docs/cellxgene_census_schema.md#data-included).

## Why do I get an `ArraySchema` error when opening the Census?

You may get this error if you are trying to open a Census data build with an old version of the Census API. Please update your Python or R Census package.

If the error persists please file a [github issue](https://github.com/chanzuckerberg/cellxgene-census/issues/new/choose).

## Why do I get an error when running `import cellxgene_census` on Databricks?

This can occur if the `cellxgene_census` Python package is installed in a Databricks notebook using `%sh pip install cellxgene_census`. This command does *not* restart the Python process after installing `cellxgene_census` and any pip package dependencies that were pre-installed by the Databricks Runtime environment but upgraded for `cellxgene_census` will not be reloaded with their new version. You may see `numba` or `pyarrow` related errors, for example.

To fix, simply install using one of the following Databricks notebook "magic" commands:

```shell
pip install -U cellxgene-census
```

or

```shell
%pip install -U cellxgene-census
```

These commands restart the Python process after installing the `cellxgene-census` package, similar to using `dbutils.library.restartPython()`. Additionally, these magic commands also ensure that the package is installed on all nodes of a multi-node cluster.

See also:

- <https://docs.databricks.com/libraries/notebooks-python-libraries.html#can-i-use-sh-pip-pip-or-pip-what-is-the-difference>
- <https://community.databricks.com/s/question/0D53f00001GHVP3CAP/whats-the-difference-between-magic-commands-pip-and-sh-pip>

Alternately, you can configure your cluster to install the `cellxgene-census` package each time it is started by adding this package to the "Libraries" tab on the cluster configuration page per these [instructions](https://docs.databricks.com/libraries/cluster-libraries.html).

## How do I connect to census from behind a proxy?

TileDB doesn't use the typical proxy environment variables and you'll need to specify these directly. That looks like:

```python
# Replace the ellipses with your proxy host and port info
config = {
    "vfs.s3.proxy_host": ..., "vfs.s3.proxy_port": ...
}

census = cellxgene_census.open_soma(tiledb_config=config)
```

It may not be obvious that a proxy is the issue. This will typically manifest as a `TileDBError` which says that a timeout was reached during a request to the s3 bucket.

You can read more about how to configure how TileDB communicates with S3 [here](https://docs.tiledb.com/main/how-to/backends/s3#aws-security-credentials).



# Section: articles-2024-20240404-categoricals

# Census supports categoricals for cell metadata

*Published:* *April 4th, 2024*

*By:* *[Emanuele Bezzi](mailto:ebezzi@chanzuckerberg.com)* & [Pablo Garcia-Nieto](mailto:pgarcia-nieto@chanzuckerberg.com)

Starting with the `2024-04-01` Census build, a subset of the columns in the `obs` dataframe are now categorical instead of strings.

Overall users will observe a smaller memory footprint when loading Census data into memory. 🚀

However, this may break some existing pipelines as explained below.

## Potential breaking changes

For **Python users**, note that Pandas will encode these columns as `pandas.Categorical`  for which some downstream operations may need to be adapted. See [this link](https://pandas.pydata.org/docs/user_guide/categorical.html#operations) for more details. In particular:

> Series methods like Series.value_counts() will use all categories, even if some categories are not present in the data

and

> DataFrame methods like sum, groupby, pivot, value_counts also show “unused” categories when observed=False, which is the default.

For **R users**, note that these columns will be encoded as `factor` and similarly downstream operations may need to be adapted. See [this link](https://r4ds.had.co.nz/factors.html) for more details.

For **Python and R users** interfacing with `arrow`, these columns will be encoded as `dictionary`, see more details for R in [this link](https://arrow.apache.org/docs/r/reference/dictionary.html) and Python in [this link](https://arrow.apache.org/docs/python/generated/pyarrow.dictionary.html).

## Identifying the `obs` columns encoded as categorical

Users can always check the the type of each cell metadata variable by inspecting the schema of `obs`. Categoricals will be shown as `dictionary`.

In Python:

```python
import cellxgene_census
census = cellxgene_census.open_soma(census_version="latest")
census["census_data"]["homo_sapiens"].obs.schema

# soma_joinid: int64 not null
# dataset_id: dictionary<values=string, indices=int16, ordered=0> not null
# assay: dictionary<values=string, indices=int8, ordered=0> not null
# assay_ontology_term_id: dictionary<values=string, indices=int8, ordered=0> not null
# cell_type: dictionary<values=string, indices=int16, ordered=0> not null
# cell_type_ontology_term_id: dictionary<values=string, indices=int16, ordered=0> not null
# development_stage: dictionary<values=string, indices=int16, ordered=0> not null
# development_stage_ontology_term_id: dictionary<values=string, indices=int16, 
# [OUTPUT TRUNCATED]
```

In R:

```r
library("cellxgene.census")
census = open_soma(census_version="latest")
census$get("census_data")$get("homo_sapiens")$obs$schema()

# Schema
# soma_joinid: int64 not null
# dataset_id: dictionary<values=string, indices=int16> not null
# assay: dictionary<values=string, indices=int8> not null
# assay_ontology_term_id: dictionary<values=string, indices=int8> not null
# cell_type: dictionary<values=string, indices=int16> not null
# cell_type_ontology_term_id: dictionary<values=string, indices=int16> not null
# development_stage: dictionary<values=string, indices=int16> not null
# development_stage_ontology_term_id: dictionary<values=string, indices=int16> not null
# [OUTPUT TRUNCATED]
```



# Section: cellxgene_census_docsite_installation

# Installation

## Requirements

The Census API requires a Linux or MacOS system with:

- Python 3.10 to Python 3.12. Or R, supported versions TBD.
- Recommended: >16 GB of memory.
- Recommended: >5 Mbps internet connection.
- Recommended: for increased performance use the API through a AWS-EC2 instance from the region `us-west-2`. The Census data builds are hosted in a AWS-S3 bucket in that region.

## Python

(Optional) In your working directory, make and activate a virtual environment or conda environment. For example:

```shell
python -m venv ./venv
source ./venv/bin/activate
```

Install the `cellxgene-census` package via pip:

```shell
pip install -U cellxgene-census
```

There are also "experimental" add-on modules that are less stable than the main API, and may have more complex dependencies. To install these,

```shell
pip install -U cellxgene-census[experimental]
```

If installing in a Databricks notebook environment, use `%pip install`. Do not use `%sh pip install`. See the [FAQ](cellxgene_census_docsite_FAQ.md#why-do-i-get-an-error-when-running-import-cellxgene-census-on-databricks).

## R

If installing from **Ubuntu**, you may need to install the following libraries via `apt install`,  `libxml2-dev` `libssl-dev` `libcurl4-openssl-dev`. In addition you must have `cmake` v3.21 or greater.

If installing from **MacOS**, you will need to install the [developer tools `Xcode`](https://apps.apple.com/us/app/xcode/id497799835?mt=12).

**Windows** is not supported.

From an R session, first install `tiledb` from R-Universe, the latest release in CRAN is not yet available.

```r
install.packages(
  "cellxgene.census",
  repos=c('https://chanzuckerberg.r-universe.dev', 'https://cloud.r-project.org')
)
```

To be able to export Census data to `Seurat` or `SingleCellExperiment` you also need to install their respective packages.

```r
# Seurat
install.packages("Seurat")

# SingleCellExperiment
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("SingleCellExperiment")
```



# Section: articles-2024-20240709-pytorch

# First stable iteration of Census (SOMA) PyTorch loaders

*Published:* *July 11th, 2024*

*Updated:* *July 19th, 2024*. Figure 3 has been improved for readability.

*By:* *[Emanuele Bezzi](mailto:ebezzi@chanzuckerberg.com), [Pablo Garcia-Nieto](mailto:pgarcia-nieto@chanzuckerberg.com), [Prathap Sridharan](mailto:psridharan@chanzuckerberg.com), [Ryan Williams](mailto:ryan.williams@tiledb.com)*

The Census team is excited to share the release of Census PyTorch loaders that work out-of-the-box for memory-efficient training across any slice of the >70M cells in Census.

In 2023, we released a beta version of the loaders and we have observed interest from users to utilize them with Census or their own data. For example [Wolf et al.](https://lamin.ai/blog/arrayloader-benchmarks) performed comparisons across different training approaches and found our loaders to be ideal for *uncached* training of Census data, albeit with some caveats.

We have continued the development of the loaders in collaboration with our partners at TileDB, and we are happy to announce this release as the first stable iteration. We hope the loaders can accelerate the development of large-scale models of single-cell data by leveraging the following main features:

- **Out-of-the-box training on all or any slice of Census data.**
- **Efficient memory usage with out-of-core training.**
- **Calibrated shuffling of observations (cells).**
- **Cloud-based or local data access.**
- **Increased training speed.**
- **Custom data encoders.**

Keep on reading for usage and more details on the main loader features.

## Census PyTorch loaders usage

The loaders are ready to use for PyTorch modeling via the specialized Data Pipe [`ExperimentDataPipe`](https://chanzuckerberg.github.io/cellxgene-census/_autosummary/cellxgene_census.experimental.ml.pytorch.ExperimentDataPipe.html#cellxgene_census.experimental.ml.pytorch.ExperimentDataPipe), which takes advantage of the out-of-core data access TileDB-SOMA offers.

Please follow the [Training a PyTorch Model](https://chanzuckerberg.github.io/cellxgene-census/notebooks/experimental/pytorch.html) tutorial for a full reproducible example to train a logistic regression on cell type labels.

In short, the following shows you how to initialize the loader to train a model on a small subset of cells. First, you can initialize a `ExperimentDataPipe` to train a model on tongue cells as follows:

```python
import cellxgene_census.experimental.ml as census_ml
import cellxgene_census
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
)
```

Then you can perform any PyTorch operations and training.

```python
# Splitting training and test sets
train_datapipe, test_datapipe = experiment_datapipe.random_split(weights={"train": 0.8, "test": 0.2}, seed=1)

# Creating data loader
experiment_dataloader = census_ml.experiment_dataloader(train_datapipe)

# Training a PyTorch model
device = torch.device("cuda") if torch.cuda.is_available() else torch.device("cpu")
model = MODEL().to(device)
model.train()
```

## Census PyTorch loaders main features

### Out-of-the-box training on all or any slice of Census data

Since the `ExperimentDataPipe` inherits from the [PyTorch Iterable-style DataPipe](https://pytorch.org/data/main/torchdata.datapipes.iter.html) it can be readily used with PyTorch models.

The single-cell expression data is encoded in numerical tensors, and for supervised training the cell metadata can be automatically transformed with a default encoder, or with custom user-defined encoders (see below).

### Efficient memory usage with out-of-core training

Thanks to the underlying backend of Census — TileDB-SOMA — the PyTorch loaders take advantage of incremental data materialization of fixed and small size to keep memory usage constant throughout training.

In addition, data is eagerly fetched while batches go through training so that compute is never idle or waiting for data to be loaded. This feature is particularly useful when fetching Census data directly from the cloud.

Memory usage is defined by the parameters `soma_chunk_size` and `shuffle_chunk_count` - see below for a full description on how these should be tuned.

### Calibrated shuffling of observations (cells)

Shuffling along efficient out-of-core data fetching is a challenge. In general, increasing randomness of shuffling leads to slower data fetching.

In the first iteration of the loaders, shuffling was done through large blocks of data of user-defined size. This shuffling strategy led to non-random distribution of observations per training batch, becasue Census has a non-random data structure (observations from the same datasets are adjacent to one another) thus training loss was unstable (Figure 1).

**Now we have implemented a scatter-gather approach**, whereby multiple chunks of data are fetched randomly from Census, then a number of chunks are concatenated into a block and all observations within the block are randomly shuffled. Adjusting the size and number of chunks per block leads to well-calibrated shuffling with stable training loss (Figure 2) while maintaining efficient data fetching (Figure 3).

The balance between memory usage, efficiency, and level of randomness can be adjusted with the parameters `soma_chunk_size` and `shuffle_chunk_count`. Increasing `shuffle_chunk_count` will improve randomness, as more scattered chunks will be collected before the pool is randomized. Increasing `soma_chunk_size` will improve I/O efficiency while decreasing it will improve memory usage. We recommend a default of `soma_chunk_size=64, shuffle_chunk_count=2000` as we determined this configuration yields a good balance.

```{figure} ./20240709-pytorch-fig-loss-before.png
:alt: Census PyTorch loaders shuffling
:align: center
:figwidth: 80%

**Figure 1. Training loss was unstable with the previous shuffling strategy**. Based on a trial scVI run on 64K Census cells.
```

```{figure} ./20240709-pytorch-fig-loss-after.png
:alt: Census PyTorch loaders callibrated shuffling
:align: center
:figwidth: 80%

**Figure 2. Training loss is well-calibrated with the current scatter-gather shuffling strategy.** Based on a trial scVI run on 250K Census cells.
```

### Increased training speed

We have made improvements to the loaders to reduce the amount of data transformations required from data fetching to model training. One such important change is to encode the expression data as a dense matrix immediately after the data is retrieved from disk/cloud.

In our benchmarks, we found that densifying data increases training speed while maintaining relatively constant memory usage (Figure 3). For this reason, we have disabled the intermediate data processing in sparse format unless Torch Sparse Tensors are requested via the `ExperimentDataPipe` parameter `return_sparse_X`.

```{figure} ./20240709-pytorch-fig-benchmark.png
:alt: Census PyTorch loaders benchmark
:align: center
:figwidth: 80%

**Figure 3. Benchmark of memory usage and speed of data processing during modeling, default parameters lead to ≈2,500 samples/sec with 27GB of memory use.** The benchmark was done processing 4M cells out of a 10M-cell Census, with data streamed from the cloud (S3). "Method" indicates the expression matrix encoding: circles are dense ("np.array", now the default behavior) and squares are sparse ("scipy.csr"). Size indicates the total number of cells per processing block (max cells materialized at any given time) and color is the number of individual randomly grabbed chunks composing a processing block; higher chunks per block lead to better shuffling. Data was fetched until modeling step, but no model was trained.
```

We repeated the benchmark in Figure 3 in different conditions encompassing varying number of total cells and multiple epochs, please [follow this link for the full benchmark report and code.](https://github.com/ryan-williams/arrayloader-benchmarks).

When comparing dense vs sparse processing in an end-to-end training exercise with scVI, we also observed slight increased speed with the dense approach and comparable memory usage to sparse processing (Figure 4). However in this full training example the differences were less substantial, highlighting that other model-specific factors during the training phase will contribute to memory and speed performance.

```{figure} ./20240709-pytorch-fig-scvi.png
:alt: Census scVI PyTorch run
:align: center
:figwidth: 80%

**Figure 4. Trial scVI training run with default parameters of the Census Pytorch loaders, highlighting increased speed of dense vs sparse data processing.** Training was done on 5684805 mouse cells for 1 epoch on a g4dn.16xlarge EC2 machine.
```

### Custom data encoders

For maximum flexibility, users can provide custom encoders for the cell metadata enabling custom transformations or interactions between different metadata variables.

To use custom encoders you need to instantiate the desired encoder via the [Encoder](https://chanzuckerberg.github.io/cellxgene-census/_autosummary/cellxgene_census.experimental.ml.encoders.Encoder.html#cellxgene_census.experimental.ml.encoders.Encoder) class and pass it to the `encoders` parameter of the `ExperimentDataPipe`.



# Section: cellxgene_census_schema

# CZ CELLxGENE Discover Census Schema

**Version**: 2.1.0

**Last edited**: June, 2024.

The key words "MUST", "MUST NOT", "REQUIRED", "SHALL", "SHALL NOT", "SHOULD", "SHOULD NOT", "RECOMMENDED", "NOT RECOMMENDED" "MAY", and "OPTIONAL" in this document are to be interpreted as described in [BCP 14](https://tools.ietf.org/html/bcp14), [RFC2119](https://www.rfc-editor.org/rfc/rfc2119.txt), and [RFC8174](https://www.rfc-editor.org/rfc/rfc8174.txt) when, and only when, they appear in all capitals, as shown here.

## Census overview

The CZ CELLxGENE Discover Census, hereafter referred as Census, is a versioned data object and API for most of the single-cell data hosted at [CZ CELLxGENE Discover](https://cellxgene.cziscience.com/). To learn more about the Census visit the `chanzuckerberg/cellxgene-census` [github repository](https://github.com/chanzuckerberg/cellxgene-census)

To better understand this document the reader should be familiar with the [CELLxGENE dataset schema](https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/5.1.0/schema.md) and [SOMA](https://github.com/single-cell-data/SOMA/blob/main/abstract_specification.md).

## Definitions

The following terms are used throughout this document:

* adata – generic variable name that refers to an [`AnnData`](https://anndata.readthedocs.io/) object.
* CELLxGENE dataset schema – the data schema for h5ad files served by CELLxGENE Discover, for this Census schema: [CELLxGENE dataset schema version is 5.1.0](https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/5.1.0/schema.md)
* census\_obj – the Census root object, a SOMACollection.
* Census data release – a versioned Census object deposited in a public bucket and accessible by APIs.
* tissue – original tissue annotation.
* tissue\_general – high-level mapping of a tissue, e.g. "Heart" is the tissue_general of "Heart left ventricle" .

## Census Schema versioning

The Census Schema follows [Semver](https://semver.org/) for its versioning:

* Major: any schema changes that make the Census incompatible with the Census API or SOMA API. Examples:
  * Column deletion in Census `obs`
  * Addition of new modality
* Minor: schema additions that are compatible with public API(s) and SOMA. Examples:
  * New column to Census `obs` is added
  * tissue/tissue_general mapping changes
* Patch: schema fixes. Examples:
  * Editorial schema changes

Changes MUST be documented in the schema [Changelog](#changelog) at the end of this document.

Census data releases are versioned separately from the schema.

## Schema

### Data included

All datasets included in the Census MUST be of [CELLxGENE dataset schema version 5.1.0](https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/5.1.0/schema.md). The following data constraints are imposed on top of the CELLxGENE dataset schema.

#### Species

The Census MUST only contain observations (cells) with an  [`organism_ontology_term_id`](https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/5.1.0/schema.md#organism_ontology_term_id) value of either "NCBITaxon:10090" for *Mus musculus* or "NCBITaxon:9606" for *Homo sapiens* MUST be included.

The Census MUST only contain features (genes) with a [`feature_reference`](https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/5.1.0/schema.md#feature_reference) value of either "NCBITaxon:10090" for *Mus musculus* or "NCBITaxon:9606" for *Homo sapiens* MUST be included

#### Multi-species data constraints

Per the CELLxGENE dataset schema, [multi-species datasets MAY contain observations (cells) of a given organism and features (genes) of a different one](https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/5.1.0/schema.md#general-requirements), as defined in [`organism_ontology_term_id`](https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/5.1.0/schema.md#organism_ontology_term_id) and [`feature_reference`](https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/5.1.0/schema.md#feature_reference) respectively.

For any given multi-species dataset, observation and features from the dataset are included in the Census as defined by the following:

* Where a dataset includes observations and features from a single species, all observations and features from the dataset are included in the Census.
* Where a dataset includes observations from a single species `S`, and includes features from multiple species *including* the species `S`, all dataset observations and all features from `S` will be included in the Census.
* Where a dataset includes features from a single species `S`, and observations from multiple species *including* the species `S`, all dataset features and all observations from species `S` are included in the Census.
* Where a species has observations *AND* features from multiple species, the dataset will be excluded from the Census.

The table below shows all possible combinations of organisms for both observations and features, assuming a Census comprised of Homo sapiens and Mus musculus. For each combination, inclusion criteria for the Census is provided.

<table>
<thead>
  <tr>
    <th>Observations (cells) from</th>
    <th>Features (genes) from</th>
    <th>Inclusion criteria</th>
  </tr>
</thead>
<tbody>
  <tr>
    <td>"NCBITaxon:9606" for <i>Homo sapiens</i></td>
    <td>"NCBITaxon:9606" for <i>Homo sapiens</i></td>
    <td>All observations and all features are included.</td>
  </tr>
  <tr>
    <td>"NCBITaxon:10090" for <i>Mus musculus</i></td>
    <td>"NCBITaxon:10090" for <i>Mus musculus</i></td>
    <td>All observations and all features are included.</td>
  </tr>
  <tr>
    <td>"NCBITaxon:9606" for <i>Homo sapiens</i> <b>AND</b> "NCBITaxon:10090" for <i>Mus musculus</i></td>
    <td>"NCBITaxon:9606" for Homo sapiens</td>
    <td>The Census MUST only contain observations from "NCBITaxon:9606" for <i>Homo sapiens</i>. All features MUST be included.</td>
  </tr>
  <tr>
    <td>"NCBITaxon:9606" for <i>Homo sapiens</i> <b>AND</b> "NCBITaxon:10090" for <i>Mus musculus</i></td>
    <td>"NCBITaxon:10090" for <i>Mus musculus</i></td>
    <td>The Census MUST only contain observations from "NCBITaxon:10090" for <i>Mus musculus</i>. All features MUST be included.</td>
  </tr>
  <tr>
    <td>"NCBITaxon:9606" for <i>Homo sapiens</i></td>
    <td>"NCBITaxon:9606" for <i>Homo sapiens</i> <b>AND</b> "NCBITaxon:10090" for <i>Mus musculus</i></td>
    <td>All observations MUST be included. The Census MUST only contain features from "NCBITaxon:9606" for <i>Homo sapiens</i>.</td>
  </tr>
  <tr>
    <td>"NCBITaxon:10090" for <i>Mus musculus</i></td>
    <td>"NCBITaxon:9606" for <i>Homo sapiens</i> <b>AND</b> "NCBITaxon:10090" for <i>Mus musculus</i></td>
    <td>All observations MUST be included. The Census MUST only contain features from "NCBITaxon:10090" for <i>Mus musculus</i>.</td>
  </tr>
  <tr>
    <td>"NCBITaxon:9606" for <i>Homo sapiens</i> <b>AND</b> "NCBITaxon:10090" for <i>Mus musculus</i></td>
    <td>"NCBITaxon:9606" for <i>Homo sapiens</i> <b>AND</b> "NCBITaxon:10090" for <i>Mus musculus</i></td>
    <td>All observations and features MUST NOT be included.</td>
  </tr>
</tbody>
</table>

#### Assays

Assays are defined in the CELLxGENE dataset schema in [`assay_ontology_term_id`](https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/5.1.0/schema.md#assay_ontology_term_id).

The Census MUST include all cells from the list of [accepted assays](./census_accepted_assays.csv).

These assays were selected with the following criteria:

> Only children "EFO:0002772" or "EFO:0010183" are shown as this is a constraint imposed by the CELLxGENE dataset schema >3.0.0.
>
> * Must measure gene expression via RNA sequencing.
> * Can be done at the single-cell level.
> * May include nascent or elongating RNA data.
> * May be targeted to specific genes in an assay-specific manner.
> * Doesn't measure other non-RNA molecules concurrently.
> * Doesn’t measure spatial information.
> * Doesn’t require author metadata for correct interpretability (e.g. perturbation-based technologies).
> * Doesn’t intend to primarily measure RNA structure, RNA fusions, RNA modifications, or RNA interactions.
> * Doesn’t intend to primarily measure non-mRNA (e.g. tRNA, rRNA, small RNAs).
> * Doesn’t intend to primarily measure viral RNA.
> * Doesn’t intend to primarily measure introns.
> * Doesn’t do ribosome profiling.

##### Full-gene sequencing assays

From the list of accepted assays, this list of [full-gene sequencing assays](./census_accepted_assays_full_gene.csv) are those that when used at the single-cell level will always perform full-gene sequencing.

These data need to be normalized by gene length for downstream analysis.

#### Data matrix types

Per the CELLxGENE dataset schema, [all RNA assays MUST include UMI or read counts](https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/5.1.0/schema.md#x-matrix-layers). Author-normalized data layers [as defined in the CELLxGENE dataset schema](https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/5.1.0/schema.md#x-matrix-layers) MUST NOT be included in the Census.

#### Sample types

Only observations (cells) from primary tissue MUST be included in the Census. Thus, ONLY those observations with a [`tissue_type`](https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/5.1.0/schema.md#tissue_type) value equal to "tissue" MUST be included; other values of `tissue_type` MUST NOT be included.

#### Repeated data

When a cell is represented multiple times in CELLxGENE Discover, only one is marked as the primary cell. This is defined in the CELLxGENE dataset schema under [`is_primary_data`](https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/5.1.0/schema.md#is_primary_data). This information MUST be included in the Census cell metadata to enable queries that retrieve datasets (see cell metadata below), and all cells MUST be included in the Census.

### Data encoding and organization

The Census MUST be encoded as a `SOMACollection` which will be referenced  as `census_obj` in the following sections. `census_obj`  MUST have two keys `"census_info"` and `"census_data"` whose contents are defined in the sections below.

#### Census information `census_obj["census_info"]` - `SOMACollection`

A series of summary and metadata tables MUST be included in this `SOMACollection`:

##### Census metadata – `census_obj​​["census_info"]["summary"]` – `SOMADataFrame`

Census metadata MUST be stored as a `SOMADataFrame` with two columns:

<table>
<thead>
  <tr>
    <th>Column</th>
    <th>Encoding</th>
    <th>Description</th>
  </tr>
</thead>
<tbody>
  <tr>
    <td>label</td>
    <td>string</td>
    <td>Human readable label of metadata variable</td>
  </tr>
  <tr>
    <td>value </td>
    <td>string</td>
    <td>Value associated to metadata variable</td>
  </tr>
</tbody>
</table>

This `SOMADataFrame` MUST have the following rows:

1. Census schema version:
   1. label: `"census_schema_version"`
   2. value: Semver schema version.
2. Census build date:
   1. label: `"census_build_date"`
   2. value: The date this Census was built in ISO 8601 date format
3. Dataset schema version:
   1. label: `"dataset_schema_version"`
   2. value: The CELLxGENE Discover schema version of the source H5AD files.
4. Total number of cells included in this Census build:
   1. label: `"total_cell_count"`
   2. value: Cell count
5. Unique number of cells included in this Census build (is_primary_data == True)
   1. label: `"unique_cell_count"`
   2. value: Cell count
6. Number of human donors included in this Census build. Donors are guaranteed to be unique within datasets, not across all Census.
   1. label: `"number_donors_homo_sapiens"`
   2. value: Donor count
7. Number of mouse donors included in this Census build. Donors are guaranteed to be unique within datasets, not across all Census.
   1. label: `"number_donors_mus_musculus"`
   2. value: Donor count

An example of this `SOMADataFrame` is shown below:

<table>
<thead>
  <tr>
    <th>label</th>
    <th>value</th>
  </tr>
</thead>
<tbody>
  <tr>
    <td>census_schema_version</td>
    <td>2.0.0</td>
  </tr>
  <tr>
    <td>census_build_date</td>
    <td>2024-04-01</td>
  </tr>
  <tr>
    <td>dataset_schema_version </td>
    <td>5.1.0</td>
  </tr>
  <tr>
    <td>total_cell_count</td>
    <td>10000</td>
  </tr>
  <tr>
    <td>unique_cell_count</td>
    <td>1000</td>
  </tr>
  <tr>
    <td>number_donors_homo_sapiens</td>
    <td>100</td>
  </tr>
  <tr>
    <td>number_donors_mus_musculus</td>
    <td>100</td>
  </tr>
</tbody>
</table>

#### Census table of CELLxGENE Discover datasets – `census_obj["census_info"]["datasets"]` – `SOMADataFrame`

All datasets used to build the Census MUST be included in a table modeled as a `SOMADataFrame`. Each row MUST correspond to an individual dataset with the following columns:

<table>
<thead>
  <tr>
    <th>Column</th>
    <th>Encoding</th>
    <th>Description</th>
  </tr>
</thead>
<tbody>
  <tr>
    <td>citation</td>
    <td>string</td>
    <td>As defined in the CELLxGENE schema.</td>
  </tr>
  <tr>
    <td>collection_id</td>
    <td>string</td>
    <td rowspan="6">As defined in CELLxGENE Discover <a href="https://api.cellxgene.cziscience.com/curation/ui/">data schema</a> (see &quot;Schemas&quot; section for field definitions)".</td>
  </tr>
  <tr>
    <td>collection_name</td>
    <td>string</td>
  </tr>
  <tr>
    <td>collection_doi</td>
    <td>string</td>
  </tr>
  <tr>
    <td>collection_doi_label</td>
    <td>string</td>
  </tr>
  <tr>
    <td>dataset_id</td>
    <td>string</td>
  </tr>
  <tr>
    <td>dataset_title</td>
    <td>string</td>
  </tr>
  <tr>
    <td>dataset_h5ad_path</td>
    <td>string</td>
    <td>Relative path to the source h5ad file in the Census storage bucket.</td>
  </tr>
  <tr>
    <td>dataset_total_cell_count</td>
    <td>int</td>
    <td>Total number of cells from the dataset included in the Census.</td>
  </tr>
  <tr>
    <td>dataset_version_id</td>
    <td>string</td>
    <td>As defined in CELLxGENE Discover <a href="https://api.cellxgene.cziscience.com/curation/ui/">data schema</a> (see &quot;Schemas&quot; section for field definitions)".</td>
  </tr>
</tbody>
</table>

#### Census summary cell counts  – `census_obj["census_info"]["summary_cell_counts"]` – `SOMADataframe`

Summary cell counts grouped by organism and relevant cell metadata MUST be modeled as a `SOMADataFrame` in `census_obj["census_info"]["summary_cell_counts"]`. Each row of MUST correspond to a combination of organism and metadata variables with the following columns:

<table>
<thead>
  <tr>
    <th>Column</th>
    <th>Encoding</th>
    <th>Description</th>
  </tr>
</thead>
<tbody>
  <tr>
    <td>organism</td>
    <td>string</td>
    <td>Organism label as shown in NCBITaxon  <code>"Homo sapiens"</code> or <code>"Mus musculus"</code></td>
  </tr>
  <tr>
    <td>category</td>
    <td>string</td>
    <td>Cell metadata used for grouping, one of the following:
        <ul>
          <li><code>all</code> (special case, no grouping)</li>
          <li><code>cell_type</code></li>
          <li><code>assay</code></li>
          <li><code>tissue</code></li>
          <li><code>tissue_general</code> (high-level mapping of a tissue)</li>
          <li><code>disease</code></li>
          <li><code>self_reported_ethnicity</code></li>
          <li><code>sex</code></li>
          <li><code>suspension_type</code></li>
        </ul>
  </tr>
  <tr>
    <td>label</td>
    <td>string</td>
    <td>Label associated to instance of metadata (e.g. <code>"lung"</code> if <code>category</code> is <code>"tissue"</code>). <code>"na"</code> if none.</td>
  </tr>
  <tr>
    <td>ontology_term_id</td>
    <td>string</td>
    <td>ID associated to instance of metadata (e.g. <code>"UBERON:0002048"</code> if category is <code>"tissue"</code>). <code>"na"</code> if none.</td>
  </tr>
  <tr>
    <td>total_cell_count</td>
    <td>int</td>
    <td>Total number of cell counts for the combination of values of all other fields above.</td>
  </tr>
  <tr>
    <td>unique_cell_count</td>
    <td>int</td>
    <td>Unique number of cells for the combination of values of all other fields above. Unique number of cells refers to the cell count, for this group, when <code><a href="https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/5.1.0/schema.md#is_primary_data">is_primary_data == True</a></code> </td>
  </tr>
</tbody>
</table>

Example of this `SOMADataFrame`:

<table>
<thead>
  <tr>
    <th>organism</th>
    <th>category</th>
    <th>label</th>
    <th>ontology_term_id</th>
    <th>total_cell_count</th>
    <th>unique_cell_count</th>
  </tr>
</thead>
<tbody>
  <tr>
    <td>[Homo sapiens|Mus musculus]</td>
    <td>all</td>
    <td>na</td>
    <td>na</td>
    <td>x</td>
    <td>x</td>
  </tr>
  <tr>
    <td>[Homo sapiens|Mus musculus]</td>
    <td>cell_type</td>
    <td>cell_type_a</td>
    <td>CL:XXXXX</td>
    <td>x</td>
    <td>x</td>
  </tr>
  <tr>
    <td>[Homo sapiens|Mus musculus]</td>
    <td>cell_type</td>
    <td>…</td>
    <td>…</td>
    <td>x</td>
    <td>x</td>
  </tr>
  <tr>
    <td>[Homo sapiens|Mus musculus]</td>
    <td>cell_type</td>
    <td>cell_type_N</td>
    <td>CL:XXXXX</td>
    <td>x</td>
    <td>x</td>
  </tr>
  <tr>
    <td>[Homo sapiens|Mus musculus]</td>
    <td>assay</td>
    <td>assay_a</td>
    <td>EFO:XXXXX</td>
    <td>x</td>
    <td>x</td>
  </tr>
  <tr>
    <td>[Homo sapiens|Mus musculus]</td>
    <td>assay</td>
    <td>…</td>
    <td>…</td>
    <td>x</td>
    <td>x</td>
  </tr>
  <tr>
    <td>[Homo sapiens|Mus musculus]</td>
    <td>assay</td>
    <td>assay_N</td>
    <td>EFO:XXXXX</td>
    <td>x</td>
    <td>x</td>
  </tr>
  <tr>
    <td>[Homo sapiens|Mus musculus]</td>
    <td>tissue</td>
    <td>tissue_a</td>
    <td>UBERON:XXXXX</td>
    <td>x</td>
    <td>x</td>
  </tr>
  <tr>
    <td>[Homo sapiens|Mus musculus]</td>
    <td>tissue</td>
    <td>…</td>
    <td>…</td>
    <td>x</td>
    <td>x</td>
  </tr>
  <tr>
    <td>[Homo sapiens|Mus musculus]</td>
    <td>tissue</td>
    <td>tissue_N</td>
    <td>UBERON:XXXXX</td>
    <td>x</td>
    <td>x</td>
  </tr>
  <tr>
    <td>[Homo sapiens|Mus musculus]</td>
    <td>tissue_general</td>
    <td>tissue_general_a</td>
    <td>UBERON:XXXXX</td>
    <td>x</td>
    <td>x</td>
  </tr>
  <tr>
    <td>[Homo sapiens|Mus musculus]</td>
    <td>tissue_general</td>
    <td>…</td>
    <td>…</td>
    <td>x</td>
    <td>x</td>
  </tr>
  <tr>
    <td>[Homo sapiens|Mus musculus]</td>
    <td>tissue_general</td>
    <td>tissue_general_N</td>
    <td>UBERON:XXXXX</td>
    <td>x</td>
    <td>x</td>
  </tr>
  <tr>
    <td>[Homo sapiens|Mus musculus]</td>
    <td>disease</td>
    <td>disease_a</td>
    <td>MONDO:XXXXX</td>
    <td>x</td>
    <td>x</td>
  </tr>
  <tr>
    <td>[Homo sapiens|Mus musculus]</td>
    <td>disease</td>
    <td>…</td>
    <td>…</td>
    <td>x</td>
    <td>x</td>
  </tr>
  <tr>
    <td>[Homo sapiens|Mus musculus]</td>
    <td>disease</td>
    <td>disease_N</td>
    <td>MONDO:XXXXX</td>
    <td>x</td>
    <td>x</td>
  </tr>
  <tr>
    <td>[Homo sapiens|Mus musculus]</td>
    <td>self_reported_ethnicity</td>
    <td>self_reported_ethnicity_a</td>
    <td>HANCESTRO:XXXXX</td>
    <td>x</td>
    <td>x</td>
  </tr>
  <tr>
    <td>[Homo sapiens|Mus musculus]</td>
    <td>self_reported_ethnicity</td>
    <td>…</td>
    <td>…</td>
    <td>x</td>
    <td>x</td>
  </tr>
  <tr>
    <td>[Homo sapiens|Mus musculus]</td>
    <td>self_reported_ethnicity</td>
    <td>self_reported_ethnicity_N</td>
    <td>HANCESTRO:XXXXX</td>
    <td>x</td>
    <td>x</td>
  </tr>
  <tr>
    <td>[Homo sapiens|Mus musculus]</td>
    <td>sex</td>
    <td>sex_a</td>
    <td>PATO:XXXXX</td>
    <td>x</td>
    <td>x</td>
  </tr>
  <tr>
    <td>[Homo sapiens|Mus musculus]</td>
    <td>sex</td>
    <td>…</td>
    <td>…</td>
    <td>x</td>
    <td>x</td>
  </tr>
  <tr>
    <td>[Homo sapiens|Mus musculus]</td>
    <td>sex</td>
    <td>sex_N</td>
    <td>PATO:XXXXX</td>
    <td>x</td>
    <td>x</td>
  </tr>
  <tr>
    <td>[Homo sapiens|Mus musculus]</td>
    <td>suspension_type</td>
    <td>suspension_type_a</td>
    <td>na</td>
    <td>x</td>
    <td>x</td>
  </tr>
  <tr>
    <td>[Homo sapiens|Mus musculus]</td>
    <td>suspension_type</td>
    <td>…</td>
    <td>…</td>
    <td>x</td>
    <td>x</td>
  </tr>
  <tr>
    <td>[Homo sapiens|Mus musculus]</td>
    <td>suspension_type</td>
    <td>suspension_type_N</td>
    <td>na</td>
    <td>x</td>
    <td>x</td>
  </tr>
</tbody>
</table>

#### Census table of organisms  – `census_obj["census_info"]["organisms"]` – `SOMADataframe`

Information about organisms whose cells are included in the Census MUST be included in a table modeled as a `SOMADataFrame`. Each row MUST correspond to an individual organism with the following columns:

<table>
<thead>
  <tr>
    <th>Column</th>
    <th>Encoding</th>
    <th>Description</th>
  </tr>
</thead>
<tbody>
  <tr>
    <td>organism_ontology_term_id</td>
    <td>string</td>
    <td>As defined in the CELLxGENE dataset schema.</td>
  </tr>
  <tr>
    <td>organism_label</td>
    <td>string</td>
    <td>Human-readable label as given by the ontology.</td>
  </tr>
  <tr>
    <td>organism</td>
    <td>string</td>
    <td>Machine-friendly label used to name the SOMA Experiments, see below  <a href="#census-data--census_objcensus_dataorganism--somaexperiment">Census Data section.</a></td>
  </tr>
</tbody>
</table>

An example of this `SOMADataFrame` is shown below:

<table>
<thead>
  <tr>
    <th>organism_ontology_term_id</th>
    <th>organism_label</th>
    <th>organism</th>
  </tr>
</thead>
<tbody>
  <tr>
    <td>NCBITaxon:9606</td>
    <td>Homo sapiens</td>
    <td>homo_sapiens</td>
  </tr>
  <tr>
    <td>NCBITaxon:10090</td>
    <td>Mus musculus</td>
    <td>mus_musculus</td>
  </tr>
</tbody>
</table>

### Census Data – `census_obj["census_data"][organism]` – `SOMAExperiment`

Data for *Homo sapiens* MUST be stored as a `SOMAExperiment` in `census_obj["homo_sapiens"]`.

Data for *Mus musculus* MUST be stored as a `SOMAExperiment` in `census_obj["mus_musculus"]`.

For each organism the `SOMAExperiment` MUST contain the following:

* Cell metadata – `census_obj["census_data"][organism].obs` – `SOMADataFrame`
* Data  –  `census_obj["census_data"][organism].ms` – `SOMACollection`. This `SOMACollection` MUST only contain one `SOMAMeasurement` in `census_obj["census_data"][organism].ms["RNA"]` with the following:
  * Matrix  data –  `census_obj["census_data"][organism].ms["RNA"].X` – `SOMACollection`. It MUST contain exactly one layer:
    * Count matrix – `census_obj["census_data"][organism].ms["RNA"].X["raw"]` – `SOMASparseNDArray`
    * Normalized count matrix – `census_obj["census_data"][organism].ms["RNA"].X["normalized"]` – `SOMASparseNDArray`
  * Feature metadata – `census_obj["census_data"][organism].ms["RNA"].var` – `SOMAIndexedDataFrame`
  * Feature dataset presence matrix – `census_obj["census_data"][organism].ms["RNA"]["feature_dataset_presence_matrix"]` – `SOMASparseNDArray`

#### Matrix Data, count (raw) matrix – `census_obj["census_data"][organism].ms["RNA"].X["raw"]` – `SOMASparseNDArray`

Per the CELLxGENE dataset schema, [all RNA assays MUST include UMI or read counts](https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/5.1.0/schema.md#x-matrix-layers). These counts MUST be encoded as `float32` in this `SOMASparseNDArray` with a fill value of zero (0), and no explicitly stored zero values.

#### Matrix Data, normalized count matrix – `census_obj["census_data"][organism].ms["RNA"].X["normalized"]` – `SOMASparseNDArray`

This is an experimental data artifact - it may be removed at any time.

A library-sized normalized layer, containing a normalized variant of the count (raw) matrix.
For [full-gene sequencing assays](#full-gene-sequencing-assays), given a value `X[i,j]` in the counts (raw) matrix, library-size normalized values are defined
as `normalized[i,j] = (X[i,j] / var[j].feature_length) / sum(X[i, ] / var.feature_length[j])`.
For all other assays, for a value `X[i,j]` in the counts (raw) matrix, library-size normalized values are defined
as `normalized[i,j] = X[i,j] / sum(X[i, ])`.

#### Feature metadata – `census_obj["census_data"][organism].ms["RNA"].var` – `SOMADataFrame`

The Census MUST only contain features with a [`feature_biotype`](https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/5.1.0/schema.md#feature_biotype) value of "gene".

The [gene references are pinned](https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/5.1.0/schema.md#required-gene-annotations) as defined in the CELLxGENE dataset schema.

The following columns MUST be included:

<table>
<thead>
  <tr>
    <th>Column</th>
    <th>Encoding</th>
    <th>Description</th>
  </tr>
</thead>
<tbody>
  <tr>
    <td>feature_id</td>
    <td>str</td>
    <td>Index of <code>adata.var</code> as defined in CELLxGENE dataset schema</td>
  </tr>
  <tr>
    <td>feature_name</td>
    <td>str</td>
    <td>As defined in CELLxGENE dataset schema</td>
  </tr>
  <tr>
    <td>feature_length</td>
    <td>int</td>
    <td>As defined in CELLxGENE dataset schema</a>.</td>
  </tr>
  <tr>
    <td>nnz</td>
    <td>int64</td>
    <td>For this feature, the number of non-zero values present in the `X['raw']` counts (raw) matrix.</td>
  </tr>
  <tr>
    <td>n_measured_obs</td>
    <td>int64</td>
    <td>For this feature, the number of observations present in the source H5AD (sum of feature presence matrix).</td>
  </tr>
</tbody>
</table>

#### Feature dataset presence matrix – `census_obj["census_data"][organism].ms["RNA"]["feature_dataset_presence_matrix"]` – `SOMASparseNDArray`

In some datasets, there are features not included in the source data. To clarify the difference between features that were not included and features that were not measured, for each `SOMAExperiment` the Census MUST include a presence matrix encoded as a `SOMASparseNDArray`.

For all features included in the Census, the dataset presence matrix MUST indicate what features are included in each dataset of the Census. This information MUST be encoded as a boolean matrix, `True` indicates the feature was included in the dataset, `False` otherwise. This is a two-dimensional matrix and it MUST be `N x M` where `N` is the number of datasets in the `SOMAExperiment` and `M` is the number of features. The matrix is indexed by the `soma_joinid` value of  `census_obj["census_info"]["datasets"]` and `census_obj["census_data"][organism].ms["RNA"].var`.

If the feature has at least one cell with a value greater than zero in the count data matrix X in the dataset of origin, the value MUST be `True`; otherwise, it MUST be `False`.

An example of this matrix is shown below:

<table>
<thead>
  <tr>
    <th></th>
    <th>Feature_1</th>
    <th>…</th>
    <th>Feature_M</th>
  </tr>
</thead>
<tbody>
  <tr>
    <td><b>Dataset_soma_joinid_1</b></td>
    <td>[<code>True</code>|<code>False</code>]</td>
    <td>…</td>
    <td>[<code>True</code>|<code>False</code>]</td>
  </tr>
  <tr>
    <td>…</td>
    <td>…</td>
    <td>…</td>
    <td>…</td>
  </tr>
  <tr>
    <td><b> Dataset_soma_joinid_N</b></td>
    <td>[<code>True</code>|<code>False</code>]</td>
    <td>…</td>
    <td>[<code>True</code>|<code>False</code>]</td>
  </tr>
</tbody>
</table>

#### Cell metadata – `census_obj["census_data"][organism].obs` – `SOMADataFrame`

Cell metadata MUST be encoded as a `SOMADataFrame` with the following columns:

<table>
<thead>
  <tr>
    <th>Column</th>
    <th>Encoding</th>
    <th>Description</th>
  </tr>
</thead>
<tbody>
  <tr>
    <td>dataset_id</td>
    <td>string</td>
    <td>CELLxGENE dataset ID</td>
  </tr>
  <tr>
    <td>tissue_general_ontology_term_id</td>
    <td>string</td>
    <td>High-level tissue UBERON ID as implemented <a href="https://github.com/chanzuckerberg/single-cell-data-portal/blob/9b94ccb0a2e0a8f6182b213aa4852c491f6f6aff/backend/wmg/data/tissue_mapper.py">here</a></td>
  </tr>
  <tr>
    <td>tissue_general</td>
    <td>string</td>
    <td>High-level tissue label as implemented <a href="https://github.com/chanzuckerberg/single-cell-data-portal/blob/9b94ccb0a2e0a8f6182b213aa4852c491f6f6aff/backend/wmg/data/tissue_mapper.py">here</a></td>
  </tr>
  <tr>
    <td>assay_ontology_term_id</td>
    <td colspan="2" rowspan="19">As defined in CELLxGENE dataset schema</td>
  </tr>
  <tr>
    <td>assay</td>
  </tr>
  <tr>
    <td>cell_type_ontology_term_id</td>
  </tr>
  <tr>
    <td>cell_type</td>
  </tr>
  <tr>
    <td>development_stage_ontology_term_id</td>
  </tr>
  <tr>
    <td>development_stage</td>
  </tr>
  <tr>
    <td>disease_ontology_term_id</td>
  </tr>
  <tr>
    <td>disease</td>
  </tr>
  <tr>
    <td>donor_id</td>
  </tr>
  <tr>
    <td>is_primary_data</td>
  </tr>
  <tr>
    <td>observation_joinid</td>
  </tr>
  <tr>
    <td>self_reported_ethnicity_ontology_term_id</td>
  </tr>
  <tr>
    <td>self_reported_ethnicity</td>
  </tr>
  <tr>
    <td>sex_ontology_term_id</td>
  </tr>
  <tr>
    <td>sex</td>
  </tr>
  <tr>
    <td>suspension_type</td>
  </tr>
  <tr>
    <td>tissue_ontology_term_id</td>
  </tr>
  <tr>
    <td>tissue</td>
  </tr>
  <tr>
    <td>tissue_type</td>
  </tr>
  <tr>
    <td>nnz</td>
    <td>int64</td>
    <td>For this observation, the number of non-zero measurements in the `X['raw']` counts (raw) matrix.</td>
  </tr>
  <tr>
    <td>n_measured_vars</td>
    <td>int64</td>
    <td>For this observation, the number of features present in the source H5AD (sum of feature presence matrix).</td>
  </tr>
  <tr>
    <td>raw_sum</td>
    <td>float32</td>
    <td>For this observation, the sum of the `X['raw']` counts (raw) matrix values.</td>
  </tr>
  <tr>
    <td>raw_mean_nnz</td>
    <td>float32</td>
    <td>For this observation, the mean of the `X['raw']` counts (raw) matrix values. Zeroes are excluded from the calculation.</td>
  </tr>
  <tr>
    <td>raw_variance_nnz</td>
    <td>float32</td>
    <td>For this observation, the variance of the `X['raw']` counts (raw) matrix values. Zeroes are excluded from the calculation.</td>
  </tr>
</tbody>
</table>

## Changelog

### Version 2.1.0

* Update to require [CELLxGENE schema version 5.1.0](https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/5.1.0/schema.md)
* Adds `collection_doi_label` to "Census table of CELLxGENE Discover datasets – `census_obj["census_info"]["datasets"]`"

### Version 2.0.1

* Update accepted assays for Census based on guidance from curators.

### Version 2.0.0

* Update to require [CELLxGENE schema version 5.0.0](https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/5.0.0/schema.md)
* Expanded list of assays included in the Census.
* Expanded the list of assays defined as full-gene sequencing assays, which have special `normalized` layer handling.
* Clarified handling of datasets which are multi-species on the obs or var axis.

### Version 1.3.0

* Update to require [CELLxGENE schema version 4.0.0](https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/4.0.0/schema.md)
* Adds `citation` to "Census table of CELLxGENE Discover datasets – `census_obj["census_info"]["datasets"]`"
* Adds `observation_joinid` and `tissue_type` to `obs` dataframe

### Version 1.2.0

* Update to require [CELLxGENE schema version 3.1.0](https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/3.1.0/schema.md)

### Version 1.1.0

* Adds `dataset_version_id` to "Census table of CELLxGENE Discover datasets – `census_obj["census_info"]["datasets"]`"
* Add `X["normalized"]` layer
* Add `nnz` and `n_measured_obs` columns to `ms["RNA"].var` dataframe
* Add `nnz`, `n_measured_vars`, `raw_sum`, `raw_mean_nnz` and `raw_variance_nnz` columns to `obs` dataframe

### Version 1.0.0

* Updates text to reflect official name: CZ CELLxGENE Discover Census.
* Updates `census["census_info"]["summary"]` to reflect official name in the column `label`:
  * From `"cell_census_build_date"` to `"census_build_date"`.
  * From `"cell_census_schema_version"` to `"census_schema_version"`.
* Adds the following row to `census["census_info"]["summary"]`:
  * `"dataset_schema_version"`

### Version 0.1.1

* Adds clarifying text for "Feature Dataset Presence Matrix"

### Version 0.1.0

* The "Dataset Presence Matrix" was renamed to "Feature Dataset Presence Matrix" and moved from  `census_obj["census_data"][organism].ms["RNA"].varp["dataset_presence_matrix"]`  to `census_obj["census_data"][organism].ms["RNA"]["feature_dataset_presence_matrix"]`.
* Editorial: changes all double quotes in the schema to ASCII quotes 0x22.

### Version 0.0.1

* Initial Census schema is published.



# Section: cellxgene_census-README

# CZ CELLxGENE Discover Census

The `cellxgene_census` package provides an API to facilitate the use of the CZ CELLxGENE Discover Census. For more information about the API and the project visit the [chanzuckerberg/cellxgene-census GitHub repo](https://github.com/chanzuckerberg/cellxgene-census/).

## For More Help

For more help, please file a issue on the repo, or contact us at <soma@chanzuckerberg.com>.

If you believe you have found a security issue, we would appreciate notification. Please send email to <security@chanzuckerberg.com>.

## Development Environment Setup

- Create a virtual environment using `venv` or `conda`
- `cd` to the root of this repository
- `pip install -e api/python/cellxgene_census`
- To install dependencies needed to work on the [experimental](./src/cellxgene_census/experimental/) portion of the API:
  `pip install -e 'api/python/cellxgene_census[experimental]'`.
- `pip install jupyterlab`
- **Test it!** Either open up a new `jupyter` notebook or the `python` interpreter and run this code:

```python
import cellxgene_census

with cellxgene_census.open_soma() as census:

    cell_metadata = cellxgene_census.get_obs(
        census,
        "homo_sapiens",
        value_filter = "sex == 'female' and cell_type in ['microglial cell', 'neuron']",
        column_names = ["assay", "cell_type", "tissue", "tissue_general", "suspension_type", "disease"]
    )
    cell_metadata
```

The output is a `pandas.DataFrame` with over 600K cells meeting our query criteria and the selected columns:

```python

The "stable" release is currently 2023-12-15. Specify 'census_version="2023-12-15"' in future calls to open_soma() to ensure data consistency.

                assay        cell_type                 tissue tissue_general suspension_type disease     sex
0        Smart-seq v4  microglial cell  middle temporal gyrus          brain         nucleus  normal  female
1        Smart-seq v4  microglial cell  middle temporal gyrus          brain         nucleus  normal  female
2        Smart-seq v4  microglial cell  middle temporal gyrus          brain         nucleus  normal  female
3        Smart-seq v4  microglial cell  middle temporal gyrus          brain         nucleus  normal  female
4        Smart-seq v4  microglial cell  middle temporal gyrus          brain         nucleus  normal  female
...               ...              ...                    ...            ...             ...     ...     ...
607636  microwell-seq           neuron          adrenal gland  adrenal gland            cell  normal  female
607637  microwell-seq           neuron          adrenal gland  adrenal gland            cell  normal  female
607638  microwell-seq           neuron          adrenal gland  adrenal gland            cell  normal  female
607639  microwell-seq           neuron          adrenal gland  adrenal gland            cell  normal  female
607640  microwell-seq           neuron          adrenal gland  adrenal gland            cell  normal  female

[607641 rows x 7 columns]

```

- Learn more about the Census API by going through the tutorials in the [notebooks](../notebooks/)



# Section: articles-2023-20230808-r_api_release

# R package `cellxgene.census` V1 is out!

*Published:* *August 7th, 2023*

*By:* *[Pablo Garcia-Nieto](mailto:pgarcia-nieto@chanzuckerberg.com)*

The Census team is pleased to announce the release of the R package `cellxgene.census`. 🎉 🎉

This has been long coming since our Python release back in May. Now, from R, computational biologists can access the Census data which is the largest standardized aggregation of single-cell data, composed of >33M cells and >60K genes.

With `cellxgene.census` in a few seconds users can access any slice of Census data using cell or gene filters across hundreds of datasets. The data can be fetched in an iterative fashion for bigger-than-memory slices of data, or quickly exported to basic R structures, and [Seurat](https://satijalab.org/seurat/) or [SingleCellExperiment](https://bioconductor.org/packages/release/bioc/html/SingleCellExperiment.html) for downstream analysis.

![image](20230808-r_api_release.svg)

## Installation and usage

Users can install `cellxgene.census` and its dependencies following the [installation instructions](../../cellxgene_census_docsite_installation.md).

To learn more about the package please make sure to check out the following resources:

* [Quick start guide.](../../cellxgene_census_docsite_quick_start.md)
* [R reference docs and tutorials.](https://chanzuckerberg.github.io/cellxgene-census/r/index.html)
* [Querying and slicing data tutorial.](https://chanzuckerberg.github.io/cellxgene-census/r/articles/census_query_extract.html)

## Census R package is made possible by `tiledbsoma`

The `cellxgene.census` package relies on [TileDB-SOMA](https://github.com/single-cell-data/TileDB-SOMA) R's package `tiledbsoma` for all of its data access capabilities as shown in the next section.

CZI and TileDB have worked closely on the development of `tiledbsoma` and recently upgraded it from beta to its first stable version. Release notes can be found [here](https://github.com/single-cell-data/TileDB-SOMA/releases/tag/1.4.0).

## Efficient access to single-cell data for >33M cells from R

Census hosts ever-growing [data releases](../../cellxgene_census_docsite_data_release_info.md) from CZ CELLxGENE Discover, representing the largest aggregation of standardized single-cell data.

Census data are accompanied by cell and gene metadata that have been standardized on ontologies across all datasets hosted in CZ CELLxGENE Discover. For example all cell types and tissues have been mapped to a value of the CL and UBERON ontologies, respectively. You can find more about the data in the [Census data and schema](../../cellxgene_census_docsite_schema.md) page.

With the `cellxgene.census` R package, researchers can have access to all of these data and metadata directly from an R session with the following capabilities:

### Easy-to-use handles to the cloud-hosted Census data

From R users can get a handle to the data by opening the Census.

```r
library("cellxgene.census")

census <- open_soma()

# Your work!

census$close()
```

### Querying and reading single-cell metadata from Census

Following our [Census data and schema](../../cellxgene_census_docsite_schema.md), users can navigate and query Census data and metadata by using any combination of gene and cell filters.

For example, reading a slice of the human cell metadata for more than 300K cells with Microglial cells or Neurons from female donors:

```r
library("cellxgene.census")

census <- open_soma()

# Open obs SOMADataFrame
cell_metadata <- census$get("census_data")$get("homo_sapiens")$get("obs")

# Read as an iterator of Arrow Tables
cell_metadata <- cell_metadata$read(
   value_filter = "sex == 'female' & cell_type %in% c('microglial cell', 'neuron')",
   column_names = c("assay", "cell_type", "sex", "tissue", "tissue_general", "suspension_type", "disease")
)

# Retrieve all metadata at once
cell_metadata <- cell_metadata$concat()

# Convert to R tibble (dataframe)
cell_metadata <- as.data.frame(cell_metadata)

census$close()
```

### Exporting Census slices to `Seurat` and `SingleCellExperiment`

Similarly, users can query both the single-cell data along with its metadata and export them to  `Seurat` or `SingleCellExperiment` objects for downstream analysis:

```r
library("cellxgene.census")

census <- open_soma()

organism <- "Homo sapiens"
gene_filter <- "feature_id %in% c('ENSG00000107317', 'ENSG00000106034')"
cell_filter <-  "cell_type == 'leptomeningeal cell'"
cell_columns <- c("assay", "cell_type", "tissue", "tissue_general", "suspension_type", "disease")

# Get Seurat object
seurat_obj <- get_seurat(
  census = census,
  organism = organism,
  var_value_filter = gene_filter,
  obs_value_filter = cell_filter,
  obs_column_names = cell_columns
)

# Get SingleCellExperiment object
sce_obj <- get_single_cell_experiment(
  census = census,
  organism = organism,
  var_value_filter = gene_filter,
  obs_value_filter = cell_filter,
  obs_column_names = cell_columns
)

census$close()
```

### Streaming data incrementally in chunks

Sometimes Census queries can be too large to be loaded in memory. TileDB-SOMA allows users to query Census data in an incremental fashion using iterators.

To find out more about iterable-based queries you can check out the following resources:

* [Memory-efficient queries from R.](../../cellxgene_census_docsite_quick_start.md#id2)
* [The SOMA objects overview from TileDB-SOMA.](https://single-cell-data.github.io/TileDB-SOMA/articles/soma-objects.html)



# Section: cellxgene_census_docsite_quick_start

# Quick start

This page provides details to start using the Census. Click [here](examples.rst) for more detailed Python tutorials (R vignettes coming soon).

**Contents:**

1. [Installation](#installation).
2. [Python quick start](python-quick-start).
3. [R quick start](r-quick-start).

## Installation

Install the Census API by following [these instructions.](cellxgene_census_docsite_installation.md)

## Python quick start

Below are 3 examples of common operations you can do with the Census. As a reminder, the reference documentation for the API can be accessed via `help()`:

```python
import cellxgene_census

help(cellxgene_census)
help(cellxgene_census.get_anndata)
# etc
```

### Querying a slice of cell metadata

The following reads the cell metadata and filters `female` cells of cell type `microglial cell` or `neuron`, and selects the columns `assay`, `cell_type`, `tissue`, `tissue_general`, `suspension_type`, and `disease`.

```python
import cellxgene_census

with cellxgene_census.open_soma() as census:

    # Reads SOMADataFrame as a slice
    cell_metadata = census["census_data"]["homo_sapiens"].obs.read(
        value_filter = "sex == 'female' and cell_type in ['microglial cell', 'neuron']",
        column_names = ["assay", "cell_type", "tissue", "tissue_general", "suspension_type", "disease"]
    )

    # Concatenates results to pyarrow.Table
    cell_metadata = cell_metadata.concat()

    # Converts to pandas.DataFrame
    cell_metadata = cell_metadata.to_pandas()

    print(cell_metadata)
```

The output is a `pandas.DataFrame` with over 300K cells meeting our query criteria and the selected columns.

```bash
The "stable" release is currently 2023-07-25. Specify 'census_version="2023-07-25"' in future calls to open_soma() to ensure data consistency.
                assay        cell_type         tissue tissue_general suspension_type disease     sex
0           10x 3' v3  microglial cell            eye            eye            cell  normal  female
1           10x 3' v3  microglial cell            eye            eye            cell  normal  female
2           10x 3' v3  microglial cell            eye            eye            cell  normal  female
3           10x 3' v3  microglial cell            eye            eye            cell  normal  female
4           10x 3' v3  microglial cell            eye            eye            cell  normal  female
...               ...              ...            ...            ...             ...     ...     ...
379219  microwell-seq           neuron  adrenal gland  adrenal gland            cell  normal  female
379220  microwell-seq           neuron  adrenal gland  adrenal gland            cell  normal  female
379221  microwell-seq           neuron  adrenal gland  adrenal gland            cell  normal  female
379222  microwell-seq           neuron  adrenal gland  adrenal gland            cell  normal  female
379223  microwell-seq           neuron  adrenal gland  adrenal gland            cell  normal  female

[379224 rows x 7 columns]
```

### Obtaining a slice as AnnData

The following creates an `anndata.AnnData` object on-demand with the same cell filtering criteria as above and filtering only the genes `ENSG00000161798`, `ENSG00000188229`.

```python
import cellxgene_census

with cellxgene_census.open_soma() as census:
    adata = cellxgene_census.get_anndata(
        census = census,
        organism = "Homo sapiens",
        var_value_filter = "feature_id in ['ENSG00000161798', 'ENSG00000188229']",
        obs_value_filter = "sex == 'female' and cell_type in ['microglial cell', 'neuron']",
        column_names = {"obs": ["assay", "cell_type", "tissue", "tissue_general", "suspension_type", "disease"]},
    )

    print(adata)
```

The output with about 300K cells and 2 genes can be now used for downstream analysis using [scanpy](https://scanpy.readthedocs.io/en/stable/).

``` bash
AnnData object with n_obs × n_vars = 379224 × 2
    obs: 'assay', 'cell_type', 'tissue', 'tissue_general', 'suspension_type', 'disease', 'sex'
    var: 'soma_joinid', 'feature_id', 'feature_name', 'feature_length'
```

### Memory-efficient queries

This example provides a demonstration to access the data for larger-than-memory operations using **TileDB-SOMA** operations.

First we initiate a lazy-evaluation query to access all brain and male cells from human. This query needs to be closed — `query.close()` — or called in a context manager — `with ...`.

```python
import cellxgene_census
import tiledbsoma

with cellxgene_census.open_soma() as census:

    human = census["census_data"]["homo_sapiens"]
    query = human.axis_query(
       measurement_name = "RNA",
       obs_query = tiledbsoma.AxisQuery(
           value_filter = "tissue == 'brain' and sex == 'male'"
       )
    )

    # Continued below

```

Now we can iterate over the matrix count, as well as the cell and gene metadata. For example, to iterate over the matrix count, we can get an iterator and perform operations for each iteration.

```python
    # Continued from above

    iterator = query.X("raw").tables()

    # Get an iterative slice as pyarrow.Table
    raw_slice = next (iterator)
    ...
```

And you can now perform operations on each iteration slice. As with any any Python iterator this logic can be wrapped around a `for` loop.

And you must close the query.

```python
    # Continued from above
    query.close()
```

## R quick start

Below are 3 examples of common operations you can do with the Census. As a reminder, the reference documentation for the API can be accessed via `?`:

```r
library("cellxgene.census")

?cellxgene.census::get_seurat
```

### Querying a slice of cell metadata

The following reads the cell metadata and filters `female` cells of cell type `microglial cell` or `neuron`, and selects the columns `assay`, `cell_type`, `tissue`, `tissue_general`, `suspension_type`, and `disease`.

The `cellxgene.census` package uses [R6](https://r6.r-lib.org/articles/Introduction.html) classes and we recommend you to get familiar with their usage.

```r
library("cellxgene.census")

census <- open_soma()

# Open obs SOMADataFrame
cell_metadata <-  census$get("census_data")$get("homo_sapiens")$get("obs")

# Read as Arrow Table
cell_metadata <-  cell_metadata$read(
   value_filter = "sex == 'female' & cell_type %in% c('microglial cell', 'neuron')",
   column_names = c("assay", "cell_type", "sex", "tissue", "tissue_general", "suspension_type", "disease")
)

# Concatenates results to an Arrow Table
cell_metadata <-  cell_metadata$concat()

# Convert to R tibble (dataframe)
cell_metadata <-  as.data.frame(cell_metadata)

print(cell_metadata)

census$close()
```

The output is a `tibble` with over 300K cells meeting our query criteria and the selected columns.

```bash
# A tibble: 379,224 × 7
   assay     cell_type       sex   tissue tissue_general suspension_type disease
   <chr>     <chr>           <chr> <chr>  <chr>          <chr>           <chr>
 1 10x 3' v3 microglial cell fema… eye    eye            cell            normal
 2 10x 3' v3 microglial cell fema… eye    eye            cell            normal
 3 10x 3' v3 microglial cell fema… eye    eye            cell            normal
 4 10x 3' v3 microglial cell fema… eye    eye            cell            normal
 5 10x 3' v3 microglial cell fema… eye    eye            cell            normal
 6 10x 3' v3 microglial cell fema… eye    eye            cell            normal
 7 10x 3' v3 microglial cell fema… eye    eye            cell            normal
 8 10x 3' v3 microglial cell fema… eye    eye            cell            normal
 9 10x 3' v3 microglial cell fema… eye    eye            cell            normal
10 10x 3' v3 microglial cell fema… eye    eye            cell            normal
# ℹ 379,214 more rows
# ℹ Use `print(n = ...)` to see more rows
```

### Obtaining a slice as a `Seurat` or `SingleCellExperiment` object

The following creates a Seurat object on-demand with a smaller set of cells and filtering only the genes `ENSG00000161798`, `ENSG00000188229`.

```r
library("cellxgene.census")
library("Seurat")

census <-  open_soma()

organism <-  "Homo sapiens"
gene_filter <-  "feature_id %in% c('ENSG00000107317', 'ENSG00000106034')"
cell_filter <-   "cell_type == 'sympathetic neuron'"
cell_columns <-  c("assay", "cell_type", "tissue", "tissue_general", "suspension_type", "disease")

seurat_obj <-  get_seurat(
   census = census,
   organism = organism,
   var_value_filter = gene_filter,
   obs_value_filter = cell_filter,
   obs_column_names = cell_columns
)

print(seurat_obj)
```

The output with over 4K cells and 2 genes can be now used for downstream analysis using [Seurat](https://satijalab.org/seurat/).

```shell
An object of class Seurat
2 features across 4744 samples within 1 assay
Active assay: RNA (2 features, 0 variable features)
```

Similarly a `SingleCellExperiment` object can be created.

```r
library("SingleCellExperiment")

sce_obj <-  get_single_cell_experiment(
   census = census,
   organism = organism,
   var_value_filter = gene_filter,
   obs_value_filter = cell_filter,
   obs_column_names = cell_columns
)

print(sce_obj)
```

The output with over 4K cells and 2 genes can be now used for downstream analysis using the [Bioconductor ecosystem](https://bioconductor.org/packages/release/bioc/html/SingleCellExperiment.html).

```shell
class: SingleCellExperiment
dim: 2 4744
metadata(0):
assays(1): counts
rownames(2): ENSG00000106034 ENSG00000107317
rowData names(2): feature_name feature_length
colnames(4744): obs48350835 obs48351829 ... obs52469564 obs52470190
colData names(6): assay cell_type ... suspension_type disease
reducedDimNames(0):
mainExpName: RNA
altExpNames(0):
```

### Memory-efficient queries

This example provides a demonstration to access the data for larger-than-memory operations using **TileDB-SOMA** operations.

First we initiate a lazy-evaluation query to access all brain and male cells from human. This query needs to be closed — `query$close()`.

```r
library("cellxgene.census")
library("tiledbsoma")

human <-  census$get("census_data")$get("homo_sapiens")
query <-  human$axis_query(
  measurement_name = "RNA",
  obs_query = SOMAAxisQuery$new(
    value_filter = "tissue == 'brain' & sex == 'male'"
  )
)

# Continued below

```

Now we can iterate over the matrix count, as well as the cell and gene metadata. For example, to iterate over the matrix count, we can get an iterator and perform operations for each iteration.

```r
# Continued from above

iterator <-  query$X("raw")$tables()
# For sparse matrices use query$X("raw")$sparse_matrix()

# Get an iterative slice as an Arrow Table
raw_slice <-  iterator$read_next()

#...
```

And you can now perform operations on each iteration slice. This logic can be wrapped around a `while()` loop and checking the iteration state by monitoring the logical output of `iterator$read_complete()`

And you must close the query and census.

```r
# Continued from above
query.close()
census.close()
```



# Section: articles-2024-20240710_embedding_metrics_dec_2023_lts

# Benchmarks of single-cell Census models

*Published:* *July 11th, 2024*

*Updated:* *July 19th, 2024*. Model names fixed.

*By:* *[Emanuele Bezzi](mailto:ebezzi@chanzuckerberg.com), [Pablo Garcia-Nieto](mailto:pgarcia-nieto@chanzuckerberg.com)*

In 2023, the Census team released a series of cells embeddings (available at the [Census Model page](https://cellxgene.cziscience.com/census-models)) compatible with the [Census LTS version `census_version="2023-12-15"`](https://chanzuckerberg.github.io/cellxgene-census/cellxgene_census_docsite_data_release_info.html#lts-2023-12-15), so that users can access and download for any slice of Census data.

These embeddings were generated via different large-scale models; in this article we present the results of light benchmarking of them. We hope that these benchmarks provide an initial picture to users on, 1) the strength of biological signal captured by these embeddings and, 2) the level of batch correction they exert.

We advise our users to consider these benchmarks as first-pass information and we recommend further benchmarking for a more comprehensive view of the embeddings and for task-oriented applications.

The benchmarks were run on the following embeddings:

- scVI latent spaces from a model trained on all Census data.
- Fine-tuned Geneformer.
- scGPT.
- Universal Cell Embeddings (UCE).

For more details on each model please see the [Census Model page](https://cellxgene.cziscience.com/census-models).

## Accessing the embeddings included in the benchmark

Please the [Census Model page](https://cellxgene.cziscience.com/census-models) for full details. Shortly, you can see the embeddings available for the Census LTS version `census_version="2023-12-15"` using the Census API as follows.

```python
import cellxgene_census.experimental.get_all_available_embeddings
cellxgene_census.experimental.get_all_available_embeddings(census_version="2023-12-15")
```

With the exception of NMF factors, all other human embeddings were included in the benchmarks below. If you would want to access the embeddings for any slice of data you can utilize the parameter `obs_embeddings` from  the`get_anndata()` method of the Census API, for example:

```python
import cellxgene_census
census = cellxgene_census.open_soma(census_version="2023-12-15")
adata = cellxgene_census.get_anndata(
    census,
    organism = "homo_sapiens",
    measurement_name = "RNA",
    obs_value_filter = "tissue_general == 'central nervous system'",
    obs_embeddings = ["scvi"]
)
```

## Benchmarks of Census Embeddings

### About the benchmarks

We executed a series of benchmarks falling into two general types: one to assess the level of biological signal contained in the embeddings, and the second to measure the level of correction for batch effects. In general, the utility of the embeddings increases as a function of these two set of benchmarks.

For each of the type, the benchmarks can be further subdivided by their "mode". A series of metrics assess the embedding space, and the others assess the capacity of the embeddings to predict labels.

The table below shows a breakdown of the benchmarks we used in this report.

<table class="custom-table">
  <thead>
      <tr>
        <th>Type</th>
        <th>Mode</th>
        <th>Metric</th>
        <th>Description</th>
      </tr>
  </thead>
  <tbody>
      <tr>
        <td rowspan="6">Bio-conservation</td>
        <td rowspan="3">Embedding<br>Space</td>
        <td><code>leiden_nmi</code></td>
        <td>Normalized Mutual Information of biological labels and leiden clusters. Described in <a href="https://scib-metrics.readthedocs.io/en/stable/generated/scib_metrics.nmi_ari_cluster_labels_leiden.html">Luecken et al.</a> and implemented in <a href="https://scib-metrics.readthedocs.io/en/stable/generated/scib_metrics.nmi_ari_cluster_labels_leiden.html">scib-metrics.</a></td>
      </tr>
      <tr>
        <td><code>leiden_ari</code></td>
        <td>Adjusted Rand Index of biological labels and leiden clusters. Described in <a href="https://scib-metrics.readthedocs.io/en/stable/generated/scib_metrics.nmi_ari_cluster_labels_leiden.html">Luecken et al.</a> and implemented in <a href="https://scib-metrics.readthedocs.io/en/stable/generated/scib_metrics.nmi_ari_cluster_labels_leiden.html">scib-metrics.</a></td>
      </tr>
      <tr>
        <td><code>silhouette_label</code></td>
        <td>Silhouette score with respect to biological labels. Described in <a href="https://scib-metrics.readthedocs.io/en/stable/generated/scib_metrics.nmi_ari_cluster_labels_leiden.html">Luecken et al.</a> and implemented in <a href="https://scib-metrics.readthedocs.io/en/stable/generated/scib_metrics.silhouette_label.html">scib-metrics.</a></td>
      </tr>
      <tr>
           <td rowspan="3">Label<br>Classifier</td>
        <td><code>classifier_svm</code></td>
        <td>Accuracy of biological label prediction using a SVM (60/40 train/test split). Implemented <a href="https://github.com/chanzuckerberg/cellxgene-census/blob/f44637ba33567400820407f4f7b9984e52966156/tools/models/metrics/run-scib.py#L36">here</a>.</td>
      </tr>
      <tr>
<td><code>classifier_forest</code></td>
        <td>Accuracy of biological label prediction using a Random Forest classifier (60/40 train/test split). Implemented <a href="https://github.com/chanzuckerberg/cellxgene-census/blob/f44637ba33567400820407f4f7b9984e52966156/tools/models/metrics/run-scib.py#L39">here</a>.</td>
      </tr>
      <tr>
<td><code>classifier_lr</code></td>
        <td>Accuracy of biological label prediction using a Logistic regression classifier (60/40 train/test split). Implemented <a href="https://github.com/chanzuckerberg/cellxgene-census/blob/f44637ba33567400820407f4f7b9984e52966156/tools/models/metrics/run-scib.py#L39">here</a>.</td>
      </tr>
      <tr>
        <td rowspan="5">Batch-correction</td>
        <td rowspan="2">Embedding<br>Space</td>
        <td><code>silhouette_batch</code></td>
        <td>1- silhouette score with respect to biological labels. Described in <a href="https://scib-metrics.readthedocs.io/en/stable/generated/scib_metrics.nmi_ari_cluster_labels_leiden.html">Luecken et al.</a> and implemented in <a href="https://scib-metrics.readthedocs.io/en/stable/generated/scib_metrics.nmi_ari_cluster_labels_leiden.html">scib-metrics.</a></td>
      </tr>
      <tr>
        <td><code>entropy</code></td>
        <td>Average of neighborhood entropy of batch labels per cell. Implemented <a href="https://github.com/chanzuckerberg/cellxgene-census/blob/f44637ba33567400820407f4f7b9984e52966156/tools/models/metrics/run-scib.py#L86">here</a>.</td>
      </tr>
      <tr>
           <td rowspan="3">Label<br>Classifier</td>
        <td><code>classifier_svm</code></td>
        <td>1 - accuracy of batch label prediction using a SVM (60/40 train/test split). Implemented <a href="https://github.com/chanzuckerberg/cellxgene-census/blob/f44637ba33567400820407f4f7b9984e52966156/tools/models/metrics/run-scib.py#L45">here</a>.</td>
      </tr>
      <tr>
<td><code>classifier_forest</code></td>
        <td>1 - accuracy of batch label prediction using a Random Forest classifier (60/40 train/test split). Implemented <a href="https://github.com/chanzuckerberg/cellxgene-census/blob/f44637ba33567400820407f4f7b9984e52966156/tools/models/metrics/run-scib.py#L48">here</a>.</td>
      </tr>
      <tr>
<td><code>classifier_lr</code></td>
        <td>1 - accuracy of batch label prediction using a Logistic regression classifier (60/40 train/test split). Implemented <a href="https://github.com/chanzuckerberg/cellxgene-census/blob/f44637ba33567400820407f4f7b9984e52966156/tools/models/metrics/run-scib.py#L42">here</a>.</td>
      </tr>
  </tbody>
</table>

**Table 1:** List of benchmarks.

### Benchmark results

As reminder the benchmarks were run on the following embeddings:

- scVI latent spaces from a model trained on all Census data.
- Fine-tuned Geneformer.
- scGPT.
- Universal Cell Embeddings (UCE).

#### Summary

The following are averages for all the metrics shown in the following sections.

```{figure} ./20240710_metrics_0_summary.png
:alt: Bio-conservation single-cell Census benchmark
:align: center
:figwidth: 90%

**Figure 1. Summary of all benchmarks.** Numerical averages across the metric types and modes from all bio- and batch-labels across the tissues in this report.
```

#### Bio-conservation

The bio-conservation metrics were run the in following biological labels in a Census cells from Adipose Tissue and Spinal Cord:

- Cell subclass: a higher definition of a cell type with maximum of 73 unique labels, as defined on the CELLxGENE collection page.
- Cell class: an even higher definition of a cell type with a maximum of 22 unique labels, also defined on the CELLxGENE collection page.

```{figure} ./20240710_metrics_1_bio_emb.png
:alt: Bio-conservation single-cell Census benchmark
:align: center
:figwidth: 90%

**Figure 2. Bio-conservation metrics on the embedding space.** Higher values signify better performance, max value for all metrics is 1.
```

```{figure} ./20240710_metrics_2_bio_classifier.png
:alt: Bio-conservation single-cell Census benchmark
:align: center
:figwidth: 90%

**Figure 3. Bio-conservation metrics based on label classifiers.** Values represent label prediction accuracy. Higher values signify better performance, max value for all metrics is 1.
```

#### Batch-correction

The batch-correction metrics were run the in following batch labels in a Census cells from Adipose Tissue and Spinal Cord:

- Assay: the sequencing technology.
- Dataset: the dataset from which the cell originated from.
- Suspension type: cell vs nucleus.
- Batch: the concatenation of values for all of the above.

```{figure} ./20240710_metrics_3_batch_emb.png
:alt: Batch-correction single-cell Census benchmark
:align: center
:figwidth: 90%

**Figure 4. Batch-correction metrics on the embedding space.** Higher values signify better performance, max value for `silhouette_batch` is 1, `entropy` values should only be compared within the tissue/label combination and not across. 
```

```{figure} ./20240710_metrics_4_batch_classifier.png
:alt: Batch-correction single-cell Census benchmark
:align: center
:figwidth: 90%

**Figure 5. Batch-correction metrics based on label classifiers.** Values represent **1 - label prediction accuracy**. In theory higher values signify better performance indicating that prediction of batch labels is not accurate. However foundation models may be designed to learn *all* information including technical variation, please refer to the original publications of the models to learn more about them. 
```

## Source data

All data was obtained from the Census API, to fetch the data used in this report you can execute the following in Python. To get the cell subclass and cell class please refer to the [CellxGene Ontology Guide API](https://github.com/chanzuckerberg/cellxgene-ontology-guide/tree/main).

```python
import cellxgene_census

val_filters = {
   "adipose": "tissue_general == 'adipose tissue' and is_primary_data == True",
   "spinal": "tissue_general == 'spinal cord' and is_primary_data == True",
}

embedding_names = ["geneformer", "scgpt", "scvi", "uce"]
embedding_names = ["scvi"]
column_names = {
   "obs": ["cell_type_ontology_term_id", "cell_type", "assay", "suspension_type", "dataset_id", "soma_joinid"]
}

census = cellxgene_census.open_soma(census_version="2023-12-15")

adatas = []
for tissue in val_filters:
    adatas.append(
       cellxgene_census.get_anndata(
           census,
           organism="homo_sapiens",
           measurement_name="RNA",
           obs_value_filter= val_filters[tissue],
           obs_embeddings=embedding_names,
           column_names=column_names,
        )
    )
```

### Batch label counts

The following shows the batch label counts per tissue:

#### Adipose tissue

<table class="custom-table">
  <thead>
      <tr>
        <th>Type</th>
        <th>Label</th>
        <th>Count</th>
      </tr>
  </thead>
  <tbody>
      <tr>
          <td rowspan="4">Assay</td>
          <td>10x 3' v3</td>
          <td>91947</td>
      </tr>
      <tr>
          <td>10x 5' transcription profiling</td>
          <td>2121</td>
      </tr>
      <tr>
          <td>microwell-seq</td>
          <td>5916</td>
      </tr>
      <tr>
          <td>Smart-seq2</td>
          <td>651</td>
      </tr>
      <tr>
          <td rowspan="2">Suspension Type</td>
          <td>nucleus</td>
          <td>72335</td>
      </tr>
      <tr>
          <td>cell</td>
          <td>23756</td>
      </tr>
      <tr>
          <td rowspan="4">Dataset</td>
          <td>9d8e5dca-03a3-457d-b7fb-844c75735c83</td>
          <td>72335</td>
      </tr>
      <tr>
          <td>53d208b0-2cfd-4366-9866-c3c6114081bc</td>
          <td>20263</td>
      </tr>
      <tr>
          <td>5af90777-6760-4003-9dba-8f945fec6fdf</td>
          <td>2121</td>
      </tr>
      <tr>
          <td>2adb1f8a-a6b1-4909-8ee8-484814e2d4bf</td>
          <td>1372</td>
      </tr>
   </tbody>
</table>

#### Spinal cord

<table class="custom-table">
  <thead>
      <tr>
        <th>Type</th>
        <th>Label</th>
        <th>Count</th>
      </tr>
  </thead>
  <tbody>
      <tr>
          <td rowspan="2">Assay</td>
          <td>10x 3' v3</td>
          <td>43840</td>
      </tr>
      <tr>
          <td>microwell-seq</td>
          <td>5916</td>
      </tr>
      <tr>
          <td rowspan="2">Suspension Type</td>
          <td>nucleus</td>
          <td>43840</td>
      </tr>
      <tr>
          <td>cell</td>
          <td>5916</td>
      </tr>
      <tr>
          <td rowspan="3">Dataset</td>
          <td>090da8ea-46e8-40df-bffc-1f78e1538d27</td>
          <td>24190</td>
      </tr>
      <tr>
          <td>c05e6940-729c-47bd-a2a6-6ce3730c4919</td>
          <td>19650</td>
      </tr>
      <tr>
          <td>2adb1f8a-a6b1-4909-8ee8-484814e2d4bf</td>
          <td>5916</td>
      </tr>
   </tbody>
</table>


