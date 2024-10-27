
# Section: notebooks-api_demo-census_summary_cell_counts

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




# Section: cellxgene_census-tests-experimental-pp-test_hvg

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
import pytest
import scanpy as sc
import tiledbsoma as soma

import cellxgene_census
from cellxgene_census.experimental.pp import (
    get_highly_variable_genes,
    highly_variable_genes,
)


@pytest.mark.experimental
@pytest.mark.live_corpus
@pytest.mark.parametrize(
    "experiment_name,obs_value_filter",
    [
        (
            "mus_musculus",
            'is_primary_data == True and tissue_general == "liver"',
        ),
        pytest.param(
            "mus_musculus",
            'is_primary_data == True and tissue_general == "skin of body"',
            marks=pytest.mark.expensive,
        ),
        pytest.param(
            "mus_musculus",
            'is_primary_data == True and tissue_general in ["heart", "lung"]',
            marks=pytest.mark.expensive,
        ),
        pytest.param(
            "mus_musculus",
            'is_primary_data == True and assay == "Smart-seq2"',
            marks=pytest.mark.expensive,
        ),
    ],
)
@pytest.mark.parametrize("n_top_genes", (50, 500))
@pytest.mark.parametrize(
    "batch_key",
    (
        None,
        "dataset_id",
        ["suspension_type", "assay_ontology_term_id"],
        pytest.param(
            ("suspension_type", "assay_ontology_term_id"),
            marks=pytest.mark.expensive,
        ),
    ),
)
@pytest.mark.parametrize(
    "span",
    (
        pytest.param(None, marks=pytest.mark.expensive),
        0.5,
    ),
)
@pytest.mark.parametrize(
    "version",
    (
        "latest",
        pytest.param("stable", marks=pytest.mark.expensive),
    ),
)
def test_hvg_vs_scanpy(
    n_top_genes: int,
    obs_value_filter: str,
    version: str,
    experiment_name: str,
    batch_key: str | tuple[str] | list[str] | None,
    span: float,
    small_mem_context: soma.SOMATileDBContext,
) -> None:
    """Compare results with ScanPy on a couple of simple tests."""

    kwargs: dict[str, Any] = {
        "n_top_genes": n_top_genes,
        "batch_key": batch_key,
        "flavor": "seurat_v3",
    }
    if span is not None:
        kwargs["span"] = span

    with cellxgene_census.open_soma(census_version=version, context=small_mem_context) as census:
        # Get the highly variable genes
        with census["census_data"][experiment_name].axis_query(
            measurement_name="RNA",
            obs_query=soma.AxisQuery(value_filter=obs_value_filter),
        ) as query:
            hvg = highly_variable_genes(query, **kwargs)
            adata = query.to_anndata(X_name="raw")

    if isinstance(batch_key, list) or isinstance(batch_key, tuple):
        # ScanPy only accepts a single column for a batch key, so create it and use it
        assert "the_batch_key" not in adata.obs.columns
        adata.obs["the_batch_key"] = (
            adata.obs[list(batch_key)].astype(str)[batch_key[0]].str.cat(adata.obs[list(batch_key[1:])])
        )
        kwargs["batch_key"] = "the_batch_key"

    try:
        scanpy_hvg = sc.pp.highly_variable_genes(adata, inplace=False, **kwargs)
    except (ZeroDivisionError, ValueError):
        # There are test cases where ScanPy will fail, rendering this "compare vs scanpy"
        # test moot. The known cases involve overly partitioned data that results in batches
        # with a very small number of samples (which manifest as a divide by zero error).
        # In these known cases, go ahead and perform the HVG (above), but skip the compare
        # assertions below.
        pytest.skip("ScanPy generated an error, likely due to batches with 1 sample")

    scanpy_hvg.index.name = "soma_joinid"
    scanpy_hvg.index = scanpy_hvg.index.astype(int)
    assert len(scanpy_hvg) == len(hvg)
    # Since scanpy 1.10, there is an extra column in the dataframe that our hvg implementation
    # does not support.  Drop it for comparison.
    assert set(scanpy_hvg.drop("gene_name", axis=1, errors="ignore").keys()) == set(hvg.keys())

    assert (hvg.index == scanpy_hvg.index).all()
    assert np.allclose(
        hvg.means.to_numpy(),
        scanpy_hvg.means.to_numpy(),
        atol=1e-5,
        rtol=1e-2,
        equal_nan=True,
    )
    assert np.allclose(
        hvg.variances.to_numpy(),
        scanpy_hvg.variances.to_numpy(),
        atol=1e-5,
        rtol=1e-2,
        equal_nan=True,
    )
    assert np.allclose(
        hvg.variances_norm.to_numpy(),
        scanpy_hvg.variances_norm.to_numpy(),
        atol=1e-5,
        rtol=1e-2,
        equal_nan=True,
    )

    # Online calculation of normalized variance will differ slightly from ScanPy's calculation,
    # so look for rank of HVGs to be close, but not identical.  Don't worry about the non-HVGs
    # (which will differ more as you get to the long tail).  This test just looks for the average
    # rank distance to be small.
    assert (
        (scanpy_hvg[scanpy_hvg.highly_variable].highly_variable_rank - hvg[hvg.highly_variable].highly_variable_rank)
        .abs()
        .sum()
        / n_top_genes
    ) < 0.05

    # Ranking will also have some noise, so check that ranking is close in the highly variable subset
    scanpy_rank = scanpy_hvg.highly_variable_rank.copy()
    hvg_rank = hvg.highly_variable_rank.copy()
    hvg_rank[pd.isna(hvg_rank)] = n_top_genes
    scanpy_rank[pd.isna(scanpy_rank)] = n_top_genes
    rank_diff = (hvg_rank - scanpy_rank)[hvg.highly_variable]
    # +/- 5 in ranking, choosen arbitrarily
    assert rank_diff.min() >= -5 and rank_diff.max() <= 5

    if "highly_variable_nbatches" in scanpy_hvg.keys() or "highly_variable_nbatches" in hvg.keys():
        # Also subject to noise, so look for "close" match
        nbatches_diff = hvg.highly_variable_nbatches - scanpy_hvg.highly_variable_nbatches
        assert nbatches_diff.min() >= -2 and nbatches_diff.max() <= 2

    assert (hvg.highly_variable == scanpy_hvg.highly_variable).all()


@pytest.mark.experimental
@pytest.mark.live_corpus
@pytest.mark.parametrize(
    "experiment_name,organism,obs_value_filter,batch_key,obs_coords",
    [
        (
            "mus_musculus",
            "Mus musculus",
            'tissue_general == "liver" and is_primary_data == True',
            None,
            None,
        ),
        (
            "mus_musculus",
            "Mus musculus",
            'is_primary_data == True and tissue_general == "heart"',
            "dataset_id",
            None,
        ),
        (
            "mus_musculus",
            "Mus musculus",
            'is_primary_data == True and tissue_general == "heart"',
            ["suspension_type", "assay_ontology_term_id"],
            None,
        ),
        pytest.param(
            "mus_musculus",
            "Mus musculus",
            "is_primary_data == True",
            "dataset_id",
            slice(750_000, 1_000_000),
            marks=pytest.mark.expensive,
        ),
        pytest.param(
            "homo_sapiens",
            "Homo sapiens",
            "is_primary_data == True",
            "dataset_id",
            slice(1_000_000, 4_000_000),
            marks=pytest.mark.expensive,
        ),
    ],
)
def test_get_highly_variable_genes(
    organism: str,
    experiment_name: str,
    obs_value_filter: str,
    batch_key: str,
    small_mem_context: soma.SOMATileDBContext,
    obs_coords: slice | None,
) -> None:
    with cellxgene_census.open_soma(census_version="2023-12-15", context=small_mem_context) as census:
        hvg = get_highly_variable_genes(
            census,
            organism=organism,
            obs_value_filter=obs_value_filter,
            n_top_genes=1000,
            batch_key=batch_key,
            obs_coords=obs_coords,
        )
        n_vars = census["census_data"][experiment_name].ms["RNA"].var.count

    assert isinstance(hvg, pd.DataFrame)
    assert len(hvg) == n_vars
    assert len(hvg[hvg.highly_variable]) == 1000


@pytest.mark.experimental
def test_hvg_error_cases(small_mem_context: soma.SOMATileDBContext) -> None:
    with cellxgene_census.open_soma(census_version="stable", context=small_mem_context) as census:
        with census["census_data"]["mus_musculus"].axis_query(measurement_name="RNA") as query:
            # Only flavor="seurat_v3" is supported
            with pytest.raises(ValueError):
                highly_variable_genes(query, flavor="oopsie")  # type: ignore[arg-type]


@pytest.mark.experimental
@pytest.mark.live_corpus
def test_max_loess_jitter_error(small_mem_context: soma.SOMATileDBContext) -> None:
    with cellxgene_census.open_soma(census_version="stable", context=small_mem_context) as census:
        with pytest.raises(ValueError):
            get_highly_variable_genes(
                census,
                organism="mus_musculus",
                obs_value_filter='is_primary_data == True and tissue_general == "heart"',
                batch_key="cell_type",
                max_loess_jitter=0.0,
            )


@pytest.mark.experimental
@pytest.mark.live_corpus
@pytest.mark.parametrize(
    "batch_key",
    [
        None,
        "suspension_type",
        ["assay_ontology_term_id", "suspension_type"],
        ["dataset_id", "assay_ontology_term_id", "suspension_type", "donor_id"],
    ],
)
def test_hvg_user_defined_batch_key_func(
    small_mem_context: soma.SOMATileDBContext,
    batch_key: str | list[str] | None,
) -> None:
    if batch_key is None:

        def batch_key_func(srs: pd.Series[Any]) -> str:
            raise AssertionError("should never be called without a batch key")

    else:
        if isinstance(batch_key, str):
            keys = set([batch_key])  # noqa: C405
        else:
            keys = set(batch_key)

        def batch_key_func(srs: pd.Series[Any]) -> str:
            assert set(srs.keys()) == keys
            return "batch0"

    with cellxgene_census.open_soma(census_version="latest", context=small_mem_context) as census:
        # Get the highly variable genes
        with census["census_data"]["mus_musculus"].axis_query(
            measurement_name="RNA",
            obs_query=soma.AxisQuery(coords=(slice(75000),)),
        ) as query:
            hvg = highly_variable_genes(
                query,
                batch_key=batch_key,
                batch_key_func=batch_key_func,
                n_top_genes=100,
            )

            assert len(hvg[hvg.highly_variable]) == 100



# Section: cellxgene_census-tests-test_get_helpers

import pytest
import scipy.sparse
import tiledbsoma as soma

import cellxgene_census
from cellxgene_census._experiment import _get_experiment


@pytest.mark.live_corpus
def test_get_experiment(census: soma.Collection) -> None:
    mouse_uri = census["census_data"]["mus_musculus"].uri
    human_uri = census["census_data"]["homo_sapiens"].uri

    assert _get_experiment(census, "mus musculus").uri == mouse_uri
    assert _get_experiment(census, "Mus musculus").uri == mouse_uri
    assert _get_experiment(census, "mus_musculus").uri == mouse_uri

    assert _get_experiment(census, "homo sapiens").uri == human_uri
    assert _get_experiment(census, "Homo sapiens").uri == human_uri
    assert _get_experiment(census, "homo_sapiens").uri == human_uri

    with pytest.raises(ValueError):
        _get_experiment(census, "no such critter")


@pytest.mark.live_corpus
@pytest.mark.parametrize("organism", ["homo_sapiens", "mus_musculus"])
def test_get_presence_matrix(organism: str, census: soma.Collection) -> None:
    census_datasets = census["census_info"]["datasets"].read().concat().to_pandas()

    pm = cellxgene_census.get_presence_matrix(census, organism)
    assert isinstance(pm, scipy.sparse.csr_matrix)
    assert pm.shape[0] == len(census_datasets)
    assert pm.shape[1] == len(
        census["census_data"][organism].ms["RNA"].var.read(column_names=["soma_joinid"]).concat().to_pandas()
    )

    census.close()



# Section: cellxgene_census-src-cellxgene_census-_presence_matrix

# Copyright (c) 2022, Chan Zuckerberg Initiative
#
# Licensed under the MIT License.

"""Presence matrix methods.

Methods to retrieve the feature dataset presence matrix.
"""

import tiledbsoma as soma
from scipy import sparse

from ._experiment import _get_experiment


def get_presence_matrix(
    census: soma.Collection,
    organism: str,
    measurement_name: str = "RNA",
) -> sparse.csr_matrix:
    """Read the feature dataset presence matrix and return as a :class:`scipy.sparse.csr_array`. The
    returned sparse matrix is indexed on the first dimension by the dataset ``soma_joinid`` values,
    and on the second dimension by the ``var`` :class:`pandas.DataFrame` ``soma_joinid`` values.

    Args:
        census:
            The census from which to read the presence matrix.
        organism:
            The organism to query, usually one of ``"Homo sapiens"`` or ``"Mus musculus"``.
        measurement_name:
            The measurement object to query. Deafults to ``"RNA"``.

    Returns:
        A :class:`scipy.sparse.csr_array` object containing the presence matrix.

    Raises:
        ValueError: if the organism cannot be found.

    Lifecycle:
        maturing

    Examples:
        >>> get_presence_matrix(census, "Homo sapiens", "RNA")
        <321x60554 sparse array of type '<class 'numpy.uint8'>'
        with 6441269 stored elements in Compressed Sparse Row format>
    """
    exp = _get_experiment(census, organism)
    presence = exp.ms[measurement_name]["feature_dataset_presence_matrix"]
    return presence.read((slice(None),)).coos().concat().to_scipy().tocsr()



# Section: cellxgene_census-src-cellxgene_census-_util

import urllib.parse

import tiledbsoma as soma

USER_AGENT_ENVVAR = "CELLXGENE_CENSUS_USERAGENT"
"""Environment variable used to add more information into the user-agent."""


def _uri_join(base: str, url: str) -> str:
    """Like urllib.parse.urljoin, but doesn't get confused by s3://."""
    p_url = urllib.parse.urlparse(url)
    if p_url.netloc:
        return url

    p_base = urllib.parse.urlparse(base)
    path = urllib.parse.urljoin(p_base.path, p_url.path)
    parts = [
        p_base.scheme,
        p_base.netloc,
        path,
        p_url.params,
        p_url.query,
        p_url.fragment,
    ]
    return urllib.parse.urlunparse(parts)


def _extract_census_version(census: soma.Collection) -> str:
    """Extract the Census version from the given Census object."""
    try:
        version: str = urllib.parse.urlparse(census.uri).path.split("/")[2]
    except (KeyError, IndexError):
        raise ValueError("Unable to extract Census version.") from None

    return version


def _user_agent() -> str:
    import os

    import cellxgene_census

    if env_specifier := os.environ.get(USER_AGENT_ENVVAR, None):
        return f"cellxgene-census-python/{cellxgene_census.__version__} {env_specifier}"
    else:
        return f"cellxgene-census-python/{cellxgene_census.__version__}"



# Section: notebooks-api_demo-census_gget_demo

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




# Section: cellxgene_census-src-cellxgene_census-experimental-ml-huggingface-geneformer_tokenizer

import pickle
from collections.abc import Sequence
from typing import Any

import numpy as np
import numpy.typing as npt
import scipy
import tiledbsoma

from .cell_dataset_builder import CellDatasetBuilder


class GeneformerTokenizer(CellDatasetBuilder):
    """Generate a Hugging Face `Dataset` containing Geneformer token sequences for each
    cell in CELLxGENE Census ExperimentAxisQuery results (human).

    This class requires the Geneformer package to be installed separately with:
    `pip install git+https://huggingface.co/ctheodoris/Geneformer@eb038a6`

    Example usage:

    ```
    import cellxgene_census
    import tiledbsoma
    from cellxgene_census.experimental.ml.huggingface import GeneformerTokenizer

    with cellxgene_census.open_soma(census_version="latest") as census:
        with GeneformerTokenizer(
            census["census_data"]["homo_sapiens"],
            # set obs_query to define some subset of Census cells:
            obs_query=tiledbsoma.AxisQuery(value_filter="is_primary_data == True and tissue_general == 'tongue'"),
            obs_column_names=(
                "soma_joinid",
                "cell_type_ontology_term_id",
            ),
        ) as tokenizer:
            dataset = tokenizer.build()
    ```

    Dataset item contents:
    - `input_ids`: Geneformer token sequence for the cell
    - `length`: Length of the token sequence
    - and the specified `obs_column_names` (cell metadata from the experiment obs dataframe)
    """

    obs_column_names: set[str]
    max_input_tokens: int
    special_token: bool

    # Newer versions of Geneformer has a consolidated gene list (gene_mapping_file), meaning the
    # counts for one or more Census genes are to be summed to get the count for one Geneformer
    # gene. model_gene_map is a sparse binary matrix to map a cell vector (or multi-cell matrix) of
    # Census gene counts onto Geneformer gene counts. model_gene_map[i,j] is 1 iff the i'th Census
    # gene count contributes to the j'th Geneformer gene count.
    model_gene_map: scipy.sparse.coo_matrix
    model_gene_tokens: npt.NDArray[np.int64]  # Geneformer token for each column of model_gene_map
    model_gene_medians: npt.NDArray[np.float64]  # float for each column of model_gene_map
    model_cls_token: np.int64 | None = None
    model_eos_token: np.int64 | None = None

    def __init__(
        self,
        experiment: tiledbsoma.Experiment,
        *,
        obs_column_names: Sequence[str] | None = None,
        obs_attributes: Sequence[str] | None = None,
        max_input_tokens: int = 2048,
        special_token: bool = False,
        token_dictionary_file: str = "",
        gene_median_file: str = "",
        gene_mapping_file: str = "",
        **kwargs: Any,
    ) -> None:
        """Initialize GeneformerTokenizer.

        Args:
        - `experiment`: Census Experiment to query
        - `obs_query`: obs AxisQuery defining the set of Census cells to process (default all)
        - `obs_column_names`: obs dataframe columns (cell metadata) to propagate into attributes
           of each Dataset item
        - `max_input_tokens`: maximum length of Geneformer input token sequence (default 2048)
        - `special_token`: whether to affix separator tokens to the sequence (default False)
        - `token_dictionary_file`, `gene_median_file`: pickle files supplying the mapping of
          Ensembl human gene IDs onto Geneformer token numbers and median expression values.
          By default, these will be loaded from the Geneformer package.
        - `gene_mapping_file`: optional pickle file with mapping for Census gene IDs to model's
        """
        if obs_attributes:  # old name of obs_column_names
            obs_column_names = obs_attributes

        self.max_input_tokens = max_input_tokens
        self.special_token = special_token
        self.obs_column_names = set(obs_column_names) if obs_column_names else set()
        self._load_geneformer_data(experiment, token_dictionary_file, gene_median_file, gene_mapping_file)
        super().__init__(
            experiment,
            measurement_name="RNA",
            layer_name="raw",
            **kwargs,
        )

    def _load_geneformer_data(
        self,
        experiment: tiledbsoma.Experiment,
        token_dictionary_file: str,
        gene_median_file: str,
        gene_mapping_file: str,
    ) -> None:
        """Load (1) the experiment's genes dataframe and (2) Geneformer's static data
        files for gene tokens and median expression; then, intersect them to compute
        self.model_gene_{ids,tokens,medians}.
        """
        # TODO: this work could be reused for all queries on this experiment

        genes_df = (
            experiment.ms["RNA"]
            .var.read(column_names=["soma_joinid", "feature_id"])
            .concat()
            .to_pandas()
            .set_index("soma_joinid")
        )

        if not (token_dictionary_file and gene_median_file):
            try:
                import geneformer
            except ImportError:
                # pyproject.toml can't express Geneformer git+https dependency
                raise ImportError(
                    "Please install Geneformer with: "
                    "pip install git+https://huggingface.co/ctheodoris/Geneformer@eb038a6"
                ) from None
            if not token_dictionary_file:
                token_dictionary_file = geneformer.tokenizer.TOKEN_DICTIONARY_FILE
            if not gene_median_file:
                gene_median_file = geneformer.tokenizer.GENE_MEDIAN_FILE
        with open(token_dictionary_file, "rb") as f:
            gene_token_dict = pickle.load(f)
        with open(gene_median_file, "rb") as f:
            gene_median_dict = pickle.load(f)

        gene_mapping = None
        if gene_mapping_file:
            with open(gene_mapping_file, "rb") as f:
                gene_mapping = pickle.load(f)

        # compute model_gene_{ids,tokens,medians} by joining genes_df with Geneformer's
        # dicts
        map_data = []
        map_i = []
        map_j = []
        model_gene_id_by_ensg: dict[str, int] = {}
        model_gene_count = 0
        model_gene_tokens: list[np.int64] = []
        model_gene_medians: list[np.float64] = []
        for gene_id, row in genes_df.iterrows():
            ensg = row["feature_id"]  # ENSG... gene id, which keys Geneformer's dicts
            if gene_mapping is not None:
                ensg = gene_mapping.get(ensg, ensg)
            if ensg in gene_token_dict:
                if ensg not in model_gene_id_by_ensg:
                    model_gene_id_by_ensg[ensg] = model_gene_count
                    model_gene_count += 1
                    model_gene_tokens.append(gene_token_dict[ensg])
                    model_gene_medians.append(gene_median_dict[ensg])
                map_data.append(1)
                map_i.append(gene_id)
                map_j.append(model_gene_id_by_ensg[ensg])

        self.model_gene_map = scipy.sparse.coo_matrix(
            (map_data, (map_i, map_j)), shape=(genes_df.index.max() + 1, model_gene_count), dtype=bool
        )
        self.model_gene_tokens = np.array(model_gene_tokens, dtype=np.int64)
        self.model_gene_medians = np.array(model_gene_medians, dtype=np.float64)

        assert len(np.unique(self.model_gene_tokens)) == len(self.model_gene_tokens)
        assert np.all(self.model_gene_medians > 0)
        # Geneformer models protein-coding and miRNA genes, so the intersection should
        # be north of 18K.
        assert (
            model_gene_count > 18_000
        ), f"Mismatch between Census gene IDs and Geneformer token dicts (only {model_gene_count} common genes)"

        # Precompute a vector by which we'll multiply each cell's expression vector.
        # The denominator normalizes by Geneformer's median expression values.
        # The numerator 10K factor follows Geneformer's tokenizer; theoretically it doesn't affect
        # affect the rank order, but is probably intended to help with numerical precision.
        self.model_gene_medians_factor = 10_000.0 / self.model_gene_medians

        if self.special_token:
            self.model_cls_token = gene_token_dict["<cls>"]
            self.model_eos_token = gene_token_dict["<eos>"]

    def __enter__(self) -> "GeneformerTokenizer":
        super().__enter__()
        # On context entry, load the necessary cell metadata (obs_df)
        obs_column_names = list(self.obs_column_names)
        if "soma_joinid" not in self.obs_column_names:
            obs_column_names.append("soma_joinid")
        self.obs_df = self.obs(column_names=obs_column_names).concat().to_pandas().set_index("soma_joinid")
        return self

    def cell_item(self, cell_joinid: int, cell_Xrow: scipy.sparse.csr_matrix) -> dict[str, Any]:
        """Given the expression vector for one cell, compute the Dataset item providing
        the Geneformer inputs (token sequence and metadata).
        """
        # Apply model_gene_map to cell_Xrow and normalize with row sum & gene medians.
        # Notice we divide by the total count of the complete row (not only of the projected
        # values); this follows Geneformer's internal tokenizer.
        model_expr = (cell_Xrow * self.model_gene_map).multiply(self.model_gene_medians_factor / cell_Xrow.sum())
        assert isinstance(model_expr, scipy.sparse.coo_matrix), type(model_expr)
        assert model_expr.shape == (1, self.model_gene_map.shape[1])

        # figure the resulting tokens in descending order of model_expr
        # (use sparse model_expr.{col,data} to naturally exclude undetected genes)
        token_order = model_expr.col[np.argsort(-model_expr.data)[: self.max_input_tokens]]
        input_ids = self.model_gene_tokens[token_order]

        if self.special_token:
            # affix special tokens, dropping the last two gene tokens if necessary
            if len(input_ids) == self.max_input_tokens:
                input_ids = input_ids[:-1]
            assert self.model_cls_token is not None
            input_ids = np.insert(input_ids, 0, self.model_cls_token)
            if len(input_ids) == self.max_input_tokens:
                input_ids = input_ids[:-1]
            assert self.model_eos_token is not None
            input_ids = np.append(input_ids, self.model_eos_token)

        ans = {"input_ids": input_ids, "length": len(input_ids)}
        # add the requested obs attributes
        for attr in self.obs_column_names:
            if attr != "soma_joinid":
                ans[attr] = self.obs_df.at[cell_joinid, attr]
            else:
                ans["soma_joinid"] = cell_joinid
        return ans



# Section: cellxgene_census-src-cellxgene_census-experimental-ml-encoders

import abc
import functools

import numpy.typing as npt
import pandas as pd
from sklearn.preprocessing import LabelEncoder as SklearnLabelEncoder


class Encoder(abc.ABC):
    """Base class for ``obs`` encoders.

    To define a custom encoder, five methods must be implemented:

    - ``fit``: defines how the encoder will be fitted to the data.
    - ``transform``: defines how the encoder will be applied to the data
      in order to create an ``obs`` tensor.
    - ``inverse_transform``: defines how to decode the encoded values back
      to the original values.
    - ``name``: The name of the encoder. This will be used as the key in the
      dictionary of encoders. Each encoder passed to a :class:`.pytorch.ExperimentDataPipe` must have a unique name.
    - ``columns``: List of columns in ``obs`` that the encoder will be applied to.

    See the implementation of :class:`LabelEncoder` for an example.
    """

    @abc.abstractmethod
    def fit(self, obs: pd.DataFrame) -> None:
        """Fit the encoder with obs."""
        pass

    @abc.abstractmethod
    def transform(self, df: pd.DataFrame) -> pd.DataFrame:
        """Transform the obs :class:`pandas.DataFrame` into a :class:`pandas.DataFrame` of encoded values."""
        pass

    @abc.abstractmethod
    def inverse_transform(self, encoded_values: npt.ArrayLike) -> npt.ArrayLike:
        """Inverse transform the encoded values back to the original values."""
        pass

    @property
    @abc.abstractmethod
    def name(self) -> str:
        """Name of the encoder."""
        pass

    @property
    @abc.abstractmethod
    def columns(self) -> list[str]:
        """Columns in ``obs`` that the encoder will be applied to."""
        pass


class LabelEncoder(Encoder):
    """Default encoder based on :class:`sklearn.preprocessing.LabelEncoder`."""

    def __init__(self, col: str) -> None:
        self._encoder = SklearnLabelEncoder()
        self.col = col

    def fit(self, obs: pd.DataFrame) -> None:
        """Fit the encoder with ``obs``."""
        self._encoder.fit(obs[self.col].unique())

    def transform(self, df: pd.DataFrame) -> pd.DataFrame:
        """Transform the obs :class:`pandas.DataFrame` into a :class:`pandas.DataFrame` of encoded values."""
        return self._encoder.transform(df[self.col])  # type: ignore

    def inverse_transform(self, encoded_values: npt.ArrayLike) -> npt.ArrayLike:
        """Inverse transform the encoded values back to the original values."""
        return self._encoder.inverse_transform(encoded_values)  # type: ignore

    @property
    def name(self) -> str:
        """Name of the encoder."""
        return self.col

    @property
    def columns(self) -> list[str]:
        """Columns in ``obs`` that the encoder will be applied to."""
        return [self.col]

    @property
    def classes_(self):  # type: ignore
        """Classes of the encoder."""
        return self._encoder.classes_


class BatchEncoder(Encoder):
    """An encoder that concatenates and encodes several ``obs`` columns."""

    def __init__(self, cols: list[str], name: str = "batch"):
        self.cols = cols
        from sklearn.preprocessing import LabelEncoder

        self._name = name
        self._encoder = LabelEncoder()

    def _join_cols(self, df: pd.DataFrame):  # type: ignore
        return functools.reduce(lambda a, b: a + b, [df[c].astype(str) for c in self.cols])

    def transform(self, df: pd.DataFrame) -> pd.DataFrame:
        """Transform the obs :class:`pandas.DataFrame` into a :class:`pandas.DataFrame` of encoded values."""
        arr = self._join_cols(df)
        return self._encoder.transform(arr)  # type: ignore

    def inverse_transform(self, encoded_values: npt.ArrayLike) -> npt.ArrayLike:
        """Inverse transform the encoded values back to the original values."""
        return self._encoder.inverse_transform(encoded_values)  # type: ignore

    def fit(self, obs: pd.DataFrame) -> None:
        """Fit the encoder with ``obs``."""
        arr = self._join_cols(obs)
        self._encoder.fit(arr.unique())

    @property
    def columns(self) -> list[str]:
        """Columns in ``obs`` that the encoder will be applied to."""
        return self.cols

    @property
    def name(self) -> str:
        """Name of the encoder."""
        return self._name

    @property
    def classes_(self):  # type: ignore
        """Classes of the encoder."""
        return self._encoder.classes_



# Section: notebooks-api_demo-census_embedding

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




# Section: cellxgene_census-tests-experimental-pp-test_stats

from typing import Any

import numpy as np
import numpy.ma as ma
import pandas as pd
import pytest
import tiledbsoma as soma
from scipy import sparse

import cellxgene_census
from cellxgene_census.experimental import pp


def var(X: sparse.csc_matrix | sparse.csr_matrix, axis: int = 0, ddof: int = 1) -> Any:
    """
    Variance of a sparse matrix calculated as mean(X**2) - mean(X)**2
    with Bessel's correction applied for unbiased estimate
    """
    X_squared = X.copy()
    X_squared.data **= 2
    n = X.shape[axis]
    return ((X_squared.sum(axis=axis).A1 / n) - np.square(X.sum(axis=axis).A1 / n)) * (n / (n - ddof))


@pytest.mark.experimental
@pytest.mark.live_corpus
@pytest.mark.parametrize("axis", [0, 1])
@pytest.mark.parametrize("calc_mean,calc_variance", [(True, True), (True, False), (False, True)])
@pytest.mark.parametrize(
    "experiment_name,obs_value_filter,obs_coords",
    [
        ("mus_musculus", 'tissue_general == "liver" and is_primary_data == True', ()),
        ("mus_musculus", 'is_primary_data == True and tissue_general == "heart"', ()),
        pytest.param(
            "mus_musculus",
            "is_primary_data == True",
            (slice(0, 400_000),),
            marks=pytest.mark.expensive,
        ),
        pytest.param(
            "homo_sapiens",
            "is_primary_data == True",
            (slice(0, 400_000),),
            marks=pytest.mark.expensive,
        ),
    ],
)
def test_mean_variance(
    experiment_name: str,
    obs_value_filter: str,
    axis: int,
    calc_mean: bool,
    calc_variance: bool,
    small_mem_context: soma.SOMATileDBContext,
    obs_coords: tuple[None, slice],
) -> None:
    with cellxgene_census.open_soma(census_version="latest", context=small_mem_context) as census:
        with census["census_data"][experiment_name].axis_query(
            measurement_name="RNA",
            obs_query=soma.AxisQuery(value_filter=obs_value_filter, coords=obs_coords),
        ) as query:
            mean_variance = pp.mean_variance(
                query,
                calculate_mean=calc_mean,
                calculate_variance=calc_variance,
                axis=axis,
            )
            assert isinstance(mean_variance, pd.DataFrame)
            if calc_mean:
                assert "mean" in mean_variance
                assert mean_variance["mean"].dtype == np.float64

            if calc_variance:
                assert "variance" in mean_variance
                assert mean_variance["variance"].dtype == np.float64

            if not calc_mean:
                assert "mean" not in mean_variance
            if not calc_variance:
                assert "variance" not in mean_variance

            assert mean_variance.index.name == "soma_joinid"
            if axis == 0:
                assert np.array_equal(mean_variance.index, query.var_joinids())
            else:
                assert np.array_equal(mean_variance.index, query.obs_joinids())

            table = query.X("raw").tables().concat()
            data = table["soma_data"].to_numpy()

            dim_0 = query.indexer.by_obs(table["soma_dim_0"])
            dim_1 = query.indexer.by_var(table["soma_dim_1"])
            coo = sparse.coo_matrix((data, (dim_0, dim_1)), shape=(query.n_obs, query.n_vars))

            if calc_mean:
                mean = coo.mean(axis=axis)
                if axis == 1:
                    mean = mean.T
                assert np.allclose(mean, mean_variance["mean"], atol=1e-5, rtol=1e-2)

            if calc_variance:
                variance = var(coo, axis=axis)
                assert np.allclose(variance, mean_variance["variance"], atol=1e-5, rtol=1e-2)


@pytest.mark.experimental
@pytest.mark.live_corpus
@pytest.mark.parametrize("axis", [0, 1])
@pytest.mark.parametrize("calc_mean,calc_variance", [(True, True)])
@pytest.mark.parametrize(
    "experiment_name,obs_coords",
    [
        ("mus_musculus", (slice(0, 1_000),)),
    ],
)
def test_mean_variance_nnz_only(
    experiment_name: str,
    axis: int,
    calc_mean: bool,
    calc_variance: bool,
    small_mem_context: soma.SOMATileDBContext,
    obs_coords: tuple[None, slice],
) -> None:
    # Note: since this test requires materializing the matrix in memory to compute the mean/variance,
    # we're going to use a coord slice based approach. This will ensure the matrix can fit in memory.
    with cellxgene_census.open_soma(census_version="latest", context=small_mem_context) as census:
        with census["census_data"][experiment_name].axis_query(
            measurement_name="RNA", obs_query=soma.AxisQuery(coords=obs_coords)
        ) as query:
            mean_variance = pp.mean_variance(
                query,
                calculate_mean=calc_mean,
                calculate_variance=calc_variance,
                axis=axis,
                nnz_only=True,
                ddof=0,
            )

            table = query.X("raw").tables().concat()
            data = table["soma_data"].to_numpy()

            dim_0 = query.indexer.by_obs(table["soma_dim_0"])
            dim_1 = query.indexer.by_var(table["soma_dim_1"])
            coo = sparse.coo_matrix((data, (dim_0, dim_1)), shape=(query.n_obs, query.n_vars))

            dense = coo.toarray()

            mask = np.ones(coo.shape)
            r, c = coo.nonzero()
            for x, y in zip(r, c):
                mask[x, y] = 0
            masked = ma.masked_array(dense, mask=mask)  # type: ignore[no-untyped-call, var-annotated]

            if calc_mean:
                mean = masked.mean(axis=axis)  # type: ignore[no-untyped-call]
                assert np.allclose(mean, mean_variance["mean"], atol=1e-5, rtol=1e-1, equal_nan=True)

            if calc_variance:
                variance = masked.var(axis=axis, ddof=0)  # type: ignore[no-untyped-call]
                va = mean_variance["variance"].to_numpy()
                assert np.allclose(variance, va, atol=1e-5, rtol=1e-2, equal_nan=True)


@pytest.mark.experimental
def test_mean_variance_no_flags() -> None:
    with pytest.raises(ValueError):
        pp.mean_variance(soma.AxisQuery(), calculate_mean=False, calculate_variance=False)


@pytest.mark.parametrize("experiment_name", ["mus_musculus"])
def test_mean_variance_empty_query(experiment_name: str, small_mem_context: soma.SOMATileDBContext) -> None:
    with cellxgene_census.open_soma(census_version="latest", context=small_mem_context) as census:
        with census["census_data"][experiment_name].axis_query(
            measurement_name="RNA",
            obs_query=soma.AxisQuery(value_filter='tissue_general == "foo"'),
        ) as query:
            with pytest.raises(ValueError):
                pp.mean_variance(query, calculate_mean=True, calculate_variance=True)


@pytest.mark.experimental
def test_mean_variance_wrong_axis() -> None:
    with pytest.raises(ValueError):
        pp.mean_variance(soma.AxisQuery(), calculate_mean=True, calculate_variance=True, axis=2)



# Section: cellxgene_census-tests-test_directory

from typing import Any

import pytest
import requests_mock as rm
import s3fs

import cellxgene_census
from cellxgene_census._release_directory import (
    CELL_CENSUS_MIRRORS_DIRECTORY_URL,
    CELL_CENSUS_RELEASE_DIRECTORY_URL,
)

# This test fixture contains 3 releases: 1 "latest" and 2 "LTS". Of the "LTS" releases, one is aliased to "stable"
# and one is "retracted", and both are aliased with "V#" aliases. The ordering of the releases is
# explicitly set to verify that the directory is sorted correctly (i.e. they are not in the desired order here).
# There is also a "dangling" tag, to verify that we handle this case correctly.
#
# TODO: Break this into multiple fixtures to test different scenarios
DIRECTORY_JSON = {
    "2022-10-01": {
        "release_date": "2022-10-30",
        "release_build": "2022-10-01",
        "flags": {"lts": True},
        "soma": {
            "uri": "s3://cellxgene-data-public/cell-census/2022-10-01/soma/",
            "s3_region": "us-west-2",
        },
        "h5ads": {
            "uri": "s3://cellxgene-data-public/cell-census/2022-10-01/h5ads/",
            "s3_region": "us-west-2",
        },
    },
    "2022-09-01": {
        "release_date": "2022-09-30",
        "release_build": "2022-09-01",
        "flags": {"lts": True, "retracted": True},
        "retraction": {
            "date": "2022-11-15",
            "reason": "mistakes happen",
            "info_permalink": "http://cellxgene.com/census/apologies",
        },
        "soma": {
            "uri": "s3://cellxgene-data-public/cell-census/2022-09-01/soma/",
            "s3_region": "us-west-2",
        },
        "h5ads": {
            "uri": "s3://cellxgene-data-public/cell-census/2022-09-01/h5ads/",
            "s3_region": "us-west-2",
        },
    },
    # Ordered the latest release to be last, to verify it is explicitly sorted
    "2022-11-01": {
        "release_date": "2022-11-30",
        "release_build": "2022-11-01",
        "soma": {
            "uri": "s3://cellxgene-data-public/cell-census/2022-11-01/soma/",
            "s3_region": "us-west-2",
        },
        "h5ads": {
            "uri": "s3://cellxgene-data-public/cell-census/2022-11-01/h5ads/",
            "s3_region": "us-west-2",
        },
    },
    # An explicitly dangling tag, to confirm we handle correct
    # Underscore indicates expected failure to test below
    "_dangling": "no-such-tag",
    # Aliases placed at bottom, to verify these are explicitly sorted to the top
    "stable": "V2",
    "latest": "2022-11-01",
    "V2": "2022-10-01",
    "V1": "2022-09-01",
}

MIRRORS_JSON = {
    "default": "AWS-S3-us-west-2",
    "AWS-S3-us-west-2": {
        "provider": "S3",
        "base_uri": "s3://cellxgene-data-public/",
        "region": "us-west-2",
    },
}


@pytest.fixture
def directory_mock(requests_mock: rm.Mocker) -> Any:
    return requests_mock.get(CELL_CENSUS_RELEASE_DIRECTORY_URL, json=DIRECTORY_JSON)


@pytest.fixture
def mirrors_mock(requests_mock: rm.Mocker) -> Any:
    return requests_mock.get(CELL_CENSUS_MIRRORS_DIRECTORY_URL, json=MIRRORS_JSON)


def test_get_census_version_directory(directory_mock: Any) -> None:
    directory = cellxgene_census.get_census_version_directory()

    assert isinstance(directory, dict)
    assert len(directory) > 0
    assert all(isinstance(k, str) for k in directory.keys())
    assert all(isinstance(v, dict) for v in directory.values())

    assert "_dangling" not in directory

    assert directory["2022-11-01"] == DIRECTORY_JSON["2022-11-01"]
    assert directory["2022-10-01"] == DIRECTORY_JSON["2022-10-01"]
    assert "2022-09-01" not in directory  # retracted excluded by default

    assert directory["latest"] == DIRECTORY_JSON["2022-11-01"]
    assert directory["stable"] == DIRECTORY_JSON["2022-10-01"]
    assert directory["V2"] == DIRECTORY_JSON["2022-10-01"]
    assert "V1" not in directory  # retracted excluded by default

    for tag in directory:
        assert directory[tag] == cellxgene_census.get_census_version_description(tag)

    # Verify that the directory is sorted according to this criteria:
    # 1. Aliases first
    # 2. Non aliases after, in reverse order
    dir_list = list(directory)
    assert dir_list == ["stable", "latest", "V2", "2022-11-01", "2022-10-01"]


def test_get_census_version_directory__lts_only(directory_mock: Any) -> None:
    directory = cellxgene_census.get_census_version_directory(lts=True)

    assert directory.keys() == {"stable", "V2", "2022-10-01"}


def test_get_census_version_directory__exclude_lts(directory_mock: Any) -> None:
    directory = cellxgene_census.get_census_version_directory(lts=False)

    assert directory.keys() == {"latest", "2022-11-01"}


def test_get_census_version_directory__include_retracted(directory_mock: Any) -> None:
    directory = cellxgene_census.get_census_version_directory(retracted=None)

    assert "V1" in directory
    assert "2022-09-01" in directory


def test_get_census_version_directory__retraction_info(directory_mock: Any) -> None:
    directory = cellxgene_census.get_census_version_directory(retracted=True)

    assert directory["2022-09-01"]["retraction"] == {
        "date": "2022-11-15",
        "reason": "mistakes happen",
        "info_permalink": "http://cellxgene.com/census/apologies",
    }

    assert directory["V1"]["retraction"] == {
        "date": "2022-11-15",
        "reason": "mistakes happen",
        "info_permalink": "http://cellxgene.com/census/apologies",
    }


def test_get_census_version_description_errors() -> None:
    with pytest.raises(ValueError):
        cellxgene_census.get_census_version_description(census_version="no/such/version/exists")


def test_get_census_mirrors_directory(mirrors_mock: Any) -> None:
    directory = cellxgene_census.get_census_mirror_directory()
    assert "default" not in directory
    assert "AWS-S3-us-west-2" in directory
    assert directory["AWS-S3-us-west-2"] == MIRRORS_JSON["AWS-S3-us-west-2"]


@pytest.mark.live_corpus
def test_live_directory_contents() -> None:
    # Sanity check that all directory contents are usable. This uses the
    # live directory, so it _could_ start failing without a code change.
    # But given the purpose of this package, that seems like a reasonable
    # tradeoff, as the data directory should never be "corrupt" or there
    # is widespread impact on users.

    fs = s3fs.S3FileSystem(anon=True, cache_regions=True)

    directory = cellxgene_census.get_census_version_directory()
    assert "latest" in directory

    for version, version_description in directory.items():
        with cellxgene_census.open_soma(census_version=version) as census:
            assert census is not None

        assert fs.exists(version_description["soma"]["uri"])
        assert fs.exists(version_description["h5ads"]["uri"])


def test_census_version_types() -> None:
    """Do a little bit of runtime type checking on the results of census version functions.

    Part of solving: https://github.com/chanzuckerberg/cellxgene-census/issues/1204
    """
    from cellxgene_census._release_directory import CensusVersionDescription

    directory = cellxgene_census.get_census_version_directory()
    for k, v in directory.items():
        assert set(v).issubset(CensusVersionDescription.__annotations__)
        desc = cellxgene_census.get_census_version_description(k)
        assert set(desc).issubset(CensusVersionDescription.__annotations__)



# Section: cellxgene_census-src-cellxgene_census-experimental-ml-__init__

"""An API to facilitate use of PyTorch ML training with data from the CZI Science CELLxGENE Census."""

from .encoders import BatchEncoder, Encoder, LabelEncoder
from .pytorch import ExperimentDataPipe, Stats, experiment_dataloader

__all__ = [
    "Stats",
    "ExperimentDataPipe",
    "experiment_dataloader",
    "Encoder",
    "LabelEncoder",
    "BatchEncoder",
]



# Section: cellxgene_census-tests-conftest

import multiprocessing

import pytest
import tiledbsoma as soma

TEST_MARKERS_SKIPPED_BY_DEFAULT = ["expensive", "experimental"]

# tiledb will complain if this isn't set and a process is spawned. May cause segfaults on the proxy test if this isn't set.
multiprocessing.set_start_method("spawn", force=True)


def pytest_addoption(parser: pytest.Parser) -> None:
    for test_option in TEST_MARKERS_SKIPPED_BY_DEFAULT:
        parser.addoption(
            f"--{test_option}",
            action="store_true",
            dest=test_option,
            default=False,
            help=f"enable '{test_option}' decorated tests",
        )

    # Add option to set the census_version (not set by default)
    parser.addoption("--census_version", action="store", default="stable")


def pytest_generate_tests(metafunc: pytest.Metafunc) -> None:
    """This is called for every test"""

    # Configure census_version if used
    census_version = metafunc.config.option.census_version
    if "census_version" in metafunc.fixturenames:
        metafunc.parametrize("census_version", [census_version])


def pytest_configure(config: pytest.Config) -> None:
    """
    Exclude tests marked with any of the TEST_MARKERS_SKIPPED_BY_DEFAULT values, unless the corresponding explicit
    flag is specified by the user.
    """
    excluded_markexprs = []

    for test_option in TEST_MARKERS_SKIPPED_BY_DEFAULT:
        if not vars(config.option).get(test_option, False):
            excluded_markexprs.append(test_option)

    if config.option.markexpr and excluded_markexprs:
        config.option.markexpr += " and "
    config.option.markexpr += " and ".join([f"not {m}" for m in excluded_markexprs])


@pytest.fixture
def small_mem_context() -> soma.SOMATileDBContext:
    """used to keep memory usage smaller for GHA runners."""
    from cellxgene_census import get_default_soma_context

    return get_default_soma_context(tiledb_config={"soma.init_buffer_bytes": 32 * 1024**2})


@pytest.fixture(scope="session")
def census() -> soma.Collection:
    import cellxgene_census

    return cellxgene_census.open_soma(census_version="latest")


@pytest.fixture(scope="session")
def lts_census() -> soma.Collection:
    import cellxgene_census

    return cellxgene_census.open_soma(census_version="stable")


@pytest.fixture(scope="session")
def dec_lts_census() -> soma.Collection:
    """Fixture for the 2023-12-15 LTS Census."""
    import cellxgene_census

    return cellxgene_census.open_soma(census_version="2023-12-15")



# Section: notebooks-analysis_demo-comp_bio_census_info

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




# Section: cellxgene_census-tests-experimental-test_embeddings_search

import json
from typing import Any

import anndata as ad
import numpy as np
import pytest

import cellxgene_census
from cellxgene_census.experimental import (
    NeighborObs,
    find_nearest_obs,
    predict_obs_metadata,
)


@pytest.mark.experimental
@pytest.mark.live_corpus
def test_embeddings_search(true_neighbors: dict[str, Any], query_result: NeighborObs) -> None:
    # check result shapes
    rslt = query_result
    assert isinstance(rslt.neighbor_ids, np.ndarray)
    assert rslt.neighbor_ids.dtype == np.uint64
    assert rslt.neighbor_ids.shape == (len(true_neighbors), TRUE_NEAREST_NEIGHBORS_K)
    assert isinstance(rslt.distances, np.ndarray)
    assert rslt.distances.dtype == np.float32
    assert rslt.distances.shape == (len(true_neighbors), TRUE_NEAREST_NEIGHBORS_K)

    # compute Jaccard index for the true neighbors & those returned by the query
    true_ids = set()
    for ns in true_neighbors.values():
        true_ids |= {n["cell_id"] for n in ns}
    rslt_ids = set(rslt.neighbor_ids.flatten())
    jaccard = len(true_ids & rslt_ids) / len(true_ids | rslt_ids)

    # Jaccard threshold: the search is approximate, and we set nprobe to a low value to speed up
    # the test
    assert jaccard >= 0.92

    return


@pytest.mark.experimental
@pytest.mark.live_corpus
@pytest.mark.parametrize("n_neighbors", [5, 7, 20])
def test_embedding_search_n_neighbors(query_anndata: ad.AnnData, n_neighbors: int) -> None:
    columns = ["cell_type"]
    result = find_nearest_obs(
        TRUE_NEAREST_NEIGHBORS_EMBEDDING_NAME,
        TRUE_NEAREST_NEIGHBORS_ORGANISM,
        TRUE_NEAREST_NEIGHBORS_CENSUS_VERSION,
        query_anndata,
        k=n_neighbors,
        nprobe=25,
    )

    # Check that the correct number of neighbors is being returned
    assert result.neighbor_ids.shape[1] == n_neighbors
    # Check that this step works
    _ = predict_obs_metadata(TRUE_NEAREST_NEIGHBORS_ORGANISM, TRUE_NEAREST_NEIGHBORS_CENSUS_VERSION, result, columns)


@pytest.mark.experimental
@pytest.mark.live_corpus
def test_embeddings_search_errors(query_anndata: ad.AnnData) -> None:
    # bogus embedding name
    with pytest.raises(ValueError, match="No embeddings found"):
        find_nearest_obs(
            "bogus123", TRUE_NEAREST_NEIGHBORS_ORGANISM, TRUE_NEAREST_NEIGHBORS_CENSUS_VERSION, query_anndata
        )
    # query anndata missing the embedding layer
    bogus_ad = query_anndata.copy()
    bogus_ad.obsm.pop(TRUE_NEAREST_NEIGHBORS_EMBEDDING_NAME)
    with pytest.raises(ValueError, match="Query does not have"):
        find_nearest_obs(
            TRUE_NEAREST_NEIGHBORS_EMBEDDING_NAME,
            TRUE_NEAREST_NEIGHBORS_ORGANISM,
            TRUE_NEAREST_NEIGHBORS_CENSUS_VERSION,
            bogus_ad,
        )
    # embedding layer has wrong number of features
    bogus_ad = query_anndata.copy()
    bogus_ad.obsm[TRUE_NEAREST_NEIGHBORS_EMBEDDING_NAME] = np.zeros((len(bogus_ad), 42))
    with pytest.raises(ValueError, match="features, expected"):
        find_nearest_obs(
            TRUE_NEAREST_NEIGHBORS_EMBEDDING_NAME,
            TRUE_NEAREST_NEIGHBORS_ORGANISM,
            TRUE_NEAREST_NEIGHBORS_CENSUS_VERSION,
            bogus_ad,
        )
    return


@pytest.mark.experimental
@pytest.mark.live_corpus
def test_predict_obs_metadata(query_anndata: ad.AnnData, query_result: NeighborObs) -> None:
    columns = ["cell_type", "tissue_general"]

    with cellxgene_census.open_soma(census_version=TRUE_NEAREST_NEIGHBORS_CENSUS_VERSION) as census:
        truth_df = (
            census["census_data"][TRUE_NEAREST_NEIGHBORS_ORGANISM]
            .obs.read(coords=(query_anndata.obs["soma_joinid"].values,), column_names=columns)
            .concat()
            .to_pandas()
        )

    pred_df = predict_obs_metadata(
        TRUE_NEAREST_NEIGHBORS_ORGANISM, TRUE_NEAREST_NEIGHBORS_CENSUS_VERSION, query_result, columns
    )
    assert len(pred_df) == len(query_anndata.obs)

    for col in columns:
        assert col in pred_df.columns
        assert f"{col}_confidence" in pred_df.columns
        accuracy = (pred_df[col] == truth_df[col]).mean()
        assert accuracy > 0.75, f"Accuracy for {col} is {accuracy}"


@pytest.fixture(scope="module")
def true_neighbors() -> dict[int, list[dict[str, Any]]]:
    ans = {}
    for line in TRUE_NEAREST_NEIGHBORS_JSON.strip().split("\n"):
        example = json.loads(line)
        ans[example["cell_id"]] = example["neighbors"][:TRUE_NEAREST_NEIGHBORS_K]
    return ans


@pytest.fixture(scope="module")
def query_anndata(true_neighbors: dict[str, Any]) -> ad.AnnData:
    with cellxgene_census.open_soma(census_version=TRUE_NEAREST_NEIGHBORS_CENSUS_VERSION) as census:
        return cellxgene_census.get_anndata(
            census,
            TRUE_NEAREST_NEIGHBORS_ORGANISM,
            obs_coords=sorted(true_neighbors.keys()),
            obs_embeddings=[TRUE_NEAREST_NEIGHBORS_EMBEDDING_NAME],
        )


@pytest.fixture(scope="module")
def query_result(query_anndata: ad.AnnData) -> NeighborObs:
    return find_nearest_obs(
        TRUE_NEAREST_NEIGHBORS_EMBEDDING_NAME,
        TRUE_NEAREST_NEIGHBORS_ORGANISM,
        TRUE_NEAREST_NEIGHBORS_CENSUS_VERSION,
        query_anndata,
        k=TRUE_NEAREST_NEIGHBORS_K,
        nprobe=25,
    )


"""
The following "truth" data were generated by randomly sampling 100 Census primary cells and their
scVI embeddings (CxG-czi-5). Then, we made an exhaustive pass through the embeddings array to find
their 25 nearest neighbors by Euclidean distance.
"""
TRUE_NEAREST_NEIGHBORS_CENSUS_VERSION = "2023-12-15"
TRUE_NEAREST_NEIGHBORS_EMBEDDING_NAME = "scvi"
TRUE_NEAREST_NEIGHBORS_ORGANISM = "homo_sapiens"
TRUE_NEAREST_NEIGHBORS_K = 10
TRUE_NEAREST_NEIGHBORS_JSON = """
{"cell_id": 48884266, "neighbors": [{"distance": 0.0, "cell_id": 48884266}, {"distance": 0.36870312719526344, "cell_id": 48882582}, {"distance": 0.45194134685206555, "cell_id": 48904243}, {"distance": 0.4689120465432001, "cell_id": 48883973}, {"distance": 0.4701634704715092, "cell_id": 48883894}, {"distance": 0.4772844749267185, "cell_id": 48881548}, {"distance": 0.47747182817763884, "cell_id": 48883215}, {"distance": 0.5171140271570219, "cell_id": 48882290}, {"distance": 0.5343526701049458, "cell_id": 48905593}, {"distance": 0.5393365808038064, "cell_id": 48896713}, {"distance": 0.5409478030410546, "cell_id": 48906104}, {"distance": 0.5461267613033821, "cell_id": 48905835}, {"distance": 0.5572477334843036, "cell_id": 48896520}, {"distance": 0.5611744042351648, "cell_id": 48885279}, {"distance": 0.5652364617866981, "cell_id": 48882507}, {"distance": 0.5657781269186272, "cell_id": 48907142}, {"distance": 0.5709338004974971, "cell_id": 48599354}, {"distance": 0.5713592066534472, "cell_id": 48904840}, {"distance": 0.5791306505071526, "cell_id": 48905823}, {"distance": 0.5796388325518227, "cell_id": 48907988}, {"distance": 0.5827178977498736, "cell_id": 48906968}, {"distance": 0.5842219992223383, "cell_id": 48905256}, {"distance": 0.5884837832056955, "cell_id": 48907109}, {"distance": 0.5886426284651249, "cell_id": 48610323}, {"distance": 0.5891020805091439, "cell_id": 48906367}]}
{"cell_id": 3285087, "neighbors": [{"distance": 0.0, "cell_id": 3285087}, {"distance": 0.0, "cell_id": 3723278}, {"distance": 0.6446024157279429, "cell_id": 20236001}, {"distance": 0.6446024157279429, "cell_id": 23353329}, {"distance": 0.6446024157279429, "cell_id": 25246595}, {"distance": 0.6643407808380998, "cell_id": 3285053}, {"distance": 0.6643407808380998, "cell_id": 3723267}, {"distance": 0.6716110836369217, "cell_id": 3285058}, {"distance": 0.6716110836369217, "cell_id": 3723484}, {"distance": 0.6989909074224284, "cell_id": 20248724}, {"distance": 0.6989909074224284, "cell_id": 23335331}, {"distance": 0.6989909074224284, "cell_id": 29081778}, {"distance": 0.6995074309524025, "cell_id": 20252505}, {"distance": 0.6995074309524025, "cell_id": 23338866}, {"distance": 0.6995074309524025, "cell_id": 29083186}, {"distance": 0.7114577326312008, "cell_id": 20246681}, {"distance": 0.7114577326312008, "cell_id": 23344825}, {"distance": 0.7114577326312008, "cell_id": 29080630}, {"distance": 0.7287329371320457, "cell_id": 3285105}, {"distance": 0.7287329371320457, "cell_id": 3723290}, {"distance": 0.7380269501365168, "cell_id": 20253327}, {"distance": 0.7380269501365168, "cell_id": 23339993}, {"distance": 0.7380269501365168, "cell_id": 29080521}, {"distance": 0.7415650428518051, "cell_id": 20248884}, {"distance": 0.7415650428518051, "cell_id": 23335308}]}
{"cell_id": 44010765, "neighbors": [{"distance": 0.0, "cell_id": 44010765}, {"distance": 0.4570227796736593, "cell_id": 44126562}, {"distance": 0.48991132196560166, "cell_id": 45608463}, {"distance": 0.48991132196560166, "cell_id": 53050823}, {"distance": 0.49183571311954527, "cell_id": 15171731}, {"distance": 0.5129020377570104, "cell_id": 43877340}, {"distance": 0.5137430792855204, "cell_id": 15369399}, {"distance": 0.5143573085759233, "cell_id": 44402470}, {"distance": 0.5165583253624815, "cell_id": 50959996}, {"distance": 0.5310796622875035, "cell_id": 44062141}, {"distance": 0.5419616322548864, "cell_id": 50923746}, {"distance": 0.5452629175422571, "cell_id": 49985581}, {"distance": 0.549922152296617, "cell_id": 44080104}, {"distance": 0.5509699626225824, "cell_id": 43566024}, {"distance": 0.5620417941913957, "cell_id": 44016658}, {"distance": 0.5647218554055395, "cell_id": 44167855}, {"distance": 0.5669505839214761, "cell_id": 50226648}, {"distance": 0.5777917192400728, "cell_id": 16612512}, {"distance": 0.5878045199961383, "cell_id": 14256077}, {"distance": 0.5968481390701431, "cell_id": 43993358}, {"distance": 0.5970031957669568, "cell_id": 43232018}, {"distance": 0.597117836540424, "cell_id": 36301462}, {"distance": 0.599696952866418, "cell_id": 49383668}, {"distance": 0.5998371978218984, "cell_id": 49305572}, {"distance": 0.6020293495008047, "cell_id": 44083732}]}
{"cell_id": 54593156, "neighbors": [{"distance": 0.0, "cell_id": 54593156}, {"distance": 0.40684991419796035, "cell_id": 39702052}, {"distance": 0.4091920837433136, "cell_id": 54191080}, {"distance": 0.42898864616400934, "cell_id": 43867753}, {"distance": 0.4457620619504269, "cell_id": 53842427}, {"distance": 0.4462893109725141, "cell_id": 61354565}, {"distance": 0.454819493086827, "cell_id": 54233878}, {"distance": 0.46508571459284154, "cell_id": 48000680}, {"distance": 0.47447094231821035, "cell_id": 49335786}, {"distance": 0.4873961122332809, "cell_id": 48493624}, {"distance": 0.4886436932873957, "cell_id": 54194143}, {"distance": 0.48872276461918474, "cell_id": 45034245}, {"distance": 0.4909224675592821, "cell_id": 54365563}, {"distance": 0.4946262368544072, "cell_id": 46642468}, {"distance": 0.5050385955392035, "cell_id": 54186608}, {"distance": 0.5052535533766086, "cell_id": 46431610}, {"distance": 0.505342932652932, "cell_id": 39633433}, {"distance": 0.5102944171747317, "cell_id": 61333555}, {"distance": 0.5118471434841758, "cell_id": 53789514}, {"distance": 0.5147027601535609, "cell_id": 49902772}, {"distance": 0.519211855553082, "cell_id": 46485691}, {"distance": 0.5193234805436283, "cell_id": 54383558}, {"distance": 0.5207687400583412, "cell_id": 50065190}, {"distance": 0.5223131852274471, "cell_id": 53782932}, {"distance": 0.5225511137222412, "cell_id": 54204543}]}
{"cell_id": 49454989, "neighbors": [{"distance": 0.0, "cell_id": 49454989}, {"distance": 0.5651309620287789, "cell_id": 49196338}, {"distance": 0.5790658189328616, "cell_id": 61555029}, {"distance": 0.6358119018825575, "cell_id": 44261450}, {"distance": 0.6452191637719351, "cell_id": 35375516}, {"distance": 0.6483208161115366, "cell_id": 3819425}, {"distance": 0.6627657433551711, "cell_id": 35407515}, {"distance": 0.6681152419436311, "cell_id": 49783022}, {"distance": 0.6928590812438838, "cell_id": 49835917}, {"distance": 0.6929768400694352, "cell_id": 13072273}, {"distance": 0.6929768400694352, "cell_id": 13285326}, {"distance": 0.6942676543274021, "cell_id": 50172165}, {"distance": 0.7081249562541901, "cell_id": 49959583}, {"distance": 0.7117484089464857, "cell_id": 49611789}, {"distance": 0.7134029912676838, "cell_id": 49243741}, {"distance": 0.7207746269057522, "cell_id": 50081883}, {"distance": 0.725391822656998, "cell_id": 35771755}, {"distance": 0.725391822656998, "cell_id": 36110319}, {"distance": 0.7303628933089745, "cell_id": 61570231}, {"distance": 0.7364214270178401, "cell_id": 50289965}, {"distance": 0.7365543235391223, "cell_id": 61429320}, {"distance": 0.7372731447756735, "cell_id": 46851442}, {"distance": 0.7423110707448584, "cell_id": 51328234}, {"distance": 0.750998436167621, "cell_id": 13058884}, {"distance": 0.750998436167621, "cell_id": 13277399}]}
{"cell_id": 26103963, "neighbors": [{"distance": 0.0, "cell_id": 10781646}, {"distance": 0.0, "cell_id": 23881353}, {"distance": 0.0, "cell_id": 26103963}, {"distance": 0.619592659769135, "cell_id": 17736061}, {"distance": 0.619592659769135, "cell_id": 42189319}, {"distance": 0.6489919810741452, "cell_id": 6324230}, {"distance": 0.6489919810741452, "cell_id": 7388047}, {"distance": 0.6732222479959671, "cell_id": 23395}, {"distance": 0.6911732860546104, "cell_id": 1631450}, {"distance": 0.6911732860546104, "cell_id": 2510944}, {"distance": 0.6982662351259253, "cell_id": 10791296}, {"distance": 0.6982662351259253, "cell_id": 23877929}, {"distance": 0.6982662351259253, "cell_id": 24150236}, {"distance": 0.7182665789398123, "cell_id": 1636210}, {"distance": 0.7182665789398123, "cell_id": 2515704}, {"distance": 0.7520928971675964, "cell_id": 137816}, {"distance": 0.7555261454814952, "cell_id": 5719829}, {"distance": 0.7555261454814952, "cell_id": 7396419}, {"distance": 0.7585939977882471, "cell_id": 10791482}, {"distance": 0.7585939977882471, "cell_id": 23789797}, {"distance": 0.7585939977882471, "cell_id": 23867139}, {"distance": 0.7671337844431296, "cell_id": 1603547}, {"distance": 0.7671337844431296, "cell_id": 2251352}, {"distance": 0.7804735398208412, "cell_id": 10769175}, {"distance": 0.7804735398208412, "cell_id": 23899965}]}
{"cell_id": 48403678, "neighbors": [{"distance": 0.0, "cell_id": 48403678}, {"distance": 0.35710639363870667, "cell_id": 45162688}, {"distance": 0.551491405855107, "cell_id": 48216524}, {"distance": 0.5650419558705774, "cell_id": 49143947}, {"distance": 0.7438835677911204, "cell_id": 53988429}, {"distance": 0.7706764170797038, "cell_id": 8881701}, {"distance": 0.8075520431790514, "cell_id": 48147826}, {"distance": 0.8168826090755377, "cell_id": 43802274}, {"distance": 0.8424177203693685, "cell_id": 8881230}, {"distance": 0.8474836807859087, "cell_id": 46862802}, {"distance": 0.8551173717229851, "cell_id": 50827187}, {"distance": 0.8594960922337449, "cell_id": 48265069}, {"distance": 0.863042103716869, "cell_id": 48292616}, {"distance": 0.8672058039273095, "cell_id": 50027013}, {"distance": 0.8672988948905219, "cell_id": 54462677}, {"distance": 0.8763203554833326, "cell_id": 8883234}, {"distance": 0.8833765105762164, "cell_id": 36629863}, {"distance": 0.8863749961068852, "cell_id": 8883048}, {"distance": 0.8871093423506627, "cell_id": 49294629}, {"distance": 0.8877135497885864, "cell_id": 43922548}, {"distance": 0.8885346967652048, "cell_id": 3800394}, {"distance": 0.8992515425956116, "cell_id": 48263226}, {"distance": 0.9041501832203231, "cell_id": 48264984}, {"distance": 0.9050695821671393, "cell_id": 7674974}, {"distance": 0.9056283323559379, "cell_id": 8886159}]}
{"cell_id": 53784183, "neighbors": [{"distance": 0.0, "cell_id": 53784183}, {"distance": 0.5445907191962778, "cell_id": 54257572}, {"distance": 0.6490806955896615, "cell_id": 54260572}, {"distance": 0.6829816533618465, "cell_id": 54206161}, {"distance": 0.6848687435822737, "cell_id": 54261678}, {"distance": 0.7115714355716817, "cell_id": 36443600}, {"distance": 0.7162098549228432, "cell_id": 54179903}, {"distance": 0.732181182320559, "cell_id": 54148657}, {"distance": 0.7442103631247078, "cell_id": 54170399}, {"distance": 0.7487969768486321, "cell_id": 34191743}, {"distance": 0.7487969768486321, "cell_id": 35077854}, {"distance": 0.7521840774454865, "cell_id": 54246667}, {"distance": 0.7640901619851383, "cell_id": 42625768}, {"distance": 0.7640901619851383, "cell_id": 42748843}, {"distance": 0.7646129549030768, "cell_id": 48280904}, {"distance": 0.7646527181299525, "cell_id": 54172016}, {"distance": 0.7651392228681514, "cell_id": 43198687}, {"distance": 0.7676610259375556, "cell_id": 54186247}, {"distance": 0.7682068996363804, "cell_id": 48537074}, {"distance": 0.7728258502966032, "cell_id": 54185718}, {"distance": 0.776029939276966, "cell_id": 36374979}, {"distance": 0.7789787599151069, "cell_id": 39729014}, {"distance": 0.7824193883156972, "cell_id": 43190102}, {"distance": 0.7889891536753498, "cell_id": 54341661}, {"distance": 0.7901072651884996, "cell_id": 54205898}]}
{"cell_id": 50925377, "neighbors": [{"distance": 0.0, "cell_id": 50925377}, {"distance": 0.8102789787877616, "cell_id": 46676049}, {"distance": 0.8725316969621637, "cell_id": 37147644}, {"distance": 0.8804761854918758, "cell_id": 8144627}, {"distance": 0.8896732454365659, "cell_id": 37241385}, {"distance": 0.9444021330571174, "cell_id": 37241629}, {"distance": 0.9478117223022737, "cell_id": 35773477}, {"distance": 0.9478117223022737, "cell_id": 36112041}, {"distance": 0.9633258720098636, "cell_id": 36621610}, {"distance": 0.9771731407037374, "cell_id": 37191068}, {"distance": 0.9786214196277997, "cell_id": 50920501}, {"distance": 0.9797288258551252, "cell_id": 37228287}, {"distance": 0.9838083344160935, "cell_id": 19487254}, {"distance": 0.9870097518300361, "cell_id": 37136940}, {"distance": 0.9877711895812153, "cell_id": 44665386}, {"distance": 0.9983874691237056, "cell_id": 37235939}, {"distance": 1.0063540238140556, "cell_id": 37143245}, {"distance": 1.0122235761964478, "cell_id": 40773960}, {"distance": 1.0350769437154923, "cell_id": 35782675}, {"distance": 1.0350769437154923, "cell_id": 36121239}, {"distance": 1.0380189703490637, "cell_id": 49777017}, {"distance": 1.0540378672178927, "cell_id": 37141340}, {"distance": 1.0568020641310703, "cell_id": 46565056}, {"distance": 1.0609501779707267, "cell_id": 35773367}, {"distance": 1.0609501779707267, "cell_id": 36111931}]}
{"cell_id": 3084688, "neighbors": [{"distance": 0.0, "cell_id": 3084688}, {"distance": 0.7004110180981762, "cell_id": 3045107}, {"distance": 0.7209157619405852, "cell_id": 3050020}, {"distance": 0.7244874167608257, "cell_id": 3136619}, {"distance": 0.7283725965865487, "cell_id": 3072913}, {"distance": 0.7394850462879035, "cell_id": 3053411}, {"distance": 0.7426367056349401, "cell_id": 3015193}, {"distance": 0.7435231586901163, "cell_id": 3084230}, {"distance": 0.7573959361040019, "cell_id": 3089508}, {"distance": 0.7624460114018342, "cell_id": 3091946}, {"distance": 0.7674003389110673, "cell_id": 3049992}, {"distance": 0.7765415489004466, "cell_id": 3068411}, {"distance": 0.7775612968287066, "cell_id": 3083200}, {"distance": 0.7813509527829179, "cell_id": 3020181}, {"distance": 0.7879243324653826, "cell_id": 54857975}, {"distance": 0.7885444994272643, "cell_id": 3073908}, {"distance": 0.7990845129107136, "cell_id": 3015933}, {"distance": 0.8105623415059661, "cell_id": 3086576}, {"distance": 0.8120944713836518, "cell_id": 3071279}, {"distance": 0.8228635153334141, "cell_id": 3049688}, {"distance": 0.8233118505626431, "cell_id": 3015279}, {"distance": 0.8303040323096381, "cell_id": 3088788}, {"distance": 0.8322192712986768, "cell_id": 3132721}, {"distance": 0.8347996516384751, "cell_id": 3081066}, {"distance": 0.8370634997598123, "cell_id": 3073339}]}
{"cell_id": 50259815, "neighbors": [{"distance": 0.0, "cell_id": 50259815}, {"distance": 0.4665665133806972, "cell_id": 49716912}, {"distance": 0.47175960548017654, "cell_id": 42499964}, {"distance": 0.47175960548017654, "cell_id": 42840466}, {"distance": 0.4809471829443771, "cell_id": 43995824}, {"distance": 0.4865533623453324, "cell_id": 44225269}, {"distance": 0.48939609199324396, "cell_id": 54710259}, {"distance": 0.494949990752749, "cell_id": 49342277}, {"distance": 0.502590448347771, "cell_id": 50943424}, {"distance": 0.5127993501353153, "cell_id": 49246822}, {"distance": 0.5267394793981474, "cell_id": 49704118}, {"distance": 0.5280888635951967, "cell_id": 50932555}, {"distance": 0.5298905971322815, "cell_id": 50871035}, {"distance": 0.5327786068714386, "cell_id": 48467427}, {"distance": 0.5335596432556128, "cell_id": 53525132}, {"distance": 0.5358108247730187, "cell_id": 50103449}, {"distance": 0.5360831233691946, "cell_id": 49156482}, {"distance": 0.5390616845009851, "cell_id": 20901404}, {"distance": 0.5433052052385646, "cell_id": 44066172}, {"distance": 0.5483657552007595, "cell_id": 49671597}, {"distance": 0.5491996317549287, "cell_id": 50920539}, {"distance": 0.5500794812712274, "cell_id": 50799304}, {"distance": 0.5509039012381665, "cell_id": 20901918}, {"distance": 0.551881849949525, "cell_id": 44450562}, {"distance": 0.5519225863942564, "cell_id": 16092408}]}
{"cell_id": 42255869, "neighbors": [{"distance": 0.0, "cell_id": 42255869}, {"distance": 0.6211045753609936, "cell_id": 33112638}, {"distance": 0.7598779625893921, "cell_id": 30565883}, {"distance": 0.7598779625893921, "cell_id": 31199441}, {"distance": 0.7746631137785531, "cell_id": 42270997}, {"distance": 0.7772834838460019, "cell_id": 42231794}, {"distance": 0.7778702853227468, "cell_id": 42271557}, {"distance": 0.8122148058827116, "cell_id": 9305231}, {"distance": 0.8122148058827116, "cell_id": 12666939}, {"distance": 0.8367783326280193, "cell_id": 42249522}, {"distance": 0.8476033289907359, "cell_id": 42265359}, {"distance": 0.8534433713351914, "cell_id": 42268649}, {"distance": 0.856168158993878, "cell_id": 30355817}, {"distance": 0.856168158993878, "cell_id": 31131217}, {"distance": 0.8574751750653401, "cell_id": 42250696}, {"distance": 0.8614885298331174, "cell_id": 42298862}, {"distance": 0.8707161943273007, "cell_id": 30355371}, {"distance": 0.8707161943273007, "cell_id": 31130771}, {"distance": 0.8812699414467022, "cell_id": 42255695}, {"distance": 0.8870692837122129, "cell_id": 42271290}, {"distance": 0.889481745697294, "cell_id": 42255290}, {"distance": 0.8961367978071566, "cell_id": 42255497}, {"distance": 0.9065199466702347, "cell_id": 42226055}, {"distance": 0.9098096327807442, "cell_id": 9303787}, {"distance": 0.9098096327807442, "cell_id": 12665495}]}
{"cell_id": 23031496, "neighbors": [{"distance": 0.0, "cell_id": 23031496}, {"distance": 0.7701697667351343, "cell_id": 23032661}, {"distance": 0.939318873987308, "cell_id": 23032024}, {"distance": 1.0283950189372664, "cell_id": 23032130}, {"distance": 1.0506817926835303, "cell_id": 23031259}, {"distance": 1.0515165389695211, "cell_id": 23032925}, {"distance": 1.0805400202158109, "cell_id": 23031012}, {"distance": 1.0902227849365007, "cell_id": 23031672}, {"distance": 1.0930595828645366, "cell_id": 23031823}, {"distance": 1.1152254932177157, "cell_id": 23032404}, {"distance": 1.1566972737787893, "cell_id": 23037053}, {"distance": 1.2017185947020752, "cell_id": 23033423}, {"distance": 1.202108283266181, "cell_id": 23031250}, {"distance": 1.212471755362381, "cell_id": 23035392}, {"distance": 1.2134093202549279, "cell_id": 23031526}, {"distance": 1.2251893385306178, "cell_id": 23032873}, {"distance": 1.2282231330893842, "cell_id": 23032914}, {"distance": 1.2343098205854557, "cell_id": 23032365}, {"distance": 1.2357748354372158, "cell_id": 23031756}, {"distance": 1.2456350904517675, "cell_id": 23032658}, {"distance": 1.2513380609578546, "cell_id": 36269617}, {"distance": 1.253091242844294, "cell_id": 23032138}, {"distance": 1.2587268991231146, "cell_id": 23035313}, {"distance": 1.2628407312490992, "cell_id": 23032921}, {"distance": 1.2658010921113947, "cell_id": 23033364}]}
{"cell_id": 58670169, "neighbors": [{"distance": 0.0, "cell_id": 58670169}, {"distance": 0.0, "cell_id": 59889126}, {"distance": 0.5195104379169908, "cell_id": 58630410}, {"distance": 0.5745081948526636, "cell_id": 58775450}, {"distance": 0.628179654321301, "cell_id": 58779138}, {"distance": 0.6441330992312582, "cell_id": 58670070}, {"distance": 0.6953855639611485, "cell_id": 58736904}, {"distance": 0.6959969547743052, "cell_id": 58688194}, {"distance": 0.6980793956260268, "cell_id": 58733482}, {"distance": 0.7081327142537599, "cell_id": 58653842}, {"distance": 0.7096506258576525, "cell_id": 58711054}, {"distance": 0.7108936818741253, "cell_id": 58727140}, {"distance": 0.7113185023249674, "cell_id": 58638607}, {"distance": 0.7159694417652266, "cell_id": 58633595}, {"distance": 0.7187356721768446, "cell_id": 58709483}, {"distance": 0.7187356721768446, "cell_id": 59898900}, {"distance": 0.7252637432103544, "cell_id": 58680265}, {"distance": 0.7276595238664447, "cell_id": 58615242}, {"distance": 0.7303753206934293, "cell_id": 58698427}, {"distance": 0.7370530332329078, "cell_id": 58704050}, {"distance": 0.7403738092249024, "cell_id": 58711579}, {"distance": 0.7404540390461167, "cell_id": 58602188}, {"distance": 0.7419821772487912, "cell_id": 58707165}, {"distance": 0.7419821772487912, "cell_id": 59898305}, {"distance": 0.7464845205911308, "cell_id": 58689813}]}
{"cell_id": 35197267, "neighbors": [{"distance": 0.0, "cell_id": 35197267}, {"distance": 0.5198027842015553, "cell_id": 44099935}, {"distance": 0.5489923509017703, "cell_id": 35197730}, {"distance": 0.5568397272169866, "cell_id": 35315778}, {"distance": 0.5897740855500684, "cell_id": 35199449}, {"distance": 0.597191445274552, "cell_id": 51324697}, {"distance": 0.6080782342392664, "cell_id": 44921364}, {"distance": 0.6176488256972791, "cell_id": 44104682}, {"distance": 0.6198292137833155, "cell_id": 17239549}, {"distance": 0.6198292137833155, "cell_id": 17294081}, {"distance": 0.6198778827377847, "cell_id": 44327135}, {"distance": 0.6227627669795563, "cell_id": 49194800}, {"distance": 0.6324143487659527, "cell_id": 35844217}, {"distance": 0.6328770571758433, "cell_id": 15666491}, {"distance": 0.644878126602745, "cell_id": 16102971}, {"distance": 0.652758243251866, "cell_id": 33667099}, {"distance": 0.652758243251866, "cell_id": 34559395}, {"distance": 0.6590647852563508, "cell_id": 44060255}, {"distance": 0.6604888096397385, "cell_id": 48519451}, {"distance": 0.6616987707746845, "cell_id": 44026421}, {"distance": 0.6686269054047427, "cell_id": 44386358}, {"distance": 0.6692690364488023, "cell_id": 33653676}, {"distance": 0.6692690364488023, "cell_id": 34545972}, {"distance": 0.6697380665292966, "cell_id": 16746763}, {"distance": 0.6719086927535023, "cell_id": 44450575}]}
{"cell_id": 13409496, "neighbors": [{"distance": 0.0, "cell_id": 13409496}, {"distance": 0.510045869493309, "cell_id": 14090901}, {"distance": 0.635275228096588, "cell_id": 45796658}, {"distance": 0.635275228096588, "cell_id": 53238994}, {"distance": 0.6376851873800936, "cell_id": 45760254}, {"distance": 0.6376851873800936, "cell_id": 53202593}, {"distance": 0.6507549156083132, "cell_id": 14175666}, {"distance": 0.6559347330444714, "cell_id": 14657105}, {"distance": 0.6585403120176514, "cell_id": 35337342}, {"distance": 0.6812754680017791, "cell_id": 15466186}, {"distance": 0.6864043240797892, "cell_id": 13661393}, {"distance": 0.6891794014772352, "cell_id": 13498332}, {"distance": 0.6939239883191547, "cell_id": 14157498}, {"distance": 0.6990654012299233, "cell_id": 14688445}, {"distance": 0.7011192440962893, "cell_id": 44687159}, {"distance": 0.7031724871922469, "cell_id": 35337610}, {"distance": 0.719727147034362, "cell_id": 15635567}, {"distance": 0.7222417722494934, "cell_id": 14022257}, {"distance": 0.7306575395245148, "cell_id": 14154874}, {"distance": 0.7379577776574411, "cell_id": 13682610}, {"distance": 0.7382308628965903, "cell_id": 35353476}, {"distance": 0.7419769652783939, "cell_id": 14689715}, {"distance": 0.7423341028306105, "cell_id": 14228775}, {"distance": 0.7460670313886428, "cell_id": 14328884}, {"distance": 0.7501700108047921, "cell_id": 15563651}]}
{"cell_id": 48900700, "neighbors": [{"distance": 0.0, "cell_id": 48900700}, {"distance": 0.5195598583533968, "cell_id": 48901762}, {"distance": 0.6405607949127488, "cell_id": 39840193}, {"distance": 0.678020526674776, "cell_id": 48633989}, {"distance": 0.7165814325200984, "cell_id": 48901804}, {"distance": 0.7236669022483474, "cell_id": 48901325}, {"distance": 0.7258579785907489, "cell_id": 62618895}, {"distance": 0.7425526663695166, "cell_id": 48903607}, {"distance": 0.7452216596881912, "cell_id": 39257942}, {"distance": 0.7452216596881912, "cell_id": 39506549}, {"distance": 0.7471917385929299, "cell_id": 39315642}, {"distance": 0.7471917385929299, "cell_id": 39567666}, {"distance": 0.7549166918029145, "cell_id": 39828310}, {"distance": 0.7597668516935967, "cell_id": 48901101}, {"distance": 0.7622777843679364, "cell_id": 39826141}, {"distance": 0.7635720904642186, "cell_id": 48901730}, {"distance": 0.7663821926883784, "cell_id": 39829961}, {"distance": 0.7672460931149864, "cell_id": 48901926}, {"distance": 0.7714480708101916, "cell_id": 48900568}, {"distance": 0.7742748781460069, "cell_id": 48901106}, {"distance": 0.7772777812984732, "cell_id": 48901946}, {"distance": 0.777833054503658, "cell_id": 48901162}, {"distance": 0.7788236459958885, "cell_id": 39324828}, {"distance": 0.7788236459958885, "cell_id": 39511459}, {"distance": 0.7817598120270497, "cell_id": 30556022}]}
{"cell_id": 39323340, "neighbors": [{"distance": 0.0, "cell_id": 39323340}, {"distance": 0.0, "cell_id": 39511101}, {"distance": 0.3005629533417228, "cell_id": 30438221}, {"distance": 0.3005629533417228, "cell_id": 31414519}, {"distance": 0.3005629533417228, "cell_id": 31915579}, {"distance": 0.3005629533417228, "cell_id": 32380912}, {"distance": 0.32438595877982856, "cell_id": 39329150}, {"distance": 0.32438595877982856, "cell_id": 39512533}, {"distance": 0.4101281341092939, "cell_id": 39329674}, {"distance": 0.42967780597092486, "cell_id": 30437598}, {"distance": 0.42967780597092486, "cell_id": 31413940}, {"distance": 0.42967780597092486, "cell_id": 31914990}, {"distance": 0.42967780597092486, "cell_id": 32380326}, {"distance": 0.43027609958143526, "cell_id": 39318045}, {"distance": 0.43027609958143526, "cell_id": 39510542}, {"distance": 0.4404446411974253, "cell_id": 39327908}, {"distance": 0.4404446411974253, "cell_id": 39512227}, {"distance": 0.4600006965059356, "cell_id": 30436613}, {"distance": 0.4600006965059356, "cell_id": 31412986}, {"distance": 0.4600006965059356, "cell_id": 31710598}, {"distance": 0.4600006965059356, "cell_id": 31914029}, {"distance": 0.4600006965059356, "cell_id": 32379368}, {"distance": 0.4688550578095917, "cell_id": 30919755}, {"distance": 0.4688550578095917, "cell_id": 31502980}, {"distance": 0.4688550578095917, "cell_id": 31725803}]}
{"cell_id": 37084953, "neighbors": [{"distance": 0.0, "cell_id": 37010884}, {"distance": 0.0, "cell_id": 37084953}, {"distance": 0.6359168140070552, "cell_id": 37019584}, {"distance": 0.6359168140070552, "cell_id": 37093556}, {"distance": 0.7122135977773104, "cell_id": 37019190}, {"distance": 0.7122135977773104, "cell_id": 37093173}, {"distance": 0.7504135368517487, "cell_id": 37020569}, {"distance": 0.7504135368517487, "cell_id": 37094512}, {"distance": 0.7678247187128978, "cell_id": 37019168}, {"distance": 0.7678247187128978, "cell_id": 37093151}, {"distance": 0.7794252093134094, "cell_id": 37017656}, {"distance": 0.7794252093134094, "cell_id": 37091685}, {"distance": 0.7835250403193181, "cell_id": 37010244}, {"distance": 0.7835250403193181, "cell_id": 37084313}, {"distance": 0.7946705374984682, "cell_id": 37020314}, {"distance": 0.7946705374984682, "cell_id": 37094263}, {"distance": 0.8080659293522148, "cell_id": 36990110}, {"distance": 0.8080659293522148, "cell_id": 37064179}, {"distance": 0.8236093755976267, "cell_id": 37017107}, {"distance": 0.8236093755976267, "cell_id": 37091153}, {"distance": 0.8257567553555338, "cell_id": 36990590}, {"distance": 0.8257567553555338, "cell_id": 37064659}, {"distance": 0.8307403551565357, "cell_id": 37000897}, {"distance": 0.8307403551565357, "cell_id": 37074966}, {"distance": 0.8381027459380952, "cell_id": 36992821}]}
{"cell_id": 16210617, "neighbors": [{"distance": 0.0, "cell_id": 16210617}, {"distance": 0.4393396945609946, "cell_id": 16178566}, {"distance": 0.45543238933547, "cell_id": 43449079}, {"distance": 0.491873486502687, "cell_id": 51155780}, {"distance": 0.49816902585945877, "cell_id": 51161564}, {"distance": 0.5040169360703042, "cell_id": 43380077}, {"distance": 0.5189412809569013, "cell_id": 48369847}, {"distance": 0.5214193247005523, "cell_id": 16663348}, {"distance": 0.532814846908902, "cell_id": 44163204}, {"distance": 0.5410817657140069, "cell_id": 43369869}, {"distance": 0.5424447193680219, "cell_id": 43904683}, {"distance": 0.54694259495202, "cell_id": 44929363}, {"distance": 0.5530808698863597, "cell_id": 49943463}, {"distance": 0.561639524749243, "cell_id": 44150753}, {"distance": 0.5653215663275172, "cell_id": 44422089}, {"distance": 0.568019310702453, "cell_id": 46653972}, {"distance": 0.577148746778592, "cell_id": 44083046}, {"distance": 0.5786609517994973, "cell_id": 44466749}, {"distance": 0.5833834342269341, "cell_id": 43764717}, {"distance": 0.5837269186014171, "cell_id": 51205098}, {"distance": 0.5842195313738732, "cell_id": 51300278}, {"distance": 0.5853297140995162, "cell_id": 51164114}, {"distance": 0.5868829810849411, "cell_id": 51145637}, {"distance": 0.5898386744235523, "cell_id": 43563506}, {"distance": 0.5966612728731171, "cell_id": 12834140}]}
{"cell_id": 18089003, "neighbors": [{"distance": 0.0, "cell_id": 18089003}, {"distance": 0.0, "cell_id": 41278049}, {"distance": 0.6826274977146373, "cell_id": 5599504}, {"distance": 0.6826274977146373, "cell_id": 7184989}, {"distance": 0.7078568342130533, "cell_id": 18192387}, {"distance": 0.7078568342130533, "cell_id": 41283853}, {"distance": 0.721103844403711, "cell_id": 18471587}, {"distance": 0.721103844403711, "cell_id": 41300047}, {"distance": 0.7307796896690513, "cell_id": 6063926}, {"distance": 0.7307796896690513, "cell_id": 7168564}, {"distance": 0.7314857652288597, "cell_id": 17795334}, {"distance": 0.7314857652288597, "cell_id": 41277105}, {"distance": 0.7331857130747755, "cell_id": 18654242}, {"distance": 0.7331857130747755, "cell_id": 41279946}, {"distance": 0.7375886019386961, "cell_id": 18585326}, {"distance": 0.7375886019386961, "cell_id": 41282585}, {"distance": 0.7456091054521012, "cell_id": 4637272}, {"distance": 0.7485100253937493, "cell_id": 6216179}, {"distance": 0.7485100253937493, "cell_id": 7181931}, {"distance": 0.7563425338725165, "cell_id": 18018386}, {"distance": 0.7563425338725165, "cell_id": 41285903}, {"distance": 0.7635755242302971, "cell_id": 5336685}, {"distance": 0.7635755242302971, "cell_id": 7184018}, {"distance": 0.7732184113240901, "cell_id": 18214765}, {"distance": 0.7732184113240901, "cell_id": 41276721}]}
{"cell_id": 56635519, "neighbors": [{"distance": 0.0, "cell_id": 56635519}, {"distance": 0.0, "cell_id": 59387692}, {"distance": 0.23060652080341715, "cell_id": 57960344}, {"distance": 0.24386801753154894, "cell_id": 57885759}, {"distance": 0.24386801753154894, "cell_id": 59696003}, {"distance": 0.2580662568238009, "cell_id": 58062936}, {"distance": 0.2629030139447843, "cell_id": 56545137}, {"distance": 0.26408247411725916, "cell_id": 57250126}, {"distance": 0.26408247411725916, "cell_id": 59539356}, {"distance": 0.26896626752693764, "cell_id": 57758835}, {"distance": 0.27118742251194716, "cell_id": 56686217}, {"distance": 0.27200487906062776, "cell_id": 57307362}, {"distance": 0.27200487906062776, "cell_id": 59553419}, {"distance": 0.27866696624437115, "cell_id": 57417452}, {"distance": 0.27866696624437115, "cell_id": 59580406}, {"distance": 0.28383811286575794, "cell_id": 56519897}, {"distance": 0.2864090765052567, "cell_id": 57523099}, {"distance": 0.2875157241488775, "cell_id": 57508396}, {"distance": 0.294659336994404, "cell_id": 56581303}, {"distance": 0.295739142091147, "cell_id": 57346829}, {"distance": 0.2959143502054007, "cell_id": 57651045}, {"distance": 0.3013194980905566, "cell_id": 57590875}, {"distance": 0.30625114878969373, "cell_id": 56953713}, {"distance": 0.3072022397109126, "cell_id": 57671193}, {"distance": 0.3072022397109126, "cell_id": 59643212}]}
{"cell_id": 5986458, "neighbors": [{"distance": 0.0, "cell_id": 4895766}, {"distance": 0.0, "cell_id": 5986458}, {"distance": 0.5896779634230227, "cell_id": 588645}, {"distance": 0.5896779634230227, "cell_id": 2169538}, {"distance": 0.6487781256043479, "cell_id": 592045}, {"distance": 0.6487781256043479, "cell_id": 1659544}, {"distance": 0.6488541826514141, "cell_id": 4897968}, {"distance": 0.6488541826514141, "cell_id": 6026713}, {"distance": 0.6589238496588009, "cell_id": 587978}, {"distance": 0.6589238496588009, "cell_id": 2168871}, {"distance": 0.6872851143394534, "cell_id": 18481495}, {"distance": 0.6872851143394534, "cell_id": 41326403}, {"distance": 0.6889204624752643, "cell_id": 561135}, {"distance": 0.6889204624752643, "cell_id": 1803459}, {"distance": 0.7002973570795211, "cell_id": 587717}, {"distance": 0.7002973570795211, "cell_id": 2168610}, {"distance": 0.7021069404410026, "cell_id": 171402}, {"distance": 0.7021069404410026, "cell_id": 536112}, {"distance": 0.7045976367410797, "cell_id": 588818}, {"distance": 0.7045976367410797, "cell_id": 2169711}, {"distance": 0.7181881271991363, "cell_id": 587932}, {"distance": 0.7181881271991363, "cell_id": 2168825}, {"distance": 0.7214925670978747, "cell_id": 419381}, {"distance": 0.7214925670978747, "cell_id": 555079}, {"distance": 0.7314593852995386, "cell_id": 549231}]}
{"cell_id": 52064871, "neighbors": [{"distance": 0.0, "cell_id": 8599841}, {"distance": 0.0, "cell_id": 52064871}, {"distance": 0.0, "cell_id": 52496544}, {"distance": 0.8772568693748752, "cell_id": 7987439}, {"distance": 1.4161158901216995, "cell_id": 8952432}, {"distance": 1.4161158901216995, "cell_id": 12375515}, {"distance": 1.4305157958079975, "cell_id": 14894385}, {"distance": 1.5257372639575615, "cell_id": 14181927}, {"distance": 1.539206178398257, "cell_id": 18900594}, {"distance": 1.578794697760199, "cell_id": 30596039}, {"distance": 1.578794697760199, "cell_id": 31211984}, {"distance": 1.579524994215551, "cell_id": 7987396}, {"distance": 1.5887482877862082, "cell_id": 8609390}, {"distance": 1.5887482877862082, "cell_id": 52074420}, {"distance": 1.5887482877862082, "cell_id": 52499290}, {"distance": 1.620108693712842, "cell_id": 14244178}, {"distance": 1.637337315813669, "cell_id": 16820306}, {"distance": 1.6431949816354363, "cell_id": 8616020}, {"distance": 1.6431949816354363, "cell_id": 52081050}, {"distance": 1.6431949816354363, "cell_id": 52500372}, {"distance": 1.643473636299484, "cell_id": 8333898}, {"distance": 1.643473636299484, "cell_id": 51981404}, {"distance": 1.643473636299484, "cell_id": 52493246}, {"distance": 1.6799439665145453, "cell_id": 21416721}, {"distance": 1.7044119705356915, "cell_id": 13550463}]}
{"cell_id": 37446172, "neighbors": [{"distance": 0.0, "cell_id": 37446172}, {"distance": 0.49523000538780376, "cell_id": 37702097}, {"distance": 0.49523000538780376, "cell_id": 38694310}, {"distance": 0.5109076507493496, "cell_id": 37684731}, {"distance": 0.5109076507493496, "cell_id": 38642866}, {"distance": 0.5248157151860108, "cell_id": 2733307}, {"distance": 0.5286864343649069, "cell_id": 40301950}, {"distance": 0.5310129800879919, "cell_id": 4053609}, {"distance": 0.5430228380142398, "cell_id": 2773212}, {"distance": 0.5484419329922569, "cell_id": 37593545}, {"distance": 0.5484419329922569, "cell_id": 38357832}, {"distance": 0.5622986115474647, "cell_id": 37443822}, {"distance": 0.5674772642857941, "cell_id": 47226542}, {"distance": 0.5674772642857941, "cell_id": 47531194}, {"distance": 0.5685584905248262, "cell_id": 37647608}, {"distance": 0.5685584905248262, "cell_id": 38521486}, {"distance": 0.5719752175135308, "cell_id": 37275694}, {"distance": 0.5748152274856402, "cell_id": 37308656}, {"distance": 0.578456578453342, "cell_id": 37322970}, {"distance": 0.5847165391202409, "cell_id": 37264861}, {"distance": 0.5847844533895662, "cell_id": 37305293}, {"distance": 0.5887064713081598, "cell_id": 37489331}, {"distance": 0.5887064713081598, "cell_id": 38052927}, {"distance": 0.591024509863269, "cell_id": 37737231}, {"distance": 0.591024509863269, "cell_id": 38796584}]}
{"cell_id": 34645682, "neighbors": [{"distance": 0.0, "cell_id": 33753386}, {"distance": 0.0, "cell_id": 34645682}, {"distance": 0.7017146756715354, "cell_id": 46709943}, {"distance": 0.732632511648059, "cell_id": 48528434}, {"distance": 0.7664481715887002, "cell_id": 61481788}, {"distance": 0.8110620579866086, "cell_id": 44512495}, {"distance": 0.8111247895165715, "cell_id": 21202883}, {"distance": 0.8245483655402699, "cell_id": 46251821}, {"distance": 0.8422265673925369, "cell_id": 16827779}, {"distance": 0.8428521603399897, "cell_id": 62823081}, {"distance": 0.8569441392602393, "cell_id": 44460542}, {"distance": 0.8598681303394557, "cell_id": 43626643}, {"distance": 0.8615336071194224, "cell_id": 43478665}, {"distance": 0.8650526540696472, "cell_id": 49051350}, {"distance": 0.8735912742451409, "cell_id": 44608972}, {"distance": 0.8748581041366331, "cell_id": 43471753}, {"distance": 0.8767503161325088, "cell_id": 61482096}, {"distance": 0.8818404865875288, "cell_id": 46309525}, {"distance": 0.8887048693716855, "cell_id": 50939359}, {"distance": 0.8893193977191719, "cell_id": 43875483}, {"distance": 0.890330356232218, "cell_id": 8633574}, {"distance": 0.890330356232218, "cell_id": 51801601}, {"distance": 0.890330356232218, "cell_id": 52291535}, {"distance": 0.8906700003275675, "cell_id": 33728779}, {"distance": 0.8906700003275675, "cell_id": 34621075}]}
{"cell_id": 53518778, "neighbors": [{"distance": 0.0, "cell_id": 53518778}, {"distance": 0.5269604474744547, "cell_id": 16561628}, {"distance": 0.540712078216867, "cell_id": 16347734}, {"distance": 0.549569995195764, "cell_id": 48109568}, {"distance": 0.5509759488802001, "cell_id": 48109575}, {"distance": 0.5625702080379268, "cell_id": 47561964}, {"distance": 0.5625702080379268, "cell_id": 61925749}, {"distance": 0.5773356431633125, "cell_id": 53518147}, {"distance": 0.5838006731074896, "cell_id": 53915486}, {"distance": 0.5843611297166441, "cell_id": 48110404}, {"distance": 0.5883360757457692, "cell_id": 53629229}, {"distance": 0.5929931355667969, "cell_id": 54655294}, {"distance": 0.5933479212084047, "cell_id": 16347780}, {"distance": 0.5939866738139721, "cell_id": 61997467}, {"distance": 0.6010117722227635, "cell_id": 16508948}, {"distance": 0.6047430492452776, "cell_id": 48157592}, {"distance": 0.6122844492348588, "cell_id": 53894850}, {"distance": 0.6167037094310007, "cell_id": 54679693}, {"distance": 0.6223203516910496, "cell_id": 53537466}, {"distance": 0.624395356992363, "cell_id": 16528080}, {"distance": 0.6269632424042586, "cell_id": 53707651}, {"distance": 0.6302146572740464, "cell_id": 48109372}, {"distance": 0.630317506455148, "cell_id": 54359839}, {"distance": 0.6307726627878651, "cell_id": 53617387}, {"distance": 0.632842996050233, "cell_id": 48108917}]}
{"cell_id": 35966857, "neighbors": [{"distance": 0.0, "cell_id": 35628293}, {"distance": 0.0, "cell_id": 35966857}, {"distance": 0.8272174498198284, "cell_id": 35522001}, {"distance": 0.8272174498198284, "cell_id": 35860565}, {"distance": 0.8419511303786827, "cell_id": 40700195}, {"distance": 0.8419511303786827, "cell_id": 40971074}, {"distance": 0.9143757108354325, "cell_id": 35629904}, {"distance": 0.9143757108354325, "cell_id": 35968468}, {"distance": 0.9332159233592723, "cell_id": 40695806}, {"distance": 0.9332159233592723, "cell_id": 40967356}, {"distance": 0.9540981682881319, "cell_id": 35630323}, {"distance": 0.9540981682881319, "cell_id": 35968887}, {"distance": 0.9801283796145477, "cell_id": 40692687}, {"distance": 0.9801283796145477, "cell_id": 40964699}, {"distance": 1.0094364582588642, "cell_id": 35630215}, {"distance": 1.0094364582588642, "cell_id": 35968779}, {"distance": 1.0230118680493312, "cell_id": 40711072}, {"distance": 1.0230118680493312, "cell_id": 40980305}, {"distance": 1.0336240431729826, "cell_id": 22634641}, {"distance": 1.042881212580461, "cell_id": 40655258}, {"distance": 1.042881212580461, "cell_id": 40932992}, {"distance": 1.0540580901892223, "cell_id": 40695624}, {"distance": 1.0540580901892223, "cell_id": 40967201}, {"distance": 1.0693565263366824, "cell_id": 22634697}, {"distance": 1.0712763073666225, "cell_id": 22634633}]}
{"cell_id": 54154148, "neighbors": [{"distance": 0.0, "cell_id": 54154148}, {"distance": 0.4478351769608494, "cell_id": 54079253}, {"distance": 0.45929296045836243, "cell_id": 54134320}, {"distance": 0.45938080817937876, "cell_id": 54109749}, {"distance": 0.5034229940904263, "cell_id": 54111974}, {"distance": 0.5087423706480712, "cell_id": 54382416}, {"distance": 0.534873263093529, "cell_id": 54128833}, {"distance": 0.5350289541170985, "cell_id": 54134366}, {"distance": 0.5362534901468372, "cell_id": 54153011}, {"distance": 0.5420022197188572, "cell_id": 54137114}, {"distance": 0.542035157796359, "cell_id": 54526466}, {"distance": 0.5534053861655035, "cell_id": 54757469}, {"distance": 0.5586695053136892, "cell_id": 54087652}, {"distance": 0.5640045335719731, "cell_id": 54380765}, {"distance": 0.56407239525692, "cell_id": 54111652}, {"distance": 0.573673606888336, "cell_id": 54382507}, {"distance": 0.5778976841590258, "cell_id": 54798557}, {"distance": 0.5788655620579118, "cell_id": 54087213}, {"distance": 0.5839106668380172, "cell_id": 54108649}, {"distance": 0.5857584663797756, "cell_id": 54151508}, {"distance": 0.5873431164804641, "cell_id": 54114066}, {"distance": 0.5894031472718699, "cell_id": 54144845}, {"distance": 0.5920024564221235, "cell_id": 54110062}, {"distance": 0.5944656508229328, "cell_id": 54111254}, {"distance": 0.5960463596580917, "cell_id": 54115984}]}
{"cell_id": 51840785, "neighbors": [{"distance": 0.0, "cell_id": 8708799}, {"distance": 0.0, "cell_id": 51840785}, {"distance": 0.0, "cell_id": 52317447}, {"distance": 0.43328526508975523, "cell_id": 30127283}, {"distance": 0.5477685545343625, "cell_id": 47674669}, {"distance": 0.5477685545343625, "cell_id": 47734407}, {"distance": 0.5795856247184947, "cell_id": 13168369}, {"distance": 0.5795856247184947, "cell_id": 13346836}, {"distance": 0.5799205416395177, "cell_id": 48287257}, {"distance": 0.5904376769156618, "cell_id": 35466151}, {"distance": 0.5932660636941421, "cell_id": 7848150}, {"distance": 0.5987651196923862, "cell_id": 18907216}, {"distance": 0.6082035709455041, "cell_id": 13154497}, {"distance": 0.6082035709455041, "cell_id": 13339398}, {"distance": 0.6141201795057931, "cell_id": 33912057}, {"distance": 0.6141201795057931, "cell_id": 34804353}, {"distance": 0.6198025149902414, "cell_id": 13896069}, {"distance": 0.623166546495952, "cell_id": 14431555}, {"distance": 0.6363340402576686, "cell_id": 36495167}, {"distance": 0.6409424999721731, "cell_id": 3797550}, {"distance": 0.6487698944682706, "cell_id": 13168141}, {"distance": 0.6487698944682706, "cell_id": 13346608}, {"distance": 0.6638572002576121, "cell_id": 39678324}, {"distance": 0.6638770039448543, "cell_id": 61891360}, {"distance": 0.6638770039448543, "cell_id": 62281511}]}
{"cell_id": 15912605, "neighbors": [{"distance": 0.0, "cell_id": 15912605}, {"distance": 0.4766138639958217, "cell_id": 48405290}, {"distance": 0.520554197439616, "cell_id": 15785776}, {"distance": 0.5217709383464493, "cell_id": 48007209}, {"distance": 0.5365955514483315, "cell_id": 45164523}, {"distance": 0.5431941502506691, "cell_id": 45008041}, {"distance": 0.5447409165729471, "cell_id": 45144013}, {"distance": 0.5559913630308797, "cell_id": 16575266}, {"distance": 0.5562505262827276, "cell_id": 45159734}, {"distance": 0.5566106083072555, "cell_id": 15909391}, {"distance": 0.5579797416375115, "cell_id": 15918214}, {"distance": 0.5614376012105968, "cell_id": 15814875}, {"distance": 0.5615430971091727, "cell_id": 15747560}, {"distance": 0.570679058213769, "cell_id": 48393561}, {"distance": 0.5731844180095674, "cell_id": 48386855}, {"distance": 0.5742215521640215, "cell_id": 16596063}, {"distance": 0.5748008774518668, "cell_id": 46404867}, {"distance": 0.5788804809550151, "cell_id": 46359307}, {"distance": 0.5807776934418497, "cell_id": 15806455}, {"distance": 0.5815155305491535, "cell_id": 45165182}, {"distance": 0.5839469837044045, "cell_id": 15928656}, {"distance": 0.5849399725618516, "cell_id": 15788865}, {"distance": 0.5851227401481093, "cell_id": 15752814}, {"distance": 0.5890413962214553, "cell_id": 45052495}, {"distance": 0.5897295566004632, "cell_id": 15710893}]}
{"cell_id": 49658246, "neighbors": [{"distance": 0.0, "cell_id": 49658246}, {"distance": 0.4904767847829311, "cell_id": 49461551}, {"distance": 0.49596596693487766, "cell_id": 43414584}, {"distance": 0.5147908602855835, "cell_id": 44502819}, {"distance": 0.5300285241084579, "cell_id": 43808827}, {"distance": 0.5522913717274639, "cell_id": 51357414}, {"distance": 0.5568261480018297, "cell_id": 51297526}, {"distance": 0.5583456880917697, "cell_id": 50843558}, {"distance": 0.5590025960543499, "cell_id": 43600746}, {"distance": 0.5678868140654274, "cell_id": 62041199}, {"distance": 0.5755737606728137, "cell_id": 49537518}, {"distance": 0.5804249830403644, "cell_id": 43848214}, {"distance": 0.5843418213658171, "cell_id": 50266556}, {"distance": 0.5844250763107033, "cell_id": 49144887}, {"distance": 0.5859170884210482, "cell_id": 54036075}, {"distance": 0.588054880398494, "cell_id": 50297614}, {"distance": 0.5897267451659528, "cell_id": 44215717}, {"distance": 0.5905720196868144, "cell_id": 49702422}, {"distance": 0.590989632349303, "cell_id": 49567790}, {"distance": 0.592742323584074, "cell_id": 43490429}, {"distance": 0.5971804992314743, "cell_id": 54111806}, {"distance": 0.601433220653065, "cell_id": 50152136}, {"distance": 0.6036783565849106, "cell_id": 49319218}, {"distance": 0.6053040016503823, "cell_id": 53745611}, {"distance": 0.6111210543277551, "cell_id": 54517824}]}
{"cell_id": 4317164, "neighbors": [{"distance": 0.0, "cell_id": 4317164}, {"distance": 0.46702833612516403, "cell_id": 4334225}, {"distance": 0.6351511143750957, "cell_id": 4244518}, {"distance": 0.6418557927183001, "cell_id": 4341796}, {"distance": 0.6647284897398255, "cell_id": 4315197}, {"distance": 0.673365502069879, "cell_id": 4357299}, {"distance": 0.6811906409062776, "cell_id": 4343184}, {"distance": 0.6826823695716288, "cell_id": 4318031}, {"distance": 0.6839221513016241, "cell_id": 4315866}, {"distance": 0.691061004708622, "cell_id": 4317449}, {"distance": 0.692996952784642, "cell_id": 4341673}, {"distance": 0.6977512714004374, "cell_id": 4318019}, {"distance": 0.6977557177265226, "cell_id": 4342652}, {"distance": 0.7036569391378299, "cell_id": 4315354}, {"distance": 0.7069138350101347, "cell_id": 4243882}, {"distance": 0.7199990475688601, "cell_id": 4318097}, {"distance": 0.720829072366339, "cell_id": 4314671}, {"distance": 0.7245093031663518, "cell_id": 4315962}, {"distance": 0.740394745702159, "cell_id": 4314973}, {"distance": 0.7413472315427225, "cell_id": 4318658}, {"distance": 0.7460083948515464, "cell_id": 4355627}, {"distance": 0.7489964191579077, "cell_id": 4250431}, {"distance": 0.7498947559945428, "cell_id": 4355317}, {"distance": 0.7551400848190399, "cell_id": 4342283}, {"distance": 0.7625457022418821, "cell_id": 4342835}]}
{"cell_id": 44034707, "neighbors": [{"distance": 0.0, "cell_id": 44034707}, {"distance": 0.7800307715150112, "cell_id": 35194148}, {"distance": 0.8094016602771666, "cell_id": 44325392}, {"distance": 0.8112879874079513, "cell_id": 46596671}, {"distance": 0.819728528568567, "cell_id": 33830149}, {"distance": 0.819728528568567, "cell_id": 34722445}, {"distance": 0.827602415167194, "cell_id": 49793590}, {"distance": 0.8434253627529301, "cell_id": 21280593}, {"distance": 0.8481973710063165, "cell_id": 20979738}, {"distance": 0.8585721341535882, "cell_id": 21245730}, {"distance": 0.8627173717146395, "cell_id": 49206126}, {"distance": 0.864766154448192, "cell_id": 20954545}, {"distance": 0.8745840517901007, "cell_id": 22308427}, {"distance": 0.8745840517901007, "cell_id": 22340735}, {"distance": 0.8825662807576243, "cell_id": 21288836}, {"distance": 0.8922702315985815, "cell_id": 21263464}, {"distance": 0.8927489255716217, "cell_id": 21268940}, {"distance": 0.9043484559879504, "cell_id": 21273008}, {"distance": 0.9610752775644821, "cell_id": 21272254}, {"distance": 0.9612460953979909, "cell_id": 49107045}, {"distance": 0.9626487144602288, "cell_id": 46142221}, {"distance": 0.9686106071978234, "cell_id": 21244879}, {"distance": 0.969071159847528, "cell_id": 21283651}, {"distance": 0.972835245749066, "cell_id": 33551536}, {"distance": 0.972835245749066, "cell_id": 34443832}]}
{"cell_id": 58373666, "neighbors": [{"distance": 0.0, "cell_id": 58373666}, {"distance": 0.6240690131119087, "cell_id": 54910465}, {"distance": 0.6240690131119087, "cell_id": 58962922}, {"distance": 0.6332818404068007, "cell_id": 58427317}, {"distance": 0.6332818404068007, "cell_id": 59829495}, {"distance": 0.6404174347531458, "cell_id": 54945753}, {"distance": 0.65299444102401, "cell_id": 58335837}, {"distance": 0.65299444102401, "cell_id": 59806718}, {"distance": 0.6569460433667005, "cell_id": 54907390}, {"distance": 0.6670011101664637, "cell_id": 58443545}, {"distance": 0.6780490103416653, "cell_id": 58379245}, {"distance": 0.6887116492916537, "cell_id": 58450273}, {"distance": 0.6926367134581352, "cell_id": 54940131}, {"distance": 0.6926367134581352, "cell_id": 58970190}, {"distance": 0.6940022766493847, "cell_id": 58917740}, {"distance": 0.6943501312388421, "cell_id": 58439231}, {"distance": 0.6960646513278816, "cell_id": 58424452}, {"distance": 0.6960646513278816, "cell_id": 59828774}, {"distance": 0.699982149616709, "cell_id": 47100953}, {"distance": 0.699982149616709, "cell_id": 47405605}, {"distance": 0.7115533886332152, "cell_id": 58480710}, {"distance": 0.7130866366018418, "cell_id": 58535212}, {"distance": 0.720794001126334, "cell_id": 58393643}, {"distance": 0.728255863106289, "cell_id": 58385974}, {"distance": 0.7305608350789272, "cell_id": 54928949}]}
{"cell_id": 16239395, "neighbors": [{"distance": 0.0, "cell_id": 16239395}, {"distance": 0.4500096794328034, "cell_id": 16136473}, {"distance": 0.465456361731817, "cell_id": 16289379}, {"distance": 0.48901883328605855, "cell_id": 16156274}, {"distance": 0.49964053381229073, "cell_id": 16282009}, {"distance": 0.5032728845029087, "cell_id": 51282832}, {"distance": 0.5279911839310052, "cell_id": 16225961}, {"distance": 0.5341833079246434, "cell_id": 16210237}, {"distance": 0.5375350251942752, "cell_id": 16053927}, {"distance": 0.5486769081058682, "cell_id": 16247865}, {"distance": 0.549270770509114, "cell_id": 54425217}, {"distance": 0.5526639473441983, "cell_id": 15863648}, {"distance": 0.553056013480457, "cell_id": 16172481}, {"distance": 0.5552710720488763, "cell_id": 46829754}, {"distance": 0.5574070093653704, "cell_id": 49216161}, {"distance": 0.5619433655868701, "cell_id": 39665422}, {"distance": 0.5645842888680087, "cell_id": 51275140}, {"distance": 0.5666030966537565, "cell_id": 16246017}, {"distance": 0.5680996207372421, "cell_id": 16183695}, {"distance": 0.5736962363520435, "cell_id": 39729388}, {"distance": 0.575322689344505, "cell_id": 46692245}, {"distance": 0.5801539743270763, "cell_id": 15722365}, {"distance": 0.5801885079621569, "cell_id": 16069233}, {"distance": 0.5805399354862696, "cell_id": 61291832}, {"distance": 0.5877571769880094, "cell_id": 45093267}]}
{"cell_id": 46963648, "neighbors": [{"distance": 0.0, "cell_id": 46963648}, {"distance": 0.6999416575067836, "cell_id": 50098225}, {"distance": 0.7021525309132568, "cell_id": 46429084}, {"distance": 0.8384830774230889, "cell_id": 49172202}, {"distance": 0.8411579966156489, "cell_id": 49732833}, {"distance": 0.8419951817233379, "cell_id": 49176020}, {"distance": 0.8496556110292034, "cell_id": 46593741}, {"distance": 0.861402185123392, "cell_id": 50039566}, {"distance": 0.8651852721582705, "cell_id": 46796105}, {"distance": 0.8774329471653868, "cell_id": 50015824}, {"distance": 0.8922652694053707, "cell_id": 62830768}, {"distance": 0.8938805475346893, "cell_id": 49587255}, {"distance": 0.9035041898769194, "cell_id": 49162688}, {"distance": 0.9076073099240455, "cell_id": 49227766}, {"distance": 0.9080644318659475, "cell_id": 49865088}, {"distance": 0.9132082830676402, "cell_id": 46699333}, {"distance": 0.9137180708549526, "cell_id": 62011001}, {"distance": 0.9174439230368272, "cell_id": 62296481}, {"distance": 0.9178398381607951, "cell_id": 53785600}, {"distance": 0.9221594752875616, "cell_id": 46919616}, {"distance": 0.925282226804348, "cell_id": 49928434}, {"distance": 0.9257731511285741, "cell_id": 62830800}, {"distance": 0.9263358263740713, "cell_id": 50253354}, {"distance": 0.9301483908715144, "cell_id": 49600924}, {"distance": 0.9305793942477888, "cell_id": 49850906}]}
{"cell_id": 3116983, "neighbors": [{"distance": 0.0, "cell_id": 3116983}, {"distance": 0.0, "cell_id": 3145879}, {"distance": 0.34416422692309906, "cell_id": 16522338}, {"distance": 0.37863150969695664, "cell_id": 16440457}, {"distance": 0.3812450481930334, "cell_id": 15959550}, {"distance": 0.3819227814635634, "cell_id": 16566957}, {"distance": 0.388372350263063, "cell_id": 16430430}, {"distance": 0.39689741845308646, "cell_id": 16444425}, {"distance": 0.40776138332185324, "cell_id": 16585744}, {"distance": 0.4103411784343723, "cell_id": 16370801}, {"distance": 0.4111234047945965, "cell_id": 16347183}, {"distance": 0.4151597277755343, "cell_id": 16408912}, {"distance": 0.42593327177377405, "cell_id": 16449060}, {"distance": 0.43338782114415925, "cell_id": 16608979}, {"distance": 0.4336451755761541, "cell_id": 16419239}, {"distance": 0.43632108087584165, "cell_id": 16410792}, {"distance": 0.43659622156926337, "cell_id": 16083141}, {"distance": 0.4459626507084446, "cell_id": 16404186}, {"distance": 0.44860367516344163, "cell_id": 16699875}, {"distance": 0.449997609253039, "cell_id": 16348935}, {"distance": 0.45169821157314377, "cell_id": 16365622}, {"distance": 0.4535226319726998, "cell_id": 16718595}, {"distance": 0.4536960751732839, "cell_id": 16487676}, {"distance": 0.45391419173621356, "cell_id": 16414664}, {"distance": 0.4554685107155715, "cell_id": 16514518}]}
{"cell_id": 5361934, "neighbors": [{"distance": 0.0, "cell_id": 5361934}, {"distance": 0.0, "cell_id": 7000010}, {"distance": 0.4197257286452684, "cell_id": 17722971}, {"distance": 0.4197257286452684, "cell_id": 41845495}, {"distance": 0.5268411179985593, "cell_id": 18628511}, {"distance": 0.5268411179985593, "cell_id": 41841968}, {"distance": 0.5271435975525807, "cell_id": 5611836}, {"distance": 0.5271435975525807, "cell_id": 6957346}, {"distance": 0.5415456297870054, "cell_id": 18363268}, {"distance": 0.5415456297870054, "cell_id": 41767636}, {"distance": 0.5503425568051248, "cell_id": 18794895}, {"distance": 0.5503425568051248, "cell_id": 41817322}, {"distance": 0.5541950529975986, "cell_id": 18154785}, {"distance": 0.5541950529975986, "cell_id": 41774019}, {"distance": 0.5618325926309553, "cell_id": 10396134}, {"distance": 0.5618325926309553, "cell_id": 24311990}, {"distance": 0.5618325926309553, "cell_id": 24936467}, {"distance": 0.5803149555042596, "cell_id": 18103556}, {"distance": 0.5803149555042596, "cell_id": 41801283}, {"distance": 0.5983329907971183, "cell_id": 18552151}, {"distance": 0.5983329907971183, "cell_id": 41852414}, {"distance": 0.6054292777073591, "cell_id": 18616708}, {"distance": 0.6054292777073591, "cell_id": 41802313}, {"distance": 0.6066675332389525, "cell_id": 313124}, {"distance": 0.6066675332389525, "cell_id": 399553}]}
{"cell_id": 44315709, "neighbors": [{"distance": 0.0, "cell_id": 44315709}, {"distance": 0.32382232482076134, "cell_id": 50194290}, {"distance": 0.35359579652411277, "cell_id": 43757649}, {"distance": 0.36105556840142533, "cell_id": 49322022}, {"distance": 0.36979371904520403, "cell_id": 49057036}, {"distance": 0.37176888357250887, "cell_id": 49274118}, {"distance": 0.37466483647441945, "cell_id": 45018409}, {"distance": 0.3776370307228959, "cell_id": 49274254}, {"distance": 0.37910657622076666, "cell_id": 49929711}, {"distance": 0.3806772324805423, "cell_id": 49452155}, {"distance": 0.3807080426274856, "cell_id": 50185996}, {"distance": 0.38909397723032874, "cell_id": 44844650}, {"distance": 0.3908252128509743, "cell_id": 49772963}, {"distance": 0.3934773095100418, "cell_id": 54460216}, {"distance": 0.3944778058851891, "cell_id": 49203258}, {"distance": 0.3983746777440161, "cell_id": 49733982}, {"distance": 0.39976118116748954, "cell_id": 4713179}, {"distance": 0.40016610851852336, "cell_id": 49989014}, {"distance": 0.401894357224822, "cell_id": 43182567}, {"distance": 0.4031796403222918, "cell_id": 50242913}, {"distance": 0.4038282791545876, "cell_id": 49293855}, {"distance": 0.40410017613552995, "cell_id": 49754533}, {"distance": 0.404592589134686, "cell_id": 4704853}, {"distance": 0.4070293068725464, "cell_id": 49525279}, {"distance": 0.40898122733476655, "cell_id": 50008108}]}
{"cell_id": 4171108, "neighbors": [{"distance": 0.0, "cell_id": 4171108}, {"distance": 0.3876584215092317, "cell_id": 4162909}, {"distance": 0.392994821990619, "cell_id": 4168885}, {"distance": 0.4281757044593793, "cell_id": 4162936}, {"distance": 0.43141498761712843, "cell_id": 4161215}, {"distance": 0.43763909255180067, "cell_id": 4170608}, {"distance": 0.443481631136373, "cell_id": 4162063}, {"distance": 0.44525588176389724, "cell_id": 4161096}, {"distance": 0.4481236912521408, "cell_id": 4163144}, {"distance": 0.45401503904944734, "cell_id": 4168083}, {"distance": 0.4551421694065302, "cell_id": 4164674}, {"distance": 0.45638723952261673, "cell_id": 4167643}, {"distance": 0.4571851201167885, "cell_id": 4159385}, {"distance": 0.4639898978507569, "cell_id": 4157926}, {"distance": 0.47015441960794263, "cell_id": 4164582}, {"distance": 0.47192235335292215, "cell_id": 4161736}, {"distance": 0.4744699458923746, "cell_id": 4164933}, {"distance": 0.47669738257941624, "cell_id": 4158308}, {"distance": 0.47779967187612404, "cell_id": 4165421}, {"distance": 0.48183029183912285, "cell_id": 4167565}, {"distance": 0.4829712897592237, "cell_id": 4166679}, {"distance": 0.4843809816113329, "cell_id": 4164089}, {"distance": 0.4876832168560228, "cell_id": 4168433}, {"distance": 0.4909144909333049, "cell_id": 4165586}, {"distance": 0.4928681368216281, "cell_id": 4157519}]}
{"cell_id": 2868231, "neighbors": [{"distance": 0.0, "cell_id": 2868231}, {"distance": 0.32317331550681616, "cell_id": 60627466}, {"distance": 0.32317331550681616, "cell_id": 60886901}, {"distance": 0.3355778232105798, "cell_id": 2901348}, {"distance": 0.35196252813800527, "cell_id": 37267179}, {"distance": 0.35555920420632675, "cell_id": 37530679}, {"distance": 0.35555920420632675, "cell_id": 38141079}, {"distance": 0.35555920420632675, "cell_id": 60570288}, {"distance": 0.35555920420632675, "cell_id": 60731643}, {"distance": 0.3635667037937195, "cell_id": 37737212}, {"distance": 0.3635667037937195, "cell_id": 38796556}, {"distance": 0.3635667037937195, "cell_id": 60658155}, {"distance": 0.3635667037937195, "cell_id": 60978958}, {"distance": 0.3642652305801184, "cell_id": 2841460}, {"distance": 0.3653629517047247, "cell_id": 60753403}, {"distance": 0.36591900214676637, "cell_id": 37308703}, {"distance": 0.3760188851639122, "cell_id": 37537011}, {"distance": 0.3760188851639122, "cell_id": 38153523}, {"distance": 0.37708263806922976, "cell_id": 37459176}, {"distance": 0.37708263806922976, "cell_id": 37986912}, {"distance": 0.37916723739671965, "cell_id": 60581345}, {"distance": 0.38197166516744163, "cell_id": 60618108}, {"distance": 0.38197166516744163, "cell_id": 60861028}, {"distance": 0.3822981104306201, "cell_id": 37514373}, {"distance": 0.3822981104306201, "cell_id": 38097493}]}
{"cell_id": 58750419, "neighbors": [{"distance": 0.0, "cell_id": 58750419}, {"distance": 0.0, "cell_id": 59908978}, {"distance": 0.40396137215356387, "cell_id": 58711819}, {"distance": 0.4567686655548558, "cell_id": 58584230}, {"distance": 0.4567686655548558, "cell_id": 59868177}, {"distance": 0.4589838835172835, "cell_id": 58700911}, {"distance": 0.4589838835172835, "cell_id": 59896748}, {"distance": 0.4827247942008069, "cell_id": 58801160}, {"distance": 0.4827247942008069, "cell_id": 59921549}, {"distance": 0.48997378154187426, "cell_id": 58654561}, {"distance": 0.4992601087839751, "cell_id": 58775917}, {"distance": 0.5082938606595735, "cell_id": 58607300}, {"distance": 0.5144610849805823, "cell_id": 58615505}, {"distance": 0.5196712841400687, "cell_id": 58646337}, {"distance": 0.5203525177106222, "cell_id": 58690230}, {"distance": 0.5220563204295209, "cell_id": 58775123}, {"distance": 0.5258586981155171, "cell_id": 58690932}, {"distance": 0.5258586981155171, "cell_id": 59894310}, {"distance": 0.5312698105774705, "cell_id": 58770821}, {"distance": 0.5312698105774705, "cell_id": 59914037}, {"distance": 0.5327783686644778, "cell_id": 58281007}, {"distance": 0.5329956526925081, "cell_id": 58701334}, {"distance": 0.5329956526925081, "cell_id": 59896861}, {"distance": 0.5407044733358726, "cell_id": 58813002}, {"distance": 0.5465625891331506, "cell_id": 58713734}]}
{"cell_id": 55075544, "neighbors": [{"distance": 0.0, "cell_id": 55075544}, {"distance": 0.5284643428074185, "cell_id": 55056925}, {"distance": 0.5431581908415779, "cell_id": 55308351}, {"distance": 0.5431581908415779, "cell_id": 59060971}, {"distance": 0.5508158984973974, "cell_id": 55079281}, {"distance": 0.55919855659575, "cell_id": 55144972}, {"distance": 0.5613985166637556, "cell_id": 55141983}, {"distance": 0.5616101710617316, "cell_id": 55008146}, {"distance": 0.565099123299701, "cell_id": 55087441}, {"distance": 0.5986706347239505, "cell_id": 55102732}, {"distance": 0.5986706347239505, "cell_id": 59010333}, {"distance": 0.6148181612117699, "cell_id": 55035980}, {"distance": 0.6148181612117699, "cell_id": 58993932}, {"distance": 0.6329158942794391, "cell_id": 55064759}, {"distance": 0.6330670388139746, "cell_id": 55066532}, {"distance": 0.6330670388139746, "cell_id": 59001399}, {"distance": 0.634902832845149, "cell_id": 55107021}, {"distance": 0.63546567274548, "cell_id": 55128481}, {"distance": 0.6363516205966235, "cell_id": 55012703}, {"distance": 0.6363516205966235, "cell_id": 58988145}, {"distance": 0.6403125067553184, "cell_id": 55064094}, {"distance": 0.6422110594165181, "cell_id": 55225578}, {"distance": 0.6433087066001669, "cell_id": 55110767}, {"distance": 0.6469120454067767, "cell_id": 55256235}, {"distance": 0.6515664538555447, "cell_id": 55324490}]}
{"cell_id": 13856700, "neighbors": [{"distance": 0.0, "cell_id": 13856700}, {"distance": 0.7166486234071366, "cell_id": 14184123}, {"distance": 0.7257699952113488, "cell_id": 14872466}, {"distance": 0.8642485648009367, "cell_id": 44265417}, {"distance": 0.9180594029867643, "cell_id": 44565900}, {"distance": 0.9236499551514376, "cell_id": 13834958}, {"distance": 0.9284363220465486, "cell_id": 15421765}, {"distance": 0.9427200621507823, "cell_id": 15053041}, {"distance": 0.9440195995336694, "cell_id": 50927083}, {"distance": 0.9442137494851758, "cell_id": 15379369}, {"distance": 0.9500105902739892, "cell_id": 14624967}, {"distance": 0.9543230190548303, "cell_id": 54798398}, {"distance": 0.9614422767897421, "cell_id": 13841606}, {"distance": 0.9704117580507822, "cell_id": 14197926}, {"distance": 0.9740844101486165, "cell_id": 15509000}, {"distance": 0.9755574519245503, "cell_id": 14522192}, {"distance": 0.9762449369222445, "cell_id": 13506868}, {"distance": 0.9772573098074864, "cell_id": 50762020}, {"distance": 0.9781214352224111, "cell_id": 50959573}, {"distance": 0.9796720880409798, "cell_id": 54575557}, {"distance": 0.9811955811907181, "cell_id": 15047459}, {"distance": 0.9883407600149445, "cell_id": 14005506}, {"distance": 0.992070021487939, "cell_id": 14623985}, {"distance": 0.9936680936399894, "cell_id": 14510473}, {"distance": 0.9958923370541647, "cell_id": 13637013}]}
{"cell_id": 29873345, "neighbors": [{"distance": 0.0, "cell_id": 29873345}, {"distance": 0.42686908876403984, "cell_id": 30872789}, {"distance": 0.42686908876403984, "cell_id": 32201262}, {"distance": 0.42686908876403984, "cell_id": 32504313}, {"distance": 0.4276090621927132, "cell_id": 29871991}, {"distance": 0.48587694231714484, "cell_id": 30339290}, {"distance": 0.48587694231714484, "cell_id": 31852402}, {"distance": 0.48587694231714484, "cell_id": 32352330}, {"distance": 0.5509608317678277, "cell_id": 30872596}, {"distance": 0.5509608317678277, "cell_id": 32201074}, {"distance": 0.5509608317678277, "cell_id": 32504227}, {"distance": 0.5531428858727105, "cell_id": 30823082}, {"distance": 0.5531428858727105, "cell_id": 32172253}, {"distance": 0.5531428858727105, "cell_id": 32491763}, {"distance": 0.5572862821313707, "cell_id": 30194444}, {"distance": 0.5572862821313707, "cell_id": 31755246}, {"distance": 0.5572862821313707, "cell_id": 32329493}, {"distance": 0.5646630282359211, "cell_id": 29867337}, {"distance": 0.5683659217215372, "cell_id": 29861519}, {"distance": 0.5718661878525638, "cell_id": 29865231}, {"distance": 0.5868564852127212, "cell_id": 29873551}, {"distance": 0.6038501949491918, "cell_id": 29872606}, {"distance": 0.6088977568591192, "cell_id": 30757058}, {"distance": 0.6088977568591192, "cell_id": 32120208}, {"distance": 0.6088977568591192, "cell_id": 32470754}]}
{"cell_id": 16619580, "neighbors": [{"distance": 0.0, "cell_id": 16619580}, {"distance": 0.30883027334021584, "cell_id": 42377607}, {"distance": 0.30883027334021584, "cell_id": 42651711}, {"distance": 0.4423917957797046, "cell_id": 54435563}, {"distance": 0.4580786045194769, "cell_id": 48566507}, {"distance": 0.45903623273454425, "cell_id": 16618839}, {"distance": 0.46556519050485395, "cell_id": 45003215}, {"distance": 0.47682521565088165, "cell_id": 16684349}, {"distance": 0.4796474948475742, "cell_id": 16348529}, {"distance": 0.48229513209871794, "cell_id": 54443809}, {"distance": 0.49481922117411514, "cell_id": 48555092}, {"distance": 0.5044889627708168, "cell_id": 16699483}, {"distance": 0.5072129756745662, "cell_id": 16662327}, {"distance": 0.5076326060555972, "cell_id": 16573639}, {"distance": 0.5091362127073888, "cell_id": 44467682}, {"distance": 0.5131306222265775, "cell_id": 16589497}, {"distance": 0.5142871839228158, "cell_id": 45330624}, {"distance": 0.5178325260597286, "cell_id": 16692456}, {"distance": 0.5193105780861084, "cell_id": 54533018}, {"distance": 0.523122251183007, "cell_id": 61458066}, {"distance": 0.5267520484904253, "cell_id": 16661033}, {"distance": 0.5312561068129917, "cell_id": 54436297}, {"distance": 0.5349844583924048, "cell_id": 54443667}, {"distance": 0.5359092438159219, "cell_id": 45004465}, {"distance": 0.5371463809284523, "cell_id": 16672121}]}
{"cell_id": 26752022, "neighbors": [{"distance": 0.0, "cell_id": 19778407}, {"distance": 0.0, "cell_id": 26752022}, {"distance": 0.0, "cell_id": 27475353}, {"distance": 0.3383034336923114, "cell_id": 19743456}, {"distance": 0.3383034336923114, "cell_id": 26752591}, {"distance": 0.3383034336923114, "cell_id": 27459206}, {"distance": 0.3726042100677344, "cell_id": 19742505}, {"distance": 0.3726042100677344, "cell_id": 27041921}, {"distance": 0.3726042100677344, "cell_id": 27516133}, {"distance": 0.41487471632939993, "cell_id": 19732903}, {"distance": 0.41487471632939993, "cell_id": 26846870}, {"distance": 0.41487471632939993, "cell_id": 27512136}, {"distance": 0.4211106672446933, "cell_id": 19748386}, {"distance": 0.4211106672446933, "cell_id": 26248945}, {"distance": 0.4211106672446933, "cell_id": 27465105}, {"distance": 0.4496040526863822, "cell_id": 19741479}, {"distance": 0.4496040526863822, "cell_id": 26750426}, {"distance": 0.4496040526863822, "cell_id": 27516591}, {"distance": 0.46421654404543966, "cell_id": 19763201}, {"distance": 0.46421654404543966, "cell_id": 27041552}, {"distance": 0.46421654404543966, "cell_id": 27479869}, {"distance": 0.46885697647834634, "cell_id": 19767732}, {"distance": 0.46885697647834634, "cell_id": 26750810}, {"distance": 0.46885697647834634, "cell_id": 27478450}, {"distance": 0.4730471867916987, "cell_id": 19775307}]}
{"cell_id": 46456148, "neighbors": [{"distance": 0.0, "cell_id": 46456148}, {"distance": 0.42434043190923104, "cell_id": 50924760}, {"distance": 0.42653021032336563, "cell_id": 15927875}, {"distance": 0.43101702058137215, "cell_id": 45208295}, {"distance": 0.434028318575565, "cell_id": 50894597}, {"distance": 0.44922251910113237, "cell_id": 49481514}, {"distance": 0.453568694845602, "cell_id": 50099181}, {"distance": 0.45546190580729357, "cell_id": 50946896}, {"distance": 0.4564892526307414, "cell_id": 49776137}, {"distance": 0.4627424247165942, "cell_id": 49499318}, {"distance": 0.47046313574458004, "cell_id": 45280229}, {"distance": 0.47383224451438727, "cell_id": 48152297}, {"distance": 0.47543169857698975, "cell_id": 49642608}, {"distance": 0.4784149017082414, "cell_id": 49425308}, {"distance": 0.48333771841217155, "cell_id": 46628011}, {"distance": 0.4855443186598186, "cell_id": 61310874}, {"distance": 0.48559240226853695, "cell_id": 50041037}, {"distance": 0.48751762266758636, "cell_id": 53599068}, {"distance": 0.4881396165551369, "cell_id": 49322538}, {"distance": 0.4905358290726284, "cell_id": 49670618}, {"distance": 0.4930710209025863, "cell_id": 46889851}, {"distance": 0.49825317274784575, "cell_id": 61220957}, {"distance": 0.4985188239956207, "cell_id": 53579708}, {"distance": 0.5001958657365527, "cell_id": 46575056}, {"distance": 0.5049767244001178, "cell_id": 49361178}]}
{"cell_id": 40064168, "neighbors": [{"distance": 0.0, "cell_id": 40064168}, {"distance": 0.0, "cell_id": 40091365}, {"distance": 0.08202870482265975, "cell_id": 45381856}, {"distance": 0.30415687107373063, "cell_id": 30716263}, {"distance": 0.30415687107373063, "cell_id": 31251368}, {"distance": 0.44561962235553776, "cell_id": 40064656}, {"distance": 0.44561962235553776, "cell_id": 40091853}, {"distance": 0.4654023632268249, "cell_id": 45382344}, {"distance": 0.49952181446127175, "cell_id": 30716290}, {"distance": 0.49952181446127175, "cell_id": 31251395}, {"distance": 0.5215347230852251, "cell_id": 45383924}, {"distance": 0.5326127045923251, "cell_id": 40066236}, {"distance": 0.5326127045923251, "cell_id": 40093433}, {"distance": 0.6241929900831042, "cell_id": 40062299}, {"distance": 0.6241929900831042, "cell_id": 40089496}, {"distance": 0.6303689054024076, "cell_id": 45379987}, {"distance": 0.6409416380967023, "cell_id": 58470536}, {"distance": 0.6409416380967023, "cell_id": 59840317}, {"distance": 0.6605322661059395, "cell_id": 45381212}, {"distance": 0.6675170955570631, "cell_id": 40063524}, {"distance": 0.6675170955570631, "cell_id": 40090721}, {"distance": 0.6731922704433464, "cell_id": 45380038}, {"distance": 0.6829011178156872, "cell_id": 40062350}, {"distance": 0.6829011178156872, "cell_id": 40089547}, {"distance": 0.6957745019094352, "cell_id": 45379708}]}
{"cell_id": 58165710, "neighbors": [{"distance": 0.0, "cell_id": 58165710}, {"distance": 0.4806728006618061, "cell_id": 55775827}, {"distance": 0.5559324203814531, "cell_id": 55364443}, {"distance": 0.5709767593900313, "cell_id": 55449427}, {"distance": 0.580868471600996, "cell_id": 55927353}, {"distance": 0.5880687574557268, "cell_id": 55731551}, {"distance": 0.5880687574557268, "cell_id": 59164976}, {"distance": 0.6012608953673726, "cell_id": 56181606}, {"distance": 0.6067909533653608, "cell_id": 55483344}, {"distance": 0.6067909533653608, "cell_id": 59103989}, {"distance": 0.610844729768711, "cell_id": 58061575}, {"distance": 0.6118744202562694, "cell_id": 57407465}, {"distance": 0.6141783093680154, "cell_id": 55401545}, {"distance": 0.62174136688353, "cell_id": 55514257}, {"distance": 0.6240476935523013, "cell_id": 56287486}, {"distance": 0.6276634881300556, "cell_id": 56262174}, {"distance": 0.6288142875589956, "cell_id": 56875663}, {"distance": 0.6288142875589956, "cell_id": 59446987}, {"distance": 0.6325433943570552, "cell_id": 55660185}, {"distance": 0.6339541977606589, "cell_id": 3886749}, {"distance": 0.639385262324112, "cell_id": 57481119}, {"distance": 0.639385262324112, "cell_id": 59596112}, {"distance": 0.6427183158186945, "cell_id": 57816360}, {"distance": 0.6427183158186945, "cell_id": 59678922}, {"distance": 0.655602255860601, "cell_id": 55747699}]}
{"cell_id": 47425285, "neighbors": [{"distance": 0.0, "cell_id": 47120633}, {"distance": 0.0, "cell_id": 47425285}, {"distance": 0.5003177524966312, "cell_id": 47154668}, {"distance": 0.5003177524966312, "cell_id": 47459320}, {"distance": 0.5741509119839365, "cell_id": 47145201}, {"distance": 0.5741509119839365, "cell_id": 47449853}, {"distance": 0.5763000469643556, "cell_id": 47148181}, {"distance": 0.5763000469643556, "cell_id": 47452833}, {"distance": 0.5885120208095231, "cell_id": 47149071}, {"distance": 0.5885120208095231, "cell_id": 47453723}, {"distance": 0.5961338389798876, "cell_id": 47149850}, {"distance": 0.5961338389798876, "cell_id": 47454502}, {"distance": 0.6235110355544113, "cell_id": 47119451}, {"distance": 0.6235110355544113, "cell_id": 47424103}, {"distance": 0.6253152296975557, "cell_id": 47154571}, {"distance": 0.6253152296975557, "cell_id": 47459223}, {"distance": 0.6261088536297642, "cell_id": 47148477}, {"distance": 0.6261088536297642, "cell_id": 47453129}, {"distance": 0.6293647844663992, "cell_id": 47147236}, {"distance": 0.6293647844663992, "cell_id": 47451888}, {"distance": 0.6353437593010459, "cell_id": 47146102}, {"distance": 0.6353437593010459, "cell_id": 47450754}, {"distance": 0.6419867592865189, "cell_id": 47116574}, {"distance": 0.6419867592865189, "cell_id": 47421226}, {"distance": 0.6614330041528984, "cell_id": 47146001}]}
{"cell_id": 48280980, "neighbors": [{"distance": 0.0, "cell_id": 48280980}, {"distance": 0.4456995437366307, "cell_id": 54184626}, {"distance": 0.45594743264273996, "cell_id": 54250472}, {"distance": 0.46824323041118504, "cell_id": 54254652}, {"distance": 0.47538380752681525, "cell_id": 54689791}, {"distance": 0.4864847505361041, "cell_id": 42620063}, {"distance": 0.4864847505361041, "cell_id": 42734537}, {"distance": 0.5129916527694796, "cell_id": 54257043}, {"distance": 0.515032581775649, "cell_id": 50830552}, {"distance": 0.5269309221609293, "cell_id": 54159427}, {"distance": 0.539930822103502, "cell_id": 54202310}, {"distance": 0.5454522697087998, "cell_id": 54268561}, {"distance": 0.5461231093276248, "cell_id": 54702962}, {"distance": 0.5509562627573765, "cell_id": 54677106}, {"distance": 0.5558393279091275, "cell_id": 54321286}, {"distance": 0.5678057256565526, "cell_id": 51125617}, {"distance": 0.5688405376369626, "cell_id": 12414831}, {"distance": 0.5688405376369626, "cell_id": 13023692}, {"distance": 0.5768249290340235, "cell_id": 54117800}, {"distance": 0.5832678625609646, "cell_id": 54174051}, {"distance": 0.5867525185274008, "cell_id": 54262374}, {"distance": 0.5895647486089821, "cell_id": 54267491}, {"distance": 0.5943867705943173, "cell_id": 54266571}, {"distance": 0.5975290729083598, "cell_id": 54256177}, {"distance": 0.5984697228997605, "cell_id": 44991487}]}
{"cell_id": 51128384, "neighbors": [{"distance": 0.0, "cell_id": 51128384}, {"distance": 0.5330822797110618, "cell_id": 51086159}, {"distance": 0.5426499198681455, "cell_id": 47588286}, {"distance": 0.5610754280039923, "cell_id": 46506021}, {"distance": 0.5996035160108195, "cell_id": 62305334}, {"distance": 0.6115943913162039, "cell_id": 61916177}, {"distance": 0.6159542005077933, "cell_id": 51215911}, {"distance": 0.626640812226715, "cell_id": 51172042}, {"distance": 0.626764359185783, "cell_id": 61293851}, {"distance": 0.6270258366944148, "cell_id": 51297309}, {"distance": 0.6421507473718634, "cell_id": 51187057}, {"distance": 0.6472245541323476, "cell_id": 51265467}, {"distance": 0.6499746444645556, "cell_id": 51008277}, {"distance": 0.6514586232464279, "cell_id": 50833346}, {"distance": 0.6577445374580887, "cell_id": 15315203}, {"distance": 0.6609281075487169, "cell_id": 61273162}, {"distance": 0.6710929476924035, "cell_id": 51297397}, {"distance": 0.6865427602652984, "cell_id": 47576596}, {"distance": 0.6867673343510245, "cell_id": 51085284}, {"distance": 0.6916483615973952, "cell_id": 51099012}, {"distance": 0.6962877549703703, "cell_id": 52738605}, {"distance": 0.6962877549703703, "cell_id": 61846318}, {"distance": 0.7014993548047893, "cell_id": 51286062}, {"distance": 0.7019647348004947, "cell_id": 51023301}, {"distance": 0.7028882318736206, "cell_id": 52826144}]}
{"cell_id": 56299610, "neighbors": [{"distance": 0.0, "cell_id": 56299610}, {"distance": 0.22412013407684495, "cell_id": 55946670}, {"distance": 0.300888691646403, "cell_id": 55366950}, {"distance": 0.3237013234262007, "cell_id": 55523787}, {"distance": 0.3522598220138071, "cell_id": 55655512}, {"distance": 0.35955874739842647, "cell_id": 55443005}, {"distance": 0.36275712696715146, "cell_id": 56027177}, {"distance": 0.36275712696715146, "cell_id": 59237819}, {"distance": 0.3671914645906413, "cell_id": 55910385}, {"distance": 0.3715538060332994, "cell_id": 55863747}, {"distance": 0.3715538060332994, "cell_id": 59197638}, {"distance": 0.3717375043375579, "cell_id": 55770800}, {"distance": 0.3744741824153991, "cell_id": 55869275}, {"distance": 0.3744741824153991, "cell_id": 59199007}, {"distance": 0.3763809029045532, "cell_id": 55703498}, {"distance": 0.3763809029045532, "cell_id": 59158077}, {"distance": 0.37835828750125755, "cell_id": 55799413}, {"distance": 0.3813813279853934, "cell_id": 56078960}, {"distance": 0.39404201104564535, "cell_id": 55348484}, {"distance": 0.396746788944026, "cell_id": 56412225}, {"distance": 0.39931845203931304, "cell_id": 55963217}, {"distance": 0.39931845203931304, "cell_id": 59222190}, {"distance": 0.3996192012653339, "cell_id": 55524226}, {"distance": 0.4013815041507816, "cell_id": 56073185}, {"distance": 0.40175898041935887, "cell_id": 56041811}]}
{"cell_id": 54104146, "neighbors": [{"distance": 0.0, "cell_id": 54104146}, {"distance": 0.27940849646228877, "cell_id": 54093260}, {"distance": 0.38434375522393355, "cell_id": 54086746}, {"distance": 0.38920031662713267, "cell_id": 54090737}, {"distance": 0.4049608552624881, "cell_id": 54077476}, {"distance": 0.406028195916218, "cell_id": 54149887}, {"distance": 0.4083587258052707, "cell_id": 54094793}, {"distance": 0.4234785802140456, "cell_id": 54086759}, {"distance": 0.4299080394582393, "cell_id": 54108401}, {"distance": 0.4327972220654157, "cell_id": 54098734}, {"distance": 0.43472510665298286, "cell_id": 54141672}, {"distance": 0.43541149877178287, "cell_id": 54097104}, {"distance": 0.4382239740599068, "cell_id": 54081302}, {"distance": 0.4386180314040186, "cell_id": 54101601}, {"distance": 0.4506042211953963, "cell_id": 54097048}, {"distance": 0.4513352336073684, "cell_id": 54079800}, {"distance": 0.4531954588345778, "cell_id": 44959589}, {"distance": 0.45409848278662174, "cell_id": 54080991}, {"distance": 0.45578101103540114, "cell_id": 54144252}, {"distance": 0.4560989505126964, "cell_id": 54082105}, {"distance": 0.45696468608126073, "cell_id": 54092646}, {"distance": 0.45813416703910026, "cell_id": 54097897}, {"distance": 0.4601689823687653, "cell_id": 54099457}, {"distance": 0.460548884684908, "cell_id": 53685492}, {"distance": 0.4608673192075211, "cell_id": 54098379}]}
{"cell_id": 17653009, "neighbors": [{"distance": 0.0, "cell_id": 17653009}, {"distance": 0.0, "cell_id": 41980593}, {"distance": 0.40116537635598865, "cell_id": 5599131}, {"distance": 0.40116537635598865, "cell_id": 7099994}, {"distance": 0.4315947149018875, "cell_id": 821116}, {"distance": 0.4315947149018875, "cell_id": 1951028}, {"distance": 0.4952405172901337, "cell_id": 4972902}, {"distance": 0.4952405172901337, "cell_id": 6285277}, {"distance": 0.5082345656671484, "cell_id": 6209524}, {"distance": 0.5082345656671484, "cell_id": 7055984}, {"distance": 0.5141221834782705, "cell_id": 5644001}, {"distance": 0.5141221834782705, "cell_id": 7050738}, {"distance": 0.5450678313761339, "cell_id": 17728247}, {"distance": 0.5450678313761339, "cell_id": 42015915}, {"distance": 0.5616949279685202, "cell_id": 17660512}, {"distance": 0.5616949279685202, "cell_id": 41960869}, {"distance": 0.5629143440994738, "cell_id": 4464198}, {"distance": 0.5668135730022962, "cell_id": 6497690}, {"distance": 0.5668135730022962, "cell_id": 7042055}, {"distance": 0.5696860110500822, "cell_id": 60044696}, {"distance": 0.5760552820726375, "cell_id": 18339087}, {"distance": 0.5760552820726375, "cell_id": 42003453}, {"distance": 0.5778428430993064, "cell_id": 1261585}, {"distance": 0.5778428430993064, "cell_id": 2236358}, {"distance": 0.5786527469865699, "cell_id": 18560037}]}
{"cell_id": 50018097, "neighbors": [{"distance": 0.0, "cell_id": 50018097}, {"distance": 0.36754692163823643, "cell_id": 49632141}, {"distance": 0.3704397881493938, "cell_id": 50015532}, {"distance": 0.37499431925177673, "cell_id": 50824122}, {"distance": 0.37636855618038595, "cell_id": 53783740}, {"distance": 0.3793436524014941, "cell_id": 54594189}, {"distance": 0.3853050457046409, "cell_id": 50257594}, {"distance": 0.38946371170188876, "cell_id": 50091789}, {"distance": 0.39753011838969476, "cell_id": 49509917}, {"distance": 0.4000804465336348, "cell_id": 49215567}, {"distance": 0.4035365782245917, "cell_id": 49548618}, {"distance": 0.40354945224604016, "cell_id": 50145925}, {"distance": 0.4055827722624526, "cell_id": 49780389}, {"distance": 0.40566057406285055, "cell_id": 49739249}, {"distance": 0.4059219485875922, "cell_id": 49441041}, {"distance": 0.407733278273729, "cell_id": 49689855}, {"distance": 0.41439183224266635, "cell_id": 50090245}, {"distance": 0.41653368385132217, "cell_id": 50242790}, {"distance": 0.4175140832891177, "cell_id": 49471861}, {"distance": 0.42174947672133056, "cell_id": 45224464}, {"distance": 0.4223197480204262, "cell_id": 53772166}, {"distance": 0.42288765487486524, "cell_id": 50281225}, {"distance": 0.42519366473361003, "cell_id": 54522672}, {"distance": 0.425609036567199, "cell_id": 43913943}, {"distance": 0.42568287006678623, "cell_id": 49851431}]}
{"cell_id": 27293760, "neighbors": [{"distance": 0.0, "cell_id": 19902709}, {"distance": 0.0, "cell_id": 27293760}, {"distance": 0.0, "cell_id": 27854485}, {"distance": 0.2806913162075483, "cell_id": 19977186}, {"distance": 0.2806913162075483, "cell_id": 26496131}, {"distance": 0.2806913162075483, "cell_id": 27828027}, {"distance": 0.2850345750973007, "cell_id": 19928500}, {"distance": 0.2850345750973007, "cell_id": 25496666}, {"distance": 0.2850345750973007, "cell_id": 27873912}, {"distance": 0.2967248995878862, "cell_id": 20008775}, {"distance": 0.2967248995878862, "cell_id": 23132401}, {"distance": 0.2967248995878862, "cell_id": 27909967}, {"distance": 0.3126019638261545, "cell_id": 20020896}, {"distance": 0.3126019638261545, "cell_id": 25499726}, {"distance": 0.3126019638261545, "cell_id": 27650268}, {"distance": 0.3205172146595023, "cell_id": 19940936}, {"distance": 0.3205172146595023, "cell_id": 25420904}, {"distance": 0.3205172146595023, "cell_id": 27808127}, {"distance": 0.33598843657958305, "cell_id": 19932037}, {"distance": 0.33598843657958305, "cell_id": 24519638}, {"distance": 0.33598843657958305, "cell_id": 27882978}, {"distance": 0.3418852145057793, "cell_id": 20012343}, {"distance": 0.3418852145057793, "cell_id": 23132892}, {"distance": 0.3418852145057793, "cell_id": 27911667}, {"distance": 0.34190826499802446, "cell_id": 3517251}]}
{"cell_id": 35350506, "neighbors": [{"distance": 0.0, "cell_id": 35350506}, {"distance": 0.48008920461984567, "cell_id": 33885564}, {"distance": 0.48008920461984567, "cell_id": 34777860}, {"distance": 0.495308127624261, "cell_id": 13946853}, {"distance": 0.5691503065187835, "cell_id": 33408312}, {"distance": 0.5691503065187835, "cell_id": 34300608}, {"distance": 0.5921755875695267, "cell_id": 33669210}, {"distance": 0.5921755875695267, "cell_id": 34561506}, {"distance": 0.6197480546364661, "cell_id": 15253305}, {"distance": 0.6285133851189402, "cell_id": 13774473}, {"distance": 0.6546436477746054, "cell_id": 13871287}, {"distance": 0.6635906358606688, "cell_id": 33817679}, {"distance": 0.6635906358606688, "cell_id": 34709975}, {"distance": 0.6692518218427969, "cell_id": 14651542}, {"distance": 0.6812035963965825, "cell_id": 33566835}, {"distance": 0.6812035963965825, "cell_id": 34459131}, {"distance": 0.6832097288215079, "cell_id": 14563474}, {"distance": 0.6883913658289997, "cell_id": 33613953}, {"distance": 0.6883913658289997, "cell_id": 34506249}, {"distance": 0.6884897075112273, "cell_id": 14663719}, {"distance": 0.69282010741597, "cell_id": 33836980}, {"distance": 0.69282010741597, "cell_id": 34729276}, {"distance": 0.7063158977881594, "cell_id": 14445980}, {"distance": 0.7129292737262718, "cell_id": 14703236}, {"distance": 0.7156959325967142, "cell_id": 33405925}]}
{"cell_id": 51201364, "neighbors": [{"distance": 0.0, "cell_id": 51201364}, {"distance": 0.546385146694375, "cell_id": 54371338}, {"distance": 0.5825560446111693, "cell_id": 52669331}, {"distance": 0.5825560446111693, "cell_id": 61733407}, {"distance": 0.5979648526217534, "cell_id": 48168457}, {"distance": 0.6082144490723214, "cell_id": 50824513}, {"distance": 0.6168550460386721, "cell_id": 8857963}, {"distance": 0.6261872056761898, "cell_id": 52819040}, {"distance": 0.6261872056761898, "cell_id": 61798942}, {"distance": 0.6339386043632361, "cell_id": 54493182}, {"distance": 0.638125685404849, "cell_id": 54650229}, {"distance": 0.6537641782337245, "cell_id": 53655920}, {"distance": 0.6580530139779404, "cell_id": 48333170}, {"distance": 0.6677913644511533, "cell_id": 53658466}, {"distance": 0.6836230661459413, "cell_id": 54328681}, {"distance": 0.6871568953229212, "cell_id": 49349389}, {"distance": 0.6871760974876511, "cell_id": 54725082}, {"distance": 0.6919440689371278, "cell_id": 54350088}, {"distance": 0.6935046809476965, "cell_id": 53874273}, {"distance": 0.6937405025635931, "cell_id": 53869839}, {"distance": 0.6956056880683645, "cell_id": 47984551}, {"distance": 0.6992606494434854, "cell_id": 47993105}, {"distance": 0.7001578252257674, "cell_id": 48355075}, {"distance": 0.7028059244696616, "cell_id": 48465919}, {"distance": 0.7034203851383669, "cell_id": 53645990}]}
{"cell_id": 14010505, "neighbors": [{"distance": 0.0, "cell_id": 14010505}, {"distance": 0.6329976410331438, "cell_id": 13893641}, {"distance": 0.8053142758561789, "cell_id": 14390736}, {"distance": 0.8127119933895754, "cell_id": 14923308}, {"distance": 0.835858480931876, "cell_id": 13450505}, {"distance": 0.8420886490311866, "cell_id": 44734876}, {"distance": 0.8653012482323283, "cell_id": 15360542}, {"distance": 0.8906221167734372, "cell_id": 54639963}, {"distance": 0.9093304431166791, "cell_id": 14581869}, {"distance": 0.9461964187876893, "cell_id": 13795873}, {"distance": 0.9575175342483262, "cell_id": 14333790}, {"distance": 0.9685917050657841, "cell_id": 34130508}, {"distance": 0.9685917050657841, "cell_id": 35016619}, {"distance": 0.9759440183582794, "cell_id": 54638695}, {"distance": 0.981070246384781, "cell_id": 34129644}, {"distance": 0.981070246384781, "cell_id": 35015755}, {"distance": 0.9831914483000109, "cell_id": 13596566}, {"distance": 0.9863156426191163, "cell_id": 13936777}, {"distance": 0.9940797632126886, "cell_id": 44861385}, {"distance": 0.9985994654168099, "cell_id": 34130022}, {"distance": 0.9985994654168099, "cell_id": 35016133}, {"distance": 1.0100802002783884, "cell_id": 15469731}, {"distance": 1.0193066601687115, "cell_id": 14514501}, {"distance": 1.024361901723799, "cell_id": 14710053}, {"distance": 1.0303804028294177, "cell_id": 15090370}]}
{"cell_id": 37169151, "neighbors": [{"distance": 0.0, "cell_id": 37169151}, {"distance": 0.4599849691349167, "cell_id": 37151963}, {"distance": 0.46350171842676147, "cell_id": 37188118}, {"distance": 0.4658971139515102, "cell_id": 37177277}, {"distance": 0.4688305028885864, "cell_id": 37165719}, {"distance": 0.4751572030198083, "cell_id": 37160852}, {"distance": 0.48902930075452555, "cell_id": 37232743}, {"distance": 0.49882146571253344, "cell_id": 37221454}, {"distance": 0.5020534959706934, "cell_id": 37234032}, {"distance": 0.5052850528845427, "cell_id": 37166960}, {"distance": 0.5081682100339229, "cell_id": 37234257}, {"distance": 0.5114539814489655, "cell_id": 37167409}, {"distance": 0.5187830415547633, "cell_id": 37167153}, {"distance": 0.5246126394052516, "cell_id": 37190208}, {"distance": 0.5288051021858071, "cell_id": 37217522}, {"distance": 0.5341825101192426, "cell_id": 37167825}, {"distance": 0.5349698677239916, "cell_id": 37168846}, {"distance": 0.5415335403987952, "cell_id": 37197498}, {"distance": 0.5419647552646896, "cell_id": 37222700}, {"distance": 0.5420367701855109, "cell_id": 37141731}, {"distance": 0.5461842675402659, "cell_id": 37170656}, {"distance": 0.5565972914741438, "cell_id": 37169471}, {"distance": 0.5635654289460922, "cell_id": 37175657}, {"distance": 0.57104918531472, "cell_id": 37157158}, {"distance": 0.5737112005222215, "cell_id": 37234049}]}
{"cell_id": 39982196, "neighbors": [{"distance": 0.0, "cell_id": 39982196}, {"distance": 0.2865501296075379, "cell_id": 39958960}, {"distance": 0.30427292621466734, "cell_id": 39982561}, {"distance": 0.3197858176260557, "cell_id": 39982375}, {"distance": 0.32589288334679917, "cell_id": 39960367}, {"distance": 0.3340863777994668, "cell_id": 39996953}, {"distance": 0.33859355638429384, "cell_id": 39985557}, {"distance": 0.34111104277910265, "cell_id": 39963709}, {"distance": 0.34465154420640537, "cell_id": 39980438}, {"distance": 0.3498627056597131, "cell_id": 39980288}, {"distance": 0.3510941431171497, "cell_id": 39963056}, {"distance": 0.35241959994114347, "cell_id": 39960571}, {"distance": 0.35998456692055986, "cell_id": 39983267}, {"distance": 0.3674454946634386, "cell_id": 39950553}, {"distance": 0.3715461400418227, "cell_id": 39980207}, {"distance": 0.371740325049836, "cell_id": 39970970}, {"distance": 0.37454082852399406, "cell_id": 39962120}, {"distance": 0.378802832351994, "cell_id": 39985565}, {"distance": 0.38093399203252276, "cell_id": 39998160}, {"distance": 0.38511927334711477, "cell_id": 39963900}, {"distance": 0.38530156631747636, "cell_id": 39965092}, {"distance": 0.38848608755894026, "cell_id": 39964114}, {"distance": 0.3888410715904984, "cell_id": 39982500}, {"distance": 0.389799315262115, "cell_id": 39976802}, {"distance": 0.3919357425010873, "cell_id": 39963893}]}
{"cell_id": 24077801, "neighbors": [{"distance": 0.0, "cell_id": 11597098}, {"distance": 0.0, "cell_id": 24077801}, {"distance": 0.0, "cell_id": 29635699}, {"distance": 0.5975005687332477, "cell_id": 11515206}, {"distance": 0.5975005687332477, "cell_id": 24073750}, {"distance": 0.5975005687332477, "cell_id": 29687034}, {"distance": 0.6605310799624007, "cell_id": 17961414}, {"distance": 0.6605310799624007, "cell_id": 41401104}, {"distance": 0.673469201896234, "cell_id": 11547266}, {"distance": 0.673469201896234, "cell_id": 24071961}, {"distance": 0.673469201896234, "cell_id": 29668151}, {"distance": 0.6801117980299222, "cell_id": 1184826}, {"distance": 0.6801117980299222, "cell_id": 2474235}, {"distance": 0.7240121506963076, "cell_id": 11585164}, {"distance": 0.7240121506963076, "cell_id": 25448884}, {"distance": 0.7240121506963076, "cell_id": 29622633}, {"distance": 0.7893375779335525, "cell_id": 11583733}, {"distance": 0.7893375779335525, "cell_id": 23214874}, {"distance": 0.7893375779335525, "cell_id": 29621244}, {"distance": 0.7952838412350582, "cell_id": 18331891}, {"distance": 0.7952838412350582, "cell_id": 41535529}, {"distance": 0.8249956004509416, "cell_id": 11550932}, {"distance": 0.8249956004509416, "cell_id": 24072995}, {"distance": 0.8249956004509416, "cell_id": 29676837}, {"distance": 0.8357971314956492, "cell_id": 18635772}]}
{"cell_id": 50234544, "neighbors": [{"distance": 0.0, "cell_id": 50234544}, {"distance": 0.5655563880465526, "cell_id": 49142097}, {"distance": 0.5952097376915361, "cell_id": 50052391}, {"distance": 0.6164298274718124, "cell_id": 33854105}, {"distance": 0.6164298274718124, "cell_id": 34746401}, {"distance": 0.6330381133971346, "cell_id": 44493277}, {"distance": 0.6589978177904394, "cell_id": 44177154}, {"distance": 0.6607934196131711, "cell_id": 44053263}, {"distance": 0.6637304022465164, "cell_id": 43441351}, {"distance": 0.6674545127129248, "cell_id": 49721926}, {"distance": 0.673849509804518, "cell_id": 49798607}, {"distance": 0.6746212289353238, "cell_id": 44243741}, {"distance": 0.6928201449276907, "cell_id": 50260828}, {"distance": 0.6934631788898128, "cell_id": 51168046}, {"distance": 0.6945534821910823, "cell_id": 44259592}, {"distance": 0.7020500463969497, "cell_id": 49954101}, {"distance": 0.7029999926203502, "cell_id": 49423902}, {"distance": 0.7085811333323424, "cell_id": 43632104}, {"distance": 0.7097151795664314, "cell_id": 44053310}, {"distance": 0.7123358108506523, "cell_id": 50173086}, {"distance": 0.7157439669187486, "cell_id": 52663714}, {"distance": 0.7157439669187486, "cell_id": 61729401}, {"distance": 0.7160596698031082, "cell_id": 49109397}, {"distance": 0.7162906747997148, "cell_id": 44269381}, {"distance": 0.7164648610737627, "cell_id": 52725454}]}
{"cell_id": 49749755, "neighbors": [{"distance": 0.0, "cell_id": 49749755}, {"distance": 0.3380231704942982, "cell_id": 49993895}, {"distance": 0.34681492871667435, "cell_id": 44249097}, {"distance": 0.35429853440247344, "cell_id": 44332389}, {"distance": 0.3579279110209156, "cell_id": 49169858}, {"distance": 0.36692381586019757, "cell_id": 49074558}, {"distance": 0.36944780099762314, "cell_id": 51347927}, {"distance": 0.3742720681250952, "cell_id": 43919971}, {"distance": 0.37479651470898007, "cell_id": 43768050}, {"distance": 0.3848720082100776, "cell_id": 43662711}, {"distance": 0.389265721257698, "cell_id": 44331276}, {"distance": 0.3919688881046273, "cell_id": 49630071}, {"distance": 0.3940527027639847, "cell_id": 44386403}, {"distance": 0.394703372480576, "cell_id": 43886125}, {"distance": 0.3990806363182422, "cell_id": 44263943}, {"distance": 0.4029295840125807, "cell_id": 49687458}, {"distance": 0.40348313902683786, "cell_id": 43869591}, {"distance": 0.4042796679684894, "cell_id": 50263229}, {"distance": 0.40440957639488795, "cell_id": 44238741}, {"distance": 0.4053822781031607, "cell_id": 44369049}, {"distance": 0.406578821204347, "cell_id": 49056663}, {"distance": 0.4068721268707113, "cell_id": 44571726}, {"distance": 0.4078477352865377, "cell_id": 43756368}, {"distance": 0.40859722869141113, "cell_id": 43758905}, {"distance": 0.40933459723447313, "cell_id": 43314590}]}
{"cell_id": 12752827, "neighbors": [{"distance": 0.0, "cell_id": 9003437}, {"distance": 0.0, "cell_id": 12752827}, {"distance": 0.6433910875480036, "cell_id": 2593684}, {"distance": 0.6776581126603077, "cell_id": 60409070}, {"distance": 0.6787683842447221, "cell_id": 61158037}, {"distance": 0.7069044486186049, "cell_id": 33573177}, {"distance": 0.7069044486186049, "cell_id": 34465473}, {"distance": 0.7450733865551518, "cell_id": 2558401}, {"distance": 0.7714001414713244, "cell_id": 61096225}, {"distance": 0.7830438812460745, "cell_id": 60370479}, {"distance": 0.783951710875547, "cell_id": 14221222}, {"distance": 0.7969377460879057, "cell_id": 17178441}, {"distance": 0.7969377460879057, "cell_id": 17638046}, {"distance": 0.8028941757376127, "cell_id": 60107682}, {"distance": 0.8162901427397666, "cell_id": 35416550}, {"distance": 0.8181030591113335, "cell_id": 60389712}, {"distance": 0.8181030591113335, "cell_id": 61122303}, {"distance": 0.8183169636025638, "cell_id": 9003924}, {"distance": 0.8183169636025638, "cell_id": 12890046}, {"distance": 0.8271436561799906, "cell_id": 14787352}, {"distance": 0.8535604169343344, "cell_id": 14168289}, {"distance": 0.8785266294278546, "cell_id": 13654531}, {"distance": 0.8801943645829281, "cell_id": 9010496}, {"distance": 0.8801943645829281, "cell_id": 12896618}, {"distance": 0.8913986483537518, "cell_id": 16887967}]}
{"cell_id": 30516539, "neighbors": [{"distance": 0.0, "cell_id": 30516539}, {"distance": 0.0, "cell_id": 31172674}, {"distance": 0.41566914113948433, "cell_id": 30517617}, {"distance": 0.41566914113948433, "cell_id": 31173515}, {"distance": 0.4707566719648022, "cell_id": 30516071}, {"distance": 0.4707566719648022, "cell_id": 31172332}, {"distance": 0.47135937319938465, "cell_id": 30517674}, {"distance": 0.47135937319938465, "cell_id": 31173560}, {"distance": 0.5584377958408098, "cell_id": 30517185}, {"distance": 0.5584377958408098, "cell_id": 31173170}, {"distance": 0.5830515178040379, "cell_id": 30518176}, {"distance": 0.5830515178040379, "cell_id": 31173946}, {"distance": 0.5925655216579877, "cell_id": 21402125}, {"distance": 0.613907591085554, "cell_id": 32993745}, {"distance": 0.6190827429850119, "cell_id": 21402458}, {"distance": 0.6208746763365063, "cell_id": 30516458}, {"distance": 0.6208746763365063, "cell_id": 31172610}, {"distance": 0.6280227773034258, "cell_id": 21404048}, {"distance": 0.628675482460288, "cell_id": 30515199}, {"distance": 0.628675482460288, "cell_id": 31171850}, {"distance": 0.6293399582118641, "cell_id": 47629281}, {"distance": 0.6293399582118641, "cell_id": 47689019}, {"distance": 0.6293399582118641, "cell_id": 47767102}, {"distance": 0.653819184882349, "cell_id": 30514730}, {"distance": 0.653819184882349, "cell_id": 31171637}]}
{"cell_id": 47485287, "neighbors": [{"distance": 0.0, "cell_id": 47180635}, {"distance": 0.0, "cell_id": 47485287}, {"distance": 0.523709926411177, "cell_id": 38450079}, {"distance": 0.523709926411177, "cell_id": 39089954}, {"distance": 0.5453497765624706, "cell_id": 47183766}, {"distance": 0.5453497765624706, "cell_id": 47488418}, {"distance": 0.5819165324883797, "cell_id": 38752913}, {"distance": 0.5819165324883797, "cell_id": 39133232}, {"distance": 0.5907388876117149, "cell_id": 47181494}, {"distance": 0.5907388876117149, "cell_id": 47486146}, {"distance": 0.6057509739008637, "cell_id": 47182735}, {"distance": 0.6057509739008637, "cell_id": 47487387}, {"distance": 0.6078531315165722, "cell_id": 47188272}, {"distance": 0.6078531315165722, "cell_id": 47492924}, {"distance": 0.6209832187884063, "cell_id": 47114209}, {"distance": 0.6209832187884063, "cell_id": 47418861}, {"distance": 0.6372672456883507, "cell_id": 38250353}, {"distance": 0.6372672456883507, "cell_id": 39061709}, {"distance": 0.6408599313401857, "cell_id": 38741728}, {"distance": 0.6408599313401857, "cell_id": 39131236}, {"distance": 0.6513489691612362, "cell_id": 47216117}, {"distance": 0.6513489691612362, "cell_id": 47520769}, {"distance": 0.6524148842371077, "cell_id": 38009557}, {"distance": 0.6524148842371077, "cell_id": 39034765}, {"distance": 0.6557918482638597, "cell_id": 60443335}]}
{"cell_id": 3122583, "neighbors": [{"distance": 0.0, "cell_id": 3122583}, {"distance": 0.9108277854443003, "cell_id": 3016227}, {"distance": 0.9615809863558402, "cell_id": 3117132}, {"distance": 0.9655932701024997, "cell_id": 3080623}, {"distance": 0.9733684844578359, "cell_id": 3121961}, {"distance": 0.9791779943012819, "cell_id": 3072644}, {"distance": 0.992886790306095, "cell_id": 61488809}, {"distance": 1.0045515140478427, "cell_id": 3116932}, {"distance": 1.0119132059500882, "cell_id": 3017052}, {"distance": 1.015043604072504, "cell_id": 3115164}, {"distance": 1.0458843021261282, "cell_id": 3114970}, {"distance": 1.0759790014745445, "cell_id": 53972210}, {"distance": 1.1084088259735176, "cell_id": 3120449}, {"distance": 1.1102860601412083, "cell_id": 3122660}, {"distance": 1.1235962956242465, "cell_id": 3122922}, {"distance": 1.1261132924840553, "cell_id": 13803715}, {"distance": 1.129148661267022, "cell_id": 3049714}, {"distance": 1.1318531628052746, "cell_id": 3116842}, {"distance": 1.1423500700545308, "cell_id": 3117358}, {"distance": 1.1423792341083399, "cell_id": 4436187}, {"distance": 1.1473659596277488, "cell_id": 7691485}, {"distance": 1.1496833248299316, "cell_id": 3121550}, {"distance": 1.1505831044896686, "cell_id": 3117883}, {"distance": 1.1524181402014004, "cell_id": 14608088}, {"distance": 1.1531603525633463, "cell_id": 21406115}]}
{"cell_id": 34836168, "neighbors": [{"distance": 0.0, "cell_id": 33943872}, {"distance": 0.0, "cell_id": 34836168}, {"distance": 0.6180509616150808, "cell_id": 33807102}, {"distance": 0.6180509616150808, "cell_id": 34699398}, {"distance": 0.7317199930603043, "cell_id": 7575972}, {"distance": 0.7820595509973453, "cell_id": 33892857}, {"distance": 0.7820595509973453, "cell_id": 34785153}, {"distance": 0.8029348949301385, "cell_id": 62825133}, {"distance": 0.8192282639108087, "cell_id": 33883592}, {"distance": 0.8192282639108087, "cell_id": 34775888}, {"distance": 0.8549758357513663, "cell_id": 35769319}, {"distance": 0.8549758357513663, "cell_id": 36107883}, {"distance": 0.8588295310632842, "cell_id": 22863797}, {"distance": 0.8719350495527721, "cell_id": 21045400}, {"distance": 0.8771195473015725, "cell_id": 35432249}, {"distance": 0.8778593165313339, "cell_id": 21033747}, {"distance": 0.8781986862957667, "cell_id": 35488216}, {"distance": 0.8799595881386719, "cell_id": 33980804}, {"distance": 0.8799595881386719, "cell_id": 34873100}, {"distance": 0.8881795182949832, "cell_id": 14452597}, {"distance": 0.8977949118028128, "cell_id": 62429626}, {"distance": 0.9040656066720045, "cell_id": 33856849}, {"distance": 0.9040656066720045, "cell_id": 34749145}, {"distance": 0.9080985301126623, "cell_id": 33860851}, {"distance": 0.9080985301126623, "cell_id": 34753147}]}
{"cell_id": 39843725, "neighbors": [{"distance": 0.0, "cell_id": 39843725}, {"distance": 0.5425541874377785, "cell_id": 39829614}, {"distance": 0.6155405995469558, "cell_id": 39835354}, {"distance": 0.6347611124678681, "cell_id": 39838625}, {"distance": 0.6436519454480402, "cell_id": 39813710}, {"distance": 0.6520285371825211, "cell_id": 39832934}, {"distance": 0.6707472189938177, "cell_id": 39844734}, {"distance": 0.6715967399407269, "cell_id": 39860434}, {"distance": 0.7066275407263424, "cell_id": 39819758}, {"distance": 0.7094335787310105, "cell_id": 39838206}, {"distance": 0.7318203689846055, "cell_id": 39857019}, {"distance": 0.7434333664696098, "cell_id": 39861594}, {"distance": 0.7525782108452947, "cell_id": 39378481}, {"distance": 0.7687932621389604, "cell_id": 39855124}, {"distance": 0.7800102900418417, "cell_id": 39798600}, {"distance": 0.7827441236251161, "cell_id": 39859479}, {"distance": 0.7899185965757809, "cell_id": 39840381}, {"distance": 0.7916424758084617, "cell_id": 39793258}, {"distance": 0.7924431897496382, "cell_id": 39796635}, {"distance": 0.7956688619526957, "cell_id": 39856964}, {"distance": 0.80356077371784, "cell_id": 39844928}, {"distance": 0.8047276710085874, "cell_id": 39860385}, {"distance": 0.8057032255174735, "cell_id": 39857447}, {"distance": 0.8211724205319783, "cell_id": 39858922}, {"distance": 0.8295541493332211, "cell_id": 39860494}]}
{"cell_id": 48037203, "neighbors": [{"distance": 0.0, "cell_id": 48037203}, {"distance": 0.4481360430555764, "cell_id": 47990012}, {"distance": 0.4773333866332414, "cell_id": 50888392}, {"distance": 0.5877821943659385, "cell_id": 34117633}, {"distance": 0.5877821943659385, "cell_id": 35003744}, {"distance": 0.5989156165540918, "cell_id": 48534659}, {"distance": 0.59953344584726, "cell_id": 50895908}, {"distance": 0.6146538111697667, "cell_id": 35318409}, {"distance": 0.6343164055857541, "cell_id": 21996808}, {"distance": 0.6393863320274892, "cell_id": 21938563}, {"distance": 0.6422634983770035, "cell_id": 47943801}, {"distance": 0.6457949954812123, "cell_id": 48094314}, {"distance": 0.6541371692668376, "cell_id": 61547278}, {"distance": 0.6579070520417265, "cell_id": 42565339}, {"distance": 0.6579070520417265, "cell_id": 42935815}, {"distance": 0.6691819888961551, "cell_id": 16695279}, {"distance": 0.6752130336442234, "cell_id": 16661398}, {"distance": 0.6794515579873925, "cell_id": 42551354}, {"distance": 0.6794515579873925, "cell_id": 42914654}, {"distance": 0.6803508536917046, "cell_id": 16616355}, {"distance": 0.6825305206773216, "cell_id": 51055013}, {"distance": 0.6843153248740308, "cell_id": 51018794}, {"distance": 0.6843434697551768, "cell_id": 36454280}, {"distance": 0.687568905013323, "cell_id": 42543858}, {"distance": 0.687568905013323, "cell_id": 42903811}]}
{"cell_id": 55711574, "neighbors": [{"distance": 0.0, "cell_id": 55711574}, {"distance": 0.19397145057249116, "cell_id": 56227817}, {"distance": 0.21112432672457324, "cell_id": 56424837}, {"distance": 0.23261964806197527, "cell_id": 56002392}, {"distance": 0.23443214357644107, "cell_id": 56098861}, {"distance": 0.23941030507195948, "cell_id": 55734096}, {"distance": 0.23941030507195948, "cell_id": 59165550}, {"distance": 0.2661298719233359, "cell_id": 56063481}, {"distance": 0.27482578874766783, "cell_id": 55701231}, {"distance": 0.27791155112666804, "cell_id": 55793428}, {"distance": 0.2796592936303854, "cell_id": 55726082}, {"distance": 0.2838133029575541, "cell_id": 56014470}, {"distance": 0.2841938730570972, "cell_id": 56231291}, {"distance": 0.2841938730570972, "cell_id": 59288259}, {"distance": 0.28612572980855105, "cell_id": 56187409}, {"distance": 0.2875859049467252, "cell_id": 56192881}, {"distance": 0.2875859049467252, "cell_id": 59278869}, {"distance": 0.28840908951802235, "cell_id": 57265606}, {"distance": 0.29142643410497737, "cell_id": 56277061}, {"distance": 0.2923621397339766, "cell_id": 56058524}, {"distance": 0.2935465479531829, "cell_id": 55727951}, {"distance": 0.2935465479531829, "cell_id": 59164099}, {"distance": 0.29384588738243533, "cell_id": 56103433}, {"distance": 0.29384588738243533, "cell_id": 59256741}, {"distance": 0.2987187904208402, "cell_id": 55936835}]}
{"cell_id": 22694370, "neighbors": [{"distance": 0.0, "cell_id": 22694370}, {"distance": 0.557692954085267, "cell_id": 22757709}, {"distance": 0.6052461751576552, "cell_id": 22768618}, {"distance": 0.666421967787795, "cell_id": 22722545}, {"distance": 0.6817451681668132, "cell_id": 22684320}, {"distance": 0.7293446715095586, "cell_id": 22697458}, {"distance": 0.7345733346715426, "cell_id": 22743082}, {"distance": 0.7512047166433969, "cell_id": 22730523}, {"distance": 0.7530043122477853, "cell_id": 22768152}, {"distance": 0.7639339258247645, "cell_id": 39285649}, {"distance": 0.7853846288342277, "cell_id": 22687172}, {"distance": 0.8042621486550235, "cell_id": 22723131}, {"distance": 0.8292428536267364, "cell_id": 22718253}, {"distance": 0.8323579049996445, "cell_id": 22761197}, {"distance": 0.8381941933987316, "cell_id": 22762046}, {"distance": 0.8391465376922579, "cell_id": 22705873}, {"distance": 0.8521733760745478, "cell_id": 22684037}, {"distance": 0.8665162018281972, "cell_id": 22713641}, {"distance": 0.8699887447088167, "cell_id": 22737949}, {"distance": 0.8784986287891393, "cell_id": 31042473}, {"distance": 0.8784986287891393, "cell_id": 31540883}, {"distance": 0.8784986287891393, "cell_id": 32319757}, {"distance": 0.8784986287891393, "cell_id": 32559493}, {"distance": 0.879049494990862, "cell_id": 22756166}, {"distance": 0.8807332134201368, "cell_id": 22723138}]}
{"cell_id": 48445734, "neighbors": [{"distance": 0.0, "cell_id": 48445734}, {"distance": 0.14266422715523286, "cell_id": 45305816}, {"distance": 0.5103016879757348, "cell_id": 46146237}, {"distance": 0.5165437376670656, "cell_id": 46379099}, {"distance": 0.5295489243892332, "cell_id": 46448221}, {"distance": 0.5681443663038256, "cell_id": 46505409}, {"distance": 0.577503510364136, "cell_id": 46435370}, {"distance": 0.610300043256313, "cell_id": 46476881}, {"distance": 0.6259917999775867, "cell_id": 46902230}, {"distance": 0.6342038698564192, "cell_id": 48446844}, {"distance": 0.6429828846255478, "cell_id": 46130119}, {"distance": 0.6452759194584379, "cell_id": 46617159}, {"distance": 0.6472810786828043, "cell_id": 44654897}, {"distance": 0.6494842032865813, "cell_id": 45300314}, {"distance": 0.6517753721351695, "cell_id": 50761705}, {"distance": 0.6588125868370622, "cell_id": 44644650}, {"distance": 0.6621440698136403, "cell_id": 54069476}, {"distance": 0.6624490239508671, "cell_id": 46134101}, {"distance": 0.6695089732974265, "cell_id": 45307420}, {"distance": 0.6700804342030138, "cell_id": 46505908}, {"distance": 0.6708658796371132, "cell_id": 54856558}, {"distance": 0.670935229114127, "cell_id": 46313992}, {"distance": 0.6840310704577047, "cell_id": 46519424}, {"distance": 0.685593373761544, "cell_id": 46535047}, {"distance": 0.6879392268080973, "cell_id": 16106862}]}
{"cell_id": 17299529, "neighbors": [{"distance": 0.0, "cell_id": 16892184}, {"distance": 0.0, "cell_id": 17299529}, {"distance": 0.6625251000631652, "cell_id": 40008953}, {"distance": 0.6938208519421697, "cell_id": 17013534}, {"distance": 0.6938208519421697, "cell_id": 17442965}, {"distance": 0.696438237198191, "cell_id": 17009571}, {"distance": 0.696438237198191, "cell_id": 17438251}, {"distance": 0.7126493856664006, "cell_id": 17131493}, {"distance": 0.7126493856664006, "cell_id": 17582436}, {"distance": 0.7361141162939713, "cell_id": 17006631}, {"distance": 0.7361141162939713, "cell_id": 17434809}, {"distance": 0.7392250535905255, "cell_id": 8971282}, {"distance": 0.7392250535905255, "cell_id": 12720672}, {"distance": 0.7401436776039485, "cell_id": 9000927}, {"distance": 0.7401436776039485, "cell_id": 12750317}, {"distance": 0.7425103759176996, "cell_id": 17110330}, {"distance": 0.7425103759176996, "cell_id": 17557393}, {"distance": 0.7427181471407206, "cell_id": 14861345}, {"distance": 0.764238667186175, "cell_id": 16854391}, {"distance": 0.764238667186175, "cell_id": 17254762}, {"distance": 0.7819353072606996, "cell_id": 16883958}, {"distance": 0.7819353072606996, "cell_id": 17289820}, {"distance": 0.7956682992212771, "cell_id": 16939139}, {"distance": 0.7956682992212771, "cell_id": 17355109}, {"distance": 0.8086082863172516, "cell_id": 30097809}]}
{"cell_id": 62573433, "neighbors": [{"distance": 0.0, "cell_id": 62573433}, {"distance": 0.9077001397825986, "cell_id": 62647462}, {"distance": 0.9300370580540218, "cell_id": 62550398}, {"distance": 0.9575100397784202, "cell_id": 35638951}, {"distance": 0.9575100397784202, "cell_id": 35977515}, {"distance": 1.0572288515694381, "cell_id": 62347613}, {"distance": 1.068779877278541, "cell_id": 62647746}, {"distance": 1.1124334950551138, "cell_id": 14717534}, {"distance": 1.115674921170161, "cell_id": 62348527}, {"distance": 1.129007373043334, "cell_id": 54830425}, {"distance": 1.1299623422705019, "cell_id": 62661673}, {"distance": 1.1380173060747383, "cell_id": 54837845}, {"distance": 1.1465753890453332, "cell_id": 35604788}, {"distance": 1.1465753890453332, "cell_id": 35943352}, {"distance": 1.1561825572403484, "cell_id": 17613546}, {"distance": 1.1607612801558709, "cell_id": 35604752}, {"distance": 1.1607612801558709, "cell_id": 35943316}, {"distance": 1.1746461174868859, "cell_id": 62899555}, {"distance": 1.1780865492202117, "cell_id": 62549435}, {"distance": 1.178348571064403, "cell_id": 62550529}, {"distance": 1.17847963617445, "cell_id": 16664640}, {"distance": 1.185403386287579, "cell_id": 35597586}, {"distance": 1.185403386287579, "cell_id": 35936150}, {"distance": 1.1942100470781125, "cell_id": 21776360}, {"distance": 1.1942100470781125, "cell_id": 22009830}]}
{"cell_id": 23167487, "neighbors": [{"distance": 0.0, "cell_id": 10446424}, {"distance": 0.0, "cell_id": 23167487}, {"distance": 0.0, "cell_id": 24992856}, {"distance": 1.070144334189743, "cell_id": 5455740}, {"distance": 1.070144334189743, "cell_id": 6992876}, {"distance": 1.2086444526584714, "cell_id": 10446086}, {"distance": 1.2086444526584714, "cell_id": 24992815}, {"distance": 1.2086444526584714, "cell_id": 28019252}, {"distance": 1.2119927523535812, "cell_id": 6441492}, {"distance": 1.2119927523535812, "cell_id": 7028968}, {"distance": 1.2196439847592229, "cell_id": 10449107}, {"distance": 1.2196439847592229, "cell_id": 24990309}, {"distance": 1.2196439847592229, "cell_id": 27392891}, {"distance": 1.2726567333453047, "cell_id": 10469090}, {"distance": 1.2726567333453047, "cell_id": 23513738}, {"distance": 1.2726567333453047, "cell_id": 24998607}, {"distance": 1.3057979542738523, "cell_id": 10465497}, {"distance": 1.3057979542738523, "cell_id": 25001829}, {"distance": 1.3057979542738523, "cell_id": 28000379}, {"distance": 1.3396206977661, "cell_id": 18372117}, {"distance": 1.3396206977661, "cell_id": 41848842}, {"distance": 1.3402844007663792, "cell_id": 10469427}, {"distance": 1.3402844007663792, "cell_id": 24774931}, {"distance": 1.3402844007663792, "cell_id": 25000669}, {"distance": 1.3478345108887053, "cell_id": 10440362}]}
{"cell_id": 26534000, "neighbors": [{"distance": 0.0, "cell_id": 20127812}, {"distance": 0.0, "cell_id": 26534000}, {"distance": 0.0, "cell_id": 27709669}, {"distance": 0.2660717750815026, "cell_id": 20113505}, {"distance": 0.2660717750815026, "cell_id": 26759111}, {"distance": 0.2660717750815026, "cell_id": 27695386}, {"distance": 0.26672197462873976, "cell_id": 19941773}, {"distance": 0.26672197462873976, "cell_id": 26110611}, {"distance": 0.26672197462873976, "cell_id": 27806932}, {"distance": 0.3216025075699551, "cell_id": 19790476}, {"distance": 0.3216025075699551, "cell_id": 23082950}, {"distance": 0.3216025075699551, "cell_id": 27602514}, {"distance": 0.34244937590100166, "cell_id": 19804782}, {"distance": 0.34244937590100166, "cell_id": 25278307}, {"distance": 0.34244937590100166, "cell_id": 27592252}, {"distance": 0.34314654553842, "cell_id": 20101431}, {"distance": 0.34314654553842, "cell_id": 25253708}, {"distance": 0.34314654553842, "cell_id": 27741181}, {"distance": 0.3622225834926773, "cell_id": 19830600}, {"distance": 0.3622225834926773, "cell_id": 25533553}, {"distance": 0.3622225834926773, "cell_id": 27618094}, {"distance": 0.36337448621378365, "cell_id": 20031421}, {"distance": 0.36337448621378365, "cell_id": 24128693}, {"distance": 0.36337448621378365, "cell_id": 27659343}, {"distance": 0.3704427987934151, "cell_id": 20106105}]}
{"cell_id": 44659388, "neighbors": [{"distance": 0.0, "cell_id": 44659388}, {"distance": 0.5502821583111552, "cell_id": 43148333}, {"distance": 0.7032255999520876, "cell_id": 44655544}, {"distance": 0.7062250747998498, "cell_id": 49893313}, {"distance": 0.7148094666356357, "cell_id": 8887511}, {"distance": 0.7202867606393827, "cell_id": 49731363}, {"distance": 0.7500826286868536, "cell_id": 39927600}, {"distance": 0.7592829717685915, "cell_id": 49312808}, {"distance": 0.7642116946164461, "cell_id": 51181532}, {"distance": 0.7837366841354377, "cell_id": 35410150}, {"distance": 0.7891523455047735, "cell_id": 53662562}, {"distance": 0.7952506216400985, "cell_id": 50984505}, {"distance": 0.8078846957835171, "cell_id": 54046910}, {"distance": 0.8087207559839766, "cell_id": 49242767}, {"distance": 0.8122862279260427, "cell_id": 44631750}, {"distance": 0.8224488682980691, "cell_id": 13570905}, {"distance": 0.8245627293081746, "cell_id": 51278495}, {"distance": 0.8294177918834793, "cell_id": 51149969}, {"distance": 0.8310536597933251, "cell_id": 50853457}, {"distance": 0.8311139100482459, "cell_id": 49326779}, {"distance": 0.8336491383341722, "cell_id": 50063936}, {"distance": 0.8338038325048883, "cell_id": 33390027}, {"distance": 0.8338038325048883, "cell_id": 34282323}, {"distance": 0.8352146396804676, "cell_id": 43158328}, {"distance": 0.8356937996608627, "cell_id": 50768253}]}
{"cell_id": 19305601, "neighbors": [{"distance": 0.0, "cell_id": 19305601}, {"distance": 0.4840845988046707, "cell_id": 62703305}, {"distance": 0.5166404086110352, "cell_id": 33225773}, {"distance": 0.522260956371888, "cell_id": 33114815}, {"distance": 0.5258143916090464, "cell_id": 39224166}, {"distance": 0.5703620210607236, "cell_id": 19297169}, {"distance": 0.5713669851409543, "cell_id": 62760515}, {"distance": 0.5906809795572281, "cell_id": 62760669}, {"distance": 0.5945336044914911, "cell_id": 19369551}, {"distance": 0.5955004765259543, "cell_id": 19355648}, {"distance": 0.5959358299561397, "cell_id": 19355388}, {"distance": 0.6002697109534365, "cell_id": 19372456}, {"distance": 0.6102235260897606, "cell_id": 18983860}, {"distance": 0.611282508275862, "cell_id": 19369459}, {"distance": 0.6173141064836054, "cell_id": 30114200}, {"distance": 0.6246227078753855, "cell_id": 19350432}, {"distance": 0.6247059743368177, "cell_id": 32986984}, {"distance": 0.6255116451565327, "cell_id": 7716100}, {"distance": 0.6285577234002349, "cell_id": 19356474}, {"distance": 0.628764018294095, "cell_id": 62761238}, {"distance": 0.6323479808277176, "cell_id": 18969245}, {"distance": 0.6329443723151437, "cell_id": 19371142}, {"distance": 0.6351559788663683, "cell_id": 19371270}, {"distance": 0.6353059649248224, "cell_id": 19305413}, {"distance": 0.6366393112019486, "cell_id": 33109170}]}
{"cell_id": 28496906, "neighbors": [{"distance": 0.0, "cell_id": 11452800}, {"distance": 0.0, "cell_id": 28496906}, {"distance": 0.0, "cell_id": 29774030}, {"distance": 0.43121253152005296, "cell_id": 11442213}, {"distance": 0.43121253152005296, "cell_id": 26154435}, {"distance": 0.43121253152005296, "cell_id": 29784128}, {"distance": 0.5113800366188728, "cell_id": 11435106}, {"distance": 0.5113800366188728, "cell_id": 23623003}, {"distance": 0.5113800366188728, "cell_id": 29763506}, {"distance": 0.5321562781685168, "cell_id": 11444973}, {"distance": 0.5321562781685168, "cell_id": 28498222}, {"distance": 0.5321562781685168, "cell_id": 29790559}, {"distance": 0.5331849622766083, "cell_id": 11555673}, {"distance": 0.5331849622766083, "cell_id": 26418302}, {"distance": 0.5331849622766083, "cell_id": 29651387}, {"distance": 0.5497468582374582, "cell_id": 11448284}, {"distance": 0.5497468582374582, "cell_id": 27341216}, {"distance": 0.5497468582374582, "cell_id": 29776030}, {"distance": 0.5584001110432851, "cell_id": 11505954}, {"distance": 0.5584001110432851, "cell_id": 23803688}, {"distance": 0.5584001110432851, "cell_id": 29759393}, {"distance": 0.5663025379288813, "cell_id": 11432983}, {"distance": 0.5663025379288813, "cell_id": 28496727}, {"distance": 0.5663025379288813, "cell_id": 29770595}, {"distance": 0.5774471248457805, "cell_id": 11567658}]}
{"cell_id": 50744667, "neighbors": [{"distance": 0.0, "cell_id": 50744667}, {"distance": 0.5188969948522565, "cell_id": 13462806}, {"distance": 0.9224211994986825, "cell_id": 50743946}, {"distance": 0.9836853940363364, "cell_id": 40351643}, {"distance": 1.0226479877963544, "cell_id": 50743899}, {"distance": 1.0353626890184005, "cell_id": 8753937}, {"distance": 1.0353626890184005, "cell_id": 8789557}, {"distance": 1.0655075032371562, "cell_id": 36893800}, {"distance": 1.070413554998173, "cell_id": 50744021}, {"distance": 1.0928246332355074, "cell_id": 37296750}, {"distance": 1.0929946782805957, "cell_id": 36934830}, {"distance": 1.1073224033044071, "cell_id": 13959485}, {"distance": 1.1283591368452974, "cell_id": 50743961}, {"distance": 1.135116370386319, "cell_id": 33007517}, {"distance": 1.1393691910534722, "cell_id": 62150578}, {"distance": 1.1397864429957236, "cell_id": 30888879}, {"distance": 1.1397864429957236, "cell_id": 31315450}, {"distance": 1.1507158678790455, "cell_id": 3217419}, {"distance": 1.1529318058035718, "cell_id": 20691045}, {"distance": 1.153075844694609, "cell_id": 50737724}, {"distance": 1.1668954562929472, "cell_id": 50756667}, {"distance": 1.1689932505329488, "cell_id": 50744008}, {"distance": 1.1734907841793636, "cell_id": 13404331}, {"distance": 1.1740573702109505, "cell_id": 50744784}, {"distance": 1.1749764475037656, "cell_id": 13621563}]}
{"cell_id": 36517874, "neighbors": [{"distance": 0.0, "cell_id": 36517874}, {"distance": 1.1530777582064413, "cell_id": 36205883}, {"distance": 1.211874117524942, "cell_id": 54847311}, {"distance": 1.2138541744645335, "cell_id": 35348727}, {"distance": 1.214192385892637, "cell_id": 45623462}, {"distance": 1.214192385892637, "cell_id": 53065820}, {"distance": 1.2284433136928963, "cell_id": 35749571}, {"distance": 1.2284433136928963, "cell_id": 36088135}, {"distance": 1.2415779797070594, "cell_id": 36525764}, {"distance": 1.245824653389873, "cell_id": 14356693}, {"distance": 1.2545935156419026, "cell_id": 13618570}, {"distance": 1.2600715608709967, "cell_id": 36534261}, {"distance": 1.2638713650924225, "cell_id": 15635982}, {"distance": 1.267242065046192, "cell_id": 45454483}, {"distance": 1.267242065046192, "cell_id": 52896860}, {"distance": 1.2724680097324337, "cell_id": 30647935}, {"distance": 1.2724680097324337, "cell_id": 32062225}, {"distance": 1.2724680097324337, "cell_id": 32687183}, {"distance": 1.2825115536617522, "cell_id": 14165402}, {"distance": 1.2867277742743024, "cell_id": 54846777}, {"distance": 1.2952720305720666, "cell_id": 15553710}, {"distance": 1.295778491650074, "cell_id": 35637096}, {"distance": 1.295778491650074, "cell_id": 35975660}, {"distance": 1.2978979920735896, "cell_id": 36211746}, {"distance": 1.2990194840717584, "cell_id": 54847360}]}
{"cell_id": 48695688, "neighbors": [{"distance": 0.0, "cell_id": 48695688}, {"distance": 0.3791479813668783, "cell_id": 48946645}, {"distance": 0.4384995072870023, "cell_id": 48742468}, {"distance": 0.4446545650223146, "cell_id": 48599973}, {"distance": 0.4629629151719755, "cell_id": 48700720}, {"distance": 0.4890045864773436, "cell_id": 48699828}, {"distance": 0.49925806105338355, "cell_id": 48744319}, {"distance": 0.4997153020133161, "cell_id": 48894958}, {"distance": 0.5006996816530115, "cell_id": 48702577}, {"distance": 0.5100868466527533, "cell_id": 48695216}, {"distance": 0.5113411795275185, "cell_id": 48911627}, {"distance": 0.5123043248642823, "cell_id": 48703000}, {"distance": 0.5147627659687511, "cell_id": 48954185}, {"distance": 0.5158309690138877, "cell_id": 48741185}, {"distance": 0.5231453589755483, "cell_id": 48696661}, {"distance": 0.523709956394548, "cell_id": 48748232}, {"distance": 0.5335615777818543, "cell_id": 48740871}, {"distance": 0.5337717374435174, "cell_id": 48702810}, {"distance": 0.5399563990145667, "cell_id": 48939186}, {"distance": 0.5417999971223761, "cell_id": 48938825}, {"distance": 0.5419779786454277, "cell_id": 48696364}, {"distance": 0.543532314645402, "cell_id": 48747611}, {"distance": 0.5469262355830463, "cell_id": 48912456}, {"distance": 0.5526002274207652, "cell_id": 57190905}, {"distance": 0.5538494094455918, "cell_id": 48746559}]}
{"cell_id": 60813275, "neighbors": [{"distance": 0.0, "cell_id": 60452528}, {"distance": 0.0, "cell_id": 60813275}, {"distance": 0.4784588678866514, "cell_id": 14821415}, {"distance": 0.6036198832743747, "cell_id": 38588310}, {"distance": 0.6036198832743747, "cell_id": 39109596}, {"distance": 0.6475748712112381, "cell_id": 60433907}, {"distance": 0.6475748712112381, "cell_id": 60950885}, {"distance": 0.656708006418707, "cell_id": 21378710}, {"distance": 0.6786835576211582, "cell_id": 38312331}, {"distance": 0.6786835576211582, "cell_id": 39071025}, {"distance": 0.6896347292690185, "cell_id": 60885080}, {"distance": 0.6968484545627439, "cell_id": 2703582}, {"distance": 0.7011205136343381, "cell_id": 2617948}, {"distance": 0.7044681185696495, "cell_id": 60363059}, {"distance": 0.7098117006694632, "cell_id": 38857489}, {"distance": 0.7098117006694632, "cell_id": 39141335}, {"distance": 0.7200703373882384, "cell_id": 21381649}, {"distance": 0.720873662883106, "cell_id": 60455714}, {"distance": 0.720873662883106, "cell_id": 60849246}, {"distance": 0.7219850492730128, "cell_id": 61089762}, {"distance": 0.7258691764282861, "cell_id": 21381896}, {"distance": 0.7264617719089121, "cell_id": 60457937}, {"distance": 0.7277207608273721, "cell_id": 60361256}, {"distance": 0.7277207608273721, "cell_id": 61083760}, {"distance": 0.7281036118155642, "cell_id": 61125297}]}
{"cell_id": 19655252, "neighbors": [{"distance": 0.0, "cell_id": 19655252}, {"distance": 0.6408510031825089, "cell_id": 19655191}, {"distance": 0.6807057987486006, "cell_id": 19670414}, {"distance": 0.7349936437242249, "cell_id": 19646905}, {"distance": 0.7399658967810823, "cell_id": 19670402}, {"distance": 0.7459772919110166, "cell_id": 20974983}, {"distance": 0.7569016964054043, "cell_id": 19664631}, {"distance": 0.7611430466396377, "cell_id": 20966496}, {"distance": 0.7686682394462221, "cell_id": 19655319}, {"distance": 0.7817270547722077, "cell_id": 19655227}, {"distance": 0.781995553960839, "cell_id": 20966070}, {"distance": 0.7853262742394614, "cell_id": 19655440}, {"distance": 0.7879382557709795, "cell_id": 19674867}, {"distance": 0.7888212590141144, "cell_id": 19670425}, {"distance": 0.7912224010954692, "cell_id": 20961890}, {"distance": 0.7922496208864932, "cell_id": 19672328}, {"distance": 0.7924991102010166, "cell_id": 19655215}, {"distance": 0.7927104845336848, "cell_id": 19663464}, {"distance": 0.7981011144794102, "cell_id": 19667663}, {"distance": 0.8055695176766731, "cell_id": 19650135}, {"distance": 0.8059614358493479, "cell_id": 19670476}, {"distance": 0.8102115776525487, "cell_id": 19655190}, {"distance": 0.8166812845563338, "cell_id": 20966754}, {"distance": 0.8224062675680186, "cell_id": 19655643}, {"distance": 0.8287063358029574, "cell_id": 19664098}]}
{"cell_id": 51210730, "neighbors": [{"distance": 0.0, "cell_id": 51210730}, {"distance": 0.4790708278935541, "cell_id": 51088388}, {"distance": 0.6065592205177877, "cell_id": 45029131}, {"distance": 0.6084043573941836, "cell_id": 22723436}, {"distance": 0.6205580787282589, "cell_id": 46815035}, {"distance": 0.6653011530321951, "cell_id": 51187764}, {"distance": 0.6695241350867374, "cell_id": 50050487}, {"distance": 0.6774195082451118, "cell_id": 49812238}, {"distance": 0.680082661194638, "cell_id": 51303732}, {"distance": 0.6859206440600827, "cell_id": 53775296}, {"distance": 0.6901902161475535, "cell_id": 44287719}, {"distance": 0.708249016267652, "cell_id": 51093140}, {"distance": 0.7130685190418454, "cell_id": 22733697}, {"distance": 0.713610800070503, "cell_id": 44191702}, {"distance": 0.715700661542474, "cell_id": 46662188}, {"distance": 0.7171231362532456, "cell_id": 12578155}, {"distance": 0.7171231362532456, "cell_id": 12977778}, {"distance": 0.7173154130408083, "cell_id": 35489637}, {"distance": 0.7180911490243059, "cell_id": 49874876}, {"distance": 0.718205498676463, "cell_id": 45676723}, {"distance": 0.718205498676463, "cell_id": 53119074}, {"distance": 0.7189438638579081, "cell_id": 8447292}, {"distance": 0.7189438638579081, "cell_id": 51653472}, {"distance": 0.7189438638579081, "cell_id": 52216416}, {"distance": 0.7191792593746424, "cell_id": 35424644}]}
{"cell_id": 62713500, "neighbors": [{"distance": 0.0, "cell_id": 62713500}, {"distance": 0.5369514480353094, "cell_id": 62711591}, {"distance": 0.6021829964612607, "cell_id": 62707364}, {"distance": 0.6358859992769841, "cell_id": 62714277}, {"distance": 0.6365235672777158, "cell_id": 62716337}, {"distance": 0.6464027903395765, "cell_id": 62718576}, {"distance": 0.6503877864806157, "cell_id": 62717175}, {"distance": 0.6649456717494066, "cell_id": 62717322}, {"distance": 0.6687284970504356, "cell_id": 62716431}, {"distance": 0.6692722150665886, "cell_id": 62620175}, {"distance": 0.6704517311374469, "cell_id": 62718189}, {"distance": 0.6762758364050578, "cell_id": 62720261}, {"distance": 0.686745915269256, "cell_id": 62633205}, {"distance": 0.6876808555557091, "cell_id": 62713320}, {"distance": 0.6886310097872755, "cell_id": 62712191}, {"distance": 0.6889758440517073, "cell_id": 62621433}, {"distance": 0.6961374416046356, "cell_id": 62709265}, {"distance": 0.705346915490458, "cell_id": 62731456}, {"distance": 0.7093420743904721, "cell_id": 62714820}, {"distance": 0.7201313793712686, "cell_id": 62712294}, {"distance": 0.7206676569936631, "cell_id": 62703539}, {"distance": 0.7219678722822523, "cell_id": 62622155}, {"distance": 0.7262395268071373, "cell_id": 62714042}, {"distance": 0.7274575518949645, "cell_id": 62705073}, {"distance": 0.728599535237919, "cell_id": 62624913}]}
{"cell_id": 14521844, "neighbors": [{"distance": 0.0, "cell_id": 14521844}, {"distance": 0.4359932091259622, "cell_id": 45991350}, {"distance": 0.4359932091259622, "cell_id": 53433672}, {"distance": 0.6468088703291106, "cell_id": 45449897}, {"distance": 0.6468088703291106, "cell_id": 52892275}, {"distance": 0.6592875014930719, "cell_id": 14315997}, {"distance": 0.698665638280336, "cell_id": 717253}, {"distance": 0.698665638280336, "cell_id": 1847165}, {"distance": 0.7048770438857611, "cell_id": 60180183}, {"distance": 0.7048770438857611, "cell_id": 60258116}, {"distance": 0.7302952463917286, "cell_id": 15572496}, {"distance": 0.7368327844061339, "cell_id": 15223637}, {"distance": 0.7408468055680458, "cell_id": 7488183}, {"distance": 0.7408468055680458, "cell_id": 7609387}, {"distance": 0.7492252985445539, "cell_id": 15359600}, {"distance": 0.7513163241595278, "cell_id": 21114557}, {"distance": 0.7546938454237285, "cell_id": 44877170}, {"distance": 0.7676241390447742, "cell_id": 13685232}, {"distance": 0.7689299455779279, "cell_id": 14595758}, {"distance": 0.7785812919140477, "cell_id": 14583483}, {"distance": 0.7862640928381238, "cell_id": 14810117}, {"distance": 0.7866922426429679, "cell_id": 15314229}, {"distance": 0.7914585342540689, "cell_id": 60173092}, {"distance": 0.7914585342540689, "cell_id": 60246859}, {"distance": 0.8132341518097352, "cell_id": 42312387}]}
{"cell_id": 61516728, "neighbors": [{"distance": 0.0, "cell_id": 61516728}, {"distance": 0.3288066889054141, "cell_id": 4730060}, {"distance": 0.40634554050503574, "cell_id": 54140806}, {"distance": 0.4108483945169508, "cell_id": 62038398}, {"distance": 0.4392674767578177, "cell_id": 4717396}, {"distance": 0.44666229144647557, "cell_id": 33549321}, {"distance": 0.44666229144647557, "cell_id": 34441617}, {"distance": 0.4473651388256785, "cell_id": 16098469}, {"distance": 0.44771335396552747, "cell_id": 4691064}, {"distance": 0.4529140867504498, "cell_id": 4731728}, {"distance": 0.46963862092991293, "cell_id": 4693640}, {"distance": 0.4705377030597628, "cell_id": 4728793}, {"distance": 0.4739011957366877, "cell_id": 54129985}, {"distance": 0.4743273219757344, "cell_id": 4686091}, {"distance": 0.4797630370134442, "cell_id": 15966487}, {"distance": 0.48396241924212136, "cell_id": 4720867}, {"distance": 0.4845753775487582, "cell_id": 4705437}, {"distance": 0.4882510664484403, "cell_id": 43007744}, {"distance": 0.49418768244854183, "cell_id": 15971269}, {"distance": 0.4954307301427093, "cell_id": 4713038}, {"distance": 0.49735237520114034, "cell_id": 49721981}, {"distance": 0.5030580360945108, "cell_id": 16080000}, {"distance": 0.5038010562119387, "cell_id": 4710816}, {"distance": 0.5051929338241935, "cell_id": 4685910}, {"distance": 0.5129379132526222, "cell_id": 61561938}]}
{"cell_id": 61427707, "neighbors": [{"distance": 0.0, "cell_id": 61427707}, {"distance": 0.6144649207655538, "cell_id": 16445418}, {"distance": 0.646616889996576, "cell_id": 16553729}, {"distance": 0.7461048480456235, "cell_id": 16596272}, {"distance": 0.7750174222617595, "cell_id": 16634084}, {"distance": 0.7750209411754649, "cell_id": 51051692}, {"distance": 0.7820762792180753, "cell_id": 16703129}, {"distance": 0.7864632516790016, "cell_id": 16678250}, {"distance": 0.8244627892189405, "cell_id": 16528060}, {"distance": 0.8288154962008026, "cell_id": 16634393}, {"distance": 0.8419775960913605, "cell_id": 61459551}, {"distance": 0.8566124028270066, "cell_id": 16666464}, {"distance": 0.8667868011056283, "cell_id": 44982793}, {"distance": 0.871987008625847, "cell_id": 16688420}, {"distance": 0.8746874573776824, "cell_id": 50997365}, {"distance": 0.8936524896238134, "cell_id": 16703390}, {"distance": 0.8937156133375643, "cell_id": 42595476}, {"distance": 0.8937156133375643, "cell_id": 42662151}, {"distance": 0.8975506445328327, "cell_id": 16437613}, {"distance": 0.8987834193452605, "cell_id": 16673434}, {"distance": 0.9054131254342005, "cell_id": 16474232}, {"distance": 0.9096838705706686, "cell_id": 16728463}, {"distance": 0.9107304351805902, "cell_id": 50995410}, {"distance": 0.911125831603279, "cell_id": 50930272}, {"distance": 0.911745422540887, "cell_id": 16719991}]}
{"cell_id": 23813076, "neighbors": [{"distance": 0.0, "cell_id": 20262953}, {"distance": 0.0, "cell_id": 23290779}, {"distance": 0.0, "cell_id": 23813076}, {"distance": 0.4338640210459501, "cell_id": 20274708}, {"distance": 0.4338640210459501, "cell_id": 23304809}, {"distance": 0.4338640210459501, "cell_id": 28017310}, {"distance": 0.47437975949736877, "cell_id": 20260661}, {"distance": 0.47437975949736877, "cell_id": 23296234}, {"distance": 0.47437975949736877, "cell_id": 26539901}, {"distance": 0.4779407484578746, "cell_id": 20242002}, {"distance": 0.4779407484578746, "cell_id": 23343990}, {"distance": 0.4779407484578746, "cell_id": 28977650}, {"distance": 0.5175152539455917, "cell_id": 20241394}, {"distance": 0.5175152539455917, "cell_id": 23344604}, {"distance": 0.5175152539455917, "cell_id": 24533338}, {"distance": 0.5250850514398399, "cell_id": 20190844}, {"distance": 0.5250850514398399, "cell_id": 23320698}, {"distance": 0.5250850514398399, "cell_id": 29176030}, {"distance": 0.5411594976097907, "cell_id": 20258885}, {"distance": 0.5411594976097907, "cell_id": 23294812}, {"distance": 0.5411594976097907, "cell_id": 26115070}, {"distance": 0.5427233883906821, "cell_id": 20269114}, {"distance": 0.5427233883906821, "cell_id": 23293146}, {"distance": 0.5427233883906821, "cell_id": 26976237}, {"distance": 0.5624900033224448, "cell_id": 20242936}]}
{"cell_id": 5254683, "neighbors": [{"distance": 0.0, "cell_id": 5046196}, {"distance": 0.0, "cell_id": 5254683}, {"distance": 0.3799534034823139, "cell_id": 5040587}, {"distance": 0.3799534034823139, "cell_id": 5872339}, {"distance": 0.3999514837859882, "cell_id": 3336226}, {"distance": 0.3999514837859882, "cell_id": 3423140}, {"distance": 0.44074845003104507, "cell_id": 5110924}, {"distance": 0.44074845003104507, "cell_id": 6334266}, {"distance": 0.44981107325761893, "cell_id": 5122034}, {"distance": 0.44981107325761893, "cell_id": 6114451}, {"distance": 0.46992864559243125, "cell_id": 18465753}, {"distance": 0.46992864559243125, "cell_id": 41181465}, {"distance": 0.4805631046645542, "cell_id": 5141713}, {"distance": 0.4805631046645542, "cell_id": 5301465}, {"distance": 0.4879051504579435, "cell_id": 5102954}, {"distance": 0.4879051504579435, "cell_id": 5485673}, {"distance": 0.5009210904922264, "cell_id": 5081349}, {"distance": 0.5009210904922264, "cell_id": 5558699}, {"distance": 0.5021257049236438, "cell_id": 5052677}, {"distance": 0.5021257049236438, "cell_id": 6042699}, {"distance": 0.526170329019862, "cell_id": 5073632}, {"distance": 0.526170329019862, "cell_id": 6210519}, {"distance": 0.5268533612914597, "cell_id": 3335214}, {"distance": 0.5268533612914597, "cell_id": 3422968}, {"distance": 0.5291550581025147, "cell_id": 5080322}]}
{"cell_id": 34710195, "neighbors": [{"distance": 0.0, "cell_id": 33817899}, {"distance": 0.0, "cell_id": 34710195}, {"distance": 0.48189921126737956, "cell_id": 43390954}, {"distance": 0.4879399007350691, "cell_id": 53667562}, {"distance": 0.5054604666655181, "cell_id": 53866735}, {"distance": 0.5130735609642998, "cell_id": 12593081}, {"distance": 0.5130735609642998, "cell_id": 12992704}, {"distance": 0.5159252318976161, "cell_id": 53517154}, {"distance": 0.547669129319318, "cell_id": 45038134}, {"distance": 0.5523441880973903, "cell_id": 44142337}, {"distance": 0.5543826117784394, "cell_id": 16303281}, {"distance": 0.5675362573770159, "cell_id": 52750094}, {"distance": 0.5675362573770159, "cell_id": 61775686}, {"distance": 0.5717852320916406, "cell_id": 51084638}, {"distance": 0.5799731243415738, "cell_id": 15196196}, {"distance": 0.590980767034079, "cell_id": 50913522}, {"distance": 0.6052205527510841, "cell_id": 16312937}, {"distance": 0.6064640888181764, "cell_id": 39736631}, {"distance": 0.6102648803582389, "cell_id": 33669519}, {"distance": 0.6102648803582389, "cell_id": 34561815}, {"distance": 0.6154782518844387, "cell_id": 51152460}, {"distance": 0.616437792007889, "cell_id": 48370110}, {"distance": 0.6200490807930268, "cell_id": 44477773}, {"distance": 0.620121634005225, "cell_id": 51076915}, {"distance": 0.627329696584522, "cell_id": 47949414}]}
{"cell_id": 7802948, "neighbors": [{"distance": 0.0, "cell_id": 7541724}, {"distance": 0.0, "cell_id": 7802948}, {"distance": 0.5451089212603693, "cell_id": 7534742}, {"distance": 0.5451089212603693, "cell_id": 7792225}, {"distance": 0.7633855554912773, "cell_id": 7459811}, {"distance": 0.7633855554912773, "cell_id": 7536458}, {"distance": 0.7633855554912773, "cell_id": 7794825}, {"distance": 0.7809397180065489, "cell_id": 7534435}, {"distance": 0.7809397180065489, "cell_id": 7791768}, {"distance": 0.8076439105514693, "cell_id": 7454413}, {"distance": 0.8076439105514693, "cell_id": 7542063}, {"distance": 0.8076439105514693, "cell_id": 7803506}, {"distance": 0.8202525515287687, "cell_id": 7542277}, {"distance": 0.8202525515287687, "cell_id": 7803820}, {"distance": 0.8211321345285852, "cell_id": 7451945}, {"distance": 0.8211321345285852, "cell_id": 7542700}, {"distance": 0.8211321345285852, "cell_id": 7804501}, {"distance": 0.8223040976076752, "cell_id": 7542111}, {"distance": 0.8223040976076752, "cell_id": 7803573}, {"distance": 0.836553518268298, "cell_id": 7535356}, {"distance": 0.836553518268298, "cell_id": 7793143}, {"distance": 0.8473930859766599, "cell_id": 7453925}, {"distance": 0.8473930859766599, "cell_id": 7532425}, {"distance": 0.8473930859766599, "cell_id": 7788790}, {"distance": 0.8577959016868808, "cell_id": 7537541}]}
{"cell_id": 18284166, "neighbors": [{"distance": 0.0, "cell_id": 18284166}, {"distance": 0.0, "cell_id": 42144748}, {"distance": 0.3326035822458916, "cell_id": 20451256}, {"distance": 0.3326035822458916, "cell_id": 23177298}, {"distance": 0.3326035822458916, "cell_id": 28731395}, {"distance": 0.43826445324080965, "cell_id": 18574135}, {"distance": 0.43826445324080965, "cell_id": 42139586}, {"distance": 0.48130913595508384, "cell_id": 17951147}, {"distance": 0.48130913595508384, "cell_id": 42154093}, {"distance": 0.4995179890879576, "cell_id": 6433031}, {"distance": 0.4995179890879576, "cell_id": 6907347}, {"distance": 0.5114650594445956, "cell_id": 18225670}, {"distance": 0.5114650594445956, "cell_id": 42135181}, {"distance": 0.5166222840856423, "cell_id": 6323368}, {"distance": 0.5166222840856423, "cell_id": 6916137}, {"distance": 0.5217361148585429, "cell_id": 18465232}, {"distance": 0.5217361148585429, "cell_id": 42142867}, {"distance": 0.5233659146979132, "cell_id": 3284009}, {"distance": 0.5233659146979132, "cell_id": 3590573}, {"distance": 0.5312164960898458, "cell_id": 6379922}, {"distance": 0.5312164960898458, "cell_id": 6918437}, {"distance": 0.5421211363480042, "cell_id": 18762359}, {"distance": 0.5421211363480042, "cell_id": 42142529}, {"distance": 0.569803157021861, "cell_id": 17875142}, {"distance": 0.569803157021861, "cell_id": 42158783}]}
{"cell_id": 57719463, "neighbors": [{"distance": 0.0, "cell_id": 57719463}, {"distance": 0.2590955690445496, "cell_id": 57659954}, {"distance": 0.2788042900637496, "cell_id": 57302135}, {"distance": 0.2799455909554843, "cell_id": 58046836}, {"distance": 0.2867939791581781, "cell_id": 56589281}, {"distance": 0.2867939791581781, "cell_id": 59376249}, {"distance": 0.2870699100962709, "cell_id": 58006277}, {"distance": 0.2870699100962709, "cell_id": 59725499}, {"distance": 0.2878982997191029, "cell_id": 56581682}, {"distance": 0.29187919392866013, "cell_id": 58010636}, {"distance": 0.29234357713600245, "cell_id": 56430389}, {"distance": 0.3038371473810605, "cell_id": 56643368}, {"distance": 0.3038371473810605, "cell_id": 59389685}, {"distance": 0.304555805259517, "cell_id": 57367848}, {"distance": 0.3086665163092873, "cell_id": 57395192}, {"distance": 0.3086665163092873, "cell_id": 59575001}, {"distance": 0.30989699734370274, "cell_id": 57292507}, {"distance": 0.30989699734370274, "cell_id": 59549779}, {"distance": 0.31158770564306154, "cell_id": 57142943}, {"distance": 0.31422531137870807, "cell_id": 57599434}, {"distance": 0.31422531137870807, "cell_id": 59625386}, {"distance": 0.3168025464874288, "cell_id": 57178244}, {"distance": 0.3168025464874288, "cell_id": 59521693}, {"distance": 0.31858290146389673, "cell_id": 57190616}, {"distance": 0.31858290146389673, "cell_id": 59524779}]}
"""



# Section: cellxgene_census-src-cellxgene_census-_experiment

# Copyright (c) 2022, Chan Zuckerberg Initiative
#
# Licensed under the MIT License.

"""Experiments handler.

Contains methods to retrieve SOMA Experiments.
"""

import re

import tiledbsoma as soma


def _get_experiment_name(organism: str) -> str:
    """Given an organism name, return the experiment name."""
    # lower/snake case the organism name to find the experiment name
    return re.sub(r"[ ]+", "_", organism).lower()


def _get_experiment(census: soma.Collection, organism: str) -> soma.Experiment:
    """Given a census :class:`tiledbsoma.Collection`, return the experiment for the named organism.
    Organism matching is somewhat flexible, attempting to map from human-friendly
    names to the underlying collection element name.

    Args:
        census:
            The census.
        organism:
            The organism name, e.g., ``"Homo sapiens"``.

    Returns:
        An :class:`tiledbsoma.Experiment` object with the requested experiment.

    Raises:
        ValueError: if unable to find the specified organism.

    Lifecycle:
        maturing

    Examples:
        >>> human = get_experiment(census, "homo sapiens")

        >>> human = get_experiment(census, "Homo sapiens")

        >>> human = get_experiment(census, "homo_sapiens")
    """
    exp_name = _get_experiment_name(organism)

    if exp_name not in census["census_data"]:
        raise ValueError(f"Unknown organism {organism} - does not exist")
    exp = census["census_data"][exp_name]
    if exp.soma_type != "SOMAExperiment":
        raise ValueError(f"Unknown organism {organism} - not a SOMA Experiment")

    return exp



# Section: notebooks-experimental-pca

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




# Section: cellxgene_census-src-cellxgene_census-experimental-pp-_online

import numba
import numpy as np
import numpy.typing as npt


class MeanVarianceAccumulator:
    """Online mean/variance for n_variables over n_samples, where the samples are
    divided into n_batches (n_batches << n_samples). Accumulates each batch separately.

    Batches implemented using Chan's parallel adaptation of Welford's online algorithm.

    References:
        https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
    and
        Knuth, Art of Computer Programming, volume II
    """

    def __init__(
        self,
        n_batches: int,
        n_samples: npt.NDArray[np.int64],
        n_variables: int,
        ddof: int = 1,
        nnz_only: bool = False,
    ):
        if n_samples.sum() <= 0:
            raise ValueError("No samples provided - can't calculate mean or variance.")

        self.nnz_only = nnz_only
        self.ddof = ddof
        self.n_batches = n_batches
        self.n_samples = n_samples
        self.n = np.zeros((n_batches, n_variables), dtype=np.int32)
        self.u = np.zeros((n_batches, n_variables), dtype=np.float64)
        self.M2 = np.zeros((n_batches, n_variables), dtype=np.float64)

        if self.nnz_only and self.n_batches > 1:
            raise ValueError("nnz_only not implemented for n_batches > 1")

    def update(
        self,
        var_vec: npt.NDArray[np.int64],
        val_vec: npt.NDArray[np.float32],
        batch_vec: npt.NDArray[np.int64] | None = None,
    ) -> None:
        if self.n_batches == 1:
            assert batch_vec is None
            _mbomv_update_single_batch(var_vec, val_vec, self.n, self.u, self.M2)
        else:
            assert batch_vec is not None
            _mbomv_update_by_batch(batch_vec, var_vec, val_vec, self.n, self.u, self.M2)

    def finalize(
        self,
    ) -> tuple[
        npt.NDArray[np.float64],
        npt.NDArray[np.float64],
        npt.NDArray[np.float64],
        npt.NDArray[np.float64],
    ]:
        # correct each batch to account for sparsity.
        # if nnz_only, the correction is not needed as we only do mean/average over nonzero values
        if not self.nnz_only:
            _mbomv_sparse_correct_batches(self.n_batches, self.n_samples, self.n, self.u, self.M2)

        # compute u, var for each batch
        batches_u = self.u

        # Note: if N-ddof is less than or equal to 0, we will return Inf - this is consistent
        # with the numpy.var behavior.
        with np.errstate(divide="ignore", invalid="ignore"):
            batches_var = (self.M2.T / np.maximum(0, (self.n_samples - self.ddof))).T

        # accum all batches using Chan's
        all_u, all_M2 = _mbomv_combine_batches(self.n_batches, self.n_samples, self.u, self.M2)

        with np.errstate(divide="ignore"):
            if self.nnz_only:
                all_var = all_M2 / np.maximum(0, self.n - self.ddof)
                all_var = all_var[0]
            else:
                all_var = all_M2 / np.maximum(0, (self.n_samples.sum() - self.ddof))

        return batches_u, batches_var, all_u, all_var


class MeanAccumulator:
    def __init__(self, n_samples: int, n_variables: int, nnz_only: bool = False):
        if n_samples <= 0:
            raise ValueError("No samples provided - can't calculate mean.")

        if n_variables <= 0:
            raise ValueError("No variables provided - can't calculate mean.")

        self.u = np.zeros(n_variables, dtype=np.float64)
        self.n_samples = n_samples
        self.nnz_only = nnz_only
        # If we want to exclude zeros, we need to keep track of the denominator
        self.n = np.zeros(n_variables)

    def update(self, var_vec: npt.NDArray[np.int64], val_vec: npt.NDArray[np.float32]) -> None:
        if self.nnz_only:
            _update_mean_and_n_vectors(self.u, self.n, var_vec, val_vec)
        else:
            _update_mean_vector(self.u, var_vec, val_vec)

    def finalize(self) -> npt.NDArray[np.float64]:
        if self.nnz_only:
            return self.u / self.n
        else:
            return self.u / self.n_samples


class CountsAccumulator:
    def __init__(self, n_batches: int, n_variables: int, clip_val: npt.NDArray[np.float64]):
        self.n_batches = n_batches
        self.n_variables = n_variables
        self.clip_val = clip_val
        self.counts_sum = np.zeros((n_batches, n_variables), dtype=np.float64)  # clipped
        self.squared_counts_sum = np.zeros((n_batches, n_variables), dtype=np.float64)  # clipped

    def update(
        self,
        var_vec: npt.NDArray[np.int64],
        val_vec: npt.NDArray[np.float32],
        batch_vec: npt.NDArray[np.int64] | None = None,
    ) -> None:
        if self.n_batches == 1:
            assert batch_vec is None
            _accum_clipped_counts(
                self.counts_sum[0],
                self.squared_counts_sum[0],
                var_vec,
                val_vec,
                self.clip_val[0],
            )
        else:
            assert batch_vec is not None
            _accum_clipped_counts_by_batch(
                self.counts_sum,
                self.squared_counts_sum,
                batch_vec,
                var_vec,
                val_vec,
                self.clip_val,
            )

    def finalize(self) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
        return self.counts_sum, self.squared_counts_sum


"""
Private performance related top-level functions.

Down the road, it would be nice to hide via @jitclass, but it is still experimental and
there appear to be performance issues.
"""


@numba.jit(
    [
        numba.void(
            numba.int64[:],
            numba.types.Array(numba.int64, 1, "C", readonly=True),
            numba.types.Array(numba.float32, 1, "C", readonly=True),
            numba.int32[:, :],
            numba.float64[:, :],
            numba.float64[:, :],
        ),
        numba.void(
            numba.int64[:],
            numba.types.Array(numba.int32, 1, "C", readonly=True),
            numba.types.Array(numba.float32, 1, "C", readonly=True),
            numba.int32[:, :],
            numba.float64[:, :],
            numba.float64[:, :],
        ),
    ],
    nopython=True,
    nogil=True,
)  # type: ignore[misc]  # See https://github.com/numba/numba/issues/7424
def _mbomv_update_by_batch(
    batch_vec: npt.NDArray[np.int64],
    var_vec: npt.NDArray[np.int64],
    val_vec: npt.NDArray[np.float32],
    n: npt.NDArray[np.int32],
    u: npt.NDArray[np.float64],
    M2: npt.NDArray[np.float64],
) -> None:
    """Incrementally accumulate mean and sum of square of distance from mean using
    Welford's online method.
    """
    for batch, col, val in zip(batch_vec, var_vec, val_vec):
        u_prev = u[batch, col]
        M2_prev = M2[batch, col]
        n[batch, col] += 1
        u[batch, col] = u_prev + (val - u_prev) / n[batch, col]
        M2[batch, col] = M2_prev + (val - u_prev) * (val - u[batch, col])


@numba.jit(
    [
        numba.void(
            numba.types.Array(numba.int32, 1, "C", readonly=True),
            numba.types.Array(numba.float32, 1, "C", readonly=True),
            numba.int32[:, :],
            numba.float64[:, :],
            numba.float64[:, :],
        ),
        numba.void(
            numba.types.Array(numba.int64, 1, "C", readonly=True),
            numba.types.Array(numba.float32, 1, "C", readonly=True),
            numba.int32[:, :],
            numba.float64[:, :],
            numba.float64[:, :],
        ),
    ],
    nopython=True,
    nogil=True,
)  # type: ignore[misc]  # See https://github.com/numba/numba/issues/7424
def _mbomv_update_single_batch(
    var_vec: npt.NDArray[np.int64],
    val_vec: npt.NDArray[np.float32],
    n: npt.NDArray[np.int32],
    u: npt.NDArray[np.float64],
    M2: npt.NDArray[np.float64],
) -> None:
    """Incrementally accumulate mean and sum of square of distance from mean using
    Welford's online method.
    """
    for col, val in zip(var_vec, val_vec):
        u_prev = u[0, col]
        M2_prev = M2[0, col]
        n[0, col] += 1
        u[0, col] = u_prev + (val - u_prev) / n[0, col]
        M2[0, col] = M2_prev + (val - u_prev) * (val - u[0, col])


@numba.jit(
    numba.void(
        numba.int64,
        numba.int64[:],
        numba.int32[:, :],
        numba.float64[:, :],
        numba.float64[:, :],
    ),
    nopython=True,
    nogil=True,
)  # type: ignore[misc]  # See https://github.com/numba/numba/issues/7424
def _mbomv_sparse_correct_batches(
    n_batches: int,
    n_samples: npt.NDArray[np.int64],
    n: npt.NDArray[np.int32],
    u: npt.NDArray[np.float64],
    M2: npt.NDArray[np.float64],
) -> None:
    """Finalize incremental accumulators to account for missing elements (due to sparse
    input). Non-sparse and sparse combined using Chan's parallel adaptation of Welford's.
    The code assumes the sparse elements are all zero.
    """
    for batch in range(n_batches):
        n_b = n_samples[batch] - n[batch]
        delta = -u[batch]  # assumes u_b == 0
        _u = (n[batch] * u[batch]) / n_samples[batch]
        _M2 = M2[batch] + delta**2 * n[batch] * n_b / n_samples[batch]  # assumes M2_b == 0
        u[batch] = _u
        M2[batch] = _M2
        n[batch] = n_samples[batch]


@numba.jit(
    numba.types.Tuple((numba.float64[:], numba.float64[:]))(
        numba.int64, numba.int64[:], numba.float64[:, :], numba.float64[:, :]
    ),
    nopython=True,
    nogil=True,
)  # type: ignore[misc]  # See https://github.com/numba/numba/issues/7424
def _mbomv_combine_batches(
    n_batches: int,
    n_samples: npt.NDArray[np.int64],
    u: npt.NDArray[np.float64],
    M2: npt.NDArray[np.float64],
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """Combine all batches using Chan's parallel adaptation of Welford's.

    Returns tuple of (u, M2).
    """
    # initialize with first batch that contains samples
    for first_batch in range(0, n_batches):
        if n_samples[first_batch] > 0:
            acc_n = n_samples[first_batch]
            acc_u = u[first_batch].copy()
            acc_M2 = M2[first_batch].copy()
            break

    # TODO: does not handle case where there is no data, i.e., n_samples.sum() == 0

    for batch in range(first_batch + 1, n_batches):
        # ignore batches with no data
        if n_samples[batch] == 0:
            continue

        n = acc_n + n_samples[batch]
        delta = u[batch] - acc_u
        _u = (acc_n * acc_u + n_samples[batch] * u[batch]) / n
        _M2 = acc_M2 + M2[batch] + delta**2 * acc_n * n_samples[batch] / n
        # TODO: reduce memory allocs?
        acc_n = n
        acc_u = _u
        acc_M2 = _M2

    return acc_u, acc_M2


@numba.jit(
    [
        numba.void(
            numba.float64[:],
            numba.float64[:],
            numba.types.Array(numba.int32, 1, "C", readonly=True),
            numba.types.Array(numba.float32, 1, "C", readonly=True),
            numba.float64[:],
        ),
        numba.void(
            numba.float64[:],
            numba.float64[:],
            numba.types.Array(numba.int64, 1, "C", readonly=True),
            numba.types.Array(numba.float32, 1, "C", readonly=True),
            numba.float64[:],
        ),
    ],
    nopython=True,
    nogil=True,
)  # type: ignore[misc]  # See https://github.com/numba/numba/issues/7424
def _accum_clipped_counts(
    counts_sum: npt.NDArray[np.float64],
    squared_counts_sum: npt.NDArray[np.float64],
    var_dim: npt.NDArray[np.int64],
    data: npt.NDArray[np.float64],
    clip_val: npt.NDArray[np.float64],
) -> None:
    for col, val in zip(var_dim, data):
        if val > clip_val[col]:
            val = clip_val[col]
        counts_sum[col] += val
        squared_counts_sum[col] += val**2


@numba.jit(
    [
        numba.void(
            numba.float64[:, :],
            numba.float64[:, :],
            numba.int64[:],
            numba.types.Array(numba.int32, 1, "C", readonly=True),
            numba.types.Array(numba.float32, 1, "C", readonly=True),
            numba.float64[:, :],
        ),
        numba.void(
            numba.float64[:, :],
            numba.float64[:, :],
            numba.int64[:],
            numba.types.Array(numba.int64, 1, "C", readonly=True),
            numba.types.Array(numba.float32, 1, "C", readonly=True),
            numba.float64[:, :],
        ),
    ],
    nopython=True,
    nogil=True,
)  # type: ignore[misc]  # See https://github.com/numba/numba/issues/7424
def _accum_clipped_counts_by_batch(
    counts_sum: npt.NDArray[np.float64],
    squared_counts_sum: npt.NDArray[np.float64],
    batch: npt.NDArray[np.int64],
    var_dim: npt.NDArray[np.int64],
    data: npt.NDArray[np.float32],
    clip_val: npt.NDArray[np.float64],
) -> None:
    for bid, col, val in zip(batch, var_dim, data):
        if val > clip_val[bid, col]:
            val = clip_val[bid, col]
        counts_sum[bid, col] += val
        squared_counts_sum[bid, col] += val**2


@numba.jit(
    [
        numba.void(
            numba.float64[:],
            numba.types.Array(numba.int32, 1, "C", readonly=True),
            numba.types.Array(numba.float32, 1, "C", readonly=True),
        ),
        numba.void(
            numba.float64[:],
            numba.types.Array(numba.int64, 1, "C", readonly=True),
            numba.types.Array(numba.float32, 1, "C", readonly=True),
        ),
    ],
    nopython=True,
    nogil=True,
)  # type: ignore[misc]  # See https://github.com/numba/numba/issues/7424
def _update_mean_vector(
    u: npt.NDArray[np.float64],
    var_vec: npt.NDArray[np.int64],
    val_vec: npt.NDArray[np.float32],
) -> None:
    for col, val in zip(var_vec, val_vec):
        u[col] = u[col] + val


@numba.jit(
    [
        numba.void(
            numba.float64[:],
            numba.float64[:],
            numba.types.Array(numba.int32, 1, "C", readonly=True),
            numba.types.Array(numba.float32, 1, "C", readonly=True),
        ),
        numba.void(
            numba.float64[:],
            numba.float64[:],
            numba.types.Array(numba.int64, 1, "C", readonly=True),
            numba.types.Array(numba.float32, 1, "C", readonly=True),
        ),
    ],
    nopython=True,
    nogil=True,
)  # type: ignore[misc]  # See https://github.com/numba/numba/issues/7424
def _update_mean_and_n_vectors(
    u: npt.NDArray[np.float64],
    n: npt.NDArray[np.float64],
    var_vec: npt.NDArray[np.int64],
    val_vec: npt.NDArray[np.float32],
) -> None:
    for col, val in zip(var_vec, val_vec):
        u[col] = u[col] + val
        n[col] = n[col] + 1



# Section: cellxgene_census-tests-test_lts_compat

"""
Compatibility tests between the installed verison of cellxgene-census and
a named LTS release. Primarilly intended to be driven by a periodic GHA.

Where there are known and accepted incompatibilities, use `pytest.skip`
to codify them.
"""

from __future__ import annotations

from collections import deque
from collections.abc import Iterator, Sequence
from typing import Literal, TypeAlias, get_args

import pyarrow as pa
import pytest
import tiledbsoma as soma

import cellxgene_census

SOMATypeNames = Literal[
    "SOMACollection",
    "SOMAExperiment",
    "SOMAMeasurement",
    "SOMADataFrame",
    "SOMASparseNDArray",
    "SOMADenseNDArray",
]
CollectionTypeNames = ["SOMACollection", "SOMAExperiment", "SOMAMeasurement"]

SOMATypes: TypeAlias = (
    soma.Collection | soma.DataFrame | soma.SparseNDArray | soma.DenseNDArray | soma.Experiment | soma.Measurement
)


def walk_census(
    census: soma.Collection, filter_types: Sequence[SOMATypeNames] | None = None
) -> Iterator[tuple[str, SOMATypes]]:
    assert census.soma_type == "SOMACollection"
    filter_types = filter_types or get_args(SOMATypeNames)
    items_to_check = deque([("census", census)])
    while items_to_check:
        key, val = items_to_check.popleft()
        if val.soma_type in CollectionTypeNames:
            items_to_check.extend(val.items())

        if val.soma_type not in filter_types:
            continue

        yield key, val


@pytest.mark.lts_compat_check
def test_open(census_version: str) -> None:
    """
    Verify we can open and walk the collections, get metadata and read schema on non-collections
    """

    with cellxgene_census.open_soma(census_version=census_version) as census:
        for name, item in walk_census(census):
            assert name
            assert list(item.metadata)
            if item.soma_type not in CollectionTypeNames:
                assert isinstance(item.schema, pa.Schema)


@pytest.mark.lts_compat_check
def test_read_dataframe(census_version: str) -> None:
    """
    Verify we can read at least one row of dataframes
    """
    with cellxgene_census.open_soma(census_version=census_version) as census:
        for name, sdf in walk_census(census, filter_types=["SOMADataFrame"]):
            assert name
            # the Census should have no zero-length DataFrames
            assert len(sdf)
            df = sdf.read(coords=([0],)).concat().to_pandas()
            assert len(df) == 1
            assert df.shape == (1, len(sdf.keys()))


@pytest.mark.lts_compat_check
def test_read_arrays(census_version: str) -> None:
    """
    Verify we can read from NDArray
    """
    with cellxgene_census.open_soma(census_version=census_version) as census:
        for name, sarr in walk_census(census, filter_types=["SOMASparseNDArray", "SOMADenseNDArray"]):
            assert name
            assert isinstance(sarr.shape, tuple)

            # There are currently no Census schema versions using DenseNDArray

            assert sarr.soma_type == "SOMASparseNDArray"
            assert sarr.nnz
            tbl = sarr.read(coords=(slice(100),)).tables().concat()
            assert len(tbl) > 0



# Section: cellxgene_census-src-cellxgene_census-experimental-__init__

"""Experimental API for the CELLxGENE Discover Census."""

from ._embedding import (
    get_all_available_embeddings,
    get_all_census_versions_with_embedding,
    get_embedding,
    get_embedding_metadata,
    get_embedding_metadata_by_name,
)
from ._embedding_search import NeighborObs, find_nearest_obs, predict_obs_metadata

__all__ = [
    "get_embedding",
    "get_embedding_metadata",
    "get_embedding_metadata_by_name",
    "get_all_available_embeddings",
    "get_all_census_versions_with_embedding",
    "find_nearest_obs",
    "NeighborObs",
    "predict_obs_metadata",
]



# Section: notebooks-experimental-highly_variable_genes

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




# Section: cellxgene_census-src-cellxgene_census-experimental-_embedding

# Copyright (c) 2022, Chan Zuckerberg Initiative
#
# Licensed under the MIT License.

"""Methods to support simplifed access to community contributed embeddings."""

from __future__ import annotations

import json
import warnings
from typing import Any, cast

import numpy as np
import numpy.typing as npt
import pandas as pd
import pyarrow as pa
import requests
import tiledbsoma as soma

from cellxgene_census._util import _user_agent

from .._open import get_default_soma_context, open_soma
from .._release_directory import (
    CensusVersionDescription,
    CensusVersionName,
    get_census_version_description,
    get_census_version_directory,
)

CELL_CENSUS_EMBEDDINGS_MANIFEST_URL = "https://contrib.cellxgene.cziscience.com/contrib/cell-census/contributions.json"


def get_embedding_metadata(embedding_uri: str, context: soma.options.SOMATileDBContext | None = None) -> dict[str, Any]:
    """Read embedding metadata and return as a Python dict.

    Args:
        embedding_uri:
            The embedding URI.
        context:
            A custom :class:`tiledbsoma.SOMATileDBContext` which will be used to open the SOMA object. Optional,
            defaults to ``None``.

    Returns:
        A Python dictionary containing metadata describing the embedding.

    Examples:
        >>> get_experiment_metadata(uri)

    """
    # Allow the user to override context for exceptional cases (e.g. the aws region)
    context = context or get_default_soma_context()

    with soma.open(embedding_uri, context=context) as E:
        # read embedding metadata and decode the JSON-encoded string
        embedding_metadata = json.loads(E.metadata["CxG_embedding_info"])
        assert isinstance(embedding_metadata, dict)

    return cast(dict[str, Any], embedding_metadata)


def _get_embedding(
    census: soma.Collection,
    census_directory: dict[CensusVersionName, CensusVersionDescription],
    census_version: str,
    embedding_uri: str,
    obs_soma_joinids: npt.NDArray[np.int64] | pa.Array,
    context: soma.options.SOMATileDBContext | None = None,
) -> npt.NDArray[np.float32]:
    """Private. Like get_embedding, but accepts a Census object and a Census directory."""
    if isinstance(obs_soma_joinids, pa.Array | pa.ChunkedArray | pd.Series):
        obs_soma_joinids = obs_soma_joinids.to_numpy()
    assert isinstance(obs_soma_joinids, np.ndarray)
    if obs_soma_joinids.dtype != np.int64:
        raise TypeError("obs_soma_joinids must be array of int64")

    # Allow the user to override context for exceptional cases (e.g. the aws region)
    context = context or get_default_soma_context()

    # Attempt to resolve census version aliases
    resolved_census_version = census_directory.get(census_version, None)

    with soma.open(embedding_uri, context=context) as E:
        embedding_metadata = json.loads(E.metadata["CxG_embedding_info"])

        if resolved_census_version is None:
            warnings.warn(
                "Unable to determine Census version - skipping validation of Census and embedding version.",
                stacklevel=1,
            )
        elif resolved_census_version != census_directory.get(embedding_metadata["census_version"], None):
            raise ValueError("Census and embedding mismatch - census_version not equal")

        with open_soma(census_version=census_version, context=context) as census:
            experiment_name = embedding_metadata["experiment_name"]
            if experiment_name not in census["census_data"]:
                raise ValueError("Census and embedding mismatch - experiment_name does not exist")
            measurement_name = embedding_metadata["measurement_name"]
            if measurement_name not in census["census_data"][experiment_name].ms:
                raise ValueError("Census and embedding mismatch - measurement_name does not exist")

        embedding_shape = (len(obs_soma_joinids), E.shape[1])
        embedding = np.full(embedding_shape, np.nan, dtype=np.float32, order="C")

        obs_indexer = soma.IntIndexer(obs_soma_joinids, context=E.context)
        for tbl in E.read(coords=(obs_soma_joinids,)).tables():
            obs_idx = obs_indexer.get_indexer(tbl.column("soma_dim_0").to_numpy())
            feat_idx = tbl.column("soma_dim_1").to_numpy()
            emb = tbl.column("soma_data")

            indices = obs_idx * E.shape[1] + feat_idx
            np.put(embedding.reshape(-1), indices, emb)

    return embedding


def get_embedding(
    census_version: str,
    embedding_uri: str,
    obs_soma_joinids: npt.NDArray[np.int64] | pa.Array,
    context: soma.options.SOMATileDBContext | None = None,
) -> npt.NDArray[np.float32]:
    """Read cell (obs) embeddings and return as a dense :class:`numpy.ndarray`. Any cells without
    an embedding will return NaN values.

    Args:
        census_version:
            The Census version tag, e.g., ``"2023-12-15"``. Used to verify that the contents of
            the embedding contain embedded cells from the same Census version.
        embedding_uri:
            The URI containing the embedding data.
        obs_soma_joinids:
            The slice of the embedding to fetch and return.
        context:
            A custom :class:`tiledbsoma.SOMATileDBContext` which will be used to open the SOMA object.
            Optional, defaults to ``None``.

    Returns:
        A :class:`numpy.ndarray` containing the embeddings. Embeddings are positionally
        indexed by the ``obs_soma_joinids``. In other words, the cell identified by
        ``obs_soma_joinids[i]`` corresponds to the ``ith`` position in the returned
        :class:`numpy.ndarray`.

    Raises:
        ValueError: if the Census and embedding are mismatched.

    Lifecycle:
        experimental

    Examples:
        >>> obs_somaids_to_fetch = np.array([10,11], dtype=np.int64)
        >>> emb = cellxgene_census.experimental.get_embedding('2023-12-15', embedding_uri, obs_somaids_to_fetch)
        >>> emb.shape
        (2, 200)
        >>> emb[:, 0:4]
        array([[ 0.02954102,  1.0390625 , -0.14550781, -0.40820312],
            [-0.00224304,  1.265625  ,  0.05883789, -0.7890625 ]],
            dtype=float32)

    """
    census_directory = get_census_version_directory()

    with open_soma(census_version=census_version, context=context) as census:
        return _get_embedding(
            census, census_directory, census_version, embedding_uri, obs_soma_joinids, context=context
        )


def get_embedding_metadata_by_name(
    embedding_name: str, organism: str, census_version: str, embedding_type: str | None = "obs_embedding"
) -> dict[str, Any]:
    """Return metadata for a specific embedding. If more embeddings match the query parameters,
    the most recent one will be returned.

    Args:
        embedding_name:
            The name of the embedding, e.g. "scvi".
        organism:
            The organism for which the embedding is associated.
        census_version:
            The Census version tag, e.g., ``"2023-12-15"``.
        embedding_type:
            Either "obs_embedding" or "var_embedding". Defaults to "obs_embedding".

    Returns:
        A dictionary containing metadata describing the embedding.

    Raises:
        ValueError: if no embeddings are found for the specified query parameters.

    """
    census_version_description = get_census_version_description(census_version)
    resolved_census_version = census_version_description["release_build"]

    response = requests.get(CELL_CENSUS_EMBEDDINGS_MANIFEST_URL, headers={"User-Agent": _user_agent()})
    response.raise_for_status()

    manifest = cast(dict[str, dict[str, Any]], response.json())
    embeddings = []
    for _, obj in manifest.items():
        if (
            obj["embedding_name"] == embedding_name
            and obj["experiment_name"] == organism
            and obj["data_type"] == embedding_type
            and obj["census_version"] == resolved_census_version
        ):
            embeddings.append(obj)

    if len(embeddings) == 0:
        raise ValueError(
            f"No embeddings found for {embedding_name}, {organism}, {resolved_census_version}, {embedding_type}"
        )

    return sorted(embeddings, key=lambda x: x["submission_date"])[-1]


def get_all_available_embeddings(census_version: str) -> list[dict[str, Any]]:
    """Return a dictionary of all available embeddings for a given Census version.

    Args:
        census_version:
            The Census version tag, e.g., ``"2023-12-15"``.

    Returns:
        A list of dictionaries, each containing metadata describing an available embedding.

    Examples:
        >>> get_all_available_embeddings('2023-12-15')
        [{
            'experiment_name': 'experiment_1',
            'measurement_name': 'RNA',
            'organism': "homo_sapiens",
            'census_version': '2023-12-15',
            'n_embeddings': 1000,
            'n_features': 200,
            'uri': 's3://bucket/embedding_1'
        }]

    """
    # Validate census_version
    census_version_description = get_census_version_description(census_version)

    response = requests.get(CELL_CENSUS_EMBEDDINGS_MANIFEST_URL, headers={"User-Agent": _user_agent()})
    response.raise_for_status()

    embeddings = []
    manifest = response.json()
    for _, obj in manifest.items():
        if obj["census_version"] == census_version_description["release_build"]:
            embeddings.append(obj)

    return embeddings


def get_all_census_versions_with_embedding(
    embedding_name: str, organism: str, embedding_type: str | None = "obs_embedding"
) -> list[str]:
    """Get a list of all census versions that contain a specific embedding.

    Args:
        embedding_name:
            The name of the embedding, e.g. "scvi".
        organism:
            The organism for which the embedding is associated.
        embedding_type:
            The type of embedding. Defaults to "obs_embedding".

    Returns:
        A list of census versions that contain the specified embedding.
    """
    response = requests.get(CELL_CENSUS_EMBEDDINGS_MANIFEST_URL, headers={"User-Agent": _user_agent()})
    response.raise_for_status()

    manifest = response.json()
    return sorted(
        {
            obj["census_version"]
            for obj in manifest.values()
            if obj["embedding_name"] == embedding_name
            and obj["experiment_name"] == organism
            and obj["data_type"] == embedding_type
        }
    )



# Section: notebooks-analysis_demo-comp_bio_geneformer_prediction

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




# Section: cellxgene_census-tests-test_user_agent

# mypy: ignore-errors
from __future__ import annotations

import json
import os
from collections.abc import Callable
from functools import partial
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
import proxy
import pytest
import requests
from urllib3.exceptions import InsecureRequestWarning

if TYPE_CHECKING:
    from _pytest.tmpdir import TempPathFactory

import cellxgene_census

# We are forcing the requests to be insecure so we can intercept them.
pytestmark = pytest.mark.filterwarnings("ignore::urllib3.exceptions.InsecureRequestWarning")


class ProxyInstance:
    def __init__(self, proxy_obj: proxy.Proxy, logpth: Path):
        self.proxy = proxy_obj
        self.logpth = logpth

    @property
    def port(self) -> int:
        return self.proxy.flags.port


@pytest.fixture(scope="session")
def ca_certificates(tmp_path_factory: TempPathFactory) -> tuple[Path, Path, Path]:
    # Adapted from https://github.com/abhinavsingh/proxy.py/blob/a7077cf8db3bb66a6667a9d968a401e8f805e092/Makefile#L68C1-L82C49
    # TODO: Figure out if we can remove this. Currently seems neccesary for intercepting tiledb s3 requests
    cert_dir = tmp_path_factory.mktemp("ca-certificates")
    KEY_FILE = cert_dir / "ca-key.pem"
    CERT_FILE = cert_dir / "ca-cert.pem"
    SIGNING_KEY_FILE = cert_dir / "ca-signing-key.pem"
    assert proxy.common.pki.gen_private_key(key_path=KEY_FILE, password="proxy.py")
    assert proxy.common.pki.remove_passphrase(key_in_path=KEY_FILE, password="proxy.py", key_out_path=KEY_FILE)
    assert proxy.common.pki.gen_public_key(
        public_key_path=CERT_FILE, private_key_path=KEY_FILE, private_key_password="proxy.py", subject="/CN=localhost"
    )
    assert proxy.common.pki.gen_private_key(key_path=SIGNING_KEY_FILE, password="proxy.py")
    assert proxy.common.pki.remove_passphrase(
        key_in_path=SIGNING_KEY_FILE, password="proxy.py", key_out_path=SIGNING_KEY_FILE
    )
    return (KEY_FILE, CERT_FILE, SIGNING_KEY_FILE)


@pytest.fixture(scope="session")
def proxy_server(
    tmp_path_factory: TempPathFactory,
    ca_certificates: tuple[Path, Path, Path],
):
    import cellxgene_census

    tmp_path = tmp_path_factory.mktemp("proxy_logs")
    # proxy.py can override passed ca-key-file and ca-cert-file with cached ones. So we create a fresh cache for each proxy server
    cert_cache_dir = tmp_path_factory.mktemp("certificates_cache")
    proxy_log_file = tmp_path / "proxy.log"
    request_log_file = tmp_path / "proxy_requests.log"
    key_file, cert_file, signing_keyfile = ca_certificates
    assert all(p.is_file() for p in (key_file, cert_file, signing_keyfile))

    # Adapted from TestCase setup from proxy.py: https://github.com/abhinavsingh/proxy.py/blob/develop/proxy/testing/test_case.py#L23
    PROXY_PY_STARTUP_FLAGS = [
        "--num-workers",
        "1",
        "--num-acceptors",
        "1",
        "--hostname",
        "127.0.0.1",
        "--port",
        "0",
        "--plugin",
        "cellxgene_census._testing.logger_proxy.RequestLoggerPlugin",
        "--ca-key-file",
        str(key_file),
        "--ca-cert-file",
        str(cert_file),
        "--ca-signing-key-file",
        str(signing_keyfile),
        "--ca-cert-dir",
        str(cert_cache_dir),
        "--log-file",
        str(proxy_log_file),
        "--request-log-file",
        str(request_log_file),
    ]
    proxy_obj = proxy.Proxy(PROXY_PY_STARTUP_FLAGS)
    with proxy_obj:
        assert proxy_obj.acceptors
        proxy.TestCase.wait_for_server(proxy_obj.flags.port)
        proxy_instance = ProxyInstance(proxy_obj, request_log_file)

        # Now that proxy is set up, set relevant environment variables/ constants to make all request making libraries use proxy
        with pytest.MonkeyPatch.context() as mp:
            # Both requests and s3fs use these environment variables:
            mp.setenv("HTTP_PROXY", f"http://localhost:{proxy_obj.flags.port}")
            mp.setenv("HTTPS_PROXY", f"http://localhost:{proxy_obj.flags.port}")

            # s3fs
            mp.setattr(
                cellxgene_census._open,
                "DEFAULT_S3FS_KWARGS",
                {
                    "anon": True,
                    "cache_regions": True,
                    "use_ssl": False,  # So we can inspect the requests on the proxy
                },
            )

            # requests
            mp.setattr(requests, "get", partial(requests.request, "get", verify=False))

            # tiledb
            tiledb_config = cellxgene_census._open.DEFAULT_TILEDB_CONFIGURATION.copy()
            tiledb_config["vfs.s3.proxy_host"] = "localhost"
            tiledb_config["vfs.s3.proxy_port"] = str(proxy_instance.port)
            tiledb_config["vfs.s3.verify_ssl"] = "false"
            mp.setattr(
                cellxgene_census._open,
                "DEFAULT_TILEDB_CONFIGURATION",
                tiledb_config,
            )

            yield proxy_instance


@pytest.fixture
def test_specific_useragent() -> str:
    """Sets custom user agent addendum for every test so they can be uniqueley identified."""
    current_test_name = os.environ["PYTEST_CURRENT_TEST"]
    with pytest.MonkeyPatch.context() as mp:
        mp.setenv("CELLXGENE_CENSUS_USERAGENT", current_test_name)
        yield current_test_name


@pytest.fixture()
def collect_proxy_requests(proxy_server: ProxyInstance):
    """Test specific fixture exposing the proxy server.

    While a proxy server is started for every test session, this fixture
    captures only the output that is written for a specific tests. These
    logged requests can be checked to make sure have the correct headers.
    """
    # If logs have already been written, count how many
    if proxy_server.logpth.is_file():
        with proxy_server.logpth.open("r") as f:
            prev_lines = len(f.readlines())
    else:
        prev_lines = 0

    def _proxy_requests():
        # For each new log written by the test, check that the correct headers were written
        with proxy_server.logpth.open("r") as f:
            records = [json.loads(line) for line in f.readlines()]
        records = records[prev_lines:]
        return records

    # Run test
    yield _proxy_requests


@pytest.fixture(scope="session")
def small_dataset_id() -> str:
    with cellxgene_census.open_soma(census_version="latest") as census:
        census_datasets = census["census_info"]["datasets"].read().concat().to_pandas()

    small_dataset = census_datasets.nsmallest(1, "dataset_total_cell_count").iloc[0]
    assert isinstance(small_dataset.dataset_id, str)
    return small_dataset.dataset_id


def check_proxy_records(records: list[dict], *, custom_user_agent: None | str = None, min_records: int = 1) -> None:
    # Check that there aren't two CONNECT requests in a row
    prev_was_connect = False
    for record in records:
        was_connect = record["method"] == "CONNECT"
        if prev_was_connect and was_connect:
            raise AssertionError(
                "Recieved multiple connect requests in a row. Some calls aren't being intercepted by the proxy."
            )

    # Check that headers were set correctly on intercepted requests
    n_records = 0
    for record in records:
        if record["method"] == "CONNECT":
            continue
        n_records += 1
        headers = record["headers"]
        user_agent = headers["user-agent"]
        assert "cellxgene-census-python" in user_agent
        assert cellxgene_census.__version__ in user_agent
        if custom_user_agent:
            assert custom_user_agent in user_agent
    assert n_records >= min_records, f"Fewer than min_records ({min_records}) were found."


def test_proxy_fixture(collect_proxy_requests: Callable[[], list[dict]]):
    """Test that our proxy testing setup is working as expected."""
    # Should just be downloading a json
    with pytest.warns(InsecureRequestWarning):
        _ = cellxgene_census.get_census_version_directory()

    records = collect_proxy_requests()

    # Expecting a CONNECT request followed by a GET request
    assert len(records) == 2
    assert records[0]["method"] == "CONNECT"
    assert records[1]["method"] == "GET"
    assert records[1]["headers"]["host"] == "census.cellxgene.cziscience.com"
    assert "cellxgene-census-python" in records[1]["headers"]["user-agent"]


def test_download_w_proxy_fixture(
    small_dataset_id: str,
    collect_proxy_requests: Callable[[], list[dict]],
    tmp_path: Path,
    test_specific_useragent: str,
):
    # Use of collect_proxy_requests forces test to use a proxy and will check headers of requests made via that proxy
    adata_path = tmp_path / "adata.h5ad"
    cellxgene_census.download_source_h5ad(small_dataset_id, adata_path.as_posix(), census_version="latest")

    records = collect_proxy_requests()
    check_proxy_records(
        records,
        custom_user_agent=test_specific_useragent,
        min_records=3,  # Should request at least a json and the download
    )


def test_query_w_proxy_fixture(collect_proxy_requests: Callable[[], list[dict]]):
    with cellxgene_census.open_soma(census_version="stable") as census:
        _ = cellxgene_census.get_obs(census, "Mus musculus", coords=slice(100, 300))

    records = collect_proxy_requests()
    check_proxy_records(
        records,
        min_records=5,  # some metadata requests, then a lot of request from tiledb
    )


def test_embedding_headers(collect_proxy_requests: Callable[[], list[dict]]):
    import cellxgene_census.experimental

    CENSUS_VERSION = "2023-12-15"

    embeddings_metadata = cellxgene_census.experimental.get_all_available_embeddings(CENSUS_VERSION)
    metadata = embeddings_metadata[0]
    embedding_uri = (
        f"s3://cellxgene-contrib-public/contrib/cell-census/soma/{metadata['census_version']}/{metadata['id']}"
    )
    _ = cellxgene_census.experimental.get_embedding(
        CENSUS_VERSION,
        embedding_uri=embedding_uri,
        obs_soma_joinids=np.arange(100),
    )

    check_proxy_records(collect_proxy_requests())


def test_dataloader_headers(collect_proxy_requests) -> None:
    import cellxgene_census
    from cellxgene_census.experimental.ml.pytorch import ExperimentDataPipe

    soma_experiment = cellxgene_census.open_soma(census_version="latest")["census_data"]["homo_sapiens"]
    dp = ExperimentDataPipe(
        soma_experiment,
        measurement_name="RNA",
        X_name="raw",
        obs_column_names=["cell_type"],
        shuffle=False,
    )
    _ = next(iter(dp))

    records = collect_proxy_requests()
    check_proxy_records(records, min_records=5)



# Section: notebooks-experimental-pytorch

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







# Section: cellxgene_census-src-cellxgene_census-experimental-pp-_highly_variable_genes

from __future__ import annotations

import os
from collections.abc import Callable, Sequence
from concurrent import futures
from typing import Any, Literal, cast

import numpy as np
import pandas as pd
import tiledbsoma as soma
from somacore.options import SparseDFCoord

from ..._experiment import _get_experiment
from ..util._eager_iter import _EagerIterator
from ._online import CountsAccumulator, MeanVarianceAccumulator

"""
Acknowledgements: ScanPy highly variable genes implementation (scanpy.pp.highly_variable_genes), in turn
based upon the original implementation in Seurat V3.

Ref:
* https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.highly_variable_genes.html#scanpy.pp.highly_variable_genes
* github.com/scverse/scanpy

Notes:
* Occasionally, skmis.loess will fail with a ValueError similar to
    `ValueError: b'reciprocal condition number  2.0467e-15\n'`.
  or
    `ValueError: b'There are other near singularities as well. 0.090619\n'`
  This is likely caused by an excess of all-zero valued counts in a given batch or
  other low-entropy data. Ref:
    * https://github.com/scverse/scanpy/issues/1504
    * https://github.com/has2k1/scikit-misc/issues/9
    * https://discourse.scverse.org/t/error-in-highly-variable-gene-selection/276/9

  It seems possible to work around by retrying failures with addition of noise/jitter.
"""


def _get_batch_index(
    query: soma.ExperimentAxisQuery,
    batch_key: str | Sequence[str],
    batch_key_func: Callable[..., Any] | None = None,
) -> pd.Series[Any]:
    """Return categorical series representing the batch key, with codes that index the key."""
    if isinstance(batch_key, str):
        batch_key = [batch_key]
    batch_key = list(batch_key)
    assert isinstance(batch_key, list) and len(batch_key) > 0

    obs: pd.DataFrame = (
        query.obs(column_names=["soma_joinid"] + batch_key).concat().to_pandas().set_index("soma_joinid")[batch_key]
    )

    batch_series: pd.Series[Any]
    if batch_key_func is not None:
        # apply user lambda
        batch_series = obs.apply(batch_key_func, axis=1, result_type="reduce")
    elif len(batch_key) > 1:
        # if multiple keys, stringify and concat
        obs = obs.astype(str)
        batch_series = obs[cast(str, batch_key[0])]
        batch_series = batch_series.str.cat(obs[batch_key[1:]])
    else:
        # if a single key, just use it
        batch_series = obs[cast(str, batch_key[0])]

    batch_series = batch_series.astype("category")
    return batch_series.cat.remove_unused_categories()


def _highly_variable_genes_seurat_v3(
    query: soma.ExperimentAxisQuery,
    batch_key: str | Sequence[str] | None = None,
    n_top_genes: int = 1_000,
    layer: str = "raw",
    span: float = 0.3,
    max_loess_jitter: float = 1e-6,
    batch_key_func: Callable[..., Any] | None = None,
) -> pd.DataFrame:
    try:
        import skmisc.loess
    except ImportError as e:
        raise ImportError("Please install skmisc package via `pip install --user scikit-misc") from e

    batch_indexer = None
    if batch_key is not None:
        batch_index = _get_batch_index(query, batch_key, batch_key_func)
        n_batches = len(batch_index.cat.categories)
        n_samples = batch_index.value_counts().loc[batch_index.cat.categories.to_numpy()].to_numpy()
        if n_batches > 1:
            batch_indexer = soma.IntIndexer(batch_index.index.to_numpy(), context=query.experiment.context).get_indexer
            batch_codes = batch_index.cat.codes.to_numpy().astype(np.int64)
    else:
        n_batches = 1
        n_samples = np.array([query.n_obs], dtype=np.int64)

    assert n_batches == len(n_samples)
    assert query.n_obs == n_samples.sum()
    assert all(n_samples > 0)
    assert (n_batches > 1) == bool(batch_indexer)

    max_workers = (os.cpu_count() or 4) + 2
    var_indexer = query.indexer

    with futures.ThreadPoolExecutor(max_workers=max_workers) as pool:
        mvn = MeanVarianceAccumulator(n_batches, n_samples, query.n_vars)
        for arrow_tbl in _EagerIterator(query.X(layer).tables(), pool=pool):
            data = arrow_tbl["soma_data"].to_numpy()
            if batch_indexer:
                _batch_take_at_future = pool.submit(batch_indexer, arrow_tbl["soma_dim_0"])
                var_dim = var_indexer.by_var(arrow_tbl["soma_dim_1"])
                _batch_vec = batch_codes[_batch_take_at_future.result()]
                mvn.update(var_dim, data, _batch_vec)
            else:
                var_dim = var_indexer.by_var(arrow_tbl["soma_dim_1"])
                mvn.update(var_dim, data)

        batches_u, batches_var, all_u, all_var = mvn.finalize()
        del mvn

    var_df = pd.DataFrame(
        index=pd.Index(data=query.var_joinids(), name="soma_joinid"),
        data={
            "means": all_u,
            "variances": all_var,
        },
    )

    # Calculate per-batch clip_val and reg_std
    estimated_variances = np.empty((query.n_vars,), dtype=np.float64)
    reg_std = np.zeros((n_batches, query.n_vars), dtype=np.float64)
    clip_val = np.zeros((n_batches, query.n_vars), dtype=np.float64)
    for batch in range(n_batches):
        estimated_variances.fill(0)
        u = batches_u[batch]
        v = batches_var[batch]
        N = n_samples[batch]

        not_const = v > 0
        if N == 1 or not not_const.any():
            reg_std[batch].fill(1)
            clip_val[batch, :] = u
            continue

        y = np.log10(v[not_const])
        x = np.log10(u[not_const])

        jitter_magnitude: float = 0
        while True:
            try:
                # Attempt to resolve low entropy loess errors by adding jitter and retrying
                # See: https://github.com/has2k1/scikit-misc/issues/9
                if jitter_magnitude != 0:
                    _x = x + np.random.default_rng().uniform(-jitter_magnitude, jitter_magnitude, x.shape[0])
                else:
                    _x = x

                model = skmisc.loess.loess(_x, y, span=span, degree=2)
                model.fit()
                estimated_variances[not_const] = model.outputs.fitted_values
                break

            except ValueError:
                jitter_magnitude = 1e-18 if jitter_magnitude == 0 else jitter_magnitude * 10.0
                if jitter_magnitude < max_loess_jitter:
                    continue
                raise

        reg_std[batch] = np.sqrt(10**estimated_variances)
        vmax = np.sqrt(N)
        clip_val[batch] = reg_std[batch] * vmax + u

    del estimated_variances

    # Read counts again, clip and save sum of counts and sum of counts squared
    with futures.ThreadPoolExecutor(max_workers=max_workers) as pool:
        acc = CountsAccumulator(n_batches, query.n_vars, clip_val)
        for arrow_tbl in _EagerIterator(query.X(layer).tables(), pool=pool):
            data = arrow_tbl["soma_data"].to_numpy()
            if batch_indexer:
                _batch_take_at_future = pool.submit(batch_indexer, arrow_tbl["soma_dim_0"])
                var_dim = var_indexer.by_var(arrow_tbl["soma_dim_1"])
                _batch_vec = batch_codes[_batch_take_at_future.result()]
                acc.update(var_dim, data, _batch_vec)
            else:
                var_dim = var_indexer.by_var(arrow_tbl["soma_dim_1"])
                acc.update(var_dim, data)

        counts_sum, squared_counts_sum = acc.finalize()
        # Don't raise python errors for 0/0, etc, just generate Inf/NaN
        with np.errstate(divide="ignore", invalid="ignore"):
            norm_gene_vars = (1 / ((n_samples - 1) * np.square(reg_std.T))).T * (
                (n_samples * np.square(batches_u.T)).T + squared_counts_sum - 2 * counts_sum * batches_u
            )
            norm_gene_vars[np.isnan(norm_gene_vars)] = 0
        del acc, counts_sum, squared_counts_sum

    ranked_norm_gene_vars = np.argsort(np.argsort(-norm_gene_vars, axis=1), axis=1)
    ranked_norm_gene_vars = ranked_norm_gene_vars.astype(np.float32)
    num_batches_high_var = np.sum((ranked_norm_gene_vars < n_top_genes).astype(int), axis=0)
    ranked_norm_gene_vars[ranked_norm_gene_vars >= n_top_genes] = np.nan
    ma_ranked = np.ma.masked_invalid(ranked_norm_gene_vars)  # type: ignore[no-untyped-call]
    median_ranked = np.ma.median(ma_ranked, axis=0).filled(np.nan)  # type: ignore[no-untyped-call]

    var_df = var_df.assign(
        highly_variable_nbatches=pd.Series(num_batches_high_var, index=var_df.index),
        highly_variable_rank=pd.Series(median_ranked, index=var_df.index),
        variances_norm=pd.Series(np.mean(norm_gene_vars, axis=0), index=var_df.index),
    )

    sorted_index = (
        var_df[["highly_variable_rank", "highly_variable_nbatches"]]
        .sort_values(
            ["highly_variable_rank", "highly_variable_nbatches"],
            ascending=[True, False],
            na_position="last",
        )
        .index
    )
    var_df["highly_variable"] = False
    var_df.loc[sorted_index[: int(n_top_genes)], "highly_variable"] = True
    if batch_key is None:
        var_df = var_df.drop(columns=["highly_variable_nbatches"])
    return var_df


def highly_variable_genes(
    query: soma.ExperimentAxisQuery,
    n_top_genes: int = 1_000,
    layer: str = "raw",
    flavor: Literal["seurat_v3"] = "seurat_v3",
    span: float = 0.3,
    batch_key: str | Sequence[str] | None = None,
    max_loess_jitter: float = 1e-6,
    batch_key_func: Callable[..., Any] | None = None,
) -> pd.DataFrame:
    """Identify and annotate highly variable genes contained in the query results.
    The API is modelled on ScanPy `scanpy.pp.highly_variable_genes` API.
    Results returned will mimic ScanPy results. The only `flavor` available
    is the Seurat V3 method, which assumes count data in the X layer.

    See
    https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.highly_variable_genes.html#scanpy.pp.highly_variable_genes
    for more information on this method.

    Args:
        query:
            A :class:`tiledbsoma.ExperimentAxisQuery`, specifying the ``obs``/``var`` selection over which genes are
            annotated.
        n_top_genes:
            Number of genes to rank.
        layer:
            X layer used, e.g., ``"raw"``.
        flavor:
            Method used to annotate genes. Must be ``"seurat_v3"``.
        span:
            If ``flavor="seurat_v3"``, the fraction of obs/cells used to estimate the LOESS variance model fit.
        batch_key:
            If specified, gene selection will be done by batch and combined. Specify the obs column name, or list of
            column names, identifying the batches. If not specified, all gene selection is done as a single batch.
            If multiple batch keys are specified, and no batch_key_func is specified, the batch key will be generated by
            converting values to string and concatenating them.
        max_lowess_jitter:
            The maximum jitter to add to data in case of LOESS failure (can occur when dataset has low entry counts.)
        batch_key_func:
            Optional function to create a user-defined batch key. Function will be called once per row in the obs
            dataframe. Function will receive a single argument: a :class:`pandas.Series` containing values specified in
            the``batch_key`` argument.

    Returns:
        A :class:`pandas.DataFrame` containing annotations for all ``var`` values specified by the ``query`` argument.
        Annotations are identical to those produced by :func:`scanpy.pp.highly_variable_genes`.

    Raises:
        ValueError: if the flavor parameter is not ``"seurat_v3"``.


    Examples:
        Fetch :class:`pandas.DataFrame` containing var annotations for the query selection, using ``"dataset_id"`` as
        ``batch_key``.

        >>> hvg = highly_variable_genes(query, batch_key="dataset_id")

        Fetch highly variable genes, using the concatenation of ``"dataset_id"`` and ``"donor_id"`` as ``batch_key``:

        >>> hvg = highly_variable_genes(query, batch_key=["dataset_id", "donor_id"])

        Fetch highly variable genes, with a user-defined ``batch_key_func``:

        >>> hvg = highly_variable_genes(
                query,
                batch_key="donor_id",
                batch_key_func=lambda s: return "batch0" if s.donor_id == "99" else "batch1"
            )

    Lifecycle:
        experimental
    """
    if flavor != "seurat_v3":
        raise ValueError('`flavor` must be "seurat_v3"')

    return _highly_variable_genes_seurat_v3(
        query,
        n_top_genes=n_top_genes,
        layer=layer,
        span=span,
        batch_key=batch_key,
        batch_key_func=batch_key_func,
        max_loess_jitter=max_loess_jitter,
    )


def get_highly_variable_genes(
    census: soma.Collection,
    organism: str,
    measurement_name: str = "RNA",
    X_name: str = "raw",
    obs_value_filter: str | None = None,
    obs_coords: SparseDFCoord | None = None,
    var_value_filter: str | None = None,
    var_coords: SparseDFCoord | None = None,
    n_top_genes: int = 1_000,
    flavor: Literal["seurat_v3"] = "seurat_v3",
    span: float = 0.3,
    batch_key: str | Sequence[str] | None = None,
    max_loess_jitter: float = 1e-6,
    batch_key_func: Callable[..., Any] | None = None,
) -> pd.DataFrame:
    """Convience wrapper around :class:`tiledbsoma.Experiment` query and
    :func:`cellxgene_census.experimental.pp.highly_variable_genes` function, to build and execute a query, and annotate
    the query result genes (``var`` dataframe) based upon variability.

    Args:
        census:
            The Census object, usually returned by :func:`open_soma`.
        organism:
            The organism to query, usually one of ``"Homo sapiens"`` or ``"Mus musculus"``.
        measurement_name:
            The measurement object to query. Defaults to ``"RNA"``.
        X_name:
            The ``X`` layer to query. Defaults to ``"raw"``.
        obs_value_filter:
            Value filter for the ``obs`` metadata. Value is a filter query written in the
            SOMA ``value_filter`` syntax.
        obs_coords:
            Coordinates for the ``obs`` axis, which is indexed by the ``soma_joinid`` value.
            May be an ``int``, a list of ``int``, or a slice. The default, ``None``, selects all.
        var_value_filter:
            Value filter for the ``var`` metadata. Value is a filter query written in the
            SOMA ``value_filter`` syntax.
        var_coords:
            Coordinates for the ``var`` axis, which is indexed by the ``soma_joinid`` value.
            May be an ``int``, a list of ``int``, or a slice. The default, ``None``, selects all.
        n_top_genes:
            Number of genes to rank.
        flavor:
            Method used to annotate genes. Must be ``"seurat_v3"``.
        span:
            If ``flavor="seurat_v3"``, the fraction of obs/cells used to estimate the LOESS variance model fit.
        batch_key:
            If specified, gene selection will be done by batch and combined. Specify the obs column name, or list of
            column names, identifying the batches. If not specified, all gene selection is done as a single batch.
            If multiple batch keys are specified, and no batch_key_func is specified, the batch key will be generated by
            converting values to string and concatenating them.
        max_lowess_jitter:
            The maximum jitter to add to data in case of LOESS failure (can occur when dataset has low entry counts.)
        batch_key_func:
            Optional function to create a user-defined batch key. Function will be called once per row in the obs
            dataframe. Function will receive a single argument: a :class:`pandas.Series` containing values specified in
            the ``batch_key`` argument.

    Returns:
        :class:`pandas.DataFrame` containing annotations for all ``var`` values specified by the query.

    Raises:
        ValueError: if the flavor paramater is not ``"seurat_v3"``.

    See Also:
        :func:`cellxgene_census.experimental.pp.highly_variable_genes`

    Examples:
        Fetch a :class:`pandas.DataFrame` containing var annotations for a subset of the cells matching the
        ``obs_value_filter``:

        >>> hvg = get_highly_variable_genes(
                census,
                organism="Mus musculus",
                obs_value_filter="is_primary_data == True and tissue_general == 'lung'",
                n_top_genes = 500
            )

        Fetch an :class:`anndata.AnnData` with top 500 genes:

        >>> with cellxgene_census.open_soma(census_version="stable") as census:
                organism = "mus_musculus"
                obs_value_filter = "is_primary_data == True and tissue_general == 'lung'"
                # Get the highly variable genes
                hvg = cellxgene_census.experimental.pp.get_highly_variable_genes(
                    census,
                    organism=organism,
                    obs_value_filter=obs_value_filter,
                    n_top_genes = 500
                )
                # Fetch AnnData - all cells matching obs_value_filter, just the HVGs
                hvg_soma_ids = hvg[hvg.highly_variable].index.values
                adata = cellxgene_census.get_anndata(
                    census, organism=organism, obs_value_filter=obs_value_filter, var_coords=hvg_soma_ids
                )

    Lifecycle:
        experimental

    """
    exp = _get_experiment(census, organism)
    obs_coords = (slice(None),) if obs_coords is None else (obs_coords,)
    var_coords = (slice(None),) if var_coords is None else (var_coords,)
    with exp.axis_query(
        measurement_name,
        obs_query=soma.AxisQuery(value_filter=obs_value_filter, coords=obs_coords),
        var_query=soma.AxisQuery(value_filter=var_value_filter, coords=var_coords),
    ) as query:
        return highly_variable_genes(
            query,
            n_top_genes=n_top_genes,
            layer=X_name,
            flavor=flavor,
            span=span,
            batch_key=batch_key,
            batch_key_func=batch_key_func,
            max_loess_jitter=max_loess_jitter,
        )



# Section: notebooks-analysis_demo-comp_bio_explore_and_load_lung_data

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




# Section: cellxgene_census-tests-experimental-ml-test_pytorch

import pathlib
from collections.abc import Callable, Sequence
from unittest.mock import patch

import numpy as np
import pyarrow as pa
import pytest
import tiledbsoma as soma
from scipy import sparse
from scipy.sparse import coo_matrix, spmatrix
from somacore import AxisQuery
from tiledbsoma import Experiment, _factory
from tiledbsoma._collection import CollectionBase

# conditionally import torch, as it will not be available in all test environments
try:
    from torch import Tensor, float32
    from torch.utils.data._utils.worker import WorkerInfo

    from cellxgene_census.experimental.ml.encoders import BatchEncoder, LabelEncoder
    from cellxgene_census.experimental.ml.pytorch import (
        ExperimentDataPipe,
        experiment_dataloader,
        list_split,
    )
except ImportError:
    # this should only occur when not running `experimental`-marked tests
    pass


def pytorch_x_value_gen(obs_range: range, var_range: range) -> spmatrix:
    occupied_shape = (
        obs_range.stop - obs_range.start,
        var_range.stop - var_range.start,
    )
    checkerboard_of_ones = coo_matrix(np.indices(occupied_shape).sum(axis=0) % 2)
    checkerboard_of_ones.row += obs_range.start
    checkerboard_of_ones.col += var_range.start
    return checkerboard_of_ones


def pytorch_seq_x_value_gen(obs_range: range, var_range: range) -> spmatrix:
    """A sparse matrix where the values of each col are the obs_range values. Useful for checking the
    X values are being returned in the correct order."""
    data = np.vstack([list(obs_range)] * len(var_range)).flatten()
    rows = np.vstack([list(obs_range)] * len(var_range)).flatten()
    cols = np.column_stack([list(var_range)] * len(obs_range)).flatten()
    return coo_matrix((data, (rows, cols)))


@pytest.fixture
def X_layer_names() -> list[str]:
    return ["raw"]


@pytest.fixture
def obsp_layer_names() -> list[str] | None:
    return None


@pytest.fixture
def varp_layer_names() -> list[str] | None:
    return None


def add_dataframe(coll: CollectionBase, key: str, value_range: range) -> None:
    df = coll.add_new_dataframe(
        key,
        schema=pa.schema(
            [
                ("soma_joinid", pa.int64()),
                ("label", pa.large_string()),
                ("label2", pa.large_string()),
            ]
        ),
        index_column_names=["soma_joinid"],
    )
    df.write(
        pa.Table.from_pydict(
            {
                "soma_joinid": list(value_range),
                "label": [str(i) for i in value_range],
                "label2": ["c" for i in value_range],
            }
        )
    )


def add_sparse_array(
    coll: CollectionBase,
    key: str,
    obs_range: range,
    var_range: range,
    value_gen: Callable[[range, range], spmatrix],
) -> None:
    a = coll.add_new_sparse_ndarray(key, type=pa.float32(), shape=(obs_range.stop, var_range.stop))
    tensor = pa.SparseCOOTensor.from_scipy(value_gen(obs_range, var_range))
    a.write(tensor)


@pytest.fixture(scope="function")
def soma_experiment(
    tmp_path: pathlib.Path,
    obs_range: int | range,
    var_range: int | range,
    X_value_gen: Callable[[range, range], sparse.spmatrix],
    obsp_layer_names: Sequence[str],
    varp_layer_names: Sequence[str],
) -> soma.Experiment:
    with soma.Experiment.create((tmp_path / "exp").as_posix()) as exp:
        if isinstance(obs_range, int):
            obs_range = range(obs_range)
        if isinstance(var_range, int):
            var_range = range(var_range)

        add_dataframe(exp, "obs", obs_range)
        ms = exp.add_new_collection("ms")
        rna = ms.add_new_collection("RNA", soma.Measurement)
        add_dataframe(rna, "var", var_range)
        rna_x = rna.add_new_collection("X", soma.Collection)
        add_sparse_array(rna_x, "raw", obs_range, var_range, X_value_gen)

        if obsp_layer_names:
            obsp = rna.add_new_collection("obsp")
            for obsp_layer_name in obsp_layer_names:
                add_sparse_array(obsp, obsp_layer_name, obs_range, var_range, X_value_gen)

        if varp_layer_names:
            varp = rna.add_new_collection("varp")
            for varp_layer_name in varp_layer_names:
                add_sparse_array(varp, varp_layer_name, obs_range, var_range, X_value_gen)
    return _factory.open((tmp_path / "exp").as_posix())


@pytest.mark.experimental
# noinspection PyTestParametrized
@pytest.mark.parametrize(
    "obs_range,var_range,X_value_gen,use_eager_fetch",
    [(6, 3, pytorch_x_value_gen, use_eager_fetch) for use_eager_fetch in (True, False)],
)
def test_non_batched(soma_experiment: Experiment, use_eager_fetch: bool) -> None:
    exp_data_pipe = ExperimentDataPipe(
        soma_experiment,
        measurement_name="RNA",
        X_name="raw",
        obs_column_names=["label"],
        shuffle=False,
        use_eager_fetch=use_eager_fetch,
    )
    row_iter = iter(exp_data_pipe)

    row = next(row_iter)
    assert row[0].int().tolist() == [0, 1, 0]
    assert row[1].tolist() == [0]


@pytest.mark.experimental
# noinspection PyTestParametrized
@pytest.mark.parametrize(
    "obs_range,var_range,X_value_gen,use_eager_fetch",
    [(6, 3, pytorch_x_value_gen, use_eager_fetch) for use_eager_fetch in (True, False)],
)
@pytest.mark.parametrize("return_sparse_X", [True, False])
def test_uneven_soma_and_result_batches(
    soma_experiment: Experiment, use_eager_fetch: bool, return_sparse_X: bool
) -> None:
    """This is checking that batches are correctly created when they require fetching multiple chunks.

    This was added due to failures in _ObsAndXIterator.__next__.
    """
    exp_data_pipe = ExperimentDataPipe(
        soma_experiment,
        measurement_name="RNA",
        X_name="raw",
        obs_column_names=["label"],
        shuffle=False,
        batch_size=3,
        soma_chunk_size=2,
        return_sparse_X=return_sparse_X,
        use_eager_fetch=use_eager_fetch,
    )
    row_iter = iter(exp_data_pipe)

    row = next(row_iter)
    X_batch = row[0].to_dense() if return_sparse_X else row[0]
    assert X_batch.int()[0].tolist() == [0, 1, 0]
    assert row[1].tolist() == [[0], [1], [2]]


@pytest.mark.experimental
# noinspection PyTestParametrized,DuplicatedCode
@pytest.mark.parametrize(
    "obs_range,var_range,X_value_gen,use_eager_fetch",
    [(6, 3, pytorch_x_value_gen, use_eager_fetch) for use_eager_fetch in (True, False)],
)
def test_batching__all_batches_full_size(soma_experiment: Experiment, use_eager_fetch: bool) -> None:
    exp_data_pipe = ExperimentDataPipe(
        soma_experiment,
        measurement_name="RNA",
        X_name="raw",
        obs_column_names=["label"],
        batch_size=3,
        shuffle=False,
        use_eager_fetch=use_eager_fetch,
    )
    batch_iter = iter(exp_data_pipe)

    batch = next(batch_iter)
    assert batch[0].int().tolist() == [[0, 1, 0], [1, 0, 1], [0, 1, 0]]
    assert batch[1].tolist() == [[0], [1], [2]]

    batch = next(batch_iter)
    assert batch[0].int().tolist() == [[1, 0, 1], [0, 1, 0], [1, 0, 1]]
    assert batch[1].tolist() == [[3], [4], [5]]

    with pytest.raises(StopIteration):
        next(batch_iter)


@pytest.mark.experimental
# noinspection PyTestParametrized,DuplicatedCode
@pytest.mark.parametrize(
    "obs_range,var_range,X_value_gen,use_eager_fetch",
    [(range(100_000_000, 100_000_003), 3, pytorch_x_value_gen, use_eager_fetch) for use_eager_fetch in (True, False)],
)
def test_unique_soma_joinids(soma_experiment: Experiment, use_eager_fetch: bool) -> None:
    exp_data_pipe = ExperimentDataPipe(
        soma_experiment,
        measurement_name="RNA",
        X_name="raw",
        obs_column_names=["label"],
        batch_size=3,
        use_eager_fetch=use_eager_fetch,
    )

    soma_joinids = np.concatenate([batch[1][:, 0].numpy() for batch in exp_data_pipe])

    assert len(np.unique(soma_joinids)) == len(soma_joinids)


@pytest.mark.experimental
# noinspection PyTestParametrized
@pytest.mark.parametrize(
    "obs_range,var_range,X_value_gen,use_eager_fetch",
    [(5, 3, pytorch_x_value_gen, use_eager_fetch) for use_eager_fetch in (True, False)],
)
def test_batching__partial_final_batch_size(soma_experiment: Experiment, use_eager_fetch: bool) -> None:
    exp_data_pipe = ExperimentDataPipe(
        soma_experiment,
        measurement_name="RNA",
        X_name="raw",
        obs_column_names=["label"],
        batch_size=3,
        shuffle=False,
        use_eager_fetch=use_eager_fetch,
    )
    batch_iter = iter(exp_data_pipe)

    next(batch_iter)
    batch = next(batch_iter)
    assert batch[0].int().tolist() == [[1, 0, 1], [0, 1, 0]]

    with pytest.raises(StopIteration):
        next(batch_iter)


@pytest.mark.experimental
# noinspection PyTestParametrized,DuplicatedCode
@pytest.mark.parametrize(
    "obs_range,var_range,X_value_gen,use_eager_fetch",
    [(3, 3, pytorch_x_value_gen, use_eager_fetch) for use_eager_fetch in (True, False)],
)
def test_batching__exactly_one_batch(soma_experiment: Experiment, use_eager_fetch: bool) -> None:
    exp_data_pipe = ExperimentDataPipe(
        soma_experiment,
        measurement_name="RNA",
        X_name="raw",
        obs_column_names=["label"],
        batch_size=3,
        shuffle=False,
        use_eager_fetch=use_eager_fetch,
    )
    batch_iter = iter(exp_data_pipe)

    batch = next(batch_iter)
    assert batch[0].int().tolist() == [[0, 1, 0], [1, 0, 1], [0, 1, 0]]
    assert batch[1].tolist() == [[0], [1], [2]]

    with pytest.raises(StopIteration):
        next(batch_iter)


@pytest.mark.experimental
# noinspection PyTestParametrized
@pytest.mark.parametrize(
    "obs_range,var_range,X_value_gen,use_eager_fetch",
    [(6, 3, pytorch_x_value_gen, use_eager_fetch) for use_eager_fetch in (True, False)],
)
def test_batching__empty_query_result(soma_experiment: Experiment, use_eager_fetch: bool) -> None:
    exp_data_pipe = ExperimentDataPipe(
        soma_experiment,
        measurement_name="RNA",
        X_name="raw",
        obs_query=AxisQuery(coords=([],)),
        obs_column_names=["label"],
        batch_size=3,
        use_eager_fetch=use_eager_fetch,
    )
    batch_iter = iter(exp_data_pipe)

    with pytest.raises(StopIteration):
        next(batch_iter)


@pytest.mark.experimental
# noinspection PyTestParametrized
@pytest.mark.parametrize(
    "obs_range,var_range,X_value_gen,use_eager_fetch",
    [(6, 3, pytorch_x_value_gen, use_eager_fetch) for use_eager_fetch in (True, False)],
)
def test_sparse_output__non_batched(soma_experiment: Experiment, use_eager_fetch: bool) -> None:
    exp_data_pipe = ExperimentDataPipe(
        soma_experiment,
        measurement_name="RNA",
        X_name="raw",
        obs_column_names=["label"],
        return_sparse_X=True,
        shuffle=False,
        use_eager_fetch=use_eager_fetch,
    )
    batch_iter = iter(exp_data_pipe)

    batch = next(batch_iter)
    assert isinstance(batch[1], Tensor)
    assert batch[0].to_dense().tolist() == [0, 1, 0]


@pytest.mark.experimental
# noinspection PyTestParametrized
@pytest.mark.parametrize(
    "obs_range,var_range,X_value_gen,use_eager_fetch",
    [(6, 3, pytorch_x_value_gen, use_eager_fetch) for use_eager_fetch in (True, False)],
)
def test_sparse_output__batched(soma_experiment: Experiment, use_eager_fetch: bool) -> None:
    exp_data_pipe = ExperimentDataPipe(
        soma_experiment,
        measurement_name="RNA",
        X_name="raw",
        obs_column_names=["label"],
        batch_size=3,
        return_sparse_X=True,
        shuffle=False,
        use_eager_fetch=use_eager_fetch,
    )
    batch_iter = iter(exp_data_pipe)

    batch = next(batch_iter)
    assert isinstance(batch[1], Tensor)
    assert batch[0].to_dense().tolist() == [[0, 1, 0], [1, 0, 1], [0, 1, 0]]


@pytest.mark.experimental
# noinspection PyTestParametrized,DuplicatedCode
@pytest.mark.parametrize(
    "obs_range,var_range,X_value_gen,use_eager_fetch",
    [(10, 1, pytorch_x_value_gen, use_eager_fetch) for use_eager_fetch in (True, False)],
)
def test_batching__partial_soma_batches_are_concatenated(soma_experiment: Experiment, use_eager_fetch: bool) -> None:
    exp_data_pipe = ExperimentDataPipe(
        soma_experiment,
        measurement_name="RNA",
        X_name="raw",
        obs_column_names=["label"],
        batch_size=3,
        # set SOMA batch read size such that PyTorch batches will span the tail and head of two SOMA batches
        soma_chunk_size=4,
        use_eager_fetch=use_eager_fetch,
    )

    full_result = list(exp_data_pipe)

    assert [len(batch[0]) for batch in full_result] == [3, 3, 3, 1]


@pytest.mark.experimental
# noinspection PyTestParametrized
@pytest.mark.parametrize("obs_range,var_range,X_value_gen", [(3, 3, pytorch_x_value_gen)])
def test_default_encoders_implicit(soma_experiment: Experiment) -> None:
    exp_data_pipe = ExperimentDataPipe(
        soma_experiment,
        measurement_name="RNA",
        X_name="raw",
        obs_column_names=["label"],
        shuffle=False,
        batch_size=3,
    )
    batch_iter = iter(exp_data_pipe)

    batch = next(batch_iter)
    assert isinstance(batch[1], Tensor)
    assert batch[0].to_dense().tolist() == [[0, 1, 0], [1, 0, 1], [0, 1, 0]]

    labels_encoded = batch[1]

    labels_decoded = exp_data_pipe.obs_encoders["label"].inverse_transform(labels_encoded)
    assert labels_decoded.tolist() == ["0", "1", "2"]  # type: ignore


@pytest.mark.experimental
# noinspection PyTestParametrized
@pytest.mark.parametrize("obs_range,var_range,X_value_gen", [(3, 3, pytorch_x_value_gen)])
def test_default_encoders_explicit(soma_experiment: Experiment) -> None:
    exp_data_pipe = ExperimentDataPipe(
        soma_experiment,
        measurement_name="RNA",
        X_name="raw",
        encoders=[LabelEncoder("label")],
        shuffle=False,
        batch_size=3,
    )
    batch_iter = iter(exp_data_pipe)

    batch = next(batch_iter)
    assert isinstance(batch[1], Tensor)

    labels_encoded = batch[1]

    labels_decoded = exp_data_pipe.obs_encoders["label"].inverse_transform(labels_encoded)
    assert labels_decoded.tolist() == ["0", "1", "2"]  # type: ignore


@pytest.mark.experimental
# noinspection PyTestParametrized
@pytest.mark.parametrize("obs_range,var_range,X_value_gen", [(3, 3, pytorch_x_value_gen)])
def test_batch_encoder(soma_experiment: Experiment) -> None:
    exp_data_pipe = ExperimentDataPipe(
        soma_experiment,
        measurement_name="RNA",
        X_name="raw",
        encoders=[BatchEncoder(["label", "label2"])],
        shuffle=False,
        batch_size=3,
    )
    batch_iter = iter(exp_data_pipe)

    batch = next(batch_iter)
    assert isinstance(batch[1], Tensor)

    labels_encoded = batch[1]

    labels_decoded = exp_data_pipe.obs_encoders["batch"].inverse_transform(labels_encoded)
    assert labels_decoded.tolist() == ["0c", "1c", "2c"]  # type: ignore


@pytest.mark.experimental
# noinspection PyTestParametrized
@pytest.mark.parametrize("obs_range,var_range,X_value_gen", [(3, 3, pytorch_x_value_gen)])
def test_custom_encoders_fail_if_duplicate(soma_experiment: Experiment) -> None:
    with pytest.raises(ValueError):
        ExperimentDataPipe(
            soma_experiment,
            measurement_name="RNA",
            X_name="raw",
            encoders=[LabelEncoder("label"), LabelEncoder("label")],
            shuffle=False,
            batch_size=3,
        )


@pytest.mark.experimental
# noinspection PyTestParametrized
@pytest.mark.parametrize("obs_range,var_range,X_value_gen", [(3, 3, pytorch_x_value_gen)])
def test_custom_encoders_fail_if_columns_defined(soma_experiment: Experiment) -> None:
    with pytest.raises(ValueError, match="Cannot specify both `obs_column_names` and `encoders`"):
        ExperimentDataPipe(
            soma_experiment,
            measurement_name="RNA",
            X_name="raw",
            obs_column_names=["label"],
            encoders=[LabelEncoder("label")],
            shuffle=False,
            batch_size=3,
        )


@pytest.mark.experimental
# noinspection PyTestParametrized
@pytest.mark.parametrize("obs_range,var_range,X_value_gen", [(6, 3, pytorch_x_value_gen)])
def test_multiprocessing__returns_full_result(soma_experiment: Experiment) -> None:
    """Tests the ExperimentDataPipe provides all data, as collected from multiple processes that are managed by a
    PyTorch DataLoader with multiple workers configured."""

    dp = ExperimentDataPipe(
        soma_experiment,
        measurement_name="RNA",
        X_name="raw",
        obs_column_names=["label"],
        soma_chunk_size=3,  # two chunks, one per worker
    )
    # Note we're testing the ExperimentDataPipe via a DataLoader, since this is what sets up the multiprocessing
    dl = experiment_dataloader(dp, num_workers=2)

    full_result = list(iter(dl))

    soma_joinids = [t[1][0].item() for t in full_result]
    assert sorted(soma_joinids) == list(range(6))


@pytest.mark.experimental
# noinspection PyTestParametrized
@pytest.mark.parametrize("obs_range,var_range,X_value_gen", [(6, 3, pytorch_x_value_gen)])
def test_distributed__returns_data_partition_for_rank(
    soma_experiment: Experiment,
) -> None:
    """Tests pytorch._partition_obs_joinids() behavior in a simulated PyTorch distributed processing mode,
    using mocks to avoid having to do real PyTorch distributed setup."""

    with (
        patch("cellxgene_census.experimental.ml.pytorch.dist.is_initialized") as mock_dist_is_initialized,
        patch("cellxgene_census.experimental.ml.pytorch.dist.get_rank") as mock_dist_get_rank,
        patch("cellxgene_census.experimental.ml.pytorch.dist.get_world_size") as mock_dist_get_world_size,
    ):
        mock_dist_is_initialized.return_value = True
        mock_dist_get_rank.return_value = 1
        mock_dist_get_world_size.return_value = 3

        dp = ExperimentDataPipe(
            soma_experiment,
            measurement_name="RNA",
            X_name="raw",
            encoders=[LabelEncoder("soma_joinid"), LabelEncoder("label")],
            soma_chunk_size=2,
            shuffle=False,
        )
        full_result = list(iter(dp))

        soma_joinids = [t[1][0].item() for t in full_result]

        # Of the 6 obs rows, the PyTorch process of rank 1 should get [2, 3]
        # (rank 0 gets [0, 1], rank 2 gets [4, 5])
        assert sorted(soma_joinids) == [2, 3]


@pytest.mark.experimental
# noinspection PyTestParametrized
@pytest.mark.parametrize("obs_range,var_range,X_value_gen", [(12, 3, pytorch_x_value_gen)])
def test_distributed_and_multiprocessing__returns_data_partition_for_rank(
    soma_experiment: Experiment,
) -> None:
    """Tests pytorch._partition_obs_joinids() behavior in a simulated PyTorch distributed processing mode and
    DataLoader multiprocessing mode, using mocks to avoid having to do distributed pytorch
    setup or real DataLoader multiprocessing."""

    with (
        patch("torch.utils.data.get_worker_info") as mock_get_worker_info,
        patch("cellxgene_census.experimental.ml.pytorch.dist.is_initialized") as mock_dist_is_initialized,
        patch("cellxgene_census.experimental.ml.pytorch.dist.get_rank") as mock_dist_get_rank,
        patch("cellxgene_census.experimental.ml.pytorch.dist.get_world_size") as mock_dist_get_world_size,
    ):
        mock_get_worker_info.return_value = WorkerInfo(id=1, num_workers=2, seed=1234)
        mock_dist_is_initialized.return_value = True
        mock_dist_get_rank.return_value = 1
        mock_dist_get_world_size.return_value = 3

        dp = ExperimentDataPipe(
            soma_experiment,
            measurement_name="RNA",
            X_name="raw",
            encoders=[LabelEncoder("soma_joinid"), LabelEncoder("label")],
            soma_chunk_size=2,
            shuffle=False,
        )

        full_result = list(iter(dp))

        soma_joinids = [t[1][0].item() for t in full_result]

        # Of the 12 obs rows, the PyTorch process of rank 1 should get [4..7], and then within that partition,
        # the 2nd DataLoader process should get the second half of the rank's partition, which is just [6, 7]
        # (rank 0 gets [0..3], rank 2 gets [8..11])
        assert sorted(soma_joinids) == [6, 7]


@pytest.mark.experimental
# noinspection PyTestParametrized,DuplicatedCode
@pytest.mark.parametrize(
    "obs_range,var_range,X_value_gen,use_eager_fetch",
    [(3, 3, pytorch_x_value_gen, use_eager_fetch) for use_eager_fetch in (True, False)],
)
def test_experiment_dataloader__non_batched(soma_experiment: Experiment, use_eager_fetch: bool) -> None:
    dp = ExperimentDataPipe(
        soma_experiment,
        measurement_name="RNA",
        X_name="raw",
        encoders=[LabelEncoder("soma_joinid"), LabelEncoder("label")],
        shuffle=False,
        use_eager_fetch=use_eager_fetch,
    )
    dl = experiment_dataloader(dp)
    torch_data = [row for row in dl]  # noqa: C416

    row = torch_data[0]
    assert row[0].to_dense().tolist() == [0, 1, 0]
    assert row[1].tolist() == [0, 0]


@pytest.mark.experimental
# noinspection PyTestParametrized,DuplicatedCode
@pytest.mark.parametrize(
    "obs_range,var_range,X_value_gen,use_eager_fetch",
    [(6, 3, pytorch_x_value_gen, use_eager_fetch) for use_eager_fetch in (True, False)],
)
def test_experiment_dataloader__batched(soma_experiment: Experiment, use_eager_fetch: bool) -> None:
    dp = ExperimentDataPipe(
        soma_experiment,
        measurement_name="RNA",
        X_name="raw",
        encoders=[LabelEncoder("soma_joinid"), LabelEncoder("label")],
        batch_size=3,
        shuffle=False,
        use_eager_fetch=use_eager_fetch,
    )
    dl = experiment_dataloader(dp)
    torch_data = [row for row in dl]  # noqa: C416

    batch = torch_data[0]
    assert batch[0].to_dense().tolist() == [[0, 1, 0], [1, 0, 1], [0, 1, 0]]
    assert batch[1].tolist() == [[0, 0], [1, 1], [2, 2]]


@pytest.mark.experimental
# noinspection PyTestParametrized,DuplicatedCode
@pytest.mark.parametrize(
    "obs_range,var_range,X_value_gen,use_eager_fetch",
    [(10, 3, pytorch_x_value_gen, use_eager_fetch) for use_eager_fetch in (True, False)],
)
def test_experiment_dataloader__batched_length(soma_experiment: Experiment, use_eager_fetch: bool) -> None:
    dp = ExperimentDataPipe(
        soma_experiment,
        measurement_name="RNA",
        X_name="raw",
        obs_column_names=["label"],
        batch_size=3,
        shuffle=False,
        use_eager_fetch=use_eager_fetch,
    )
    dl = experiment_dataloader(dp)
    assert len(dl) == len(list(dl))


@pytest.mark.experimental
# noinspection PyTestParametrized,DuplicatedCode
@pytest.mark.parametrize(
    "obs_range,var_range,X_value_gen,use_eager_fetch",
    [(6, 3, pytorch_x_value_gen, use_eager_fetch) for use_eager_fetch in (True, False)],
)
def test__X_tensor_dtype_matches_X_matrix(soma_experiment: Experiment, use_eager_fetch: bool) -> None:
    dp = ExperimentDataPipe(
        soma_experiment,
        measurement_name="RNA",
        X_name="raw",
        obs_column_names=["label"],
        batch_size=3,
        use_eager_fetch=use_eager_fetch,
    )
    torch_data = next(iter(dp))

    assert torch_data[0].dtype == float32


@pytest.mark.experimental
# noinspection PyTestParametrized,DuplicatedCode
@pytest.mark.parametrize("obs_range,var_range,X_value_gen", [(10, 1, pytorch_x_value_gen)])
def test__pytorch_splitting(soma_experiment: Experiment) -> None:
    dp = ExperimentDataPipe(
        soma_experiment,
        measurement_name="RNA",
        X_name="raw",
        obs_column_names=["label"],
    )
    dp_train, dp_test = dp.random_split(weights={"train": 0.7, "test": 0.3}, seed=1234)
    dl = experiment_dataloader(dp_train)

    all_rows = list(iter(dl))
    assert len(all_rows) == 7


@pytest.mark.experimental
# noinspection PyTestParametrized,DuplicatedCode
@pytest.mark.parametrize("obs_range,var_range,X_value_gen", [(16, 1, pytorch_seq_x_value_gen)])
def test__shuffle(soma_experiment: Experiment) -> None:
    dp = ExperimentDataPipe(
        soma_experiment,
        measurement_name="RNA",
        X_name="raw",
        encoders=[LabelEncoder("soma_joinid"), LabelEncoder("label")],
        shuffle=True,
    )

    all_rows = list(iter(dp))

    soma_joinids = [row[1][0].item() for row in all_rows]
    X_values = [row[0][0].item() for row in all_rows]

    # same elements
    assert set(soma_joinids) == set(range(16))
    # not ordered! (...with a `1/16!` probability of being ordered)
    assert soma_joinids != list(range(16))
    # randomizes X in same order as obs
    # note: X values were explicitly set to match obs_joinids to allow for this simple assertion
    assert X_values == soma_joinids


@pytest.mark.experimental
@pytest.mark.skip(reason="Not implemented")
def test_experiment_dataloader__multiprocess_sparse_matrix__fails() -> None:
    pass


@pytest.mark.experimental
@pytest.mark.skip(reason="Not implemented")
def test_experiment_dataloader__multiprocess_dense_matrix__ok() -> None:
    pass


@pytest.mark.experimental
def test_experiment_dataloader__unsupported_params__fails() -> None:
    with patch("cellxgene_census.experimental.ml.pytorch.ExperimentDataPipe") as dummy_exp_data_pipe:
        with pytest.raises(ValueError):
            experiment_dataloader(dummy_exp_data_pipe, shuffle=True)
        with pytest.raises(ValueError):
            experiment_dataloader(dummy_exp_data_pipe, batch_size=3)
        with pytest.raises(ValueError):
            experiment_dataloader(dummy_exp_data_pipe, batch_sampler=[])
        with pytest.raises(ValueError):
            experiment_dataloader(dummy_exp_data_pipe, sampler=[])
        with pytest.raises(ValueError):
            experiment_dataloader(dummy_exp_data_pipe, collate_fn=lambda x: x)


@pytest.mark.experimental
def test_list_split() -> None:
    data = list(range(10))
    chunks = list_split(data, 3)
    assert len(chunks) == 4
    assert len(chunks[0]) == 3
    assert len(chunks[1]) == 3
    assert len(chunks[2]) == 3
    assert len(chunks[3]) == 1



# Section: conf

# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'cellxgene-census'
copyright = '2022, Chan Zuckerberg Initiative'
author = 'Chan Zuckerberg Initiative'

from packaging.version import Version
import git

repo = git.Repo(search_parent_directories=True)
all_versions = sorted([Version(t.name) for t in repo.tags])

version = str(all_versions[-1])

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc', 
    "nbsphinx", 
    "sphinx.ext.intersphinx", 
    'sphinx.ext.napoleon', 
    'sphinx.ext.autosummary', 
    'myst_parser'
]

autosummary_generate = True

napoleon_custom_sections = ["Lifecycle"]

tiledb_version = "latest"

intersphinx_mapping = {
    "tiledbsoma-py": (
        "https://tiledbsoma.readthedocs.io/en/%s/"
        % tiledb_version,
        None,
    ),
    'python': ('https://docs.python.org/3', None),
    'numpy': ('http://docs.scipy.org/doc/numpy', None),
    'scipy': ('http://docs.scipy.org/doc/scipy/reference', None),
    'anndata': ('https://anndata.readthedocs.io/en/latest/', None),
    'pandas': ('https://pandas.pydata.org/docs/', None),
    'scanpy': ('https://scanpy.readthedocs.io/en/stable/', None),
    'torch': ('https://pytorch.org/docs/stable/', None),
    'torchdata': ('https://pytorch.org/data/beta/', None),
    'sklearn': ('http://scikit-learn.org/stable', None),
}

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

source_suffix = ['.rst', '.md']



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

# html_theme = 'alabaster'
# Inject custom css files in `/_static/css/*`
html_static_path = ['_static']

html_theme = "sphinx_rtd_theme"

html_js_files = [
    ('https://plausible.io/js/script.js', {"data-domain": "chanzuckerberg.github.io/cellxgene-census", "defer": "defer"}),
]

def setup(app):
    app.add_css_file("css/custom.css")



# Section: notebooks-api_demo-census_dataset_presence

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




# Section: notebooks-api_demo-census_query_extract

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




# Section: notebooks-experimental-mean_variance

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




# Section: cellxgene_census-tests-test_get_anndata

from typing import Any, Literal

import numpy as np
import pandas as pd
import pytest
import tiledbsoma as soma

import cellxgene_census


@pytest.mark.live_corpus
def test_get_anndata_value_filter(census: soma.Collection) -> None:
    ad = cellxgene_census.get_anndata(
        census,
        organism="Mus musculus",
        obs_value_filter="tissue_general == 'vasculature'",
        var_value_filter="feature_name in ['Gm53058', '0610010K14Rik']",
        obs_column_names=[
            "soma_joinid",
            "cell_type",
            "tissue",
            "tissue_general",
            "assay",
        ],
        var_column_names=["soma_joinid", "feature_id", "feature_name", "feature_length"],
    )

    assert ad is not None
    assert ad.n_vars == 2
    assert ad.n_obs > 0
    assert (ad.obs.tissue_general == "vasculature").all()
    assert set(ad.var.feature_name) == {"Gm53058", "0610010K14Rik"}


@pytest.mark.live_corpus
def test_get_anndata_coords(census: soma.Collection) -> None:
    ad = cellxgene_census.get_anndata(
        census,
        organism="Mus musculus",
        obs_coords=slice(1000),
        var_coords=slice(2000),
    )

    assert ad is not None
    assert ad.n_vars == 2001
    assert ad.n_obs == 1001


@pytest.mark.live_corpus
def test_get_anndata_allows_missing_obs_or_var_filter(census: soma.Collection) -> None:
    # This test is slightly sensitive to the live data, in that it assumes the
    # existance of certain cell tissue labels and gene feature ids.
    mouse = census["census_data"]["mus_musculus"]

    adata = cellxgene_census.get_anndata(census, organism="Mus musculus", obs_value_filter="tissue == 'aorta'")
    assert adata.n_obs == len(mouse.obs.read(value_filter="tissue == 'aorta'", column_names=["soma_joinid"]).concat())
    assert adata.n_vars == len(mouse.ms["RNA"].var.read(column_names=["soma_joinid"]).concat())

    adata = cellxgene_census.get_anndata(
        census,
        organism="Mus musculus",
        obs_coords=slice(10000),
        var_value_filter="feature_id == 'ENSMUSG00000069581'",
    )
    assert adata.n_obs == 10001
    assert adata.n_vars == 1


@pytest.mark.live_corpus
@pytest.mark.parametrize("layer", ["raw", "normalized"])
def test_get_anndata_x_layer(census: soma.Collection, layer: str) -> None:
    ad = cellxgene_census.get_anndata(
        census,
        organism="Homo sapiens",
        X_name=layer,
        obs_coords=slice(100),
        var_coords=slice(200),
    )

    assert ad.X.shape == (101, 201)
    assert len(ad.layers) == 0


@pytest.mark.live_corpus
@pytest.mark.parametrize("layers", [["raw", "normalized"], ["normalized", "raw"]])
def test_get_anndata_two_layers(census: soma.Collection, layers: list[str]) -> None:
    ad_primary_layer_in_X = cellxgene_census.get_anndata(
        census,
        organism="Homo sapiens",
        X_name=layers[0],
        obs_coords=slice(100),
        var_coords=slice(200),
    )

    ad_secondary_layer_in_X = cellxgene_census.get_anndata(
        census,
        organism="Homo sapiens",
        X_name=layers[1],
        obs_coords=slice(100),
        var_coords=slice(200),
    )

    ad_multiple_layers = cellxgene_census.get_anndata(
        census,
        organism="Homo sapiens",
        X_name=layers[0],
        X_layers=[layers[1]],
        obs_coords=slice(100),
        var_coords=slice(200),
    )

    assert layers[1] in ad_multiple_layers.layers
    assert ad_multiple_layers.X.shape == (101, 201)
    assert ad_multiple_layers.layers[layers[1]].shape == (101, 201)

    # Assert that matrices of multilayer anndata are equal to one-layer-at-time anndatas
    assert np.array_equal(ad_multiple_layers.X.data, ad_primary_layer_in_X.X.data)
    assert np.array_equal(ad_multiple_layers.layers[layers[1]].data, ad_secondary_layer_in_X.X.data)


@pytest.mark.live_corpus
def test_get_anndata_wrong_layer_names(census: soma.Collection) -> None:
    with pytest.raises(ValueError) as raise_info:
        cellxgene_census.get_anndata(
            census,
            organism="Homo sapiens",
            X_name="this_layer_name_is_bad",
            obs_coords=slice(100),
            var_coords=slice(200),
        )

        assert raise_info.value.args[0] == "Unknown X layer name"

    with pytest.raises(ValueError) as raise_info:
        cellxgene_census.get_anndata(
            census,
            organism="Homo sapiens",
            X_name="raw",
            X_layers=["this_layer_name_is_bad"],
            obs_coords=slice(100),
            var_coords=slice(200),
        )

        assert raise_info.value.args[0] == "Unknown X layer name"


@pytest.mark.live_corpus
@pytest.mark.parametrize("obsm_layer", ["scvi", "geneformer"])
def test_get_anndata_obsm_one_layer(dec_lts_census: soma.Collection, obsm_layer: str) -> None:
    # NOTE: This test only works on the 2023-12-15 LTS Census, since in newer releases
    # the embeddings aren't distributed via the `obsm_layer` parameter.
    ad = cellxgene_census.get_anndata(
        dec_lts_census,
        organism="Homo sapiens",
        X_name="raw",
        obs_coords=slice(100),
        var_coords=slice(200),
        obsm_layers=[obsm_layer],
    )

    assert len(ad.obsm.keys()) == 1
    assert obsm_layer in ad.obsm.keys()
    assert ad.obsm[obsm_layer].shape[0] == 101


@pytest.mark.live_corpus
@pytest.mark.parametrize("obsm_layers", [["scvi", "geneformer"]])
def test_get_anndata_obsm_two_layers(dec_lts_census: soma.Collection, obsm_layers: list[str]) -> None:
    # NOTE: This test only works on the 2023-12-15 LTS Census, since in newer releases
    # the embeddings aren't distributed via the `obsm_layer` parameter.
    ad = cellxgene_census.get_anndata(
        dec_lts_census,
        organism="Homo sapiens",
        X_name="raw",
        obs_coords=slice(100),
        var_coords=slice(200),
        obsm_layers=obsm_layers,
    )

    assert len(ad.obsm.keys()) == 2
    for obsm_layer in obsm_layers:
        assert obsm_layer in ad.obsm.keys()
        assert ad.obsm[obsm_layer].shape[0] == 101


@pytest.mark.live_corpus
@pytest.mark.parametrize("obs_embeddings", [["scvi", "geneformer"]])
def test_get_anndata_obs_embeddings(lts_census: soma.Collection, obs_embeddings: list[str]) -> None:
    # NOTE: when the next LTS gets released (>2023-12-15), embeddings may or may not be available,
    # so this test could require adjustments.
    ad = cellxgene_census.get_anndata(
        lts_census,
        organism="Homo sapiens",
        X_name="raw",
        obs_coords=slice(100),
        var_coords=slice(200),
        obs_embeddings=obs_embeddings,
    )

    assert len(ad.obsm.keys()) == 2
    assert len(ad.varm.keys()) == 0
    for obsm_layer in obs_embeddings:
        assert obsm_layer in ad.obsm.keys()
        assert ad.obsm[obsm_layer].shape[0] == 101


@pytest.mark.live_corpus
@pytest.mark.parametrize("var_embeddings", [["nmf"]])
def test_get_anndata_var_embeddings(dec_lts_census: soma.Collection, var_embeddings: list[str]) -> None:
    # NOTE: this test only works on the 2023-12-15 LTS Census, since var embeddings
    # aren't available in the newer releases.

    ad = cellxgene_census.get_anndata(
        dec_lts_census,
        organism="Homo sapiens",
        X_name="raw",
        obs_coords=slice(100),
        var_coords=slice(200),
        var_embeddings=var_embeddings,
    )

    assert len(ad.obsm.keys()) == 0
    assert len(ad.varm.keys()) == 1
    for varm_layers in var_embeddings:
        assert varm_layers in ad.varm.keys()
        assert ad.varm[varm_layers].shape[0] == 201


@pytest.mark.live_corpus
def test_get_anndata_obsm_layers_and_add_obs_embedding_fails(lts_census: soma.Collection) -> None:
    """Fails if both `obsm_layers` and `obs_embeddings` are specified."""
    with pytest.raises(ValueError):
        cellxgene_census.get_anndata(
            lts_census,
            organism="Homo sapiens",
            X_name="raw",
            obs_coords=slice(100),
            var_coords=slice(200),
            obsm_layers=["scvi"],
            obs_embeddings=["scvi"],
        )


@pytest.mark.live_corpus
def test_deprecated_column_api(census: soma.Collection) -> None:
    """Testing for previous `column_names` argument.

    See: https://github.com/chanzuckerberg/cellxgene-census/issues/1035
    """
    ad_curr = cellxgene_census.get_anndata(
        census,
        organism="Mus musculus",
        obs_value_filter="tissue_general == 'vasculature'",
        var_value_filter="feature_name in ['Gm53058', '0610010K14Rik']",
        obs_column_names=[
            "soma_joinid",
            "cell_type",
            "tissue",
            "tissue_general",
            "assay",
        ],
        var_column_names=["soma_joinid", "feature_id", "feature_name", "feature_length"],
    )
    with pytest.warns(FutureWarning):
        ad_prev = cellxgene_census.get_anndata(
            census,
            organism="Mus musculus",
            obs_value_filter="tissue_general == 'vasculature'",
            var_value_filter="feature_name in ['Gm53058', '0610010K14Rik']",
            column_names={
                "obs": [
                    "soma_joinid",
                    "cell_type",
                    "tissue",
                    "tissue_general",
                    "assay",
                ],
                "var": ["soma_joinid", "feature_id", "feature_name", "feature_length"],
            },
        )
    with pytest.raises(
        ValueError, match=r"Both the deprecated 'column_names' argument and its replacements were used."
    ):
        cellxgene_census.get_anndata(
            census,
            organism="Mus musculus",
            obs_value_filter="tissue_general == 'vasculature'",
            var_value_filter="feature_name in ['Gm53058', '0610010K14Rik']",
            obs_column_names=[
                "soma_joinid",
                "cell_type",
            ],
            column_names={
                "obs": [
                    "soma_joinid",
                    "cell_type",
                ],
            },
        )
    pd.testing.assert_frame_equal(ad_curr.obs, ad_prev.obs)
    pd.testing.assert_frame_equal(ad_curr.var, ad_prev.var)


def _map_to_get_anndata_args(query: dict[str, Any], axis: Literal["obs", "var"]) -> dict[str, Any]:
    """Helper to map arguments of get_obs/ get_var to get_anndata."""
    result = {}
    if "coords" in query:
        result[f"{axis}_coords"] = query["coords"]
    if "value_filter" in query:
        result[f"{axis}_value_filter"] = query["value_filter"]
    if "column_names" in query:
        result["column_names"] = {axis: query["column_names"]}
    return result


@pytest.mark.live_corpus
@pytest.mark.parametrize(
    "query",
    [
        pytest.param(
            {
                "coords": slice(100),
                "column_names": [
                    "soma_joinid",
                    "cell_type",
                    "tissue",
                    "tissue_general",
                    "assay",
                ],
            },
            id="coords+column-names",
        ),
        pytest.param({"coords": slice(100, 300)}, id="coords"),
        pytest.param({"value_filter": "tissue_general == 'vasculature'"}, id="value_filter"),
    ],
)
def test_get_obs(lts_census: soma.Collection, query: dict[str, Any]) -> None:
    adata_obs = cellxgene_census.get_anndata(
        lts_census, organism="Mus musculus", **_map_to_get_anndata_args(query, "obs")
    ).obs
    only_obs = cellxgene_census.get_obs(lts_census, "Mus musculus", **query)
    # account for a difference:
    only_obs.index = only_obs.index.astype(str)

    pd.testing.assert_frame_equal(adata_obs, only_obs)


@pytest.mark.live_corpus
@pytest.mark.parametrize(
    "query",
    [
        pytest.param(
            {
                "coords": slice(100),
                "column_names": ["soma_joinid", "feature_id", "feature_name", "feature_length"],
            },
            id="coords+column-names",
        ),
        pytest.param({"coords": slice(100, 300)}, id="coords"),
        pytest.param({"value_filter": "feature_name in ['Gm53058', '0610010K14Rik']"}, id="value_filter"),
    ],
)
def test_get_var(lts_census: soma.Collection, query: dict[str, Any]) -> None:
    adata_var = cellxgene_census.get_anndata(
        lts_census, organism="Mus musculus", obs_coords=slice(0), **_map_to_get_anndata_args(query, "var")
    ).var
    only_var = cellxgene_census.get_var(lts_census, "Mus musculus", **query)
    # AnnData instantiation converts the index to string, so we match that behaviour for comparisons sake
    only_var.index = only_var.index.astype(str)

    pd.testing.assert_frame_equal(adata_var, only_var)



# Section: cellxgene_census-src-cellxgene_census-experimental-pp-__init__

"""API to facilitate preprocessing of SOMA datasets."""

from ._highly_variable_genes import get_highly_variable_genes, highly_variable_genes
from ._stats import mean_variance

__all__ = [
    "get_highly_variable_genes",
    "highly_variable_genes",
    "mean_variance",
]



# Section: notebooks-analysis_demo-comp_bio_normalizing_full_gene_sequencing

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



# Section: notebooks-analysis_demo-comp_bio_embedding_exploration

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




# Section: cellxgene_census-tests-__init__




# Section: cellxgene_census-tests-test_acceptance

"""
Acceptance tests for the Census.

NOTE: those marked `expensive` are not run in the CI as they are, well, expensive...

Several of them will not run to completion except on VERY large hosts.

Intended use:  periodically do a manual run, including the expensive tests, on an
appropriately large host.

See README.md for historical data.
"""

from collections.abc import Iterator
from typing import Any

import pyarrow as pa
import pytest
import tiledbsoma as soma

import cellxgene_census
from cellxgene_census._open import DEFAULT_TILEDB_CONFIGURATION


def make_context(census_version: str, config: dict[str, Any] | None = None) -> soma.SOMATileDBContext:
    config = config or {}
    version = cellxgene_census.get_census_version_description(census_version)
    s3_region = version["soma"].get("s3_region", "us-west-2")
    config.update({"vfs.s3.region": s3_region})
    config.update({"vfs.s3.no_sign_request": "true"})
    return soma.options.SOMATileDBContext().replace(**{"tiledb_config": config})


@pytest.mark.live_corpus
@pytest.mark.parametrize("organism", ["homo_sapiens", "mus_musculus"])
def test_load_axes(organism: str) -> None:
    """Verify axes can be loaded into a Pandas DataFrame"""
    census = cellxgene_census.open_soma(census_version="latest")

    # use subset of columns for speed
    obs_df = (
        census["census_data"][organism]
        .obs.read(column_names=["soma_joinid", "cell_type", "tissue"])
        .concat()
        .to_pandas()
    )
    assert len(obs_df)
    del obs_df

    var_df = census["census_data"][organism].ms["RNA"].var.read().concat().to_pandas()
    assert len(var_df)
    del var_df


def table_iter_is_ok(tbl_iter: Iterator[pa.Table], stop_after: int | None = 2) -> bool:
    """
    Utility that verifies that the value is an iterator of pa.Table.

    Will only call __next__ as many times as the `stop_after` param specifies,
    or will read until end of iteration of it is None.
    """
    assert isinstance(tbl_iter, Iterator)
    for n, tbl in enumerate(tbl_iter):
        # keep things speedy by quitting early if stop_after specified
        if stop_after is not None and n > stop_after:
            break
        assert isinstance(tbl, pa.Table)
        assert len(tbl)

    return True


@pytest.mark.live_corpus
@pytest.mark.parametrize("organism", ["homo_sapiens", "mus_musculus"])
@pytest.mark.parametrize(
    ("stop_after", "ctx_config"),
    [
        pytest.param(2, None),
        pytest.param(None, DEFAULT_TILEDB_CONFIGURATION, marks=pytest.mark.expensive),
    ],
)
def test_incremental_read_obs(organism: str, stop_after: int | None, ctx_config: dict[str, Any] | None) -> None:
    """Verify that obs, var and X[raw] can be read incrementally, i.e., in chunks"""

    # ctx_config=None open census with a small (default) TileDB buffer size, which reduces
    # memory use, and makes it feasible to run in a GHA.
    ctx_config = ctx_config or {}
    context = make_context("latest", ctx_config)
    with cellxgene_census.open_soma(census_version="latest", context=context) as census:
        assert table_iter_is_ok(
            census["census_data"][organism].obs.read(column_names=["soma_joinid", "tissue"]),
            stop_after=stop_after,
        )


@pytest.mark.live_corpus
@pytest.mark.parametrize("organism", ["homo_sapiens", "mus_musculus"])
@pytest.mark.parametrize(
    ("stop_after", "ctx_config"),
    [
        pytest.param(2, None),
        pytest.param(None, DEFAULT_TILEDB_CONFIGURATION, marks=pytest.mark.expensive),
    ],
)
def test_incremental_read_var(organism: str, stop_after: int | None, ctx_config: dict[str, Any] | None) -> None:
    """Verify that var can be read incrementally, i.e., in chunks"""

    # ctx_config=None open census with a small (default) TileDB buffer size, which reduces
    # memory use, and makes it feasible to run in a GHA.
    ctx_config = ctx_config or {}
    context = make_context("latest", ctx_config)
    with cellxgene_census.open_soma(census_version="latest", context=context) as census:
        assert table_iter_is_ok(
            census["census_data"][organism].ms["RNA"].var.read(column_names=["soma_joinid", "feature_id"]),
            stop_after=stop_after,
        )


@pytest.mark.live_corpus
@pytest.mark.parametrize("organism", ["homo_sapiens", "mus_musculus"])
@pytest.mark.parametrize(
    ("stop_after", "ctx_config", "coords"),
    [
        pytest.param(2, None, (None, None)),
        pytest.param(
            None,
            DEFAULT_TILEDB_CONFIGURATION,
            (slice(0, 500_000), slice(None)),
            marks=pytest.mark.expensive,
        ),
        pytest.param(
            None,
            {"soma.init_buffer_bytes": 4 * 1024**3},
            (slice(0, 500_000), slice(None)),
            marks=pytest.mark.expensive,
        ),
        pytest.param(
            None,
            {"soma.init_buffer_bytes": 4 * 1024**3},
            (slice(0, 500_000), slice(0, 1_000)),
            marks=pytest.mark.expensive,
        ),
    ],
)
def test_incremental_read_X(
    organism: str,
    stop_after: int | None,
    ctx_config: dict[str, Any] | None,
    coords: tuple[slice, slice] | None,
) -> None:
    """Verify that obs, var and X[raw] can be read incrementally, i.e., in chunks"""

    ctx_config = ctx_config or {}
    context = make_context("latest", ctx_config)
    with cellxgene_census.open_soma(census_version="latest", context=context) as census:
        assert table_iter_is_ok(
            census["census_data"][organism].ms["RNA"].X["raw"].read(coords=coords).tables(),
            stop_after=stop_after,
        )


@pytest.mark.live_corpus
@pytest.mark.parametrize("organism", ["homo_sapiens", "mus_musculus"])
@pytest.mark.parametrize(
    "obs_value_filter",
    ["tissue=='aorta'", pytest.param("tissue=='brain'", marks=pytest.mark.expensive)],
)
@pytest.mark.parametrize("stop_after", [2, pytest.param(None, marks=pytest.mark.expensive)])
def test_incremental_query(organism: str, obs_value_filter: str, stop_after: int | None) -> None:
    """Verify incremental read of query result."""
    # use default TileDB configuration
    with cellxgene_census.open_soma(census_version="latest") as census:
        with census["census_data"][organism].axis_query(
            measurement_name="RNA",
            obs_query=soma.AxisQuery(value_filter=obs_value_filter),
        ) as query:
            assert table_iter_is_ok(query.obs(), stop_after=stop_after)
            assert table_iter_is_ok(query.var(), stop_after=stop_after)
            assert table_iter_is_ok(query.X("raw").tables(), stop_after=stop_after)


@pytest.mark.live_corpus
@pytest.mark.parametrize("organism", ["homo_sapiens", "mus_musculus"])
@pytest.mark.parametrize(
    ("obs_value_filter", "obs_coords", "ctx_config"),
    [
        # small query, should be runable in CI
        pytest.param("tissue=='aorta'", None, DEFAULT_TILEDB_CONFIGURATION),
        # 10K cells, also small enough to run in CI
        pytest.param(None, slice(0, 10_000), DEFAULT_TILEDB_CONFIGURATION, id="First 10K cells"),
        # 100K cells, standard buffer size
        pytest.param(
            None,
            slice(0, 100_000),
            DEFAULT_TILEDB_CONFIGURATION,
            marks=pytest.mark.expensive,
            id="First 100K cells",
        ),
        # 250K cells, standard buffer size
        pytest.param(
            None,
            slice(0, 250_000),
            DEFAULT_TILEDB_CONFIGURATION,
            marks=pytest.mark.expensive,
            id="First 250K cells",
        ),
        # 500K cells, standard buffer size
        pytest.param(
            None,
            slice(0, 500_000),
            DEFAULT_TILEDB_CONFIGURATION,
            marks=pytest.mark.expensive,
            id="First 500K cells",
        ),
        # 750K cells, standard buffer size
        pytest.param(
            None,
            slice(0, 750_000),
            DEFAULT_TILEDB_CONFIGURATION,
            marks=pytest.mark.expensive,
            id="First 750K cells",
        ),
        # 1M cells, standard buffer size
        pytest.param(
            None,
            slice(0, 1_000_000),
            DEFAULT_TILEDB_CONFIGURATION,
            marks=pytest.mark.expensive,
            id="First 1M cells",
        ),
        # very common tissue, with standard buffer size
        pytest.param(
            "tissue_general=='brain'",
            slice(0, 1_000_000),
            DEFAULT_TILEDB_CONFIGURATION,
            marks=pytest.mark.expensive,
        ),
        # very common cell type, with standard buffer size
        pytest.param(
            "cell_type=='neuron'",
            slice(0, 1_000_000),
            DEFAULT_TILEDB_CONFIGURATION,
            marks=pytest.mark.expensive,
        ),
        # many primary cells, with big buffer size
        pytest.param(
            "is_primary_data==True",
            slice(0, 1_000_000),
            {"soma.init_buffer_bytes": 4 * 1024**3},
            marks=pytest.mark.expensive,
        ),
        #
        # a large enchilada, if not the whole thing, with a big buffer size
        pytest.param(
            None,
            slice(0, 1_000_000),
            {"soma.init_buffer_bytes": 4 * 1024**3},
            marks=pytest.mark.expensive,
        ),
    ],
)
def test_get_anndata(
    organism: str,
    obs_value_filter: str | None,
    obs_coords: slice | None,
    ctx_config: dict[str, Any] | None,
) -> None:
    """Verify query and read into AnnData"""
    ctx_config = ctx_config or {}
    context = make_context("latest", ctx_config)
    with cellxgene_census.open_soma(census_version="latest", context=context) as census:
        ad = cellxgene_census.get_anndata(census, organism, obs_value_filter=obs_value_filter, obs_coords=obs_coords)
        assert ad is not None

        # sanity checks
        with census["census_data"][organism].axis_query(
            measurement_name="RNA",
            obs_query=soma.AxisQuery(value_filter=obs_value_filter, coords=(obs_coords,)),
        ) as query:
            assert ad.n_obs == query.n_obs
            assert ad.n_vars == query.n_vars



# Section: cellxgene_census-tests-test_open

import os
import pathlib
import re
import time
from typing import TYPE_CHECKING
from unittest.mock import ANY, patch

import anndata
import numpy as np
import pytest
import requests_mock as rm
import tiledb
import tiledbsoma as soma

import cellxgene_census
from cellxgene_census import get_default_soma_context
from cellxgene_census._open import DEFAULT_TILEDB_CONFIGURATION
from cellxgene_census._release_directory import (
    CELL_CENSUS_MIRRORS_DIRECTORY_URL,
    CELL_CENSUS_RELEASE_DIRECTORY_URL,
    CensusLocator,
)

if TYPE_CHECKING:
    # You're not supposed to import this, but mypy demands it
    pass


@pytest.mark.live_corpus
def test_open_soma_stable() -> None:
    # There should _always_ be a 'stable'
    with cellxgene_census.open_soma(census_version="stable") as census:
        assert census is not None
        assert isinstance(census, soma.Collection)

    # and it should be the latest, until the first "stable" build is available
    with cellxgene_census.open_soma() as default_census:
        assert default_census.uri == census.uri
        for k, v in DEFAULT_TILEDB_CONFIGURATION.items():
            assert census.context.tiledb_config[k] == str(v)


@pytest.fixture(scope="module")
def latest_locator() -> CensusLocator:
    return cellxgene_census.get_census_version_description("latest")["soma"]


@pytest.mark.live_corpus
def test_open_soma_latest(latest_locator: CensusLocator) -> None:
    with cellxgene_census.open_soma(census_version="latest") as census:
        # There should _always_ be a 'latest'
        assert census is not None

        # It should always be a SOMA Collection
        assert isinstance(census, soma.Collection)

        # Verify that open_soma() actually opened "latest"
        assert census.uri == latest_locator["uri"]


@pytest.mark.live_corpus
def test_open_soma_with_customized_tiledb_config(latest_locator: CensusLocator) -> None:
    soma_init_buffer_bytes = "221000"
    tiledb_config = {
        "soma.init_buffer_bytes": soma_init_buffer_bytes,
        "vfs.s3.region": latest_locator.get("s3_region"),
    }
    with cellxgene_census.open_soma(uri=latest_locator["uri"], tiledb_config=tiledb_config) as census:
        assert census.uri == latest_locator["uri"]
        # Verify that user-provided custom config is passed through correctly
        assert census.context.tiledb_config["soma.init_buffer_bytes"] == soma_init_buffer_bytes


@pytest.mark.live_corpus
def test_open_soma_with_customized_plain_soma_context(
    latest_locator: CensusLocator,
) -> None:
    soma_init_buffer_bytes = "221000"
    timestamp_ms = int(time.time() * 1000) - 10  # don't use exactly current time, as that is the default
    cfg = {
        "timestamp": timestamp_ms,
        "tiledb_config": {
            "soma.init_buffer_bytes": soma_init_buffer_bytes,
            # The below settings are required to access the Census, but otherwise not material to the test.
            # By virtue of the Census opening successfully, we know these settings are being applied.
            "vfs.s3.region": latest_locator.get("s3_region"),
            "vfs.s3.no_sign_request": "true",
        },
    }
    context = soma.SOMATileDBContext().replace(**cfg)
    with cellxgene_census.open_soma(uri=latest_locator["uri"], context=context) as census:
        # Verify that the user-provided config settings are set correctly in the TileDB context object.
        assert census.context.tiledb_config["soma.init_buffer_bytes"] == soma_init_buffer_bytes
        assert census.context.timestamp_ms == timestamp_ms


@pytest.mark.live_corpus
def test_open_soma_with_customized_default_soma_context(
    latest_locator: CensusLocator,
) -> None:
    soma_init_buffer_bytes = "221000"

    timestamp_ms = int(time.time() * 1000) - 10  # don't use exactly current time, as that is the default
    custom_context = get_default_soma_context().replace(
        tiledb_config={"soma.init_buffer_bytes": soma_init_buffer_bytes},
        timestamp=timestamp_ms,
    )

    with cellxgene_census.open_soma(census_version="latest", context=custom_context) as census:
        # Verify the non-overriden soma context defaults are set correctly in the TileDB context object.
        assert census.context.tiledb_config["vfs.s3.no_sign_request"] == "true"
        assert census.context.tiledb_config["vfs.s3.region"] == latest_locator.get("s3_region")
        assert census.context.tiledb_config["py.init_buffer_bytes"] == f"{1 * 1024 ** 3}"

        # Verify that the user-overridden config settings are set correctly in the TileDB context object.
        assert census.context.tiledb_config["soma.init_buffer_bytes"] == soma_init_buffer_bytes
        assert census.context.timestamp_ms == timestamp_ms


def test_open_soma_uri_with_custom_s3_region() -> None:
    assert get_default_soma_context().tiledb_config["vfs.s3.region"] != "region-1", "test pre-condition"

    with patch("cellxgene_census._open.soma.open") as m:
        cellxgene_census.open_soma(
            uri="s3://bucket/cell-census/2022-11-01/soma/",
            tiledb_config={"vfs.s3.region": "region-1"},
        )

        m.assert_called_once_with(
            "s3://bucket/cell-census/2022-11-01/soma/",
            mode="r",
            soma_type=soma.Collection,
            context=ANY,
        )
        assert m.call_args[1]["context"].tiledb_config["vfs.s3.region"] == "region-1"


def test_open_soma_census_version_always_uses_mirror_s3_region(
    requests_mock: rm.Mocker,
) -> None:
    assert get_default_soma_context().tiledb_config["vfs.s3.region"] != "mirror-region-1", "test pre-condition"

    mock_mirrors = {
        "default": "test-mirror",
        "test-mirror": {
            "provider": "S3",
            "base_uri": "s3://mirror-bucket/",
            "region": "mirror-region-1",
        },
    }
    requests_mock.get(CELL_CENSUS_MIRRORS_DIRECTORY_URL, json=mock_mirrors)

    dir = {
        "latest": "2022-11-01",
        "2022-11-01": {
            "release_date": "2022-11-30",
            "soma": {
                "relative_uri": "/cell-census/2022-11-01/soma/",
            },
        },
    }
    requests_mock.get(CELL_CENSUS_RELEASE_DIRECTORY_URL, json=dir)

    # Verify that the mirror's S3 region is used, overriding the default
    with patch("cellxgene_census._open.soma.open") as m:
        cellxgene_census.open_soma(census_version="latest")

        m.assert_called_once_with(
            "s3://mirror-bucket/cell-census/2022-11-01/soma/",
            mode="r",
            soma_type=soma.Collection,
            context=ANY,
        )
        assert m.call_args[1]["context"].tiledb_config["vfs.s3.region"] == "mirror-region-1"

    # Verify that the mirror's S3 region is used, overriding even a user-provided region
    with patch("cellxgene_census._open.soma.open") as m:
        cellxgene_census.open_soma(census_version="latest", tiledb_config={"vfs.s3.region": "region-2"})

        m.assert_called_once_with(
            "s3://mirror-bucket/cell-census/2022-11-01/soma/",
            mode="r",
            soma_type=soma.Collection,
            context=ANY,
        )
        assert m.call_args[1]["context"].tiledb_config["vfs.s3.region"] == "mirror-region-1"


def test_open_soma_invalid_args() -> None:
    with pytest.raises(
        ValueError,
        match=re.escape("Must specify either a census version or an explicit URI."),
    ):
        cellxgene_census.open_soma(census_version=None)

    with pytest.raises(
        ValueError,
        match=re.escape("Only one of tiledb_config and context can be specified."),
    ):
        cellxgene_census.open_soma(tiledb_config={}, context=soma.SOMATileDBContext())


def test_open_soma_errors(requests_mock: rm.Mocker) -> None:
    requests_mock.get(CELL_CENSUS_RELEASE_DIRECTORY_URL, json={})
    requests_mock.real_http = True
    with pytest.raises(
        ValueError,
        match=re.escape(
            'The "does-not-exist" Census version is not valid. Use get_census_version_directory() to retrieve available versions.'
        ),
    ):
        cellxgene_census.open_soma(census_version="does-not-exist")


def test_open_soma_uses_correct_mirror(requests_mock: rm.Mocker) -> None:
    mock_mirrors = {
        "default": "test-mirror",
        "test-mirror": {
            "provider": "S3",
            "base_uri": "s3://mirror-bucket-1/",
            "region": "region-1",
        },
        "test-mirror-2": {
            "provider": "S3",
            "base_uri": "s3://mirror-bucket-2/",
            "region": "region-2",
        },
    }
    requests_mock.get(CELL_CENSUS_MIRRORS_DIRECTORY_URL, json=mock_mirrors)

    dir = {
        "stable": "2022-11-01",
        "2022-11-01": {
            "release_date": "2022-11-30",
            "release_build": "2022-11-01",
            "soma": {
                "uri": "s3://ignored-bucket/cell-census/2022-11-01/soma/",
                "relative_uri": "/cell-census/2022-11-01/soma/",
                "s3_region": "ignored",
            },
            "h5ads": {
                "uri": "s3://ignored-bucket/cell-census/2022-11-01/h5ads/",
                "relative_uri": "/cell-census/2022-11-01/soma/",
                "s3_region": "ignored",
            },
        },
    }

    requests_mock.get(CELL_CENSUS_RELEASE_DIRECTORY_URL, json=dir)

    # Verify that the default mirror is used if no mirror is specified
    with patch("cellxgene_census._open._open_soma") as m:
        cellxgene_census.open_soma()
        m.assert_called_once_with(
            {
                "uri": "s3://mirror-bucket-1/cell-census/2022-11-01/soma/",
                "region": "region-1",
                "provider": "S3",
            },
            None,
        )

    # Verify that the correct mirror is used if a mirror parameter is specified
    with patch("cellxgene_census._open._open_soma") as m:
        cellxgene_census.open_soma(mirror="test-mirror-2")
        m.assert_called_once_with(
            {
                "uri": "s3://mirror-bucket-2/cell-census/2022-11-01/soma/",
                "region": "region-2",
                "provider": "S3",
            },
            None,
        )

    # Verify that an error is raised if a non existing mirror is specified
    with patch("cellxgene_census._open._open_soma") as m:
        with pytest.raises(
            ValueError,
            match=re.escape("Mirror not found."),
        ):
            cellxgene_census.open_soma(mirror="bogus-mirror")


def test_open_soma_rejects_non_s3_mirror(requests_mock: rm.Mocker) -> None:
    mock_mirrors = {
        "default": "test-mirror",
        "test-mirror": {"provider": "GCS", "base_uri": "gcs://mirror-bucket-1/"},
    }
    requests_mock.real_http = True
    requests_mock.get(CELL_CENSUS_MIRRORS_DIRECTORY_URL, json=mock_mirrors)

    with pytest.raises(
        ValueError,
        match=re.escape(
            "Unsupported mirror provider: GCS. Try upgrading the cellxgene-census package to the latest version."
        ),
    ):
        cellxgene_census.open_soma()


def test_open_soma_works_if_no_relative_uri_specified(requests_mock: rm.Mocker) -> None:
    requests_mock.real_http = True
    """
    This test ensures that the Census works even if the relative_uri is not specified in the directory.
    This ensures backwards compatibility with the v1 route.
    """

    dir = {
        "stable": "2022-11-01",
        "2022-11-01": {
            "release_date": "2022-11-30",
            "release_build": "2022-11-01",
            "soma": {
                "uri": "s3://bucket-from-absolute-uri/cell-census/2022-11-01/soma/",
                "s3_region": "us-west-2",
            },
            "h5ads": {
                "uri": "s3://bucket-from-absolute-uri/cell-census/2022-11-01/h5ads/",
                "s3_region": "us-west-2",
            },
        },
    }

    requests_mock.get(CELL_CENSUS_RELEASE_DIRECTORY_URL, json=dir)
    with patch("cellxgene_census._open._open_soma") as m:
        cellxgene_census.open_soma()
        m.assert_called_once_with(
            {
                "uri": "s3://bucket-from-absolute-uri/cell-census/2022-11-01/soma/",
                "region": "us-west-2",
                "provider": "S3",
            },
            None,
        )


def test_open_soma_defaults_to_stable(requests_mock: rm.Mocker) -> None:
    requests_mock.real_http = True
    directory_with_stable = {
        "stable": "2022-10-01",
        "2022-10-01": {
            "release_date": "2022-10-30",
            "release_build": "2022-10-01",
            "soma": {
                "uri": "s3://cellxgene-census-public-us-west-2/cell-census/2022-10-01/soma/",
                "relative_uri": "/cell-census/2022-10-01/soma/",
                "s3_region": "us-west-2",
            },
            "h5ads": {
                "uri": "s3://cellxgene-census-public-us-west-2/cell-census/2022-10-01/h5ads/",
                "relative_uri": "/cell-census/2022-10-01/soma/",
                "s3_region": "us-west-2",
            },
        },
    }

    requests_mock.get(CELL_CENSUS_RELEASE_DIRECTORY_URL, json=directory_with_stable)
    with patch("cellxgene_census._open._open_soma") as m:
        cellxgene_census.open_soma()
        m.assert_called_once_with(
            {
                "uri": "s3://cellxgene-census-public-us-west-2/cell-census/2022-10-01/soma/",
                "region": "us-west-2",
                "provider": "S3",
            },
            None,
        )


@pytest.mark.live_corpus
def test_get_source_h5ad_uri() -> None:
    with cellxgene_census.open_soma(census_version="latest") as census:
        census_datasets = census["census_info"]["datasets"].read().concat().to_pandas()

    rng = np.random.default_rng()
    for idx in rng.choice(np.arange(len(census_datasets)), size=3, replace=False):
        a_dataset = census_datasets.iloc[idx]
        locator = cellxgene_census.get_source_h5ad_uri(a_dataset.dataset_id, census_version="latest")
        assert isinstance(locator, dict)
        assert "uri" in locator
        assert locator["uri"].endswith(a_dataset.dataset_h5ad_path)


def test_get_source_h5ad_uri_errors() -> None:
    with pytest.raises(KeyError):
        cellxgene_census.get_source_h5ad_uri(dataset_id="no/such/id")


@pytest.fixture
def small_dataset_id() -> str:
    with cellxgene_census.open_soma(census_version="latest") as census:
        census_datasets = census["census_info"]["datasets"].read().concat().to_pandas()

    small_dataset = census_datasets.nsmallest(1, "dataset_total_cell_count").iloc[0]
    assert isinstance(small_dataset.dataset_id, str)
    return small_dataset.dataset_id


@pytest.mark.live_corpus
def test_download_source_h5ad(tmp_path: pathlib.Path, small_dataset_id: str) -> None:
    adata_path = tmp_path / "adata.h5ad"
    cellxgene_census.download_source_h5ad(small_dataset_id, adata_path.as_posix(), census_version="latest")
    assert adata_path.exists() and adata_path.is_file()
    ad = anndata.read_h5ad(adata_path.as_posix())
    assert ad is not None


def test_download_source_h5ad_errors(tmp_path: pathlib.Path, small_dataset_id: str) -> None:
    existing_file = tmp_path / "existing_file.h5ad"
    existing_file.touch()
    assert existing_file.exists()

    with pytest.raises(ValueError):
        cellxgene_census.download_source_h5ad(small_dataset_id, existing_file.as_posix(), census_version="latest")

    with pytest.raises(ValueError):
        cellxgene_census.download_source_h5ad(small_dataset_id, "/tmp/dirname/", census_version="latest")


@pytest.mark.parametrize("progress_bar", [True, False])
def test_download_h5ad_progress_bar(  # type: ignore[no-untyped-def]
    capsys,
    tmp_path: pathlib.Path,
    small_dataset_id: str,
    progress_bar: bool,
) -> None:
    adata_path = tmp_path / "adata.h5ad"
    _ = capsys.readouterr()  # Clearing any previously captured output
    cellxgene_census.download_source_h5ad(
        small_dataset_id, adata_path.as_posix(), census_version="latest", progress_bar=progress_bar
    )
    captured = capsys.readouterr()
    if progress_bar:
        assert ("Downloading" in captured.err) and ("100%" in captured.err)
    else:
        assert captured.err == ""


@pytest.mark.live_corpus
def test_opening_census_without_anon_access_fails_with_bogus_creds() -> None:
    os.environ["AWS_ACCESS_KEY_ID"] = "fake_id"
    os.environ["AWS_SECRET_ACCESS_KEY"] = "fake_key"
    # Passing an empty context
    with pytest.raises(
        (tiledb.TileDBError, soma.DoesNotExistError),
        match=r"does not exist",
    ):
        cellxgene_census.open_soma(census_version="latest", context=soma.SOMATileDBContext())


@pytest.mark.live_corpus
def test_can_open_with_anonymous_access() -> None:
    """
    With anonymous access, `open_soma` must be able to access the census even with bogus credentials
    """
    os.environ["AWS_ACCESS_KEY_ID"] = "fake_id"
    os.environ["AWS_SECRET_ACCESS_KEY"] = "fake_key"
    with cellxgene_census.open_soma(census_version="latest") as census:
        assert census is not None
        assert isinstance(census, soma.Collection)


def test_get_default_soma_context_tiledb_config_overrides() -> None:
    context = get_default_soma_context(
        tiledb_config={
            "nondefault.config.option": "true",
            "vfs.s3.no_sign_request": "false",
        }
    )
    assert context.tiledb_config["nondefault.config.option"] == "true", "adds new option"
    assert context.tiledb_config["vfs.s3.no_sign_request"] == "false", "overrides existing default"
    assert context.tiledb_config["vfs.s3.region"] == "us-west-2", "keeps existing default"



# Section: cellxgene_census-src-cellxgene_census-_get_anndata

# Copyright (c) 2022, Chan Zuckerberg Initiative
#
# Licensed under the MIT License.

"""Get slice as AnnData.

Methods to retrieve slices of the census as AnnData objects.
"""

from collections.abc import Sequence
from typing import Literal
from warnings import warn

import anndata
import pandas as pd
import tiledbsoma as soma
from somacore.options import SparseDFCoord

from ._experiment import _get_experiment, _get_experiment_name
from ._release_directory import get_census_version_directory
from ._util import _extract_census_version, _uri_join

CENSUS_EMBEDDINGS_LOCATION_BASE_URI = "s3://cellxgene-contrib-public/contrib/cell-census/soma/"


def get_anndata(
    census: soma.Collection,
    organism: str,
    measurement_name: str = "RNA",
    X_name: str = "raw",
    X_layers: Sequence[str] | None = (),
    obsm_layers: Sequence[str] | None = (),
    obsp_layers: Sequence[str] | None = (),
    varm_layers: Sequence[str] | None = (),
    varp_layers: Sequence[str] | None = (),
    obs_value_filter: str | None = None,
    obs_coords: SparseDFCoord | None = None,
    var_value_filter: str | None = None,
    var_coords: SparseDFCoord | None = None,
    column_names: soma.AxisColumnNames | None = None,
    obs_embeddings: Sequence[str] | None = (),
    var_embeddings: Sequence[str] | None = (),
    obs_column_names: Sequence[str] | None = None,
    var_column_names: Sequence[str] | None = None,
) -> anndata.AnnData:
    """Convenience wrapper around :class:`tiledbsoma.Experiment` query, to build and execute a query,
    and return it as an :class:`anndata.AnnData` object.

    Args:
        census:
            The census object, usually returned by :func:`open_soma`.
        organism:
            The organism to query, usually one of ``"Homo sapiens`` or ``"Mus musculus"``.
        measurement_name:
            The measurement object to query. Defaults to ``"RNA"``.
        X_name:
            The ``X`` layer to query. Defaults to ``"raw"``.
        X_layers:
            Additional layers to add to :attr:`anndata.AnnData.layers`.
        obs_value_filter:
            Value filter for the ``obs`` metadata. Value is a filter query written in the
            SOMA ``value_filter`` syntax.
        obs_coords:
            Coordinates for the ``obs`` axis, which is indexed by the ``soma_joinid`` value.
            May be an ``int``, a list of ``int``, or a slice. The default, ``None``, selects all.
        var_value_filter:
            Value filter for the ``var`` metadata. Value is a filter query written in the
            SOMA ``value_filter`` syntax.
        var_coords:
            Coordinates for the ``var`` axis, which is indexed by the ``soma_joinid`` value.
            May be an ``int``, a list of ``int``, or a slice. The default, ``None``, selects all.
        obsm_layers:
            Additional obsm layers to read and return in the ``obsm`` slot.
        obsp_layers:
            Additional obsp layers to read and return in the ``obsp`` slot.
        varm_layers:
            Additional varm layers to read and return in the ``varm`` slot.
        varp_layers:
            Additional varp layers to read and return in the ``varp`` slot.
        obs_embeddings:
            Additional embeddings to be returned as part of the ``obsm`` slot.
            Use :func:`get_all_available_embeddings` to retrieve available embeddings
            for this Census version and organism.
        var_embeddings:
            Additional embeddings to be returned as part of the ``varm`` slot.
            Use :func:`get_all_available_embeddings` to retrieve available embeddings
            for this Census version and organism.
        obs_column_names:
            Columns to fetch for ``obs`` dataframe.
        var_column_names:
            Columns to fetch for ``var`` dataframe.

    Returns:
        An :class:`anndata.AnnData` object containing the census slice.

    Lifecycle:
        experimental

    Examples:
        >>> get_anndata(census, "Mus musculus", obs_value_filter="tissue_general in ['brain', 'lung']")

        >>> get_anndata(census, "Homo sapiens", obs_column_names=["tissue"])

        >>> get_anndata(census, "Homo sapiens", obs_coords=slice(0, 1000))
    """
    exp = _get_experiment(census, organism)
    obs_coords = (slice(None),) if obs_coords is None else (obs_coords,)
    var_coords = (slice(None),) if var_coords is None else (var_coords,)

    if obsm_layers and obs_embeddings and set(obsm_layers) & set(obs_embeddings):
        raise ValueError("Cannot request both `obsm_layers` and `obs_embeddings` for the same embedding name")

    if varm_layers and var_embeddings and set(varm_layers) & set(var_embeddings):
        raise ValueError("Cannot request both `varm_layers` and `var_embeddings` for the same embedding name")

    # Backwards compat for old column_names argument
    if column_names is not None:
        if obs_column_names is not None or var_column_names is not None:
            raise ValueError(
                "Both the deprecated 'column_names' argument and its replacements were used. Please use 'obs_column_names' and 'var_column_names' only."
            )
        else:
            warn(
                "The argument `column_names` is deprecated and will be removed in a future release. Please use `obs_column_names` and `var_column_names` instead.",
                FutureWarning,
                stacklevel=2,
            )
        if "obs" in column_names:
            obs_column_names = column_names["obs"]
        if "var" in column_names:
            var_column_names = column_names["var"]

    with exp.axis_query(
        measurement_name,
        obs_query=soma.AxisQuery(value_filter=obs_value_filter, coords=obs_coords),
        var_query=soma.AxisQuery(value_filter=var_value_filter, coords=var_coords),
    ) as query:
        adata = query.to_anndata(
            X_name=X_name,
            column_names={"obs": obs_column_names, "var": var_column_names},
            X_layers=X_layers,
            obsm_layers=obsm_layers,
            varm_layers=varm_layers,
            obsp_layers=obsp_layers,
            varp_layers=varp_layers,
        )

        # If obs_embeddings or var_embeddings are defined, inject them in the appropriate slot
        if obs_embeddings or var_embeddings:
            from .experimental._embedding import _get_embedding, get_embedding_metadata_by_name

            census_version = _extract_census_version(census)
            experiment_name = _get_experiment_name(organism)
            census_directory = get_census_version_directory()

            if obs_embeddings:
                obs_soma_joinids = query.obs_joinids()
                for emb in obs_embeddings:
                    emb_metadata = get_embedding_metadata_by_name(emb, experiment_name, census_version, "obs_embedding")
                    uri = _uri_join(CENSUS_EMBEDDINGS_LOCATION_BASE_URI, f"{census_version}/{emb_metadata['id']}")
                    embedding = _get_embedding(census, census_directory, census_version, uri, obs_soma_joinids)
                    adata.obsm[emb] = embedding

            if var_embeddings:
                var_soma_joinids = query.var_joinids()
                for emb in var_embeddings:
                    emb_metadata = get_embedding_metadata_by_name(emb, experiment_name, census_version, "var_embedding")
                    uri = _uri_join(CENSUS_EMBEDDINGS_LOCATION_BASE_URI, f"{census_version}/{emb_metadata['id']}")
                    embedding = _get_embedding(census, census_directory, census_version, uri, var_soma_joinids)
                    adata.varm[emb] = embedding

        return adata


def _get_axis_metadata(
    census: soma.Collection,
    axis: Literal["obs", "var"],
    organism: str,
    *,
    value_filter: str | None = None,
    coords: SparseDFCoord | None = slice(None),
    column_names: Sequence[str] | None = None,
) -> pd.DataFrame:
    exp = _get_experiment(census, organism)
    coords = (slice(None),) if coords is None else (coords,)
    if axis == "obs":
        df = exp.obs
    elif axis == "var":
        df = exp.ms["RNA"].var
    else:
        raise ValueError(f"axis should be either 'obs' or 'var', but '{axis}' was passed")
    result: pd.DataFrame = (
        df.read(coords=coords, column_names=column_names, value_filter=value_filter).concat().to_pandas()
    )
    return result


def get_obs(
    census: soma.Collection,
    organism: str,
    *,
    value_filter: str | None = None,
    coords: SparseDFCoord | None = slice(None),
    column_names: Sequence[str] | None = None,
) -> pd.DataFrame:
    """Get the observation metadata for a query on the census.

    Args:
        census:
            The census object, usually returned by :func:`open_soma`.
        organism:
            The organism to query, usually one of ``"Homo sapiens`` or ``"Mus musculus"``
        value_filter:
            Value filter for the ``obs`` metadata. Value is a filter query written in the
            SOMA ``value_filter`` syntax.
        coords:
            Coordinates for the ``obs`` axis, which is indexed by the ``soma_joinid`` value.
            May be an ``int``, a list of ``int``, or a slice. The default, ``None``, selects all.
        column_names:
            Columns to fetch.

    Returns:
        A :class:`pandas.DataFrame` object containing metadata for the queried slice.
    """
    return _get_axis_metadata(
        census, "obs", organism, value_filter=value_filter, coords=coords, column_names=column_names
    )


def get_var(
    census: soma.Collection,
    organism: str,
    *,
    value_filter: str | None = None,
    coords: SparseDFCoord | None = slice(None),
    column_names: Sequence[str] | None = None,
) -> pd.DataFrame:
    """Get the variable metadata for a query on the census.

    Args:
        census:
            The census object, usually returned by :func:`open_soma`.
        organism:
            The organism to query, usually one of ``"Homo sapiens`` or ``"Mus musculus"``
        value_filter:
            Value filter for the ``var`` metadata. Value is a filter query written in the
            SOMA ``value_filter`` syntax.
        coords:
            Coordinates for the ``var`` axis, which is indexed by the ``soma_joinid`` value.
            May be an ``int``, a list of ``int``, or a slice. The default, ``None``, selects all.
        column_names:
            Columns to fetch.

    Returns:
        A :class:`pandas.DataFrame` object containing metadata for the queried slice.
    """
    return _get_axis_metadata(
        census, "var", organism, value_filter=value_filter, coords=coords, column_names=column_names
    )



# Section: cellxgene_census-src-cellxgene_census-experimental-ml-huggingface-cell_dataset_builder

import uuid
from abc import ABC, abstractmethod
from collections.abc import Generator
from typing import Any

import scipy.sparse
from datasets import Dataset
from tiledbsoma import Experiment, ExperimentAxisQuery


class CellDatasetBuilder(ExperimentAxisQuery[Experiment], ABC):  # type: ignore
    """Abstract base class for methods to process CELLxGENE Census ExperimentAxisQuery
    results into a Hugging Face Dataset in which each item represents one cell.
    Subclasses implement the `cell_item()` method to process each row of an X layer
    into a Dataset item, and may also override `__init__()` and context `__enter__()`
    to perform any necessary preprocessing.

    The base class inherits ExperimentAxisQuery, so typical usage would be:

    ```
    import cellxgene_census
    import tiledbsoma
    from cellxgene_census.experimental.ml import GeneformerTokenizer

    with cellxgene_census.open_soma() as census:
        with SubclassOfCellDatasetBuilder(
            census["census_data"]["homo_sapiens"],
            obs_query=tilebsoma.AxisQuery(...),  # define some subset of Census cells
            ... # other ExperimentAxisQuery parameters e.g. var_query
        ) as builder:
            dataset = builder.build()
    ```
    """

    def __init__(
        self,
        experiment: Experiment,
        measurement_name: str = "RNA",
        layer_name: str = "raw",
        *,
        block_size: int | None = None,
        **kwargs: Any,
    ):
        """Initialize the CellDatasetBuilder to process the results of a Census
        ExperimentAxisQuery.

        - `experiment`: Census Experiment to be queried.
        - `measurement_name`: Measurement in the experiment, default "RNA".
        - `layer_name`: Name of the X layer to process, default "raw".
        - `block_size`: Number of cells to process in-memory at once. If unspecified,
           `tiledbsoma.SparseNDArrayRead.blockwise()` will select a default.
        - `kwargs`: passed through to `ExperimentAxisQuery()`, especially `obs_query`
           and `var_query`.
        """
        super().__init__(experiment, measurement_name, **kwargs)
        self.layer_name = layer_name
        self.block_size = block_size

    def build(self, from_generator_kwargs: dict[str, Any] | None = None) -> Dataset:
        """Build the dataset from query results.

        - `from_generator_kwargs`: kwargs passed through to `Dataset.from_generator()`
        """

        def gen() -> Generator[dict[str, Any], None, None]:
            for Xblock, (block_cell_joinids, _) in (
                self.X(self.layer_name).blockwise(axis=0, reindex_disable_on_axis=[1], size=self.block_size).scipy()
            ):
                assert isinstance(Xblock, scipy.sparse.csr_matrix)
                assert Xblock.shape[0] == len(block_cell_joinids)
                for i, cell_joinid in enumerate(block_cell_joinids):
                    yield self.cell_item(cell_joinid, Xblock.getrow(i))

        return Dataset.from_generator(_DatasetGeneratorPickleHack(gen), **(from_generator_kwargs or {}))

    @abstractmethod
    def cell_item(self, cell_joinid: int, Xrow: scipy.sparse.csr_matrix) -> dict[str, Any]:
        """Abstract method to process the X row for one cell into a Dataset item.

        - `cell_joinid`: The cell `soma_joinid`.
        - `Xrow`: The `X` row for this cell. This csr_matrix has a single row 0, equal
          to the `cell_joinid` row of the full `X` layer matrix.
        """
        ...


class _DatasetGeneratorPickleHack:
    """SEE: https://github.com/huggingface/datasets/issues/6194."""

    def __init__(self, generator: Any, generator_id: str | None = None) -> None:
        self.generator = generator
        self.generator_id = generator_id if generator_id is not None else str(uuid.uuid4())

    def __call__(self, *args: Any, **kwargs: Any) -> Any:
        return self.generator(*args, **kwargs)

    def __reduce__(self) -> Any:
        return (_DatasetGeneratorPickleHack_raise, (self.generator_id,))


def _DatasetGeneratorPickleHack_raise(*args: Any, **kwargs: Any) -> None:
    raise AssertionError("cannot actually unpickle _DatasetGeneratorPickleHack!")



# Section: cellxgene_census-tests-experimental-test_embeddings

from functools import partial

import pytest
import requests_mock as rm

import cellxgene_census
from cellxgene_census.experimental import (
    get_all_available_embeddings,
    get_all_census_versions_with_embedding,
    get_embedding_metadata_by_name,
)
from cellxgene_census.experimental._embedding import CELL_CENSUS_EMBEDDINGS_MANIFEST_URL


def test_get_embedding_metadata_by_name(requests_mock: rm.Mocker) -> None:
    mock_embeddings = {
        "embedding-id-1": {
            "id": "embedding-id-1",
            "embedding_name": "emb_1",
            "title": "Embedding 1",
            "description": "First embedding",
            "experiment_name": "homo_sapiens",
            "data_type": "obs_embedding",
            "census_version": "2023-12-15",
            "submission_date": "2023-11-15",
        },
        "embedding-id-2": {
            "id": "embedding-id-2",
            "embedding_name": "emb_1",
            "title": "Embedding 2",
            "description": "Second embedding",
            "experiment_name": "homo_sapiens",
            "data_type": "obs_embedding",
            "census_version": "2023-12-15",
            "submission_date": "2023-12-31",
        },
        "embedding-id-3": {
            "id": "embedding-id-3",
            "embedding_name": "emb_3",
            "title": "Embedding 3",
            "description": "Third embedding",
            "experiment_name": "homo_sapiens",
            "data_type": "obs_embedding",
            "census_version": "2023-12-15",
            "submission_date": "2023-11-15",
        },
    }
    requests_mock.real_http = True
    requests_mock.get(CELL_CENSUS_EMBEDDINGS_MANIFEST_URL, json=mock_embeddings)

    embedding = get_embedding_metadata_by_name(
        "emb_1", organism="homo_sapiens", census_version="2023-12-15", embedding_type="obs_embedding"
    )
    assert embedding is not None
    assert embedding["id"] == "embedding-id-2"  # most recent version
    assert embedding == mock_embeddings["embedding-id-2"]

    embedding = get_embedding_metadata_by_name(
        "emb_3", organism="homo_sapiens", census_version="2023-12-15", embedding_type="obs_embedding"
    )
    assert embedding is not None
    assert embedding["id"] == "embedding-id-3"
    assert embedding == mock_embeddings["embedding-id-3"]

    with pytest.raises(ValueError):
        get_embedding_metadata_by_name(
            "emb_2", organism="homo_sapiens", census_version="2023-12-15", embedding_type="obs_embedding"
        )
        get_embedding_metadata_by_name(
            "emb_1", organism="mus_musculus", census_version="2023-12-15", embedding_type="obs_embedding"
        )
        get_embedding_metadata_by_name(
            "emb_1", organism="homo_sapiens", census_version="2023-10-15", embedding_type="obs_embedding"
        )
        get_embedding_metadata_by_name(
            "emb_1", organism="mus_musculus", census_version="2023-12-15", embedding_type="var_embedding"
        )


def test_get_embedding_by_name_w_version_aliases() -> None:
    """https://github.com/chanzuckerberg/cellxgene-census/issues/1202"""
    # Only testing "stable" as "latest" doesn't have embeddings
    version = "stable"
    resolved_version = cellxgene_census.get_census_version_description(version)["release_build"]

    metadata = get_all_available_embeddings(version)[0]

    _get_metadata = partial(
        get_embedding_metadata_by_name,
        embedding_name=metadata["embedding_name"],
        organism=metadata["experiment_name"],
        embedding_type=metadata["data_type"],
    )

    w_alias = _get_metadata(census_version=version)
    w_resolved = _get_metadata(census_version=resolved_version)

    assert w_resolved == w_alias
    assert metadata == w_alias


def test_get_all_available_embeddings(requests_mock: rm.Mocker) -> None:
    mock_embeddings = {
        "embedding-id-1": {
            "id": "embedding-id-1",
            "embedding_name": "emb_1",
            "title": "Embedding 1",
            "description": "First embedding",
            "experiment_name": "homo_sapiens",
            "measurement_name": "RNA",
            "n_embeddings": 1000,
            "n_features": 200,
            "data_type": "obs_embedding",
            "census_version": "2023-12-15",
        },
        "embedding-id-2": {
            "id": "embedding-id-2",
            "embedding_name": "emb_2",
            "title": "Embedding 2",
            "description": "Second embedding",
            "experiment_name": "homo_sapiens",
            "measurement_name": "RNA",
            "n_embeddings": 1000,
            "n_features": 200,
            "data_type": "obs_embedding",
            "census_version": "2023-12-15",
        },
    }
    requests_mock.real_http = True
    requests_mock.get(CELL_CENSUS_EMBEDDINGS_MANIFEST_URL, json=mock_embeddings)

    embeddings = get_all_available_embeddings("2023-12-15")
    assert embeddings is not None
    assert len(embeddings) == 2

    # Query for a version of the census that doesn't have embeddings
    embeddings = get_all_available_embeddings("2023-05-15")
    assert len(embeddings) == 0


def test_get_all_census_versions_with_embedding(requests_mock: rm.Mocker) -> None:
    mock_embeddings = {
        "embedding-id-1": {
            "id": "embedding-id-1",
            "embedding_name": "emb_1",
            "title": "Embedding 1",
            "description": "First embedding",
            "experiment_name": "homo_sapiens",
            "data_type": "obs_embedding",
            "census_version": "2023-12-15",
        },
        "embedding-id-2": {
            "id": "embedding-id-2",
            "embedding_name": "emb_1",
            "title": "Embedding 2",
            "description": "Second embedding",
            "experiment_name": "homo_sapiens",
            "data_type": "obs_embedding",
            "census_version": "2023-12-15",
        },
        "embedding-id-3": {
            "id": "embedding-id-3",
            "embedding_name": "emb_1",
            "title": "Embedding 3",
            "description": "Third embedding",
            "experiment_name": "mus_musculus",
            "data_type": "obs_embedding",
            "census_version": "2023-12-15",
        },
        "embedding-id-4": {
            "id": "embedding-id-4",
            "embedding_name": "emb_1",
            "title": "Embedding 4",
            "description": "Fourth embedding",
            "experiment_name": "mus_musculus",
            "data_type": "obs_embedding",
            "census_version": "2024-01-01",
        },
        "embedding-id-5": {
            "id": "embedding-id-5",
            "embedding_name": "emb_2",
            "title": "Embedding 5",
            "description": "Fifth embedding",
            "experiment_name": "mus_musculus",
            "data_type": "var_embedding",
            "census_version": "2023-12-15",
        },
    }
    requests_mock.real_http = True
    requests_mock.get(CELL_CENSUS_EMBEDDINGS_MANIFEST_URL, json=mock_embeddings)

    versions = get_all_census_versions_with_embedding("emb_1", organism="homo_sapiens", embedding_type="obs_embedding")
    assert versions == ["2023-12-15"]

    versions = get_all_census_versions_with_embedding("emb_1", organism="mus_musculus", embedding_type="obs_embedding")
    assert versions == ["2023-12-15", "2024-01-01"]

    versions = get_all_census_versions_with_embedding("emb_1", organism="mus_musculus", embedding_type="var_embedding")
    assert versions == []

    versions = get_all_census_versions_with_embedding("emb_2", organism="mus_musculus", embedding_type="var_embedding")
    assert versions == ["2023-12-15"]


@pytest.mark.parametrize("version", ["stable", "latest"])
def test_get_all_available_embeddings_w_version_aliases(version: str) -> None:
    """https://github.com/chanzuckerberg/cellxgene-census/issues/1202"""
    resolved_version = cellxgene_census.get_census_version_description(version)["release_build"]

    assert get_all_available_embeddings(version) == get_all_available_embeddings(resolved_version)


def test_get_all_available_embeddings_non_existing_version() -> None:
    false_version = "not a real version"

    with pytest.raises(ValueError, match=f"Unable to locate Census version: {false_version}"):
        get_all_available_embeddings(false_version)



# Section: cellxgene_census-src-cellxgene_census-experimental-_embedding_search

"""Nearest-neighbor search based on vector index of Census embeddings."""

from collections.abc import Sequence
from contextlib import ExitStack
from typing import Any, NamedTuple, cast

import anndata as ad
import numpy as np
import numpy.typing as npt
import pandas as pd
import tiledbsoma as soma
from scipy import sparse

from .._experiment import _get_experiment_name
from .._open import DEFAULT_TILEDB_CONFIGURATION, open_soma
from .._release_directory import CensusMirror, _get_census_mirrors
from .._util import _uri_join
from ._embedding import get_embedding_metadata_by_name


class NeighborObs(NamedTuple):
    """Results of nearest-neighbor search for Census obs embeddings."""

    distances: npt.NDArray[np.float32]
    """
    Distances to the nearest neighbors for each query obs embedding (q by k, where q is the number
    of query embeddings and k is the desired number of neighbors). The distance metric is
    implementation-dependent.
    """

    neighbor_ids: npt.NDArray[np.int64]
    """
    obs soma_joinid's of the nearest neighbors for each query embedding (q by k).
    """


def find_nearest_obs(
    embedding_name: str,
    organism: str,
    census_version: str,
    query: ad.AnnData,
    *,
    k: int = 10,
    nprobe: int = 100,
    memory_GiB: int = 4,
    mirror: str | None = None,
    embedding_metadata: dict[str, Any] | None = None,
    **kwargs: dict[str, Any],
) -> NeighborObs:
    """Search Census for similar obs (cells) based on nearest neighbors in embedding space.

    Args:
        embedding_name, organism, census_version:
            Identify the embedding to search, as in :func:`get_embedding_metadata_by_name`.
        query:
            AnnData object with an obsm layer embedding the query cells. The obsm layer name
            matches ``embedding_metadata["embedding_name"]`` (e.g. scvi, geneformer). The layer
            shape matches the number of query cells and the number of features in the embedding.
        k:
            Number of nearest neighbors to return for each query obs.
        nprobe:
            Sensitivity parameter; defaults to 100 (roughly N^0.25 where N is the number of Census
            cells) for a thorough search. Decrease for faster but less accurate search.
        memory_GiB:
            Memory budget for the search index, in gibibytes; defaults to 4 GiB.
        mirror:
            Name of the Census mirror to use for the search.
        embedding_metadata:
            The result of `get_embedding_metadata_by_name(embedding_name, organism, census_version)`.
            Supplying this saves a network request for repeated searches.
    """
    import tiledb.vector_search as vs

    if embedding_metadata is None:
        embedding_metadata = get_embedding_metadata_by_name(embedding_name, organism, census_version)
    assert embedding_metadata["embedding_name"] == embedding_name
    n_features = embedding_metadata["n_features"]

    # validate query (expected obsm layer exists with the expected dimensionality)
    if embedding_name not in query.obsm:
        raise ValueError(f"Query does not have the expected layer {embedding_name}")
    if query.obsm[embedding_name].shape[1] != n_features:
        raise ValueError(
            f"Query embedding {embedding_name} has {query.obsm[embedding_name].shape[1]} features, expected {n_features}"
        )

    # formulate index URI and run query
    resolved_index = _resolve_embedding_index(embedding_metadata, mirror=mirror)
    if not resolved_index:
        raise ValueError("No suitable embedding index found for " + embedding_name)
    index_uri, index_region = resolved_index
    config = {k: str(v) for k, v in DEFAULT_TILEDB_CONFIGURATION.items()}
    config["vfs.s3.region"] = index_region
    memory_vectors = memory_GiB * (2**30) // (4 * n_features)  # number of float32 vectors
    index = vs.ivf_flat_index.IVFFlatIndex(uri=index_uri, config=config, memory_budget=memory_vectors)
    distances, neighbor_ids = index.query(query.obsm[embedding_name], k=k, nprobe=nprobe, **kwargs)

    return NeighborObs(distances=distances, neighbor_ids=neighbor_ids)


def _resolve_embedding_index(
    embedding_metadata: dict[str, Any],
    mirror: str | None = None,
) -> tuple[str, str] | None:
    index_metadata = embedding_metadata.get("indexes", None)
    if not index_metadata:
        return None
    # TODO (future): support multiple index [types]
    assert index_metadata[0]["type"] == "IVFFlat", "Only IVFFlat index is supported (update cellxgene_census)"
    mirrors = _get_census_mirrors()
    mirror = mirror or cast(str, mirrors["default"])
    mirror_info = cast(CensusMirror, mirrors[mirror])
    uri = _uri_join(mirror_info["embeddings_base_uri"], index_metadata[0]["relative_uri"])
    return uri, cast(str, mirror_info["region"])


def predict_obs_metadata(
    organism: str,
    census_version: str,
    neighbors: NeighborObs,
    column_names: Sequence[str],
    experiment: soma.Experiment | None = None,
) -> pd.DataFrame:
    """Predict obs metadata attributes for the query cells based on the embedding nearest neighbors.

    Args:
        organism, census_version:
            Embedding information as supplied to :func:`find_nearest_obs`.
        neighbors:
            Results of a :func:`find_nearest_obs` search.
        column_names:
            Desired obs metadata column names. The current implementation is suitable for
            categorical attributes (e.g. cell_type, tissue_general).
        experiment:
            Open handle for the relevant SOMAExperiment, if available (otherwise, will be opened
            internally). e.g. ``census["census_data"]["homo_sapiens"]`` with the relevant Census
            version.

    Returns:
        Pandas DataFrame with the desired column predictions. Additionally, for each predicted
        column ``col``, an additional column ``col_confidence`` with a confidence score between 0
        and 1.
    """
    with ExitStack() as cleanup:
        if experiment is None:
            # open Census transiently
            census = cleanup.enter_context(open_soma(census_version=census_version))
            experiment = census["census_data"][_get_experiment_name(organism)]

        # fetch the desired obs metadata for all of the found neighbors
        neighbor_obs = (
            experiment.obs.read(
                coords=(neighbors.neighbor_ids.flatten(),), column_names=(["soma_joinid"] + list(column_names))
            )
            .concat()
            .to_pandas()
        ).set_index("soma_joinid")

        # step through query cells to generate prediction for each column as the plurality value
        # found among its neighbors, with a confidence score based on the simple fraction (for now)
        # TODO: something more intelligent for numeric columns! also use distances, etc.
        max_joinid = neighbor_obs.index.max()
        out: dict[str, pd.Series[Any]] = {}
        n_queries, n_neighbors = neighbors.neighbor_ids.shape
        indices = np.broadcast_to(np.arange(n_queries), (n_neighbors, n_queries)).T
        g = sparse.csr_matrix(
            (
                np.broadcast_to(1, n_queries * n_neighbors),
                (
                    indices.flatten(),
                    neighbors.neighbor_ids.astype(np.int64).flatten(),
                ),
            ),
            shape=(n_queries, max_joinid + 1),
        )
        for col in column_names:
            col_categorical = neighbor_obs[col].astype("category")
            joinid2category = sparse.coo_matrix(
                (np.broadcast_to(1, len(neighbor_obs)), (neighbor_obs.index, col_categorical.cat.codes)),
                shape=(max_joinid + 1, len(col_categorical.cat.categories)),
            )
            counts = g @ joinid2category
            rel_counts = counts / counts.sum(axis=1)
            out[col] = col_categorical.cat.categories[rel_counts.argmax(axis=1).A.flatten()].astype(object)
            out[f"{col}_confidence"] = rel_counts.max(axis=1).toarray().flatten()

    return pd.DataFrame.from_dict(out)



# Section: notebooks-api_demo-census_and_cell_guide_example

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




# Section: cellxgene_census-src-cellxgene_census-__init__

"""An API to facilitate use of the CZI Science CELLxGENE Census. The Census is a versioned container of single-cell data hosted at `CELLxGENE Discover`_.

The API is built on the `tiledbsoma` SOMA API, and provides a number of helper functions including:

    * Open a named version of the Census, for use with the SOMA API
    * Get a list of available Census versions, and for each version, a description
    * Get a slice of the Census as an AnnData, for use with ScanPy
    * Get the URI for, or directly download, underlying data in H5AD format

For more information on the API, visit the `cellxgene_census repo`_. For more information on SOMA, see the `tiledbsoma repo`_.

.. _CELLxGENE Discover:
    https://cellxgene.cziscience.com/

.. _cellxgene_census repo:
    https://github.com/chanzuckerberg/cellxgene-census/

.. _tiledbsoma repo:
    https://github.com/single-cell-data/TileDB-SOMA
"""

from importlib import metadata

from ._get_anndata import get_anndata, get_obs, get_var
from ._open import (
    download_source_h5ad,
    get_default_soma_context,
    get_source_h5ad_uri,
    open_soma,
)
from ._presence_matrix import get_presence_matrix
from ._release_directory import (
    get_census_mirror_directory,
    get_census_version_description,
    get_census_version_directory,
)

try:
    __version__ = metadata.version("cellxgene_census")
except metadata.PackageNotFoundError:
    # package is not installed
    __version__ = "0.0.0-unknown"

__all__ = [
    "download_source_h5ad",
    "get_anndata",
    "get_obs",
    "get_var",
    "get_census_version_description",
    "get_census_version_directory",
    "get_census_mirror_directory",
    "get_default_soma_context",
    "get_presence_matrix",
    "get_source_h5ad_uri",
    "open_soma",
]



# Section: cellxgene_census-src-cellxgene_census-experimental-ml-pytorch

import gc
import itertools
import logging
import os
import typing
from collections.abc import Iterator, Sequence
from contextlib import contextmanager
from datetime import timedelta
from math import ceil
from time import time
from typing import Any, TypeAlias

import numpy as np
import numpy.typing as npt
import pandas as pd
import psutil
import tiledbsoma as soma
import torch
import torchdata.datapipes.iter as pipes
from attr import define
from numpy.random import Generator
from pyarrow import Table
from scipy import sparse
from torch import Tensor
from torch import distributed as dist
from torch.utils.data import DataLoader
from torch.utils.data.dataset import Dataset

from ... import get_default_soma_context
from ..util._eager_iter import _EagerIterator
from .encoders import Encoder, LabelEncoder

pytorch_logger = logging.getLogger("cellxgene_census.experimental.pytorch")

# TODO: Rename to reflect the correct order of the Tensors within the tuple: (X, obs)
ObsAndXDatum = tuple[Tensor, Tensor]
"""Return type of ``ExperimentDataPipe`` that pairs a Tensor of ``obs`` row(s) with a Tensor of ``X`` matrix row(s).
The Tensors are rank 1 if ``batch_size`` is 1, otherwise the Tensors are rank 2."""


# "Chunk" of X data, returned by each `Method` above
ChunkX: TypeAlias = npt.NDArray[Any] | sparse.csr_matrix


@define
class _SOMAChunk:
    """Return type of ``_ObsAndXSOMAIterator`` that pairs a chunk of ``obs`` rows with the respective rows from the ``X``
    matrix.

    Lifecycle:
        experimental
    """

    obs: pd.DataFrame
    X: ChunkX
    stats: "Stats"

    def __len__(self) -> int:
        return len(self.obs)


Encoders = dict[str, Encoder]
"""A dictionary of ``Encoder``s keyed by the ``obs`` column name."""


@define
class Stats:
    """Statistics about the data retrieved by ``ExperimentDataPipe`` via SOMA API. This is useful for assessing the read
    throughput of SOMA data.

    Lifecycle:
        experimental
    """

    n_obs: int = 0
    """The total number of obs rows retrieved"""

    nnz: int = 0
    """The total number of values retrieved"""

    elapsed: float = 0
    """The total elapsed time in seconds for retrieving all batches"""

    n_soma_chunks: int = 0
    """The number of chunks retrieved"""

    def __str__(self) -> str:
        return f"{self.n_soma_chunks=}, {self.n_obs=}, {self.nnz=}, " f"elapsed={timedelta(seconds=self.elapsed)}"

    def __add__(self, other: "Stats") -> "Stats":
        self.n_obs += other.n_obs
        self.nnz += other.nnz
        self.elapsed += other.elapsed
        self.n_soma_chunks += other.n_soma_chunks
        return self


@contextmanager
def _open_experiment(
    uri: str,
    aws_region: str | None = None,
) -> soma.Experiment:
    """Internal method for opening a SOMA ``Experiment`` as a context manager."""
    context = get_default_soma_context().replace(tiledb_config={"vfs.s3.region": aws_region} if aws_region else {})

    with soma.Experiment.open(uri, context=context) as exp:
        yield exp


def _tables_to_np(
    tables: Iterator[tuple[Table, Any]], shape: tuple[int, int]
) -> typing.Generator[tuple[npt.NDArray[Any], Any, int], None, None]:
    for tbl, indices in tables:
        row_indices, col_indices, data = (x.to_numpy() for x in tbl.columns)
        nnz = len(data)
        dense_matrix = np.zeros(shape, dtype=data.dtype)
        dense_matrix[row_indices, col_indices] = data
        yield dense_matrix, indices, nnz


class _ObsAndXSOMAIterator(Iterator[_SOMAChunk]):
    """Iterates the SOMA chunks of corresponding ``obs`` and ``X`` data. This is an internal class,
    not intended for public use.
    """

    X: soma.SparseNDArray
    """A handle to the full X data of the SOMA ``Experiment``"""

    obs_joinids_chunks_iter: Iterator[npt.NDArray[np.int64]]

    var_joinids: npt.NDArray[np.int64]
    """The ``var`` joinids to be retrieved from the SOMA ``Experiment``"""

    def __init__(
        self,
        obs: soma.DataFrame,
        X: soma.SparseNDArray,
        obs_column_names: Sequence[str],
        obs_joinids_chunked: list[npt.NDArray[np.int64]],
        var_joinids: npt.NDArray[np.int64],
        shuffle_chunk_count: int | None = None,
        shuffle_rng: Generator | None = None,
        return_sparse_X: bool = False,
    ):
        self.obs = obs
        self.X = X
        self.obs_column_names = obs_column_names
        if shuffle_chunk_count is not None:
            assert shuffle_rng is not None

            # At the start of this step, `obs_joinids_chunked` is a list of one dimensional
            # numpy arrays. Each numpy array corresponds to a chunk of contiguous rows in `obs`.
            # Critically, `obs_joinids_chunked` is randomly ordered where each chunk is
            # from a random section of `obs`.
            # We then take `shuffle_chunk_count` of these in order, concatenate them into
            # a larger numpy array and shuffle this larger numpy array.
            # The result is again a list of numpy arrays.
            self.obs_joinids_chunks_iter = (
                shuffle_rng.permutation(np.concatenate(grouped_chunks))
                for grouped_chunks in list_split(obs_joinids_chunked, shuffle_chunk_count)
            )
        else:
            self.obs_joinids_chunks_iter = iter(obs_joinids_chunked)
        self.var_joinids = var_joinids
        self.shuffle_chunk_count = shuffle_chunk_count
        self.return_sparse_X = return_sparse_X

    def __next__(self) -> _SOMAChunk:
        pytorch_logger.debug("Retrieving next SOMA chunk...")
        start_time = time()

        # If no more chunks to iterate through, raise StopIteration, as all iterators do when at end
        obs_joinids_chunk = next(self.obs_joinids_chunks_iter)

        if "soma_joinid" not in self.obs_column_names:
            cols = ["soma_joinid", *self.obs_column_names]
        else:
            cols = list(self.obs_column_names)

        obs_batch = (
            self.obs.read(
                coords=(obs_joinids_chunk,),
                column_names=cols,
            )
            .concat()
            .to_pandas()
            .set_index("soma_joinid")
        )
        assert obs_batch.shape[0] == obs_joinids_chunk.shape[0]

        # handle case of empty result (first batch has 0 rows)
        if len(obs_batch) == 0:
            raise StopIteration

        # reorder obs rows to match obs_joinids_chunk ordering, which may be shuffled
        obs_batch = obs_batch.reindex(obs_joinids_chunk, copy=False)

        # note: the `blockwise` call is employed for its ability to reindex the axes of the sparse matrix,
        # but the blockwise iteration feature is not used (block_size is set to retrieve the chunk as a single block)
        blockwise_iter = self.X.read(coords=(obs_joinids_chunk, self.var_joinids)).blockwise(
            axis=0, size=len(obs_joinids_chunk), eager=False
        )

        X_batch: ChunkX
        if not self.return_sparse_X:
            res = next(_tables_to_np(blockwise_iter.tables(), shape=(obs_batch.shape[0], len(self.var_joinids))))
            X_batch, nnz = res[0], res[2]
        else:
            X_batch = next(blockwise_iter.scipy(compress=True))[0]
            nnz = X_batch.nnz

        assert obs_batch.shape[0] == X_batch.shape[0]

        end_time = time()
        stats = Stats()
        stats.n_obs += X_batch.shape[0]
        stats.nnz += nnz
        stats.elapsed += end_time - start_time
        stats.n_soma_chunks += 1

        pytorch_logger.debug(f"Retrieved SOMA chunk: {stats}")
        return _SOMAChunk(obs=obs_batch, X=X_batch, stats=stats)


def list_split(arr_list: list[Any], sublist_len: int) -> list[list[Any]]:
    """Splits a python list into a list of sublists where each sublist is of size `sublist_len`.
    TODO: Replace with `itertools.batched` when Python 3.12 becomes the minimum supported version.
    """
    i = 0
    result = []
    while i < len(arr_list):
        if (i + sublist_len) >= len(arr_list):
            result.append(arr_list[i:])
        else:
            result.append(arr_list[i : i + sublist_len])

        i += sublist_len

    return result


def run_gc() -> tuple[tuple[Any, Any, Any], tuple[Any, Any, Any], float]:  # noqa: D103
    proc = psutil.Process(os.getpid())

    pre_gc = proc.memory_full_info(), psutil.virtual_memory(), psutil.swap_memory()
    start = time()
    gc.collect()
    gc_elapsed = time() - start
    post_gc = proc.memory_full_info(), psutil.virtual_memory(), psutil.swap_memory()

    pytorch_logger.debug(f"gc:  pre={pre_gc}")
    pytorch_logger.debug(f"gc: post={post_gc}")

    return pre_gc, post_gc, gc_elapsed


class _ObsAndXIterator(Iterator[ObsAndXDatum]):
    """Iterates through a set of ``obs`` and corresponding ``X`` rows, where the rows to be returned are specified by
    the ``obs_tables_iter`` argument. For the specified ``obs` rows, the corresponding ``X`` data is loaded and
    joined together. It is returned from this iterator as 2-tuples of ``X`` and obs Tensors.

    Internally manages the retrieval of data in SOMA-sized chunks, fetching the next chunk of SOMA data as needed.
    Supports fetching the data in an eager manner, where the next SOMA chunk is fetched while the current chunk is
    being read. This is an internal class, not intended for public use.
    """

    soma_chunk_iter: Iterator[_SOMAChunk]
    """The iterator for SOMA chunks of paired obs and X data"""

    soma_chunk: _SOMAChunk | None
    """The current SOMA chunk of obs and X data"""

    i: int = -1
    """Index into current obs ``SOMA`` chunk"""

    def __init__(
        self,
        obs: soma.DataFrame,
        X: soma.SparseNDArray,
        obs_column_names: Sequence[str],
        obs_joinids_chunked: list[npt.NDArray[np.int64]],
        var_joinids: npt.NDArray[np.int64],
        batch_size: int,
        encoders: list[Encoder],
        stats: Stats,
        return_sparse_X: bool,
        use_eager_fetch: bool,
        shuffle_chunk_count: int | None = None,
        shuffle_rng: Generator | None = None,
    ) -> None:
        self.soma_chunk_iter = _ObsAndXSOMAIterator(
            obs,
            X,
            obs_column_names,
            obs_joinids_chunked,
            var_joinids,
            shuffle_chunk_count,
            shuffle_rng,
            return_sparse_X=return_sparse_X,
        )
        if use_eager_fetch:
            self.soma_chunk_iter = _EagerIterator(self.soma_chunk_iter)
        self.soma_chunk = None
        self.var_joinids = var_joinids
        self.batch_size = batch_size
        self.return_sparse_X = return_sparse_X
        self.encoders = encoders
        self.stats = stats
        self.gc_elapsed = 0.0
        self.max_process_mem_usage_bytes = 0
        self.X_dtype = X.schema[2].type.to_pandas_dtype()

    def __next__(self) -> ObsAndXDatum:
        """Read the next torch batch, possibly across multiple soma chunks."""
        obss: list[pd.DataFrame] = []
        Xs: list[ChunkX] = []
        n_obs = 0

        while n_obs < self.batch_size:
            try:
                obs_partial, X_partial = self._read_partial_torch_batch(self.batch_size - n_obs)
                n_obs += len(obs_partial)
                obss.append(obs_partial)
                Xs.append(X_partial)
            except StopIteration:
                break

        if len(Xs) == 0:  # If we ran out of data
            raise StopIteration
        else:
            if self.return_sparse_X:
                X = sparse.vstack(Xs)
            else:
                X = np.concatenate(Xs, axis=0)
            obs = pd.concat(obss, axis=0)

        obs_encoded = pd.DataFrame()

        # Add the soma_joinid to the original obs, in case that is requested by the encoders.
        obs["soma_joinid"] = obs.index

        for enc in self.encoders:
            obs_encoded[enc.name] = enc.transform(obs)

        # `to_numpy()` avoids copying the numpy array data
        obs_tensor = torch.from_numpy(obs_encoded.to_numpy())

        if not self.return_sparse_X:
            X_tensor = torch.from_numpy(X)
        else:
            coo = X.tocoo()

            X_tensor = torch.sparse_coo_tensor(
                # Note: The `np.array` seems unnecessary, but PyTorch warns bare array is "extremely slow"
                indices=torch.from_numpy(np.array([coo.row, coo.col])),
                values=coo.data,
                size=coo.shape,
            )

        if self.batch_size == 1:
            X_tensor = X_tensor[0]
            obs_tensor = obs_tensor[0]

        return X_tensor, obs_tensor

    def _read_partial_torch_batch(self, batch_size: int) -> tuple[pd.DataFrame, ChunkX]:
        """Reads a torch-size batch of data from the current SOMA chunk, returning a torch-size batch whose size may
        contain fewer rows than the requested ``batch_size``. This can happen when the remaining rows in the current
        SOMA chunk are fewer than the requested ``batch_size``.
        """
        if self.soma_chunk is None or not (0 <= self.i < len(self.soma_chunk)):
            # GC memory from previous soma_chunk
            self.soma_chunk = None
            pre_gc, _, gc_elapsed = run_gc()
            self.max_process_mem_usage_bytes = max(self.max_process_mem_usage_bytes, pre_gc[0].uss)

            self.soma_chunk: _SOMAChunk = next(self.soma_chunk_iter)
            self.stats += self.soma_chunk.stats
            self.gc_elapsed += gc_elapsed
            self.i = 0

            pytorch_logger.debug(
                f"Retrieved SOMA chunk totals: {self.stats}, gc_elapsed={timedelta(seconds=self.gc_elapsed)}"
            )

        obs_batch = self.soma_chunk.obs
        X_chunk = self.soma_chunk.X

        safe_batch_size = min(batch_size, len(obs_batch) - self.i)
        slice_ = slice(self.i, self.i + safe_batch_size)
        assert slice_.stop <= obs_batch.shape[0]

        obs_rows = obs_batch.iloc[slice_]
        assert obs_rows.index.is_unique
        assert safe_batch_size == obs_rows.shape[0]

        X_batch = X_chunk[slice_]

        assert obs_rows.shape[0] == X_batch.shape[0]

        self.i += safe_batch_size

        return obs_rows, X_batch


class ExperimentDataPipe(pipes.IterDataPipe[Dataset[ObsAndXDatum]]):  # type: ignore
    r"""An :class:`torchdata.datapipes.iter.IterDataPipe` that reads ``obs`` and ``X`` data from a
    :class:`tiledbsoma.Experiment`, based upon the specified queries along the ``obs`` and ``var`` axes. Provides an
    iterator over these data when the object is passed to Python's built-in ``iter`` function.

    >>> for batch in iter(ExperimentDataPipe(...)):
            X_batch, y_batch = batch

    The ``batch_size`` parameter controls the number of rows of ``obs`` and ``X`` data that are returned in each
    iteration. If the ``batch_size`` is 1, then each Tensor will have rank 1:

    >>> (tensor([0., 0., 0., 0., 0., 1., 0., 0., 0.]),  # X data
         tensor([2415,    0,    0], dtype=torch.int64)) # obs data, encoded

    For larger ``batch_size`` values, the returned Tensors will have rank 2:

    >>> DataLoader(..., batch_size=3, ...):
        (tensor([[0., 0., 0., 0., 0., 1., 0., 0., 0.],     # X batch
                 [0., 0., 0., 0., 0., 0., 0., 0., 0.],
                 [0., 0., 0., 0., 0., 0., 0., 0., 0.]]),
         tensor([[2415,    0,    0],                       # obs batch
                 [2416,    0,    4],
                 [2417,    0,    3]], dtype=torch.int64))

    The ``return_sparse_X`` parameter controls whether the ``X`` data is returned as a dense or sparse
    :class:`torch.Tensor`. If the model supports use of sparse :class:`torch.Tensor`\ s, this will reduce memory usage.

    The ``obs_column_names`` parameter determines the data columns that are returned in the ``obs`` Tensor. String-typed
    columns are encoded as integer values. If needed, these values can be decoded by obtaining the encoder for a given
    ``obs`` column name and calling its ``inverse_transform`` method:

    >>> exp_data_pipe.obs_encoders["<obs_attr_name>"].inverse_transform(encoded_values)

    Lifecycle:
        experimental
    """

    _initialized: bool

    _obs_joinids: npt.NDArray[np.int64] | None

    _var_joinids: npt.NDArray[np.int64] | None

    _encoders: list[Encoder]

    _stats: Stats

    _shuffle_rng: Generator | None

    # TODO: Consider adding another convenience method wrapper to construct this object whose signature is more closely
    #  aligned with get_anndata() params (i.e. "exploded" AxisQuery params).
    def __init__(
        self,
        experiment: soma.Experiment,
        measurement_name: str = "RNA",
        X_name: str = "raw",
        obs_query: soma.AxisQuery | None = None,
        var_query: soma.AxisQuery | None = None,
        obs_column_names: Sequence[str] = (),
        batch_size: int = 1,
        shuffle: bool = True,
        seed: int | None = None,
        return_sparse_X: bool = False,
        soma_chunk_size: int | None = 64,
        use_eager_fetch: bool = True,
        shuffle_chunk_count: int | None = 2000,
        encoders: list[Encoder] | None = None,
    ) -> None:
        r"""Construct a new ``ExperimentDataPipe``.

        Args:
            experiment:
                The :class:`tiledbsoma.Experiment` from which to read data.
            measurement_name:
                The name of the :class:`tiledbsoma.Measurement` to read. Defaults to ``"RNA"``.
            X_name:
                The name of the X layer to read. Defaults to ``"raw"``.
            obs_query:
                The query used to filter along the ``obs`` axis. If not specified, all ``obs`` and ``X`` data will
                be returned, which can be very large.
            var_query:
                The query used to filter along the ``var`` axis. If not specified, all ``var`` columns (genes/features)
                will be returned.
            obs_column_names:
                The names of the ``obs`` columns to return. If custom encoders are passed, this parameter must not be used,
                since the columns will be inferred automatically from the encoders.
            batch_size:
                The number of rows of ``obs`` and ``X`` data to return in each iteration. Defaults to ``1``. A value of
                ``1`` will result in :class:`torch.Tensor` of rank 1 being returns (a single row); larger values will
                result in :class:`torch.Tensor`\ s of rank 2 (multiple rows).
            shuffle:
                Whether to shuffle the ``obs`` and ``X`` data being returned. Defaults to ``True``.
                For performance reasons, shuffling is not performed globally across all rows, but rather in chunks.
                More specifically, we select ``shuffle_chunk_count`` non-contiguous chunks across all the observations
                in the query, concatenate the chunks and shuffle the associated observations.
                The randomness of the shuffling is therefore determined by the
                (``soma_chunk_size``, ``shuffle_chunk_count``) selection. The default values have been determined
                to yield a good trade-off between randomness and performance. Further tuning may be required for
                different type of models. Note that memory usage is correlated to the product
                ``soma_chunk_size * shuffle_chunk_count``.
            seed:
                The random seed used for shuffling. Defaults to ``None`` (no seed). This *must* be specified when using
                :class:`torch.nn.parallel.DistributedDataParallel` to ensure data partitions are disjoint across worker
                processes.
            return_sparse_X:
                Controls whether the ``X`` data is returned as a dense or sparse :class:`torch.Tensor`. As ``X`` data is
                very sparse, setting this to ``True`` will reduce memory usage, if the model supports use of sparse
                :class:`torch.Tensor`\ s. Defaults to ``False``, since sparse :class:`torch.Tensor`\ s are still
                experimental in PyTorch.
            soma_chunk_size:
                The number of ``obs``/``X`` rows to retrieve when reading data from SOMA. This impacts two aspects of
                this class's behavior: 1) The maximum memory utilization, with larger values providing
                better read performance, but also requiring more memory; 2) The granularity of the global shuffling
                step (see ``shuffle`` parameter for details). The default value of 64 works well in conjunction
                with the default ``shuffle_chunk_count`` value.
            use_eager_fetch:
                Fetch the next SOMA chunk of ``obs`` and ``X`` data immediately after a previously fetched SOMA chunk is made
                available for processing via the iterator. This allows network (or filesystem) requests to be made in
                parallel with client-side processing of the SOMA data, potentially improving overall performance at the
                cost of doubling memory utilization. Defaults to ``True``.
            shuffle_chunk_count:
                The number of contiguous blocks (chunks) of rows sampled to then concatenate and shuffle.
                Larger numbers correspond to more randomness per training batch.
                If ``shuffle == False``, this parameter is ignored. Defaults to ``2000``.
            encoders:
                Specify custom encoders to be used. If not specified, a LabelEncoder will be created and
                used for each column in ``obs_column_names``. If specified, only columns for which an encoder
                has been registered will be returned in the ``obs`` tensor. Each encoder needs to have a unique name.
                If this parameter is specified, the ``obs_column_names`` parameter must not be used,
                since the columns will be inferred automatically from the encoders.

        Lifecycle:
            experimental
        """
        self.exp_uri = experiment.uri
        self.aws_region = experiment.context.tiledb_config.get("vfs.s3.region")
        self.measurement_name = measurement_name
        self.layer_name = X_name
        self.obs_query = obs_query
        self.var_query = var_query
        self.obs_column_names = obs_column_names
        self.batch_size = batch_size
        self.return_sparse_X = return_sparse_X
        self.soma_chunk_size = soma_chunk_size
        self.use_eager_fetch = use_eager_fetch
        self._stats = Stats()
        self._encoders = encoders or []
        self._obs_joinids = None
        self._var_joinids = None
        self._shuffle_chunk_count = shuffle_chunk_count if shuffle else None
        self._shuffle_rng = np.random.default_rng(seed) if shuffle else None
        self._initialized = False
        self.max_process_mem_usage_bytes = 0

        if obs_column_names and encoders:
            raise ValueError(
                "Cannot specify both `obs_column_names` and `encoders`. If `encoders` are specified, columns will be inferred automatically."
            )

        if encoders:
            # Check if names are unique
            if len(encoders) != len({enc.name for enc in encoders}):
                raise ValueError("Encoders must have unique names")

            self.obs_column_names = list(dict.fromkeys(itertools.chain(*[enc.columns for enc in encoders])))

    def _init(self) -> None:
        if self._initialized:
            return

        pytorch_logger.debug("Initializing ExperimentDataPipe")

        with _open_experiment(self.exp_uri, self.aws_region) as exp:
            query = exp.axis_query(
                measurement_name=self.measurement_name,
                obs_query=self.obs_query,
                var_query=self.var_query,
            )

            # The to_numpy() call is a workaround for a possible bug in TileDB-SOMA:
            # https://github.com/single-cell-data/TileDB-SOMA/issues/1456
            self._obs_joinids = query.obs_joinids().to_numpy()
            self._var_joinids = query.var_joinids().to_numpy()

            self._encoders = self._build_obs_encoders(query)

        self._initialized = True

    @staticmethod
    def _subset_ids_to_partition(
        ids_chunked: list[npt.NDArray[np.int64]],
        partition_index: int,
        num_partitions: int,
    ) -> list[npt.NDArray[np.int64]]:
        """Returns a single partition of the obs_joinids_chunked (a 2D ndarray), based upon the current process's distributed rank and world
        size.
        """
        # subset to a single partition
        # typing does not reflect that is actually a List of 2D NDArrays
        partition_indices = np.array_split(range(len(ids_chunked)), num_partitions)
        partition = [ids_chunked[i] for i in partition_indices[partition_index]]

        if pytorch_logger.isEnabledFor(logging.DEBUG) and len(partition) > 0:
            pytorch_logger.debug(
                f"Process {os.getpid()} handling partition {partition_index + 1} of {num_partitions}, "
                f"partition_size={sum([len(chunk) for chunk in partition])}"
            )

        return partition

    @staticmethod
    def _compute_partitions(
        loader_partition: int,
        loader_partitions: int,
        dist_partition: int,
        num_dist_partitions: int,
    ) -> tuple[int, int]:
        # NOTE: Can alternately use a `worker_init_fn` to split among workers split workload
        total_partitions = num_dist_partitions * loader_partitions
        partition = dist_partition * loader_partitions + loader_partition
        return partition, total_partitions

    def __iter__(self) -> Iterator[ObsAndXDatum]:
        self._init()
        assert self._obs_joinids is not None
        assert self._var_joinids is not None

        if self.soma_chunk_size is None:
            # set soma_chunk_size to utilize ~1 GiB of RAM per SOMA chunk; assumes 95% X data sparsity, 8 bytes for the
            # X value and 8 bytes for the sparse matrix indices, and a 100% working memory overhead (2x).
            X_row_memory_size = 0.05 * len(self._var_joinids) * 8 * 3 * 2
            self.soma_chunk_size = int((1 * 1024**3) / X_row_memory_size)
        pytorch_logger.debug(f"Using {self.soma_chunk_size=}")

        if (
            self.return_sparse_X
            and torch.utils.data.get_worker_info()
            and torch.utils.data.get_worker_info().num_workers > 0
        ):
            raise NotImplementedError(
                "torch does not work with sparse tensors in multi-processing mode "
                "(see https://github.com/pytorch/pytorch/issues/20248)"
            )

        # chunk the obs joinids into batches of size soma_chunk_size
        obs_joinids_chunked = self._chunk_ids(self._obs_joinids, self.soma_chunk_size)

        # globally shuffle the chunks, if requested
        if self._shuffle_rng:
            self._shuffle_rng.shuffle(obs_joinids_chunked)

        # subset to a single partition, as needed for distributed training and multi-processing datat loading
        worker_info = torch.utils.data.get_worker_info()
        partition, partitions = self._compute_partitions(
            loader_partition=worker_info.id if worker_info else 0,
            loader_partitions=worker_info.num_workers if worker_info else 1,
            dist_partition=dist.get_rank() if dist.is_initialized() else 0,
            num_dist_partitions=dist.get_world_size() if dist.is_initialized() else 1,
        )
        obs_joinids_chunked_partition: list[npt.NDArray[np.int64]] = self._subset_ids_to_partition(
            obs_joinids_chunked, partition, partitions
        )

        with _open_experiment(self.exp_uri, self.aws_region) as exp:
            obs_and_x_iter = _ObsAndXIterator(
                obs=exp.obs,
                X=exp.ms[self.measurement_name].X[self.layer_name],
                obs_column_names=self.obs_column_names,
                obs_joinids_chunked=obs_joinids_chunked_partition,
                var_joinids=self._var_joinids,
                batch_size=self.batch_size,
                encoders=self._encoders,
                stats=self._stats,
                return_sparse_X=self.return_sparse_X,
                use_eager_fetch=self.use_eager_fetch,
                shuffle_rng=self._shuffle_rng,
                shuffle_chunk_count=self._shuffle_chunk_count,
            )

            yield from obs_and_x_iter

            self.max_process_mem_usage_bytes = obs_and_x_iter.max_process_mem_usage_bytes
            pytorch_logger.debug(
                "max process memory usage=" f"{self.max_process_mem_usage_bytes / (1024 ** 3):.3f} GiB"
            )

    @staticmethod
    def _chunk_ids(ids: npt.NDArray[np.int64], chunk_size: int) -> list[npt.NDArray[np.int64]]:
        num_chunks = max(1, ceil(len(ids) / chunk_size))
        pytorch_logger.debug(f"Shuffling {len(ids)} obs joinids into {num_chunks} chunks of {chunk_size}")
        return np.array_split(ids, num_chunks)

    def __len__(self) -> int:
        self._init()
        assert self._obs_joinids is not None

        div, rem = divmod(len(self._obs_joinids), self.batch_size)
        return div + bool(rem)

    def __getitem__(self, index: int) -> ObsAndXDatum:
        raise NotImplementedError("IterDataPipe can only be iterated")

    def _build_obs_encoders(self, query: soma.ExperimentAxisQuery) -> list[Encoder]:
        pytorch_logger.debug("Initializing encoders")

        encoders = []

        if "soma_joinid" not in self.obs_column_names:
            cols = ["soma_joinid", *self.obs_column_names]
        else:
            cols = list(self.obs_column_names)

        obs = query.obs(column_names=cols).concat().to_pandas()

        if self._encoders:
            # Fit all the custom encoders with obs
            for enc in self._encoders:
                enc.fit(obs)
                encoders.append(enc)
        else:
            # Create one LabelEncoder for each column, and fit it with obs
            for col in self.obs_column_names:
                enc = LabelEncoder(col)
                enc.fit(obs)
                encoders.append(enc)

        return encoders

    # TODO: This does not work in multiprocessing mode, as child process's stats are not collected
    def stats(self) -> Stats:
        """Get data loading stats for this :class:`cellxgene_census.experimental.ml.pytorch.ExperimentDataPipe`.

        Returns:
            The :class:`cellxgene_census.experimental.ml.pytorch.Stats` object for this
            :class:`cellxgene_census.experimental.ml.pytorch.ExperimentDataPipe`.

        Lifecycle:
            experimental
        """
        return self._stats

    @property
    def shape(self) -> tuple[int, int]:
        """Get the shape of the data that will be returned by this :class:`cellxgene_census.experimental.ml.pytorch.ExperimentDataPipe`.
        This is the number of obs (cell) and var (feature) counts in the returned data. If used in multiprocessing mode
        (i.e. :class:`torch.utils.data.DataLoader` instantiated with num_workers > 0), the obs (cell) count will reflect
        the size of the partition of the data assigned to the active process.

        Returns:
            A 2-tuple of ``int``s, for obs and var counts, respectively.

        Lifecycle:
            experimental
        """
        self._init()
        assert self._obs_joinids is not None
        assert self._var_joinids is not None

        return len(self._obs_joinids), len(self._var_joinids)

    @property
    def obs_encoders(self) -> Encoders:
        """Returns a dictionary of :class:`sklearn.preprocessing.LabelEncoder` objects, keyed on ``obs`` column names,
        which were used to encode the ``obs`` column values.

        These encoders can be used to decode the encoded values as follows:

        >>> exp_data_pipe.obs_encoders["<obs_attr_name>"].inverse_transform(encoded_values)

        Returns:
            A ``Dict[str, LabelEncoder]``, mapping column names to :class:`sklearn.preprocessing.LabelEncoder` objects.
        """
        self._init()
        assert self._encoders is not None

        return {enc.name: enc for enc in self._encoders}


# Note: must be a top-level function (and not a lambda), to play nice with multiprocessing pickling
def _collate_noop(x: Any) -> Any:
    return x


# TODO: Move into somacore.ExperimentAxisQuery
def experiment_dataloader(
    datapipe: pipes.IterDataPipe,
    num_workers: int = 0,
    **dataloader_kwargs: Any,
) -> DataLoader:
    """Factory method for :class:`torch.utils.data.DataLoader`. This method can be used to safely instantiate a
    :class:`torch.utils.data.DataLoader` that works with :class:`cellxgene_census.experimental.ml.pytorch.ExperimentDataPipe`,
    since some of the :class:`torch.utils.data.DataLoader` constructor parameters are not applicable when using a
    :class:`torchdata.datapipes.iter.IterDataPipe` (``shuffle``, ``batch_size``, ``sampler``, ``batch_sampler``,
    ``collate_fn``).

    Args:
        datapipe:
            An :class:`torchdata.datapipes.iter.IterDataPipe`, which can be an
            :class:`cellxgene_census.experimental.ml.pytorch.ExperimentDataPipe` or any other
            :class:`torchdata.datapipes.iter.IterDataPipe` that has been chained to the
            :class:`cellxgene_census.experimental.ml.pytorch.ExperimentDataPipe`.
        num_workers:
            Number of worker processes to use for data loading. If ``0``, data will be loaded in the main process.
        **dataloader_kwargs:
            Additional keyword arguments to pass to the :class:`torch.utils.data.DataLoader` constructor,
            except for ``shuffle``, ``batch_size``, ``sampler``, ``batch_sampler``, and ``collate_fn``, which are not
            supported when using :class:`cellxgene_census.experimental.ml.pytorch.ExperimentDataPipe`.

    Returns:
        A :class:`torch.utils.data.DataLoader`.

    Raises:
        ValueError: if any of the ``shuffle``, ``batch_size``, ``sampler``, ``batch_sampler``, or ``collate_fn`` params
            are passed as keyword arguments.

    Lifecycle:
        experimental
    """
    unsupported_dataloader_args = [
        "shuffle",
        "batch_size",
        "sampler",
        "batch_sampler",
        "collate_fn",
    ]
    if set(unsupported_dataloader_args).intersection(dataloader_kwargs.keys()):
        raise ValueError(f"The {','.join(unsupported_dataloader_args)} DataLoader params are not supported")

    if num_workers > 0:
        _init_multiprocessing()

    return DataLoader(
        datapipe,
        batch_size=None,  # batching is handled by our ExperimentDataPipe
        num_workers=num_workers,
        # avoid use of default collator, which adds an extra (3rd) dimension to the tensor batches
        collate_fn=_collate_noop,
        # shuffling is handled by our ExperimentDataPipe
        shuffle=False,
        **dataloader_kwargs,
    )


def _init_multiprocessing() -> None:
    """Ensures use of "spawn" for starting child processes with multiprocessing.

    Forked processes are known to be problematic:
      https://pytorch.org/docs/stable/notes/multiprocessing.html#avoiding-and-fighting-deadlocks
    Also, CUDA does not support forked child processes:
      https://pytorch.org/docs/stable/notes/multiprocessing.html#cuda-in-multiprocessing

    """
    torch.multiprocessing.set_start_method("fork", force=True)
    orig_start_method = torch.multiprocessing.get_start_method()
    if orig_start_method != "spawn":
        if orig_start_method:
            pytorch_logger.warning(
                "switching torch multiprocessing start method from "
                f'"{torch.multiprocessing.get_start_method()}" to "spawn"'
            )
        torch.multiprocessing.set_start_method("spawn", force=True)



# Section: cellxgene_census-src-cellxgene_census-_open

# Copyright (c) 2022, Chan Zuckerberg Initiative
#
# Licensed under the MIT License.

"""Open census and related datasets.

Contains methods to open publicly hosted versions of Census object and access its source datasets.
"""

import logging
import os.path
import urllib.parse
from typing import Any, get_args

import s3fs
import tiledbsoma as soma
from fsspec.callbacks import NoOpCallback, TqdmCallback

from ._release_directory import (
    CensusLocator,
    CensusMirror,
    Provider,
    ResolvedCensusLocator,
    _get_census_mirrors,
    get_census_version_description,
)
from ._util import _uri_join, _user_agent

DEFAULT_CENSUS_VERSION = "stable"

DEFAULT_S3FS_KWARGS = {
    "anon": True,
    "cache_regions": True,
}
DEFAULT_TILEDB_CONFIGURATION: dict[str, Any] = {
    # https://docs.tiledb.com/main/how-to/configuration#configuration-parameters
    "py.init_buffer_bytes": 1 * 1024**3,
    "soma.init_buffer_bytes": 1 * 1024**3,
    # S3 requests should not be signed, since we want to allow anonymous access
    "vfs.s3.no_sign_request": "true",
    "vfs.s3.region": "us-west-2",
}

api_logger = logging.getLogger("cellxgene_census")
api_logger.setLevel(logging.INFO)
api_logger.addHandler(logging.StreamHandler())


def _assert_mirror_supported(mirror: CensusMirror) -> None:
    """Verifies if the mirror is supported by this version of the census API.
    This method provides a proper error message in case an old version of the census
    tries to connect to an unsupported mirror.
    """
    if mirror["provider"] not in get_args(Provider):
        raise ValueError(
            f"Unsupported mirror provider: {mirror['provider']}. Try upgrading the cellxgene-census package to the latest version."
        )


def _resolve_census_locator(locator: CensusLocator, mirror: CensusMirror) -> ResolvedCensusLocator:
    _assert_mirror_supported(mirror)

    if locator.get("relative_uri"):
        uri = _uri_join(mirror["base_uri"], locator["relative_uri"])
        region = mirror["region"]
    else:
        uri = locator["uri"]
        region = locator.get("s3_region")
    return ResolvedCensusLocator(uri=uri, region=region, provider=mirror["provider"])


def _open_soma(
    locator: ResolvedCensusLocator,
    context: soma.options.SOMATileDBContext | None = None,
) -> soma.Collection:
    """Private. Merge config defaults and return open census as a soma Collection/context."""
    # if no user-defined context, cellxgene_census defaults take precedence over SOMA defaults
    context = context or get_default_soma_context()

    # A locator's S3 bucket and region are intrinsically coupled, so ensure that the context always uses the
    # locator's region, even if the passed-in context specifies an alternate S3 region
    if locator["provider"] == "S3":
        context = context.replace(tiledb_config={"vfs.s3.region": locator.get("region")})

    return soma.open(locator["uri"], mode="r", soma_type=soma.Collection, context=context)


def get_default_soma_context(tiledb_config: dict[str, Any] | None = None) -> soma.options.SOMATileDBContext:
    """Return a :class:`tiledbsoma.SOMATileDBContext` with sensible defaults that can be further customized by the
    user. The customized context can then be passed to :func:`cellxgene_census.open_soma` with the ``context``
    argument or to :meth:`somacore.SOMAObject.open` with the ``context`` argument, such as
    :meth:`tiledbsoma.Experiment.open`. Use the :meth:`tiledbsoma.SOMATileDBContext.replace` method on the returned
    object to customize its settings further.

    Args:
        tiledb_config:
            A dictionary of TileDB configuration parameters. If specified, the parameters will override the
            defaults. If not specified, the default configuration will be returned.

    Returns:
        A :class:`tiledbsoma.SOMATileDBContext` object with sensible defaults.

    Examples:
        To reduce the amount of memory used by TileDB-SOMA I/O operations:

        .. highlight:: python
        .. code-block:: python

            ctx = cellxgene_census.get_default_soma_context(
                tiledb_config={"py.init_buffer_bytes": 128 * 1024**2,
                               "soma.init_buffer_bytes": 128 * 1024**2})
            c = census.open_soma(uri="s3://my-private-bucket/census/soma", context=ctx)

        To access a copy of the Census located in a private bucket that is located in a different S3 region, use:

        .. highlight:: python
        .. code-block:: python

            ctx = cellxgene_census.get_default_soma_context(
                tiledb_config={"vfs.s3.no_sign_request": "false",
                               "vfs.s3.region": "us-east-1"})
            c = census.open_soma(uri="s3://my-private-bucket/census/soma", context=ctx)

    Lifecycle:
        experimental
    """
    tiledb_config = dict(
        DEFAULT_TILEDB_CONFIGURATION, **{"vfs.s3.custom_headers.User-Agent": _user_agent()}, **(tiledb_config or {})
    )
    return soma.options.SOMATileDBContext().replace(tiledb_config=tiledb_config)


def open_soma(
    *,
    census_version: str | None = DEFAULT_CENSUS_VERSION,
    mirror: str | None = None,
    uri: str | None = None,
    tiledb_config: dict[str, Any] | None = None,
    context: soma.options.SOMATileDBContext | None = None,
) -> soma.Collection:
    """Open the Census by version or URI.

    Args:
        census_version:
            The version of the Census, e.g. ``"latest"`` or ``"stable"``. Defaults to ``"stable"``.
        mirror:
            The mirror used to retrieve the Census. If not specified, a suitable mirror
            will be chosen automatically.
        uri:
            The URI containing the Census SOMA objects. If specified, will take precedence
            over ``census_version`` parameter.
        tiledb_config:
            A dictionary of TileDB configuration parameters that will be used to open the SOMA object. Optional,
            defaults to ``None``. If specified, the parameters will override the default settings specified by
            ``get_default_soma_context().tiledb_config``. Only one of the ``tiledb_config`` and ``context`` params
            can be specified.
        context:
            A custom :class:`tiledbsoma.SOMATileDBContext` that will be used to open the SOMA object.
            Optional, defaults to ``None``. Only one of the ``tiledb_config`` and ``context`` params can be specified.

    Returns:
        A :class:`tiledbsoma.Collection` object containing the top-level census.
        It can be used as a context manager, which will automatically close upon exit.

    Raises:
        ValueError: if the census cannot be found, the URI cannot be opened, neither a URI
            or a version are specified, or an invalid mirror is provided.

    See Also:
        - :func:`get_source_h5ad_uri`: Look up the location of the source H5AD.
        - :func:`get_default_soma_context`: Get a default SOMA context that can be updated with custom settings and used
          as the ``context`` argument.

    Lifecycle:
        maturing

    Examples:
        Open the default Census version, using a context manager which will automatically
        close the Census upon exit of the context.

        >>> with cellxgene_census.open_soma() as census:
                ...

        Open and close:

        >>> census = cellxgene_census.open_soma()
            ...
            census.close()

        Open a specific Census by version:

        >>> with cellxgene_census.open_soma("2022-12-31") as census:
                ...

        Open a Census by S3 URI, rather than by version.

        >>> with cellxgene_census.open_soma(uri="s3://bucket/path") as census:
                ...

        Open a Census by path (file:// URI), rather than by version.

        >>> with cellxgene_census.open_soma(uri="/tmp/census") as census:
                ...

        Open a Census using a mirror.

        >>> with cellxgene_census.open_soma(mirror="s3-us-west-2") as census:
                ...

        Open a Census with a custom TileDB configuration setting.

        >>> with cellxgene_census.open_soma(tiledb_config={"py.init_buffer_bytes": 128 * 1024**2}) as census:
                ...
    """
    if tiledb_config is not None and context is not None:
        raise ValueError("Only one of tiledb_config and context can be specified.")

    if tiledb_config is not None:
        context = get_default_soma_context(tiledb_config=tiledb_config)

    if uri is not None:
        return _open_soma({"uri": uri, "region": None, "provider": "unknown"}, context)

    if census_version is None:
        raise ValueError("Must specify either a census version or an explicit URI.")

    mirrors = _get_census_mirrors()
    selected_mirror: CensusMirror
    if mirror is not None:
        if mirror not in mirrors:
            raise ValueError("Mirror not found.")
        selected_mirror = mirrors[mirror]  # type: ignore
    else:
        selected_mirror = mirrors[mirrors["default"]]  # type: ignore

    # TODO: Consider raising exceptions instead of issuing warnings, possibly introducing a "strict" mode to control the
    #  behavior
    try:
        description = get_census_version_description(census_version)  # raises
    except ValueError:
        raise ValueError(
            f'The "{census_version}" Census version is not valid. Use get_census_version_directory() to retrieve '
            f"available versions."
        ) from None

    if description.get("flags", {}).get("retracted", False):
        api_logger.warning(
            f"The \"{census_version}\" Census version has been retracted!\n{description['retraction']}."
            f'Use "stable" or "latest", or use get_census_version_directory() to retrieve valid versions.'
        )
    elif census_version == "stable":
        api_logger.info(
            f"The \"{census_version}\" release is currently {description['release_build']}. Specify "
            f"'census_version=\"{description['release_build']}\"' in future calls to open_soma() to ensure data "
            "consistency."
        )

    locator = _resolve_census_locator(description["soma"], selected_mirror)

    return _open_soma(locator, context)


def get_source_h5ad_uri(dataset_id: str, *, census_version: str = DEFAULT_CENSUS_VERSION) -> CensusLocator:
    """Open the named version of the census, and return the URI for the ``dataset_id``. This
    does not guarantee that the H5AD exists or is accessible to the user.

    Args:
        dataset_id:
            The ``dataset_id`` of interest.
        census_version:
            The census version. Defaults to ``"stable"``.

    Returns:
        A :class:`cellxgene_census._release_directory.CensusLocator` object that contains the URI and optional S3 region
        for the source H5AD.

    Raises:
        KeyError: if either ``dataset_id`` or ``census_version`` do not exist.

    Lifecycle:
        maturing

    Examples:
        >>> cellxgene_census.get_source_h5ad_uri("cb5efdb0-f91c-4cbd-9ad4-9d4fa41c572d")
        {'uri': 's3://cellxgene-data-public/cell-census/2022-12-01/h5ads/cb5efdb0-f91c-4cbd-9ad4-9d4fa41c572d.h5ad',
        's3_region': 'us-west-2'}
    """
    description = get_census_version_description(census_version)  # raises

    # For h5ads, it makes sense to use the default mirror, since the artifacts themselves won't be mirrored
    mirrors = _get_census_mirrors()
    selected_mirror: CensusMirror = mirrors[mirrors["default"]]  # type: ignore
    census_locator = _resolve_census_locator(description["soma"], selected_mirror)

    census = _open_soma(census_locator)
    dataset = census["census_info"]["datasets"].read(value_filter=f"dataset_id == '{dataset_id}'").concat().to_pandas()
    if len(dataset) == 0:
        raise KeyError("Unknown dataset_id")

    locator = description["h5ads"].copy()
    h5ads_base_uri = locator["uri"]
    dataset_h5ad_path = dataset.dataset_h5ad_path.iloc[0]
    locator["uri"] = _uri_join(h5ads_base_uri, dataset_h5ad_path)
    return locator


def download_source_h5ad(
    dataset_id: str, to_path: str, *, census_version: str = DEFAULT_CENSUS_VERSION, progress_bar: bool = True
) -> None:
    """Download the source H5AD dataset, for the given `dataset_id`, to the user-specified
    file name.

    Args:
        dataset_id
            Fetch the source (original) H5AD associated with this ``dataset_id``.
        to_path:
            The file name where the downloaded H5AD will be written. Must not already exist.
        census_version:
            The census version name. Defaults to ``"stable"``.
        progress_bar:
            Whether to display a progress bar. Defaults to ``True``.

    Raises:
        ValueError: if the path already exists (i.e., will not overwrite an existing file), or is not a file.

    Lifecycle:
        maturing

    See Also:
        :func:`get_source_h5ad_uri`: Look up the location of the source H5AD.

    Examples:
        >>> download_source_h5ad("8e47ed12-c658-4252-b126-381df8d52a3d", to_path="/tmp/data.h5ad")
    """
    if os.path.exists(to_path):
        raise ValueError("Path exists - will not overwrite existing file.")
    if to_path.endswith("/"):
        raise ValueError("Specify to_path as a file name, not a directory name.")

    if progress_bar:
        callback = TqdmCallback(
            tqdm_kwargs={"unit": "B", "unit_scale": True, "unit_divisor": 1024, "desc": "Downloading"}
        )
    else:
        callback = NoOpCallback()

    locator = get_source_h5ad_uri(dataset_id, census_version=census_version)
    protocol = urllib.parse.urlparse(locator["uri"]).scheme
    assert protocol == "s3"

    fs = s3fs.S3FileSystem(
        config_kwargs={"user_agent": _user_agent()},
        **DEFAULT_S3FS_KWARGS,
    )
    fs.get_file(
        locator["uri"],
        to_path,
        callback=callback,
    )



# Section: cellxgene_census-src-cellxgene_census-experimental-pp-_stats

from __future__ import annotations

from collections.abc import Generator
from concurrent import futures
from typing import Any

import numpy as np
import numpy.typing as npt
import pandas as pd
import tiledbsoma as soma

from ..util._eager_iter import _EagerIterator
from ._online import MeanAccumulator, MeanVarianceAccumulator


def mean_variance(
    query: soma.ExperimentAxisQuery,
    layer: str = "raw",
    axis: int = 0,
    calculate_mean: bool = False,
    calculate_variance: bool = False,
    ddof: int = 1,
    nnz_only: bool = False,
) -> pd.DataFrame:
    """Calculate  mean and/or variance along the ``obs`` axis from query results. Calculations are done in an accumulative
    chunked fashion. For the mean and variance calculations, the total number of elements (N) is, by default, the
    corresponding dimension size: for column-wise calculations (``axis = 0``) N is number of rows, for row-wise
    calculations (``axis = 1``) N is number of columns. For metrics calculated only on nnz (explicitly stored) values of
    the sparse matrix, specify ``nnz_only=True``.

    Args:
        query:
            A :class:`tiledbsoma.ExperimentAxisQuery`, specifying the ``obs``/``var`` selection over which mean and
            variance are calculated.
        layer:
            X layer used, e.g., ``"raw"``.
        axis:
           Axis or axes along which the statistics are computed.
        calculate_mean:
            If ``True`` it calculates mean, otherwise skips calculation.
        calculate_variance:
            If ``True`` it calculates variance, otherwise skips calculation.
        ddof:
            "Delta Degrees of Freedom": the divisor used in the calculation for variance is ``N - ddof``, where ``N``
            represents the number of elements.
        nnz_only:
            If ``True`` mean and variance will only be calculated over explicitly stored values in the sparse matrix.
            Defaults to ``False``.

    Returns:
        :class:`pandas.DataFrame` indexed by the ``soma_joinid`` and with columns ``mean`` (if
        ``calculate_mean = True``), and ``variance`` (if ``calculate_variance = True``).

    Lifecycle:
        experimental
    """
    if axis not in (0, 1):
        raise ValueError("axis must be 0 or 1")

    if calculate_mean is False and calculate_variance is False:
        raise ValueError("At least one of `calculate_mean` or `calculate_variance` must be True")

    if query.n_obs == 0 or query.n_vars == 0:
        raise ValueError("The query cannot yield an empty result")

    n_dim_0 = query.n_obs if axis == 1 else query.n_vars
    n_dim_1 = query.n_vars if axis == 1 else query.n_obs

    n_batches = 1
    n_samples = np.array([n_dim_1], dtype=np.int64)

    idx = pd.Index(
        data=query.obs_joinids() if axis == 1 else query.var_joinids(),
        name="soma_joinid",
    )

    def iterate() -> Generator[tuple[npt.NDArray[np.int64], Any], None, None]:
        with futures.ThreadPoolExecutor(max_workers=1) as pool:  # Note: _EagerIterator only supports one thread
            for arrow_tbl in _EagerIterator(query.X(layer).tables(), pool=pool):
                dim = idx.get_indexer(arrow_tbl[f"soma_dim_{1-axis}"].to_numpy())
                data = arrow_tbl["soma_data"].to_numpy()
                yield dim, data

    result = pd.DataFrame(
        index=idx,
    )

    if calculate_variance:
        mvn = MeanVarianceAccumulator(n_batches, n_samples, n_dim_0, ddof=ddof, nnz_only=nnz_only)
        for dim, data in iterate():
            mvn.update(dim, data)
        _, _, all_u, all_var = mvn.finalize()
        if calculate_mean:
            result["mean"] = all_u
        result["variance"] = all_var
    else:
        mn = MeanAccumulator(n_dim_1, n_dim_0, nnz_only=nnz_only)
        for dim, data in iterate():
            mn.update(dim, data)
        all_u = mn.finalize()
        result["mean"] = all_u

    return result



# Section: cellxgene_census-tests-experimental-ml-huggingface-test_geneformer

import sys

import datasets
import pytest
import tiledbsoma
from py.path import local as Path
from scipy.stats import spearmanr

import cellxgene_census

try:
    from geneformer import TranscriptomeTokenizer

    from cellxgene_census.experimental.ml.huggingface import GeneformerTokenizer
except ImportError:
    # this should only occur when not running `experimental`-marked tests
    pass


CENSUS_VERSION_FOR_GENEFORMER_TESTS = "2023-12-15"


@pytest.mark.skip("Needs to be investigated.")
@pytest.mark.experimental
@pytest.mark.live_corpus
def test_GeneformerTokenizer_correctness(tmpdir: Path) -> None:
    """
    Test that GeneformerTokenizer produces the same token sequences as the original
    geneformer.TranscriptomeTokenizer (modulo a small tolerance on Spearman rank correlation)
    """
    # causes deterministic selection of roughly 1,000 cells:
    MODULUS = 32768
    # minimum Spearman rank correlation to consider token sequences effectively identical; this
    # allows for rare, slight differences in token sequences possibly arising from unstable sorting
    # and/or minor numerical precision differences in lowly-expressed genes.
    RHO_THRESHOLD = 0.99
    # notwithstanding RHO_THRESHOLD, we'll check that almost all token sequences are -exactly-
    # identical.
    EXACT_THRESHOLD = 0.98

    with cellxgene_census.open_soma(census_version=CENSUS_VERSION_FOR_GENEFORMER_TESTS) as census:
        human = census["census_data"]["homo_sapiens"]
        # read obs dataframe to get soma_joinids of all primary cells
        obs_df = (
            human.obs.read(column_names=["soma_joinid"], value_filter="is_primary_data == True").concat().to_pandas()
        )
        # select those with soma_joinid == 0 (mod MODULUS)
        cell_ids = [it for it in obs_df["soma_joinid"].tolist() if it % MODULUS == 0]

        # run our GeneformerTokenizer on them
        with GeneformerTokenizer(
            human,
            obs_query=tiledbsoma.AxisQuery(coords=(cell_ids,)),
        ) as tokenizer:
            test_tokens = [it["input_ids"] for it in tokenizer.build()]

        # write h5ad for use with geneformer.TranscriptomeTokenizer
        ad = cellxgene_census.get_anndata(
            census,
            "homo_sapiens",
            X_name="raw",
            obs_coords=cell_ids,
            column_names=tiledbsoma.AxisColumnNames(var=["feature_id"]),
        )
        ad.var.rename(columns={"feature_id": "ensembl_id"}, inplace=True)
        ad.obs["n_counts"] = ad.X.sum(axis=1)
        h5ad_dir = tmpdir.join("h5ad")
        h5ad_dir.mkdir()
        ad.write_h5ad(h5ad_dir.join("tokenizeme.h5ad"))
        # run geneformer.TranscriptomeTokenizer to get "true" tokenizations
        # see: https://huggingface.co/ctheodoris/Geneformer/blob/main/geneformer/tokenizer.py
        TranscriptomeTokenizer({}).tokenize_data(h5ad_dir, str(tmpdir), "tk", file_format="h5ad")
        true_tokens = [it["input_ids"] for it in datasets.load_from_disk(tmpdir.join("tk.dataset"))]

        # check GeneformerTokenizer sequences against geneformer.TranscriptomeTokenizer's
        assert len(test_tokens) == len(cell_ids)
        assert len(true_tokens) == len(cell_ids)
        identical = 0
        for i, cell_id in enumerate(cell_ids):
            assert len(test_tokens[i]) == len(true_tokens[i])
            rho, _ = spearmanr(test_tokens[i], true_tokens[i])
            if rho < RHO_THRESHOLD:
                # token sequences are too dissimilar; assert exact identity so that pytest -vv will
                # show the complete diff:
                assert (
                    test_tokens[i] == true_tokens[i]
                ), f"Discrepant token sequences for cell soma_joinid={cell_id}; Spearman rho={rho}"
            elif test_tokens[i] == true_tokens[i]:
                identical += 1
        assert identical / len(cell_ids) >= EXACT_THRESHOLD


@pytest.mark.skipif(sys.version_info < (3, 10), reason="requires python3.10 or higher")
@pytest.mark.experimental
@pytest.mark.live_corpus
def test_GeneformerTokenizer_docstring_example() -> None:
    with cellxgene_census.open_soma(census_version=CENSUS_VERSION_FOR_GENEFORMER_TESTS) as census:
        with GeneformerTokenizer(
            census["census_data"]["homo_sapiens"],
            # set obs_query to define some subset of Census cells:
            obs_query=tiledbsoma.AxisQuery(value_filter="is_primary_data == True and tissue_general == 'tongue'"),
            obs_attributes=(
                "soma_joinid",
                "cell_type_ontology_term_id",
            ),
        ) as tokenizer:
            dataset = tokenizer.build()
            assert len(dataset) == 15020
            assert sum(it.length for it in dataset.to_pandas().itertuples()) == 27798388



# Section: cellxgene_census-tests-test_util

import re

import pytest

import cellxgene_census
from cellxgene_census._util import _extract_census_version, _uri_join


def test_uri_join() -> None:
    assert _uri_join("https://foo/", "bar") == "https://foo/bar"
    assert _uri_join("https://foo/a", "bar") == "https://foo/bar"
    assert _uri_join("https://foo/a/", "bar") == "https://foo/a/bar"
    assert _uri_join("https://foo/", "a/b") == "https://foo/a/b"
    assert _uri_join("https://foo/", "a/b/") == "https://foo/a/b/"

    assert _uri_join("https://foo/?a=b#99", "a?b=c#1") == "https://foo/a?b=c#1"

    assert _uri_join("http://foo/bar/", "a") == "http://foo/bar/a"
    assert _uri_join("https://foo/bar/", "a") == "https://foo/bar/a"
    assert _uri_join("s3://foo/bar/", "a") == "s3://foo/bar/a"
    assert _uri_join("foo/bar", "a") == "foo/a"
    assert _uri_join("/foo/bar", "a") == "/foo/a"
    assert _uri_join("file://foo/bar", "a") == "file://foo/a"
    assert _uri_join("file:///foo/bar", "a") == "file:///foo/a"

    assert _uri_join("https://foo/bar", "https://a/b") == "https://a/b"


@pytest.mark.live_corpus
def test_extract_census_version() -> None:
    """Ensures that extracting the Census version from a Collection object does not break"""

    pattern = r"^\d{4}-\d{2}-\d{2}$"

    with cellxgene_census.open_soma(census_version="stable") as census:
        assert census is not None
        version = _extract_census_version(census)
        assert re.match(pattern, version)

    with cellxgene_census.open_soma(census_version="latest") as census:
        assert census is not None
        version = _extract_census_version(census)
        assert re.match(pattern, version)



# Section: notebooks-api_demo-census_compute_over_X

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




# Section: cellxgene_census-src-cellxgene_census-_testing-logger_proxy

"""This module defines a plugin class that logs each request to a logfile.

This class needs to be importable by the proxy server which runs in a separate process.
See the user agent tests for usage.
"""

import json
import traceback
from pathlib import Path

import proxy
from proxy.common.flag import flags

flags.add_argument(
    "--request-log-file",
    type=str,
    default="",
    help="Where to log the requests to.",
)


class RequestLoggerPlugin(proxy.http.proxy.HttpProxyBasePlugin):  # type: ignore
    def handle_client_request(self, request: proxy.http.parser.HttpParser) -> proxy.http.parser.HttpParser:
        # If anything fails in here, it just fails to respond
        try:
            with Path(self.flags.request_log_file).open("a") as f:
                record = {
                    "method": request.method.decode(),
                    "url": str(request._url),
                }

                if request.headers:
                    record["headers"] = {k2.decode().lower(): v.decode() for _, (k2, v) in request.headers.items()}
                f.write(f"{json.dumps(record)}\n")
        except Exception as e:
            # Making sure there is some visible output
            print(repr(e))
            traceback.print_exception(e)
            raise e
        return request



# Section: cellxgene_census-src-cellxgene_census-experimental-ml-huggingface-__init__

"""An API to facilitate using Hugging Face ML tools with the CZI Science CELLxGENE Census."""

from .cell_dataset_builder import CellDatasetBuilder
from .geneformer_tokenizer import GeneformerTokenizer

__all__ = [
    "CellDatasetBuilder",
    "GeneformerTokenizer",
]



# Section: notebooks-api_demo-census_citation_generation

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




# Section: notebooks-analysis_demo-comp_bio_scvi_model_use

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




# Section: cellxgene_census-tests-experimental-pp-test_online

import numpy as np
import numpy.ma as ma
import numpy.typing as npt
import pytest
from scipy import sparse

from cellxgene_census.experimental.pp._online import (
    CountsAccumulator,
    MeanAccumulator,
    MeanVarianceAccumulator,
)


def allclose(a: npt.NDArray[np.float64], b: npt.NDArray[np.float64]) -> bool:
    return np.allclose(a, b, atol=1e-5, rtol=1e-2, equal_nan=True)


@pytest.fixture
def matrix(m: int, n: int) -> sparse.coo_matrix:
    m = 100 * sparse.random(
        m,
        n,
        density=0.1,
        format="coo",
        dtype=np.float32,
        random_state=np.random.default_rng(),
    )
    m.row.flags.writeable = False  # type: ignore[attr-defined]
    m.col.flags.writeable = False  # type: ignore[attr-defined]
    m.data.flags.writeable = False  # type: ignore[attr-defined]
    return m


@pytest.mark.experimental
@pytest.mark.parametrize("stride", [101, 53])
@pytest.mark.parametrize("n_batches", [1, 3, 11, 101])
@pytest.mark.parametrize("m,n", [(1200, 511), (100001, 57)])
@pytest.mark.parametrize("ddof", [0, 1, 100])
@pytest.mark.parametrize("nnz_only", [True, False])
def test_meanvar(matrix: sparse.coo_matrix, n_batches: int, stride: int, ddof: int, nnz_only: bool) -> None:
    rng = np.random.default_rng()
    batches_prob = rng.random(n_batches)
    batches_prob /= batches_prob.sum()
    batches = rng.choice(n_batches, matrix.shape[0], p=batches_prob)
    batches.flags.writeable = False

    batch_id, batch_count = np.unique(batches, return_counts=True)
    n_samples = np.zeros((n_batches,), dtype=np.int64)
    n_samples[batch_id] = batch_count
    assert n_samples.sum() == matrix.shape[0]
    assert len(n_samples) == n_batches

    # nnz_only only if there is a single batch
    should_nnz_only = nnz_only and n_batches == 1

    olmv = MeanVarianceAccumulator(n_batches, n_samples, matrix.shape[1], ddof, nnz_only=should_nnz_only)
    for i in range(0, matrix.nnz, stride):
        batch_vec = batches[matrix.row[i : i + stride]] if n_batches > 1 else None
        olmv.update(matrix.col[i : i + stride], matrix.data[i : i + stride], batch_vec)
    batches_u, batches_var, all_u, all_var = olmv.finalize()

    assert isinstance(all_u, np.ndarray)
    assert isinstance(all_var, np.ndarray)
    assert isinstance(batches_u, np.ndarray)
    assert isinstance(batches_var, np.ndarray)

    assert all_u.shape == (matrix.shape[1],)
    assert all_var.shape == (matrix.shape[1],)
    assert batches_u.shape == (n_batches, matrix.shape[1])
    assert batches_var.shape == (n_batches, matrix.shape[1])

    dense = matrix.toarray()

    if should_nnz_only:
        mask = np.ones(matrix.shape)
        r, c = matrix.tolil().nonzero()
        for x, y in zip(r, c):
            mask[x, y] = 0
        masked: np.ma.MaskedArray[tuple[int, ...], np.dtype[np.float64]] = ma.masked_array(dense, mask=mask)  # type: ignore[no-untyped-call]
        mean = masked.mean(axis=0)  # type: ignore[no-untyped-call]
        assert allclose(all_u, mean)

        nv = masked.var(axis=0, ddof=ddof)  # type: ignore[no-untyped-call]
        assert allclose(all_var, nv)

    else:
        assert allclose(all_u, dense.mean(axis=0))
        assert allclose(all_var, dense.var(axis=0, ddof=ddof, dtype=np.float64))

        matrix = matrix.tocsr()  # COO not conveniently indexable
        for batch in range(n_batches):
            dense = matrix[batches == batch, :].toarray()
            assert allclose(batches_u[batch], dense.mean(axis=0))
            assert allclose(batches_var[batch], dense.var(axis=0, ddof=ddof, dtype=np.float64))


@pytest.mark.experimental
@pytest.mark.parametrize("m,n", [(1200, 511), (100001, 57)])
def test_meanvar_nnz_only_batches_fails(matrix: sparse.coo_matrix) -> None:
    n_batches = 10
    nnz_only = True
    n_samples = np.zeros((n_batches,), dtype=np.int64)
    ddof = 1
    with pytest.raises(ValueError):
        MeanVarianceAccumulator(n_batches, n_samples, matrix.shape[1], ddof, nnz_only=nnz_only)


@pytest.mark.experimental
@pytest.mark.parametrize("stride", [101, 53])
@pytest.mark.parametrize("m,n", [(1200, 511), (100001, 57)])
def test_mean(matrix: sparse.coo_matrix, stride: int) -> None:
    m_acc = MeanAccumulator(matrix.shape[0], matrix.shape[1])
    for i in range(0, matrix.nnz, stride):
        m_acc.update(matrix.col[i : i + stride], matrix.data[i : i + stride])
    u = m_acc.finalize()

    assert isinstance(u, np.ndarray)
    assert u.shape == (matrix.shape[1],)

    dense = matrix.toarray()
    assert allclose(u, dense.mean(axis=0))


@pytest.mark.experimental
@pytest.mark.parametrize("stride", [101, 53])
@pytest.mark.parametrize("n_batches", [1, 3, 11, 101])
@pytest.mark.parametrize("m,n", [(1200, 511), (100001, 57)])
def test_counts(matrix: sparse.coo_matrix, n_batches: int, stride: int) -> None:
    rng = np.random.default_rng()
    batches_prob = rng.random(n_batches)
    batches_prob /= batches_prob.sum()
    batches = rng.choice(n_batches, matrix.shape[0], p=batches_prob)
    batches.flags.writeable = False

    batch_id, batch_count = np.unique(batches, return_counts=True)
    n_samples = np.zeros((n_batches,), dtype=np.int64)
    n_samples[batch_id] = batch_count
    assert n_samples.sum() == matrix.shape[0]
    assert len(n_samples) == n_batches

    clip_val = 50 * np.random.rand(n_batches, matrix.shape[1])

    ca = CountsAccumulator(n_batches, matrix.shape[1], clip_val)
    for i in range(0, matrix.nnz, stride):
        batch_vec = batches[matrix.row[i : i + stride]] if n_batches > 1 else None
        ca.update(matrix.col[i : i + stride], matrix.data[i : i + stride], batch_vec)
    counts_sum, counts_squared_sum = ca.finalize()

    assert isinstance(counts_sum, np.ndarray)
    assert isinstance(counts_squared_sum, np.ndarray)

    assert counts_sum.shape == (n_batches, matrix.shape[1])
    assert counts_squared_sum.shape == (n_batches, matrix.shape[1])

    matrix = matrix.tocsr()  # COO not conveniently indexable
    for batch in range(n_batches):
        dense = matrix[batches == batch, :].toarray()
        assert allclose(counts_sum[batch], np.minimum(dense, clip_val[batch]).sum(axis=0))
        assert allclose(
            counts_squared_sum[batch],
            (np.minimum(dense, clip_val[batch]) ** 2).sum(axis=0),
        )


@pytest.mark.experimental
def test_mean_fails_no_variables_or_samples() -> None:
    with pytest.raises(ValueError, match=r"No samples provided - can't calculate mean."):
        MeanAccumulator(n_samples=0, n_variables=100)
    with pytest.raises(ValueError, match=r"No variables provided - can't calculate mean."):
        MeanAccumulator(n_samples=1000, n_variables=0)



# Section: notebooks-api_demo-census_embedding_search

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




# Section: notebooks-api_demo-census_datasets

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




# Section: notebooks-analysis_demo-comp_bio_summarize_axis_query

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




# Section: cellxgene_census-src-cellxgene_census-_release_directory

# Copyright (c) 2022, Chan Zuckerberg Initiative
#
# Licensed under the MIT License.

"""Versioning of Census builds.

Methods to retrieve information about versions of the publicly hosted Census object.
"""

from collections import OrderedDict
from typing import Any, Literal, cast

import requests
from typing_extensions import NotRequired, TypedDict

from cellxgene_census._util import _user_agent

"""
The following types describe the expected directory of Census builds, used
to bootstrap all data location requests.
"""
CensusVersionName = str  # census version name, e.g., "release-99", "2022-10-01-test", etc.


class CensusLocator(TypedDict):
    """A locator for a Census resource.

    Args:
        uri:
            Absolute resource URI (deprecated: only used in census < 1.6.0).
        relative_uri:
            Resource URI (relative).
        s3_region:
             If an S3 URI, has optional region (deprecated: only used in census < 1.6.0).
    """

    uri: str
    relative_uri: str
    s3_region: str | None


class CensusVersionRetraction(TypedDict):
    """A retraction of a Census version.

    Args:
        date:
            The date of retraction.
        reason:
            The reason for retraction.
        info_url:
            A permalink to more information.
        replaced_by:
            The census version that replaces this one.
    """

    date: str
    reason: str | None
    info_url: str | None
    replaced_by: str | None


ReleaseFlag = Literal["lts", "retracted"]
ReleaseFlags = dict[ReleaseFlag, bool]


class CensusVersionDescription(TypedDict):
    """A description of a Census version.

    Args:
        release_date:
            The date of the release (deprecated).
        release_build:
            Date of build.
        soma:
            SOMA objects locator.
        h5ads:
            Source H5ADs locator.
        flags:
            Flags for the release.
        retraction:
            If retracted, details of the retraction.
    """

    release_date: str | None
    release_build: str
    soma: CensusLocator
    h5ads: CensusLocator
    flags: NotRequired[ReleaseFlags]
    retraction: NotRequired[CensusVersionRetraction]


CensusDirectory = dict[CensusVersionName, CensusVersionName | CensusVersionDescription]

"""
A provider identifies a storage medium for the Census, which can either be a cloud provider or a local file.
A value of "unknown" can be specified if the provider isn't specified - the API will try to determine
the correct configuration based on the URI.
"""
Provider = Literal["S3", "file", "unknown"]


CensusMirrorName = str  # name of the mirror


class CensusMirror(TypedDict):
    """A mirror for a Census resource.

    A mirror identifies a location that can host the census artifacts. A dict of available mirrors exists in the
    ``mirrors.json`` file, and looks like this:

    .. highlight:: json
    .. code-block:: json

        {
            "default": "default-mirror",
            "default-mirror": {
                "provider": "S3",
                "base_uri": "s3://a-public-bucket/",
                "region": "us-west-2"
            }
        }

    Args:
        provider:
            Provider of the mirror.
        base_uri:
            Base URI for the mirror location, e.g. s3://cellxgene-data-public/.
        region:
            Region of the bucket or resource.
    """

    provider: Provider
    base_uri: str
    region: str | None
    embeddings_base_uri: str


CensusMirrors = dict[CensusMirrorName, CensusMirrorName | CensusMirror]


class ResolvedCensusLocator(TypedDict):
    """A resolved locator for a Census resource.

    A `ResolvedCensusLocator` represent an absolute location of a Census resource, including the provider info. It is
    obtained by resolving a relative location against a specified mirror.

    Args:
        uri:
            Resource URI (absolute).
        region:
            If an S3 URI, has optional region.
        provider:
            Provider.
    """

    uri: str
    region: str | None
    provider: str


# URL for the default top-level directory of all public data
CELL_CENSUS_RELEASE_DIRECTORY_URL = "https://census.cellxgene.cziscience.com/cellxgene-census/v1/release.json"
CELL_CENSUS_MIRRORS_DIRECTORY_URL = "https://census.cellxgene.cziscience.com/cellxgene-census/v1/mirrors.json"


def get_census_version_description(census_version: str) -> CensusVersionDescription:
    """Get release description for given Census version, from the Census release directory.

    Args:
        census_version:
            The census version name.

    Returns:
        ``CensusVersionDescription`` - a dictionary containing a description of the release.

    Raises:
        ValueError: if unknown ``census_version`` value.

    Lifecycle:
        maturing

    See Also:
        :func:`get_census_version_directory`: returns the entire directory as a dict.

    Examples:
        >>> cellxgene_census.get_census_version_description("latest")
        {'release_date': None,
        'release_build': '2022-12-01',
        'soma': {'uri': 's3://cellxgene-data-public/cell-census/2022-12-01/soma/',
        's3_region': 'us-west-2'},
        'h5ads': {'uri': 's3://cellxgene-data-public/cell-census/2022-12-01/h5ads/',
        's3_region': 'us-west-2'}}
    """
    census_directory = get_census_version_directory()
    description = census_directory.get(census_version, None)
    if description is None:
        raise ValueError(f"Unable to locate Census version: {census_version}.")
    return description


def get_census_version_directory(
    *, lts: bool | None = None, retracted: bool | None = False
) -> dict[CensusVersionName, CensusVersionDescription]:
    """Get the directory of Census versions currently available, optionally filtering by specified
    flags. If a filtering flag is not specified, Census versions will not be filtered by that flag.
    Defaults to including both "long-term stable" (LTS) and weekly Census versions, and excluding
    retracted versions.

    Args:
        lts:
            A filtering flag to either include or exclude long-term stable releases in the result.
            If ``None``, no filtering is performed based on this flag. Defaults to ``None``, which
            includes both LTS and non-LTS (weekly) versions.
        retracted:
            A filtering flag to either include or exclude retracted releases in the result. If
            ``None``, no filtering is performed based on this flag. Defaults to ``False``, which
            excludes retracted releases in the result.

    Returns:
        A dictionary that contains Census version names and their corresponding descriptions. Census
        versions are always named by their release date (``YYYY-MM-DD``) but may also have aliases.
        If an alias is specified, the Census version will appear multiple times in the dictionary,
        once under it's release date name, and again for each alias. Aliases may be: ``"stable"``,
        ``"latest"``, or ``"V#"``. The ``"stable"`` alias is used for the most recent LTS release,
        the ``"latest"`` alias is used for the most recent weekly release, and the ``"V#"`` aliases
        are used to identify LTS releases by a sequentially incrementing version number.

    Lifecycle:
        maturing

    See Also:
        :func:`get_census_version_description`: get description by ``census_version``.

    Examples:
        Get all LTS and weekly versions, but exclude retracted LTS versions:

        >>> cellxgene_census.get_census_version_directory()
            {
                'stable': {
                    'release_date': None,
                    'release_build': '2022-11-29',
                    'soma': {'uri': 's3://cellxgene-data-public/cell-census/2022-11-29/soma/',
                             's3_region': 'us-west-2'},
                    'h5ads': {'uri': 's3://cellxgene-data-public/cell-census/2022-11-29/h5ads/',
                              's3_region': 'us-west-2'},
                    'flags': {'lts': True, 'retracted': False}
                },
                'latest': {
                    'release_date': None,
                    'release_build': '2022-12-01',
                    'soma': {'uri': 's3://cellxgene-data-public/cell-census/2022-12-01/soma/',
                             's3_region': 'us-west-2'},
                    'h5ads': {'uri': 's3://cellxgene-data-public/cell-census/2022-12-01/h5ads/',
                              's3_region': 'us-west-2'},
                    'flags': {'lts': True, 'retracted': False}
                },
                'V2': {
                    'release_date': None,
                    'release_build': '2022-11-29',
                    'soma': {'uri': 's3://cellxgene-data-public/cell-census/2022-11-29/soma/',
                             's3_region': 'us-west-2'},
                    'h5ads': {'uri': 's3://cellxgene-data-public/cell-census/2022-11-29/h5ads/',
                              's3_region': 'us-west-2'},
                    'flags': {'lts': True, 'retracted': False}
                },
                '2022-12-01': {
                    'release_date': None,
                    'release_build': '2022-12-01',
                    'soma': {'uri': 's3://cellxgene-data-public/cell-census/2022-12-01/soma/',
                             's3_region': 'us-west-2'},
                    'h5ads': {'uri': 's3://cellxgene-data-public/cell-census/2022-12-01/h5ads/',
                              's3_region': 'us-west-2'},
                    'flags': {'lts': False, 'retracted': False}
                },
                '2022-11-29': {
                    'release_date': None,
                    'release_build': '2022-11-29',
                    'soma': {'uri': 's3://cellxgene-data-public/cell-census/2022-11-29/soma/',
                             's3_region': 'us-west-2'},
                    'h5ads': {'uri': 's3://cellxgene-data-public/cell-census/2022-11-29/h5ads/',
                              's3_region': 'us-west-2'},
                    'flags': {'lts': True, 'retracted': False}
                }
            }

        Get only LTS versions that are not retracted:

        >>> cellxgene_census.get_census_version_directory(lts=True)
            {
                'stable': {
                    'release_date': None,
                    'release_build': '2022-11-29',
                    'soma': {'uri': 's3://cellxgene-data-public/cell-census/2022-11-29/soma/',
                             's3_region': 'us-west-2'},
                    'h5ads': {'uri': 's3://cellxgene-data-public/cell-census/2022-11-29/h5ads/',
                              's3_region': 'us-west-2'},
                    'flags': {'lts': True, 'retracted': False}
                },
                'V2': {
                    'release_date': None,
                    'release_build': '2022-11-29',
                    'soma': {'uri': 's3://cellxgene-data-public/cell-census/2022-11-29/soma/',
                             's3_region': 'us-west-2'},
                    'h5ads': {'uri': 's3://cellxgene-data-public/cell-census/2022-11-29/h5ads/',
                              's3_region': 'us-west-2'},
                    'flags': {'lts': True, 'retracted': False}
                },
                '2022-11-29': {
                    'release_date': None,
                    'release_build': '2022-11-29',
                    'soma': {'uri': 's3://cellxgene-data-public/cell-census/2022-11-29/soma/',
                             's3_region': 'us-west-2'},
                    'h5ads': {'uri': 's3://cellxgene-data-public/cell-census/2022-11-29/h5ads/',
                              's3_region': 'us-west-2'},
                    'flags': {'lts': True, 'retracted': False}
                }
            }

        Get only retracted releases:

        >>> cellxgene_census.get_census_version_directory(retracted=True)
            {
                'V1': {
                    'release_date': None,
                    'release_build': '2022-10-15',
                    'soma': {'uri': 's3://cellxgene-data-public/cell-census/2022-10-15/soma/',
                             's3_region': 'us-west-2'},
                    'h5ads': {'uri': 's3://cellxgene-data-public/cell-census/2022-10-15/h5ads/',
                              's3_region': 'us-west-2'},
                    'flags': {'lts': True, 'retracted': True},
                    'retraction': {
                        'date': '2022-10-30',
                        'reason': 'mistakes happen',
                        'info_url': 'http://cellxgene.com/census/errata/v1',
                        'replaced_by': 'V2'
                    },
                },
                '2022-10-15': {
                    'release_date': None,
                    'release_build': '2022-10-15',
                    'soma': {'uri': 's3://cellxgene-data-public/cell-census/2022-10-15/soma/',
                             's3_region': 'us-west-2'},
                    'h5ads': {'uri': 's3://cellxgene-data-public/cell-census/2022-10-15/h5ads/',
                              's3_region': 'us-west-2'},
                    'flags': {'lts': True, 'retracted': True},
                    'retraction': {
                        'date': '2022-10-30',
                        'reason': 'mistakes happen',
                        'info_url': 'http://cellxgene.com/census/errata/v1',
                        'replaced_by': 'V2'
                    }
                }
            }
    """
    response = requests.get(CELL_CENSUS_RELEASE_DIRECTORY_URL, headers={"User-Agent": _user_agent()})
    response.raise_for_status()

    directory: dict[str, str | dict[str, Any]] = response.json()
    directory_out: CensusDirectory = {}
    aliases: set[CensusVersionName] = set()

    # Resolve all aliases for easier use
    for census_version_name in list(directory.keys()):
        # Strings are aliases for other census_version_name
        directory_value = directory[census_version_name]
        alias = None
        while isinstance(directory_value, str):
            alias = directory_value
            # resolve aliases
            if alias not in directory:
                # oops, dangling pointer -- drop original census_version_name
                directory.pop(census_version_name)
                break

            directory_value = directory[alias]

        if alias:
            aliases.add(census_version_name)

        # exclude aliases
        if not isinstance(directory_value, dict):
            continue

        # Filter fields
        directory_value = {
            k: directory_value[k] for k in CensusVersionDescription.__annotations__ if k in directory_value
        }

        # filter by release flags
        census_version_description = cast(CensusVersionDescription, directory_value)
        release_flags = cast(ReleaseFlags, {"lts": lts, "retracted": retracted})
        admitted = all(
            census_version_description.get("flags", {}).get(flag_name, False) == release_flags[flag_name]
            for flag_name, flag_value in release_flags.items()
            if flag_value is not None
        )
        if not admitted:
            continue

        directory_out[census_version_name] = census_version_description.copy()

    # Cast is safe, as we have removed all aliases
    unordered_directory = cast(dict[CensusVersionName, CensusVersionDescription], directory_out)

    # Sort by aliases and release date, descending
    aliased_releases = [(k, v) for k, v in unordered_directory.items() if k in aliases]
    concrete_releases = [(k, v) for k, v in unordered_directory.items() if k not in aliases]
    ordered_directory = OrderedDict()
    # Note: reverse sorting of aliases serendipitously orders the names we happen to use in a desirable manner:
    # "stable", "latest", "V#"). This will require a more explicit ordering if we change alias naming conventions.
    for k, v in sorted(aliased_releases, key=lambda k: k[0], reverse=True) + sorted(
        concrete_releases, key=lambda k: k[0], reverse=True
    ):
        ordered_directory[k] = v

    return ordered_directory


def get_census_mirror_directory() -> dict[CensusMirrorName, CensusMirror]:
    """Get the directory of Census mirrors currently available.

    Returns:
        A dictionary that contains mirror names and their corresponding info,
        like the provider and the region.

    Lifecycle:
        maturing
    """
    mirrors = _get_census_mirrors()
    del mirrors["default"]
    return cast(dict[CensusMirrorName, CensusMirror], mirrors)


def _get_census_mirrors() -> CensusMirrors:
    response = requests.get(CELL_CENSUS_MIRRORS_DIRECTORY_URL, headers={"User-Agent": _user_agent()})
    response.raise_for_status()
    return cast(CensusMirrors, response.json())



# Section: notebooks-analysis_demo-comp_bio_data_integration_scvi

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



# Section: notebooks-api_demo-census_access_maintained_embeddings

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




# Section: cellxgene_census-src-cellxgene_census-experimental-util-_eager_iter

import logging
import threading
from collections import deque
from collections.abc import Iterator
from concurrent import futures
from concurrent.futures import Future
from typing import TypeVar

util_logger = logging.getLogger("cellxgene_census.experimental.util")

_T = TypeVar("_T")


class _EagerIterator(Iterator[_T]):
    def __init__(
        self,
        iterator: Iterator[_T],
        pool: futures.Executor | None = None,
    ):
        super().__init__()
        self.iterator = iterator
        self._pool = pool or futures.ThreadPoolExecutor()
        self._own_pool = pool is None
        self._future: Future[_T] | None = None
        self._begin_next()

    def _begin_next(self) -> None:
        self._future = self._pool.submit(self.iterator.__next__)
        util_logger.debug("Fetching next iterator element, eagerly")

    def __next__(self) -> _T:
        try:
            assert self._future
            res = self._future.result()
            self._begin_next()
            return res
        except StopIteration:
            self._cleanup()
            raise

    def _cleanup(self) -> None:
        util_logger.debug("Cleaning up eager iterator")
        if self._own_pool:
            self._pool.shutdown()

    def __del__(self) -> None:
        # Ensure the threadpool is cleaned up in the case where the
        # iterator is not exhausted. For more information on __del__:
        # https://docs.python.org/3/reference/datamodel.html#object.__del__
        self._cleanup()
        super_del = getattr(super(), "__del__", lambda: None)
        super_del()


class _EagerBufferedIterator(Iterator[_T]):
    def __init__(
        self,
        iterator: Iterator[_T],
        max_pending: int = 1,
        pool: futures.Executor | None = None,
    ):
        super().__init__()
        self.iterator = iterator
        self.max_pending = max_pending
        self._pool = pool or futures.ThreadPoolExecutor()
        self._own_pool = pool is None
        self._pending_results: deque[futures.Future[_T]] = deque()
        self._lock = threading.Lock()
        self._begin_next()

    def __next__(self) -> _T:
        try:
            res = self._pending_results[0].result()
            self._pending_results.popleft()
            self._begin_next()
            return res
        except StopIteration:
            self._cleanup()
            raise

    def _begin_next(self) -> None:
        def _fut_done(fut: futures.Future[_T]) -> None:
            util_logger.debug("Finished fetching next iterator element, eagerly")
            if fut.exception() is None:
                self._begin_next()

        with self._lock:
            not_running = len(self._pending_results) == 0 or self._pending_results[-1].done()
            if len(self._pending_results) < self.max_pending and not_running:
                _future = self._pool.submit(self.iterator.__next__)
                util_logger.debug("Fetching next iterator element, eagerly")
                _future.add_done_callback(_fut_done)
                self._pending_results.append(_future)
            assert len(self._pending_results) <= self.max_pending

    def _cleanup(self) -> None:
        util_logger.debug("Cleaning up eager iterator")
        if self._own_pool:
            self._pool.shutdown()

    def __del__(self) -> None:
        # Ensure the threadpool is cleaned up in the case where the
        # iterator is not exhausted. For more information on __del__:
        # https://docs.python.org/3/reference/datamodel.html#object.__del__
        self._cleanup()
        super_del = getattr(super(), "__del__", lambda: None)
        super_del()



# Section: notebooks-api_demo-census_duplicated_cells

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



