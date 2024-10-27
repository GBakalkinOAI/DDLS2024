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

