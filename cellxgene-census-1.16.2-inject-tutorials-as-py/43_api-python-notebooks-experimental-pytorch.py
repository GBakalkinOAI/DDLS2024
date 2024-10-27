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




