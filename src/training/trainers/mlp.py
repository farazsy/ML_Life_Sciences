import numpy as np
import torch
import torch.nn as nn
from src.training.models.base_model import BaseModel
from torch.utils.data import DataLoader, TensorDataset


class MLP(BaseModel):
    def __init__(
        self,
        hidden_dims=(256, 128),
        epochs=100,
        batch_size=32,
        learning_rate=0.001,
        weight_decay=0.0001,
        device=None,
    ):

        self.hidden_dims = hidden_dims
        self.epochs = epochs
        self.batch_size = batch_size
        self.learning_rate = learning_rate
        self.weight_decay = weight_decay

        self.device = device or ("cuda" if torch.cuda.is_available() else "cpu")

        self.model = None
        self.loss_history = []

    def build_network(self, input_dim, output_dim):

        layers = []
        previous_dim = input_dim

        for hidden_dim in self.hidden_dims:
            layers.extend([nn.Linear(previous_dim, hidden_dim), nn.ReLU()])

            previous_dim = hidden_dim

        # Output layer: no activation for regression
        layers.append(nn.Linear(previous_dim, output_dim))

        return nn.Sequential(*layers)

    def fit(self, X, y):

        X = torch.tensor(np.asarray(X), dtype=torch.float32)

        y = torch.tensor(np.asarray(y), dtype=torch.float32)

        if y.ndim == 1:
            y = y.unsqueeze(1)

        input_dim = X.shape[1]
        output_dim = y.shape[1]

        self.model = self.build_network(input_dim, output_dim).to(self.device)

        dataset = TensorDataset(X, y)

        loader = DataLoader(dataset, batch_size=self.batch_size, shuffle=True)

        criterion = nn.MSELoss()

        optimizer = torch.optim.Adam(
            self.model.parameters(),
            lr=self.learning_rate,
            weight_decay=self.weight_decay,
        )

        self.loss_history = []

        for epoch in range(self.epochs):
            self.model.train()
            epoch_loss = 0.0

            for X_batch, y_batch in loader:
                X_batch = X_batch.to(self.device)
                y_batch = y_batch.to(self.device)

                optimizer.zero_grad()

                predictions = self.model(X_batch)

                loss = criterion(predictions, y_batch)

                loss.backward()
                optimizer.step()

                epoch_loss += loss.item() * len(X_batch)

            epoch_loss /= len(dataset)

            self.loss_history.append(epoch_loss)

        return self

    def predict(self, X):

        if self.model is None:
            raise RuntimeError("Model has not been trained yet.")

        X = torch.tensor(np.asarray(X), dtype=torch.float32).to(self.device)

        self.model.eval()

        with torch.no_grad():
            predictions = self.model(X)

        return predictions.cpu().numpy()
