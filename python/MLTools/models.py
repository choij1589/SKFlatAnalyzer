import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.nn import Sequential, Linear, ReLU, LeakyReLU, Dropout, BatchNorm1d
from torch_geometric.nn import global_mean_pool, knn_graph
from torch_geometric.nn import GraphNorm, LayerNorm
from torch_geometric.nn import MessagePassing
from torch_geometric.utils import dropout_edge

class SNN(nn.Module):
    def __init__(self, nFeatures, nClasses, num_nodes, dropout_p):
        super(SNN, self).__init__()
        # Lecun init is default for pytorch
        self.dropout_p = dropout_p
        self.bn = nn.BatchNorm1d(nFeatures)
        self.dense1 = nn.Linear(nFeatures, num_nodes, bias=True)
        self.dense2 = nn.Linear(num_nodes, num_nodes, bias=True)
        self.dense3 = nn.Linear(num_nodes, num_nodes, bias=True)
        self.dense4 = nn.Linear(num_nodes, num_nodes, bias=True)
        self.dense5 = nn.Linear(num_nodes, num_nodes, bias=True)
        self.dense6 = nn.Linear(num_nodes, nClasses, bias=True)

    def forward(self, x):
        x = F.selu(self.dense1(x))
        x = F.alpha_dropout(x, p=self.dropout_p, training=self.training)
        x = F.selu(self.dense2(x))
        x = F.alpha_dropout(x, p=self.dropout_p, training=self.training)
        x = F.selu(self.dense3(x))
        x = F.alpha_dropout(x, p=self.dropout_p, training=self.training)
        x = F.selu(self.dense4(x))
        x = F.alpha_dropout(x, p=self.dropout_p, training=self.training)
        x = F.selu(self.dense5(x))
        x = F.alpha_dropout(x, p=self.dropout_p, training=self.training)
        out = F.softmax(self.dense6(x), dim=1)

        return out

# NOTE: Using SELU activation layer in mlp block make training unstable
class EdgeConv(MessagePassing):
    def __init__(self, in_channels, out_channels, dropout_p):
        super().__init__(aggr="mean")
        self.mlp = Sequential(
                Linear(2*in_channels, out_channels), LeakyReLU(), BatchNorm1d(out_channels), Dropout(dropout_p),
                Linear(out_channels, out_channels), LeakyReLU(), BatchNorm1d(out_channels), Dropout(dropout_p),
                Linear(out_channels, out_channels), LeakyReLU(), BatchNorm1d(out_channels), Dropout(dropout_p)
                )

    def forward(self, x, edge_index, batch=None):
        return self.propagate(edge_index, x=x, batch=batch)

    def message(self, x_i, x_j):
        tmp = torch.cat([x_i, x_j - x_i], dim=1)
        return self.mlp(tmp)


class DynamicEdgeConv(EdgeConv):
    def __init__(self, in_channels, out_channels, dropout_p, training, k=4):
        super().__init__(in_channels, out_channels, dropout_p=dropout_p)
        self.shortcut = Sequential(Linear(in_channels, out_channels), BatchNorm1d(out_channels), Dropout(dropout_p))
        #self.layer_norm = LayerNorm(out_channels)
        self.dropout_p = dropout_p
        self.training = training
        self.k = k

    def forward(self, x, edge_index=None, batch=None):
        if edge_index is None:
            edge_index = knn_graph(x, self.k, batch, loop=False, flow=self.flow)
        edge_index, _ = dropout_edge(edge_index, p=self.dropout_p, training=self.training)
        out = super().forward(x, edge_index, batch=batch)
        out += self.shortcut(x)
        #out = self.layer_norm(out)
        return out

class ParticleNet(torch.nn.Module):
    def __init__(self, num_node_features, num_graph_features, num_classes, num_hidden, dropout_p):
        super(ParticleNet, self).__init__()
        self.gn0 = GraphNorm(num_node_features)
        self.conv1 = DynamicEdgeConv(num_node_features, num_hidden, dropout_p, training=self.training, k=4)
        self.conv2 = DynamicEdgeConv(num_hidden, num_hidden, dropout_p, training=self.training, k=4)
        self.conv3 = DynamicEdgeConv(num_hidden, num_hidden, dropout_p, training=self.training, k=4)
        #self.linear = Linear(3*num_hidden, num_hidden)

        self.bn0 = BatchNorm1d(num_hidden*3+num_graph_features)
        self.dense1 = Linear(num_hidden*3+num_graph_features, num_hidden)
        self.bn1 = BatchNorm1d(num_hidden)
        self.dense2 = Linear(num_hidden, num_hidden)
        self.bn2 = BatchNorm1d(num_hidden)
        self.output = Linear(num_hidden, num_classes)
        self.dropout_p = dropout_p

    def forward(self, x, edge_index, graph_input, batch=None):
        # Convolution layers
        x = self.gn0(x, batch=batch)
        conv1 = self.conv1(x, edge_index, batch=batch)
        conv2 = self.conv2(conv1, batch=batch)
        conv3 = self.conv3(conv2, batch=batch)
        x = torch.cat([conv1, conv2, conv3], dim=1)

        # readout layers
        x = global_mean_pool(x, batch=batch)
        x = torch.cat([x, graph_input], dim=1)
        x = self.bn0(x)

        # dense layers
        x = F.leaky_relu(self.dense1(x))
        x = self.bn1(x)
        x = F.dropout(x, p=self.dropout_p, training=self.training)
        x = F.leaky_relu(self.dense2(x))
        x = self.bn2(x)
        x = F.dropout(x, p=self.dropout_p, training=self.training)
        x = self.output(x)

        return F.softmax(x, dim=1)

class ParticleNetV2(torch.nn.Module):
    def __init__(self, num_node_features, num_graph_features, num_classes, num_hidden, dropout_p):
        super(ParticleNetV2, self).__init__()
        self.gn0 = GraphNorm(num_node_features)
        self.gn1 = GraphNorm(num_hidden)
        self.gn2 = GraphNorm(num_hidden)
        self.gn3 = GraphNorm(num_hidden)
        self.bn0 = BatchNorm1d(num_hidden+num_graph_features)
        self.bn1 = BatchNorm1d(num_hidden)
        self.bn2 = BatchNorm1d(num_hidden)
        self.conv1 = DynamicEdgeConv(num_node_features, num_hidden, dropout_p, training=self.training, k=4)
        self.conv2 = DynamicEdgeConv(num_hidden, num_hidden, dropout_p, training=self.training, k=4)
        self.conv3 = DynamicEdgeConv(num_hidden, num_hidden, dropout_p, training=self.training, k=4)
        self.dense1 = Linear(num_hidden+num_graph_features, num_hidden)
        self.dense2 = Linear(num_hidden, num_hidden)
        self.output = Linear(num_hidden, num_classes)
        self.dropout_p = dropout_p

    def forward(self, x, edge_index, graph_input, batch=None):
        # Convolution layers
        x = self.gn0(x, batch=batch)
        #y = self.bn0(graph_input)
        conv1 = self.conv1(x, edge_index, batch=batch)
        conv1 = self.gn1(conv1, batch=batch)
        conv2 = self.conv2(conv1, batch=batch)
        conv2 = self.gn2(conv2, batch=batch)
        conv3 = self.conv3(conv2, batch=batch)
        conv3 = self.gn3(conv3, batch=batch)
        x = conv1 + conv2 + conv3

        # readout layers
        x = global_mean_pool(x, batch=batch)
        x = torch.cat([x, graph_input], dim=1)
        x = self.bn0(x)

        # dense layers
        x = F.relu(self.dense1(x))
        x = self.bn1(x)
        x = F.dropout1d(x, p=self.dropout_p, training=self.training)
        x = F.relu(self.dense2(x))
        x = self.bn2(x)
        x = F.dropout1d(x, p=self.dropout_p, training=self.training)
        x = self.output(x)

        return F.softmax(x, dim=1)
