from torch.nn import Module, Linear, Sigmoid


class MLP(Module):

    def __init__(self, n_inputs, n_outputs):
        self.l1 = Linear(n_inputs, n_outputs)
        self.sigoid = Sigmoid()

    def forward(self, x):
        return self.sigmoid(self.l1(x))
