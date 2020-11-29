from torch.nn import Module, Linear, Softmax 


class MLP(Module):

    def __init__(self, n_inputs, n_outputs):
        super(MLP, self).__init__()
        self.l1 = Linear(n_inputs, n_outputs)
        self.softmax = Softmax(dim=1)

    def forward(self, x):
        return self.softmax(self.l1(x))
