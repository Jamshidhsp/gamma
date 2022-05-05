print('hi, this is python')
import torch
from kymatio import Scattering2D

x = torch.randn(3, 32, 32)
scattering = Scattering2D(J=1, shape=(32,32))

y = scattering(x)
print('y.shape', y.shape)

phi, psi = scattering.load_filters()

print('phi')
# print(phi)