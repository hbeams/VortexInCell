import numpy as np
import array as arr
import matplotlib.pyplot as plt
from scipy import signal
import Hockney
import copy

class PoissonSolver:
    def __init__(self):
        self.data = []

    def __init__(self, a_dx, a_M, Center = None):
        self._m_dx = a_dx;
        self._m_M = a_M;

    def solve(self, RHSArray, fullGrid = True):
        if (not fullGrid):
            inputArray = copy.deepcopy(RHSArray)
        else:
            inputArray = copy.deepcopy(RHSArray);
        Hockney.convolve(self._m_dx, self._m_M, inputArray)
        return inputArray

M = 3
gridSize = 2**M
L = 1
xCen = L/2
yCen = xCen
r0 = L/8
x_grid = np.linspace(0, L, gridSize+1)
print(x_grid.shape)
dx=np.abs(x_grid[0] - x_grid[1])
print(dx)

RHSArray = np.zeros((gridSize+1, gridSize+1))
rhsRead = np.loadtxt("rhs.txt", delimiter=" ")
for line in np.arange(0, len(rhsRead)):
        RHSArray[int(float(rhsRead[line][0])/dx), int(float(rhsRead[line][1])/dx)] = float(rhsRead[line][2])
        
#plt.imshow(RHSArray)
#plt.colorbar(orientation='vertical')
#plt.show()

PS = PoissonSolver(dx, M)
soln = PS.solve(RHSArray)

#plt.imshow(soln)
#plt.colorbar(orientation='vertical')
#plt.show()

SolnArray = np.zeros((gridSize+1, gridSize+1))
solnRead = np.loadtxt("soln.txt", delimiter=" ")
for line in np.arange(0, len(solnRead)):
        SolnArray[int(float(solnRead[line][0])/dx), int(float(solnRead[line][1])/dx)] = float(solnRead[line][2])

#plt.imshow(SolnArray)
#plt.colorbar(orientation='vertical')
#plt.show()
    
