import numpy as np 
import scipy as sp
import matplotlib.pyplot as plt 
print('Importing Hockney')
import Hockney
M = 4
gridSize = 2**M
L = 1.0
x_grid = np.linspace(0, L, gridSize+1)
dx = abs(x_grid[0] - x_grid[1])

RHSArray = np.zeros((gridSize+1, gridSize+1))
iter = 0
for indexY, Y in np.ndenumerate(x_grid):
    for indexX, X in np.ndenumerate(x_grid):
        RHSArray[indexX, indexY] = iter
        iter = iter+1
print('RHSArray[15, 16] = ', RHSArray[15, 16])
print('RHSArray[16, 15] = ', RHSArray[16, 15])
print('RHSArray[1, 0] = ', RHSArray[1, 0])
print('RHSArray[0, 1] = ', RHSArray[0, 1])
Hockney.convertToRMDA(dx, int(M), RHSArray)
