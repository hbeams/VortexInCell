import numpy as np
import array as arr
from scipy import signal


class PoissonSolver:
    # Default Constructor
    def __init__(self):
        self.data = []
        
    # 2D Poisson Greens function
    def kernel(self, x):
        if x == 0:
            return 0
        else: 
            return -1/(2*np.pi)*np.log(x) 
        
    # Assumes grid is a 1D grid and the domain is square
    def __init__(self, grid, Center = None):
        self.grid = grid
        if Center == None:
            xCen = grid[int(len(grid)/2)]
            yCen = xCen
        else:
            xCen = Center
            yCen = Center
            
        gridSize = int(len(self.grid))
        # Computes Greens function kernel for 2d Poisson on domain 
        KernelArray = np.zeros((gridSize,gridSize))
        for indexX, X in np.ndenumerate(grid):
            for indexY, Y in np.ndenumerate(grid):
                KernelArray[indexX, indexY] = self.kernel(np.sqrt((X - xCen)**2 + (Y - yCen)**2))
        # Extend domain for kernel and pad with zeros
        extendedKernel = np.zeros((2*gridSize, 2*gridSize))
        extendedKernel[0:gridSize, 0:gridSize] = KernelArray

        # Take fft of kernel and store this array
        self.fftKernel = np.fft.fft2(extendedKernel, norm = "ortho")
        
    def solve(self, RHSArray):
        gridSize = int(len(self.grid))
        # Build padded domain and write input to padded array
        extendedRHS = np.zeros((2*gridSize, 2*gridSize))
        extendedRHS[0:gridSize, 0:gridSize] = RHSArray
        # FFT rhs. Use ortho to help with normalization
        fftRHS = np.fft.fft2(extendedRHS, norm = "ortho")
        # Multiply transformed arrays 
        fftmultipliedArray = np.multiply(fftRHS, self.fftKernel)
       
        # Not ortho to get the constants correct
        solnShifted = np.fft.ifft2(fftmultipliedArray)
        # Take real values and shift transform back to domain 
        soln = np.real(solnShifted[int(gridSize/2):int(3/2*gridSize), int(gridSize/2):int(3/2*gridSize)])
        
        return soln

