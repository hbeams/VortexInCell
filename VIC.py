#!/usr/bin/env python
# coding: utf-8

# In[267]:


import numpy as np
import array as arr
from scipy import signal
from matplotlib import animation
import sys
sys.path.insert(0, "./../../../vortexincell/src/wrappedCode/")
import PoissonSolver as ps
import matplotlib.pyplot as plt
import matplotlib as mpl
import random
import math
np.set_printoptions(suppress=True)
import itertools


# In[268]:


Debugging = False


# In[269]:


##Generates a random collection of particles 
##Particles will have random positions with x values in [0,xrange] and y values in [0,yrange]
##Particle weights generated with values in [-1,1]

class points:
    def __init__(self):
        self.data = []
    
    def __init__(self, x, y, weight):
        self.x = x
        self.y = y
        self.weight = weight
    
    def generate(number, rangex, rangey):
        data = np.zeros(((number,3)))
        data[:,0] = rangex*np.random.random_sample(number)  ## x value generated
        data[:,1] = rangey*np.random.random_sample(number)  ## y value generated
        data[:,2] = 2*np.random.random_sample(number)-1     ## weight generated
      
        return(data)


# In[270]:


##Generates a one dimensional grid used in the Poisson Solver
##The grid is formed based off of the minimum and maximum particle locations
##The grid that is formed assumes that the x and y dimensions are the same



class gridform:
    def __init__(self):
        self.grid = []
    
    def __init__(self):
        self.grid = grid
        
    def form(particle, size):
        if max(particle[:,0])-min(particle[:,0])> max(particle[:,1])-min(particle[:,1]):
            rangex = math.ceil(max(particle[:,0])-math.floor(min(particle[:,0])))
        else:
            rangex = math.ceil(max(particle[:,1])-math.floor(min(particle[:,1])))
        
        Cen = rangex/2
        
        grid = np.linspace(0, rangex, size)
        grid = grid[:-1]

        
        return grid
    


# In[271]:


##Takes finite differences of a grid


class differences:
    def __init__(self):
        self.fdiffx = []
        self.bdiffx = []
        self.fdiffy = []
        self.bdiffy = []
        self.cdiffx = []
        self.cdiffy = []
        
    def __init__(self):
        self.fdiffx = fdiffx
        self.bdiffx = bdiffx
        self.fdiffy = fdiffy
        self.bdiffy = bdiffy
        self.cdiffx = cdiffx
        self.cdiffy = cdiffy
        
    ##Forwad X difference
    def FDX (grid, dx):
        fdiffx= np.zeros((len(grid[0]),(len(grid) )))


        for i in range (0, len(grid)):
            for k in range (0, len(grid[0])-1):
                fdiffx[i,k] = (grid[i,k]-grid[i,k+1])/dx

        return fdiffx
    
    ##Forward Y difference
    def FDY (grid, dy):
        fdiffy= np.zeros((((len(grid))), ((len(grid[0])))))

        for i in range (0, len(grid)-1):
            for k in range (0, len(grid[0])):
                fdiffy[i,k] = (grid[i,k]-grid[i+1,k])/dy
                
        return fdiffy
    
    ##Backward X difference
    def BDX (grid,dx):

        bdiffx= np.zeros((((len(grid))), ((len(grid[0])))))

        for i in range (0, len(grid)):
            for k in range (1, len(grid[0])):
                bdiffx[i,k] = (grid[i,k]-grid[i,k-1])/dx
        return bdiffx
    
    ##Backward Y difference
    def BDY (grid, dy):
                
        bdiffy= np.zeros((((len(grid))), ((len(grid[0])))))

        for i in range (1, len(grid)):
            for k in range (0, len(grid[0])):
                bdiffy[i,k] = (grid[i,k]-grid[i-1,k])/dy
              
        return bdiffy
    
    ##Central X difference
    def CDX (grid, dx):       
        cdiffx= np.zeros((((len(grid))),((len(grid[0])))))

        for i in range (0, len(grid)):
            for k in range (1, len(grid[0])-1):
                cdiffx[i,k] = (grid[i,k+1]-grid[i,k-1])/(2*dx)

        return cdiffx
    
    ##Central Y difference
    def CDY (grid, dy):
                
        cdiffy= np.zeros((((len(grid))),((len(grid[0])))))

        for i in range (1, len(grid)-1):
            for k in range (0, len(grid[0])):
                cdiffy[i,k] = (grid[i-1,k]-grid[i+1,k])/(2*dy)

        return cdiffy


# In[272]:


##2nd order spline VIC
class dep2nd:
    def __init__(self):
        self.index = []
        self.weight= []
        self.grid = []
        self.indices =[]
        self.domainx = []
        self.domainy = []
        self.Uind = []
        self.U = []
        self.indices = []
        self.column = []
        self.field = []
        self.newpoints = []
        self.rangexmax = []
        self.rangexmin = []
        self.rangeymax = []
        self.rangeymin = []
        
    def __init__(self):
        self.index = index
        self.weight = weight
        self.grid = grid
        self.indices = indices
        self.domainx = domainx
        self.domainy = domainy
        self.Uind = Uind
        self.U = U
        self.indices = indices
        self.column = column
        self.field = field
        self.newpoints = newpoints
        self.rangexmax = rangexmax
        self.rangexmin = rangexmin
        self.rangeymax = rangeymax
        self.rangeymin = rangeymin
        
        
    ##particles are deposited onto a grid with a second order B-spline function    
    def interptogrid(particle, number, size):
        dx = np.abs(gridform.form(particle, size)[0] - gridform.form(particle, size)[1])
        dy = dx
        rangexmax = math.ceil(max(particle[:,0]))
        rangexmin = math.floor(min(particle[:,0]))
        rangeymax = math.ceil(max(particle[:,1]))
        rangeymin = math.floor(min(particle[:,1]))
       
        
        
        index = np.zeros(((number),2))
        weight=np.zeros(((number),4))
        
        if max(particle[:,0])-min(particle[:,0])> max(particle[:,1])-min(particle[:,1]):
            rangex = math.ceil(max(particle[:,0]))-math.floor(min(particle[:,0]))
        else:
            rangex = math.ceil(max(particle[:,1]))-math.floor(min(particle[:,1]))
      
        grid = np.zeros(((int(rangex/dx)),(int(rangex/dx))))
     
    
        
        ##defines the bottom left index for a given particle
        i=0
        while i<number:
            index[i,0]= math.floor((((particle[i,0])-rangexmin)/dx))
            index[i,1]= math.floor((((particle[i,1])-rangeymin)/dy))
            i=i+1
        
        
        
        ##implements spline function to calculate the four weights at each particle
        for i in range (0,number):
            weight[i,0]=(1-abs(((((((particle[i,0]))-(index[i,0]*dx)))/dx))))  ##x left
            weight[i,1]=(1-abs((((((particle[i,0]))-(((index[i,0]+1)*dx)))/dx))))    ##x right
            weight[i,2]=(1-abs((((((particle[i,1]))-((index[i,1]*dy)))/dy))))  ##y bottom
            weight[i,3]=(1-abs(((((((particle[i,1]))-((index[i,1]+1)*dy)))/dy))))    ##y top
       
        
        ##interpolates each particle to the its index points around it
        i=0
        while i<number:
            grid[(int((index[i,0])), int((index[i,1])))]=+ grid[(int((index[i,0]))), (int((index[i,1])))] + ((((particle[i,2]))*weight[i,0]*weight[i,2]))
            grid[(int((index[i,0])), (int((index)[i,1]+1)))]=+ grid[(int((index[i,0]))), (int((index[i,1]+1)))] + (((particle[i,2])*weight[i,0]*weight[i,3]))
            grid[(int((index[i,0]))+1), (int((index)[i,1]))]=+ grid[(int((index[i,0])+1)), (int((index[i,1])))] + (((particle[i,2])*weight[i,1]*weight[i,2]))
            grid[(int((index[i,0])+1)), (int((index)[i,1]+1))]=+ grid[(int((index[i,0]+1))), (int((index[i,1]+1)))] + (((particle[i,2])*weight[i,1]*weight[i,3]))
            i=i+1
        
        ##This remaining part can generate the values at the grid points in a form that is plottable
        ##If needing to plot these values at this step, uncomment the print statement that is #print(indices, 'indices')
        
        indices= np.where(grid !=0)      
        indices=np.reshape(indices,(2, (int((np.size(indices))/2))))
        indices=indices.T
        

        column=np.zeros((((int((np.size(indices))/2),1))))

        indices=np.append(indices,column, axis=1 )
    
        i=0
        while i <len(indices):
            indices[i,2]= (grid[int(indices[i,0]), int(indices[i,1])])
            indices[i,0] = (indices[i,0])*(dx)
            indices[i,1] = (indices[i,1])*(dy)
            i=i+1
        #print(indices, 'indices')
    
     
        return grid
    
 
    ##potential of grid are calculated
    def potential(particle, number, size):
        dx = np.abs(gridform.form(particle, size)[0] - gridform.form(particle, size)[1])
        dy = dx
        mySolver = ps.PoissonSolver(gridform.form(particle, size), False)  ##initializes the Poisson solver with the 1D grid
        #print(particle, 'potparticle')
        domainx = []
        domainy = []
        gridshape= mySolver.solve((dep2nd.interptogrid(particle, number, size))).shape
        lenx= gridshape[0]
        leny = gridshape[1]
    
        
        gridx = np.zeros((lenx,leny))  ##two separate grids are needed since each grid point has both
        gridy = np.zeros((lenx,leny))  ##an x and a y value attached to it

        for i in range(0,lenx):
            domainx.append(i)
        for i in range(0,leny):
            domainy.append(i)
        
        Uind = list(itertools.product(domainx, domainy)) #indices of potentials given in array U

        
        U = np.zeros(((lenx*leny),2))

        U=np.append(U,Uind, axis=1 )
        
        
        solve = mySolver.solve((dep2nd.interptogrid(particle, number, size)))  ##solves Poisson's equation 
        
        
        if (Debugging):
            print('Poisson Solve')
            plt.imshow(solve)
            plt.colorbar(orientation = 'vertical')
            plt.show()
        
        
        
        i=0
        while i <(len(U)):
            U[i,0]= differences.CDX(solve  ,dy)[int(U[i,3]), int(U[i,2])]
            U[i,1]= differences.CDY(solve  ,dx)[int(U[i,3]), int(U[i,2])]  ##Note that CDY here does not have a negative
            i = i+1                                                         ##This is because of the orientation of the input grid
        #print(U,'U')
        i = 0
        while i<(len(U)):
            gridx[int(U[i,3]), int(U[i,2])] = U[i,0]  ##places potentials onto the grids 
            gridy[int(U[i,3]), int(U[i,2])] = U[i,1]
       
            i=i+1
#         print('Ux')
#         plt.imshow(gridx)
#         plt.colorbar(orientation = 'vertical')
#         plt.show()
#         print('Uy')
#         plt.imshow(gridy)
#         plt.colorbar(orientation = 'vertical')
#         plt.show()
        #print(gridx)
        return (gridx, gridy)  
    
    def interptopart(particle, number, size):
        dx = np.abs(gridform.form(particle, size)[0] - gridform.form(particle, size)[1])
        dy = dx
        #print(particle, 'particle')
        rangexmax = math.ceil(max(particle[:,0]))
        rangexmin = math.floor(min(particle[:,0]))
        rangeymax = math.ceil(max(particle[:,1]))
        rangeymin = math.floor(min(particle[:,1]))
       
        grid = dep2nd.potential(particle, number, size)  
        
        ##dep2nd.potential returns the gridx and gridy as a single array
        ##the following unpacks the two arrays
        gridy = np.array(grid[len(dep2nd.potential(particle, number, size))//2:])
        gridx = np.array(grid[:len(dep2nd.potential(particle, number, size))//2])
        gridx = np.reshape(gridx,(len(gridx[0]),len(gridx[0])))
        gridy = np.reshape(gridy,(len(gridy[0]), len(gridy[0])))
        
        
        
        index = np.zeros(((number),2))
        weight = np.zeros(((number, 4)))
        newpoints = np.zeros(((number,4)))
        
        for i in range (0, number):
            newpoints[i,0] = particle[i,0]
            newpoints[i,1] = particle[i,1]
            
        i=0
        while i<number:
            index[i,0]= math.floor((((particle[i,0])-math.floor(rangexmin))/dx))
            index[i,1]= math.floor((((particle[i,1])-math.floor(rangeymin))/dy))
            i=i+1
        #print(index)
        
        
        for i in range (0,number):
            weight[i,0]=(1-abs((((((index[i,0]*dx))-((particle[i,0])))/dx))))  ##x left
            weight[i,1]=(1-abs(((((((index[i,0]+1)*dx))-((particle[i,0])))/dx))))    ##x right
            weight[i,2]=(1-abs((((((index[i,1]*dy))-((particle[i,1])))/dy))))  ##y bottom
            weight[i,3]=(1-abs(((((((index[i,1]+1)*dy))-((particle[i,1])))/dy))))    ##y top
        #print(weight)
    
        ##interpolate back to the initial particles 
        i=0
        while i<number:
            newpoints[i,2] = newpoints[i,2] + ((gridx[int(index[i,0]), int(index[i,1])])*weight[i,2]*weight[i,0])
            #print(gridx[int(index[i,0]), int(index[i,1])])
            #print(newpoints)
            newpoints[i,2] = newpoints[i,2] + ((gridx[int(index[i,0]), int(index[i,1])+1])*weight[i,0]*weight[i,3])
            #print(newpoints)
            #print(gridx[int(index[i,0])+1, int(index[i,1])])
            newpoints[i,2] = newpoints[i,2] + ((gridx[int(index[i,0])+1, int(index[i,1])])*weight[i,1]*weight[i,2])
            #print(newpoints)
            newpoints[i,2] = newpoints[i,2] + ((gridx[int(index[i,0])+1, int(index[i,1])+1])*weight[i,1]*weight[i,3])
            #print(newpoints)
            newpoints[i,3] = newpoints[i,3] + ((gridy[int(index[i,0]), int(index[i,1])])*weight[i,2]*weight[i,0])
            #print(newpoints)
            newpoints[i,3] = newpoints[i,3] + ((gridy[int(index[i,0]), int(index[i,1])+1])*weight[i,3]*weight[i,0])
            #print(newpoints)
            newpoints[i,3] = newpoints[i,3] + ((gridy[int(index[i,0])+1, int(index[i,1])])*weight[i,2]*weight[i,1])
            #print(newpoints)
            newpoints[i,3] = newpoints[i,3] + ((gridy[int(index[i,0])+1, int(index[i,1])+1])*weight[i,3]*weight[i,1])
 
            
            
            
            i = i+1
        
        column=np.zeros((((int((len(newpoints))),1))))
        newpoints = np.append(newpoints, column, axis=1)
        for i in range(0, len(newpoints)):
            newpoints[i,4] = particle[i,2]

        return newpoints
        
    ##same code as remaining part of interptogrid.
    ##Can call separately here
    
    def plotvalues (particle, number, size):
        dx = np.abs(gridform.form(particle, size)[0] - gridform.form(particle, size)[1])
        dy = dx
        
        indices= np.where((dep2nd.interptogrid(particle, number, size)) !=0)      
        indices=np.reshape(indices,(2, (int((np.size(indices))/2))))
        indices=indices.T
        

        column=np.zeros((((int((np.size(indices))/2),1))))

        indices=np.append(indices,column, axis=1 )
    
        i=0
        while i <len(indices):
            indices[i,2]= (dep2nd.interptogrid(particle, number, size))[int(indices[i,0]), int(indices[i,1])]
            indices[i,0] = (indices[i,0])*(dx)
            indices[i,1] = (indices[i,1])*(dy)
            i=i+1
        #print(indices, 'indices')
        return indices
    

    
    
    


# In[273]:


##4th order spline VIC
class dep4th:
    def __init__(self):
        self.index = []
        self.weight= []
        self.grid = []
        self.indices =[]
        self.domainx = []
        self.domainy = []
        self.Uind = []
        self.U = []
        self.indices = []
        self.column = []
        self.field = []
        self.newpoints = []
        self.rangexmax = []
        self.rangexmin = []
        self.rangeymax = []
        self.rangeymin = []
        
    def __init__(self):
        self.index = index
        self.weight = weight
        self.grid = grid
        self.indices = indices
        self.domainx = domainx
        self.domainy = domainy
        self.Uind = Uind
        self.U = U
        self.indices = indices
        self.column = column
        self.field = field
        self.newpoints = newpoints
        self.rangexmax = rangexmax
        self.rangexmin = rangexmin
        self.rangeymax = rangeymax
        self.rangeymin = rangeymin
        
        
    def interptogrid(particle, number, size):
        dx = np.abs(gridform.form(particle, size)[0] - gridform.form(particle, size)[1])
        dy = dx
        rangexmax = math.ceil(max(particle[:,0]))
        rangexmin = math.floor(min(particle[:,0]))
        rangeymax = math.ceil(max(particle[:,1]))
        rangeymin = math.floor(min(particle[:,1]))
       
        
        
        index = np.zeros(((number),2))
        weight=np.zeros(((number),8))
        
        if max(particle[:,0])-min(particle[:,0])> max(particle[:,1])-min(particle[:,1]):
            rangex = math.ceil(max(particle[:,0]))-math.floor(min(particle[:,0]))
        else:
            rangex = math.ceil(max(particle[:,1]))-math.floor(min(particle[:,1]))
      
        grid = np.zeros(((int(rangex/dx)),(int(rangex/dx))))
      
      
        i=0
        while i<number:
            index[i,0]= math.floor((((particle[i,0])-rangexmin)/dx))-1
            index[i,1]= math.floor((((particle[i,1])-rangeymin)/dy))-1
            i=i+1
      
        
        for i in range (0,number):
            weight[i,0] = 1 - ((5*((abs((((index[i,0]+1)*dx)-particle[i,0])/dx))**2))/2) + ((3*((abs((((index[i,0]+1)*dx)-particle[i,0])/dx))**3))/2) ##x left
            weight[i,1] = 1 - ((5*((abs((((index[i,1]+2)*dy)-particle[i,1])/dy))**2))/2) + ((3*((abs((((index[i,1]+2)*dy)-particle[i,1])/dy))**3))/2) ##y top
            weight[i,2] = 1 - ((5*((abs((((index[i,0]+2)*dx)-particle[i,0])/dx))**2))/2) + ((3*((abs((((index[i,0]+2)*dx)-particle[i,0])/dx))**3))/2) ##x right
            weight[i,3] = 1 - ((5*((abs((((index[i,1]+1)*dy)-particle[i,1])/dy))**2))/2) + ((3*((abs((((index[i,1]+1)*dy)-particle[i,1])/dy))**3))/2) ##y bottom
            weight[i,4] = (1/2)*((2-abs((((index[i,0]+3)*dx)- particle[i,0])/dx))**2)*(1-abs((((index[i,0]+3)*dx)- particle[i,0])/dx)) ##x right right
            weight[i,5] = (1/2)*((2-abs((((index[i,0])*dx)- particle[i,0])/dx))**2)*(1-abs((((index[i,0])*dx)- particle[i,0])/dx)) ##x left left
            weight[i,6] = (1/2)*((2-abs((((index[i,1]+3)*dy)- particle[i,1])/dy))**2)*(1-abs((((index[i,1]+3)*dy)- particle[i,1])/dy)) ##y top top
            weight[i,7] = (1/2)*((2-abs((((index[i,1])*dy)- particle[i,1])/dy))**2)*(1-abs((((index[i,1])*dy)- particle[i,1])/dy)) ##y bottom bottom                                                                       
      
        i=0
        while i<number:
            grid[(int((index[i,0])), int((index[i,1])))]=+ grid[(int((index[i,0]))), (int((index[i,1])))] + ((((particle[i,2]))*weight[i,5]*weight[i,7]))
            grid[(int((index[i,0])), (int((index)[i,1]+1)))]=+ grid[(int((index[i,0]))), (int((index[i,1]+1)))] + (((particle[i,2])*weight[i,5]*weight[i,3]))
            grid[(int((index[i,0])), (int((index)[i,1]+2)))]=+ grid[(int((index[i,0]))), (int((index[i,1]+2)))] + (((particle[i,2])*weight[i,5]*weight[i,1]))
            grid[(int((index[i,0])), (int((index)[i,1]+3)))]=+ grid[(int((index[i,0]))), (int((index[i,1]+3)))] + (((particle[i,2])*weight[i,5]*weight[i,6]))
            
            grid[(int((index[i,0]))+1), (int((index)[i,1]))]=+ grid[(int((index[i,0])+1)), (int((index[i,1])))] + (((particle[i,2])*weight[i,0]*weight[i,7]))
            grid[(int((index[i,0]))+1, (int((index)[i,1]+1)))]=+ grid[(int((index[i,0])+1)), (int((index[i,1]+1)))] + (((particle[i,2])*weight[i,0]*weight[i,3]))
            grid[(int((index[i,0]))+1, (int((index)[i,1]+2)))]=+ grid[(int((index[i,0])+1)), (int((index[i,1]+2)))] + (((particle[i,2])*weight[i,0]*weight[i,1]))
            grid[(int((index[i,0]))+1, (int((index)[i,1]+3)))]=+ grid[(int((index[i,0])+1)), (int((index[i,1]+3)))] + (((particle[i,2])*weight[i,0]*weight[i,6])) 
            grid[(int((index[i,0]))+2), (int((index)[i,1]))]=+ grid[(int((index[i,0])+2)), (int((index[i,1])))] + (((particle[i,2])*weight[i,2]*weight[i,7]))
            grid[(int((index[i,0]))+2, (int((index)[i,1]+1)))]=+ grid[(int((index[i,0])+2)), (int((index[i,1]+1)))] + (((particle[i,2])*weight[i,2]*weight[i,3]))
            grid[(int((index[i,0]))+2, (int((index)[i,1]+2)))]=+ grid[(int((index[i,0])+2)), (int((index[i,1]+2)))] + (((particle[i,2])*weight[i,2]*weight[i,1]))
            grid[(int((index[i,0]))+2, (int((index)[i,1]+3)))]=+ grid[(int((index[i,0])+2)), (int((index[i,1]+3)))] + (((particle[i,2])*weight[i,2]*weight[i,6])) 
            grid[(int((index[i,0]))+3), (int((index)[i,1]))]=+ grid[(int((index[i,0])+3)), (int((index[i,1])))] + (((particle[i,2])*weight[i,4]*weight[i,7]))
            grid[(int((index[i,0]))+3, (int((index)[i,1]+1)))]=+ grid[(int((index[i,0])+3)), (int((index[i,1]+1)))] + (((particle[i,2])*weight[i,4]*weight[i,3]))
            grid[(int((index[i,0]))+3, (int((index)[i,1]+2)))]=+ grid[(int((index[i,0])+3)), (int((index[i,1]+2)))] + (((particle[i,2])*weight[i,4]*weight[i,1]))
            grid[(int((index[i,0]))+3, (int((index)[i,1]+3)))]=+ grid[(int((index[i,0])+3)), (int((index[i,1]+3)))] + (((particle[i,2])*weight[i,4]*weight[i,6]))                                                                            
            
            i=i+1
        ##This remaining part can generate the values at the grid points in a form that is plottable
        ##If needing to plot these values at this step, uncomment the print statement that is #print(indices, 'indices')
     
        indices= np.where(grid !=0)      
        indices=np.reshape(indices,(2, (int((np.size(indices))/2))))
        indices=indices.T
        

        column=np.zeros((((int((np.size(indices))/2),1))))

        indices=np.append(indices,column, axis=1 )
    
        i=0
        while i <len(indices):
            indices[i,2]= (grid[int(indices[i,0]), int(indices[i,1])])
            indices[i,0] = (indices[i,0])*(dx)
            indices[i,1] = (indices[i,1])*(dy)
            i=i+1
        #print(indices, 'indices')
    
     
        return grid
   

    
    ##particles are deposited onto a grid with a second order B-spline function    
    
    
 
    ##potential of grid are calculated
    def potential(particle, number, size):
        dx = np.abs(gridform.form(particle, size)[0] - gridform.form(particle, size)[1])
        dy = dx
        mySolver = ps.PoissonSolver(gridform.form(particle, size), False)  ##initializes the Poisson solver with the 1D grid
      
        domainx = []
        domainy = []
        gridshape= mySolver.solve((dep4th.interptogrid(particle, number, size))).shape
        lenx= gridshape[0]
        leny = gridshape[1]
    
        
        gridx = np.zeros((lenx,leny))  ##two separate grids are needed since each grid point has both
        gridy = np.zeros((lenx,leny))  ##an x and a y value attached to it

        for i in range(0,lenx):
            domainx.append(i)
        for i in range(0,leny):
            domainy.append(i)
        
        Uind = list(itertools.product(domainx, domainy)) #indices of potentials given in array U

        
        U = np.zeros(((lenx*leny),2))

        U=np.append(U,Uind, axis=1 )
        
        
        solve = mySolver.solve((dep4th.interptogrid(particle, number, size)))  ##solves Poisson's equation 
       
        
        if (Debugging):
            print('Poisson Solve')
            plt.imshow(solve)
            plt.colorbar(orientation = 'vertical')
            plt.show()
        
        
        
        i=0
        while i <(len(U)):
            U[i,0]= differences.CDX(solve  ,dy)[int(U[i,3]), int(U[i,2])]
            U[i,1]= differences.CDY(solve  ,dx)[int(U[i,3]), int(U[i,2])]  ##Note that CDY here does not have a negative
            i = i+1                                                         ##This is because of the orientation of the input grid
        
        
        i = 0
        while i<(len(U)):
            gridx[int(U[i,3]), int(U[i,2])] = U[i,0]  ##places potentials onto the grids 
            gridy[int(U[i,3]), int(U[i,2])] = U[i,1]
       
            i=i+1
#         print('Ux')
#         plt.imshow(gridx)
#         plt.colorbar(orientation = 'vertical')
#         plt.show()
#         print('Uy')
#         plt.imshow(gridy)
#         plt.colorbar(orientation = 'vertical')
#         plt.show()
        #print(gridx)
        return (gridx, gridy)  
    
    def interptopart(particle, number, size):
        dx = np.abs(gridform.form(particle, size)[0] - gridform.form(particle, size)[1])
        dy = dx
        
        rangexmax = math.ceil(max(particle[:,0]))
        rangexmin = math.floor(min(particle[:,0]))
        rangeymax = math.ceil(max(particle[:,1]))
        rangeymin = math.floor(min(particle[:,1]))
       
        grid = dep4th.potential(particle, number, size)

        gridy = np.array(grid[len(dep4th.potential(particle, number, size))//2:])
        gridx = np.array(grid[:len(dep4th.potential(particle, number, size))//2])
        gridx = np.reshape(gridx,(len(gridx[0]),len(gridx[0])))
        gridy = np.reshape(gridy,(len(gridy[0]), len(gridy[0])))
        

        index = np.zeros(((number),2))
        weight = np.zeros(((number, 8)))
        newpoints = np.zeros(((number,4)))
        
        for i in range (0, number):
            newpoints[i,0] = particle[i,0]
            newpoints[i,1] = particle[i,1]

        i=0
        while i<number:
            index[i,0]= math.floor((((particle[i,0])-rangexmin)/dx))-1
            index[i,1]= math.floor((((particle[i,1])-rangeymin)/dy))-1
            i=i+1
      
        
        for i in range (0,number):
            weight[i,0] = 1 - ((5*((abs((particle[i,0]-((index[i,0]+1)*dx))/dx))**2))/2) + ((3*((abs((particle[i,0]-((index[i,0]+1)*dx))/dx))**3))/2) ##x left
            weight[i,1] = 1 - ((5*((abs((particle[i,1]-((index[i,1]+2)*dy))/dy))**2))/2) + ((3*((abs((particle[i,1]-((index[i,1]+2)*dy))/dy))**3))/2) ##y top
            weight[i,2] = 1 - ((5*((abs((particle[i,0]-((index[i,0]+2)*dx))/dx))**2))/2) + ((3*((abs((particle[i,0]-((index[i,0]+2)*dx))/dx))**3))/2) ##x right
            weight[i,3] = 1 - ((5*((abs((particle[i,1]-((index[i,1]+1)*dy))/dy))**2))/2) + ((3*((abs((particle[i,1]-((index[i,1]+1)*dy))/dy))**3))/2) ##y bottom
            weight[i,4] = (1/2)*((2-abs((particle[i,0]-((index[i,0]+3)*dx))/dx))**2)*(1-abs((particle[i,0]-((index[i,0]+3)*dx))/dx)) ##x right right
            weight[i,5] = (1/2)*((2-abs((particle[i,0]-((index[i,0])*dx))/dx))**2)*(1-abs((particle[i,0]-((index[i,0])*dx))/dx)) ##x left left
            weight[i,6] = (1/2)*((2-abs((particle[i,1]-((index[i,1]+3)*dy))/dy))**2)*(1-abs((particle[i,1]-((index[i,1]+3)*dy))/dy)) ##y top top
            weight[i,7] = (1/2)*((2-abs((particle[i,1]-((index[i,1])*dy))/dy))**2)*(1-abs((particle[i,1]-((index[i,1])*dy))/dy)) ##y bottom bottom                                                                       
        
    
        
        i=0
        while i<number:
            newpoints[i,2] =+ newpoints[i,2] + ((gridx[(int((index[i,0]))), (int((index[i,1])))])*weight[i,5]*weight[i,7])
            
            newpoints[i,2] =+ newpoints[i,2] + ((gridx[(int((index[i,0]))), (int((index[i,1]+1)))])*weight[i,5]*weight[i,3])

            newpoints[i,2] =+ newpoints[i,2] + ((gridx[(int((index[i,0]))), (int((index[i,1]+2)))])*weight[i,5]*weight[i,1])

            newpoints[i,2] =+ newpoints[i,2] + ((gridx[(int((index[i,0]))), (int((index[i,1]+3)))])*weight[i,5]*weight[i,6])
       
            newpoints[i,2] =+ newpoints[i,2] + ((gridx[(int((index[i,0]))+1), (int((index[i,1])))])*weight[i,0]*weight[i,7])
  
            newpoints[i,2] = newpoints[i,2] + ((gridx[(int((index[i,0]))+1), (int((index[i,1]+1)))])*weight[i,0]*weight[i,3])

            newpoints[i,2] = newpoints[i,2] + ((gridx[(int((index[i,0]))+1), (int((index[i,1]+2)))])*weight[i,0]*weight[i,1])
   
            newpoints[i,2] = newpoints[i,2] + ((gridx[(int((index[i,0]))+1), (int((index[i,1]+3)))])*weight[i,0]*weight[i,6])
           
            newpoints[i,2] = newpoints[i,2] + ((gridx[(int((index[i,0]))+2), (int((index[i,1])))])*weight[i,2]*weight[i,7])
          
            newpoints[i,2] = newpoints[i,2] + ((gridx[(int((index[i,0]))+2), (int((index[i,1]+1)))])*weight[i,2]*weight[i,3])
           
            newpoints[i,2] = newpoints[i,2] + ((gridx[(int((index[i,0]))+2), (int((index[i,1]+2)))])*weight[i,2]*weight[i,1])
       
            newpoints[i,2] = newpoints[i,2] + ((gridx[(int((index[i,0]))+2), (int((index[i,1]+3)))])*weight[i,2]*weight[i,6])
     
            newpoints[i,2] = newpoints[i,2] + ((gridx[(int((index[i,0]))+3), (int((index[i,1])))])*weight[i,4]*weight[i,7])

            newpoints[i,2] = newpoints[i,2] + ((gridx[(int((index[i,0]))+3), (int((index[i,1]+1)))])*weight[i,4]*weight[i,3])

            newpoints[i,2] = newpoints[i,2] + ((gridx[(int((index[i,0]))+3), (int((index[i,1]+2)))])*weight[i,4]*weight[i,1])

            newpoints[i,2] = newpoints[i,2] + ((gridx[(int((index[i,0]))+3), (int((index[i,1]+3)))])*weight[i,4]*weight[i,6])

            newpoints[i,3] = newpoints[i,3] + ((gridy[(int((index[i,0]))), (int((index[i,1])))])*weight[i,5]*weight[i,7])
            newpoints[i,3] = newpoints[i,3] + ((gridy[(int((index[i,0]))), (int((index[i,1]+1)))])*weight[i,5]*weight[i,3])
            newpoints[i,3] = newpoints[i,3] + ((gridy[(int((index[i,0]))), (int((index[i,1]+2)))])*weight[i,5]*weight[i,1])
            newpoints[i,3] = newpoints[i,3] + ((gridy[(int((index[i,0]))), (int((index[i,1]+3)))])*weight[i,5]*weight[i,6])
            newpoints[i,3] = newpoints[i,3] + ((gridy[(int((index[i,0]))+1), (int((index[i,1])))])*weight[i,0]*weight[i,7])
            newpoints[i,3] = newpoints[i,3] + ((gridy[(int((index[i,0]))+1), (int((index[i,1]+1)))])*weight[i,0]*weight[i,3])
            newpoints[i,3] = newpoints[i,3] + ((gridy[(int((index[i,0]))+1), (int((index[i,1]+2)))])*weight[i,0]*weight[i,1])
            newpoints[i,3] = newpoints[i,3] + ((gridy[(int((index[i,0]))+1), (int((index[i,1]+3)))])*weight[i,0]*weight[i,6])
            newpoints[i,3] = newpoints[i,3] + ((gridy[(int((index[i,0]))+2), (int((index[i,1])))])*weight[i,2]*weight[i,7])
            newpoints[i,3] = newpoints[i,3] + ((gridy[(int((index[i,0]))+2), (int((index[i,1]+1)))])*weight[i,2]*weight[i,3])
            newpoints[i,3] = newpoints[i,3] + ((gridy[(int((index[i,0]))+2), (int((index[i,1]+2)))])*weight[i,2]*weight[i,1])
            newpoints[i,3] = newpoints[i,3] + ((gridy[(int((index[i,0]))+2), (int((index[i,1]+3)))])*weight[i,2]*weight[i,6])
            newpoints[i,3] = newpoints[i,3] + ((gridy[(int((index[i,0]))+3), (int((index[i,1])))])*weight[i,4]*weight[i,7])
            newpoints[i,3] = newpoints[i,3] + ((gridy[(int((index[i,0]))+3), (int((index[i,1]+1)))])*weight[i,4]*weight[i,3])
            newpoints[i,3] = newpoints[i,3] + ((gridy[(int((index[i,0]))+3), (int((index[i,1]+2)))])*weight[i,4]*weight[i,1])
            newpoints[i,3] = newpoints[i,3] + ((gridy[(int((index[i,0]))+3), (int((index[i,1]+3)))])*weight[i,4]*weight[i,6])
            i = i+1
        
        column=np.zeros((((int((len(newpoints))),1))))
        newpoints = np.append(newpoints, column, axis=1)
        for i in range(0, len(newpoints)):
            newpoints[i,4] = particle[i,2]

        return newpoints
    
    
    
    ##same code as remaining part of interptogrid.
    ##Can call separately here
    
    def plotvalues (particle, number, size):
        dx = np.abs(gridform.form(particle, size)[0] - gridform.form(particle, size)[1])
        dy = dx
        
        indices= np.where((dep4th.interptogrid(particle, number, size)) !=0)      
        indices=np.reshape(indices,(2, (int((np.size(indices))/2))))
        indices=indices.T
        

        column=np.zeros((((int((np.size(indices))/2),1))))

        indices=np.append(indices,column, axis=1 )
    
        i=0
        while i <len(indices):
            indices[i,2]= (dep4th.interptogrid(particle, number, size))[int(indices[i,0]), int(indices[i,1])]
            indices[i,0] = (indices[i,0])*(dx)
            indices[i,1] = (indices[i,1])*(dy)
            i=i+1
        #print(indices, 'indices')
        return indices
    


# In[274]:


##Test 1 parameters

# particles = np.array([[0.24,0.26]])
# number = 1
# rangex = 1
# rangey = 1
# size = (2**5)+1     ##size must be in form of (2**N)+1 where N is any integer. This is needed by the PoissonSolver


# In[289]:


#Test 2 parameters

number = 2
rangex = 1
rangey = 1
size = (2**6)+1
particles = np.array([[0.25, 0.5],[0.5, 0.5]])
dx = np.abs(gridform.form(particles, size)[0] - gridform.form(particles, size)[1])
print(dx)
column=np.zeros((((int((len(particles))),1))))
particles = np.append(particles, column, axis=1)
particles[1,2] = (1/(dx**2))
particles[0,2] = 0


# In[290]:


##Test 3 parameters

# number = 2
# rangex = 1
# rangey = 1
# size = (2**5)+1
# particles = np.array([[0.25, 0.5],[0.75, 0.5]])
# dx = np.abs(gridform.form(particles, size)[0] - gridform.form(particles, size)[1])
# print(dx)
# column=np.zeros((((int((len(particles))),1))))
# particles = np.append(particles, column, axis=1)
# particles[1,2] = (1/(dx**2))
# particles[0,2] = (1/(dx**2))


# In[291]:


##Test 4 a parameters
# dx = 0.03125    #need to know dx to begin with to define Np
# Np = 1/(dx/2)
# hp = 1/Np
# particles = []

# rangex = 1
# rangey = 1
# size = (2**5)+1
# for i in range(0, int(Np)):
#     for k in range(0, int(Np)):
#         cond1 = math.sqrt(((i*dx - 0.5)**2 + (k*dx - 0.375)**2))
#         #print(cond1)
#         cond2 = math.sqrt(((i*dx - 0.5)**2 + (k*dx - 0.625)**2))
#         #print(cond2)
#         if cond1 <= 0.12 or  cond2<= 0.12:
#             particles = np.append(particles,[i*dx, k*dx, (hp**2)/(dx**2)])
# particles = np.reshape(particles, (int((np.shape(particles)[0])/3), 3))
# number = int(len(particles))


# In[292]:


# ##Test 4 b
# dx = 0.03125
# Np = 1/(dx/2)
# hp = 1/Np
# rbdry = 0.25
# particles = []

# rangex = 1
# rangey = 1
# size = (2**5)+1

# for i in range(0, int(Np)):
#     for k in range(0, int(Np)):
#         cond1 = math.sqrt(((i*dx - 0.5)**2 + (k*dx - 0.5)**2))
#         #print(cond1)
#         if cond1 <= rbdry:
#             particles = np.append(particles,[i*dx, k*dx, ((rbdry- (((i*dx - 0.5)**2 + (k*dx - 0.5)**2)))**4)])
# #             print(((rbdry- (((i*dx - 0.5)**2 + (k*dx - 0.5)**2)))**7))
# #             print(rbdry-(((i*dx - 0.5)**2 + (k*dx - 0.5)**2)), 'n7')
# particles = np.reshape(particles, (int((np.shape(particles)[0])/3), 3))
# number = int(len(particles))


# In[293]:


plt.scatter(particles[:,0], particles[:,1], particles[:,2])


# In[294]:


##newpoints of the form [x, y, dx/dt, dy/dt, weight]
grid1D = gridform.form(particles, size)
plotdata = dep4th.plotvalues(particles, number, size)
potential = dep4th.potential(particles, number, size)   ##indexed potentials
newpoints = dep4th.interptopart(particles, number, size)
print(newpoints)


if (Debugging):
    print(grid1D, 'grid')
    print(grid1D)
    dx = np.abs(gridform.form(particle, size)[0] - gridform.form(particle, size)[1])
    dy = dx
    print(dx)
    data = dep4th.interptogrid(particles,number,size)
    print(data)
    print(len(data), 'LENGTH')
    step = PoissonSolver(grid1D, Center=None)
    fieldsolve = PoissonSolver.solve(step,data)  
    print(fieldsolve, 'fieldsolve')
    print(plotdata)
    potentialongridX = differences.CDX(fieldsolve,dx) ##x potentials on grid
    print(potentialongridX, 'px')
    potentialongridY = differences.CDY(fieldsolve,dy) ##y potentials on grid
    print(potentialongridY, 'py')
    print(differences.FDX(fieldsolve,dx), 'fdx')
    print(differences.BDX(fieldsolve,dx), 'bdx')
    print(differences.FDY(fieldsolve,dy), 'fdy')
    print(differences.BDY(fieldsolve,dy), 'bdy')
    print(potential, 'p')
    print(plotdata)



# In[295]:



fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
dx = np.abs(gridform.form(particles, size)[0] - gridform.form(particles, size)[1])
# Major ticks every 20, minor ticks every 5
major_ticks = np.arange(0, rangex+1, dx)
minor_ticks = np.arange(0, rangey+1, dx)

ax.set_xticks(major_ticks)
ax.set_xticks(minor_ticks, minor=True)
ax.set_yticks(major_ticks)
ax.set_yticks(minor_ticks, minor=True)

# And a corresponding grid
ax.grid(which='both')

# Or if you want different settings for the grids:
ax.grid(which='minor', alpha=0.2)
ax.grid(which='major', alpha=0.5)
    
    
    
    
    

plt.scatter(plotdata[np.where(plotdata[:,2]>0),0], plotdata[np.where(plotdata[:,2]>0),1], plotdata[np.where(plotdata[:,2]>0),2], 'r')
plt.scatter(plotdata[np.where(plotdata[:,2]<0),0], plotdata[np.where(plotdata[:,2]<0),1], abs(plotdata[np.where(plotdata[:,2]<0),2]), 'b')
plt.scatter(particles[:,0], particles[:,1], 10)
plt.quiver(newpoints[:,0], newpoints[:,1], newpoints[:,2], newpoints[:,3])
# plt.scatter(particles[np.where(particles[:,2]>0),0], particles[np.where(particles[:,2]>0),1], particles[np.where(particles[:,2]>0),2]/(dx*dy),'m')
# plt.scatter(particles[np.where(particles[:,2]<0),0], particles[np.where(particles[:,2]<0),1], abs(particles[np.where(particles[:,2]<0),2]/(dx*dy)), 'c')

# plt.scatter(newpoints[np.where(newpoints[:,2]>0),0], newpoints[np.where(newpoints[:,2]>0),1], newpoints[np.where(newpoints[:,2]>0),2]*200, 'k')
# plt.scatter(newpoints[np.where(newpoints[:,3]>0),0], newpoints[np.where(newpoints[:,3]>0),1], newpoints[np.where(newpoints[:,3]>0),2]*200, 'k')
# plt.scatter(newpoints[np.where(newpoints[:,2]<0),0], newpoints[np.where(newpoints[:,2]<0),1], newpoints[np.where(newpoints[:,2]<0),2]*200, 'k')
# plt.scatter(newpoints[np.where(newpoints[:,3]<0),0], newpoints[np.where(newpoints[:,3]<0),1], newpoints[np.where(newpoints[:,3]<0),2]*200, 'k')


# for i in range(0, len(newpoints)):
#     if newpoints[i,2]!=0 or newpoints[i,3]!=0:
#         plt.arrow(newpoints[i,0], newpoints[i,1], newpoints[i,2]*(50), newpoints[i,3]*(50) ,head_width=(10/50), head_length=10/50, fc='k', ec='k')
# print(plotdata)        
       


# In[296]:


class RHS:
    def __init__ (self, data):
        self.data = []
        self.data = data
    
    def computex(self, number, state):
        temp = self.data
        dx = temp[number,2]  
        return dx
    
    def computey(self, number, state):
        temp = self.data
        dy = temp[number,3]
        return dy


# In[297]:


##RK2 for 4th order spline code


class RK2:
    def __init__ (self, RHS):
        self.RHS = RHS
        
    
    def advance(number, dt, CurrentState, spline):
        finalstate = np.zeros(((len(CurrentState), 3)))
        midpos = np.zeros(((len(CurrentState),3)))
        
        
        for i in range (0, number):
            k = RHS.computex(RHS(CurrentState),i, CurrentState)*dt
            midpos[i,0] = (1/2)*k + CurrentState[i,0]
            
            if (Debugging):
                print(k, 'kx')
                
            k = RHS.computey(RHS(CurrentState),i, CurrentState)*dt
            midpos[i,1] = (1/2)*k +CurrentState[i,1]
            
            if (Debugging):
                print(k, 'ky')
                
            midpos[i,2] = newpoints[i,4]

        midstate = spline.interptopart(midpos, number, size)
        
        if (Debugging):
            print(midstate, 'midstate')
            
        for i in range(0,number):
            k = RHS.computex(RHS(midstate), i , midstate)*dt
            finalstate[i,0] = CurrentState[i,0] + k
            
            k = RHS.computey(RHS(midstate), i , midstate)*dt
            finalstate[i,1] = CurrentState[i,1] + k
            finalstate[i,2] = newpoints[i,4]

        return finalstate
        


# In[298]:


class move:
    def __init__(self, spline):
        self.spline = spline

    def forward(time, dt, number, CurrentState, spline):
        i = 0

        stuff = RK2.advance(number, dt , newpoints, spline)

        stuff = spline.interptopart(stuff, number, size)

        fulldata = stuff

        while i<time:

            stuff = RK2.advance(number, dt , stuff, spline)

            stuff = spline.interptopart(stuff, number, size)

            print(i, 'i')
            fulldata = np.vstack([fulldata, stuff] )

            i = i + dt
            print(fulldata, 'fdata')
            
        fulldata= np.delete(fulldata, slice(0, number), axis=0)
        column=np.zeros((((int((len(fulldata))),1))))
        fulldata = np.append(fulldata, column, axis=1)
        for i in range(0, len(fulldata)):
            fulldata[i,5] = math.floor(i/number)
        return fulldata


    


# In[ ]:


time = 3
dt = 1/16

#data4th = move.forward(time, dt, number, newpoints, dep4th)
data2nd = move.forward(time, dt, number, newpoints, dep2nd)



fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
plt.axis([0.1, 0.9, 0.1, 0.9])
i = 0
fig = plt.scatter(data4th[:,0], data4th[:,1], 10, data4th[:,5], cmap = 'autumn')
cbar = plt.colorbar(cmap = 'autumn', label = 'Time')
#plt.savefig('4thorder.png')


fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
plt.axis([0.1, 0.9, 0.1, 0.9])
i = 0
fig = plt.scatter(data2nd[:,0], data2nd[:,1], 10, data2nd[:,5], cmap = 'autumn')
cbar = plt.colorbar(cmap = 'autumn', label = 'Time')
#plt.savefig('2ndorder.png')



#for saving image at each timestep to make video

# plt.scatter(particles[:,0], particles[:,1], 10, 'r')
# plt.savefig("this.png", format="PNG")
# i =0
# while i<((len(fulldata))):
#     for k in range(0,90):
#         fig = plt.scatter(fulldata[k+i,0], fulldata[k+i,1], 10, 'r')
    
#     plt.savefig("that" + str(i) +".png", format="PNG")
#     print(i)
#     i = i + number


