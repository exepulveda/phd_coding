import numpy as np

class BlockModel(object):
    def __init__(self,nx,ny,nz,sx=1.0,sy=1.0,sz=1.0):
        self.nx= nx
        self.ny= ny
        self.nz= nz
        self.nxy= nx*ny
        self.n = self.nxy*nz
        self.sizex = sx
        self.sizey = sy
        self.sizez = sz
        self.volume = sx * sy * sz
        
    def get_blockid(self,i,j,k):
        return k*self.nxy + j*self.nx + i

    def get_volume(self):
        return self.volume
