import h5py
import numpy as np
import matplotlib.pyplot as plt
class Sim1D():
    def __init__(self,num,base='sod_'):
        with h5py.File(base + str(num) + '.h5','r') as file:
            f = file['/Data']
            self.time = float(f['time'][...])
            self.xm1 = f['xm1'][...]
            self.xc1 = .5*(self.xm1[1:] + self.xm1[:-1])
            self.dx1 = np.diff(self.xm1)
            nx1 = len(self.xc1)
            self.gamma = float(f['Gamma'][...])
            self.rho = f['Density'][...].reshape((nx1,))
            self.vx1 = f['Mx1'][...].reshape((nx1,)) / self.rho
            self.vx2 = f['Mx2'][...].reshape((nx1,)) / self.rho
            self.vx3 = f['Mx3'][...].reshape((nx1,)) / self.rho
            self.energy = f['Energy'][...].reshape((nx1,))
            self.ke = .5*self.rho*(self.vx1**2 + self.vx2**2 + self.vx3**2)
            self.pres = (self.energy-self.ke)*(self.gamma-1)
            self.intenergy = f['InternalEnergy'][...].reshape((nx1,))
        if self.nan_check():
            print('NaN detected!')

        return
    def nan_check(self):
        func = lambda x: np.any(np.isnan(x))
        return func(self.vx1)|func(self.rho)|func(self.pres)|func(self.energy)
    def plot(self,val='rho',func = None, fig=None,ax=None,ylbl='',**kargs):
        if ax is None:
            fig,ax=plt.subplots(figsize=(8,6))
        if func is not None:
            q = func(self)
        else:
            q = getattr(self,val)
        ax.plot(self.xc1,q,**kargs)
        ax.minorticks_on()
        ax.set_ylabel(ylbl,fontsize=20)

    def sum(self,fig=None,axes=None,**kargS):
        if axes is None:
            fig,axes = plt.subplots(2,2,figsize=(8,8))
        self.plot(val='rho',ax=axes[0,0],fig=fig,ylbl='Density')
        self.plot(val='vx1',ax=axes[0,1],fig=fig,ylbl='Velocity')
        self.plot(val='pres',ax=axes[1,0],fig=fig,ylbl='Pressure')
        self.plot(func=lambda x: x.pres/(x.gamma-1)/x.rho,ax=axes[1,1],fig=fig,ylbl='Energy')
        fig.tight_layout()
        return fig,axes




