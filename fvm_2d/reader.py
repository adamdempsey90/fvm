import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
class Sim2D():
    def __init__(self,num,base='kh_'):
        with h5py.File(base + str(num) + '.h5','r') as file:
            f = file['/Data']
            self.time = float(f['time'][...])

            self.xm1 = f['xm1'][...]
            self.xc1 = .5*(self.xm1[1:] + self.xm1[:-1])
            self.dx1 = np.diff(self.xm1)
            nx1 = len(self.xc1)

            self.xm2 = f['xm2'][...]
            self.xc2 = .5*(self.xm2[1:] + self.xm2[:-1])
            self.dx2 = np.diff(self.xm2)
            nx2 = len(self.xc2)

            self.gamma = float(f['Gamma'][...])
            self.rho = f['Density'][...].reshape(nx1+6,nx2+6)
            self.vx1 = f['Mx1'][...].reshape(nx1+6,nx2+6) / self.rho
            self.vx2 = f['Mx2'][...].reshape(nx1+6,nx2+6) / self.rho
            self.vx3 = f['Mx3'][...].reshape(nx1+6,nx2+6) / self.rho
            self.energy = f['Energy'][...].reshape(nx1+6,nx2+6)
            self.ke = .5*self.rho*(self.vx1**2 + self.vx2**2 + self.vx3**2)
            self.pres = (self.energy-self.ke)*(self.gamma-1)
            self.intenergy = f['InternalEnergy'][...].reshape((nx1+6,nx2+6))
            self.extent = (self.xm1.min(),self.xm1.max(),self.xm2.min(),
                    self.xm2.max())

        self.rho = self.rho[3:-3,3:-3]
        self.pres = self.pres[3:-3,3:-3]
        self.vx1 = self.vx1[3:-3,3:-3]
        self.vx2 = self.vx2[3:-3,3:-3]
        self.vx3 = self.vx3[3:-3,3:-3]
        self.pres = self.pres[3:-3,3:-3]
        self.ke = self.ke[3:-3,3:-3]
        self.intenergy = self.intenergy[3:-3,3:-3]
        if self.nan_check():
            print('NaN detected!')

        return
    def nan_check(self):
        func = lambda x: np.any(np.isnan(x))
        return func(self.vx1)|func(self.rho)|func(self.pres)|func(self.energy)
    def plot(self,val='rho',func = None,norm=None, fig=None,ax=None,ylbl='',
            cmap='viridis',**kargs):
        if ax is None:
            fig,ax=plt.subplots(figsize=(6,6))
        if func is not None:
            q = func(self)
        else:
            q = getattr(self,val)
        if norm is None:
            norm = colors.Normalize()
        ax.imshow(q,origin='lower',extent=self.extent,norm=norm,cmap=cmap,**kargs)
        _create_colorbar(ax,norm,cmap=cmap)
        ax.minorticks_on()
        ax.set_aspect('equal')

    def contour(self,val='rho',func = None,norm=colors.Normalize(), fig=None,ax=None,ylbl='',**kargs):
        if ax is None:
            fig,ax=plt.subplots(figsize=(6,6))
        if func is not None:
            q = func(self)
        else:
            q = getattr(self,val)
        cont = ax.contour(q.T,origin='lower',extent=self.extent,norm=norm,**kargs)
        _create_colorbar(ax,norm)
        ax.minorticks_on()
        return fig,ax,cont
    def sum(self,fig=None,axes=None,**kargS):
        if axes is None:
            fig,axes = plt.subplots(2,2,figsize=(8,8))
        self.plot(val='rho',ax=axes[0,0],fig=fig,ylbl='Density')
        self.plot(val='vx1',ax=axes[0,1],fig=fig,ylbl='Velocity')
        self.plot(val='pres',ax=axes[1,0],fig=fig,ylbl='Pressure')
        self.plot(func=lambda x: x.pres/(x.gamma-1)/x.rho,ax=axes[1,1],fig=fig,ylbl='Energy')
        fig.tight_layout()
        return fig,axes



def _create_colorbar(ax,norm,cax=None,log=False,cmap='viridis',**kargs):
    import matplotlib
    import matplotlib.cm
    from mpl_toolkits.axes_grid1 import make_axes_locatable


    labelsize = kargs.pop('labelsize',12)

    if cax is None:
        divider = make_axes_locatable(ax)
        cax = divider.append_axes('top',size='3%',pad=.05)


    cmap = matplotlib.cm.get_cmap(cmap)
    cb = matplotlib.colorbar.ColorbarBase(ax=cax,cmap=cmap,norm=norm,orientation='horizontal',**kargs)
    cb.ax.xaxis.set_ticks_position('top')
    cb.ax.xaxis.set_label_position('top')
    cb.ax.tick_params(labelsize=labelsize)

    return


def animate_contour(irange,**kargs):
    import matplotlib.animation as animation
    nt = len(irange)
    fld = Sim2D(irange[0],base='imp_')
    fig,ax,cont = fld.contour(val='rho',**kargs)


    def update(i):
        for c in cont.collections:
            c.remove()
        fld = Sim2D(i,base='imp_')
        fig,ax,cont =fld.contour(val='rho',fig=fig,ax=ax,**kargs)
        return cont

    anim = animation.FuncAnimation(fig, update, frames=nt, repeat=False)
    anim.save('implosion.mp4', writer=animation.FFMpegWriter())
