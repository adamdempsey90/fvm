import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
class Sim2D():
    def __init__(self,num,base='kh_',with_ghost=False):
        with h5py.File(base + str(num) + '.h5','r') as file:
            f = file['/Data']
            self.time = float(f['time'][...])

            self.xm1 = f['xm1'][...]

            nx1 = len(self.xm1) - 1 -6
            self.nx1 = nx1
            if not with_ghost:
                self.xm1 = self.xm1[3:-3]
            self.xc1 = .5*(self.xm1[1:] + self.xm1[:-1])
            self.dx1 = np.diff(self.xm1)

            self.xm2 = f['xm2'][...]
            nx2 = len(self.xm2) - 1 - 6
            self.nx2 = nx2
            if not with_ghost:
                self.xm2 = self.xm1[3:-3]
            self.xc2 = .5*(self.xm2[1:] + self.xm2[:-1])
            self.dx2 = np.diff(self.xm2)

            self.gamma = float(f['Gamma'][...])
            self.rho = f['Density'][...].reshape(nx2+6,nx1+6)
            self.vx1 = f['Mx1'][...].reshape(nx2+6,nx1+6) / self.rho
            self.vx2 = f['Mx2'][...].reshape(nx2+6,nx1+6) / self.rho
            self.vx3 = f['Mx3'][...].reshape(nx2+6,nx1+6) / self.rho
            self.energy = f['Energy'][...].reshape(nx2+6,nx1+6)
            self.ke = .5*self.rho*(self.vx1**2 + self.vx2**2 + self.vx3**2)
            self.intenergy = f['InternalEnergy'][...].reshape((nx2+6,nx1+6))

            try:
                self.scalar = f['Scalar1'][...].reshape(nx2 + 6, nx1 +6) / self.rho
            except:
                pass

            self.extent = (self.xm1.min(),self.xm1.max(),self.xm2.min(),
                    self.xm2.max())

        if not with_ghost:
            self.rho = self.rho[3:-3,3:-3]
            self.vx1 = self.vx1[3:-3,3:-3]
            self.vx2 = self.vx2[3:-3,3:-3]
            self.vx3 = self.vx3[3:-3,3:-3]
            self.energy = self.energy[3:-3,3:-3]
            self.ke = self.ke[3:-3,3:-3]
            self.intenergy = self.intenergy[3:-3,3:-3]
            try:
                self.scalar = self.scalar[3:-3,3:-3]
            except AttributeError:
                pass
        self.pres = (self.energy-self.ke)*(self.gamma-1)
        self.cs = np.sqrt(self.gamma*self.pres/self.rho)
        self.temp = self.intenergy*self.gamma/self.rho
        delad = 1. - 1./self.gamma
        self.entropy = np.log(self.temp * self.pres**(-delad))
        self.vort = np.gradient(self.vx2,self.dx1[0],axis=1,edge_order=2) - np.gradient(self.vx1,self.dx2[0],axis=0,edge_order=2)
        if self.nan_check():
            print('NaN detected!')

        return
    def nan_check(self):
        func = lambda x: np.any(np.isnan(x))
        return func(self.vx1)|func(self.rho)|func(self.pres)|func(self.energy)
    def plot(self,val='rho',func = None,norm=None, fig=None,ax=None,ylbl='',
            cmap='viridis',conts=None,**kargs):
        first = ax is None
        if first:
            fig,ax=plt.subplots(figsize=(8,8))
        if func is not None:
            q = func(self)
        else:
            q = getattr(self,val)
        if norm is None:
            norm = colors.Normalize()
        img = ax.imshow(q.T,origin='lower',extent=self.extent,norm=norm,cmap=cmap,**kargs)
        if first:
            cb = _create_colorbar(ax,norm,cmap=cmap)
        else:
            cb = None
        ax.minorticks_on()
        ax.set_aspect('equal')
        ax.tick_params(labelsize=20)
        ax.text(.05,.05,'$t={:.2f}$'.format(self.time),transform=ax.transAxes,fontsize=20)
        if conts is not None:
            cont = ax.contour(q,levels=conts,origin='lower',extent=self.extent,norm=norm,colors='k',**kargs)
        return fig,ax,cb,img
    def plotavg(self,val='rho',axis=1,func = None,norm=1,shift=0, fig=None,ax=None,ylbl='',**kargs):
        if ax is None:
            fig,ax=plt.subplots(figsize=(8,6))
        if func is not None:
            q = func(self).mean(axis=axis)
        else:
            q = (getattr(self,val).mean(axis=axis) - shift) / norm

        if axis == 1 or axis == -1:
            x = self.xc2
        else:
            x = self.xc1
        print(q.shape,x.shape)
        ax.plot(x,q,**kargs)
        ax.minorticks_on()
        ax.tick_params(labelsize=20)
        return fig,ax

    def contour(self,val='rho',colorbar=False,func = None,norm=colors.Normalize(), fig=None,ax=None,ylbl='',**kargs):
        if ax is None:
            fig,ax=plt.subplots(figsize=(6,6))
        if func is not None:
            q = func(self)
        else:
            q = getattr(self,val)
        levels = kargs.pop('levels',None)
        if levels is None:
            cont = ax.contour(q,origin='lower',extent=self.extent,norm=norm,**kargs)
        else:
            cont = ax.contour(q,levels=levels,origin='lower',extent=self.extent,norm=norm,**kargs)

        if colorbar:
            cb = _create_colorbar(ax,norm)
        else:
            cb = None
        ax.minorticks_on()
        return fig,ax,cb,cont
    def presconts(self,conts,streams=True,fig=None,ax=None,**kargs):
        if ax is not None:
            self.plot('pres',fig=fig,ax=ax,**kargs)
        else:
            fig,ax,_,_ = self.plot('pres',**kargs)

        self.contour('rho',levels=conts,fig=fig,ax=ax,colors='k',clrbar=False)
        xlims = ax.get_xlim()
        ylims = ax.get_ylim()
        if streams:
            ax.streamplot(self.xc1,self.xc2,self.vx1,self.vx2,color='k')
        ax.set_xlim(xlims)
        ax.set_ylim(ylims)
        return fig,ax

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

    return cb

class Animation():
    def __init__(self,base='kh_',**kargs):
        self.base = base
        self.kargs = kargs
        self.fixbar=False
        self.fig,self.ax,self.cb,self.img = Sim2D(0,base=base).plot(**kargs)
    def update(self,i):
        fld = Sim2D(i,base=self.base)
        try:
            func =self.kargs['func']
            d = func(fld)
        except KeyError:
            try:
                val = self.kargs['val']
            except KeyError:
                val = 'rho'
                self.kargs['val'] = 'rho'
            d = getattr(fld,val)
        if not self.fixbar:
            self.cb.set_clim([d.min(),d.max()])
            self.cb.draw_all()
        self.img.set_data(d)
        self.ax.texts[0].remove()
        self.ax.text(.05,.05,'$t={:.2f}$'.format(fld.time),transform=self.ax.transAxes,fontsize=20)


    def animate(self,irange,fixbar=False,fname='mov',frames=None):
        self.fixbar = fixbar
        import matplotlib.animation as animation
        frames = len(irange)
        anim = animation.FuncAnimation(self.fig, self.update, frames=frames, repeat=False)
        anim.save('{}.mp4'.format(fname), writer=animation.FFMpegWriter())



