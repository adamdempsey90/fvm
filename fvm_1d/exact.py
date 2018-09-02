import numpy as np
from scipy.optimize import fsolve

def fk_prime(ps,d,p,a,g):
    dL,dR = d[0],d[1]
    pL,pR = p[0],p[1]
    aL,aR = a[0],a[1]

    Al  = 2./(g+1) / dL
    Ar = 2./(g+1) / dR
    Bl = (g-1)/(g+1) * pL
    Br = (g-1)/(g+1) * pR
    if ps > pL:
        fL = np.sqrt(Al/(ps + Bl))*(1 - .5*(ps-pL)/(Bl +ps))
    else:
        fL = 1./(aL*dL) *  (ps/pL)**(-(g+1)/(2*g))
    if ps > pR:
        fR = np.sqrt(Ar/(ps + Br))*(1 - .5*(ps-pR)/(Br +ps))
    else:
        fR = 1./(aR*dR) *  (ps/pR)**(-(g+1)/(2*g))

    return fL,fR

def fk_star(ps,d,p,a,g):
    dL,dR = d[0],d[1]
    pL,pR = p[0],p[1]
    aL,aR = a[0],a[1]

    Al  = 2./(g+1) / dL
    Ar = 2./(g+1) / dR
    Bl = (g-1)/(g+1) * pL
    Br = (g-1)/(g+1) * pR
    if ps > pL:
        fL =  (ps-pL)*np.sqrt(Al/(ps+Bl))
    else:
        fL = 2*aL/(g-1) * ((ps/pL)**((g-1)/(2*g)) - 1.)
    if ps > pR:
        fR =  (ps-pR)*np.sqrt(Ar/(ps+Br))
    else:
        fR = 2*aR/(g-1)*((ps/pR)**((g-1)/(2*g)) -1.)

    return fL,fR


def star_region(d,p,u,g):
    dL,dR = d[0],d[1]
    pL,pR = p[0],p[1]
    uL,uR = u[0],u[1]
    aL,aR = np.sqrt(g*pL/dL), np.sqrt(g*pR/dR)

    du = uR-uL

    ps = fsolve(lambda ps: sum(fk_star(ps,d,p,(aL,aR),g)) +du, .5*(pL+pR),
            fprime = lambda ps: sum(fk_prime(ps,d,p,(aL,aR),g)),
            factor=.1)[0]
    fL,fR = fk_star(ps,d,p,(aL,aR),g)
    us = .5*(uL+uR) + .5*(fR-fL)
    return ps,us


def exact_solve(x,d,p,u,t,g):

    S = x/t
    vel = np.zeros(S.shape)
    dens = np.zeros(S.shape)
    pres = np.zeros(S.shape)

    pL,pR = p[0],p[1]
    dL,dR = d[0],d[1]
    uL,uR = u[0],u[1]
    aL,aR = np.sqrt(g*pL/dL), np.sqrt(g*pR/dR)


    ps,us = star_region(d,p,u,g)
    print("ps, us ", ps,us)

    # Left of contact
    ind_left = S < us

    if  ps > pL:
        # Shock Wave
        SL = uL - aL*np.sqrt( (g+1)/(2*g) * ps/pL + (g-1)/(2*g))
        indL = (ind_left)&(S<=SL)
        indR = (ind_left)&(S>SL)

        vel[indL] = uL
        pres[indL] = pL
        dens[indL] = dL

        vel[indR] = us
        pres[indR] = ps
        dens[indR] = dL *( ps/pL + (g-1)/(g+1))/((g-1)/(g+1) * ps/pL + 1)

    else:
        # Rarefaction
        SHL = uL-aL
        STL = us - aL*(ps/pL)**((g-1)/(2*g))
        indL = (ind_left)&(S <= SHL)
        indM = (ind_left)&(S >= SHL)&(S<=STL)
        indR = (ind_left)&(S >= STL)

        vel[indL] = uL
        pres[indL] = pL
        dens[indL] = dL

        sm = S[indM]
        vel[indM]  =2./(g+1) *(aL + (g-1)*uL/2. + sm)
        dens[indM] = dL*(2./(g+1) + (g-1)*(uL-sm)/(aL*(g+1)))**(2./(g-1))
        pres[indM] = pL *(2./(g+1) + (g-1)/(aL*(g+1)) *(uL-sm))**(2*g/(g-1))

        vel[indR] = us
        pres[indR] = ps
        dens[indR] = dL * (ps/pL)**(1./g)




    # Right of contact
    ind_right = ~ind_left

    if ps > pR:
        # Shock wave
        SR = uR + aR*np.sqrt((g+1)/(2*g) * ps/pR + (g-1)/(2*g))
        indL  = (ind_right)&(S <= SR)
        indR = (ind_right)&(S >= SR)

        pres[indL] = ps
        vel[indL] = us
        dens[indL] = dR *( ps/pR + (g-1)/(g+1))/((g-1)/(g+1) *ps/pR + 1)

        pres[indR] = pR
        vel[indR] = uR
        dens[indR] = dR
    else:
        # Rarefaction
        SHR = uR +aR
        STR =  us + aR*(ps/pR)**((g-1)/(2*g))


        indL = (ind_right)&(S <= STR)
        indM = (ind_right)&(S>=STR)&(S<=SHR)
        indR = (ind_right)&(S >= SHR)

        vel[indL] = us
        pres[indL] = ps
        dens[indL] = dR*(ps/pR)**(1./g)

        sm = S[indM]
        vel[indM] =2./(g+1) *( -aR + (g-1)*uR/2  + sm)
        pres[indM] = pR*(2./(g+1) - (g-1)/(g+1) *(uR-sm)/aR)**(2*g/(g-1))
        dens[indM] = dR*(2./(g+1) - (g-1)/(g+1)*(uR-sm)/aR)**(2./(g-1))

        vel[indR] = uR
        pres[indR] = pR
        dens[indR] = dR

    return vel,pres,dens


def plot_exact(x,vel,pres,dens,g,fig=None,axes=None):
    import matplotlib.pyplot as plt
    if axes is None:
        fig,axes = plt.subplots(2,2,figsize=(8,8))

    energ = pres/(g-1)/ dens

    axes[0,0].plot(x,dens,'-k')
    axes[0,1].plot(x,vel,'-k')
    axes[1,0].plot(x,pres,'-k')
    axes[1,1].plot(x,energ,'-k')

    axes[0,0].set_ylabel('Density',fontsize=20)
    axes[0,1].set_ylabel('Velocity',fontsize=20)
    axes[1,0].set_ylabel('Pressure',fontsize=20)
    axes[1,1].set_ylabel('Internal Energy',fontsize=20)

    for ax in axes.flatten():
        ax.minorticks_on()
    fig.tight_layout()
    return fig,axes


def toro_test_ics(num):
    tests = { 1: {'p': (1.,.1), 'd': (1.,.125), 'u': (0.,0.), 'ps': .30313, 't': .25},
            2: {'p': (.4,.4), 'd': (1.,1.), 'u': (-2.,2.), 'ps':.00189, 't': .15},
            3: {'p': (1000.,.01), 'd': (1.,1.), 'u': (0.,0.), 'ps': 460.894, 't': .012},
            4: {'p': (.01,100.), 'd': (1.,1.), 'u': (0.,0.), 'ps': 46.095, 't': .035},
            5: {'p': (460.894,46.0950), 'd': (5.99924,5.99242), 'u': (19.5975,-6.19633), 'ps': 1691.64, 't': .035},
            }
    return tests[num]

def plot_toro_test(num=1,g=1.4, nx=10000):
    x = np.linspace(-.5,.5,nx)

    states = toro_test_ics(num)

    vel,pres,dens = exact_solve(x,states['d'],states['p'],states['u'],states['t'],g)
    fig,ax=plot_exact(x+.5,vel,pres,dens,g)
    return fig,ax



def test_exact(num=1,g=1.4,nx=10000):
    fig,axes = plot_toro_test(num=num,g=g,nx=nx)
    dat = np.fromfile('test{:d}.dat'.format(num))
    x = dat[:nx]; d = dat[nx:2*nx]; u = dat[2*nx:3*nx]; p = dat[3*nx:4*nx];
    x += .5
    axes[0,0].plot(x,d,'.-b')
    axes[0,1].plot(x,u,'.-b')
    axes[1,0].plot(x,p,'.-b')
    axes[1,1].plot(x,p/d/(g-1),'.-b')


def test_star():
    g = 1.4
    for i in range(1,6):
        ics = toro_test_ics(i)
        d = ics['d']
        p = ics['p']
        u = ics['u']
        ans = ics['ps']
        ps,us = star_region(d,p,u,g)
        err = abs(ps-ans)/ans
        print('Test {:d}, Relative Error {:.3e}'.format(i,err))
