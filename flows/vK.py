import numpy as np

class vK(object):
    
    ndim = 2
    L = 2
    alpha = 2
    R0 = .35
    y0 = .3
    Tc = 1.107
    w = 8*12/np.pi
    u0 = 14/Tc

    @staticmethod
    # stream function
    def psi(X, t):
        x, y = X[...,0], X[...,1]
        return vK.f(x,y)*vK.g(x,y,t)

    # velocity
    def vel(X, t):
        x, y = X[...,0], X[...,1]
        vx, vy = vK.vx(x,y,t), vK.vy(x,y,t)
        return np.stack((vx, vy), axis=-1)

    # velocity gradient
    def u_grad(X, t):
        
        x, y = X[...,0], X[...,1]
        grad = np.zeros(np.array(x).shape + (2,2))
        grad[...,0,0] = vK.dvxdx(x,y,t)
        grad[...,0,1] = vK.dvxdy(x,y,t)
        grad[...,1,0] = vK.dvydx(x,y,t)
        grad[...,1,1] = vK.dvydy(x,y,t)
        
        return grad

    # material derivative
    def DuDt(X, t=0):
        
        x, y = X[...,0], X[...,1]
        duxdt = vK.dvxdt(x,y,t)
        duydt = vK.dvydt(x,y,t)
        
        u = vK.vel(X, t)
        nabla_u = vK.u_grad(X, t)
        DuDt = np.squeeze(np.matmul(nabla_u, np.expand_dims(u, -1)), axis=-1)
        
        DuDt[...,0] += duxdt
        DuDt[...,1] += duydt
        
        return DuDt

    # stability criterion
    def lmin(X, mu, t):
        
        grad = vK.u_grad(X, t)
        S = .5*(grad + np.swapaxes(grad, -1, -2))
        S[...,0,0] += mu
        S[...,1,1] += mu
        
        S_eig, _ = np.linalg.eig(S)
        lmda_min = S_eig.min(axis=-1)
        
        return lmda_min

    ## component functions ##
    def f(x,y):
        r = np.sqrt(x**2 + y**2)
        return 1 - np.exp(-(r-1)**2)

    def dfdx(x,y):
        r = np.sqrt(x**2 + y**2)
        return 2*x*(r-1)*np.exp(-(r-1)**2)/r

    def dfdy(x,y):
        r = np.sqrt(x**2 + y**2)
        return 2*y*(r-1)*np.exp(-(r-1)**2)/r

    def d2fdx2(x,y):
    	r = np.sqrt(x**2 + y**2)
    	ff = -2*x*(r-1)/r - x/r**2 + x/r/(r-1) + 1/x
    	return ff*vK.dfdx(x,y)

    def d2fdy2(x,y):
    	r = np.sqrt(x**2 + y**2)
    	ff = -2*y*(r-1)/r - y/r**2 + y/r/(r-1) + 1/y
    	return ff*vK.dfdy(x,y)

    def d2fdxdy(x,y):
    	r = np.sqrt(x**2 + y**2)
    	ff = -2*y*(r-1)/r - y/r**2 + y/r/(r-1)
    	return ff*vK.dfdx(x,y)

    def x1(t):
        return 1 + vK.L*(np.mod(t/vK.Tc, 1))

    def x2(t):
        return vK.x1(t - vK.Tc/2)

    def h1(t):
        return np.abs(np.sin(np.pi*t/vK.Tc))

    def dh1dt(t):
        return np.pi/vK.Tc*np.cos(np.pi/vK.Tc*np.mod(t, vK.Tc))

    def h2(t):
        return vK.h1(t - vK.Tc/2)

    def dh2dt(t):
        return vK.dh1dt(t - vK.Tc/2)

    def g1(x,y,t):
        return np.exp(-vK.R0*((x-vK.x1(t))**2 + vK.alpha**2*(y-vK.y0)**2))

    def dg1dx(x,y,t):
        return vK.g1(x,y,t)*(-2*vK.R0*(x-vK.x1(t)))

    def dg1dy(x,y,t):
        return vK.g1(x,y,t)*(-2*vK.R0*vK.alpha**2*(y-vK.y0))

    def dg1dt(x,y,t):
    	return vK.g1(x,y,t)*(2*vK.R0*(x-vK.x1(t))*vK.L/vK.Tc)

    def d2g1dx2(x,y,t):
    	ff = (-2*vK.R0*(x-vK.x1(t)))**2 - 2*vK.R0
    	return ff*vK.g1(x,y,t)

    def d2g1dy2(x,y,t):
    	ff = (-2*vK.R0*vK.alpha**2*(y-vK.y0))**2 - 2*vK.R0*vK.alpha**2
    	return ff*vK.g1(x,y,t)

    def d2g1dxdy(x,y,t):
    	ff = (-2*vK.R0*(x-vK.x1(t)))*(-2*vK.R0*vK.alpha**2*(y-vK.y0))
    	return ff*vK.g1(x,y,t)

    def d2g1dxdt(x,y,t):
    	return vK.dg1dt(x,y,t)*(-2*vK.R0*(x-vK.x1(t))) + vK.g1(x,y,t)*2*vK.R0*vK.L/vK.Tc

    def d2g1dydt(x,y,t):
    	return vK.dg1dt(x,y,t)*(-2*vK.R0*vK.alpha**2*(y-vK.y0))

    def g2(x,y,t):
        return np.exp(-vK.R0*((x-vK.x2(t))**2 + vK.alpha**2*(y+vK.y0)**2))

    def dg2dx(x,y,t):
        return vK.g2(x,y,t)*(-2*vK.R0*(x-vK.x2(t)))

    def dg2dy(x,y,t):
        return vK.g2(x,y,t)*(-2*vK.R0*vK.alpha**2*(y+vK.y0))

    def dg2dt(x,y,t):
    	return vK.g2(x,y,t)*(2*vK.R0*(x-vK.x2(t))*vK.L/vK.Tc)

    def d2g2dx2(x,y,t):
    	ff = (-2*vK.R0*(x-vK.x2(t)))**2 - 2*vK.R0
    	return ff*vK.g2(x,y,t)

    def d2g2dy2(x,y,t):
    	ff = (-2*vK.R0*vK.alpha**2*(y+vK.y0))**2 - 2*vK.R0*vK.alpha**2
    	return ff*vK.g2(x,y,t)

    def d2g2dxdy(x,y,t):
    	ff = (-2*vK.R0*(x-vK.x2(t)))*(-2*vK.R0*vK.alpha**2*(y+vK.y0))
    	return ff*vK.g2(x,y,t)

    def d2g2dxdt(x,y,t):
    	return vK.dg2dt(x,y,t)*(-2*vK.R0*(x-vK.x2(t))) + vK.g2(x,y,t)*2*vK.R0*vK.L/vK.Tc 

    def d2g2dydt(x,y,t):
    	return vK.dg2dt(x,y,t)*(-2*vK.R0*vK.alpha**2*(y+vK.y0))

    def s(x,y):
        return 1 - np.exp(-(x-1)**2/vK.alpha**2 - y**2)

    def dsdx(x,y):
        return np.exp(-(x-1)**2/vK.alpha**2 - y**2)*2*(x-1)/vK.alpha**2

    def dsdy(x,y):
        return np.exp(-(x-1)**2/vK.alpha**2 - y**2)*2*y

    def d2sdx2(x,y):
    	ff = 2*(vK.alpha**2 - 2*(x-1)**2)/vK.alpha**4
    	return ff*np.exp(-(x-1)**2/vK.alpha**2 - y**2)

    def d2sdy2(x,y):
    	ff = 2 - 4*y**2
    	return ff*np.exp(-(x-1)**2/vK.alpha**2 - y**2)

    def d2sdxdy(x,y):
    	ff = 4*(1-x)*y/vK.alpha**2
    	return ff*np.exp(-(x-1)**2/vK.alpha**2 - y**2)

    def g(x,y,t):
        return -vK.w*vK.h1(t)*vK.g1(x,y,t) + vK.w*vK.h2(t)*vK.g2(x,y,t) + vK.u0*y*vK.s(x,y)

    def dgdx(x,y,t):
        return -vK.w*vK.h1(t)*vK.dg1dx(x,y,t) + vK.w*vK.h2(t)*vK.dg2dx(x,y,t) + vK.u0*y*vK.dsdx(x,y)

    def dgdy(x,y,t):
        return -vK.w*vK.h1(t)*vK.dg1dy(x,y,t) + vK.w*vK.h2(t)*vK.dg2dy(x,y,t) + vK.u0*(vK.s(x,y) + y*vK.dsdy(x,y))

    def dgdt(x,y,t):
    	return vK.w*(-vK.h1(t)*vK.dg1dt(x,y,t) - vK.dh1dt(t)*vK.g1(x,y,t) + vK.h2(t)*vK.dg2dt(x,y,t) + vK.dh2dt(t)*vK.g2(x,y,t))

    def d2gdxdt(x,y,t):
    	return vK.w*(-vK.h1(t)*vK.d2g1dxdt(x,y,t) - vK.dh1dt(t)*vK.dg1dx(x,y,t) + vK.h2(t)*vK.d2g2dxdt(x,y,t) + vK.dh2dt(t)*vK.dg2dx(x,y,t))

    def d2gdydt(x,y,t):
    	return vK.w*(-vK.h1(t)*vK.d2g1dydt(x,y,t) - vK.dh1dt(t)*vK.dg1dy(x,y,t) + vK.h2(t)*vK.d2g2dydt(x,y,t) + vK.dh2dt(t)*vK.dg2dy(x,y,t))

    def d2gdx2(x,y,t):
    	return -vK.w*vK.h1(t)*vK.d2g1dx2(x,y,t) + vK.w*vK.h2(t)*vK.d2g2dx2(x,y,t) + vK.u0*y*vK.d2sdx2(x,y)

    def d2gdxdy(x,y,t):
    	return -vK.w*vK.h1(t)*vK.d2g1dxdy(x,y,t) + vK.w*vK.h2(t)*vK.d2g2dxdy(x,y,t) + vK.u0*(vK.dsdx(x,y) + y*vK.d2sdxdy(x,y))

    def d2gdy2(x,y,t):
    	return -vK.w*vK.h1(t)*vK.d2g1dy2(x,y,t) + vK.w*vK.h2(t)*vK.d2g2dy2(x,y,t) + vK.u0*(2*vK.dsdy(x,y) + y*vK.d2sdy2(x,y))

    def vx(x,y,t):
        return vK.dfdy(x,y)*vK.g(x,y,t) + vK.dgdy(x,y,t)*vK.f(x,y)

    def vy(x,y,t):
        return -vK.dfdx(x,y)*vK.g(x,y,t) - vK.dgdx(x,y,t)*vK.f(x,y)
        
    def dvxdt(x,y,t):
    	return vK.dfdy(x,y)*vK.dgdt(x,y,t) + vK.d2gdydt(x,y,t)*vK.f(x,y)

    def dvydt(x,y,t):
    	return -vK.dfdx(x,y)*vK.dgdt(x,y,t) - vK.d2gdxdt(x,y,t)*vK.f(x,y)

    def dvxdx(x,y,t):
    	return vK.d2fdxdy(x,y)*vK.g(x,y,t) + vK.dfdy(x,y)*vK.dgdx(x,y,t) + vK.d2gdxdy(x,y,t)*vK.f(x,y) + vK.dgdy(x,y,t)*vK.dfdx(x,y)

    def dvxdy(x,y,t):
    	return vK.d2fdy2(x,y)*vK.g(x,y,t) + vK.dfdy(x,y)*vK.dgdy(x,y,t) + vK.d2gdy2(x,y,t)*vK.f(x,y) + vK.dgdy(x,y,t)*vK.dfdy(x,y)

    def dvydx(x,y,t):
    	return -vK.d2fdx2(x,y)*vK.g(x,y,t) - vK.dfdx(x,y)*vK.dgdx(x,y,t) - vK.d2gdx2(x,y,t)*vK.f(x,y) - vK.dgdx(x,y,t)*vK.dfdx(x,y)

    def dvydy(x,y,t):
    	return -vK.dvxdx(x,y,t)