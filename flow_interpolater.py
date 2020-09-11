import numpy as np


class flow_interpolater(object):
    
    def __init__(self, interp_file):
        file = np.load(interp_file)
        self.f_ux, self.f_uy, self.f_t = file['f_ux'], file['f_uy'], file['f_t']
        self.f_duxdt, self.f_duydt = file['f_duxdt'], file['f_duydt']
        self.t_max = self.f_t[-1]
        self.dt = self.f_t[1] - self.f_t[0]
        
    def t_index(self, t, approx):
        
        assert t <= self.t_max and t >= 0, "time out of range!"    # check if t <= t_max
        # if approx is false, t needs to be found in self.f_t
        if approx is False:
            idx_list = np.where(self.f_t == t)[0]
            assert len(idx_list) is not 0, "time (%.3f) needs to be multiple of %.3f!"%(t, self.dt)
            idx = idx_list[0]
        # if approx is true, we use the CLOSEST time in self.f_t
        else:
            idx = int(t/self.dt)
            
        return idx
        
    def vel(self, X, t, approx=True):
        
        t_idx = self.t_index(t, approx)
        Xr = np.mod(X, 2*np.pi)    # assume flow domain is 2pi periodic
        x, y = Xr[...,0], Xr[...,1]
        ux = self.f_ux[t_idx].ev(x, y)
        uy = self.f_uy[t_idx].ev(x, y)
        
        return np.stack((ux, uy), axis=-1)
    
    def grad_u(self, X, t, approx=True):
        # very crude estimate of the flow gradient; assume X consists of a single point - i.e. has shape (2,)
        grad = np.zeros((2, 2))
        grad[:,0] = (self.vel(X + np.array([1e-3, 0]), t) - self.vel(X - np.array([1e-3, 0]), t))/2e-3
        grad[:,1] = (self.vel(X + np.array([0, 1e-3]), t) - self.vel(X - np.array([0, 1e-3]), t))/2e-3
        
        return grad
    
    def DuDt(self, X, t, approx=True):
        
        t_idx = self.t_index(t, approx)
        Xr = np.mod(X, 2*np.pi)    # assume flow domain is 2pi periodic
        x, y = Xr[...,0], Xr[...,1]
        duxdt = self.f_duxdt[t_idx].ev(x, y)
        duydt = self.f_duydt[t_idx].ev(x, y)
        
        return np.stack((duxdt, duydt), axis=-1)
    
    def lmin(self, X, t, mu, approx=True):
        # estimate stability using self.grad_u; assume X consists of a single point - i.e. has shape (2,)
        grad = self.grad_u(X, t, approx)
        S = .5*(grad + grad.T) + np.array([[mu, 0],[0, mu]])
        S_eig, _ = np.linalg.eig(S)
        lmda_min = S_eig.min()
    
        return lmda_min

