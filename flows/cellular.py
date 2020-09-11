import numpy as np

class cellular_flow(object):
    
    ndim = 2

    def __init__(self, B, k=1, w=np.pi, A=15):
        self.B = B
        self.k = k
        self.w = w
        self.A = A
        
    def vel(self, X, t=0):
        
        x, y = X[...,0], X[...,1]
        vx = self.A * np.sin(self.k*y) * np.cos(self.k*x + self.B*np.sin(self.w*t))
        vy = - self.A * np.sin(self.k*x + self.B*np.sin(self.w*t)) * np.cos(self.k*y)
        
        return np.stack((vx, vy), axis=-1)
    
    def u_grad(self, X, t=0):
        
        x, y = X[...,0], X[...,1]
        grad = np.zeros(np.array(x).shape + (2,2))
        grad[...,0,0] = -self.A*self.k*np.sin(self.k*y)*np.sin(self.k*x + self.B*np.sin(self.w*t))
        grad[...,0,1] = self.A*self.k*np.cos(self.k*y)*np.cos(self.k*x+ self.B*np.sin(self.w*t))
        grad[...,1,0] = -self.A*self.k*np.cos(self.k*x + self.B*np.sin(self.w*t))*np.cos(self.k*y)
        grad[...,1,1] = self.A*self.k*np.sin(self.k*x + self.B*np.sin(self.w*t))*np.sin(self.k*y)
        
        return grad
    
    def DuDt(self, X, t=0):
        
        if self.B != 0:
            x, y = X[...,0], X[...,1]
            duxdt = -self.A*self.B*self.w*np.sin(self.k*y)*np.cos(self.w*t)*np.sin(self.k*x+self.B*(np.sin(self.w*t)))
            duydt = -self.A*self.B*self.w*np.cos(self.k*y)*np.cos(self.w*t)*np.cos(self.k*x+self.B*(np.sin(self.w*t)))
        else:
            duxdt, duydt = 0, 0

        u = self.vel(X, t)
        nabla_u = self.u_grad(X, t)
        DuDt = np.squeeze(np.matmul(nabla_u, np.expand_dims(u, -1)), axis=-1)
        
        DuDt[...,0] += duxdt
        DuDt[...,1] += duydt
        
        return DuDt
    
    def Q(self, X, t=0):
        
        grad = self.u_grad(X, t)
        return grad[...,0,1]*grad[...,1,0] - grad[...,0,0]*grad[...,1,1]
    
    def lmin(self, X, mu, t=0):
        
        grad = self.u_grad(X, t)
        S = .5*(grad + np.swapaxes(grad, -1, -2))
        S[...,0,0] += mu
        S[...,1,1] += mu
        
        S_eig, _ = np.linalg.eig(S)
        lmda_min = S_eig.min(axis=-1)
        
        return lmda_min
