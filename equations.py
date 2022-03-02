
import spectral
import numpy as np
from scipy import sparse
import scipy.sparse.linalg as spla

class SoundWaves:

    def __init__(self, domain, u, p, p0):
        self.u = u
        self.p = p
        self.p0 = p0
        self.domain = domain
        self.dtype = dtype = u.dtype
        self.dudx = spectral.Field(domain, dtype=dtype)
        self.ux = spectral.Field(domain, dtype=dtype)
        self.u_RHS = spectral.Field(domain, dtype=dtype)
        self.p_RHS = spectral.Field(domain, dtype=dtype)
        
        self.problem = spectral.InitialValueProblem(domain, [u, p], [self.u_RHS, self.p_RHS],
                                                    num_BCs=2, dtype=dtype)
        
        prob = self.problem.pencils[0]
        
        self.N = N = domain.bases[0].N
        Z = np.zeros((N, N))
        
        diag = np.arange(N-1)+1
        D = sparse.diags(diag, offsets=1)
        length=domain.bases[0].interval[1]-domain.bases[0].interval[0]
        D = (2/length)*D
        
        self.D = D

        diag0 = np.ones(N)/2
        diag0[0] = 1
        diag2 = -np.ones(N-2)/2
        self.C = C = sparse.diags((diag0, diag2), offsets=(0,2))
        
        M = sparse.csr_matrix((2*N+2,2*N+2))
        M[:N, :N] = C
        M[N:2*N, N:2*N] = C
        prob.M = M
        
        # L matrix
        BC_rows = np.zeros((2, 2*N))
        i = np.arange(N)
        BC_rows[0, :N] = (-1)**i
        BC_rows[1, :N] = (+1)**i
        #BC_rows[0, N:2*N] = (-1)**i
        #BC_rows[1, N:2*N] = (+1)**i


        cols = np.zeros((2*N,2))
        cols[  N-1, 0] = 1
        cols[2*N-1, 1] = 1
        corner = np.zeros((2,2))

        Z = np.zeros((N, N))
        L = sparse.bmat([[Z, D],
                         [D, Z]])
        L = sparse.bmat([[      L,   cols],
                         [BC_rows, corner]])
        L = L.tocsr()
        prob.L = L
        self.t = 0

    def evolve(self, timestepper, dt, num_steps):
        ts = timestepper(self.problem)
        u = self.u
        ux = self.ux
        p = self.p
        p0 = self.p0
        p_RHS = self.p_RHS
        D = self.D
        BC_func = lambda t: [0, 0]

        for i in range(num_steps):
            u.require_coeff_space()
            ux.require_coeff_space()
            #ux.data = self.C @ u.data
            ux.data = self.D @ u.data
            ux.data = spla.spsolve(self.C,ux.data)
            #ux.data = self.C @ u.data
            
            p.require_coeff_space()
            p0.require_coeff_space()
            p_RHS.require_coeff_space()
            
            u.require_grid_space(scales=3/2)
            ux.require_grid_space(scales=3/2)
            p0.require_grid_space(scales=3/2)
            p_RHS.require_grid_space(scales=3/2)
            
            p_RHS.data = (1-p0.data) * ux.data
            
            p_RHS.require_coeff_space()
            p_RHS.data = self.C @ p_RHS.data
            u.require_coeff_space()
            ux.require_coeff_space()
            p0.require_coeff_space()
            
            ts.step(dt, BC_func(self.t))
            self.t += dt





class BurgersEquation:
    
    def __init__(self, domain, u, nu):
        dtype = u.dtype
        self.u = u
        self.u_RHS = spectral.Field(domain, dtype=dtype)
        self.dudx = spectral.Field(domain, dtype=dtype)
        self.problem = spectral.InitialValueProblem(domain, [u], [self.u_RHS], dtype=dtype)
        
        p = self.problem.pencils[0]
        x_basis = domain.bases[0]
        I = sparse.eye(x_basis.N)
        p.M = I
        D = x_basis.derivative_matrix(dtype)
        p.L = -nu*D@D
        
    def evolve(self, timestepper, dt, num_steps):
        ts = timestepper(self.problem)
        u = self.u
        dudx = self.dudx
        u_RHS = self.u_RHS
        for i in range(num_steps):
            dudx.require_coeff_space()
            u.require_coeff_space()
            dudx.data = u.differentiate(0)
            u.require_grid_space()
            dudx.require_grid_space()
            u_RHS.require_grid_space()
            u_RHS.data = -u.data*dudx.data
            ts.step(dt)


class KdVEquation:
    
    def __init__(self, domain, u):
        dtype = u.dtype
        self.dealias = 3/2
        self.u = u
        self.u_RHS = spectral.Field(domain, dtype=dtype)
        self.dudx = spectral.Field(domain, dtype=dtype)
        self.problem = spectral.InitialValueProblem(domain, [u], [self.u_RHS], dtype=dtype)
        
        p = self.problem.pencils[0]
        x_basis = domain.bases[0]
        I = sparse.eye(x_basis.N)
        p.M = I
        D = x_basis.derivative_matrix(dtype)
        p.L = D@D@D
        
    def evolve(self, timestepper, dt, num_steps):
        ts = timestepper(self.problem)
        u = self.u
        dudx = self.dudx
        u_RHS = self.u_RHS
        for i in range(num_steps):
            dudx.require_coeff_space()
            u.require_coeff_space()
            dudx.data = u.differentiate(0)
            u.require_grid_space(scales=self.dealias)
            dudx.require_grid_space(scales=self.dealias)
            u_RHS.require_grid_space(scales=self.dealias)
            u_RHS.data = 6*u.data*dudx.data
            ts.step(dt)


class SHEquation:

    def __init__(self, domain, u):
        dtype = u.dtype
        self.dealias = 2
        self.u = u
        self.u_RHS = spectral.Field(domain, dtype=dtype)
        self.problem = spectral.InitialValueProblem(domain, [u], [self.u_RHS], dtype=dtype)

        p = self.problem.pencils[0]
        x_basis = domain.bases[0]
        I = sparse.eye(x_basis.N)
        p.M = I
        D = x_basis.derivative_matrix(dtype)
        op = I + D@D
        p.L = op @ op + 0.3*I

    def evolve(self, timestepper, dt, num_steps):
        ts = timestepper(self.problem)
        u = self.u
        u_RHS = self.u_RHS
        for i in range(num_steps):
            u.require_coeff_space()
            u.require_grid_space(scales=self.dealias)
            u_RHS.require_grid_space(scales=self.dealias)
            u_RHS.data = 1.8*u.data**2 - u.data**3
            ts.step(dt)



