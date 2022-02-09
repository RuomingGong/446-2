
import numpy as np
import spectral
from scipy import sparse

class KdVEquation:

    def __init__(self, domain, u):
        self.dtype = u.dtype
        dtype = self.dtype
        self.u = u
        self.domain = domain
        self.dudx = spectral.Field(domain, dtype=dtype)
        self.RHS = spectral.Field(domain, dtype=dtype) # -u*dudx
        self.problem = spectral.InitialValueProblem(domain, [self.u], [self.RHS], dtype=dtype)

        p = self.problem.pencils[0]

        x_basis = domain.bases[0]
        I = sparse.eye(x_basis.N, dtype=dtype)
        p.M = I
        if self.dtype == np.complex128:
            diag = -1j*x_basis.wavenumbers(dtype)**3 
            p.L = sparse.diags(diag)
        if self.dtype == np.float64:
            upper_diag = np.zeros(x_basis.N-1)
            upper_diag[::2] = x_basis.wavenumbers(dtype)[::2]
            upper_diag = upper_diag**3
            lower_diag = -upper_diag
            diag = sparse.diags([upper_diag, lower_diag], offsets=(1,-1))
            p.L = diag

    def evolve(self, timestepper, dt, num_steps):
        ts = timestepper(self.problem)
        x_basis = self.domain.bases[0]
        u = self.u
        dudx = self.dudx
        RHS = self.RHS

        for i in range(num_steps):
            # need to calculate -u*ux and put it into RHS
            u.require_coeff_space()
            dudx.require_coeff_space()
            if self.dtype == np.complex128:
                dudx.data = 1j*x_basis.wavenumbers(self.dtype)*u.data 
            if self.dtype == np.float64:
                dudx.data[::2] = -x_basis.wavenumbers(self.dtype)[::2]*u.data[1::2]
                dudx.data[1::2] = x_basis.wavenumbers(self.dtype)[::2]*u.data[::2]
            u.require_grid_space(scales=3/2)
            dudx.require_grid_space(scales=3/2)
            RHS.require_grid_space(scales=3/2)
            RHS.data = 6*u.data * dudx.data

            # take timestep
            ts.step(dt)

class SHEquation:

    def __init__(self, domain, u):
        self.dtype = u.dtype
        dtype = self.dtype
        self.u = u
        self.domain = domain
        self.dudx = spectral.Field(domain, dtype=dtype)
        self.RHS = spectral.Field(domain, dtype=dtype) # -u*dudx
        self.problem = spectral.InitialValueProblem(domain, [self.u], [self.RHS], dtype=dtype)

        p = self.problem.pencils[0]

        x_basis = domain.bases[0]

        r=-0.3

        I = sparse.eye(x_basis.N, dtype=dtype)
        p.M = I

        diag = (1-r) + (-2)*x_basis.wavenumbers(dtype)**2 + x_basis.wavenumbers(dtype)**4
        p.L = sparse.diags(diag)

    def evolve(self, timestepper, dt, num_steps):
        ts = timestepper(self.problem)
        x_basis = self.domain.bases[0]
        u = self.u
        dudx = self.dudx
        RHS = self.RHS

        for i in range(num_steps):
            # need to calculate -u*ux and put it into RHS
            u.require_coeff_space()          
            u.require_grid_space(scales=3)
            RHS.require_grid_space(scales=3)
            RHS.data = -(u.data)**3+1.8*(u.data)**2

            # take timestep
            ts.step(dt)


