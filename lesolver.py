import numpy as np
from scipy.integrate import cumtrapz, solve_ivp

# Defining an event which causes the termination of the integration when it passes 0 in the negative direction
def zero_crossing(t, y):
    return y[0]
zero_crossing.terminal = True
zero_crossing.direction = -1

class LaneEmdenSolver():
    def __init__(self, alpha, rho0, nnn):
        self.alpha = alpha
        self.rho0 = rho0
        self.nnn = nnn
    
    def set_core_mass(self, mass_ratio):
        self.core_mass = self.bcs.m * mass_ratio
        self.rhobar = 3. * self.core_mass / (4. * np.pi * self.bcs.r**3)

    def set_kernel(self, kernel):
        self.kernel = kernel

    def update_profile(self):
        xi, theta, eta = self.evaluate()
        self.xi    = xi
        self.theta = theta
        self.eta   = eta

        self.mass = cumtrapz(4. * np.pi * alphacgs**2 * xi**2 * rho0 * theta**nnn, alphacgs * xi, initial=0.)

    def set_boundary_conditions(self, bcs):
        self.bcs = bcs

    # The modified Lane-Emden equation
    def modified_lane_emden(t, y, alpha, rho0):
        theta, eta = y
        xi         = t

        dthetadxi = -eta/(xi**2)

        u, chi, dchidu = self.kernel(alpha * xi)
        
        dfdt0 = -eta / xi**2
        dfdt1 = xi**2 * (theta**self.nnn + self.rhobar / rho0 * (chi + 1./3. * u * dchidu))

        dfdt = np.array([dfdt0, dfdt1])

        return dfdt

    def evaluate(self, alpha=self.alpha, rho0=self.rho0):
        t_span = (tmin, xi_max)
        y0     = [1., 0.] #Initial conditions for integration
        t_eval = np.linspace(tmin, xi_max, 1000)

        fun = lambda t, y: self.modified_lane_emden(t, y, alpha, rho0)

        sol = solve_ivp(fun=fun, t_span=t_span, y0=y0, t_eval=t_eval, events=zero_crossing, dense_output=True)

        xi    = np.array(sol.t)
        theta = np.array(sol.y[0])
        eta   = np.array(sol.y[1])

        return xi, eta, theta

    def get_cut_rho(self, alpha=self.alpha, rho0=self.rho0):
        xi, eta, theta = self.evaluate(alpha, rho0)

        r = alpha * xi
        icut = (np.abs(r-rcut)).argmin()
        rho = rho0 * theta[icut]**self.nnn
        rhop = rho0 * theta[icut-1]**self.nnn
        drhodr = (rho - rhop) / (r[icut] - r[icut-1])  #This could probably be improved
        return rho, drhodr

    def drhodx(self, x, name):
        h = np.sqrt(sys.float_info.epsilon) * x
        rhon2, drhon2 = self.get_cut_rho(*{name:x-2*h})
        rhon1, drhon1 = self.get_cut_rho(*{name:x-h})
        rhop1, drhop1 = self.get_cut_rho(*{name:x+h})
        rhop2, drhop2 = self.get_cut_rho(*{name:x+2*h})
        drho = (rhon2 - 8.*rhon1 + 8.*rhop1 - rhop2) / (12.*h)
        ddrho = (drhon2 - 8.*drhon1 + 8.*drhop1 - drhop2) / (12.*h)
        return drho, ddrho




# Abstract out anything above NR, i.e. eval_rk, drhodalpha etc.
# Hmmm...
