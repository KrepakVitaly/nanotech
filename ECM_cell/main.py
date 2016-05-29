__author__ = 'zorroxied'

import numpy as np
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
import matplotlib.pyplot as plt
from PyQt4 import Qt
import sys
from scipy import constants as const
from PyQt4 import QtGui
from PyQt4 import QtCore
from PyQt4.QtCore import QObject, pyqtSignal, pyqtSlot
import random
from sympy.solvers import solve
from sympy import Symbol, nsolve
from sympy import exp, sinh, ln, asinh, sqrt
import mpmath


class Cell():
    def __init__(self):
        # Simulation model parameters
        self.M_me       = 1.79e-25  # kg
        self.z          = 1.0  # n/a
        self.rho_me     = 10.49e3  # kg m-3
        self.m_r        = 0.023  # n/a
        self.alpha      = 0.3  # n/a
        self.j_0et      = 3.2e5  # A m-2
        self.DeltaG_et  = 0.6 * const.elementary_charge  # Joule
        self.j_0hop     = 1.1e11  # A m-2
        self.a          = 0.25e-9  # m
        self.DeltaG_hop = 0.32 * const.elementary_charge  # Joule
        self.DeltaG_nuc = 0.80 * const.elementary_charge  # Joule
        self.t_0nuc     = 2e-8  # s
        self.N_c        = 3.0  # n/a
        self.A_ac       = 804.25e-18  # m2
        self.A_fil      = 12.57e-18  # m2
        self.A_is       = 12.57e-18  # m2
        self.L          = 30.0e-9  # m
        self.rho_fil    = 1.7e-8  # Ohm m
        self.R_el       = 76.4  # Ohm
        self.R_S        = 1.0e6  # Ohm

        self.C = 2.7
        self.T = 300  # K

    def eta_fil(self, x, V_app, apprx=(0, 0, 0, 0)):
        m_eff = self.m_r * const.electron_mass

        mpmath.mp.dps = 20
        x0 = Symbol('x0')  # eta_fil
        x1 = Symbol('x1')  # eta_ac
        x2 = Symbol('x2')  # eta_hop
        x3 = Symbol('x3')  # V_tunnel

        f0 = const.Boltzmann * self.T / (1 - self.alpha) / const.elementary_charge / self.z * \
             ln(self.A_fil/self.A_ac*(exp(- self.alpha * const.elementary_charge * self.z / const.Boltzmann / self.T * x0) - 1) + 1) - x1# eta_ac = f(eta_fil) x1 = f(x0)
        f1 = x*2*const.Boltzmann*self.T/self.a/self.z/const.elementary_charge*\
             asinh(self.j_0et/self.j_0hop*(exp(- self.alpha * const.elementary_charge * self.z / const.Boltzmann / self.T * x0) - 1)) - x2# eta_hop = f(eta_fil)
        f2 = x1 - x0 + x2 - x3

        f3 = -V_app + ((self.C * 3 * sqrt(2 * m_eff * ((4+x3/2)*const.elementary_charge)) / 2 / x * (const.elementary_charge / const.Planck)**2 * \
             exp(- 4 * const.pi * x / const.Planck * sqrt(2 * m_eff * ((4+x3/2)*const.elementary_charge))) * self.A_fil*x3)
                       + (self.j_0et*self.A_fil*(exp(-self.alpha*const.elementary_charge*self.z*x0/const.Boltzmann/self.T) - 1))) * (self.R_el + self.R_S + self.rho_fil*(self.L - x) / self.A_fil) \
             + x3

        eta_fil, eta_ac, eta_hop, V_tunnel = nsolve((f0, f1, f2, f3), [x0, x1, x2, x3], apprx)
        eta_fil = np.real(np.complex128(eta_fil))
        eta_ac = np.real(np.complex128(eta_ac))
        eta_hop = np.real(np.complex128(eta_hop))
        V_tunnel = np.real(np.complex128(V_tunnel))
        current = ((self.C * 3 * sqrt(2 * m_eff * ((4+V_tunnel)*const.elementary_charge)) / 2 / x * (const.elementary_charge / const.Planck)**2 * \
            exp(- 4 * const.pi * x / const.Planck * sqrt(2 * m_eff * ((4+V_tunnel)*const.elementary_charge))) * self.A_fil*V_tunnel)
                       + (self.j_0et*self.A_fil*(exp(-self.alpha*const.elementary_charge*self.z*eta_fil/const.Boltzmann/self.T) - 1)))
        print(eta_fil, eta_ac, eta_hop, V_tunnel)
        # print(eta_ac - eta_fil + eta_hop - V_tunnel)
        return eta_fil, eta_ac, eta_hop, V_tunnel, current

    def tafel(self, y, V_app, apprx):
        eta_fil, eta_ac, eta_hop, V_tunnel, current_full = self.eta_fil(y, V_app, apprx)
        current = self.j_0et * (np.exp(np.float64(-self.alpha*const.elementary_charge*self.z/const.Boltzmann/self.T*eta_fil)) - 1.0)
        out = -self.M_me/self.z/const.elementary_charge/self.rho_me * current
        return out, eta_fil, eta_ac, eta_hop, V_tunnel, current_full

    def filament_growth(self, V_app, time, h):
        n_steps = np.int(time/h)

        y = np.zeros(n_steps+1, dtype=np.float64)
        y_s = np.zeros(n_steps+1, dtype=np.float64)
        eta_fil = np.zeros(n_steps+1, dtype=np.float64)
        eta_ac = np.zeros(n_steps+1, dtype=np.float64)
        eta_hop = np.zeros(n_steps+1, dtype=np.float64)
        V_tunnel = np.zeros(n_steps+1, dtype=np.float64)
        current = np.zeros(n_steps+1, dtype=np.float64)

        # plt.title("Beamscan DOA Estimation with a ULA. Spatial Spectrum")
        # plt.xlabel("Broadside Angle (degrees)")
        # plt.ylabel("Magnitude")
        y[0] = self.L
        for i in range(0, n_steps):
            # print(i)
            yi, eta_fil[i], eta_ac[i], eta_hop[i], V_tunnel[i], current[i] = self.tafel(y[i], V_app, [eta_fil[i], eta_ac[i], eta_hop[i], V_tunnel[i]])
            y_s[i+1] = y[i] + h * yi
            ysi1, _, _, _, _, _ = self.tafel(y_s[i+1], V_app, [eta_fil[i], eta_ac[i], eta_hop[i], V_tunnel[i]])
            y[i+1] = y[i] + h/2 * (yi + ysi1)

        plt.figure(0)
        plt.plot(np.asarray(range(0, n_steps+1)), y)
        plt.axis([0, n_steps, np.min(y), np.max(y)])
        plt.figure(1)
        plt.plot(np.asarray(range(0, n_steps+1)), -eta_fil)
        plt.axis([0, n_steps*10, 0, 0.16])
        plt.figure(2)
        plt.plot(np.asarray(range(0, n_steps+1)), current)
        plt.axis([0, n_steps, np.min(current), np.max(current)])
        plt.show()

    def t_nuc(self, eta):
        t_nuc = self.t_0nuc * np.exp(self.DeltaG_nuc/const.Boltzmann/self.T) * np.exp(-(self.N_c + self.alpha)*self.z * const.elementary_charge * eta/const.Boltzmann/self.T)
        # print(t_nuc)
        return t_nuc


if __name__ == "__main__":
    # app = Qt.QApplication(sys.argv)

    cell = Cell()
    V_app = 0.15
    eta_fil, eta_ac, eta_hop, V_tunnel, current = cell.eta_fil(cell.L, V_app)
    print(V_tunnel.dtype)
    print(cell.t_nuc(V_tunnel))
    cell.filament_growth(V_app, 5e-4, 1e-6)

    # sys.exit(app.exec_())
