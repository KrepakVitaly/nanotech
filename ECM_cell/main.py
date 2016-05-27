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
from sympy import Symbol
from sympy import exp, sinh, ln, asinh


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
        self.N_c        = 3  # n/a
        self.A_ac       = 804.25e-18  # m2
        self.A_fil      = 12.57e-18  # m2
        self.A_is       = 12.57e-18  # m2
        self.L          = 20e-9  # m
        self.rho_fil    = 1.7e-8  # Ohm m
        self.R_el       = 76.4  # Ohm
        self.R_S        = 1e6  # Ohm

        self.C = 2.7
        self.T = 300  # K

    def eta_fil(self, x, V_app):
        DeltaW_0 = 1 * const.elementary_charge  # Joule
        # print (const.elementary_charge)
        m_eff = self.m_r * const.electron_mass

        A = self.C * 3 * np.sqrt(2 * m_eff * DeltaW_0) / 2 / x * (const.elementary_charge / const.Planck)**2 * \
            np.exp(- 4 * const.pi * x / const.Planck * np.sqrt(2 * m_eff * DeltaW_0)) * self.A_fil
        print(A)

        B = self.j_0et * self.A_fil
        print(B)

        C = self.R_el + self.R_S + self.rho_fil * (self.L - x) / self.A_fil
        print(C)

        eta = 0

        eta_ac = (const.Boltzmann * self.T / (1-self.alpha) / const.elementary_charge / self.z * ln(self.A_fil/self.A_ac*(exp(- self.alpha * const.elementary_charge * self.z / const.Boltzmann / self.T * eta) - 1) + 1))
        print(eta_ac)

        eta_hop = (x*2*const.Boltzmann*self.T/self.a/self.z/const.elementary_charge*asinh(self.j_0et/self.j_0hop*(exp(- self.alpha * const.elementary_charge * self.z / const.Boltzmann / self.T * eta) - 1)))
        print(eta_hop)

        # g_eta_fil = \
        #     lambda eta: np.exp(- self.alpha * const.elementary_charge * self.z / const.Boltzmann / self.T * eta) - 1

        # eta = Symbol('eta')
        # eq = np.exp(- self.alpha * const.elementary_charge * self.z / const.Boltzmann / self.T * eta) - 1

        # eq = (A*C + 1)*(const.Boltzmann * self.T / (1-self.alpha) / const.elementary_charge / self.z * ln(self.A_fil/self.A_ac*(exp(- self.alpha * const.elementary_charge * self.z / const.Boltzmann / self.T * eta) - 1) + 1)) - \
        #      (A*C + 1)*eta + \
        #      (A*C + 1)*(x*2*const.Boltzmann*self.T/self.a/self.z/const.elementary_charge*asinh(self.j_0et/self.j_0hop*(exp(- self.alpha * const.elementary_charge * self.z / const.Boltzmann / self.T * eta) - 1))) +\
        #      B*(exp(- self.alpha * const.elementary_charge * self.z / const.Boltzmann / self.T * eta) - 1) - V_app

        # print(solve(eq, eta, implicit=True))


if __name__ == "__main__":
    # app = Qt.QApplication(sys.argv)

    cell = Cell()
    cell.eta_fil(cell.L, 0.125)

    # sys.exit(app.exec_())
