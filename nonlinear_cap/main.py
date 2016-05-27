import numpy as np
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
import matplotlib.pyplot as plt
from PyQt4 import Qt
import sys
from PyQt4 import QtGui
from PyQt4 import QtCore
from PyQt4.QtCore import QObject, pyqtSignal, pyqtSlot
import random


class Window(QtGui.QDialog):
    def __init__(self, parent=None):
        super(Window, self).__init__(parent)

        self.v_end = 1.
        self.v_start = -1.
        self.alpha1 = 0.0002
        self.alpha2 = 0.0001
        self.alpha_slider_coeff = 100000
        self.alpha_slider_coeff = 100000
        self.c10 = 1.
        self.c20 = 1.
        # a figure instance to plot on
        self.figure = plt.figure()

        # this is the Canvas Widget that displays the `figure`
        # it takes the `figure` instance as a parameter to __init__
        self.canvas = FigureCanvas(self.figure)

        # this is the Navigation widget
        # it takes the Canvas widget and a parent
        self.toolbar = NavigationToolbar(self.canvas, self)

        # Just some button connected to `plot` method
        self.button = QtGui.QPushButton('Plot')
        self.button.clicked.connect(self.plot)

        # Add slider
        self.alpha1_lcd = QtGui.QLabel("alpha1 = " + str(self.alpha1))
        self.alpha1_slider = QtGui.QSlider(QtCore.Qt.Horizontal, self)
        self.alpha1_slider.setRange(-100, 100)
        self.alpha1_slider.setSliderPosition(float(self.alpha1 * self.alpha_slider_coeff))
        # self.alpha1_lcd.display(float(self.alpha1 * 100))
        self.connect(self.alpha1_slider,  QtCore.SIGNAL('valueChanged(int)'),
                     self, QtCore.SLOT('update_alpha1_label(int)'))
        self.connect(self.alpha1_slider,  QtCore.SIGNAL('valueChanged(int)'), self.change_alpha1)
        alpha1_slider_layout = QtGui.QVBoxLayout()
        alpha1_slider_layout.addWidget(self.alpha1_lcd)
        alpha1_slider_layout.addWidget(self.alpha1_slider)

        self.alpha2_lcd = QtGui.QLabel("alpha2 = " + str(self.alpha2))
        self.alpha2_slider = QtGui.QSlider(QtCore.Qt.Horizontal, self)
        self.alpha2_slider.setRange(-100, 100)
        self.alpha2_slider.setSliderPosition(self.alpha2 * self.alpha_slider_coeff)
        self.connect(self.alpha2_slider,  QtCore.SIGNAL('valueChanged(int)'),
                     self, QtCore.SLOT('update_alpha2_label(int)'))
        self.connect(self.alpha2_slider,  QtCore.SIGNAL('valueChanged(int)'), self.change_alpha2)
        alpha2_slider_layout = QtGui.QVBoxLayout()
        alpha2_slider_layout.addWidget(self.alpha2_lcd)
        alpha2_slider_layout.addWidget(self.alpha2_slider)

        self.c10_lcd = QtGui.QLabel("c10 = " + str(self.c10))
        self.c10_slider = QtGui.QSlider(QtCore.Qt.Horizontal, self)
        self.c10_slider.setSliderPosition(float(self.c10 * 10))
        self.connect(self.c10_slider,  QtCore.SIGNAL('valueChanged(int)'),
                     self, QtCore.SLOT('update_c10_label(int)'))
        self.connect(self.c10_slider,  QtCore.SIGNAL('valueChanged(int)'), self.change_c10)
        c10_slider_layout = QtGui.QVBoxLayout()
        c10_slider_layout.addWidget(self.c10_lcd)
        c10_slider_layout.addWidget(self.c10_slider)

        self.c20_lcd = QtGui.QLabel("c20 = " + str(self.c20))
        self.c20_slider = QtGui.QSlider(QtCore.Qt.Horizontal, self)
        self.c20_slider.setSliderPosition(float(self.c20 * 10))
        self.connect(self.c20_slider,  QtCore.SIGNAL('valueChanged(int)'),
                     self, QtCore.SLOT('update_c20_label(int)'))
        self.connect(self.c20_slider,  QtCore.SIGNAL('valueChanged(int)'), self.change_c20)
        c20_slider_layout = QtGui.QVBoxLayout()
        c20_slider_layout.addWidget(self.c20_lcd)
        c20_slider_layout.addWidget(self.c20_slider)

        slider1_layout = QtGui.QHBoxLayout()
        slider1_layout.addLayout(alpha1_slider_layout)
        slider1_layout.addLayout(alpha2_slider_layout)

        slider2_layout = QtGui.QHBoxLayout()
        slider2_layout.addLayout(c10_slider_layout)
        slider2_layout.addLayout(c20_slider_layout)

        slider_layout = QtGui.QVBoxLayout()
        slider_layout.addLayout(slider1_layout)
        slider_layout.addLayout(slider2_layout)

        self.status = QtGui.QLabel("Status: OK")
        self.console = QtGui.QLabel("console")
        # set the layout
        layout = QtGui.QVBoxLayout()
        layout.addWidget(self.status)
        layout.addWidget(self.console)
        layout.addWidget(self.toolbar)
        layout.addWidget(self.canvas)
        layout.addWidget(self.button)
        layout.addLayout(slider_layout)
        self.setLayout(layout)

    @pyqtSlot(int)
    def update_alpha1_label(self, value):
        self.alpha1_lcd.setText("alpha1 = " + str(value/self.alpha_slider_coeff))

    @pyqtSlot(int)
    def update_alpha2_label(self, value):
        self.alpha2_lcd.setText("alpha2 = " + str(value/self.alpha_slider_coeff))

    @pyqtSlot(int)
    def update_c10_label(self, value):
        self.c10_lcd.setText("c10 = " + str(value/10.))

    @pyqtSlot(int)
    def update_c20_label(self, value):
        self.c20_lcd.setText("c20 = " + str(value/10.))

    def change_v_end(self, v_end):
        self.v_end = float(v_end/10)
        self.plot()

    def change_v_start(self, v_start):
        self.v_start = float(v_start/10)
        self.plot()

    def change_alpha1(self, alpha1):
        self.alpha1 = float(alpha1/100)
        self.plot()

    def change_alpha2(self, alpha2):
        self.alpha2 = float(alpha2/100)
        self.plot()

    def change_c10(self, c10):
        self.c10 = float(c10/10)
        self.plot()

    def change_c20(self, c20):
        self.c20 = float(c20/10)
        self.plot()

    def plot(self):
        # Все здесь написанное проверено только для дефолтных значений напряжений!!!!
        # print("Plot!")
        self.status.setText("Status: OK")
        v_start = self.v_start
        v_end = self.v_end
        alpha1 = self.alpha1
        alpha2 = self.alpha2
        c10 = self.c10
        c20 = self.c20

        # c1 = c10 * (alpha1 * v1*v1 + 1)
        # c2 = c20 * (alpha2 * v2*v2 + 1)
        # v = v1 + v2
        # if v_end == 0:
        #     print("V_end = 0")
        #     return
        # c1 * v1 = c2 * v2 = q => c = q/v
        # c = c10 * (alpha1 * v1**2 + 1) * v1 / v
        # c1(v1) * v1 = c2(v - v1) * (v - v1)
        # c10 * (alpha1 * v1**2 + 1) * v1 = c20 * (alpha2 * (v-v1)**2 + 1) * (v - v1)
        # c10 * alpha1 * v1**3 + c10 * v1 = c20 * alpha2 * (v**3 - 3*v**2*v1 + 3*v*v1**2 - v1**3) + c20*v - c20*v1
        #  0 = c20 * alpha2 * (v**3) + c20*v
        # (c10 * alpha1 + c20 * alpha2) * v1**3 + (- c20 * alpha2 * 3 * v) * v1**2 + (c10 + 3*c20 * alpha2 * v**2 + c20) * v1 + (-c20 * alpha2 * (v**3) - c20*v)
        # a3 = c10 * alpha1 + c20 * alpha2
        # a2 = - c20 * alpha2 * 3 * v
        # a1 = c10 + 3*c20 * alpha2 * v**2 + c20
        # a0 = -c20 * alpha2 * v**3 - c20 * v
        num_v_dots = 100
        a = np.zeros((4, num_v_dots,))
        roots = np.zeros((num_v_dots,), dtype=complex)
        c = np.zeros((num_v_dots,))

        for i in range(int(v_start*(num_v_dots/(v_end - v_start))), int(v_end*(num_v_dots/(v_end - v_start)))):
            c[i] = None
            roots[i] = 666j
            if i == 0:
                continue
            v = i * (v_end - v_start) / float(num_v_dots)
            a3 = c10 * alpha1 + c20 * alpha2
            a2 = -c20 * alpha2 * 3 * v
            a1 = c10 + 3*c20 * alpha2 * v**2 + c20
            a0 = -c20 * alpha2 * v**3 - c20 * v
            a[0, i] = a3
            a[1, i] = a2
            a[2, i] = a1
            a[3, i] = a0
            roots_tmp = np.roots(a[:,  i])
            c_all = (c10 * (alpha1 * np.real(roots_tmp) * np.real(roots_tmp) + 1) * np.real(roots_tmp) / float(v))
            for k in range(roots_tmp.shape[0]):
                if np.abs(np.imag(roots_tmp[k])) < 1e-10:
                    c_temp = (c10 * (alpha1 * np.real(roots_tmp[k]) * np.real(roots_tmp[k]) + 1) * np.real(roots_tmp[k]) / float(v))
                    if c_temp < 0 or np.real(roots_tmp[k]) < 0 or np.real(roots_tmp[k]) > v_end * i / float(num_v_dots):
                        pass
                    elif np.isnan(c[i]):
                        c[i] = c_temp
                        self.console.setText("console: roots_tmp = " + str(roots_tmp) + " c[" + str(i) + "] = " + str(c_all) )
                    else:
                        self.status.setText("Status: Uncertainity")
                        self.console.setText("console: roots_tmp = " + str(roots_tmp) + " c[" + str(i) + "] = " + str(c_all) )
                        # c[i] = None
                        # print ("Achtung!!")
                        break

            #c_tmp = (c10 * (alpha1 * np.real(roots_tmp) * np.real(roots_tmp) + 1) * np.real(roots_tmp) / float(v))

            #roots[i] = np.roots(a[:,  i])[2]
            #c[i] = (c10 * (alpha1 * np.real(roots[i]) * np.real(roots[i]) + 1) * np.real(roots[i]) / float(v))

        # create an axis
        ax = self.figure.add_subplot(111)
        # discards the old graph
        ax.hold(False)
        # plot data
        ax.plot(c)
        # refresh canvas
        self.canvas.draw()


if __name__ == "__main__":
    app = Qt.QApplication(sys.argv)
    main = Window()
    main.show()
    sys.exit(app.exec_())
