import csv
import sys
import math
from scipy import integrate
from scipy.stats import linregress
import numpy as np
import matplotlib.pyplot as plt

# Here is our parameters:
# D is the diameter of the telescope and F is the magnification of the telescope.
# S represents the integral of human eye sensitivity.
# T represents the transmittance which we assume it to be a gaussian distribution in visible spectrum.
# I represents the intensity of light under low light conditions which we assumed to be constant.
# Q represents quantum efficiency which we think it has a Poisson distribution in visible spectrum.
# SNR represents the signal-to-noise ratio index which we believe follows a linear distribution in visible spectrum.

class mytool1:

    # Sensitivity and wavelength relationship data path.
    data1_path = 'E:\pycharm location\pythonProject\code\CIE_sle_scotopic.csv'

    # Define initial variable.
    def __init__(self, diameter, transmittance, luminance):
        self.D = diameter
        self.F = transmittance
        self.I = luminance

    # The pupil diameter is defined, and the pupil distance and pupil size under dark light are
    # determined. It is believed that the pupil is generally 4-6mm.
    def pupil_diameter(self):
        z = self.D/self.F
        # return z
        if z>=4 and z<=7:
            print('Reasonable diameter and magnification')
            return z
        else:
            print('Unreasonable diameter and magnification')
            sys.exit()

    # read csv data.
    def read_data_csv(self):
        with open(self.data1_path, "r") as f:
            data = csv.reader(f)
            data_ls =list(data)
        return data_ls

    # Calculate human eye sensitivity.
    def calc_S(self):
        data1 = np.array(self.read_data_csv())
        x1 = np.array(data1[:, 0], dtype='float')
        y1 = np.array(data1[:, 1], dtype='float')
        S = integrate.trapz(y1, x1)/10000
        return S, y1

    # The integral function of transmittance is calculated and assumed to satisfy.
    # the Gaussian distribution in the visible region.
    def calc_T(self):
        # Define Gaussian function.
        def gaus_distribution(x):
            u = 535
            s = 100
            return math.exp(-0.5 * ((x - u) / s) ** 2) / (s * math.sqrt(2 * math.pi))
        x2 = np.linspace(400, 760, 361)
        y2 = []
        for i in range(x2.shape[0]):
            y2.append(gaus_distribution(x2[i]))
        T = integrate.trapz(y2, x2)
        return T,y2

    # The final evaluation function.
    def calc_V(self):
        V1 = math.sqrt(self.D*self.F)*self.calc_S()[0]*self.pupil_diameter()*self.calc_T()[0]*self.I
        print('Final evaluate function result:', V1)
    def calc_V1(self):
        V1 = math.sqrt(self.D*self.F)*self.calc_S()[0]*self.pupil_diameter()*self.calc_T()[0]*self.I
        return V1

    # Plot the Gaussian distribution transmission and human eye sensitivity.
    def paint(self):
        x1 = np.linspace(380, 780, 401)
        x2 = np.linspace(400, 760, 361)
        y1 = self.calc_S()[1]
        y2 = self.calc_T()[1]
        plt.subplot(2, 1, 1)
        plt.plot(x1, y1, 'r-')
        plt.ylabel("Eye sensitivity", fontsize=30)
        plt.title("Diameter "+str(self.D)+" magnification "+str(self.F), fontsize=30)
        plt.subplot(2, 1, 2)
        plt.plot(x2, np.array(y2)*1000, 'b-')
        plt.xlabel("Visible wavelength", fontsize=30)
        plt.ylabel("Mixed telescope transmittance*1000", fontsize=30)
        plt.show()

# Plot the corrected twilight factor for 9 points with pupil distance of 7.
def paint_z_7():
    y_7 = []
    x_7 = []
    for i in range(1, 10):
        a_7 = i * 7
        x_7.append(a_7)
        a_7 = mytool1(a_7, i, 50).calc_V1()
        y_7.append(a_7)
    plt.plot(x_7, y_7, marker='o')
    plt.grid(axis='x')
    plt.xticks(size=20)
    plt.yticks(size=20)
    plt.title("Exit pupil diameter is 7", fontsize=30)
    plt.xlabel("Lens diameter(mm)", fontsize=30)
    plt.ylabel("Twilight vision factor", fontsize=30)
    plt.show()

# Plot the conjugate diameter and magnification of the corrected twilight coefficient.
def paint_conjugate():
    x = []
    y = []
    m = mytool1(56, 8, 50)
    x.append(56)
    y.append(m.calc_V1())
    n = mytool1(8, 56, 50)
    x.append(8)
    y.append(n.calc_V1())
    m = mytool1(50, 10, 50)
    x.append(50)
    y.append(m.calc_V1())
    n = mytool1(10, 50, 50)
    x.append(10)
    y.append(n.calc_V1())
    m = mytool1(60, 12, 50)
    x.append(60)
    y.append(m.calc_V1())
    n = mytool1(12, 60, 50)
    x.append(12)
    y.append(n.calc_V1())
    n = mytool1(80, 20, 50)
    x.append(80)
    y.append(n.calc_V1())
    n = mytool1(20, 80, 50)
    x.append(20)
    y.append(n.calc_V1())
    n = mytool1(25, 8, 50)
    x.append(25)
    y.append(n.calc_V1())
    n = mytool1(8, 25, 50)
    x.append(8)
    y.append(n.calc_V1())
    sizes = np.array([600, 150, 600, 150, 600, 150, 600, 150, 600, 150])
    colors = np.array(["red", "red", "green", "green", "black", "black",
                       "blue", "blue", "purple", "purple"])
    plt.scatter(x, y, s=sizes, c=colors)
    plt.xticks(size=20)
    plt.yticks(size=20)
    plt.title("Conjugate diameter and magnification", fontsize=30)
    plt.xlabel("Lens diameter", fontsize=30)
    plt.ylabel("Twilight vision factor", fontsize=30)
    plt.show()

class mytool2():

    # SNR data path.
    data2_path = 'E:\pycharm location\pythonProject\code\SNR_DATA.csv'

    def __init__(self, diameter, transmittance, luminance):
        self.D = diameter
        self.F = transmittance
        self.I = luminance

    # The integral function of transmittance is calculated and assumed to satisfy the Gaussian
    # distribution in the visible region.
    def calc_T(self):
        # Define Gaussian function.
        def gaus_distribution(x):
            u = 535
            s = 100
            return math.exp(-0.5 * ((x - u) / s) ** 2) / (s * math.sqrt(2 * math.pi))
        x2 = np.linspace(400, 760, 361)
        y2 = []
        for i in range(x2.shape[0]):
            y2.append(gaus_distribution(x2[i]))
        T = integrate.trapz(y2, x2)
        return T, y2

    # Read the SNR data and calculate the fit results.
    def SNR_paras(self):
        # Calculate the slope and intercept of the SNR line fitted.
        with open(self.data2_path, "r") as f:
            data2 = csv.reader(f)
            data_ls = list(data2)
        data2 = np.array(np.array(data_ls)[1:, 1:], dtype='float')
        y = np.squeeze(data2[:, 0])
        x = np.squeeze(data2[:, 1])
        slope, intercept, r, p, e = linregress(x, y)
        return slope, intercept

    # The SNR integral function is calculated and assumed to be a linear distribution。
    def calc_SNR(self):
        # Defined linear function
        def liner_func(x):
            slope, intercept = self.SNR_paras()
            return slope*x+intercept
        x = np.linspace(400, 760, 361)
        y = []
        for i in range(x.shape[0]):
            y.append(liner_func(x[i]))
        SNR = integrate.trapz(y, x)
        return SNR, y, x

    # t1 and t2 represent the Poisson distribution bound, and l1 and l2 are the target intervals.
    def calc_Q(self):
        l1 = 400
        l2 = 760
        l = 0.25
        def calc_Poisson(j1, j2, i):
            t1 = 300
            t2 = 1500
            # The default sampling interval between points is 10 because of computational power issues.
            m = int((t2-t1)/10)
            n = int((j2-j1)/10)
            x = []
            y = []
            for k in range(0, n+1):
                k = int(k+j1/10-t1/10)
                p = (i*m)**k/math.factorial(k)*math.exp(-1*i*m)
                x.append(k*10+t1)
                y.append(p)
            return x, y
        #  When the first two parameters here are the same as t2 and t1, you should get the whole interval.
        a, b = calc_Poisson(l1, l2, l)
        return integrate.trapz(b, a), b, a

    # The final CMOS evaluation function.
    def calc_V(self):
        V2 = math.sqrt(self.D*self.F*self.calc_Q()[0]*self.calc_T()[0]*self.calc_SNR()[0])*self.I**2/10**6
        print('Final CMOS evaluate function result:', V2)
    def calc_V2(self):
        V2 = math.sqrt(self.D*self.F*self.calc_Q()[0]*self.calc_T()[0]*self.calc_SNR()[0])*self.I**2/10**6
        return V2

    # Plot SNR signal-to-noise ratio distribution and quantum efficiency distribution in the visible spectrum。
    def plot_Poisson_SNR(self):
        x1 = np.linspace(400, 760, 361)
        y1 = self.calc_SNR()[1]
        plt.plot(x1, y1, 'r-')
        plt.xticks(size=20)
        plt.yticks(size=20)
        plt.title("SNR linear distribution", fontsize=30)
        plt.xlabel("Visble wavelength", fontsize=30)
        plt.ylabel("Probability", fontsize=30)
        plt.show()
        x2 = np.linspace(400, 760, 37)
        y2 = self.calc_Q()[1]
        plt.plot(x2, y2, 'r-')
        plt.xticks(size=20)
        plt.yticks(size=20)
        plt.title("Quantum efficiency distribution", fontsize=30)
        plt.xlabel("Visble wavelength", fontsize=30)
        plt.ylabel("Probability", fontsize=30)
        plt.show()

# Plot cmos evluaiton with differnent illuminance.
def plot_cmos():
    x = np.linspace(0.25, 1000, 40)
    y = []
    D = 56
    F = 8
    for i in x:
        m = mytool2(D, F, i)
        y.append(m.calc_V2())
    plt.plot(x, y, 'r-')
    plt.xticks(size=20)
    plt.yticks(size=20)
    plt.title("Diameter "+str(D)+" magnification "+str(F), fontsize=30)
    plt.xlabel("Visble wavelength(mm)", fontsize=30)
    plt.ylabel("CMOS twilight vision factor", fontsize=30)
    plt.show()