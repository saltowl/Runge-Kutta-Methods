import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

class DifferentialEquation:

    def __init__(self, a, b, h, y0):
        self.__begin = a
        self.__end = b
        self.__step = h
        self.__y0 = y0

        self.__x = [a]
        i = 1
        while self.__x[-1] != b:
            self.__x.append(a + i * h)
            i += 1

        self.__func = lambda y, x: (2 - x**2 - y**2) / (2 + x**2 + x * y)

    def EulersMethod(self):
        y = [self.__y0]
        i = 1
        while self.__x[i-1] < self.__end:
            y.append(y[-1] + self.__step * self.__func(y[-1], self.__x[i - 1]))
            i += 1
        return y
    
    def SecondOrderMethod(self):
        y = [self.__y0]
        i = 1
        while self.__x[i-1] < self.__end:
            y.append(y[-1] + self.__step 
                     * self.__func(y[-1] + self.__step / 2 * self.__func(y[-1], self.__x[i - 1]),
                                   self.__x[i-1] + self.__step / 2)) 
            i += 1
        return y
    
    def RK4(self):
        y = [self.__y0]
        i = 1
        while self.__x[i-1] < self.__end:
            k1 = self.__step * self.__func(y[-1], self.__x[i - 1])
            k2 = self.__step * self.__func(y[-1] + k1 / 2, 
                                           self.__x[i - 1] + self.__step / 2)
            k3 = self.__step * self.__func(y[-1] + k2 / 2,
                                           self.__x[i - 1] + self.__step / 2)
            k4 = self.__step * self.__func(y[-1] + k3,
                                           self.__x[i - 1] + self.__step)
            y.append(y[-1] + (k1 + 2*k2 + 2*k3 + k4) / 6)
            i += 1
        return y

    def RealValue(self):
        return odeint(self.__func, self.__y0, self.__x)


    def Graph(self):
        plt.style.use("ggplot")
        plt.plot(self.__x, self.EulersMethod(), color = 'tomato', 
                 label = 'Method Runge-Kutta first order (Euler\'s Method)')
        plt.plot(self.__x, self.SecondOrderMethod(), color = 'royalblue', linestyle = '--', 
                 label = 'Method Runge-Kutta second order')
        plt.plot(self.__x, self.RK4(), color = 'black', linestyle = '-.', 
                 label = 'Method Runge-Kutta fourth order')
        plt.plot(self.__x, self.RealValue(), color = 'sandybrown', linestyle = '-', 
                 label = 'Real solution of equation')
        plt.plot(self.__begin, self.__y0, 'g^')
        plt.xlabel('X')
        plt.ylabel('Y')
        plt.title('A particular solution of equation')
        plt.legend()
        plt.show()

def Main():
    a = 0
    b = 1
    h = 0.1
    y0 = 0
    ude = DifferentialEquation(a, b, h, y0)
    ude.Graph()
    pass

Main()