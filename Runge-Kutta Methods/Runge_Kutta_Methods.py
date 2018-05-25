import numpy as np
from scipy. integrate import odeint
import matplotlib.pyplot as plt

class DifferentialEquation:

    def __init__(self, a, b, h, y0):
        self.__begin = a
        self.__end = b
        self.__step = h
        self.__n = (b - a) // h
        self.__y0 = y0
        self.__x = np.linspace(a, b, self.__n + 1)
        self.__func = lambda y, x: x + y

    def EulersMethod(self):
        y = [self.__y0]
        i = 1
        while i <= self.__n:
            y.append(y[-1] + self.__step * self.__func(y[-1], self.__x[i - 1]))
            i += 1
        return y
    
    def SecondOrderMethod(self):
        y = [self.__y0]
        i = 0
        while i <= self.__n:
            y.append(y[-1] + self.__step 
                     * self.__func(y[-1] + self.__step / 2 * self.__func(y[-1], self.__x[i - 1]),
                                   self.__x[i-1] + self.__step / 2)) 
            i += 1
        return y
    
    def RK4(self):
        y = [self.__y0]
        i = 0
        while i <= self.__n:
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
        plt.style.use("bmh")
        plt.xlim([self.__begin, self.__end])
        plt.plot(self.__x, EulersMethod(), color = 'tomato')
        plt.plot(self.__x, SecondOrderMethod(), color = 'blue')
        plt.plot(self.__x, RK4(), color = 'black')
        plt.plot(self.__x, RealValue(), color = 'pink')
        plt.plot(self.__begin, self.__y0, 'g^')
        plt.show()

def Main():
    ude = DifferentialEquation(a, b, h, y0)
    ude.Graph()
    pass

Main()