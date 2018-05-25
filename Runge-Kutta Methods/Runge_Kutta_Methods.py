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
		plt.plot(self.__begin, self.__y0, 'g^', label = 'The initial condition')
		plt.xlabel('X')
		plt.ylabel('Y')
		plt.title('A particular solution of equation')
		plt.legend()
		plt.show()


	def PrintResult(self):
		RK1 = self.EulersMethod()
		RK2 = self.SecondOrderMethod()
		RK4 = self.RK4()
		real = self.RealValue()

		print("_____________________________________________________________________________________")
		print('|     | Method Runge-Kutta | Method Runge-Kutta | Method Runge-Kutta |               |')
		print('|  X  |    first order     |     second order   |     fourth order   |   Real Value  |')
		print("|_____|____________________|____________________|____________________|_______________|")

		for i in range(len(self.__x)):
			print('|{0:^5}|{1:^20}|{2:^20}|{3:^20}|{4:^15}|'
				  .format(round(self.__x[i], 2), round(RK1[i], 9), round(RK2[i], 9), 
						  round(RK4[i], 9), round(real[i][0], 9)))
		
		print('|_____|____________________|____________________|____________________|_______________|')



def GetInputData():
	while True:
		bounds = input("Enter a space between the beginning and the end of the gap: ").split(' ')
		if (len(bounds) == 2):
			try:
				begin = float(bounds[0])
				end = float(bounds[1])
			except (ValueError):
				print("Please enter valid numbers! Try again.\n")
				continue
		elif (bounds[0] == ''):
			begin = 0
			end = 1
		else:
			print("The gap must be given by two real numbers!\n")
			continue
		break

	while True:
		stepString = input("Enter a step to specify a uniform grid: ")
		if (len(stepString) != 0):
			try:
				step = float(stepString)
			except (ValueError):
				print("Enter the decimals! Try again.\n")
				continue
		else:
			step = 0.1
		break

	while True:
		sourceString = input("Enter the initial condition for the solution of the Cauchy problem - y0(x): ")
		if (len(sourceString) != 0):
			try:
				source = float(sourceString)
			except (ValueError):
				print("Enter the decimals! Try again.\n")
				continue
		else:
			source = 0
		break
	
	return (begin, end, step, source)


def Main():
	data = GetInputData()
	a = data[0]
	b = data[1]
	h = data[2]
	y0 = data[3]
	print("Solution of the Cauchy problem on the interval [{0}; {1}] in increments of {2} with y({0}) = {3} ..."
	   .format(a, b, h, y0))

	ude = DifferentialEquation(a, b, h, y0)
	ude.PrintResult()
	ude.Graph()
	pass

Main()