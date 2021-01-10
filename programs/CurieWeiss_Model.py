import numpy as np
from Bloch_Model import GraphPlot
from scipy.optimize import curve_fit, newton
from matplotlib import pyplot as plt

class CurieWeiss_fitting():
	def __init__(self):
		self.normalisedMVals    = []
		self.normalisedMValErrs = []

		self.Temperatures = (10,100,200,300,400,450,500,550,600,650,700,750,800)

	#function in the form from m = B(m*c) into B(m*c)-m = 0
	def ZeroFunction(self,m,c):
		return (1/np.tanh(m*c)) - 1/(m*c) - m

	def CurieWeiss(self,T,Tc):
		#initialise rootVal and array
		rootVal    = 1
		foundRoots = []

		for Temp in T:
			#if the temperature is equal to zero, this will result in a MathException. Whereas it should be impossible to get less than 0K. In either of these cases, use root of 1
			if Temp > 0:
				#value calculated for CurieWeiss equation used in self.ZeroFunction above
				const = 3*Tc/Temp

				#attempt to use NR root solver to find the root of the equation. Upon fail, use value of 0
				try:
					if rootVal != 0:
						foundRoots.append(newton(self.ZeroFunction,rootVal,args=(const,),maxiter=10000,rtol=1e-6))
				except:
					foundRoots.append(0)
			else:
				foundRoots.append(1)

		return foundRoots

	#function for obtaining the Ms0 value
	def getMs0Val(self,DIR):
		#reads from the chosen DIR and splits into an array based on spaces
		with open(DIR,'r') as f:
			line = f.readlines()[0].split()
			#line = line.split()

			#finds the individual vals for the Ms0 and its uncertainty respectfully and returns them
			return (int(line[1]),int(line[3]))

	def getNormalisedMVals(self,Ms0DIR,MsValsDIR='../plots/MsTCurve/data.txt'):
		#ensures that the arrays are empty before appending to them in the case of the function being used more than once for the object
		self.normalisedMVals    = []
		self.normalisedMValErrs = []

		#reads from the chosen DIR and obtains the Ms values
		Ms0,Ms0Err = self.getMs0Val(Ms0DIR)
		with open(MsValsDIR,'r') as f:
			for line in f:
				line = line.split()

				#normalised M = Ms/Ms0
				Ms 	  = int(line[3])
				MsErr = int(line[5])

				MNormalised = Ms/Ms0
				MNormError  = MNormalised * np.sqrt(((Ms0Err/Ms0)**2)+((MsErr/Ms)**2))

				self.normalisedMVals.append(MNormalised)
				self.normalisedMValErrs.append(MNormError)

	def graphCW(self,Tc,filename='data'):
		#create Graph() object to fit points to CW function
		#xvals = Temperature, yvals = normalised mag
		graph = GraphPlot(self.Temperatures,self.normalisedMVals,dy=self.normalisedMValErrs,title=f'Curie Weiss Law using {filename}',x_axis='Temperature(K)',y_axis='Normalised Magnetisation',saveDIR=f'../plots/CurieWeiss/{filename}.png')

		#fit graph for CurieWeiss Law
		popt,perr = graph.plot_curve_fit(self.CurieWeiss,bounds=(400,Tc),save=True)

		#save data obtained from fit
		with open(f'../plots/CurieWeiss/{filename}.txt','w') as file:
			file.write(f'Tc= {popt[0]} +/- {perr[0]}')


def main():
	#initialises the object
	CurieWeissFit = CurieWeiss_fitting()

	#creates the plot based from the none fitted values for B
	CurieWeissFit.getNormalisedMVals('../plots/MsTCurve/Bnofit.txt')
	CurieWeissFit.graphCW(836,'data for B=0.333')

	#tests values in the range of the error found for Tc=854+/-4K
	templist = [850,851,852,853,854,855,856,857,858]
	for i in templist:
		CurieWeissFit.graphCW(i,f'data from alpha fitting Tc={i}K')

	#creates the plot based from the fitted values for B=1/3
	CurieWeissFit.getNormalisedMVals('../plots/MsTCurve/Bfit.txt')
	#CurieWeissFit.graphCW(854,'data from a fitting Ms0=7.99e5')
	CurieWeissFit.graphCW(839,'data for B=0.340')

if __name__ == '__main__':
	main()