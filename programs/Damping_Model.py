from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import numpy as np
from Bloch_Model import GraphPlot

class Damping_fitting():
	def __init__(self):
		self.constant = (2.002319*1.602177e-19/(2*9.109304e-31))*1e-7
		self.x  = []
		self.y  = []
		self.dy = []

	def Alpha(self,x,alpha0,Tc):
		return alpha0*(1-(x/(Tc*3)))

	#function for returning a string of the entered values rounded to an appropriate number of sig figs based on the uncertainty provided
	def sf(self,val,error):
		round_to = 1 - len(str(int(error)))
		return f'{int(round(val,round_to))} +/- {int(round(error,round_to))}'

	#special function for adjusting for different numbers for rounding
	def sfalpha(self,val,error):
		round_to = 5
		naughts = '0.0000'

		#if this is too high of a rounding value, increase by one
		if round(error*10**(round_to)) == 0:
			round_to += 1
			naughts += '0'

		return f'{round(val,round_to)} +/- {naughts}{round(error*10**(round_to))}'

	#function for loading the data and creating datapoints for the M against T graph
	def create_points(self,DIR='../plots/data.txt'):
		with open(DIR,'r') as f:

			#empty arrays for temp storage for finding mean values
			AlphaWiVals = []
			WiVals 		= []

			#itterates through every line in the given file
			for line in f:

				line = line.split()

				#data for frequency is stored in GHz, converted to Hz via *e9
				freq = int(line[1][2:-3])*1e9
				Dh   = int(line[2][3:])
				DhE  = int(line[4])

				#calcualtes the value for damping (alpha) at frequency using damning equation formula
				Alpha 		= self.constant * (Dh/freq)
				AlphaError	= Alpha * DhE/Dh

				Wi = 1/(AlphaError**2)
				WiVals.append(Wi)
				AlphaWiVals.append(Alpha*Wi)
				

				#if all 4 Ms values have been obtained for the 4 values for frequency
				if len(AlphaWiVals) == 4:
					#weighted mean calculation
					aMean    = sum(AlphaWiVals)/sum(WiVals)
					aMeanErr = np.sqrt(1/sum(WiVals))

					#obtain value for Temperature of the data set
					T = int((line[0])[2:-1])
					
					#appends to the datapoints and resets the temp arrays
					self.x.append(T)
					self.y.append(aMean)
					self.dy.append(aMeanErr)
					AlphaWiVals = []
					WiVals      = []

	def graphAlpha(self):
		#create Graph() object to fit points to Bloch function
		graph = GraphPlot(self.x,self.y,dy=self.dy,title='Damping function fit',x_axis='T(K)',y_axis='Damping',saveDIR='../plots/Damping Values/Damping plot.png')

		#fit graph to Bloch function with fitted value for Beta
		popt,perr = graph.plot_curve_fit(self.Alpha,[0.0201,2500],save=True)

		#save data obtained from fit
		with open('../plots/Damping Values/data.txt','w') as data:
			for i in range(len(self.x)):
				data.write(f'T={self.x[i]}K\ta={self.sfalpha(self.y[i],self.dy[i])}\n')
		#save data obtained from fit
		with open('../plots/Damping Values/fitted values.txt','w') as out:
			out.write(f'a0={round(popt[0],6)} +/- 0.00000{int(round(perr[0]*1e6))}\tTc={self.sf(popt[1],perr[1])}')

def main():
	alpha = Damping_fitting()
	alpha.create_points()
	alpha.graphAlpha()

if __name__ == '__main__':
	from Bloch_Model import GraphPlot
	main()