from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import numpy as np

class GraphPlot():
	def __init__(self,x_points,y_points,dx=None,dy=None,title=None,x_axis=None,y_axis=None,saveDIR=None):
		#Graph info
		self.title   = title
		self.xlabel  = x_axis
		self.ylabel  = y_axis
		#location to save graphs to
		self.saveDIR = saveDIR

		#plot points. Should be in array format
		self.x_points = x_points
		self.y_points = y_points
		self.x_uncert = dx
		self.y_uncert = dy

	#function for plotting graph
	def plot(self,save=False,show=False):
		plt.errorbar(self.x_points,self.y_points,xerr=self.x_uncert,yerr=self.y_uncert,fmt='bo')

		#fills axis names and title if provided
		if self.title != None:
			plt.title(self.title)
		if self.ylabel != None:
			plt.ylabel(self.ylabel)
		if self.xlabel != None:
			plt.xlabel(self.xlabel)

		#draws to graph
		plt.draw()

		#optional input for showing the plot
		if show:
			plt.show()

		#if set, will save the graph image as the title in png format if no filename is given in the DIR. Otherwise, default filename is "graphplot.png"
		if save:
			if self.saveDIR != None:
				if self.saveDIR.endswith('.png'):
					plt.savefig(self.saveDIR)
				elif self.title != None:
					plt.savefig(self.saveDIR + self.title + '.png')
				else:
					plt.savefig(self.saveDIR + 'graphplot.png')
			else:
				print('No directory has been set for the graph images')

	#function for plotting graph with fitting
	def plot_curve_fit(self,fitting_function,fitting_vals=None,bounds=None,save=False,show=False):
		#fits the plot points to the function given
		if bounds != None:
			popt,pcov = curve_fit(fitting_function,xdata=self.x_points,ydata=self.y_points,bounds=bounds)
		else:
			popt,pcov = curve_fit(fitting_function,xdata=self.x_points,ydata=self.y_points,p0=fitting_vals)
		#calculates chi squared uncertainty values
		perr = np.sqrt(np.diag(pcov))

		#plot based from best fit with standard plot points
		self.plot()
		plt.plot(self.x_points,fitting_function(self.x_points,*popt))
		plt.draw()

		#optional input for showing the plot
		if show:
			plt.show()

		#if set, will save the graph image as the title in png format if no filename is given in the DIR. Otherwise, default filename is "graphplot.png"
		if save:
			if self.saveDIR != None:
				if self.saveDIR.endswith('.png'):
					plt.savefig(self.saveDIR)
				elif self.title != None:
					plt.savefig(self.saveDIR + self.title + '.png')
				else:
					plt.savefig(self.saveDIR + 'graphplot.png')
			else:
				print('No directory has been set for the graph images')

		plt.clf()

		return popt,perr


class Bloch_fitting():
	def __init__(self):
		self.constant   = ((9.109304e-31)**2)/((1e-14)*(2.002319**2)*(1.602177e-19**2))
		self.x  = []
		self.y  = []
		self.dy = []

	#### Fitting functions:
	#function of Bloch fit with fitted B value
	def Bloch(self,x,Ms0,Tc,B):
		return Ms0*((1-(x/Tc))**B)
	#function of Bloch fit with fixed value of B=1/3
	def Bloch_fixed_B(self,x,Ms0,Tc):
		return self.Bloch(x,Ms0,Tc,1/3)

	#quick function for finding the number of significant figures based on the error and returning a string of the rounded results
	def sf(self,val,error):
		round_to = 1 - len(str(int(error)))
		return f'{int(round(val,round_to))} +/- {int(round(error,round_to))}'

	#function for loading the data and creating datapoints for the M against T graph
	def create_points(self,DIR='../plots/data.txt'):
		with open(DIR,'r') as f:

			#empty arrays for temp storage for finding mean values
			MsWiVals = []
			WiVals   = []

			#itterates through every line in the given file
			for line in f:

				line = line.split() # i=1 for f, i=5 for h0, i=7 for h0error

				#data for frequency is stored in GHz, converted to Hz via *e9
				freq = int(line[1][2:-3])*1e9
				h0   = int(line[5][3:])
				h0E  = int(line[7])

				#calcualtes the value for Ms at that frequency using Kittels formula in terms of Ms
				Ms = self.constant*(freq**2/h0) - h0

				MsErrorSqrd = (self.constant * h0E*(freq**2)/(h0**2))**2 + h0E**2
				Wi = 1/MsErrorSqrd
				WiVals.append(Wi)
				MsWiVals.append(Ms*Wi)
				

				#if all 4 Ms values have been obtained for the 4 values for frequency
				if len(MsWiVals) == 4:
					#weighted mean calculation
					MsMean    = sum(MsWiVals)/sum(WiVals)
					MsMeanErr = np.sqrt(1/sum(WiVals))

					#obtain value for Temperature of the data set
					T = int((line[0])[2:-1])
					
					#appends to the datapoints and resets the temp arrays
					self.x.append(T)
					self.y.append(MsMean)
					self.dy.append(MsMeanErr)
					MsWiVals = []
					WiVals   = []

	def graphBloch(self):
		#create Graph() object to fit points to Bloch function
		graph = GraphPlot(self.x,self.y,dy=self.dy,title='Bloch function fit (B fitted)',x_axis='T(K)',y_axis='Ms(Am-1)',saveDIR='../plots/MsTCurve/Bloch_w_fitted_B.png')

		#saves the values for the Ms and T
		with open('../plots/MsTCurve/data.txt','w') as dat:
			for i in range(len(self.x)):
				dat.write(f'T= {self.x[i]}\tMs= {self.sf(self.y[i],self.dy[i])}\n')

		#fit graph to Bloch function with fitted value for Beta
		popt,perr = graph.plot_curve_fit(self.Bloch,[800000,900,1/3],save=True)

		#save data obtained from fit
		with open('../plots/MsTCurve/Bfit.txt','w') as Bfit:
			Bfit.write(f'Ms0= {int(round(popt[0],-3))} +/- {int(round(perr[0],-3))}\t\tTc= {int(popt[1])} +/- {int(perr[1])}\t\tB= {round(popt[2],3)} +/- {round(perr[2],3)}')

		graph.title   = 'Bloch function fit (Fixed B=1/3)'
		graph.saveDIR = '../plots/MsTCurve/Bloch_w_fixed_B.png'

		#fit graph to Bloch function with fixed value for Beta
		popt,perr = graph.plot_curve_fit(self.Bloch_fixed_B,[800000,900],save=True)

		#save data obtained from fit
		with open('../plots/MsTCurve/Bnofit.txt','w') as Bnofit:
			Bnofit.write(f'Ms0= {int(round(popt[0],-3))} +/- {int(round(perr[0],-3))}\t\tTc= {int(popt[1])} +/- {int(perr[1])}')

def main():
	bloch = Bloch_fitting()
	bloch.create_points()
	bloch.graphBloch()

if __name__ == '__main__':
	main()