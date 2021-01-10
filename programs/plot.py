from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
import numpy as np

#function for plotting best fit to lorentzian curve
def lorentzian(x,y0,s,dx,x0):
	return y0 + (s*dx)/((4*((x-x0)**2))+dx**2)

#quick function for finding the number of significant figures based on the error and returning a string of the rounded results
def sf(val,error):
	round_to = 1 - len(str(int(error)))
	return f'{int(round(val,round_to))} +/- {int(round(error,round_to))}'

def main():
	#temperature range used
	temp_range = (10,100,200,300,400,450,500,550,600,650,700,750,800)

	#due to the method used for appending data to the file, this is needed to ensure the file is empty beforehand
	with open('../plots/data.txt','w') as file:
		file.write('')

	#nested for loop iterating through the temperature ranges used as well as the frequency of the voltage for the file name
	for temp in temp_range:
		#frequency from 10GHz to 40GHz for each temperature
		for freq in range(10,50,10):
			x = []
			y = []
			with open(f'../data/fmr_permalloy_thinfilm_T_K_{temp}_f_GHz_{freq}.0.txt','r') as data:
				for line in data:
					line = line.split()
					x.append(int(line[0]))
					y.append(round(float(line[1])**2,4))

			#estimations used for fitting
			peak = x[y.index(max(y))]
			dh = 6000
			s = 1e8*dh

			#fits the values to the curve and finds the appropriate errors
			popt,pcov = curve_fit(lorentzian,x,y,[min(y),s,dh,peak])
			perr = np.sqrt(np.diag(pcov))
			
			#plots the new plot
			plt.plot(x,y,'bo')
			plt.title(f'T={temp}K, f={freq}GHz')
			plt.ylabel('FMR Signal Squared')
			plt.xlabel('Magnetic Field (H)')
			plt.plot(x,lorentzian(x,*popt))
			plt.draw()
			plt.savefig(f'../plots/T={temp}K,f={freq}GHz.png')

			plt.clf()

			#saves the plot data
			with open('../plots/data.txt','a') as f:
				f.write(f'T={temp}K\tf={freq}GHz\t\t\tdH={sf(popt[2],perr[2])}\t\tH0={sf(popt[3],perr[3])}\n')

if __name__ == '__main__':
	main()