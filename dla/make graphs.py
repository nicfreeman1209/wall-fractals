import os
import csv
#import Tkinter as tk
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def make_graph(filename, data):
	x = []
	y = []
	for line in data:
		n = float(line[0])
		r = float(line[3])
		x.append(n)
		y.append(r**1.75/n)
	ax = plt.axes()
	ax.set_axisbelow(False)
	plt.grid(axis='y', which='both')
	plt.plot(x, y)
	plt.xlabel("n")
	plt.ylabel("r^1.75/n")
	plt.savefig(filename + '.png', dpi=150)	
	plt.clf()
	
for root,dirs,files in os.walk("./"):
	for filename in files:
		if filename.endswith(".csv"):
			print(filename)
			f = open(filename, 'r')
			for i in range(3):
				f.readline() #ignore first three lines
			data = csv.reader(f, delimiter=",")
			make_graph(filename, data)
			f.close()