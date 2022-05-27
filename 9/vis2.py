import matplotlib.pyplot as plt
def main():
	f = open('nonlinear.txt','r')
	l = []
	for line in f:
    		l.append(list(map(float, line.split())))
	plt.plot(l[0], l[1], label='newton')
	plt.plot(l[0], l[2], label='chordes')
	plt.title("Nolinear ODU")
	plt.grid()
	plt.legend()
	plt.show()
main()
