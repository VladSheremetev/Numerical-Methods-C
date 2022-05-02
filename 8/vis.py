import matplotlib.pyplot as plt
def main():
	f = open('y.txt','r')
	l = []
	for line in f:
    		l.append(list(map(float, line.split())))
	plt.plot(l[0], l[1], label='euler')
	plt.plot(l[0], l[2], label='euler_predict')
	plt.plot(l[0], l[3], label='runge2')
	plt.plot(l[0], l[4], label='runge4')
	plt.plot(l[0], l[5], label='adams3')
	plt.title("y(x)")
	plt.grid()
	plt.legend()
	plt.show()
main()
