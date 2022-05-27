import matplotlib.pyplot as plt
def main():
	f = open('linear.txt','r')
	l = []
	for line in f:
    		l.append(list(map(float, line.split())))
	plt.plot(l[0], l[1], label='ballistic')
	plt.plot(l[0], l[2], label='finite_substrations')
	plt.title("Linear ODU")
	plt.grid()
	plt.legend()
	plt.show()
main()
