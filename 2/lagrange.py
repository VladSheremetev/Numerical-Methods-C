import matplotlib.pyplot as plt
import pandas as pd


df_func = pd.read_csv('C:\\Users\\Владислав\\func.csv')
df_lagrang = pd.read_csv('C:\\Users\\Владислав\\lagrange.csv')
df_knots = pd.read_csv('C:\\Users\\Владислав\\knots.csv')
size_array = len(df_lagrang['x_coordinate'])
plt.xlim(df_lagrang['x_coordinate'][0], df_lagrang['x_coordinate'][size_array - 1])
plt.plot(df_lagrang['x_coordinate'], df_lagrang['y_coordinate'], "r", label="L(x)", linewidth=2)
plt.plot(df_func['x_coordinate'], df_func['y_coordinate'], "b", label="F(x)")
plt.plot(df_knots['x_coordinate'], df_knots['y_coordinate'], 'co', alpha=1, markersize=4)
plt.grid(1)
plt.legend()
plt.show()
