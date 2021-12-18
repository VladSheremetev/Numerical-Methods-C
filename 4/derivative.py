import matplotlib.pyplot as plt
import pandas as pd


df_func = pd.read_csv('C:\\Users\\Владислав\\analytical_derivative.csv')
df_lagrang = pd.read_csv('C:\\Users\\Владислав\\derivative_step_h.csv')
df_knots = pd.read_csv('C:\\Users\\Владислав\\derivative_step_half.csv')
df_runge = pd.read_csv('C:\\Users\\Владислав\\derivative_runge.csv')
size_array = len(df_lagrang['x_coordinate'])
plt.xlim(df_lagrang['x_coordinate'][0], df_lagrang['x_coordinate'][size_array - 1])
plt.plot(df_lagrang['x_coordinate'], df_lagrang['y_coordinate'], "g", label="oneh", linewidth = 3)
plt.plot(df_knots['x_coordinate'], df_knots['y_coordinate'], "r", label="halfh", linewidth = 3)
plt.plot(df_runge['x_coordinate'], df_runge['y_coordinate'], "c", label="runge", linewidth = 3)
plt.plot(df_func['x_coordinate'], df_func['y_coordinate'], "b", label="dy/dx")
plt.grid(1)
plt.legend()
plt.show()
