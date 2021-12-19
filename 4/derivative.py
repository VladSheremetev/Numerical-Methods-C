import matplotlib.pyplot as plt
import pandas as pd


df_func = pd.read_csv('C:\\Users\\Владислав\\analytical_derivative.csv')
df_lagrang = pd.read_csv('C:\\Users\\Владислав\\derivative_step_h.csv')
df_knots = pd.read_csv('C:\\Users\\Владислав\\derivative_step_half.csv')
df_runge = pd.read_csv('C:\\Users\\Владислав\\derivative_runge.csv')
size_array = len(df_lagrang['x_coordinate'])
plt.xlim(df_lagrang['x_coordinate'][0], df_lagrang['x_coordinate'][size_array - 1])
plt.plot(df_lagrang['x_coordinate'], df_lagrang['y_coordinate'], "go", label="oneh", linewidth=0.5, markersize=4)
plt.plot(df_knots['x_coordinate'], df_knots['y_coordinate'], "ro", label="halfh", linewidth=0.5, markersize=4)
plt.plot(df_runge['x_coordinate'], df_runge['y_coordinate'], "co", label="runge", linewidth=0.8, markersize=4)
plt.plot(df_func['x_coordinate'], df_func['y_coordinate'], "b", label="dy/dx")
plt.grid(1)
plt.legend()
plt.show()

