from enum import Enum
import matplotlib.pyplot as plt
import numpy as np



sub_conc = [50, 100, 200, 400, 600, 800, 1000]
initial_v = [1, 1.687, 2.57, 3.48, 3.95, 4.23, 4.4]

plt.scatter(sub_conc, initial_v)
plt.xlabel('[S]')
plt.ylabel('Vi')
plt.title('Uninhibited Enzyme Velocities')
plt.show()

# Convert the lists to numpy arrays
sub_conc_array = np.array(sub_conc)
initial_v_array = np.array(initial_v)

one_over_sub_conc = 1 / sub_conc_array
one_over_initial_v = 1 / initial_v_array

# Perform linear regression
a, b = np.polyfit(one_over_sub_conc, one_over_initial_v, 1)

# Plot the enzyme data
plt.scatter(one_over_sub_conc, one_over_initial_v, label='Uninhibited Velocity')
plt.plot(one_over_sub_conc, a * one_over_sub_conc + b, label=f'Line of Best Fit (Vmax={1/b:.2f}, Km={a/b:.2f})', color='red')
plt.xlabel('1 / Substrate Concentration')
plt.ylabel('1 / Initial Velocity')
plt.title('Uninhibited Enzyme Lineweaver-Burk')
plt.legend()
plt.grid(True)
plt.show()

def find_vmax_km(sub_conc, initial_v, apparent): #apparent is for inhib and slices 0th value as that is inhib conc
    initial_v = initial_v[1:] if apparent else initial_v
   
    sub_conc_reciprocal = [1/i for i in sub_conc]
    initial_v_reciprocal = [1/i for i in initial_v]
   
    slope = (initial_v_reciprocal[-1] - initial_v_reciprocal[0]) / (sub_conc_reciprocal[-1] - sub_conc_reciprocal[0])
   
    y_intercept = initial_v_reciprocal[0] - (slope * sub_conc_reciprocal[0])
   
    vmax = 1/y_intercept if y_intercept != 0 else 0
    km = slope * vmax
   
    return (vmax, km, slope, y_intercept)

enzyme_vmax, enzyme_km, enzyme_slope, enzyme_intercept = find_vmax_km(sub_conc, initial_v, False)

inhibitor_data = [[10, 0.93, 1.5, 2.16, 2.79, 3.05, 3.22, 3.33], 
                  [25, 0.84375, 1.285714, 1.741935, 2.117647, 2.2816901, 2.373626, 2.4324324],
                  [50, 0.72972972, 1.038461, 1.3170731, 1.52112676, 1.603960, 1.648854, 1.677018],
                  [100, 0.574468, 0.75, 0.885245, 0.972972, 1.00621, 1.02369, 1.03448],
                  [150, 0.473684, 0.58695, 0.666666, 0.71523, 0.733031, 0.742268, 0.747922]] # use excel .csv and convert to pandas DF next time

inhib_concs = [i[0] for i in inhibitor_data]
app_inhib_constants = []

# Perform linear regression and plot lines of best fit for inhibitor data
for i in range(len(inhibitor_data)):
    inhib_array = np.array(inhibitor_data[i][1:])
    inhib_array = 1 / inhib_array

    a, b = np.polyfit(one_over_sub_conc, inhib_array, 1)

    app_inhib_constants.append((1/b, a/b, a, b)) # vmax, km, slope, y-intercept

    # Plot the line of best fit for the inhibitor data
    plt.plot(one_over_sub_conc, a * one_over_sub_conc + b, label=f'Inhibitor Conc: {inhibitor_data[i][0]} (Vmax={1/b:.2f}, Km={a/b:.2f})')

plt.xlabel('1 / Substrate Concentration')
plt.ylabel('1 / Initial Velocity')
plt.title('Inhibitor(s) Lineweaver-Burk Plot')
plt.legend()
plt.grid(True)
plt.show()

one_over_app_vmax = [1/i[0] for i in app_inhib_constants]
one_over_app_km = [1/i[1] for i in app_inhib_constants]
inhib_slopes = [i[2] for i in app_inhib_constants]

# positive delta's mean they increased.
app_vmax = [i[0] for i in app_inhib_constants]
app_km = [i[1] for i in app_inhib_constants]

vmax_delta_avg = 0
km_delta_avg = 0

for i in range(len(app_vmax) - 1):
  vmax_delta = round(app_vmax[i+1], 2) - round(app_vmax[i], 2)
  km_delta = round(app_km[i+1], 2) - round(app_km[i], 2)
  
  vmax_delta_avg += vmax_delta
  km_delta_avg += km_delta

vmax_delta_avg /= len(app_vmax) - 1
km_delta_avg /= len(app_km) - 1


class InhibitorTypes(Enum):
    COMPETITIVE = 1
    NON_COMPETITIVE = 2
    UNCOMPETITIVE = 3
    MIXED = 4


#use delta's to determine inhib type based on how Vmax and Km change w/ diff inhibs
myInhibitorType = None

if vmax_delta_avg == 0 and km_delta_avg > 0:
    myInhibitorType = InhibitorTypes.COMPETITIVE
elif vmax_delta_avg < 0 and km_delta_avg == 0:
    myInhibitorType = InhibitorTypes.NON_COMPETITIVE
elif vmax_delta_avg < 0 and km_delta_avg < 0:
    myInhibitorType = InhibitorTypes.UNCOMPETITIVE
else:
    myInhibitorType = InhibitorTypes.MIXED
  

#Ki & Ki' plots

def plotKiPrimeGraph(one_over_app_vmax, inhib_concs):
  one_over_app_vmax_arr, inhib_concs_arr = np.array(one_over_app_vmax), np.array(inhib_concs)
  a, b = np.polyfit(inhib_concs_arr, one_over_app_vmax_arr, 1)

  
  plt.scatter(inhib_concs_arr, one_over_app_vmax_arr, label='Inhib Data')
  plt.plot(inhib_concs_arr, a * inhib_concs_arr + b, label=f'Line of Best Fit (Ki Prime={b/a:.2f})', color='red')
  plt.xlabel('Inhib Concentration')
  plt.ylabel('1 / Apparent Vmax')
  plt.title("Ki' Plot")
  plt.legend()
  plt.grid(True)
  plt.show()

  return b/a


def plotKiGraph(slopes, inhib_concs):
  slopes, inhib_concs_arr = np.array(inhib_slopes), np.array(inhib_concs)
  a, b = np.polyfit(inhib_concs_arr, slopes, 1)

  
  plt.scatter(inhib_concs_arr, slopes, label='Inhib Data')
  plt.plot(inhib_concs_arr, a * inhib_concs_arr + b, label=f'Line of Best Fit (Ki={b/a:.2f})', color='red')
  plt.xlabel('Inhib Concentration')
  plt.ylabel('Slope')
  plt.title("Ki Plot")
  plt.legend()
  plt.grid(True)
  plt.show()

  return b/a

if myInhibitorType == InhibitorTypes.COMPETITIVE or myInhibitorType == InhibitorTypes.NON_COMPETITIVE or myInhibitorType == InhibitorTypes.MIXED:
    kiValue = plotKiGraph(inhib_slopes, inhib_concs)
    print(f'The Ki for my {myInhibitorType.name} Inhibitor is: {kiValue:.2f}')
if myInhibitorType == InhibitorTypes.UNCOMPETITIVE or myInhibitorType == InhibitorTypes.MIXED:
    kiPrimeValue = plotKiPrimeGraph(one_over_app_vmax, inhib_concs)
    print(f"The Ki' (Ki Prime) for my {myInhibitorType.name} Inhibitor is: {kiPrimeValue:.2f}")
