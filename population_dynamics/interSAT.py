''' ****************************************************************** '''
''' ****************************************************************** '''
'''              INTERPOP, 2 layers, q = 1 : varying <k>               '''
''' ****************************************************************** '''
''' ****************************************************************** '''

# Main 1: Higher-order case yields cooling spinodal thresholds at q=1; 
# Main 2: Higher-order case smoothens out the abrupt jump for k_av < 3; 
# Main 3: Thermal and higher-order cases match hypergraph results for q=1;
# Main 4: Higher-order catches heating threshold 4-spin hypergraph for q = 1;
# Main 5: Higher-order does not work good for q<1 for the cooling.

import matplotlib.pyplot as plt
plt.style.use('classic')

import numpy as np
import random as rm
import time


''' ************************************************************************************* '''
''' ----------------- Thermal algorithm, segregated time scales 2 layers ---------------- '''
''' ************************************************************************************* '''

# Kernel functions solving the cavity fields via saddle-point
def kernel2sp(x, y):
    return (y / 2) * np.log((np.cosh((x + 1.0) / y)) / (np.cosh((x - 1.0) / y)))

# Zero-truncated Poisson distribution (for sim's only with the GCC)
import scipy.stats as sct

def truncated_Poisson(mu, min_value, size):
    temp_size = size
    while True:
        temp_size *= 2
        temp = sct.poisson.rvs(mu, size=temp_size)
        truncated = temp[temp >= min_value]
        if len(truncated) >= size:
            return truncated[:size]


# ***************** Multilayer population settings *******************

k_av = 4 
N_pop = 10**3 # population size

# Indices identifying independent and interdependent populations
net_indices1 = list(range(0, N_pop))
net_indices2 = rm.sample(net_indices1, N_pop) # randomly interdependent
int_pairs = [x for x in zip(net_indices1, net_indices2)] #interdependent pairs

# Initializing lists of degrees and nearest neighbors
deg_list1, deg_list2, nn_list1, nn_list2 = [], [], [], []

for i in range(N_pop):
    
    # Network 1 ******************************
    deg_list1.append(np.random.poisson(k_av, 1)[0])
    nn_list1.append(rm.sample(np.delete(net_indices1, i).tolist(), deg_list1[i])) 
    
    # Network 2 ******************************
    deg_list2.append(np.random.poisson(k_av, 1)[0])
    nn_list2.append(rm.sample(np.delete(net_indices2, i).tolist(), deg_list2[i]))  


# Initial conditions: random, then by annealed heating/cooling
T = 10**3 # popdyn convergence
NOIs = 10**2 # interdependent iterations
cav_array1 = k_av * np.random.random_sample((N_pop, 1))
cav_array2 = k_av * np.random.random_sample((N_pop, 1))
loc_mag1 = np.full(N_pop, 1.0)
loc_mag2 = np.full(N_pop, 1.0)


# ******************* HEATING PART **************************
start_time = time.time()

magNOI1_th_loc_heat, magNOI2_th_loc_heat = [], [] 
temps2sp = np.linspace(0.01 * k_av, 0.6 * k_av, 10**2)
err_heat = 1e-15

for temp in temps2sp:
    
    print('The heating temperatures is:', temp)
    # Generating two arrays of local temperatures
    temp_eff1 = np.full(N_pop, temp)
    temp_eff2 = np.full(N_pop, temp)
    
    for run in range(NOIs):     
        
        # ****************** NETWORK 1 *********************
        
        # Local temperatures: network 1
        for i in range(N_pop):
            
            int_i2 = int_pairs[i][1]
            if loc_mag2[int_i2] >= 1e-10:
                temp_eff1[i] = temp / loc_mag2[int_i2]
            else: 
                temp_eff1[i] = temp / (loc_mag2[int_i2] + 1e-10)
            
        # Population dynamics: network 1 
        for pditer in range(T):
            for i in range(N_pop):
                Sum1 = err_heat
                for index in nn_list1[i]:
                    Sum1 += kernel2sp(cav_array1[index][0], temp_eff1[i])
                cav_array1[i] = Sum1    
                
        # Calculating the average spins' magnetization          
        mag_array1 = []
        for i in range(N_pop):
            a = cav_array1[i]
            mag_array1.append(np.tanh(a / temp_eff1[i]))
           
        # Calculating the local magnetization of spins' neighborhoods
        for i in range(N_pop):
            loc_sum1 = 0
            for index in nn_list1[i]:
                loc_sum1 += mag_array1[index] / deg_list1[i]
            loc_mag1[i] = loc_sum1


        # ****************** NETWORK 2 *********************
        
        # Local temperatures: network 2
        for i in range(N_pop):
            
            ind2 = int_pairs[i][1]
            int_i1 = int_pairs[i][0]
            if loc_mag1[int_i1] >= 1e-10:
                temp_eff2[ind2] = temp / loc_mag1[int_i1]
            else: 
                temp_eff2[ind2] = temp / (loc_mag1[int_i1] + 1e-10)  
        
        # Population dynamics: network 2
        for pditer in range(T):                
            for i in range(N_pop):
                Sum2 = err_heat
                for index in nn_list2[i]:
                    Sum2 += kernel2sp(cav_array2[index][0], temp_eff2[i])
                cav_array2[i] = Sum2
                                
        # Calculating the average spins' magnetization          
        mag_array2 = []
        for i in range(N_pop):
            a = cav_array2[i]
            mag_array2.append(np.tanh(a / temp_eff2[i]))
           
        # Calculating the local magnetization of spins' neighborhoods
        for i in range(N_pop):
            loc_sum2 = 0
            for index in nn_list2[i]:
                loc_sum2 += mag_array2[index] / deg_list2[i]
            loc_mag2[i] = loc_sum2

    magNOI1_th_loc_heat.append(np.mean(mag_array1))
    magNOI2_th_loc_heat.append(np.mean(mag_array2))
    print(magNOI1_th_loc_heat[-1], magNOI2_th_loc_heat[-1])


# ******************* COOLING PART **************************

magNOI1_th_loc_cool, magNOI2_th_loc_cool = [], [] 
err_cool = 1e-15

for temp in temps2sp[::-1]:
    
    print('The cooling temperature is:', temp)
    # Generating two arrays of local temperatures
    temp_eff1 = np.full(N_pop, temp)
    temp_eff2 = np.full(N_pop, temp)
    
    for run in range(NOIs):
        
        
        # ****************** NETWORK 1 *********************
        
        # Local temperatures: network 1
        for i in range(N_pop):
            
            int_i2 = int_pairs[i][1]
            if loc_mag2[int_i2] >= 1e-10:
                temp_eff1[i] = temp / loc_mag2[int_i2]
            else: 
                temp_eff1[i] = temp / (loc_mag2[int_i2] + 1e-10)
            
        # Population dynamics: network 1 
        for pditer in range(T):
            for i in range(N_pop):
                Sum1 = err_cool
                for index in nn_list1[i]:
                    Sum1 += kernel2sp(cav_array1[index][0], temp_eff1[i])
                cav_array1[i] = Sum1    
                
        # Calculating the average spins' magnetization          
        mag_array1 = []
        for i in range(N_pop):
            a = cav_array1[i]
            mag_array1.append(np.tanh(a / temp_eff1[i]))
           
        # Calculating the local magnetization of spins' neighborhoods
        for i in range(N_pop):
            loc_sum1 = 0
            for index in nn_list1[i]:
                loc_sum1 += mag_array1[index] / deg_list1[i]
            loc_mag1[i] = loc_sum1


        # ****************** NETWORK 2 *********************
        
        # Local temperatures: network 2
        for i in range(N_pop):
            
            ind2 = int_pairs[i][1]
            int_i1 = int_pairs[i][0]
            if loc_mag1[int_i1] >= 1e-10:
                temp_eff2[ind2] = temp / loc_mag1[int_i1]
            else: 
                temp_eff2[ind2] = temp / (loc_mag1[int_i1] + 1e-10)
        
        # Population dynamics: network 2
        for pditer in range(T):                
            for i in range(N_pop):
                Sum2 = err_cool
                for index in nn_list2[i]:
                    Sum2 += kernel2sp(cav_array2[index][0], temp_eff2[i])
                cav_array2[i] = Sum2
                                
        # Calculating the average spins' magnetization          
        mag_array2 = []
        for i in range(N_pop):
            a = cav_array2[i]
            mag_array2.append(np.tanh(a / temp_eff2[i]))
           
        # Calculating the local magnetization of spins' neighborhoods
        for i in range(N_pop):
            loc_sum2 = 0
            for index in nn_list2[i]:
                loc_sum2 += mag_array2[index] / deg_list2[i]
            loc_mag2[i] = loc_sum2

    magNOI1_th_loc_cool.append(np.mean(mag_array1))
    magNOI2_th_loc_cool.append(np.mean(mag_array2))
    print(magNOI1_th_loc_cool[-1], magNOI2_th_loc_cool[-1])


# Checking for the heating and cooling thresholds
for i in range(5, len(temps2sp)): 
    if abs(magNOI1_th_loc_heat[i] - magNOI1_th_loc_heat[i - 1]) > 1e-1:
        break
Tsp_th_loc_heat = round((temps2sp[i] + temps2sp[i - 1]) / (2 * k_av), 2)

for i in range(5, len(temps2sp)): 
    if abs(magNOI1_th_loc_cool[i] - magNOI1_th_loc_cool[i - 1]) > 1e-1:
        break
Tsp_th_loc_cool = round((temps2sp[::-1][i] + temps2sp[::-1][i - 1]) / (2 * k_av), 2)


# Plotting the results - heating curves
plt.plot(temps2sp / k_av, magNOI1_th_loc_heat, '.', color = 'red')
plt.plot(temps2sp / k_av, magNOI2_th_loc_heat, '.', color = 'orange') 
# Plotting the results - cooling curves
plt.plot(temps2sp[::-1] / k_av, magNOI1_th_loc_cool, '.', color = 'red')
plt.plot(temps2sp[::-1] / k_av, magNOI2_th_loc_cool, '.', color = 'orange') 
plt.title('$k_{av} =$ %s, $T_{sp}^{heat} = $ %s, $T_{sp}^{cool} = $ %s' %(k_av, Tsp_th_loc_heat, Tsp_th_loc_cool))    
plt.xlim(0.0, (k_av + 0.5) / k_av)
plt.ylim(0.0, 1.0)
plt.show()

print("--- %s seconds ---" % (time.time() - start_time))



''' ****************************************************************************************** '''
''' ------------------  Thermal algorithm, mixed time scales: 2 ER graphs  ------------------- '''
''' ****************************************************************************************** '''
# Since it uses the equilibrium sol of the 2-spin model, it is equivalent to the algorithm above.
# However, in this form, the algorithm is O(N)-faster than the thermal one since here NOI = 1.

import matplotlib.pyplot as plt
plt.style.use('classic')

import numpy as np
import random as rm
import time


# Kernel functions solving the cavity fields via saddle-point
def kernel2sp(x, y):
    return (y / 2) * np.log((np.cosh((x + 1.0) / y)) / (np.cosh((x - 1.0) / y)))


# ***************** MULTILAYER POPULATION SETTINGS *******************

k_av = 4 
N_pop = 10**3 # population size

# Indices identifying independent and interdependent populations
net_indices1 = list(range(0, N_pop))
net_indices2 = rm.sample(net_indices1, N_pop) # randomly interdependent
int_pairs = [x for x in zip(net_indices1, net_indices2)] #interdependent pairs

# Initializing lists of degrees and nearest neighbors
deg_list1, deg_list2, nn_list1, nn_list2 = [], [], [], []

for i in range(N_pop):
    
    # Network 1 ******************************
    deg_list1.append(np.random.poisson(k_av, 1)[0])
    nn_list1.append(rm.sample(np.delete(net_indices1, i).tolist(), deg_list1[i])) 
    
    # Network 2 ******************************
    deg_list2.append(np.random.poisson(k_av, 1)[0])
    nn_list2.append(rm.sample(np.delete(net_indices2, i).tolist(), deg_list2[i]))  


# Initial conditions: first step is random, then by annealed heating/cooling
T = 10**1 # popdyn convergence
cav_array1 = k_av * np.random.random_sample((N_pop, 1))
cav_array2 = k_av * np.random.random_sample((N_pop, 1))
temps2sp = np.linspace(0.01 * k_av, 0.6 * k_av, 10**2)
mag_list1, mag_list2 = np.zeros((N_pop, len(temps2sp))), np.zeros((N_pop, len(temps2sp)))
loc_mag1, loc_mag2 = np.full(N_pop, 1.0), np.full(N_pop, 1.0)


# ********************** HEATING PART ************************

start_time = time.time()
mag2spin1_heat, mag2spin2_heat = [], []
err_heat = 1e-15

t = -1
for temp in temps2sp:
    t += 1
    
    temp_eff1 = np.full(N_pop, temp)
    temp_eff2 = np.full(N_pop, temp)
    print('The heating temperature is:', temp)

    for pditer in range(T):
        
        # Local temperature distributions
        for i in range(N_pop):
            
            # Network 1 *******************
            int_i2 = int_pairs[i][1]
            if loc_mag2[int_i2] >= 1e-10:
                temp_eff1[i] = temp / loc_mag2[int_i2]
            else: 
                temp_eff1[i] = temp / (loc_mag2[int_i2] + 1e-10)
            
            # Network 2 *******************
            ind2 = int_pairs[i][1]
            int_i1 = int_pairs[i][0]
            if loc_mag1[int_i1] >= 1e-10:
                temp_eff2[ind2] = temp / loc_mag1[int_i1]
            else: 
                temp_eff2[ind2] = temp / (loc_mag1[int_i1] + 1e-10)   
        
        # Calculating instantaneous cavities and local magnetizations
        for i in range(N_pop):
            
            # Network 1 *******************
            Sum1 = err_heat
            for index in nn_list1[i]: 
                Sum1 += kernel2sp(cav_array1[index][0], temp_eff1[i])
            cav_array1[i] = Sum1
            mag_list1[i][t] = np.tanh(cav_array1[i] / temp_eff1[i])
            
            # Network 2 *******************
            Sum2 = err_heat
            for index in nn_list2[i]: 
                Sum2 += kernel2sp(cav_array2[index][0], temp_eff2[i])
            cav_array2[i] = Sum2
            mag_list2[i][t] = np.tanh(cav_array2[i] / temp_eff2[i])
        
        # Calculating the local OPs at this PD's iteration
        for i in range(N_pop):
            
            # Network 1 *******************
            loc_sum1 = 0
            for index in nn_list1[i]:
                loc_sum1 += mag_list1[index][t] / deg_list1[i]
            loc_mag1[i] = loc_sum1
            
            # Network 2 *******************
            loc_sum2 = 0
            for index in nn_list2[i]:
                loc_sum2 += mag_list2[index][t] / deg_list2[i]
            loc_mag2[i] = loc_sum2
            
    mag_array1, mag_array2 = [], []
    for i in range(N_pop):
        mag_array1.append(mag_list1[i][t])
        mag_array2.append(mag_list2[i][t])

    mag2spin1_heat.append(np.mean(mag_array1))
    mag2spin2_heat.append(np.mean(mag_array2))
    print(np.mean(mag_array1), np.mean(mag_array2))
    


# ********************** COOLING PART ************************

loc_mag1, loc_mag2 = np.full(N_pop, 1.0), np.full(N_pop, 1.0)
mag2spin1_cool, mag2spin2_cool = [], []
err_cool = 1e-15

t = -1
for temp in temps2sp[::-1]:
    t += 1
    
    temp_eff1 = np.full(N_pop, temp)
    temp_eff2 = np.full(N_pop, temp)
    
    # initiating cavities to negative values to catch the transcritical point
    cav_array1 = - k_av * np.random.random_sample((N_pop, 1))
    cav_array2 = - k_av * np.random.random_sample((N_pop, 1))
    print('The cooling temperature is:', temp)

    for pditer in range(T):
        
        # Local temperatures distributions 
        for i in range(N_pop):
            
            # Network 1 *******************
            int_i2 = int_pairs[i][1]
            if loc_mag2[int_i2] >= 1e-10:
                temp_eff1[i] = temp / loc_mag2[int_i2]
            else: 
                temp_eff1[i] = temp / (loc_mag2[int_i2] + 1e-10)
            
            # Network 2 *******************
            ind2 = int_pairs[i][1]
            int_i1 = int_pairs[i][0]
            if loc_mag1[int_i1] >= 1e-10:
                temp_eff2[ind2] = temp / loc_mag1[int_i1]
            else: 
                temp_eff2[ind2] = temp / (loc_mag1[int_i1] + 1e-10)     
        
        # Calculating instantaneous cavities and local magnetizations
        for i in range(N_pop):
            
            # Network 1 *******************
            Sum1 = err_cool
            for index in nn_list1[i]: 
                Sum1 += kernel2sp(cav_array1[index][0], temp_eff1[i])
            cav_array1[i] = Sum1
            mag_list1[i][t] = np.tanh(cav_array1[i] / temp_eff1[i])
            
            # Network 2 *******************
            Sum2 = err_cool
            for index in nn_list2[i]: 
                Sum2 += kernel2sp(cav_array2[index][0], temp_eff2[i])
            cav_array2[i] = Sum2
            mag_list2[i][t] = np.tanh(cav_array2[i] / temp_eff2[i])
        
        
        # Calculating the local OPs at this PD's iteration
        for i in range(N_pop):
            
             # Network 1 *******************
            loc_sum1 = 0
            for index in nn_list1[i]:
                loc_sum1 += mag_list1[index][t] / deg_list1[i]
            loc_mag1[i] = loc_sum1
            
            # Network 2 *******************
            loc_sum2 = 0
            for index in nn_list2[i]:
                loc_sum2 += mag_list2[index][t] / deg_list2[i]
            loc_mag2[i] = loc_sum2
            
    mag_array1, mag_array2 = [], []
    for i in range(N_pop):
        mag_array1.append(mag_list1[i][t])
        mag_array2.append(mag_list2[i][t])

    mag2spin1_cool.append(np.mean(mag_array1))
    mag2spin2_cool.append(np.mean(mag_array2))
    print(np.mean(mag_array1), np.mean(mag_array2))


# Checking for the heating threshold
for i in range(5, len(temps2sp)):
    if mag2spin1_heat[i] < 1/20 or mag2spin2_heat[i] < 1/20:
        break
Tsp_mix_loc_heat = round((temps2sp[i] + temps2sp[i - 1]) / (2 * k_av), 2)

# Checking for the cooling threshold 
for i in range(5, len(temps2sp)):
    if abs(mag2spin1_heat[i]) < 1/20 or abs(mag2spin2_heat[i]) < 1/20:
        break
Tsp_mix_loc_cool = round((temps2sp[::-1][i] + temps2sp[::-1][i - 1]) / (2 * k_av), 2)


# Plotting the results for the heating curves
plt.plot(temps2sp / k_av, mag2spin1_heat, '--', color = 'red', linewidth = 1.5) 
plt.plot(temps2sp / k_av, mag2spin2_heat, '--', color = 'black', linewidth = 1.5) 
# Ploting the results for the cooling curves
plt.plot(temps2sp[::-1] / k_av, mag2spin1_cool, '--', color = 'red', linewidth = 1.5) 
plt.plot(temps2sp[::-1] / k_av, mag2spin2_cool, '--', color = 'black', linewidth = 1.5) 
plt.title(' $T_{sp}^{heat} / k_{av} =$ %s, $T_{sp}^{cool} / k_{av} =$ %s' %(Tsp_mix_loc_heat, Tsp_mix_loc_cool))
plt.xlim(0.0, (k_av + 0.5) / k_av)
plt.ylim(0.0, 1.0)
plt.show()

print("--- %s seconds ---" % (time.time() - start_time))


''' ************************************************************************ '''
''' --------------- Hyphergraph : random 3-spin interactions --------------- '''
''' ************************************************************************ '''

# Integral kernel solving the saddle-point replica symmetric equation
def kernel3sp(x1, x2, y):
    alpha_p = np.cosh((x2 + x1) / y) + np.exp(- 2 / y) * np.cosh((x2 - x1) / y)
    alpha_m = np.cosh((x2 - x1) / y) + np.exp(- 2 / y) * np.cosh((x2 + x1) / y)
    return np.log(alpha_p / alpha_m)

# Function grouping list of elements in n-uples
def grouped_list(lst, n):
    if len(lst) == 0:
        return []
    else:
        return list(zip(*[iter(lst)] * n))


# **************** Structural part ******************
k_av = 4
N_pop = 10**3
net_indices = (np.arange(N_pop))
deg_list, nn_tri_list = [], []
T = 10**1
cav_array = k_av * np.random.random_sample((N_pop, 1))
temps3sp = np.linspace(0.01 * k_av, 0.6 * k_av, 10**2)

for i in range(N_pop): 
    deg_list.append(np.random.poisson(k_av, 1)[0])
    indices_i = rm.sample(np.delete(net_indices, i).tolist(), 2 * deg_list[i])
    nn_tri_list.append(grouped_list(indices_i, 2))


start_time = time.time()
# ****************** HEATING PART **********************
mag3sp_heat = []
for temp in temps3sp:
    for runs in range(T):
        for i in range(N_pop):

            Sum = 0
            for index in nn_tri_list[i]: 
                Sum += kernel3sp(cav_array[index[1]], cav_array[index[0]], temp)
            cav_array[i] = temp * Sum / 2
    
    mag_array3sp = []
    for i in range(N_pop):
        a = cav_array[i]
        mag_array3sp.append(np.tanh(a / temp))

    mag3sp_heat.append(np.mean(mag_array3sp))


# ****************** COOLING PART **********************
mag3sp_cool = []
for temp in temps3sp[::-1]:
    for runs in range(T):
        for i in range(N_pop):

            Sum = 0
            for index in nn_tri_list[i]: 
                Sum += kernel3sp(cav_array[index[1]], cav_array[index[0]], temp)
            cav_array[i] = temp * Sum / 2
    
    mag_array3sp = []
    for i in range(N_pop):
        a = cav_array[i]
        mag_array3sp.append(np.tanh(a / temp))

    mag3sp_cool.append(np.mean(mag_array3sp))


# Searching for the heating threshold
for i in range(5, len(temps3sp)): 
    if abs(mag3sp_heat[i] - mag3sp_heat[i - 1]) > 1e-1:
        break
Tsp_hy_heat = round((temps3sp[i] + temps3sp[i - 1]) / (2 * k_av), 2)

# Searching for the cooling threshold
for i in range(5, len(temps3sp)): 
    if abs(mag3sp_cool[i] - mag3sp_cool[i - 1]) > 1e-1:
        break
Tsp_hy_cool = round((temps3sp[::-1][i] + temps3sp[::-1][i - 1]) / (2 * k_av), 2)
    

# Plotting the results
plt.plot(temps3sp / k_av, mag3sp_heat, '-', color = 'k', label = '3-spin hypergraph')
plt.plot(temps3sp[::-1] / k_av, mag3sp_cool, '-', color = 'k')
plt.legend(loc = 'best')
plt.title(' $T_{sp}^{heat} / k_{av} =$ %s, $T_{sp}^{cool} / k_{av} =$ %s' %(Tsp_hy_heat, Tsp_hy_cool))
plt.ylim(0.0, 1.0)
plt.xlim(0.0, (k_av + 0.5) / k_av)
plt.show()

print("--- %s seconds ---" % (time.time() - start_time))



''' ************************************************************************ '''
''' ---------------------  M2 - synthetic fluctuations --------------------- '''
''' ************************************************************************ '''

import matplotlib.pyplot as plt
plt.style.use('classic')

import numpy as np
import random as rm
import time


# Kernel functions solving the cavity fields via saddle-point
def kernel2sp(x, y):
    return (y / 2) * np.log((np.cosh((x + 1.0) / y)) / (np.cosh((x - 1.0) / y)))


# ***************** MULTILAYER POPULATION SETTINGS *******************

N_pop = 10**3 # population size
T_eq = 10**2 # popdyn convergence
Tsp_kav_q1_fluct = []
loc_thr = 1e-10
err_heat = 1e-20
N_runs = 10**2


# Calling data from single runs to optimize interval of temperatures
data_th_fluc = np.loadtxt('/Users/physivan/Desktop/Tsp_kav_q1_M2_fluc_1run.txt')
k_avs = data_th_fluc[3:8, 0]
Tsp_kav_1run = data_th_fluc[3:8, 1]


start_time = time.time()

ik = -1
for k_av in k_avs: 
    
    ik += 1
    # Collecting statistics over runs for better numerical thresholds
    Tsp_kav_fluc = []
    
    for run in range(N_runs):
        # Indices identifying independent and interdependent populations
        net_indices1 = list(range(0, N_pop))
        net_indices2 = rm.sample(net_indices1, N_pop) # randomly interdependent
        int_pairs = [x for x in zip(net_indices1, net_indices2)] #interdependent pairs
        
        # Initializing lists of degrees and nearest neighbors
        deg_list1, deg_list2, nn_list1, nn_list2 = [], [], [], []
        
        for i in range(N_pop):
        
            # Network 1 ******************************
            deg_list1.append(np.random.poisson(k_av, 1)[0])
            nn_list1.append(rm.sample(np.delete(net_indices1, i).tolist(), deg_list1[i])) 
        
            # Network 2 ******************************
            deg_list2.append(np.random.poisson(k_av, 1)[0])
            nn_list2.append(rm.sample(np.delete(net_indices2, i).tolist(), deg_list2[i]))  
        
        # Initial conditions: first step is random, then by annealed heating/cooling
        temps2sp = np.linspace((Tsp_kav_1run[ik] - 0.1) * k_av, (Tsp_kav_1run[ik] + 0.1) * k_av, 20)
        cav_array1 = k_av * np.random.random_sample((N_pop, 1))
        cav_array2 = k_av * np.random.random_sample((N_pop, 1))
        mag_list1, mag_list2 = np.zeros((N_pop, len(temps2sp))), np.zeros((N_pop, len(temps2sp)))
        loc_mag1, loc_mag2 = np.full(N_pop, 1.0), np.full(N_pop, 1.0)
        
        
        # ********************** Annealing part ************************
        
        mag2spin1_heat, mag2spin2_heat = [], []
        sigma = 1 / k_av
        
        t = -1
        for temp in temps2sp:
            t += 1
        
            temp_eff1 = np.full(N_pop, temp)
            temp_eff2 = np.full(N_pop, temp)
        
            for pditer in range(T_eq):
        
                # Local temperature distributions
                for i in range(N_pop):
        
                    # Network 1 *******************
                    int_i2 = int_pairs[i][1]
                    if loc_mag2[int_i2] >= loc_thr:
                        temp_eff1[i] = temp / loc_mag2[int_i2]
                    else: 
                        temp_eff1[i] = temp / (loc_mag2[int_i2] + loc_thr)
        
                    # Network 2 *******************
                    ind2 = int_pairs[i][1]
                    int_i1 = int_pairs[i][0]
                    if loc_mag1[int_i1] >= loc_thr:
                        temp_eff2[ind2] = temp / loc_mag1[int_i1]
                    else: 
                        temp_eff2[ind2] = temp / (loc_mag1[int_i1] + loc_thr)   
        
                # Calculating instantaneous cavities and local magnetizations
                for i in range(N_pop):
        
                    # Network 1 *******************
                    Sum1 = err_heat
                    for index in nn_list1[i]: 
                        Sum1 += kernel2sp(cav_array1[index][0], temp_eff1[i])
                    cav_array1[i] = Sum1 
                    mag_list1[i][t] = np.tanh(cav_array1[i] / temp_eff1[i])
        
                    # Network 2 *******************
                    Sum2 = err_heat
                    for index in nn_list2[i]: 
                        Sum2 += kernel2sp(cav_array2[index][0], temp_eff2[i])
                    cav_array2[i] = Sum2 
                    mag_list2[i][t] = np.tanh(cav_array2[i] / temp_eff2[i])
        
                # Calculating the local OPs at this PD's iteration
                for i in range(N_pop):
        
                    # Network 1 *******************
                    loc_sum1 = 0
                    for index in nn_list1[i]:
                        loc_sum1 += mag_list1[index][t] / deg_list1[i]
                    loc_mag1[i] = loc_sum1 + np.random.normal(0, sigma, 1)
        
                    # Network 2 *******************
                    loc_sum2 = 0
                    for index in nn_list2[i]:
                        loc_sum2 += mag_list2[index][t] / deg_list2[i]
                    loc_mag2[i] = loc_sum2 + np.random.normal(0, sigma, 1)
        
            mag_array1, mag_array2 = [], []
            for i in range(N_pop):
                mag_array1.append(mag_list1[i][t])
                mag_array2.append(mag_list2[i][t])
        
            mag2spin1_heat.append(np.mean(mag_array1))
            mag2spin2_heat.append(np.mean(mag_array2))
            print('k_av is:', k_av, 'The run is:', run, 'The temperature is:', temp, 'The mags are:', np.mean(mag_array1), np.mean(mag_array2))
        
            if np.isnan(mag2spin1_heat[-1]):
                cav_array1 = k_av * np.random.random_sample((N_pop, 1))
                cav_array2 = k_av * np.random.random_sample((N_pop, 1))
                mag_list1, mag_list2 = np.zeros((N_pop, len(temps2sp))), np.zeros((N_pop, len(temps2sp)))
                loc_mag1, loc_mag2 = np.full(N_pop, 1.0), np.full(N_pop, 1.0)
                
            if mag2spin1_heat[-1] < 1/10 or mag2spin2_heat[-1] < 1/10:
                break
        
        len_m = int(len(mag2spin1_heat) - 1)
        Tsp_mix_loc_fluc = temps2sp[len_m] / k_av
        Tsp_kav_fluc.append(Tsp_mix_loc_fluc)
    
    # Collecting the data and saving them in an external file
    Tsp_kav_q1_fluct.append( (k_av, np.mean(Tsp_kav_fluc), np.std(Tsp_kav_fluc)) )
    np.savetxt('/Users/physivan/Desktop/Tsp_kav_q1_M2_fluc1.txt', Tsp_kav_q1_fluct)

print("--- %s seconds ---" % (time.time() - start_time))




''' ****************************************************************** '''
''' ****************************************************************** '''
'''              INTERPOP, 3 layers, q = 1 : varying <k>               '''
''' ****************************************************************** '''
''' ****************************************************************** '''

import matplotlib.pyplot as plt
plt.style.use('classic')

import numpy as np
import random as rm
import time


''' *********************************************************************************** '''
''' ------------------- Segregated time scales (thermal) algorithm  ------------------- '''
''' *********************************************************************************** '''

# Kernel functions solving the cavity fields via saddle-point
def kernel2sp(x, y):
    return (y / 2) * np.log((np.cosh((x + 1.0) / y)) / (np.cosh((x - 1.0) / y)))


# ***************** Multilayer population settings *******************
k_av = 4 
N_pop = 10**3 # population size

# Indices identifying independent and interdependent populations
net_indices1 = list(range(0, N_pop))
net_indices2 = rm.sample(net_indices1, N_pop) # randomly interdependent
net_indices3 = rm.sample(net_indices2, N_pop) # randomly interdependent
int_pairs = [x for x in zip(net_indices1, net_indices2, net_indices3)] #interdependent pairs

# Initializing lists of degrees and nearest neighbors
deg_list1, deg_list2, deg_list3 = [], [], []
nn_list1, nn_list2, nn_list3 = [], [], []

for i in range(N_pop):
    
    # Network 1 ******************************
    deg_list1.append(np.random.poisson(k_av, 1)[0])
    nn_list1.append(rm.sample(np.delete(net_indices1, i).tolist(), deg_list1[-1])) 
    
    # Network 2 ******************************
    deg_list2.append(np.random.poisson(k_av, 1)[0])
    nn_list2.append(rm.sample(np.delete(net_indices2, i).tolist(), deg_list2[-1]))  
    
    # Network 3 ******************************
    deg_list3.append(np.random.poisson(k_av, 1)[0])
    nn_list3.append(rm.sample(np.delete(net_indices3, i).tolist(), deg_list3[-1]))  


# Initial conditions: random, then by annealed heating/cooling
T = 10**1 # popdyn convergence
NOIs = 10**1 # interdependent iterations
cav_array1 = k_av * np.random.random_sample((N_pop, 1))
cav_array2 = k_av * np.random.random_sample((N_pop, 1))
cav_array3 = k_av * np.random.random_sample((N_pop, 1))
loc_mag1 = np.full(N_pop, 1.0)
loc_mag2 = np.full(N_pop, 1.0)
loc_mag3 = np.full(N_pop, 1.0)


# ******************* HEATING PART **************************
start_time = time.time()

magNOI1_th_loc_heat, magNOI2_th_loc_heat, magNOI3_th_loc_heat = [], [], []
temps2sp_M3 = np.linspace(0.01 * k_av, 0.4 * k_av, 10**2)
err_heat = 1e-15

for temp in temps2sp_M3:
    
    print('The heating temperatures is:', temp)
    # Generating two arrays of local temperatures
    temp_eff1 = np.full(N_pop, temp)
    temp_eff2 = np.full(N_pop, temp)
    temp_eff3 = np.full(N_pop, temp)
    
    for run in range(NOIs):
        
        
        # ****************** NETWORK 1 *********************
        
        # Local temperatures: network 1
        for i in range(N_pop):
            
            int_i2 = int_pairs[i][1]
            int_i3 = int_pairs[i][2]
            if loc_mag2[int_i2] >= 1e-10 and loc_mag3[int_i3] >= 1e-15:
                temp_eff1[i] = temp / (loc_mag2[int_i2] * loc_mag3[int_i3])
            elif loc_mag2[int_i2] <= 1e-10 or loc_mag3[int_i3] <= 1e-15:
                temp_eff1[i] = temp / (loc_mag2[int_i2] * loc_mag3[int_i3] + 1e-15)
            
        # Population dynamics: network 1 
        for pditer in range(T):
            for i in range(N_pop):
                Sum1 = err_heat
                for index in nn_list1[i]:
                    Sum1 += kernel2sp(cav_array1[index][0], temp_eff1[i])
                cav_array1[i] = Sum1    
                
        # Calculating the average spins' magnetization          
        mag_array1 = []
        for i in range(N_pop):
            a = cav_array1[i]
            mag_array1.append(np.tanh(a / temp_eff1[i]))
           
        # Calculating the local magnetization of spins' neighborhoods
        for i in range(N_pop):
            loc_sum1 = 0
            for index in nn_list1[i]:
                loc_sum1 += mag_array1[index] / deg_list1[i]
            loc_mag1[i] = loc_sum1


        # ****************** NETWORK 2 *********************
        
        # Local temperatures: network 2
        for i in range(N_pop):
            
            ind2 = int_pairs[i][1]
            int_i1 = int_pairs[i][0]
            int_i3 = int_pairs[i][2]
            if loc_mag1[int_i1] >= 1e-10 and loc_mag3[int_i3] >= 1e-10:
                temp_eff2[ind2] = temp / (loc_mag1[int_i1] * loc_mag3[int_i3])
            elif loc_mag1[int_i1] <= 1e-10 or loc_mag3[int_i3] <= 1e-10:
                temp_eff2[ind2] = temp / (loc_mag1[int_i1] * loc_mag3[int_i3] + 1e-10)
        
        # Population dynamics: network 2
        for pditer in range(T):                
            for i in range(N_pop):
                Sum2 = err_heat
                for index in nn_list2[i]:
                    Sum2 += kernel2sp(cav_array2[index][0], temp_eff2[i])
                cav_array2[i] = Sum2
                                
        # Calculating the average spins' magnetization          
        mag_array2 = []
        for i in range(N_pop):
            a = cav_array2[i]
            mag_array2.append(np.tanh(a / temp_eff2[i]))
           
        # Calculating the local magnetization of spins' neighborhoods
        for i in range(N_pop):
            loc_sum2 = 0
            for index in nn_list2[i]:
                loc_sum2 += mag_array2[index] / deg_list2[i]
            loc_mag2[i] = loc_sum2
            
        
        # ****************** NETWORK 3 *********************
        
        # Local temperatures: network 3
        for i in range(N_pop):
            
            ind3 = int_pairs[i][2]
            int_i1 = int_pairs[i][0]
            int_i2 = int_pairs[i][1]
            if loc_mag1[int_i1] >= 1e-10 and loc_mag2[int_i2] >= 1e-10:
                temp_eff3[ind3] = temp / (loc_mag1[int_i1] * loc_mag2[int_i2])
            elif loc_mag1[int_i1] < 1e-10 or loc_mag2[int_i2] < 1e-10:
                temp_eff3[ind3] = temp / (loc_mag1[int_i1] * loc_mag2[int_i2] + 1e-10)
        
        # Population dynamics: network 3
        for pditer in range(T):                
            for i in range(N_pop):
                Sum3 = err_heat
                for index in nn_list3[i]:
                    Sum3 += kernel2sp(cav_array3[index][0], temp_eff3[i])
                cav_array3[i] = Sum3
                                
        # Calculating the average spins' magnetization          
        mag_array3 = []
        for i in range(N_pop):
            a = cav_array3[i]
            mag_array3.append(np.tanh(a / temp_eff3[i]))
           
        # Calculating the local magnetization of spins' neighborhoods
        for i in range(N_pop):
            loc_sum3 = 0
            for index in nn_list3[i]:
                loc_sum3 += mag_array3[index] / deg_list3[i]
            loc_mag3[i] = loc_sum3

    magNOI1_th_loc_heat.append(np.mean(mag_array1))
    magNOI2_th_loc_heat.append(np.mean(mag_array2))
    magNOI3_th_loc_heat.append(np.mean(mag_array3))
    
    print('The layers magnetizations are: ', magNOI1_th_loc_heat[-1], \
          magNOI2_th_loc_heat[-1], magNOI3_th_loc_heat[-1])


# ******************* COOLING PART **************************
magNOI1_th_loc_cool, magNOI2_th_loc_cool, magNOI3_th_loc_cool = [], [], []
loc_mag1, loc_mag2, loc_mag3 = np.full(N_pop, 1.0), np.full(N_pop, 1.0), np.full(N_pop, 1.0)
mag2spin1_cool, mag2spin2_cool, mag2spin3_cool= [], [], []
err_cool = 1e-15

for temp in temps2sp_M3[::-1]:
    
    print('The cooling temperatures is:', temp)
    # Generating two arrays of local temperatures
    temp_eff1 = np.full(N_pop, temp)
    temp_eff2 = np.full(N_pop, temp)
    temp_eff3 = np.full(N_pop, temp)
    
    # Cavities to the negative direction for catching transcritical bif'n
    cav_array1 = - k_av * np.random.random_sample((N_pop, 1))
    cav_array2 = - k_av * np.random.random_sample((N_pop, 1))
    cav_array3 = - k_av * np.random.random_sample((N_pop, 1))
    
    for run in range(NOIs):
        
        # ****************** NETWORK 1 *********************
        
        # Local temperatures: network 1
        for i in range(N_pop):
            
            int_i2 = int_pairs[i][1]
            int_i3 = int_pairs[i][2]
            if loc_mag2[int_i2] >= 1e-10 and loc_mag3[int_i3] >= 1e-15:
                temp_eff1[i] = temp / (loc_mag2[int_i2] * loc_mag3[int_i3])
            elif loc_mag2[int_i2] <= 1e-10 or loc_mag3[int_i3] <= 1e-15:
                temp_eff1[i] = temp / (loc_mag2[int_i2] * loc_mag3[int_i3] + 1e-15)
            
        # Population dynamics: network 1 
        for pditer in range(T):
            for i in range(N_pop):
                Sum1 = err_heat
                for index in nn_list1[i]:
                    Sum1 += kernel2sp(cav_array1[index][0], temp_eff1[i])
                cav_array1[i] = Sum1    
                
        # Calculating the average spins' magnetization          
        mag_array1 = []
        for i in range(N_pop):
            a = cav_array1[i]
            mag_array1.append(np.tanh(a / temp_eff1[i]))
           
        # Calculating the local magnetization of spins' neighborhoods
        for i in range(N_pop):
            loc_sum1 = 0
            for index in nn_list1[i]:
                loc_sum1 += mag_array1[index] / deg_list1[i]
            loc_mag1[i] = loc_sum1


        # ****************** NETWORK 2 *********************
        
        # Local temperatures: network 2
        for i in range(N_pop):
            
            ind2 = int_pairs[i][1]
            int_i1 = int_pairs[i][0]
            int_i3 = int_pairs[i][2]
            if loc_mag1[int_i1] >= 1e-10 and loc_mag3[int_i3] >= 1e-10:
                temp_eff2[ind2] = temp / (loc_mag1[int_i1] * loc_mag3[int_i3])
            elif loc_mag1[int_i1] <= 1e-10 or loc_mag3[int_i3] <= 1e-10:
                temp_eff2[ind2] = temp / (loc_mag1[int_i1] * loc_mag3[int_i3] + 1e-10)
        
        # Population dynamics: network 2
        for pditer in range(T):                
            for i in range(N_pop):
                Sum2 = err_heat
                for index in nn_list2[i]:
                    Sum2 += kernel2sp(cav_array2[index][0], temp_eff2[i])
                cav_array2[i] = Sum2
                                
        # Calculating the average spins' magnetization          
        mag_array2 = []
        for i in range(N_pop):
            a = cav_array2[i]
            mag_array2.append(np.tanh(a / temp_eff2[i]))
           
        # Calculating the local magnetization of spins' neighborhoods
        for i in range(N_pop):
            loc_sum2 = 0
            for index in nn_list2[i]:
                loc_sum2 += mag_array2[index] / deg_list2[i]
            loc_mag2[i] = loc_sum2
            
        
        # ****************** NETWORK 3 *********************
        
        # Local temperatures: network 3
        for i in range(N_pop):
            
            ind3 = int_pairs[i][2]
            int_i1 = int_pairs[i][0]
            int_i2 = int_pairs[i][1]
            if loc_mag1[int_i1] >= 1e-10 and loc_mag2[int_i2] >= 1e-10:
                temp_eff3[ind3] = temp / (loc_mag1[int_i1] * loc_mag2[int_i2])
            elif loc_mag1[int_i1] < 1e-10 or loc_mag2[int_i2] < 1e-10:
                temp_eff3[ind3] = temp / (loc_mag1[int_i1] * loc_mag2[int_i2] + 1e-10)
        
        # Population dynamics: network 3
        for pditer in range(T):                
            for i in range(N_pop):
                Sum3 = err_heat
                for index in nn_list3[i]:
                    Sum3 += kernel2sp(cav_array3[index][0], temp_eff3[i])
                cav_array3[i] = Sum3
                                
        # Calculating the average spins' magnetization          
        mag_array3 = []
        for i in range(N_pop):
            a = cav_array3[i]
            mag_array3.append(np.tanh(a / temp_eff3[i]))
           
        # Calculating the local magnetization of spins' neighborhoods
        for i in range(N_pop):
            loc_sum3 = 0
            for index in nn_list3[i]:
                loc_sum3 += mag_array3[index] / deg_list3[i]
            loc_mag3[i] = loc_sum3

    magNOI1_th_loc_cool.append(np.mean(mag_array1))
    magNOI2_th_loc_cool.append(np.mean(mag_array2))
    magNOI3_th_loc_cool.append(np.mean(mag_array3))
    
    print('The layers magnetizations are: ', magNOI1_th_loc_cool[-1], \
          magNOI2_th_loc_cool[-1], magNOI3_th_loc_cool[-1])


# Checking for the heating and cooling thresholds
for i in range(5, len(temps2sp_M3)): 
    if abs(magNOI1_th_loc_heat[i] - magNOI1_th_loc_heat[i - 1]) > 1e-1:
        break
Tsp_th_loc_M3_heat = round((temps2sp_M3[i] + temps2sp_M3[i - 1]) / (2 * k_av), 2)

for i in range(5, len(temps2sp_M3)): 
    if abs(magNOI1_th_loc_cool[i] - magNOI1_th_loc_cool[i - 1]) > 1e-1:
        break
Tsp_th_loc_M3_cool = round((temps2sp_M3[::-1][i] + temps2sp_M3[::-1][i - 1]) / (2 * k_av), 2)


# Plotting the results - heating curves
plt.plot(temps2sp_M3 / k_av, magNOI1_th_loc_heat, '.', color = 'red')
plt.plot(temps2sp_M3 / k_av, magNOI2_th_loc_heat, '.', color = 'orange') 
plt.plot(temps2sp_M3 / k_av, magNOI3_th_loc_heat, '.', color = 'blue') 
# Plotting the results - cooling curves
plt.plot(temps2sp_M3[::-1] / k_av, magNOI1_th_loc_cool, '.', color = 'red')
plt.plot(temps2sp_M3[::-1] / k_av, magNOI2_th_loc_cool, '.', color = 'orange') 
plt.plot(temps2sp_M3[::-1] / k_av, magNOI2_th_loc_cool, '.', color = 'blue') 
plt.title('$k_{av} =$ %s, $T_{sp}^{heat} = $ %s, $T_{sp}^{cool} = $ %s' %(k_av, Tsp_th_loc_M3_heat, Tsp_th_loc_M3_cool))    
plt.xlim(0.0, (k_av + 0.5) / k_av)
plt.ylim(0.0, 1.0)
plt.show()

print("--- %s seconds ---" % (time.time() - start_time))




''' ************************************************************************ '''
''' -----------  M3 - THERMAL (mixed) + synthetic fluctuations  ------------ '''
''' ************************************************************************ '''

import matplotlib.pyplot as plt
plt.style.use('classic')

import numpy as np
import random as rm
import time


# Kernel functions solving the cavity fields via saddle-point
def kernel2sp(x, y):
    return (y / 2) * np.log((np.cosh((x + 1.0) / y)) / (np.cosh((x - 1.0) / y)))


# ***************** MULTILAYER POPULATION SETTINGS *******************

N_pop = 10**3 # population size
T_eq = 10**2 # popdyn convergence
Tsp_kav_q1_fluct = []
loc_thr = 1e-10
err_heat = 1e-20
N_runs = 10**2

# Calling data from single runs to optimize interval of temperatures
data_th_fluc = np.loadtxt('/Users/physivan/Desktop/Tsp_kav_q1_M3_fluc_1run.txt')
k_avs = data_th_fluc[2:10, 0]
Tsp_kav_1run = data_th_fluc[2:10, 1]


start_time = time.time()

ik = -1
for k_av in k_avs: 

    ik += 1
    # Collecting statistics over runs for better numerical thresholds
    Tsp_kav_fluc = []

    for run in range(N_runs):
        # Indices identifying independent and interdependent populations
        net_indices1 = list(range(0, N_pop))
        net_indices2 = rm.sample(net_indices1, N_pop) # randomly interdependent - layer 2
        net_indices3 = rm.sample(net_indices2, N_pop) # randomly interdependent - layer 3
        int_pairs = [x for x in zip(net_indices1, net_indices2, net_indices3)] #interdependent pairs

        # Initializing lists of degrees and nearest neighbors
        deg_list1, deg_list2, deg_list3 = [], [], []
        nn_list1, nn_list2, nn_list3 = [], [], []

        for i in range(N_pop):

            # Network 1 ******************************
            deg_list1.append(np.random.poisson(k_av, 1)[0])
            nn_list1.append(rm.sample(np.delete(net_indices1, i).tolist(), deg_list1[-1])) 

            # Network 2 ******************************
            deg_list2.append(np.random.poisson(k_av, 1)[0])
            nn_list2.append(rm.sample(np.delete(net_indices2, i).tolist(), deg_list2[-1]))  

            # Network 3 ******************************
            deg_list3.append(np.random.poisson(k_av, 1)[0])
            nn_list3.append(rm.sample(np.delete(net_indices3, i).tolist(), deg_list3[-1]))


        # Initial conditions: first step is random, then by annealed heating/cooling
        temps2sp = np.linspace((Tsp_kav_1run[ik] - 0.1) * k_av, (Tsp_kav_1run[ik] + 0.1) * k_av, 20)
        cav_array1 = k_av * np.random.random_sample((N_pop, 1))
        cav_array2 = k_av * np.random.random_sample((N_pop, 1))
        cav_array3 = k_av * np.random.random_sample((N_pop, 1))
        mag_list1, mag_list2, mag_list3 = np.zeros((N_pop, len(temps2sp))), np.zeros((N_pop, len(temps2sp))), np.zeros((N_pop, len(temps2sp)))
        loc_mag1, loc_mag2, loc_mag3 = np.full(N_pop, 1.0), np.full(N_pop, 1.0), np.full(N_pop, 1.0)


        # ********************** Annealing part ************************

        mag2spin1_heat, mag2spin2_heat, mag2spin3_heat = [], [], []
        sigma = 1.15 / k_av

        t = -1
        for temp in temps2sp:
            t += 1

            temp_eff1 = np.full(N_pop, temp)
            temp_eff2 = np.full(N_pop, temp)
            temp_eff3 = np.full(N_pop, temp)

            for pditer in range(T_eq):

                # Local temperature distributions
                for i in range(N_pop):

                    # Network 1 *******************
                    int_i2 = int_pairs[i][1]
                    int_i3 = int_pairs[i][2]
                    if loc_mag2[int_i2] >= loc_thr and loc_mag3[int_i3] >= loc_thr:
                        temp_eff1[i] = temp / (loc_mag2[int_i2] * loc_mag3[int_i3])
                    elif loc_mag2[int_i2] < loc_thr or loc_mag3[int_i3] < loc_thr:
                        temp_eff1[i] = temp / (loc_mag2[int_i2] * loc_mag3[int_i3] + loc_thr)

                    # Network 2 *******************
                    ind2 = int_pairs[i][1]
                    int_i1 = int_pairs[i][0]
                    int_i3 = int_pairs[i][2]
                    if loc_mag1[int_i1] >= loc_thr and loc_mag3[int_i3] >= loc_thr:
                        temp_eff2[ind2] = temp / (loc_mag1[int_i1] * loc_mag3[int_i3])
                    elif loc_mag1[int_i1] < loc_thr or loc_mag3[int_i3] < loc_thr:
                        temp_eff2[ind2] = temp / (loc_mag1[int_i1] * loc_mag3[int_i3] + loc_thr)  

                    # Network 3 *******************
                    ind3 = int_pairs[i][2]
                    int_i1 = int_pairs[i][0]
                    int_i2 = int_pairs[i][1]
                    if loc_mag1[int_i1] >= loc_thr and loc_mag3[int_i2] >= loc_thr:
                        temp_eff3[ind3] = temp / (loc_mag1[int_i1] * loc_mag2[int_i2])
                    elif loc_mag1[int_i1] < loc_thr or loc_mag3[int_i2] < loc_thr:
                        temp_eff3[ind3] = temp / (loc_mag1[int_i1] * loc_mag2[int_i2] + loc_thr)

                # Calculating instantaneous cavities and local magnetizations
                for i in range(N_pop):

                    # Network 1 *******************
                    Sum1 = err_heat
                    for index in nn_list1[i]: 
                        Sum1 += kernel2sp(cav_array1[index][0], temp_eff1[i])
                    cav_array1[i] = Sum1 
                    mag_list1[i][t] = np.tanh(cav_array1[i] / temp_eff1[i])

                    # Network 2 *******************
                    Sum2 = err_heat
                    for index in nn_list2[i]: 
                        Sum2 += kernel2sp(cav_array2[index][0], temp_eff2[i])
                    cav_array2[i] = Sum2 
                    mag_list2[i][t] = np.tanh(cav_array2[i] / temp_eff2[i])

                    # Network 3 *******************
                    Sum3 = err_heat
                    for index in nn_list3[i]: 
                        Sum3 += kernel2sp(cav_array3[index][0], temp_eff3[i])
                    cav_array3[i] = Sum3
                    mag_list3[i][t] = np.tanh(cav_array3[i] / temp_eff3[i])

                # Calculating the local OPs at this PD's iteration
                for i in range(N_pop):

                    # Network 1 *******************
                    loc_sum1 = 0
                    for index in nn_list1[i]:
                        loc_sum1 += mag_list1[index][t] / deg_list1[i]
                    loc_mag1[i] = loc_sum1 + np.random.normal(0, sigma, 1)

                    # Network 2 *******************
                    loc_sum2 = 0
                    for index in nn_list2[i]:
                        loc_sum2 += mag_list2[index][t] / deg_list2[i]
                    loc_mag2[i] = loc_sum2 + np.random.normal(0, sigma, 1)

                    # Network 2 *******************
                    loc_sum3 = 0
                    for index in nn_list3[i]:
                        loc_sum3 += mag_list3[index][t] / deg_list3[i]
                    loc_mag3[i] = loc_sum3 + np.random.normal(0, sigma, 1)


            mag_array1, mag_array2, mag_array3 = [], [], []
            for i in range(N_pop):
                mag_array1.append(mag_list1[i][t])
                mag_array2.append(mag_list2[i][t])
                mag_array3.append(mag_list3[i][t])

            mag2spin1_heat.append(np.mean(mag_array1))
            mag2spin2_heat.append(np.mean(mag_array2))
            mag2spin3_heat.append(np.mean(mag_array3))
            print('k_av is:', k_av, 'The run is:', run, 'The temperature is:', temp, 'The mags are:', np.mean(mag_array1), np.mean(mag_array2), np.mean(mag_array3))

            if np.isnan(mag2spin1_heat[-1]):
                cav_array1 = k_av * np.random.random_sample((N_pop, 1))
                cav_array2 = k_av * np.random.random_sample((N_pop, 1))
                cav_array3 = k_av * np.random.random_sample((N_pop, 1))
                mag_list1, mag_list2, mag_list3 = np.zeros((N_pop, len(temps2sp))), np.zeros((N_pop, len(temps2sp))), np.zeros((N_pop, len(temps2sp)))
                loc_mag1, loc_mag2, loc_mag3 = np.full(N_pop, 1.0), np.full(N_pop, 1.0), np.full(N_pop, 1.0)

            if mag2spin1_heat[-1] < 1/10 or mag2spin2_heat[-1] < 1/10 or mag2spin3_heat[-1] < 1/10:
                break

        len_m = int(len(mag2spin1_heat) - 1)
        Tsp_mix_loc_fluc = temps2sp[len_m] / k_av
        Tsp_kav_fluc.append(Tsp_mix_loc_fluc)

    # Collecting the data and saving them in an external file
    Tsp_kav_q1_fluct.append( (k_av, np.mean(Tsp_kav_fluc), np.std(Tsp_kav_fluc)) )
    np.savetxt('/Users/physivan/Desktop/Tsp_kav_q1_M3_fluc2.txt', Tsp_kav_q1_fluct)

print("--- %s seconds ---" % (time.time() - start_time))



''' ************************************************************************ '''
''' --------------- Hyphergraph : random 4-spin interactions --------------- '''
''' ************************************************************************ '''

import matplotlib.pyplot as plt
plt.style.use('classic')

import numpy as np
import random as rm
import time

# Integral kernels solving the cavity fields via saddle-point
def kernel2sp(x, y):
    return (y / 2) * np.log((np.cosh((x + 1.0) / y)) / (np.cosh((x - 1.0) / y)))

def kernel4sp(x1, x2, x3, y):
    alpha4sp_p = np.cosh((+ 1 + x3) / y) * np.cosh((x2 + x1) / y) + \
        np.cosh((+ 1 - x3) / y) * np.cosh((x2 - x1) / y)
    alpha4sp_m = np.cosh((- 1 + x3) / y) * np.cosh((x2 + x1) / y) + \
        np.cosh((- 1 - x3) / y) * np.cosh((x2 - x1) / y)
    return (y / 2) * np.log(alpha4sp_p / alpha4sp_m)


# Function grouping list of elements in n-uples
def grouped_list(lst, n):
    if len(lst) == 0:
        return []
    else:
        return list(zip(*[iter(lst)] * n))
    
    
# **************** Structural part ******************
k_av = 4
N_pop = 10**3
net_indices = (np.arange(N_pop))
deg_list, nn_quad_list = [], []
T = 10**2
cav_array = k_av * np.random.random_sample((N_pop, 1))
temps4sp = np.linspace(0.01 * k_av, 0.4 * k_av, 10**2)

for i in range(N_pop): 
    deg_list.append(np.random.poisson(k_av, 1)[0])
    indices_i = rm.sample(np.delete(net_indices, i).tolist(), 3 * deg_list[-1])
    nn_quad_list.append(grouped_list(indices_i, 3))


start_time = time.time()
# ****************** HEATING PART **********************
mag4sp_heat = []

for temp in temps4sp:
    
    for runs in range(T):
        for i in range(N_pop):
                        
            Sum = 0
            for index in nn_quad_list[i]:
                Sum += kernel4sp(cav_array[index[2]], cav_array[index[1]], cav_array[index[0]], temp)
            cav_array[i] = Sum 
                    
    mag_array4sp = []
    for i in range(N_pop):
        a = cav_array[i]
        mag_array4sp.append(np.tanh(a / temp))

    mag4sp_heat.append(np.mean(mag_array4sp))
    

# ****************** COOLING PART **********************
mag4sp_cool = []

for temp in temps4sp[::-1]:
    
    for runs in range(T):
        for i in range(N_pop):
                        
            Sum = 0
            for index in nn_quad_list[i]:
                Sum += kernel4sp(cav_array[index[2]], cav_array[index[1]], cav_array[index[0]], temp)
            cav_array[i] = Sum 
                    
    mag_array4sp = []
    for i in range(N_pop):
        a = cav_array[i]
        mag_array4sp.append(np.tanh(a / temp))

    mag4sp_cool.append(np.mean(mag_array4sp))


# Searching for the heating threshold
for i in range(5, len(temps4sp)): 
    if abs(mag4sp_heat[i] - mag4sp_heat[i - 1]) > 1e-1:
        break
Tsp_hy_4sp_heat = round((temps4sp[i] + temps4sp[i - 1]) / (2 * k_av), 2)

# Searching for the cooling threshold
for i in range(5, len(temps4sp)): 
    if abs(mag4sp_cool[i] - mag4sp_cool[i - 1]) > 1e-1:
        break
Tsp_hy_4sp_cool = round((temps4sp[::-1][i] + temps4sp[::-1][i - 1]) / (2 * k_av), 2)
    

# Plotting the results
plt.plot(temps4sp / k_av, mag4sp_heat, '-', color = 'k', label = '3-spin hypergraph')
plt.plot(temps4sp[::-1] / k_av, mag4sp_cool, '-', color = 'k')
plt.legend(loc = 'best')
plt.title(' $T_{sp}^{heat} / k_{av} =$ %s, $T_{sp}^{cool} / k_{av} =$ %s' %(Tsp_hy_4sp_heat, Tsp_hy_4sp_cool))
plt.ylim(0.0, 1.0)
plt.xlim(0.0, (k_av + 0.5) / k_av)
plt.show()

print("--- %s seconds ---" % (time.time() - start_time))




'''    ****************************************************************************    '''
'''    ****************************************************************************    '''
'''    ---------  Random (p+q)-spin interactions vs Interdependent spins  ---------    '''
'''    ****************************************************************************    '''
'''    ****************************************************************************    '''

# Code solving the case of 2+q and 3+q spin models. 
# Random ER hypergraphs vs randomly interdependent networks.
# It contains 1) population dynamics for 2+q and 3+q spin models;
# 2) Interdependent population dynamics: global OP and local OP interactions
# 3) Interdependent population dynamics with synthetic fluctuations.



'''    ****************************************************************************    '''
'''    -------------  2q3 - THERMAL (local) + synthetic fluctuations --------------    '''
'''    ****************************************************************************    '''

# We adopt the mixed time scale settings instead of proper thermal.
# This is numerically equivalent to the thermal one, but O(N) faster. 

import matplotlib.pyplot as plt
plt.style.use('classic')

import numpy as np
import random as rm
import time


# Kernel functions solving the cavity fields via saddle-point
def kernel2sp(x, y):
    return (y / 2) * np.log((np.cosh((x + 1.0) / y)) / (np.cosh((x - 1.0) / y)))


# ---------------- Static parameters of the runs ----------------
k_av = 4
N_pop = 10**3
T = 10**2
err_heat, err_cool = 1e-20, 1e-20
loc_thr = 1e-10
Tsp_2q3_th_k4_fluc = []
runs = 10**0
qvals = np.arange(0.00, 1.05, 0.05)


''' ******************* Starting the script ******************* '''

start_time = time.time()

for q in qvals:
    
    # Threshold vectors and q-dependent parameters
    Tsp_2q3_heat, Tsp_2q3_cool = [], []
    N_pop_q = int(q * N_pop)
    sigma_q = 1.1 * q / k_av # Comment this line to remove fluctuations
    
    
    # ********* Collecting runs for statistics *********
    for run in range(runs):
        
        # Structure of the populations at this run        
        net_indices1 = list(range(0, N_pop))
        net_indices2 = rm.sample(net_indices1, N_pop)
        int_pairs = [x for x in zip(net_indices1, net_indices2)] 
        
        deg_list1, deg_list2, nn_list1, nn_list2 = [], [], [], []
        for i in range(N_pop):
            
            deg_list1.append(np.random.poisson(k_av, 1)[0])
            nn_list1.append(rm.sample(np.delete(net_indices1, i).tolist(), deg_list1[i])) 
            
            deg_list2.append(np.random.poisson(k_av, 1)[0])
            nn_list2.append(rm.sample(np.delete(net_indices2, i).tolist(), deg_list2[i])) 
        
        
        # **************** Heating part ******************
        temps2sp = np.arange(0.01 * k_av, 1.1 * k_av, 0.01 * k_av)
        
        cav_array1 = k_av * np.random.random_sample((N_pop, 1)) 
        cav_array2 = k_av * np.random.random_sample((N_pop, 1)) 
        
        mag_list1, mag_list2 = np.zeros((N_pop, len(temps2sp))), np.zeros((N_pop, len(temps2sp)))
        loc_mag1, loc_mag2 = np.full(N_pop, 1.0), np.full(N_pop, 1.0)
        mag2spin1_heat, mag2spin2_heat = [], []
        
        t = -1
        for temp in temps2sp:
            t += 1
            
            temp_eff1 = np.full(N_pop, temp)
            temp_eff2 = np.full(N_pop, temp)
            print('Thermal case. q =', q, 'run =', run, '- The reduced heating temperature is:', temp / k_av)
        
            for pditer in range(T):
        
                for i in range(N_pop_q):
        
                    int_i2 = int_pairs[i][1]
                    if loc_mag2[int_i2] >= loc_thr:
                        temp_eff1[i] = temp / loc_mag2[int_i2]
                    else: 
                        temp_eff1[i] = temp / (loc_mag2[int_i2] + loc_thr)
        
                    ind2 = int_pairs[i][1]
                    int_i1 = int_pairs[i][0]
                    if loc_mag1[int_i1] >= loc_thr:
                        temp_eff2[ind2] = temp / loc_mag1[int_i1]
                    else: 
                        temp_eff2[ind2] = temp / (loc_mag1[int_i1] + loc_thr)
                
                # Calculating instantaneous cavities and local magnetizations
                for i in range(N_pop):
                    
                    Sum1 = err_heat
                    for index in nn_list1[i]: 
                        Sum1 += kernel2sp(cav_array1[index][0], temp_eff1[i])
                    cav_array1[i] = Sum1
                    mag_list1[i][t] = np.tanh(cav_array1[i] / temp_eff1[i])
                    
                    Sum2 = err_heat
                    for index in nn_list2[i]: 
                        Sum2 += kernel2sp(cav_array2[index][0], temp_eff2[i])
                    cav_array2[i] = Sum2
                    mag_list2[i][t] = np.tanh(cav_array2[i] / temp_eff2[i])
                
                # Calculating the local OPs at this PD's iteration
                for i in range(N_pop):
                    
                    loc_sum1 = 0
                    for index in nn_list1[i]:
                        loc_sum1 += mag_list1[index][t] / deg_list1[i]
                    # Adding synthetic fluctuations - if not needed remove random
                    loc_mag1[i] = loc_sum1 + np.random.normal(0, sigma_q, 1)
        
                    loc_sum2 = 0
                    for index in nn_list2[i]:
                        loc_sum2 += mag_list2[index][t] / deg_list2[i]
                    # Adding synthetic fluctuations - if not needed remove random
                    loc_mag2[i] = loc_sum2 + np.random.normal(0, sigma_q, 1)
                    
            mag_array_heat1, mag_array_heat2 = [], []
            for i in range(N_pop):
                mag_array_heat1.append(mag_list1[i][t])
                mag_array_heat2.append(mag_list2[i][t])
        
            mag2spin1_heat.append(np.mean(mag_array_heat1))
            mag2spin2_heat.append(np.mean(mag_array_heat2))
            print('The magnetizations are:', mag2spin1_heat[-1], mag2spin2_heat[-1])
            
            if mag2spin1_heat[-1] < 1e-2 or mag2spin2_heat[-1] < 1e-2:
                break
        
        # Collecting the heating thresholds
        Tsp_2q3_heat.append(temp / k_av)
            
        
        # **************** Cooling part ******************
        mag2spin1_cool, mag2spin2_cool = [], []
        
        t = -1
        for temp in temps2sp[::-1]:
            t += 1
            
            temp_eff1 = np.full(N_pop, temp)
            temp_eff2 = np.full(N_pop, temp)
            print('Thermal alg: this is q =', q, 'run_q =', run, '- The reduced cooling temperature is:', temp / k_av)
            
            # Initializing the cavities down for catching the transcritical instabilities
            mag_list1, mag_list2 = np.zeros((N_pop, len(temps2sp))), np.zeros((N_pop, len(temps2sp)))
            cav_array1 = - k_av * np.random.random_sample((N_pop, 1)) 
            cav_array2 = - k_av * np.random.random_sample((N_pop, 1)) 
            
        
            for pditer in range(T):
        
                for i in range(N_pop_q):
        
                    int_i2 = int_pairs[i][1]
                    if loc_mag2[int_i2] >= loc_thr:
                        temp_eff1[i] = temp / loc_mag2[int_i2]
                    else: 
                        temp_eff1[i] = temp / (loc_mag2[int_i2] + loc_thr)
        
                    ind2 = int_pairs[i][1]
                    int_i1 = int_pairs[i][0]
                    if loc_mag1[int_i1] >= loc_thr:
                        temp_eff2[ind2] = temp / loc_mag1[int_i1]
                    else: 
                        temp_eff2[ind2] = temp / (loc_mag1[int_i1] + loc_thr)      
                
                # Calculating instantaneous cavities and local magnetizations
                for i in range(N_pop):
                    
                    Sum1 = err_cool
                    for index in nn_list1[i]: 
                        Sum1 += kernel2sp(cav_array1[index][0], temp_eff1[i])
                    cav_array1[i] = Sum1
                    mag_list1[i][t] = np.tanh(cav_array1[i] / temp_eff1[i])
                    
                    Sum2 = err_cool
                    for index in nn_list2[i]: 
                        Sum2 += kernel2sp(cav_array2[index][0], temp_eff2[i])
                    cav_array2[i] = Sum2
                    mag_list2[i][t] = np.tanh(cav_array2[i] / temp_eff2[i])
                
                
                # Calculating the local OPs at this PD's iteration
                for i in range(N_pop):
                    
                    loc_sum1 = 0
                    for index in nn_list1[i]:
                        loc_sum1 += mag_list1[index][t] / deg_list1[i]
                    # Adding synthetic fluctuations - if not needed remove random
                    loc_mag1[i] = loc_sum1 + np.random.normal(0, sigma_q, 1)
        
                    loc_sum2 = 0
                    for index in nn_list2[i]:
                        loc_sum2 += mag_list2[index][t] / deg_list2[i]
                    # Adding synthetic fluctuations - if not needed remove random
                    loc_mag2[i] = loc_sum2 + np.random.normal(0, sigma_q, 1)
                    
            mag_array_cool1, mag_array_cool2 = [], []
            for i in range(N_pop):
                mag_array_cool1.append(mag_list1[i][t])
                mag_array_cool2.append(mag_list2[i][t])
        
            mag2spin1_cool.append(np.mean(mag_array_cool1))
            mag2spin2_cool.append(np.mean(mag_array_cool2))
            print('The layers magnetizations are:', mag2spin1_cool[-1], mag2spin2_cool[-1])
                
            if abs(mag2spin1_cool[-1]) >= 1/20 or abs(mag2spin2_cool[-1]) >= 1/20:
                break
            
        Tsp_2q3_cool.append(temp / k_av)
        
    # Saving the datapoints in an external file
    Tsp_2q3_th_k4_fluc.append( (q, np.mean(Tsp_2q3_cool), np.mean(Tsp_2q3_heat)) )
    np.savetxt('/Users/physivan/Desktop/Tsp_2q3_th_k4_fluct_1run.txt', Tsp_2q3_th_k4_fluc)
    
print("--- %s seconds ---" % (time.time() - start_time))



'''    *****************************************************************************    '''
'''    --------------------------  2q3 - THERMAL (global) --------------------------    '''
'''    *****************************************************************************    '''

import matplotlib.pyplot as plt
plt.style.use('classic')

import numpy as np
import random as rm
import time


# Integral kernels solving the saddle-point replica symmetric equations
def kernel2sp(x, y):
    return (y / 2) * np.log((np.cosh((x + 1.0) / y)) / (np.cosh((x - 1.0) / y)))


# ---------------- Static parameters of the runs ----------------
k_av = 4
N_pop = 10**3
T = 10**2
err_heat = 1e-20 
Tsp_th_k4 = []
runs = 10**2

# Optimizing the range of temperature after a single run search
data_Tsp_2q3_glob = np.loadtxt('/Users/physivan/Desktop/Tsp_2q3_th_k4_glob_100runs.txt')
q_vals = data_Tsp_2q3_glob[:, 0]
Tsp_glob_heat = data_Tsp_2q3_glob[:, 1]
Tsp_glob_cool = data_Tsp_2q3_glob[:, 2]


''' ------------- Starting the script --------------- '''
start_time = time.time()

j_q = -1
for q in q_vals:

    j_q += 1
    # Arrays for calculating mean thresholds
    Tsp_heat, Tsp_cool = [], []
    M_sp_heat_jump1, M_sp_heat_jump2 = [], []

    # Setting range of temperatures
    temps2q3sp_heat = np.linspace((Tsp_glob_heat[j_q] - 0.1) * k_av, (Tsp_glob_heat[j_q] + 0.1) * k_av, 20)
    temps2q3sp_cool = np.arange(0.01 * k_av, (Tsp_glob_cool[j_q] + 0.10) * k_av, 0.01 * k_av)
    
    # Collecting more runs for statistical significance
    for run in range(runs):

        ''' --------------------------------------------------------------- '''
        ''' Interdependent solution: mixed time scales, same saddle kernels '''
        ''' --------------------------------------------------------------- '''
        # Numerically equivalent to the thermal one, but O(N) faster. 

        N_pop_q = int(q * N_pop)
        net_indices1 = list(range(0, N_pop))
        net_indices2 = rm.sample(net_indices1, N_pop)
        int_pairs = [x for x in zip(net_indices1, net_indices2)] 

        deg_list1, deg_list2, nn_list1, nn_list2 = [], [], [], []
        for i in range(N_pop):

            deg_list1.append(np.random.poisson(k_av, 1)[0])
            nn_list1.append(rm.sample(np.delete(net_indices1, i).tolist(), deg_list1[i])) 

            deg_list2.append(np.random.poisson(k_av, 1)[0])
            nn_list2.append(rm.sample(np.delete(net_indices2, i).tolist(), deg_list2[i])) 


        # **************** Heating part ******************
       
        mag2spin1_heat, mag2spin2_heat = [], []
        mag_list1, mag_list2 = np.zeros((N_pop, len(temps2q3sp_heat))), np.zeros((N_pop, len(temps2q3sp_heat)))
        cav_array1 = k_av * np.random.random_sample((N_pop, 1)) 
        cav_array2 = k_av * np.random.random_sample((N_pop, 1))
        glob_mag1, glob_mag2 = 1.0, 1.0
       
        t = -1
        for temp in temps2q3sp_heat:
            t += 1

            temp_eff1 = np.full(N_pop, temp)
            temp_eff2 = np.full(N_pop, temp)
            print('Thermal case. q =', q, 'run =', run, '- The cooling temperature is:', temp)
           
            for pditer in range(T):

                for i in range(N_pop_q):

                    if glob_mag2 >= 1e-10:
                        temp_eff1[i] = temp / glob_mag2
                    else: 
                        temp_eff1[i] = temp / (glob_mag2 + 1e-10)

                    ind2 = int_pairs[i][1]
                    if glob_mag1 >= 1e-10:
                        temp_eff2[ind2] = temp / glob_mag1
                    else: 
                        temp_eff2[ind2] = temp / (glob_mag1 + 1e-10)   

                # Calculating instantaneous cavities and local magnetizations
                for i in range(N_pop):

                    Sum1 = err_heat
                    for index in nn_list1[i]: 
                        Sum1 += kernel2sp(cav_array1[index][0], temp_eff1[i])
                    cav_array1[i] = Sum1
                    mag_list1[i][t] = np.tanh(cav_array1[i] / temp_eff1[i])

                    Sum2 = err_heat
                    for index in nn_list2[i]: 
                        Sum2 += kernel2sp(cav_array2[index][0], temp_eff2[i])
                    cav_array2[i] = Sum2
                    mag_list2[i][t] = np.tanh(cav_array2[i] / temp_eff2[i])

                mag_array1, mag_array2 = [], []
                for i in range(N_pop):
                    mag_array1.append(mag_list1[i][t])
                    mag_array2.append(mag_list2[i][t])

                glob_mag1 = np.mean(mag_array1)
                glob_mag2 = np.mean(mag_array2)
            
            mag2spin1_heat.append(glob_mag1)
            mag2spin2_heat.append(glob_mag2)
            print('The layers magnetizations are:', glob_mag1, glob_mag2)

            if abs(glob_mag1) <= 1e-3 and abs(glob_mag2) <= 1e-3:
                break
        
        # Collecting the threshold and the magnetization jumps at this run
        Tsp_heat.append(temp / k_av)
        
        for n in range(len(mag2spin1_heat)): 
           if mag2spin1_heat[n] <= 1e-2:
              M_sp_heat_jump1.append(abs(mag2spin1_heat[n] - mag2spin1_heat[n-1]))
              break
        for n in range(len(mag2spin2_heat)): 
           if mag2spin2_heat[n] <= 1e-2:
              M_sp_heat_jump2.append(abs(mag2spin2_heat[n] - mag2spin2_heat[n-1]))
              break
        
        print('The jumps at this run are:', M_sp_heat_jump1[-1], M_sp_heat_jump2[-1])
        
        # **************** Cooling part ******************
        mag2spin1_cool, mag2spin2_cool = [], []

        t = -1
        for temp in temps2q3sp_cool[::-1]:
            t += 1

            temp_eff1 = np.full(N_pop, temp)
            temp_eff2 = np.full(N_pop, temp)
            print('Thermal case. q =', q, 'run =', run, '- The cooling temperature is:', temp)

            mag_list1, mag_list2 = np.zeros((N_pop, len(temps2q3sp_cool))), np.zeros((N_pop, len(temps2q3sp_cool)))
            glob_mag1, glob_mag2 = 1.0, 1.0
            cav_array1 = - k_av * np.random.random_sample((N_pop, 1)) 
            cav_array2 = - k_av * np.random.random_sample((N_pop, 1))
           
            for pditer in range(T):

                for i in range(N_pop_q):

                    if glob_mag2 >= 1e-10:
                        temp_eff1[i] = temp / glob_mag2
                    else: 
                        temp_eff1[i] = temp / (glob_mag2 + 1e-10)

                    ind2 = int_pairs[i][1]
                    if glob_mag1 >= 1e-10:
                        temp_eff2[ind2] = temp / glob_mag1
                    else: 
                        temp_eff2[ind2] = temp / (glob_mag1 + 1e-10)   

                # Calculating instantaneous cavities and local magnetizations
                for i in range(N_pop):

                    Sum1 = err_heat
                    for index in nn_list1[i]: 
                        Sum1 += kernel2sp(cav_array1[index][0], temp_eff1[i])
                    cav_array1[i] = Sum1
                    mag_list1[i][t] = np.tanh(cav_array1[i] / temp_eff1[i])

                    Sum2 = err_heat
                    for index in nn_list2[i]: 
                        Sum2 += kernel2sp(cav_array2[index][0], temp_eff2[i])
                    cav_array2[i] = Sum2
                    mag_list2[i][t] = np.tanh(cav_array2[i] / temp_eff2[i])

                mag_array1, mag_array2 = [], []
                for i in range(N_pop):
                    mag_array1.append(mag_list1[i][t])
                    mag_array2.append(mag_list2[i][t])

                glob_mag1 = np.mean(mag_array1)
                glob_mag2 = np.mean(mag_array2)
            
            mag2spin1_cool.append(glob_mag1)
            mag2spin2_cool.append(glob_mag2)
            print('The layers magnetizations are:', glob_mag1, glob_mag2)

            if abs(glob_mag1) >= 1/10 and abs(glob_mag2) >= 1/10:
                break
        
        # Collecting the threshold and the magnetization jumps at this run
        Tsp_cool.append(temp / k_av)
        
    ''' ------------ Collecting the thresholds ------------ '''
    Tsp_th_k4.append( (q, np.mean(Tsp_heat), np.std(Tsp_heat), \
                          np.mean(Tsp_cool), np.std(Tsp_cool), \
                          np.mean(M_sp_heat_jump1), np.std(M_sp_heat_jump1), \
                          np.mean(M_sp_heat_jump2), np.std(M_sp_heat_jump2)) )

    # Saving the datapoints in an external file
    np.savetxt('/Users/physivan/Desktop/Tsp_2q3_th_k4_glob_100runs.txt', Tsp_th_k4)


print("--- %s seconds ---" % (time.time() - start_time))



''' ---------------------------------------------------------- '''
'''   Hyphergraph - Population dynamics for (2+q)-spin model   '''
''' ---------------------------------------------------------- '''

import matplotlib.pyplot as plt
plt.style.use('classic')

import numpy as np
import random as rm
import time


# Integral kernels solving the saddle-point replica symmetric equations
def kernel2sp(x, y):
    return (y / 2) * np.log((np.cosh((x + 1.0) / y)) / (np.cosh((x - 1.0) / y)))

def kernel3sp(x1, x2, y):
    alpha_p = np.cosh((x2 + x1) / y) + np.exp(- 2 / y) * np.cosh((x2 - x1) / y)
    alpha_m = np.cosh((x2 - x1) / y) + np.exp(- 2 / y) * np.cosh((x2 + x1) / y)
    return (y / 2) * np.log(alpha_p / alpha_m)

# Function grouping list of elements in n-uples
def grouped_list(lst, n):
    if len(lst) == 0:
        return []
    else:
        return list(zip(*[iter(lst)] * n))


# **************** Structure settings ******************
k_av = 4
N_pop = 10**3
q = 0.0

# Generating the list of nearest neighbors according to their interaction pattern
net_indices = (np.arange(N_pop))
nn_2sp_list, nn_3sp_list = [], []

for i in range(N_pop): 
    
    # Each node has a given, Poisson distributed degree
    deg_i = np.random.poisson(k_av, 1)[0]
    
    # A fraction 1-q of its edges are involved in pairwise interactions
    deg_2sp_i = int(np.floor(deg_i * (1 - q)))
    indices_2sp_minus_i = np.delete(net_indices, i).tolist()
    indices_2sp_i = rm.sample(indices_2sp_minus_i, deg_2sp_i)
    nn_2sp_list.append(indices_2sp_i)
    
    # The remaining fraction is partecipating to 3-spin interactions
    deg_3sp_i = deg_i - deg_2sp_i
    indices_3sp_minus_2sp_minus_i = [a for a in indices_2sp_minus_i if a not in indices_2sp_i]
    indices_3sp_i = rm.sample(indices_3sp_minus_2sp_minus_i, 2 * deg_3sp_i)
    nn_3sp_list.append(grouped_list(indices_3sp_i, 2))


start_time = time.time()

# **************** Dynamics settings ******************
T = 10**2
cav_array = k_av * np.random.random_sample((N_pop, 1))
temps2q3sp = np.linspace(0.4 * k_av, 1.1 * k_av, 10**2)


# ****************** Heating part **********************
mag2q3sp_heat = []
for temp in temps2q3sp:
    
    print('The heating temperature is:', temp)
    for pditer in range(T):
        for i in range(N_pop):

            Sum = 0
            for index2sp in nn_2sp_list[i]: 
                Sum += kernel2sp(cav_array[index2sp], temp)
            for index3sp in nn_3sp_list[i]:
                Sum += kernel3sp(cav_array[index3sp[1]], cav_array[index3sp[0]], temp)
            cav_array[i] = Sum
                
    mag_array2q3sp = []
    for i in range(N_pop):
        a = cav_array[i]
        mag_array2q3sp.append(np.tanh(a / temp))

    mag2q3sp_heat.append(np.mean(mag_array2q3sp))
    print('The magnetization is:', mag_array2q3sp[-1])
    
    # Magnetization threshold: stopping annealed heating
    if mag2q3sp_heat[-1] < 1e-2:
        break
    
Tsp_hy_2q3_heat = temp / k_av
print('The heating spinodal is at:', Tsp_hy_2q3_heat)


# ****************** Cooling part **********************
mag2q3sp_cool = []
for temp in temps2q3sp[::-1]:
    
    cav_array = - k_av * np.random.random_sample((N_pop, 1))
    print('The cooling temperature is:', temp)
    for pditer in range(T):
        for i in range(N_pop):

            Sum = 0
            for index2sp in nn_2sp_list[i]: 
                Sum += kernel2sp(cav_array[index2sp], temp)
            for index3sp in nn_3sp_list[i]:
                Sum += kernel3sp(cav_array[index3sp[1]], cav_array[index3sp[0]], temp)
            cav_array[i] = Sum
                
    mag_array2q3sp = []
    for i in range(N_pop):
        a = cav_array[i]
        mag_array2q3sp.append(np.tanh(a / temp))

    mag2q3sp_cool.append(np.mean(mag_array2q3sp))
    print('The magnetization is:', mag_array2q3sp[-1])
    
    # Magnetization threshold: stopping annealed cooling
    if abs(mag2q3sp_heat[-1]) > 1e-1:
        break
    
Tsp_hy_2q3_cool = temp / k_av
print('The cooling transcritical point is at:', Tsp_hy_2q3_cool)

    
# Plotting the results
plt.plot(temps2q3sp / k_av, mag2q3sp_heat, '-', color = 'k', label = '3-spin hypergraph')
plt.plot(temps2q3sp[::-1] / k_av, mag2q3sp_cool, '-', color = 'k')
plt.legend(loc = 'best')
plt.title(' $T_{sp}^{heat} / k_{av} =$ %s, $T_{sp}^{cool} / k_{av} =$ %s' %(Tsp_hy_2q3_heat, Tsp_hy_2q3_cool))
plt.ylim(0.0, 1.0)
plt.xlim(0.0, (k_av + 0.5) / k_av)
plt.show()

print("--- %s seconds ---" % (time.time() - start_time))




'''    *****************************************************    '''
'''    -------------  3+q - THERMAL (global)  --------------    '''
'''    *****************************************************    '''
# Also here we adopt the mixed time scales. 
# Numerically equivalent to the thermal one with segregated time scales, but O(2N) faster. 

import matplotlib.pyplot as plt
plt.style.use('classic')

import numpy as np
import random as rm
import time


# Kernel functions solving the cavity fields via saddle-point
def kernel2sp(x, y):
    return (y / 2) * np.log((np.cosh((x + 1.0) / y)) / (np.cosh((x - 1.0) / y)))


 # ---------------- Static parameters of the runs ----------------
k_av = 4
N_pop = 10**3
T = 10**2
err_heat = 1e-20 
Tsp_3q4_glob = []
runs = 10**0
mag_thr = 1e-10
qvals = np.arange(0.00, 1.05, 0.05)


''' ------------- Starting the script --------------- '''
start_time = time.time()

for q in qvals: 

    Tsp_heat, Tsp_cool = [], []

    # Collecting more runs for statistical significance
    for run in range(runs):

        ''' --------------------------------------------------------------- '''
        ''' Interdependent solution: mixed time scales, same saddle kernels '''
        ''' --------------------------------------------------------------- '''
        # Numerically equivalent to the thermal one, but O(2N) faster. 


        # **************** Structural part ****************
        N_pop_q = int(q * N_pop)
        net_indices1 = list(range(0, N_pop))
        net_indices2 = rm.sample(net_indices1, N_pop)
        net_indices3 = rm.sample(net_indices2, N_pop)
        int_pairs = [x for x in zip(net_indices1, net_indices2, net_indices3)] 

        deg_list1, deg_list2, deg_list3 = [], [], []
        nn_list1, nn_list2, nn_list3 = [], [], []

        for i in range(N_pop):

            deg_list1.append(np.random.poisson(k_av, 1)[0])
            nn_list1.append(rm.sample(np.delete(net_indices1, i).tolist(), deg_list1[-1])) 

            deg_list2.append(np.random.poisson(k_av, 1)[0])
            nn_list2.append(rm.sample(np.delete(net_indices2, i).tolist(), deg_list2[-1])) 

            deg_list3.append(np.random.poisson(k_av, 1)[0])
            nn_list3.append(rm.sample(np.delete(net_indices3, i).tolist(), deg_list3[-1])) 


        ''' **************** Heating part ****************** '''

        temps_th_3q4sp_heat = np.linspace(0.31 * k_av, 0.61 * k_av, 30)
        mag2spin1_heat, mag2spin2_heat, mag2spin3_heat = [], [], []
        mag_list1 = np.zeros((N_pop, len(temps_th_3q4sp_heat)))
        mag_list2 = np.zeros((N_pop, len(temps_th_3q4sp_heat)))
        mag_list3 = np.zeros((N_pop, len(temps_th_3q4sp_heat)))
        cav_array1 = k_av * np.random.random_sample((N_pop, 1)) 
        cav_array2 = k_av * np.random.random_sample((N_pop, 1)) 
        cav_array3 = k_av * np.random.random_sample((N_pop, 1)) 
        glob_mag1, glob_mag2, glob_mag3 = 1.0, 1.0, 1.0

        t = -1
        for temp in temps_th_3q4sp_heat:
            t += 1

            temp_eff1 = np.full(N_pop, temp)
            temp_eff2 = np.full(N_pop, temp)
            temp_eff3 = np.full(N_pop, temp)
            print('Thermal case. q =', q, 'run =', run, '- The reduced heating temperature is:', temp / k_av)

            for pditer in range(T):

                ''' ------- Generating distribution of effective temperatures --------- '''       
                
                # Networks 1, 2 and 3 are directly dependent on each other closing the loop
                for i in net_indices1[:int((1 - q) * N_pop)]:

                    # Network 1 ----------------
                    if glob_mag2 >= mag_thr:
                        temp_eff1[i] = temp / glob_mag2
                    else: 
                        temp_eff1[i] = temp / (glob_mag2 + mag_thr)

                    # Network 2 ----------------
                    ind2 = int_pairs[i][1]
                    if glob_mag3 >= mag_thr:
                        temp_eff2[ind2] = temp / glob_mag3
                    else: 
                        temp_eff2[ind2] = temp / (glob_mag3 + mag_thr) 
                        
                   # Network 3 ----------------
                    ind3 = int_pairs[i][2]
                    if glob_mag1 >= mag_thr:
                        temp_eff3[ind3] = temp / glob_mag1
                    else: 
                        temp_eff3[ind3] = temp / (glob_mag1 + mag_thr)  

                # Fraction q of triads are closed in dependency bi-directionality
                for i in net_indices1[int((1 - q) * N_pop):]:

                    # Network 1 ----------------
                    if glob_mag2 >= mag_thr and glob_mag3 >= mag_thr:
                        temp_eff1[i] = temp / (glob_mag2 * glob_mag3)
                    elif glob_mag2 < mag_thr or glob_mag3 < mag_thr:
                        temp_eff1[i] = temp / (glob_mag2 * glob_mag3 + mag_thr)

                    # Network 2 ----------------
                    ind2 = int_pairs[i][1]
                    if glob_mag1 >= mag_thr and glob_mag3 >= mag_thr:
                        temp_eff2[ind2] = temp / (glob_mag1 * glob_mag3)
                    elif glob_mag1 < mag_thr or glob_mag3 < mag_thr:
                        temp_eff2[ind2] = temp / (glob_mag1 * glob_mag3 + mag_thr)    

                    # Network 3 ----------------
                    ind3 = int_pairs[i][2]
                    if glob_mag1 >= mag_thr and glob_mag2 >= mag_thr:
                        temp_eff3[ind3] = temp / (glob_mag1 * glob_mag2)
                    elif glob_mag1 < mag_thr or glob_mag2 < mag_thr:
                        temp_eff3[ind3] = temp / (glob_mag1 * glob_mag2 + mag_thr) 

                ''' --------- Calculating cavities and global magnetizations -------- '''                
                for i in range(N_pop):

                    # Network 1 ----------------
                    Sum1 = err_heat
                    for index in nn_list1[i]: 
                        Sum1 += kernel2sp(cav_array1[index][0], temp_eff1[i])
                    cav_array1[i] = Sum1
                    mag_list1[i][t] = np.tanh(cav_array1[i] / temp_eff1[i])

                    # Network 2 ----------------
                    Sum2 = err_heat
                    for index in nn_list2[i]: 
                        Sum2 += kernel2sp(cav_array2[index][0], temp_eff2[i])
                    cav_array2[i] = Sum2
                    mag_list2[i][t] = np.tanh(cav_array2[i] / temp_eff2[i])

                    # Network 3 ----------------
                    Sum3 = err_heat
                    for index in nn_list3[i]: 
                        Sum3 += kernel2sp(cav_array3[index][0], temp_eff3[i])
                    cav_array3[i] = Sum3
                    mag_list3[i][t] = np.tanh(cav_array3[i] / temp_eff3[i])

                mag_array1, mag_array2, mag_array3 = [], [], []
                for i in range(N_pop):
                    mag_array1.append(mag_list1[i][t])
                    mag_array2.append(mag_list2[i][t])
                    mag_array3.append(mag_list3[i][t])

                glob_mag1 = np.mean(mag_array1)
                glob_mag2 = np.mean(mag_array2)
                glob_mag3 = np.mean(mag_array3)

            # -------- Appending converged magnetizations at this temperature -------- 
            mag2spin1_heat.append(glob_mag1)
            mag2spin2_heat.append(glob_mag2)
            mag2spin3_heat.append(glob_mag3)
            print('The layers magnetizations are:', glob_mag1, glob_mag2, glob_mag3)

            if abs(glob_mag1) <= 1e-1 or abs(glob_mag1) <= 1e-1 or abs(glob_mag1) <= 1e-1:
                break

        # -------- Collecting the heating threshold at this run --------
        Tsp_heat.append(temp / k_av)


        ''' **************** Cooling part ****************** '''

        temps_th_3q4sp_cool = np.linspace(0.01 * k_av, 0.61 * k_av, 60)
        mag2spin1_cool, mag2spin2_cool, mag2spin3_cool = [], [], []
        mag_list1 = np.zeros((N_pop, len(temps_th_3q4sp_cool)))
        mag_list2 = np.zeros((N_pop, len(temps_th_3q4sp_cool)))
        mag_list3 = np.zeros((N_pop, len(temps_th_3q4sp_cool)))
        cav_array1 = - k_av * np.random.random_sample((N_pop, 1)) 
        cav_array2 = - k_av * np.random.random_sample((N_pop, 1))
        cav_array3 = - k_av * np.random.random_sample((N_pop, 1))
        #glob_mag1, glob_mag2, glob_mag3 = 1.0, 1.0, 1.0

        t = -1
        for temp in temps_th_3q4sp_cool[::-1]:
            t += 1

            temp_eff1 = np.full(N_pop, temp)
            temp_eff2 = np.full(N_pop, temp)
            temp_eff3 = np.full(N_pop, temp)
            print('Thermal case. q =', q, 'run =', run, '- The reduced cooling temperature is:', temp / k_av)

            for pditer in range(T):

                ''' --------- Generating distribution of effective temperatures --------- '''
                
                # Networks 1, 2 and 3 are directly dependent on each other closing the loop
                for i in net_indices1[:int((1 - q) * N_pop)]:

                    # Network 1 ----------------
                    if glob_mag2 >= mag_thr:
                        temp_eff1[i] = temp / glob_mag2
                    else: 
                        temp_eff1[i] = temp / (glob_mag2 + mag_thr)

                    # Network 2 ----------------
                    ind2 = int_pairs[i][1]
                    if glob_mag3 >= mag_thr:
                        temp_eff2[ind2] = temp / glob_mag3
                    else: 
                        temp_eff2[ind2] = temp / (glob_mag3 + mag_thr) 
                        
                   # Network 3 ----------------
                    ind3 = int_pairs[i][2]
                    if glob_mag1 >= mag_thr:
                        temp_eff3[ind3] = temp / glob_mag1
                    else: 
                        temp_eff3[ind3] = temp / (glob_mag1 + mag_thr)  

                # Fraction q of triads are closed in dependency bi-directionality    
                for i in net_indices1[int((1 - q) * N_pop):]:

                    # Network 1 ----------------
                    if glob_mag2 >= mag_thr and glob_mag3 >= mag_thr:
                        temp_eff1[i] = temp / (glob_mag2 * glob_mag3)
                    elif glob_mag2 < mag_thr or glob_mag3 < mag_thr:
                        temp_eff1[i] = temp / (glob_mag2 * glob_mag3 + mag_thr)

                    # Network 2 ----------------
                    ind2 = int_pairs[i][1]
                    if glob_mag1 >= mag_thr and glob_mag3 >= mag_thr:
                        temp_eff2[ind2] = temp / (glob_mag1 * glob_mag3)
                    elif glob_mag1 < mag_thr or glob_mag3 < mag_thr:
                        temp_eff2[ind2] = temp / (glob_mag1 * glob_mag3 + mag_thr)    

                    # Network 3 ----------------
                    ind3 = int_pairs[i][2]
                    if glob_mag1 >= mag_thr and glob_mag2 >= mag_thr:
                        temp_eff3[ind3] = temp / (glob_mag1 * glob_mag2)
                    elif glob_mag1 < mag_thr or glob_mag2 < mag_thr:
                        temp_eff3[ind3] = temp / (glob_mag1 * glob_mag2 + mag_thr) 

                ''' --------- Calculating cavities and global magnetizations -------- '''
                for i in range(N_pop):

                    # Network 1 ----------------
                    Sum1 = err_heat
                    for index in nn_list1[i]: 
                        Sum1 += kernel2sp(cav_array1[index][0], temp_eff1[i])
                    cav_array1[i] = Sum1
                    mag_list1[i][t] = np.tanh(cav_array1[i] / temp_eff1[i])

                    # Network 2 ----------------
                    Sum2 = err_heat
                    for index in nn_list2[i]: 
                        Sum2 += kernel2sp(cav_array2[index][0], temp_eff2[i])
                    cav_array2[i] = Sum2
                    mag_list2[i][t] = np.tanh(cav_array2[i] / temp_eff2[i])

                    # Network 3 ----------------
                    Sum3 = err_heat
                    for index in nn_list3[i]: 
                        Sum3 += kernel2sp(cav_array3[index][0], temp_eff3[i])
                    cav_array3[i] = Sum3
                    mag_list3[i][t] = np.tanh(cav_array3[i] / temp_eff3[i])

                mag_array1, mag_array2, mag_array3 = [], [], []
                for i in range(N_pop):
                    mag_array1.append(mag_list1[i][t])
                    mag_array2.append(mag_list2[i][t])
                    mag_array3.append(mag_list3[i][t])

                glob_mag1 = np.mean(mag_array1)
                glob_mag2 = np.mean(mag_array2)
                glob_mag3 = np.mean(mag_array3)

            # -------- Appending converged magnetizations at this temperature -------- 
            mag2spin1_cool.append(glob_mag1)
            mag2spin2_cool.append(glob_mag2)
            mag2spin3_cool.append(glob_mag3)
            print('The layers magnetizations are:', glob_mag1, glob_mag2, glob_mag3)

            if abs(glob_mag1) >= 1e-1:
                break

        # -------- Collecting the cooling threshold at this run --------
        Tsp_cool.append(temp / k_av)

    ''' **************** Saving the thresholds statistics **************** '''

    Tsp_3q4_glob.append( (q, np.mean(Tsp_heat), np.std(Tsp_heat), np.mean(Tsp_cool), np.std(Tsp_cool)) )
    # Saving the datapoints in an external file
    np.savetxt('/Users/ivan_bonamassa/Desktop/Tsp_3q4_th_k4_glob_1run.txt', Tsp_3q4_glob)


print("--- %s seconds ---" % (time.time() - start_time))



'''    ****************************************************************************    '''
'''    -------------  3+q - THERMAL (local) + synthetic fluctuations --------------    '''
'''    ****************************************************************************    '''
# Also here we adopt the mixed time scales. 
# Numerically equivalent to the thermal one with segregated time scales, but O(2N) faster. 

import matplotlib.pyplot as plt
plt.style.use('classic')

import numpy as np
import random as rm
import time


# Kernel functions solving the cavity fields via saddle-point
def kernel2sp(x, y):
    return (y / 2) * np.log((np.cosh((x + 1.0) / y)) / (np.cosh((x - 1.0) / y)))


 # ---------------- Static parameters of the runs ----------------
k_av = 4
N_pop = 10**3
T = 10**2
err_heat = 1e-20 
mag_thr = 1e-15
Tsp_th_3q4_k4 = []
runs = 10**0
qvals = np.arange(0.00, 1.05, 0.05)
temps_th_3q4sp = np.linspace(0.05 * k_av, 0.45 * k_av, 40)


''' ------------- Starting the script --------------- '''
start_time = time.time()

for q in qvals: 

    Tsp_heat = []
    # Setting the values depending on q
    sigma_q = (1.25 + q / 4) / k_av
    
    # Collecting more runs for statistical significance
    for run in range(runs):
    
        ''' --------------------------------------------------------------- '''
        ''' Interdependent solution: mixed time scales, same saddle kernels '''
        ''' --------------------------------------------------------------- '''
        # Numerically equivalent to the thermal one, but O(2N) faster. 
        
        
        # **************** Structural part ****************
        
        N_pop_q = int(q * N_pop)
        net_indices1 = list(range(0, N_pop))
        net_indices2 = rm.sample(net_indices1, N_pop)
        net_indices3 = rm.sample(net_indices2, N_pop)
        int_pairs = [x for x in zip(net_indices1, net_indices2, net_indices3)] 

        deg_list1, deg_list2, deg_list3 = [], [], []
        nn_list1, nn_list2, nn_list3 = [], [], []
        
        for i in range(N_pop):
            
            deg_list1.append(np.random.poisson(k_av, 1)[0])
            nn_list1.append(rm.sample(np.delete(net_indices1, i).tolist(), deg_list1[-1])) 
            
            deg_list2.append(np.random.poisson(k_av, 1)[0])
            nn_list2.append(rm.sample(np.delete(net_indices2, i).tolist(), deg_list2[-1])) 
            
            deg_list3.append(np.random.poisson(k_av, 1)[0])
            nn_list3.append(rm.sample(np.delete(net_indices3, i).tolist(), deg_list3[-1])) 


        ''' **************** Heating part ****************** '''
        
        loc_mag1, loc_mag2, loc_mag3 = np.full(N_pop, 1.0), np.full(N_pop, 1.0), np.full(N_pop, 1.0)
        mag2spin1_heat, mag2spin2_heat, mag2spin3_heat = [], [], []
        mag_list1 = np.zeros((N_pop, len(temps_th_3q4sp)))
        mag_list2 = np.zeros((N_pop, len(temps_th_3q4sp)))
        mag_list3 = np.zeros((N_pop, len(temps_th_3q4sp)))
        cav_array1 = k_av * np.random.random_sample((N_pop, 1)) 
        cav_array2 = k_av * np.random.random_sample((N_pop, 1)) 
        cav_array3 = k_av * np.random.random_sample((N_pop, 1)) 

        t = -1
        for temp in temps_th_3q4sp:
            t += 1
            
            temp_eff1 = np.full(N_pop, temp)
            temp_eff2 = np.full(N_pop, temp)
            temp_eff3 = np.full(N_pop, temp)
            print('Thermal case. q =', q, 'run =', run, '- The reduced heating temperature is:', temp / k_av)
        
            for pditer in range(T):
                
                ''' Generating distribution of effective temperatures --------- '''
                
                # Network 1, Network 2 and Network 3 are directly interdependent to form a closed loop of triples
                for i in net_indices1[:int(N_pop * (1 - q))]:
                    
                    # Network 1 ----------------
                    int_i2 = int_pairs[i][1]
                    if loc_mag2[int_i2] >= mag_thr:
                        temp_eff1[i] = temp / loc_mag2[int_i2]
                    else: 
                         temp_eff1[i] = temp / (loc_mag2[int_i2] + mag_thr)
                    
                    # Network 2 ----------------
                    int_i2 = int_pairs[i][1]
                    int_i3 = int_pairs[i][2]
                    if loc_mag3[int_i3] >= mag_thr:
                        temp_eff2[int_i2] = temp / loc_mag3[int_i3]
                    else: 
                        temp_eff2[int_i2] = temp / (loc_mag3[int_i3] + mag_thr)      
                        
                    # Network 3 ----------------
                    int_i3 = int_pairs[i][2]
                    int_i1 = int_pairs[i][0]
                    if loc_mag1[int_i1] >= mag_thr:
                        temp_eff3[int_i3] = temp / loc_mag1[int_i1]
                    else: 
                        temp_eff3[int_i3] = temp / (loc_mag1[int_i1] + mag_thr) 
                        
                # Partially completing the triples with bi-directional dependencies     
                for i in net_indices1[int(N_pop * (1 - q)):]:
                    
                    # Network 1 ----------------
                    int_i1 = int_pairs[i][0]
                    int_i2 = int_pairs[i][1]
                    int_i3 = int_pairs[i][2]
                    if loc_mag2[int_i2] >= mag_thr and loc_mag3[int_i3] >= mag_thr:
                        temp_eff1[int_i1] = temp / (loc_mag2[int_i2] * loc_mag3[int_i3])
                    elif loc_mag2[int_i2] < mag_thr or loc_mag3[int_i3] < mag_thr:
                        temp_eff1[int_i1] = temp / (loc_mag2[int_i2] * loc_mag3[int_i3] + mag_thr)
                        
                    # Network 2 ----------------
                    int_i1 = int_pairs[i][0]
                    int_i2 = int_pairs[i][1]
                    int_i3 = int_pairs[i][2]
                    if loc_mag1[int_i1] >= mag_thr and loc_mag3[int_i3] >= mag_thr:
                        temp_eff2[int_i2] = temp / (loc_mag1[int_i1] * loc_mag3[int_i3])
                    elif loc_mag1[int_i1] < mag_thr or loc_mag3[int_i3] < mag_thr:
                        temp_eff2[int_i2] = temp / (loc_mag1[int_i1] * loc_mag3[int_i3] + mag_thr) 
                    
                    # Network 3 ----------------
                    int_i1 = int_pairs[i][0]
                    int_i2 = int_pairs[i][1]
                    int_i3 = int_pairs[i][2]
                    if loc_mag1[int_i1] >= mag_thr and loc_mag2[int_i2] >= mag_thr:
                        temp_eff3[int_i3] = temp / (loc_mag1[int_i1] * loc_mag2[int_i2])
                    elif loc_mag1[int_i1] < mag_thr or loc_mag2[int_i2] < mag_thr:
                        temp_eff3[int_i3] = temp / (loc_mag1[int_i1] * loc_mag2[int_i2] + mag_thr)      
                
                # Calculating cavities and local magnetizations
                for i in range(N_pop):
                    
                    # Network 1 ----------------
                    Sum1 = err_heat
                    for index in nn_list1[i]: 
                        Sum1 += kernel2sp(cav_array1[index][0], temp_eff1[i])
                    cav_array1[i] = Sum1
                    mag_list1[i][t] = np.tanh(cav_array1[i] / temp_eff1[i])
                    
                    # Network 2 ----------------
                    Sum2 = err_heat
                    for index in nn_list2[i]: 
                        Sum2 += kernel2sp(cav_array2[index][0], temp_eff2[i])
                    cav_array2[i] = Sum2
                    mag_list2[i][t] = np.tanh(cav_array2[i] / temp_eff2[i])
                    
                    # Network 3 ----------------
                    Sum3 = err_heat
                    for index in nn_list3[i]: 
                        Sum3 += kernel2sp(cav_array3[index][0], temp_eff3[i])
                    cav_array3[i] = Sum3
                    mag_list3[i][t] = np.tanh(cav_array3[i] / temp_eff3[i])
                
                # Calculating the local OPs at this PD's iteration
                for i in range(N_pop):
                    
                    # Network 1 ----------------
                    loc_sum1 = 0
                    for index in nn_list1[i]:
                        loc_sum1 += mag_list1[index][t] / deg_list1[i]
                    loc_mag1[i] = loc_sum1 + np.random.normal(0, sigma_q, 1)
                    
                    # Network 2 ----------------
                    loc_sum2 = 0
                    for index in nn_list2[i]:
                        loc_sum2 += mag_list2[index][t] / deg_list2[i]
                    loc_mag2[i] = loc_sum2 + np.random.normal(0, sigma_q, 1)
                    
                    # Network 3 ----------------
                    loc_sum3 = 0
                    for index in nn_list3[i]:
                        loc_sum3 += mag_list3[index][t] / deg_list3[i]
                    loc_mag3[i] = loc_sum3 + np.random.normal(0, sigma_q, 1)
                    
            mag_array1, mag_array2, mag_array3 = [], [], []
            for i in range(N_pop):
                mag_array1.append(mag_list1[i][t])
                mag_array2.append(mag_list2[i][t])
                mag_array3.append(mag_list3[i][t])
        
            mag2spin1_heat.append(np.mean(mag_array1))
            mag2spin2_heat.append(np.mean(mag_array2))
            mag2spin3_heat.append(np.mean(mag_array3))
            print('The magnetizations are: ', np.mean(mag_array1), np.mean(mag_array2), np.mean(mag_array3))
            
            if mag2spin1_heat[-1] < 1e-1 or mag2spin2_heat[-1] < 1e-1 or mag2spin3_heat[-1] < 1e-1:
                break

        # -------- Collecting the heating threshold at this run --------
        Tsp_heat.append(temp / k_av)
        
    ''' **************** Saving the thresholds statistics **************** '''

    Tsp_th_3q4_k4.append( (q, np.mean(Tsp_heat), np.std(Tsp_heat), 0.01, 0.00) )
    # Saving the datapoints in an external file
    np.savetxt('/Users/ivan_bonamassa/Desktop/Tsp_3q4_th_k4_fluc_1run.txt', Tsp_th_3q4_k4)
            

print("--- %s seconds ---" % (time.time() - start_time))



''' ---------------------------------------------------------- '''
'''   Hyphergraph - Population dynamics for (3+q)-spin model   '''
''' ---------------------------------------------------------- '''
import matplotlib.pyplot as plt
plt.style.use('classic')

import numpy as np
import random as rm
import time


# Integral kernels solving the saddle-point replica symmetric equations
def kernel2sp(x, y):
    return (y / 2) * np.log((np.cosh((x + 1.0) / y)) / (np.cosh((x - 1.0) / y)))

def kernel3sp(x1, x2, y):
    alpha_p = np.cosh((x2 + x1) / y) + np.exp(- 2 / y) * np.cosh((x2 - x1) / y)
    alpha_m = np.cosh((x2 - x1) / y) + np.exp(- 2 / y) * np.cosh((x2 + x1) / y)
    return (y / 2) * np.log(alpha_p / alpha_m)

def kernel4sp(x1, x2, x3, y):
    alpha4sp_p = np.cosh((+ 1 + x3) / y) * np.cosh((x2 + x1) / y) + \
        np.cosh((+ 1 - x3) / y) * np.cosh((x2 - x1) / y)
    alpha4sp_m = np.cosh((- 1 + x3) / y) * np.cosh((x2 + x1) / y) + \
        np.cosh((- 1 - x3) / y) * np.cosh((x2 - x1) / y)
    return (y / 2) * np.log(alpha4sp_p / alpha4sp_m)


# Function grouping list of elements in n-uples
def grouped_list(lst, n):
    if len(lst) == 0:
        return []
    else:
        return list(zip(*[iter(lst)] * n))


# **************** Structure settings ******************
k_av = 4
N_pop = 10**3
q = 0.0

# Generating the list of nearest neighbors according to their interaction pattern
net_indices = (np.arange(N_pop))
nn_3sp_list, nn_4sp_list = [], []

for i in range(N_pop): 
    
    # Each node has a given, Poisson distributed degree
    deg_i = np.random.poisson(k_av, 1)[0]
    
    # A fraction 1-q of its edges are involved in 3-spin interactions
    deg_3sp_i = int(np.floor(deg_i * (1 - q)))
    indices_3sp_minus_i = np.delete(net_indices, i).tolist()
    indices_3sp_i = rm.sample(indices_3sp_minus_i, 2 * deg_3sp_i)
    nn_3sp_list.append(grouped_list(indices_3sp_i, 2))
    
    # The remaining fraction is partecipating to 4-spin interactions
    deg_4sp_i = deg_i - deg_3sp_i
    indices_4sp_minus_3sp_minus_i = [a for a in indices_3sp_minus_i if a not in indices_3sp_i]
    indices_4sp_i = rm.sample(indices_4sp_minus_3sp_minus_i, 3 * deg_4sp_i)
    nn_4sp_list.append(grouped_list(indices_4sp_i, 3))


start_time = time.time()

# **************** Dynamics settings ******************
T = 10**2
cav_array = k_av * np.random.random_sample((N_pop, 1))
temps3q4sp = np.linspace(0.25 * k_av, 0.5 * k_av, 10**2)


# ****************** Heating part **********************
mag3q4sp_heat = []
for temp in temps3q4sp:
    
    print('The heating temperature is:', temp)
    for pditer in range(T):
        for i in range(N_pop):

            Sum = 0
            for index3sp in nn_3sp_list[i]:
                Sum += kernel3sp(cav_array[index3sp[1]], cav_array[index3sp[0]], temp)
            for index4sp in nn_4sp_list[i]:
                Sum += kernel4sp(cav_array[index4sp[2]], cav_array[index4sp[1]], cav_array[index4sp[0]], temp)
            cav_array[i] = Sum
                
    mag_array3q4sp = []
    for i in range(N_pop):
        a = cav_array[i]
        mag_array3q4sp.append(np.tanh(a / temp))

    mag3q4sp_heat.append(np.mean(mag_array3q4sp))
    print('The magnetization is:', mag_array3q4sp[-1])

# ****************** Cooling part **********************
mag3q4sp_cool = []
for temp in temps3q4sp[::-1]:
    
    print('The cooling temperature is:', temp)
    for pditer in range(T):
        for i in range(N_pop):

            Sum = 0
            for index3sp in nn_3sp_list[i]:
                Sum += kernel3sp(cav_array[index3sp[1]], cav_array[index3sp[0]], temp)
            for index4sp in nn_4sp_list[i]:
                Sum += kernel4sp(cav_array[index4sp[2]], cav_array[index4sp[1]], cav_array[index4sp[0]], temp)
            cav_array[i] = Sum
                
    mag_array3q4sp = []
    for i in range(N_pop):
        a = cav_array[i]
        mag_array3q4sp.append(np.tanh(a / temp))

    mag3q4sp_cool.append(np.mean(mag_array3q4sp))
    print('The magnetization is:', mag_array3q4sp[-1])


# Searching for the heating spinodal
for i in range(5, len(temps3q4sp)): 
    if mag3q4sp_heat[i] < 1/2:
        break
Tsp_hy_3q4_heat = round((temps3q4sp[i] + temps3q4sp[i - 1]) / (2 * k_av), 2)
print('The heating spinodal is at:', Tsp_hy_3q4_heat)

# Searching for the cooling threshold
for i in range(5, len(temps3q4sp)): 
    if mag3q4sp_cool[i] > 1/2:
        break
Tsp_hy_3q4_cool = round((temps3q4sp[::-1][i] + temps3q4sp[::-1][i - 1]) / (2 * k_av), 2)
print('The cooling spinodal is at:', Tsp_hy_3q4_cool)

    
# Plotting the results
plt.plot(temps3q4sp / k_av, mag3q4sp_heat, '-', color = 'k', label = '3-spin hypergraph')
plt.plot(temps3q4sp[::-1] / k_av, mag3q4sp_cool, '-', color = 'k')
plt.legend(loc = 'best')
plt.title(' $T_{sp}^{heat} / k_{av} =$ %s, $T_{sp}^{cool} / k_{av} =$ %s' %(Tsp_hy_3q4_heat, Tsp_hy_3q4_cool))
plt.ylim(0.0, 1.0)
plt.xlim(0.0, (k_av + 0.5) / k_av)
plt.show()

print("--- %s seconds ---" % (time.time() - start_time))

    

''' ----------------------------------------------------------------- '''
'''     Fully MF sol. : 2+q and 3+q Ising on complete hypergraphs     '''
''' ----------------------------------------------------------------- '''

# 2+q case ------------------------
accuracy = 1e-10
temp = np.linspace(1e-5, 1.0, 10**4)
qvals_mf = np.linspace(0.0, 1.0, 10**2)

def full_mf_hy_2q3(q, T, x):
    return np.tanh((q * x ** 2 + (1 - q) * x) / T)

Tsp_2q3_hy = []

for q in qvals_mf:
    m_2q3 = []
    for T in temp:
        m1 = 1.0
        error = 1.0

        # Loop until error is small enough
        while error > accuracy:
            m2 = m1
            m1 = full_mf_hy_2q3(q, T, m1)
            error = abs(m1 - m2)
        m_2q3.append(m1)

    for i in range(len(m_2q3)):
        if abs(m_2q3[i]) < 1e-3:
            break
    Tsp_2q3_hy.append(temp[i])

Tsp_2q3_list_hy = []
for i in range(len(qvals_mf)):
    Tsp_2q3_list_hy.append((qvals_mf[i], Tsp_2q3_hy[i]))

np.savetxt('/Users/physivan/Desktop/Tsp_2q3_MF_hy.txt', Tsp_2q3_list_hy)


# 3+q case --------------------------
accuracy = 1e-10
temp = np.linspace(0.4, 0.6, 10**4)
qvals_mf = np.linspace(0.0, 1.0, 10**2)

def full_mf_hy_3q4(q, T, x):
    return np.tanh((q * x ** 3 + (1 - q) * x ** 2) / T)

Tsp_3q4_hy = []

for q in qvals_mf:
    m_3q4 = []
    for T in temp:
        m1 = 1.0
        error = 1.0

        # Loop until error is small enough
        while error > accuracy:
            m2 = m1
            m1 = full_mf_hy_3q4(q, T, m1)
            error = abs(m1 - m2)
        m_3q4.append(m1)

    for i in range(len(m_3q4)):
        if abs(m_3q4[i]) < 1/2:
            break
    Tsp_3q4_hy.append(temp[i])

Tsp_3q4_list_hy = []
for i in range(len(qvals_mf)):
    Tsp_3q4_list_hy.append((qvals_mf[i], Tsp_3q4_hy[i]))

np.savetxt('/Users/physivan/Desktop/Tsp_3q4_MF_hy.txt', Tsp_3q4_list_hy)



''' ----------------------------------------------------------------- '''
'''      Fully MF sol. : 2+q and 3+q interdependent Ising model       '''
''' ----------------------------------------------------------------- '''

# 2+q case, i.e. 2 partially interdependent layers
accuracy = 1e-10
temp = np.linspace(1e-5, 1.0, 10**4)
qvals_mf = np.linspace(0.0, 1.0, 10**2)

def full_mf_2q3(q, T, x):
    return q * np.tanh(x ** 2 / T) + (1 - q) * np.tanh(x / T)

Tsp_2q3_int = []

for q in qvals_mf:
    m_2q3 = []
    for T in temp:
        m1 = 1.0
        error = 1.0

        # Loop until error is small enough
        while error > accuracy:
            m2 = m1
            m1 = full_mf_2q3(q, T, m1)
            error = abs(m1 - m2)
        m_2q3.append(m1)

    for i in range(len(m_2q3)):
        if abs(m_2q3[i]) < 1e-3:
            break
    Tsp_2q3_int.append(temp[i])

Tsp_2q3_list_int = []
for i in range(len(qvals_mf)):
    Tsp_2q3_list_int.append((qvals_mf[i], Tsp_2q3_int[i]))

np.savetxt('/Users/physivan/Desktop/Tsp_2q3_MF_thermal.txt', Tsp_2q3_list_int)


# 3+q case, i.e. 2 interdependent layers and 1 partially interdependent one
accuracy = 1e-10
temp = np.linspace(0.4, 0.6, 10**4)
qvals_mf = np.linspace(0.0, 1.0, 10**2)

def full_mf_3q4(q, T, x):
    return q * np.tanh(x ** 3 / T) + (1 - q) * np.tanh(x ** 2 / T)

Tsp_3q4_int = []

for q in qvals_mf:
    m_3q4 = []
    for T in temp:
        m1 = 1.0
        error = 1.0

        # Loop until error is small enough
        while error > accuracy:
            m2 = m1
            m1 = full_mf_3q4(q, T, m1)
            error = abs(m1 - m2)
        m_3q4.append(m1)

    for i in range(len(m_3q4)):
        if abs(m_3q4[i]) < 1/2:
            break
    Tsp_3q4_int.append(temp[i])

Tsp_3q4_list_int = []
for i in range(len(qvals_mf)):
    Tsp_3q4_list_int.append((qvals_mf[i], Tsp_3q4_int[i]))

np.savetxt('/Users/physivan/Desktop/Tsp_3q4_MF_thermal.txt', Tsp_3q4_list_int)







