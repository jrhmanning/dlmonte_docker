# This script performs multiple DL_MONTE simulatins to calculate the mean and uncertainty in the energy and 
# number of molecules in the system reflected in the CONTROL, CONFIG and FIELD files in the current directory 
# - which correspond to a Lennard-Jones fluid above the critical temperature. Simulations are performed until
# the energy is determined to within an uncertainty of 0.2 energy units.

import logging

import dlmontepython.simtask.dlmonteinterface as interface
import dlmontepython.simtask.measurement as measurement
import dlmontepython.simtask.analysis as analysis
import dlmontepython.simtask.task as task

# Set up the logger, which determines the nature of information output by the machinery in the 'task' package.
# The code below results in logging information being output to stdout

handler = logging.StreamHandler()

# Replacing 'logging.INFO' with 'logging.DEBUGGING' below results in more information being output by the
# logger. Using 'logging.WARNING' results in less information being output: only 'warnings'

measurement.logger.setLevel(logging.INFO)
measurement.logger.addHandler(handler)

# Set up the relevant TaskInterface object: which tells the low-level machinery in the 'task' package
# which code will be used to perform the simulations, and how to perform various tasks specific to that
# code, e.g. extracting the energy from output files created by the code.
# In this case we use DL_MONTE to perform our simulations; thus the TaskInterface object we will use
# is in fact a DLMonteInterface object (DLMonteInterface is a subclass of TaskInterface). The line
# below sets up a DL_MONTE-specific interface. Note that the interface must know the location of the 
# DL_MONTE executable - which is specified as the argument to the DLMonteInterface constructor.

interface = interface.DLMonteInterface("DLMONTE-SRL.X")

# Set up a list of 'observables' to track and analyse. Observables must be Observable objects, and the nature
# of Observable objects may vary between simulation codes. For DL_MONTE only observables corresponding to variables
# output periodically in YAMLDATA are currently supported. For a variable 'foo' specified in the YAMLDATA file
# the corresponding Observable object is returned by the command 'task.Observable( ("foo",) )'. Note the essential
# comma after "foo"! For a variable in YAMLDATA which is an array (e.g., 'nmol'), the observable corresponding to 
# the nth element in the array is returned by the command 'task.Observable( ("foo",n-1) )'. See below: 'energy_obs'
# corresponds to the 'energy' variable in YAMLDATA, and 'nmol_obs' corresponds to the 1st element in the 'nmol'
# array in YAMLDATA (which in fact is the number of molecules belonging to the 1st molecular species)

energy_obs = task.Observable( ("energy",) )
nmol_obs = task.Observable( ("nmol",0) )
observables = [ energy_obs, nmol_obs ]

# Set up a dictionary containing the threshold precisions to which to determine the observables. Here we only
# specify a threshold precision for the energy observable: the simulations will terminate when the energy is
# determined to an uncertainty less than 0.2 (energy units).

precisions = { nmol_obs : 10 }

# Set up the relevant Measurement object - which will actually perform the simulations and data analysis.
# Note that all the simulations and output files pertaining to analysis will be created in the directory
# 'fixedprecision' (via the 'outputdir' argument), and that we will impose a threshold precision on the
# energy as implied by the dictionary 'precisions' described above. Furthermore, we specify that no more 
# than 20 simulations will ever be performed (via the 'maxsims' argument)

m = measurement.Measurement(interface, observables, maxsims=20, precisions=precisions, outputdir='fixedprecision_sweep')

# Run the task

#m.run()

# Once the task is complete in the 'fixedprecision' directory there will be directories 'sim_1', 'sim_2', etc.
# containing the input and output files corresponding to each simulation. Note that 'sim_2' is a resumption
# of 'sim_1', 'sim_3' is a resumption of 'sim_2', etc. Furthermore the files 'energy_converge.dat' and
# 'nmol_0_converge.dat' contain, respectively, the mean and uncertainty in the energy and number of molecules
# after each simulation. Note that the uncertainty decreases as expected after each simulation. Note also that
# no more simulations were performed once the uncertainty in the energy was determined to within 0.2 (units of
# energy).

from copy import deepcopy
def isotherm_func(pressure, template):
    new_meas = deepcopy(template)
    new_meas.interface.amend_input_parameter(new_meas.inputdir, "molchempot", pressure)
    new_meas.outputdir = str(pressure)
    new_meas.run()

    output = (new_meas.mean[nmol_obs], new_meas.stderr[nmol_obs])
    return output


def bin_search(objective, low=(1e-5,-6.9,0), high=(1-1e-5, 11.5,0),res=5):
    mid_x = round(low[0]+((high[0]-low[0])/2), res)
    obj_low = low[1]
    obj_high = high[1]
    obj_mid, err_mid = objective(mid_x, m)
    mid = (mid_x, obj_mid, err_mid)
    low_interval = (obj_mid-obj_low)/(obj_high-obj_low)
    high_interval = (obj_high-obj_mid)/(obj_high-obj_low)
    return (mid, high), (low, mid)

def bin_sweep(max_iters, min_x_spacing, min_y_spacing, low = 1e-5, high = 0.1):
    iteration=0
    min_x = low
    min_y, err_min = isotherm_func(min_x, m)
    
    max_x = high
    max_y, err_max = isotherm_func(max_x, m)
    
    interval = ((min_x, min_y, err_min), (max_x, max_y, err_max))
    interval_buffer = []
    points = [min_x, max_x]
    values = [(min_x, min_y, err_min), (max_x, max_y, err_max)]

    def completion_tests(half_interval, 
                         min_x_spacing=min_x_spacing, 
                         min_y_spacing=min_y_spacing, 
                         max_x=max_x, 
                         max_y=max_y, 
                         min_x=min_x, 
                         min_y=min_y):
        """Tests if the interval should be further subdivided"""
        test1 = ((half_interval[1][0]-half_interval[0][0])/(max_x-min_x))<min_x_spacing # if the minimum x spacing has been exceeded
        test2 = ((half_interval[1][1]-half_interval[0][1])/(max_y-min_y))<min_y_spacing # if the minimum y spacing has been exceeded
        return [test1, test2]
    
    while iteration<max_iters:
        working = bin_search(isotherm_func, *interval)
        iteration +=1

        points.append(working[0][0][0])
        values.append(working[0][0])
        if not any(completion_tests(working[1])):
            interval_buffer.append(working[1])
        if any(completion_tests(working[0])):
            break
        interval = working[0]
    

    while len(interval_buffer)>0:
        interval = interval_buffer[-1]
        interval_buffer.pop(-1)
        while iteration<max_iters:
            working = bin_search(isotherm_func, *interval)
            iteration +=1
            points.append(working[0][0][0])
            values.append(working[0][0])
            if not any(completion_tests(working[1])):
                interval_buffer.append(working[1])
            if any(completion_tests(working[0])):
                break
            interval = working[0]     
    return values 

if __name__ =="__main__":
    output = bin_sweep(50, 0.001, 0.1, 1e-5, 0.1)
    with open('test_sweep.dat', 'w') as f:
        for x in output:
            f.write('{0}, {1}, {2}\n'.format(x[0], x[1], x[2]))
    
