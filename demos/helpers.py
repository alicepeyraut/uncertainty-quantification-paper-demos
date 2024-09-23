import numpy


def write_porosity(porosity_field = [], n_cells = 0, filepath = "./"):
    print("filepath", filepath)
    with open(filepath, "w") as file:
        file.write('<?xml version="1.0" encoding="UTF-8"?>\n')
        file.write('<dolfin xmlns:dolfin="http://fenicsproject.org">\n')
        file.write('  <mesh_function type="double" dim="'+str(3)+'" size="'+str(n_cells)+'">\n')
        for k_cell in range(n_cells):
            file.write('    <entity index="'+str(k_cell)+'" value="'+str(porosity_field[k_cell])+'"/>\n')
        file.write('  </mesh_function>\n')
        file.write('</dolfin>\n')
        file.close()


def checking_if_processes_converged(processes):
    if not processes:  # Check if the list is empty
        print("Warning: No processes found.")
        return True 
    over = True
    for process in processes:
        if process.poll() is None:  # If any process is still running
            over = False  # Not all processes are done
            break 
    return(over)


def checking_if_converged(params_opt = {}, nb_processes_converged = 0, tol = 1e-3):
    converged =False
    crit_max = []
    if nb_processes_converged == 0:
        nb_processes_converged = 1 
    for list_ in params_opt:
        l = list_[1:]
        if abs(numpy.std(l)) > 45:
            diverge = True
            return(converged, None)
        else:
            if numpy.percentile(l[:-nb_processes_converged], 50) == 0:
                median = abs((numpy.percentile(l[:],50)-numpy.percentile(l[:-nb_processes_converged],50)))
            else:
                median = abs((numpy.percentile(l[:],50)-numpy.percentile(l[:-nb_processes_converged],50))/numpy.percentile(l[:-nb_processes_converged],50))
            if numpy.percentile(l[:-nb_processes_converged],25) == 0:
                q1 = abs((numpy.percentile(l[:],25)-numpy.percentile(l[:-nb_processes_converged],25)))
            else:
                q1 = abs((numpy.percentile(l[:],25)-numpy.percentile(l[:-nb_processes_converged],25))/numpy.percentile(l[:-nb_processes_converged],25))
            if numpy.percentile(l[:-nb_processes_converged],25) == 0:
                q3 = abs((numpy.percentile(l[:],75)-numpy.percentile(l[:-nb_processes_converged],75)))
            else:
                q3 = abs((numpy.percentile(l[:],75)-numpy.percentile(l[:-nb_processes_converged],75))/numpy.percentile(l[:-nb_processes_converged], 75))
            ##### to modify, but tempor fix with problem with e-16 values as it is ==> std 1 while it is between 2 numerical 0
            if numpy.std(l[:-nb_processes_converged])==0:
                std_plus = abs(numpy.percentile(l[:],50)+numpy.std(l[:])-numpy.percentile(l[:-nb_processes_converged],50)-numpy.std(l[:-nb_processes_converged]))
                std_minus = abs(numpy.percentile(l[:],50)-numpy.std(l[:])-numpy.percentile(l[:-nb_processes_converged],50)-numpy.std(l[:-nb_processes_converged]))
            else:
                std_plus = abs((numpy.percentile(l[:],50) + numpy.std(l[:]) - numpy.std(l[:-nb_processes_converged])-numpy.percentile(l[:-nb_processes_converged],50))/(numpy.std(l[:-nb_processes_converged])+numpy.percentile(l[:-nb_processes_converged],50)))
                std_minus = abs((numpy.percentile(l[:],50) - numpy.std(l[:]) + numpy.std(l[:-nb_processes_converged])-numpy.percentile(l[:-nb_processes_converged],50))/(-numpy.std(l[:-nb_processes_converged])+numpy.percentile(l[:-nb_processes_converged],50)))
            std = max(std_plus, std_minus)
            # print("average for list", l, numpy.average(l[:]), numpy.average(l[:-1]))
            if numpy.average(l[:-nb_processes_converged]) !=0 :
                average = abs((numpy.average(l[:])-numpy.average(l[:-nb_processes_converged]))/numpy.average(l[:-nb_processes_converged]))
            else:
                average = abs((numpy.average(l[:])-numpy.average(l[:-nb_processes_converged])))
            crit_max.append(median)
            crit_max.append(q1)
            crit_max.append(q3)
            crit_max.append(std)
            crit_max.append(average)   
        
        criteria = max(crit_max)
    if criteria < float(tol):
        converged = True
    return(converged, criteria)
