from ainet_validation import OptAiNetUtils
from multiprocessing import Process, Pipe
from itertools import izip
from time import time, strftime, localtime
import subprocess

import numpy as np


def spawn(f):
    def fun(pipe, x):
        pipe.send(f(x))
        pipe.close()
    return fun

def parmap(f, X):
    pipe = [Pipe() for x in X]
    proc = [Process(target = spawn(f), args=(c, x)) for x,(p, c) in izip(X, pipe)]
    [p.start() for p in proc]
    [p.join() for p in proc]
    return [p.recv() for (p, c) in pipe]

def write(Ab, fit):
    f = open("antibodies", "w") 
    f.write("Ab\n\n")
    f.write(' '.join(str(e) for e in Ab))
    f.write("\n")
    f.write("fit\n\n")
    f.write(' '.join(str(e) for e in fit))

class OptAiNet6:
    def __init__(self, nr_of_attr, _Ab, nr_of_Ab, nr_of_clones, gen, beta, tsup, tgen, pinc, x_bounds, y_bounds, nr_of_centers, dim):
        self.nr_of_attr = nr_of_attr
        self.nr_of_Ab = nr_of_Ab
        self.nr_of_clones = nr_of_clones
        self.gen = gen  
        self.beta = beta
        self.tsup = tsup  
        self.tgen = tgen
        self.nr_of_sup = 0
        self.sum_eucl_dist = 0
        self.x_min = x_bounds[0]
        self.x_max = x_bounds[1]
        self.y_min = y_bounds[0]
        self.y_max = y_bounds[1]
        self.nr_of_centers = nr_of_centers
        self.dim = dim
        self.fit = []
        self.opt_ainet_utils = OptAiNetUtils()
                
        if len(_Ab) != 0:
            self.Ab = _Ab
            # fitness of initial population of antibodies
            ids = np.array([[i for i in range(len(self.Ab))]])
            Ab_t = np.insert(self.Ab, 0, ids, axis = 1)
            
            if __name__ == '__main__':                    
                self.fit = parmap(self.fitness, Ab_t)
                   
            self.print_ab_fit(self.Ab, self.fit)                        
        else:
            self.Ab = np.empty([nr_of_Ab, nr_of_attr])
           
            while True:
                # initial population of antibodies Ab
                self.init_antibodies() 
                # fitness of initial population of antibodies
                ids = np.array([[i for i in range(len(self.Ab))]])
                Ab_t = np.insert(self.Ab, 0, ids, axis = 1)
                
                if __name__ == '__main__':                    
                    self.fit = parmap(self.fitness, Ab_t)
                   
                self.print_ab_fit(self.Ab, self.fit)
                # validate if antibodies have fitness different zero
                # we need a maximum fit greater than zero because it will be used as a denominator
                # fitness lesser or equals 1.0 may represent false positives band gaps
                if max(self.fit) > 0.0:
                    break
                
    def init_antibodies(self):
        # for the problem of photonic crystal, we have a set of centers and related radii of a unit cell
        # centers and radii have different ranges in initialization, so we need to set them separately
        # center is represented in 2D
        for i in range(self.nr_of_Ab):
            # all antibodies must respect the bounds of domain
            while True:
                centers = np.random.uniform(self.x_min, self.x_max, [self.nr_of_centers * self.dim])
                radii = np.random.uniform(0, self.x_max, [self.nr_of_centers])
                self.Ab[i] = np.concatenate((centers, radii))               
                
                is_inside_bounds = self.opt_ainet_utils.is_inside_bounds(self.Ab[i], \
                                                                         self.x_min, \
                                                                         self.x_max, \
                                                                         self.y_min, \
                                                                         self.y_max, \
                                                                         self.nr_of_centers, \
                                                                         self.dim)
                if (is_inside_bounds):                    
                    break
        
    def fitness(self, Ab):
        current_bands = []       
        tid = int(Ab[0])
        filename = "Pbg" + str(self.nr_of_centers) + "/photonic_crystal" + str(self.nr_of_centers) + "_" + str(tid) + ".ctl"     
        lines_centers, lines_radius = self.opt_ainet_utils.grep(filename, "center", "radius") 
        list_of_centers, list_of_radii = self.opt_ainet_utils.split_attr(Ab[1:], (self.nr_of_centers * self.dim), self.nr_of_centers)  
        is_rod = True
            
        for j in range(self.nr_of_centers):
            if (list_of_radii[j] <= 0.0):
                is_rod = False
         
        if (is_rod):
            # before run MPB, write the corresponding file
            self.opt_ainet_utils.write(filename, lines_centers, lines_radius, list_of_centers, list_of_radii)
            # call MPB
            p_mpb = subprocess.Popen(['mpb', filename], stdout = subprocess.PIPE)
            out_mpb, mpb_err = p_mpb.communicate()
                                
            # search for a 'Gap' word in mpb's output
            if "Gap" in out_mpb:
                v_out = out_mpb.split("Gap")
                #print v_out
                                         
                for j in range(len(v_out) - 1):
                    current_str = v_out[j + 1] 
                    begining = current_str.find(", ") + 2
                    end = current_str.find("%")
                    current_band = float(current_str[begining:end])
                    current_bands.append(current_band)
            else:
                current_bands.append(0.)   
                     
            # get the max band           
            band = max(current_bands)                 
        else:
            band = 0
                
        return band
    
    def print_ab_fit(self, Ab, fit):
        for i in range(len(Ab)):
            print "cell " + str(i) + ": " + str(Ab[i]) + " - fit: " + str(fit[i])
         
    def clone_mut_select(self, Ab, fit):
        # clones of memory cells
        Abc = np.empty([len(Ab), self.nr_of_clones + 1, self.nr_of_attr])
        # best clones for each cell
        C = np.empty([len(Ab), self.nr_of_attr])
        # fitness of clones from each cell
        fit_of_clones = np.empty([self.nr_of_clones + 1])
        # best clone's fitness for each cell
        best_fit_of_clones = np.empty([len(Ab)])
        # label clones with fitness equals zero
        mask = np.ones(len(Ab), dtype = bool)
        # number of 'bad' mutations
        nr_of_mut = 0
        # threshold of the number of 'bad' mutations
        tr = 50
        
        for i in range(len(Ab)):
            # first "clone" is the current antibody
            Abc[i][0] = Ab[i]
            
            print "----------------------------------------------------------------"
            print "cloning original antibody number: ", i, "...\n"         
            # counter of clones    
            k = 1
            
            while (k <= self.nr_of_clones):                
                _beta = self.beta
                # all clones must respect the bounds of domain
                while True:     
                    # change mutation function when a clone doesn't respect the boundary limits
                    if (nr_of_mut > tr):                        
                        _beta += tr
                        tr = tr + 50
                    
                    mut = np.random.randn(self.nr_of_attr) / _beta * np.exp(-(fit[i] / max(fit)))
                    #print "mut: ", mut
                    # print mut
                    # validate if the clones respect the bounds of solution
                    Ab_mut = Ab[i] + mut
                    is_inside_bounds = self.opt_ainet_utils.is_inside_bounds(Ab_mut, \
                                                                             self.x_min, \
                                                                             self.x_max, \
                                                                             self.y_min, \
                                                                             self.y_max, \
                                                                             self.nr_of_centers, \
                                                                             self.dim)
                      
                    nr_of_mut += 1
                                      
                    if is_inside_bounds:
                        Abc[i][k] = Ab_mut   
                        break                                   
                
                k = k + 1
                                
            # fitness of current cell and its clones
            fit_of_clones[0] = fit[i]
            ids = np.array([[j for j in range(len(Abc[i][1:]))]])
            Ab_t = np.insert(Abc[i][1:], 0, ids, axis = 1)
               
            if __name__ == '__main__':
                fit_of_clones[1:] = parmap(self.fitness, Ab_t)        
            
            self.print_ab_fit(Abc[i], fit_of_clones)
            print "----------------------------------------------------------------"
            
            if max(fit_of_clones > 0.):
                # id of the best fit
                best_fit_of_clones_id = np.argmax(fit_of_clones)
                # store the highest fitness
                best_fit_of_clones[i] = fit_of_clones[best_fit_of_clones_id]                
                # store the cell with the highest fitness
                C[i] = Abc[i][best_fit_of_clones_id]
                mask[i] = True
            else:
                mask[i] = False
                    
        C = C[mask]
        best_fit_of_clones = best_fit_of_clones[mask]
        
        return (C, best_fit_of_clones)
    
    def supress(self, Ab, fit):
        # indicate which antibodies must remain
        mask = np.ones(len(Ab), dtype = bool)
        sum_euclid_dist = 0
        nr_of_dist = 0
                
        for i in range(len(Ab) - 1):
            k = i + 1     
                            
            while k < len(Ab):       
                sum_square = np.sum((Ab[i] - Ab[k]) ** 2)                                 
                # compute euclidian distance    
                eucl_dist = np.sqrt(sum_square)
                sum_euclid_dist += eucl_dist   
                nr_of_dist += 1 
                print "-------\nEuclidean distance between cell", i , " e ", k, ": ", eucl_dist, "\n-------"
                
                if eucl_dist < self.tsup:
                    if fit[i] <= fit[k]:
                        mask[i] = False
                    else:
                        mask[k] = False
                
                k = k + 1
    
        print "mean of Euclidean distance is: ", (sum_euclid_dist / nr_of_dist)
        self.sum_eucl_dist += (sum_euclid_dist / nr_of_dist)
        Ab1 = Ab[mask]
               
        return (Ab1, mask)
    
    def search(self):
        fit_after_supress = []       
        memory_length = 0
        avg_fit_old = np.mean(self.fit)
        Ab_new = self.Ab
                
        for i in range(self.gen):
            print "----------------------------------------------------------------"
            print "Generation ", i + 1
            print "----------------------------------------------------------------"
            C, fit_C = self.clone_mut_select(Ab_new, self.fit)
            Ab_new = C
            self.fit = fit_C
            avg_fit = np.mean(self.fit)   
            
            print "----------------------------------------------------------------"
            print "mean fitness of generation ", i + 1, " after cloning process: ", avg_fit
            print "----------------------------------------------------------------"
            
            
            if abs(1 - (avg_fit_old / avg_fit)) < self.tgen:    
                self.nr_of_sup += 1             
                print "----------------------------------------------------------------"
                print "Suppression Process number: ", self.nr_of_sup
                print "----------------------------------------------------------------"
                print "length of current network cell: ", len(Ab_new)
                    
                Ab_new, mask = self.supress(Ab_new, self.fit)
                    
                print "length of current network cell after suppression: ", len(Ab_new)
                print "memory length of previous suppression process (", self.nr_of_sup - 1, "): ", memory_length
                print "----------------------------------------------------------------"
                print "mean fitness of generation ", i + 1, " after suppression process: ", sum(self.fit[mask]) / len(self.fit[mask])
                print "----------------------------------------------------------------"
                    
                if memory_length == len(Ab_new):                    
                    if (int(sum(self.fit[mask])) == int(sum(fit_after_supress))) and (i > 450):
                        print "----------------------------------------------------------------"
                        print "stabilized at generation: ", i + 1
                        print "mean of last suppression fitness: ", sum(fit_after_supress) / len(fit_after_supress)
                        print "mean of last fitness: ", sum(self.fit[mask]) / len(self.fit[mask])
                        print "----------------------------------------------------------------"
                        break
                       
                memory_length = len(Ab_new)
                fit_after_supress = self.fit[mask]            
                # older antibodies are joined to a number 'ninc' of new incomers
                #self.ninc = int(round(pinc * memory_length))
                self.ninc = int(round(pinc * self.nr_of_Ab))
                Ab_new_aux = np.empty([len(Ab_new) + self.ninc, self.nr_of_attr])   
                # copy current memory
                Ab_new_aux[0:len(Ab_new), 0:self.nr_of_attr] = Ab_new
                             
                for j in range(self.ninc):
                    while True:
                        centers = np.random.uniform(self.x_min, self.x_max, [self.nr_of_centers * self.dim])
                        radii = np.random.uniform(0, self.x_max, [self.nr_of_centers])
                        Ab_new_aux[len(Ab_new) + j] = np.concatenate((centers, radii))                                
                        is_inside_bounds = self.opt_ainet_utils.is_inside_bounds(Ab_new_aux[len(Ab_new) + j], \
                                                                                    self.x_min, \
                                                                                    self.x_max, \
                                                                                    self.y_min, \
                                                                                    self.y_max, \
                                                                                    self.nr_of_centers, \
                                                                                    self.dim)
                
                        if (is_inside_bounds):
                            break
                                              
                Ab_new = Ab_new_aux                   
                # fitness of remaining antibodies
                self.fit = fit_after_supress
                ids = np.array([[k for k in range(len(Ab_new[len(fit_after_supress):]))]])
                Ab_t = np.insert(Ab_new[len(fit_after_supress):], 0, ids, axis = 1)
                # fitness of incomers antibodies
                if __name__ == '__main__':
                    new_fit = parmap(self.fitness, Ab_t)
                                       
                self.fit = np.append(self.fit, new_fit)
                avg_fit = np.mean(self.fit)
                self.print_ab_fit(Ab_new, self.fit)
                
            avg_fit_old = avg_fit
        
        self.Ab = Ab_new
        print "----------------------------------------------------------------"
        print "suppression process was executed ", self.nr_of_sup, " times."
        if self.nr_of_sup > 0:
            print "Euclidean mean during optimization process was: ", self.sum_eucl_dist / self.nr_of_sup
        print "----------------------------------------------------------------"

 
# _Ab = [[-0.14096478, -0.07743338,  0.14285484,  0.20504023,  0.34335043,  0.21850687], \
#        [0.01974602,  0.1720212,   0.15539254, -0.15287318,  0.19448198,  0.25977011], \
#        [0.00938384, -0.03623394,  0.25962439, -0.46875238,  0.21405318,  0.00303153], \
#        [0.33513698, -0.38779207,  0.03963977, -0.17961856,  0.06515107,  0.29420122], \
#        [-0.20059415, -0.07934438,  0.38963888, -0.16390214,  0.13721494,  0.10668681], \
#        [0.098838, -0.096287, -0.216067, -0.403265, 0.210949, 0.074313], \
#        [-0.182345, 0.075337, 0.163385, 0.083506, 0.219030, 0.045193], \
#        [-0.119950, 0.281852, -0.235096, 0.084128, 0.115204, 0.160295], \
#        [0.080398, 0.017353, -0.170057, 0.196013, 0.146961, 0.187795], \
#        [-0.229475, -0.166715, 0.030983, 0.248050, 0.138090, 0.248451]
#        ] 

# _Ab = [[0.28667414, -0.21824958, -0.25922606, -0.11682518, 0.21632138, 0.1247642, 0.12247624, 0.20736425, 0.21181645], \
#        [-0.23687193, -0.32313312, -0.20690953, 0.11754077, 0.35803307, -0.10385738, 0.13303399, 0.13269204, 0.14154781], \
#        [-0.17902514, -0.24293117, 0.15824563, 0.13079739, 0.41458493, 0.33868087, 0.18005933, 0.26464033, 0.06706045], \
#        [-0.31013553,  0.2120338, 0.00316904, 0.13382456, -0.21835986,-0.29682483, 0.11394526, 0.21376525, 0.19394304], \
#        [ 0.26079079,  0.26013298, 0.13573794, 0.105514, -0.36689479,  0.10750183, 0.23878677, 0.23337446, 0.13239466], \
#        [2.41743765e-01, 3.60568698e-01, -2.59895837e-01, -3.62518238e-01, 1.92964811e-01, -1.10393442e-04, 1.39056063e-01, 1.37313457e-01, 4.73628964e-02], \
#        [-0.27227112, -0.13636523,  0.29239915, -0.36206272, -0.22426166,  0.10605274, 0.11751039,  0.07050917,  0.17773942], \
#        [-0.06378043, 0.07362566, -0.22193323, -0.38736242, -0.18076471, -0.19641157, 0.2346229,   0.02049827,  0.18381986], \
#        [ 0.16271414, 0.05075439,  0.01696329, 0.13868631, 0.29888251, -0.27496444, 0.26593269,  0.13858041,  0.16695358], \
#        [ 0.19689689, 0.12897751, 0.22173313, 0.37919621, 0.24991474, 0.13411761, 0.21434954, 0.02398293, 0.06555841], \
#        [-0.14017874, -0.18409463,  0.11549039, -0.26757377, 0.10079597, 0.09800457, 0.27034218, 0.02449276, 0.37352072], \
#        [-0.36897179, -0.36121485, -0.03719754, -0.40785559, 0.24820732, -0.13162016, 0.05788878, 0.06061621, 0.11244832], \
#        [-0.24408482,  0.27187135, -0.17881624, 0.10086204, 0.03785726, 0.15115122, 0.08148387, 0.203962,  0.08809029]]
  
nr_of_antibodies = 5
print "initial number of antibodies: ", nr_of_antibodies
nr_of_clones = 5
print "number of clones: ", nr_of_clones
gen = 2
print "number of generations: ", gen
beta = 100
print "beta values: ", beta
tsup = 0.5
print "supression's threshold: ", tsup 
tgen = 0.001
print "stabilization value's threshold", tgen
pinc = 0.4
print "percent of new incomers: ", pinc
x_bounds = [-0.5, 0.5]
print "bounds in x-axis: ", x_bounds 
y_bounds = x_bounds
print "bounds in y-axis: ", y_bounds
nr_of_centers = 2
print "number of rods: ", nr_of_centers
dimension = 2
nr_of_attributes = (nr_of_centers * dimension) + nr_of_centers
print "number of attributes: ", nr_of_attributes
print "-----------------------------------------"

_Ab = []
#_Ab = np.asarray(_Ab)
# sum_euclid_dist = 0
# nr_of_dist = 0
#  
# for i in range(len(_Ab) - 1):
#     k = i + 1
#     while k < len(_Ab):       
#         sum_square = np.sum((_Ab[i] - _Ab[k]) ** 2)                                 
#         # compute euclidian distance    
#         eucl_dist = np.sqrt(sum_square)
#         sum_euclid_dist += eucl_dist   
#         nr_of_dist += 1 
#         print "-------\nEuclidean distance between cell", i , " e ", k, ": ", eucl_dist, "\n-------"
#         k = k + 1

print strftime("Started on: %a, %d %b %Y %H:%M:%S", localtime())
print "-----------------------------------------"
start = time()
opt_ainet = OptAiNet6(nr_of_attributes, _Ab, nr_of_antibodies, nr_of_clones, gen, beta, tsup, tgen, pinc, x_bounds, y_bounds, nr_of_centers, dimension)
opt_ainet.search()
print "-----------------------------------------"
print "\nElapsed runtime: %f" % (time() - start)
print "-----------------------------------------"
print strftime("Finished on: %a, %d %b %Y %H:%M:%S", localtime())
# write final results
#write(opt_ainet.Ab, opt_ainet.fit)
