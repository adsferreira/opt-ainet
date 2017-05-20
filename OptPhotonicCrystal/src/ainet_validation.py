import numpy as np

class OptAiNetUtils:
    def join_attr(self, matrix1, matrix2):
        nr_of_attr_mat1 = len(matrix1[0])
        nr_of_attr_mat2 = len(matrix2[0])
        nr_of_attr = nr_of_attr_mat1 + nr_of_attr_mat2
        matrix = [[float(0.) for _ in range(nr_of_attr)] for _ in range(len(matrix1))]
        
        for i in range(len(matrix1)):
            for j in range(nr_of_attr):
                if (j < (nr_of_attr_mat1)):
                    matrix[i][j] = matrix1[i][j]
                else:
                    matrix[i][j] = matrix2[i][j - nr_of_attr_mat1]
        
        return matrix
    
    def split_attr(self, matrix, nr_of_attr_mat1, nr_of_attr_mat2):
        #nr_of_el = len(matrix)
        #print matrix
        matrix1 = matrix[0:nr_of_attr_mat1]
        #matrix1 = matrix[0:nr_of_el, 0:nr_of_attr_mat1]
#         print "\n\n"
#         print matrix1
        matrix2 = matrix[nr_of_attr_mat1:]
        #matrix2 = matrix[0:nr_of_el, nr_of_attr_mat1:]
#         print "\n\n"
#         print matrix2
        
#         nr_of_attr = nr_of_attr_mat1 + nr_of_attr_mat2
#         matrix1 = [[float(0.) for _ in range(nr_of_attr_mat1)] for _ in range(len(matrix))]
#         matrix2 = [[float(0.) for _ in range(nr_of_attr_mat2)] for _ in range(len(matrix))]
#         
#         for i in range(len(matrix)):
#             for j in range(nr_of_attr):
#                 if (j < nr_of_attr_mat1):
#                     matrix1[i][j] = matrix[i][j]
#                 else:
#                     matrix2[i][j - nr_of_attr_mat1] = matrix[i][j]
        
        return (matrix1, matrix2)

    def grep(self, filename, pattern_center, pattern_radius):
        pos_p_center = []
        pos_p_radius = []
         
        for n,line in enumerate(open(filename)):
            if pattern_center in line:
                pos_p_center.append(n)
            if pattern_radius in line:
                pos_p_radius.append(n)
     
        return (pos_p_center, pos_p_radius)
 
    def write(self, filename, lines_centers, lines_radius, centers, radii):
        unit_cell_centers = np.reshape(centers, (len(radii), 2))
        
        with open(filename, 'r') as c_file:
            # read a list of lines into data
            data = c_file.readlines()
     
            # access the lines that contains the center parameter and set it
            for i in range(len(lines_centers)):
                data[lines_centers[i]] = '                       (center ' + str(unit_cell_centers[i][0]) + ' ' + \
                                                                             str(unit_cell_centers[i][1]) + ' ' + \
                                                                             '0.0) \n'
     
            # access the lines that contains the radius parameter and set it
            for i in range(len(lines_radius)):
                data[lines_radius[i]] = '                       (radius ' + str(radii[i]) + ') \n'                                                                 
        
        # and write everything back
        with open(filename, 'w') as c_file:
            c_file.writelines(data)
            
    def are_rods_into_cell(self, centers, radii, x_min, x_max, y_min, y_max, dim):
        for i in range(len(radii)):
            current_x = centers[i * dim]
            current_y = centers[i * dim + 1]
            current_radius = radii[i]
            
            # validate left and right bounds
            if (current_x - current_radius < x_min) or (current_x + current_radius > x_max):
                return False
            # validate bottom and top bounds
            if (current_y - current_radius < y_min) or (current_y + current_radius > y_max):
                return False
                
        return True
        
    def is_inside_bounds(self, Ab, x_min, x_max, y_min, y_max, nr_of_centers, dim):
        for i in range(nr_of_centers):
            if (Ab[i * dim] < x_min) or (Ab[i * dim] > x_max):
                return False
            if (Ab[(i * dim) + 1] < y_min) or (Ab[(i * dim) + 1] > y_max):
                return False 
        
        # split centers and radii data
        centers = Ab[0:(nr_of_centers * dim)]
        radii = Ab[(nr_of_centers * dim):]
        # validate if rods are inside the unit cell
        if self.are_rods_into_cell(centers, radii, x_min, x_max, y_min, y_max, dim):
            return True
        else:
            return False
                
           
            
                
        

