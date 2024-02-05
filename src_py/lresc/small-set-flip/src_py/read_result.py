import os.path
import glob


column_names_for_simulation_description = "algo,dv,dc,nv,nc,code_id,p_phys,p_mask"
column_names_for_simulation_results = "no_test,no_success,no_stop"
column_names_for_extras = "p_log"
# DO not use separator = ','
separator = "\t\t"
first_line = column_names_for_simulation_description + separator + column_names_for_simulation_results + separator + column_names_for_extras + '\n'



class Result:
    def __init__(self,algo,dv,dc,nv,nc,code_id,p_phys,p_mask,no_test,no_success,no_stop):
        if not (type(algo) == int or type(dv) == int or\
                type(dc) == int or type(nv) == int or\
                type(nc) == int or type(code_id) == str or\
                type(p_phys) == float or type(p_mask) == float or type(no_test) == int or\
                type(no_success) == int or type(no_stop) == int):
            raise NameError('Bad result format')
        self.algo = algo  # int
        self.dv = dv  # int
        self.dc = dc  # int
        self.nv = nv  # int
        self.nc = nc  # int
        self.code_id = code_id  # str
        self.p_phys = p_phys  # float
        self.p_mask = p_mask
        self.no_test = no_test  # int
        self.no_success = no_success  # int
        self.no_stop = no_stop  # int

# r is an object of the class Result
# A line is a string
def res_to_line(r):
    p_log = r.no_success / r.no_test
    line = str(r.algo) + ',' + str(r.dv) + ',' + str(r.dc) + ',' + str(r.nv) + ',' + str(r.nc) + ',' + r.code_id +',' + str(r.p_phys) + ',' + str(r.p_mask)
    line += separator + str(r.no_test) + ',' + str(r.no_success) + ',' + str(r.no_stop) + ','
    line += separator + str(p_log) + '\n'
    return line

def line_to_res(line):
    tmp = line.strip('\n').split(separator)
    tmp0 = tmp[0].split(',')
    tmp1 = tmp[1].split(',')
    r = Result(int(tmp0[0]),int(tmp0[1]),int(tmp0[2]),int(tmp0[3]),int(tmp0[4]),tmp0[5],float(tmp0[6]),float(tmp0[7]),int(tmp1[0]),int(tmp1[1]),int(tmp1[2]))
    return r

# Creates a file whose lines are stored in lines_list
def create_file(file_name, lines_list):
    lines_list.sort()
    file = open(file_name, 'w')
    file.write(first_line)
    for line in lines_list:
        file.write(line)
    file.close()


    


# Returns lines which match one of the possibilities.
# e.g: algo = [], dv = [4], dc = [5], nv = [40,50] returns lines such that dv = 4 and dc = 5 and (nv = 40 or dc = 50)
def file_to_res_list(file_name,algo=[],dv=[],dc=[],nv=[],nc=[],code_id=[],p_phys=[],p_mask=[]):
    file = open(file_name, 'r')
    if file.readline() != first_line:
        file.close()
        raise NameError('Bad file format')
    line_list = [line_to_res(line) for line in file]
    file.close()
    if len(algo) != 0:
        line_list = [r for r in line_list if r.algo in algo]
    if len(dv) != 0:
        line_list = [r for r in line_list if r.dv in dv]
    if len(dc) != 0:
        line_list = [r for r in line_list if r.dc in dc]
    if len(nv) != 0:
        line_list = [r for r in line_list if r.nv in nv]
    if len(nc) != 0:
        line_list = [r for r in line_list if r.nc in nc]
    if len(code_id) != 0:
        line_list = [r for r in line_list if r.code_id in code_id]
    if len(p_phys) != 0:
        line_list = [r for r in line_list if r.p_phys in p_phys]
    if len(p_mask) != 0:
        line_list = [r for r in line_list if r.p_mask in p_mask]
    return line_list


# r1 and r2 are objects of the class Result
# This function tries to combine the results r1 and r2.
# It is possible to combine r1 and r2 when r1.algo == r2. algo and r1.dv == r2.dv and r1.dc == r2.dv and r1.nv == r2.nv and r1.nc == r2.nc and r1.code_id == r2.code_id and r1.p_phys == r2.p_phys
# It returns None when it 
def combine_res(r1,r2):
    if r1.algo == r2.algo and r1.dv == r2.dv and\
       r1.dc == r2.dc and r1.nv == r2.nv and\
       r1.nc == r2.nc and r1.code_id == r2.code_id and\
       r1.p_phys == r2.p_phys and r1.p_mask == r2.p_mask:
        no_test = r1.no_test + r2.no_test
        no_success = r1.no_success + r2.no_success
        no_stop = r1.no_stop + r2.no_stop
        return Result(r1.algo,r1.dv,r1.dc,r1.nv,r1.nc,r1.code_id,r1.p_phys,r1.p_mask,no_test,no_success,no_stop)
    return None



############# To store the results during simulations #############
# This function adds the result r in the list 'res_list'
# In place function
def add_new_res(res_list, r):    
    done = False
    i = 0
    while not done and i < len(res_list):
        r2 = res_list[i]
        r_new = combine_res(r,r2)
        if r_new == None:
            r_new = r2
        else:
            done = True
        res_list[i] = r_new
        i = i+1
    if not done:
        res_list.append(r)


# This function adds the result r in the file 'file_name'
# We create a tmp file and then rename it. If we had modified directly the file 'file_name' then we lose the result when there is Ctrl-C or a timeout during this function.
def save_new_res(file_name, new_res_list):
    if not os.path.exists(file_name):
        create_file(file_name,[])
    tmp_file_name = file_name + "tmp"
    file = open(file_name, 'r')
    if file.readline() != first_line:
        file.close()
        raise NameError('Bad file format')

    res_list = [line_to_res(line) for line in file]
    file.close()
    for r in new_res_list:
        add_new_res(res_list, r)
    new_lines_list = [res_to_line(r) for r in res_list]
    create_file(tmp_file_name,new_lines_list)
    os.replace(tmp_file_name,file_name)
###########################################################################


############# To create the summary file #############
# BE CAREFULL: this function erases the file summary_file_name
# This function merges the files in directory_name which name finishes by result_extension into the file summary_file_name
def create_summary_file(file_name_list, summary_file_name):
    file_name_list = [f for f in file_name_list if f != summary_file_name]
    res_list = []
    for file_name in file_name_list:
        file = open(file_name,'r')
        if file.readline() != first_line:
            file.close()
            raise NameError('Bad file format')
        for line in file:
            add_new_res(res_list, line_to_res(line))
        file.close()
    line_list = [res_to_line(r) for r in res_list]
    create_file(summary_file_name,line_list)
