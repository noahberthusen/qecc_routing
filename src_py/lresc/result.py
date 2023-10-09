import os.path

column_names_for_simulation = "t,m,l,a1,a2,a3,b1,b2,b3,emb_m,emb_l,aa,bb,p_phys,p_mask,no_test,no_success,p_log"
first_line = column_names_for_simulation + '\n'

class Result:
    def __init__(self,t,m,l,a1,a2,a3,b1,b2,b3,emb_m,emb_l,aa,bb,p_phys,p_mask,no_test,no_success):
        if not (type(t) == int or type(m) == int or type(l) == int or\
                type(a1) == int or type(a2) == int or type(a3) == int or type(b1) == int or type(b2) == int or type(b3) == int or\
                type(emb_m) == int or type(emb_l) == int or type(aa) == int or type(bb) == int or\
                type(p_phys) == float or type(p_mask) == float or type(no_test) == int or type(no_success) == int):
            raise NameError('Bad result format')
        self.t = t
        self.m = m
        self.l = l
        self.a1 = a1
        self.a2 = a2
        self.a3 = a3
        self.b1 = b1
        self.b2 = b2
        self.b3 = b3
        self.emb_m = emb_m
        self.emb_l = emb_l
        self.aa = aa
        self.bb = bb
        self.p_phys = p_phys
        self.p_mask = p_mask
        self.no_test = no_test
        self.no_success = no_success

def res_to_line(r):
    p_log = r.no_success / r.no_test
    line = str(r.t) + ',' + str(r.m) + ',' + str(r.l) + ',' +\
           str(r.a1) + ',' + str(r.a2) + ',' + str(r.a3) + ',' + str(r.b1) + ',' + str(r.b2) + ',' + str(r.b3) + ',' +\
           str(r.emb_m) + ',' + str(r.emb_l) + ',' + str(r.aa) + ',' + str(r.bb) + ',' +\
           str(r.p_phys) + ',' + str(r.p_mask) + "," + str(r.no_test) + ',' + str(r.no_success) + ',' + str(p_log) + '\n'
    return line

def line_to_res(line):
    tmp = line.strip('\n').split(',')
    r = Result(int(tmp[0]),int(tmp[1]),int(tmp[2]),
               int(tmp[3]),int(tmp[4]),int(tmp[5]),int(tmp[6]),int(tmp[7]),int(tmp[8]),
               int(tmp[9]),int(tmp[10]),int(tmp[11]),int(tmp[12]),
               float(tmp[13]),float(tmp[14]),int(tmp[15]),int(tmp[16]))
    return r

# Creates a file whose lines are stored in lines_list
def create_file(file_name, lines_list):
    lines_list.sort()
    file = open(file_name, 'w')
    file.write(first_line)
    for line in lines_list:
        file.write(line)
    file.close()


# r1 and r2 are objects of the class Result
# This function tries to combine the results r1 and r2.
# It is possible to combine r1 and r2 when r1.algo == r2. algo and r1.dv == r2.dv and r1.dc == r2.dv and r1.nv == r2.nv and r1.nc == r2.nc and r1.code_id == r2.code_id and r1.p_phys == r2.p_phys
# It returns None when it
def combine_res(r1,r2):
    if (r1.t == r2.t and r1.m == r2.m and r1.l == r2.l and\
        r1.a1 == r2.a1 and r1.a2 == r2.a2 and r1.a3 == r2.a3 and r1.b1 == r2.b1 and r1.b2 == r2.b2 and r1.b3 == r2.b3 and\
        r1.emb_m == r2.emb_m and r1.emb_l == r2.emb_l and r1.aa == r2.aa and r1.bb == r2.bb and\
        r1.p_phys == r2.p_phys and r1.p_mask == r2.p_mask):
        no_test = r1.no_test + r2.no_test
        no_success = r1.no_success + r2.no_success
        return Result(r1.t,r1.m,r1.l,
                      r1.a1,r1.a2,r1.a3,r1.b1,r1.b2,r1.b3,
                      r1.emb_m,r1.emb_l,r1.aa,r1.bb,
                      r1.p_phys,r1.p_mask,no_test,no_success)
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
    tmp_file_name = file_name + ".tmp"
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
