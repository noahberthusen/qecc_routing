import glob
# from decoder import Classical_code


################# Classical code #####################
################# Copy/pasted from decoder.pyx #####################

def is_regular(bit_nbhd, dv):
    """
    Input: 
      - bit_nbhd is defined like in the class Classical_code
      - dv is an integer
    Output:
      Return True when each bit has exactely 'dv' different neighbours
    """
    for i in range(len(bit_nbhd)):
        L = bit_nbhd[i]
        if len(set(L)) < dv:
            return False
    return True

class Classical_code:
    """
    - n is the number of bits
    - m is the number of checknodes
    - bit_nbhd is a size n list of lists. bit_nbhd[no_bit] is the list
    of checknodes which involve the bit number no_bit.
    - check_nbhd is a size m list of lists. check_nbhd[no_check] is the list
    of bits which are involved in the check number no_check.
    """
    def __init__(self, n, m, bit_nbhd, check_nbhd, dv, dc, id, not_regular = False):
        if not not_regular and (not is_regular(bit_nbhd, dv) or not is_regular(check_nbhd, dc)):
            raise NameError('Cannot create non-regular classical code')
        self.n = n
        self.m = m
        # self.bit_nbhd = [list(set(nbhd)) for nbhd in bit_nbhd]
        # self.check_nbhd = [list(set(nbhd)) for nbhd in check_nbhd]
        self.bit_nbhd = bit_nbhd
        self.check_nbhd = check_nbhd
        self.dv = dv
        self.dc = dc
        self.id = id



def write_ccode(file_name, ccode):
    f = open(file_name, 'w+')
    f.write('n,' + str(ccode.n) + '\n')
    f.write('m,' + str(ccode.m) + '\n')
    f.write('dv,' + str(ccode.dv) + '\n')
    f.write('dc,' + str(ccode.dc) + '\n')
    f.write('id,' + str(ccode.id) + '\n')
    f.write('bit_nbhd\n')
    for nbhd in ccode.bit_nbhd:
        for nbr in nbhd:
            f.write(str(nbr))
            f.write(',')
        f.write('\n')
    f.write('check_nbhd\n')
    for nbhd in ccode.check_nbhd:
        for nbr in nbhd:
            f.write(str(nbr))
            f.write(',')
        f.write('\n')
        
    f.close()


# This function returns the list of classical codes in the files 'file_name_list' such that ccode.n in n_list, ccode.m in m_list ...
# If n_list = [] then we don't use the parameter n_list
def read_ccode(file_name_list, n_list, m_list, dv_list, dc_list, id_list):
    res = []
    for file_name in file_name_list:
        file = open(file_name,'r')
        tmp_n = file.readline().strip(',\n').split(',')
        while tmp_n != ['']:
            tmp_m = file.readline().strip(',\n').split(',')
            tmp_dv = file.readline().strip(',\n').split(',')
            tmp_dc = file.readline().strip(',\n').split(',')
            tmp_id = file.readline().strip(',\n').split(',')
            if tmp_n[0] != 'n' or tmp_m[0] != 'm' or tmp_dv[0] != 'dv' or tmp_dc[0] != 'dc' or\
               file.readline() != "bit_nbhd\n":
                print(tmp_n[0], tmp_m[0], tmp_dv[0], tmp_dc[0])
                raise NameError('Bad file format1')
            n = int(tmp_n[1])
            m = int(tmp_m[1])
            dv = int(tmp_dv[1])
            dc = int(tmp_dc[1])
            id = tmp_id[1]
            bit_nbhd = []
            check_nbhd = []
            for j in range(n):
                nbhd = [int(c) for c in file.readline().strip(',\n').split(',')]
                bit_nbhd.append(nbhd)
            if file.readline() != "check_nbhd\n":
                raise NameError('Bad file format2')
            for j in range(m):
                nbhd = [int(v) for v in file.readline().strip(',\n').split(',')]
                check_nbhd.append(nbhd)
            if (n_list == [] or n in n_list) and\
               (m_list == [] or m in m_list) and\
               (dv_list == [] or dv in dv_list) and\
               (dc_list == [] or dc in dc_list) and\
               (id_list == [] or id in id_list):
                res.append(Classical_code(n,m,bit_nbhd,check_nbhd,dv,dc,id))
            tmp_n = file.readline().strip(',\n').split(',')
        file.close()
    return res


        




    
