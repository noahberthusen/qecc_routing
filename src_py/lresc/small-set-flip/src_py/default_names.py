import uuid

# Directory where the classical codes are stored
ccode_dir = "../ccode"
# Directory where the results are stored
res_dir = "../results"

###################################################
# When we generate classical codes
###################################################
# Name of the file where we store the codes for given n,m,dv,dc
def ccode_file_name(n,m,dv,dc,swap):
    if swap == 0:
        return ccode_dir + '/' + str(n) + '_' + str(m) + '_' + str(dv) + '_' + str(dc) + '.code'
    elif swap == 1:
        return ccode_dir + '/swap_' + str(n) + '_' + str(m) + '_' + str(dv) + '_' + str(dc) + '.code'
    elif swap == 2:
        return ccode_dir + '/swap2_' + str(n) + '_' + str(m) + '_' + str(dv) + '_' + str(dc) + '.code'
    else:
        return ccode_dir + '/swap3_' + str(n) + '_' + str(m) + '_' + str(dv) + '_' + str(dc) + '.code'
    
###################################################
# When we run simulations
###################################################
# File where the result is stored
id = str(uuid.uuid4())[:8]
res_file_name = res_dir + "/py_laptop_" + id + ".res"
cluster_size_file_name = res_dir + "/cluster_size_" + id + ".res"
error_size_file_name = res_dir + "/error_size_" + id + ".res"

###################################################
# When we create a summary of the results
###################################################
# File which summarizes the results
summary_res_file = "all_results.res"
