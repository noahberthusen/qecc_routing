import stim
import numpy as np

def lr_bell_pair(paths, p):  # THIS FUNCTION IS WRONG FOR ODD PATH LENGTHS, doesn't matter. Only even
    c = stim.Circuit()

    for path in paths:
        size = len(path)
        c.append("H", path[:-1][::2])
        c.append("DEPOLARIZE1", path[:-1][(size%2)::2], p/10)
    c.append("TICK")

    for path in paths:
        size = len(path)
        c.append("CNOT", path[:size-(size%2)])
        c.append("DEPOLARIZE2", path[:size-(size%2)], p)
    c.append("TICK")

    for path in paths:
        size = len(path)
        c.append("CNOT", path[1:size-1+(size%2)])
        c.append("DEPOLARIZE2", path[1:size-1+(size%2)], p)
    c.append("TICK")

    for path in paths:
        c.append("H", path[:-1][1::2])
        c.append("DEPOLARIZE1", path[:-1][2::2], p/10)
    c.append("TICK")

    for path in paths:
        c.append("X_ERROR", path[1:-1], p)
        c.append("MR", path[1:-1])
        c.append("X_ERROR", path[1:-1], p/10)
    c.append("TICK")

    for j, path in enumerate(paths):
        tot_len = sum([len(p[1:-1]) for p in paths[j:]])
        size = len(path)
        for i in range(1 + (size%2), size-1, 2):
            c.append("CZ", [stim.target_rec(-tot_len+i-1), path[0]])
        c.append("DEPOLARIZE1", path[0], p/10)
        for i in range(2 - (size%2), size-1, 2):
            c.append("CX", [stim.target_rec(-tot_len+i-1), path[-1]])
        c.append("DEPOLARIZE1", path[-1], p/10)
    c.append("TICK")

    return c

def lr_CNOT_bell(paths, p):
    c = stim.Circuit()

    for path in paths:
        c.append("CNOT", [path[0], path[2][0], path[2][1], path[1]])
        c.append("DEPOLARIZE2", [path[0], path[2][0], path[2][1], path[1]], p)
    c.append("TICK")

    for path in paths:
        c.append("X_ERROR", [path[2][0], path[2][1]], p)
        c.append("MR", path[2][0])
        c.append("MRX", path[2][1])
        c.append("X_ERROR", [path[2][0], path[2][1]], p/10)
    c.append("TICK")

    for j, path in enumerate(paths[::-1]):
        c.append("CX", [stim.target_rec(-2*j-2), path[1]])
        c.append("CZ", [stim.target_rec(-2*j-1), path[0]])
        c.append("DEPOLARIZE1", [path[0], path[1]], p/10)
    c.append("TICK")
    return c

def lr_CNOT_no_bell(paths):
    # path[0] is control, path[-1] is target. Reverse path to get reverse CNOT
    c = stim.Circuit()

    for path in paths:
        size = len(path)
        c.append("H", path[:-1][2::2])
        c.append("DEPOLARIZE1", path[:-1][2::2], 0.001)
    c.append("TICK")

    for path in paths:
        size = len(path)
        c.append("CNOT", path[:size-(size%2)])
        c.append("DEPOLARIZE2", path[:size-(size%2)], 0.001)
    c.append("TICK")

    for path in paths:
        size = len(path)
        c.append("CNOT", path[1:size-1+(size%2)])
        c.append("DEPOLARIZE2", path[1:size-1+(size%2)], 0.001)
    c.append("TICK")

    for path in paths:
        c.append("H", path[:-1][1::2])
        c.append("DEPOLARIZE1", path[:-1][2::2], 0.001)
    c.append("TICK")

    for path in paths:
        c.append("X_ERROR", path[1:-1], 0.001)
        c.append("MR", path[1:-1])
        c.append("X_ERROR", path[1:-1], 0.001)
    c.append("TICK")

    for j, path in enumerate(paths):
        tot_len = sum([len(p[1:-1]) for p in paths[j:]])
        size = len(path)
        for i in range(1 + (size%2), size-1, 2):
            c.append("CZ", [stim.target_rec(-tot_len+i-1), path[0]])
        c.append("DEPOLARIZE1", path[0], 0.001)
        for i in range(2 - (size%2), size-1, 2):
            c.append("CX", [stim.target_rec(-tot_len+i-1), path[-1]])
        c.append("DEPOLARIZE1", path[-1], 0.001)
    c.append("TICK")

    return c, sum([len(p[1:-1]) for p in paths])


def test_bell(l, p):
    bell_circuit = lr_bell_pair([np.arange(l)], p)
    bell_circuit2 = lr_bell_pair([[i+l for i in np.arange(l)]], p)
    # bell_circuit3 = lr_bell_pair([[i+2*l for i in np.arange(l)]], p)

    control = 100
    target = 101

    c_0 = 0
    c_1 = 0
    c_01 = 0
    c_11 = 0

    num_iters = 1
    for i in range(num_iters):
        s = stim.TableauSimulator()
        s.do(bell_circuit)
        s.do(bell_circuit2)
        # s.do(bell_circuit3)


        # s.do(stim.Circuit(
        # f"""
        # # SQRT_X {2*l} {3*l-1}
        # CX {2*l} 0 {3*l-1} {l-1}
        # # DEPOLARIZE2(0.001) {2*l} 0 {3*l-1} {l-1}
        # H {2*l} {3*l-1}
        # # X 0
        # M {2*l} {3*l-1}
        # """
        # ))
        # res = s.current_measurement_record()
        # if res[-1] != res[-2]: continue


        s.do(stim.Circuit(
        f"""
        CX 0 {l} {l-1} {2*l-1}
        DEPOLARIZE2({p}) 0 {l} {l-1} {2*l-1}
        X_ERROR({p}) {l} {2*l-1}
        MR {l} {2*l-1}
        """
        ))
        res = s.current_measurement_record()

        if res[-1] != res[-2]:
            c_0 += 1

            s.h(control)
            s.do(lr_CNOT_bell([[control, target, (0, l-1)]], p))
            res = s.measure_many(control, target)

            if res[0] != res[1]: c_01 += 1
        else:
            c_1 += 1

            # if (np.random.random() < 0.5): s.x(control)
            s.h(control)
            s.do(lr_CNOT_bell([[control, target, (0, l-1)]], p))
            res = s.measure_many(control, target)
            # res = s.measure_many(0, l-1)

            if res[0] != res[1]: c_11 += 1

    res = (l, c_01/c_0, c_11/c_1, c_1/num_iters)
    return res

arr = [test_bell(i, 0.001) for i in range(52,82,2)]