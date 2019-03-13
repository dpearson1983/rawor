import numpy as np
import rawor
import time

V_box = 1024.0**3
num_parts = 499362
num_rans = 4993620
num_shells = 32
r_max = 32.0
r_min = 0.0
DD = np.int_([41388, 145384, 175586, 215284, 248902, 302862, 361566, 421164, 482874, 551938, 623613, 695893,
            774264, 856986, 942491, 1036535, 1129877, 1230247, 1336693, 1449855, 1562763, 1685596, 1809689,
            1940755, 2079489, 2225490, 2376902, 2531352, 2688256, 2855352, 3032008, 3205774])
print(DD)

my_predictor = rawor.rawor(num_parts, num_rans, num_shells, V_box, r_max, r_min)
print(my_predictor.get_num_parts())
print(my_predictor.get_num_rans())
print(my_predictor.get_num_shells())
print(my_predictor.get_V_box())
print(my_predictor.get_r_max())
print(my_predictor.get_r_min())
start = time.time()
RRR = my_predictor.get_RRR()
DRR = my_predictor.get_DRR()
DDR = my_predictor.get_DDR(DD)
end = time.time()
test = np.zeros(len(RRR))

print(RRR)
print(DRR)
print(DDR)
print(test)
print(RRR[0])
print(RRR[1])
print(RRR[2])
print(RRR[3])
print(RRR[4])

ZZZ = RRR*test
print(ZZZ)

print("Time to predict randoms: " + str(end - start) + " s")
