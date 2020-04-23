import os
import struct
import numpy as np
from scipy.linalg import svd
import matplotlib.pyplot as plt

fname = "dreggn"
# node size is 536 bytes or half that in uint16 (as the payload data is)
node_size = 536//2
header_size = 96//2
motion_samples=64

nodes = np.fromfile(fname, dtype='<i8')[1]//2
data = np.fromfile(fname, dtype='<f2')

payload_size = 512*motion_samples*2

payload = data[header_size:nodes]
num_blocks=payload.size/payload_size


# analysis of time dependent data:
# time = payload.reshape(num_blocks, 2, 512, -1)[...,1,...,...]
# only look at one block:
time = payload.reshape(num_blocks, 2, 512, -1)[20,1,...,...]
for i in range(0,time.shape[0]):
  plt.plot(time[i,...])
plt.show()


# block compression via svd:
# # reshape into 1D arrays per block (one for temperature, one for density)
# blocks = payload.reshape(num_blocks, 2, -1)
# density = blocks[...,1,...]
# # svd of all blocks:
# U, s, Vt = svd(density, full_matrices=False)
# # s are the sorted singular values, s[0] is the largest one
# # Vt[0,...] is the first eigenvector
# # for i in range(0,density.shape[0]):
# #   plt.plot(density[i,...])
# # plt.show()
# # 
# # for i in range(0,Vt.shape[0]):
# #   plt.plot(density[i,...])
# # plt.show()
# 
# # now project the data to the first couple of eigenvectors,
# # and write the file back
# coeffs=200
# compressed_density=np.array(density, copy=True)
# compressed_density[...]=0
# 
# for i in range(0,density.shape[0]):
#   for c in range(0,coeffs):
#     compressed_density[i,...] += np.dot(density[i,...], Vt[c,...]) * Vt[c,...]
# 
# # open file for random access without truncation
# f = open(fname, "rb+")
# f.seek(header_size*2, 0)
# for i in range(0,density.shape[0]):
#   pos = f.seek(payload_size, 1) # jump over temperature in bytes
#   written = f.write(compressed_density[i,...].data)
# f.close()
  

