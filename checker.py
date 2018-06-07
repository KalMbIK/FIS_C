import numpy as np
import scipy
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import gmres
import csv

def read_mtx_toCOO(fileName):
    I = []
    J = []
    V = []
    with file(fileName, 'r') as csvfile:
        reader = csv.reader(csvfile, delimiter=' ')
        reader.next()
        row = reader.next()
        nx = int(row[0])
        ny = int(row[1])
        nz = int(row[2])
        for row in reader:
            I.append(int(row[0])-1)
            J.append(int(row[1])-1)
            if row[2] is not '':
                V.append(float(row[2]))
            else:
                V.append(float(row[3]))
    I = np.array(I)
    J = np.array(J)
    V = np.array(V)
    coo = coo_matrix((V, (I,J)))
    return [nx, ny, nz, I, J, V, coo]

def read_mtx_toCOO_sym(fileName):
    I = []
    J = []
    V = []
    with file(fileName, 'r') as ffile:
        rows = file.readlines(ffile)
        i = 0
        for row in rows:
            if row[0] is not '%':
                numbrs = row.split()
                i += 1
                if i is 1:
                    nx = int(numbrs[0])
                    ny = int(numbrs[1])
                    nz = int(numbrs[2])
                else:
                    I.append(int(numbrs[0])-1)
                    J.append(int(numbrs[1])-1)
                    V.append(float(numbrs[2]))
    I = np.array(I)
    J = np.array(J)
    V = np.array(V)
    coo = coo_matrix((V, (I,J)))
    return [nx, ny, nz, I, J, V, coo]


def save_npVec_toFile(fileName, vector):
    vector_strings = [str(el)+'\n' for el in vector]
    vector_strings.insert(0,str(len(vector_strings))+'\n')
    with file(fileName,'w') as saveFile:
        saveFile.writelines(vector_strings)

myMatrix = 'data/myMatrix.mtx'
orsirr = 'data/orsirr_1.mtx'
s3rmt3m3 = 'data/s3rmt3m3.mtx'
nx, ny, nz, I, J, V, coo = read_mtx_toCOO_sym(orsirr)
# x_star = [i for i in range(nx)]
x_star = np.ones(nx)
A = coo.todense()
A = np.array(A)
A = A+np.transpose(A)-np.diag(A.diagonal())
b = A.dot(x_star)
b = coo.dot(x_star)
print b
save_npVec_toFile('checkCOO.txt', b)
# mK = 400
# x = gmres(coo,b,np.zeros_like(b),1e-8,mK)
# print np.linalg.norm(x[0]-x_star)
# print x
