
import numpy as np
np.random.seed(2)

def prune(R,A,B):
	n = A.shape[0]
	for i in range(n):
		for j in range(n):
			if R[j,i] == 1:
				# now i,j are potential 
				for k in range(n):
					if A[i,k] == 1:
						#now k is a neighbor of i
						possible = 0
						for l in range(n):
							if B[j,l] == 1:
								# now l is a neighbor
								if R[l,k] == 1:
									possible = 1
						if possible == 0:
							R[i,j] = 0
							continue 

	return R

def prune2(R,A,B):	
	for i in range(n):
		found = 0
		for j in range(n):
			if R[j,i] == 1:
				found += 1
				ind = j
		if found == 1:
			R[ind,:] = 0
			R[ind,i] = 1

	for i in range(n):
		found = 0
		for j in range(n):
			if R[i,j] == 1:
				found += 1
				ind = j
		if found == 1:
			R[:,ind] = 0
			R[i,ind] = 1

	return R

def create_one_mask(v, w):
    sp2 = np.outer(v, np.ones(n)) - np.outer(np.ones(n), w)
    #print sp2
    return(np.abs(sp2).T < 1e-2 ).astype(int)#.astype(float)


n = 10

for i in range(1):
	perm = np.random.permutation(n)
	P = np.zeros([n,n])
	for i in range(n):
		P[i,perm[i]] = 1


	A = np.random.rand(n,n)
	A = np.round((A+A.T)/2)
	B = np.dot(P,np.dot(A,P.T))
	B = A

	degrees_A = np.dot(A,np.ones(n))
	degrees_B = np.dot(B,np.ones(n))

	R = create_one_mask(degrees_A,degrees_B)
	#print A
	#print R
	print ((R))
	Rp = prune(R,A,B)
	#print Rp
	print ((Rp))
	Rpp = prune2(Rp,A,B)
	print ((Rpp))