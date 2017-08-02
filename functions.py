import sys
from math import log

def initialize_theta(theta, P, R, S, alpha, beta):
  for i,j in P:
    #theta[(i,j)] = float(R[(i,j)] + alpha) / (S[(i,j)] + alpha + beta)
    theta[(i,j)] = float(R[(i,j)]) / S[(i,j)]
    if theta[(i,j)] < 0.0 or theta[(i,j)] > 1.0:
      print "ERROR: theta[(%s,%s)] = %.30f in initialization" % (i, j, theta[(i,j)])
      sys.exit()

def calc_likelihood(Q, theta, M, N, Z, alpha, beta):
  l = 0.0
  bottom = 10e-30
  for i,j in Q:
    a = max(pow(theta[(i,j)], M[(i,j)] + alpha), bottom)
    # Below I had M instead of N. Serious bug.
    b = max(pow(1-theta[(i,j)], N[(i,j)] + Z[(i,j)] + beta), bottom)
    p = log(a)
    q = log(b)
    l += (p + q)
  return l

def initialize_C(C, nodes, edges):
  for e in edges.keys():
    C[e] = {}
    m = nodes[edges[e][0]]
    n = nodes[edges[e][1]]
    for i in m:
      for j in n:
        C[e][(i,j)] = 1.0

def update_C(C, theta):
  for e in C.keys():
    # First the denominator, P(m <-> n)
    prod = 1.0
    for k,l in C[e]:
      prod *= (1.0 - theta[(k,l)])
    #p = round(1.0 - prod, 4)
    p = 1.0 - prod
    for i,j in C[e]:
      #q =  round(theta[(i,j)], 4) / p
      q =  round(theta[(i,j)] / p, 4)
      if q < 0.0 or q > 1.0:
        print "ERROR: C[%d][(%s,%s)] =\t%.30f" % (e, i, j, q) 
	print "  theta[(%s,%s)] =\t%.30f" % (i, j, theta[(i,j)])
	print "  p =\t\t\t%.30f" % p
	sys.exit()
      C[e][(i,j)] = q

def initialize_MN(P, M, N):
  for i,j in P:
    M[(i,j)] = 0.0
    N[(i,j)] = 0.0

def populate_MN(C, M, N):
  for e in C.keys():
    for i,j in C[e]:
      M[(i,j)] += C[e][(i,j)]
      if i != j:
        M[(j,i)] += C[e][(i,j)]
      N[(i,j)] += 1 - C[e][(i,j)]
      if i != j:
        N[(j,i)] += 1 - C[e][(i,j)]

def update_MN(C, M, N):
  # Zero them out
  for e in C.keys():
    for i,j in C[e]:
      M[(i,j)] = 0.0
      M[(j,i)] = 0.0
      N[(i,j)] = 0.0
      N[(j,i)] = 0.0
  # Now count them
  for e in C.keys():
    for i,j in C[e]:
      M[(i,j)] += C[e][(i,j)]
      if i != j:
        M[(j,i)] += C[e][(i,j)]
      N[(i,j)] += (1 - C[e][(i,j)])
      if i != j:
        N[(j,i)] += (1 - C[e][(i,j)])
  # Now do an integrity check. Make sure no counts are below zero. There was
  # a bug with raising negative numbers to fractional powers.
  for e in C.keys():
    for i,j in C[e]:
      if M[(i,j)] < 0.0:
        print "ERROR: M[(%s,%s)] = %f" % (i, j, M[(i,j)])
	sys.exit()
      if M[(j,i)] < 0.0:
        print "ERROR: M[(%s,%s)] = %f" % (j, i, M[(j,i)])
	sys.exit()
      if N[(i,j)] < 0.0:
        print "ERROR: N[(%s,%s)] = %f" % (i, j, N[(i,j)])
	sys.exit()
      if N[(j,i)] < 0.0:
        print "ERROR: N[(%s,%s)] = %f" % (j, i, M[(j,i)])
	sys.exit()

def update_theta(Q, theta, M, N, Z, alpha, beta):
  for i,j in Q:
    theta[(i,j)] = (M[(i,j)] + alpha) / (M[(i,j)] + N[(i,j)] + Z[(i,j)] + alpha + beta)
    theta[(j,i)] = (M[(i,j)] + alpha) / (M[(i,j)] + N[(i,j)] + Z[(i,j)] + alpha + beta)
    #theta[(i,j)] = (M[(i,j)] + alpha) / (M[(i,j)] + N[(i,j)]  + beta)
    #theta[(j,i)] = (M[(i,j)] + alpha) / (M[(i,j)] + N[(i,j)]  + beta)
    if theta[(i,j)] < 0.0 or theta[(i,j)] > 1.0:
      print "ERROR: theta[(%s,%s)] = %.30f in update" % (i, j, theta[(i,j)])
      sys.exit()
    

def check_symmetry(Q, list):
  for l in list:
    for i,j in Q:
      if l[(i,j)] != l[(j,i)]:
        print "ERROR: assymetry"
        sys.exit()

def calc_Eij(i, j, ddi2edge, edges, C, theta, p_random, E):
  for e in ddi2edge[(i,j)]:
    # Numerator
    nprod = 1.0
    for k,l in C[e]:
      nprod *= (1 - theta[(k,l)])
    p = max(1 - nprod, p_random)
    # Denominator
    dprod = 1.0
    for k,l in C[e]:
      if (k,l) != (i,j) and (k,l) != (j,i):
        dprod *= (1 - theta[(k,l)])
    q = max(1 - dprod, p_random)
    E[(i,j)] += log(p / q)
    if i != j:
      E[(j,i)] += log(p / q)


# FUNCTIONS FOR DEALING WITH AN EXCLUDED EDGE

def calc_Eij_null(i, j, ddi2edge, edges, C, theta, theta_null, p_random, E):
  for e in ddi2edge[(i,j)]:
    # Numerator
    nprod = 1.0
    for k,l in C[e]:
      nprod *= (1 - theta[(k,l)])
    p = max(1 - nprod, p_random)
    # Denominator
    dprod = 1.0
    for k,l in C[e]:
      #if (k,l) != (i,j) and (k,l) != (j,i):
      dprod *= (1 - theta_null[(k,l)])
    q = max(1 - dprod, p_random)
    E[(i,j)] += log(p / q)
    if i != j:
      E[(j,i)] += log(p / q)
