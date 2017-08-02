################################################################################
# This program measures the evidence for domain interactions in protein 
# interaction data. All of the lookup was done by a previous program and 
# dumped using the pickle module
################################################################################

import sys
#from MySQLdb import *
from functions import *
from time import time
from math import exp
from re import search
from pickle import load
import psyco
psyco.full()


TRUE = 1
FALSE = 0

# PROCESS COMMAND LINE
if len(sys.argv) < 3 or len(sys.argv) > 4:
  print "ARGS: <input prefix> <output prefix>"
  sys.exit()
elif len(sys.argv) == 3:
  list_provided = FALSE
elif len(sys.argv) == 4:
  list_provided = TRUE
  dom_pair_file = sys.argv[3]
input_name = sys.argv[1]
output_name = sys.argv[2]

# OPEN SOME OUTPUT FILES
out = open(output_name + '.out', 'w')
logfile = open(output_name + '.log', 'w')

# DEFINE ALL DATA STRUCTURES
nodes = {}        # Nodes, each being a protein. node_key -> [domain0, domain1, ... , domainN]
edges = {}        # The edges. Each is a protein-protein interaction. edge_key -> (n1, n2) 
P = []            # List list of potentially interacting domain pairs. For lookup mainly.
Q = []            # Nonredundant list of potentially interacting domain pairs. For iterating and for output.
pairs = []        # A list read from a file, if provided. Otherwise Q.
ddi2edge = {}     # Edges dependent on a given domain interaction. (i,j) -> [e0, e1,...,eN]
R = {}            # Number of reported protein pairs interacting through i<->j
S = {}            # Number of possible protein pairs interacting through i<->j
theta = {}        # Pr(i<->j)
C = {}            # edge_key -> [(i0, j0) -> c0, (i1, j1) -> c1, ..., (iN, jN) -> cN]. 1 if domains interact, 0 if not
M = {}            # C[(m,n)][(i,j)] summed over all m<->n dependent on i<->j. The i<->j pairs that interact in interacting proteins.
N = {}            # 1 - C[(m,n)][(i,j)] summed all m<->n dependent on i<->j. The i<->j pairs that do not interact in interacting proteins.
Z = {}            # Z = S - R
E = {}            # Total evidence for the domain interaction i<->j

# SOME PARAMETERS. MODIFY AS NEEDED
alpha = 1.0
beta = 1.0
p_random = .001  # A hack
NUM_ITERS = 100
t0 = 0.0
t1 = 0.0

# FILL UP THE DATA STRUCTURES
f = open(input_name + '_ddi2edge', 'r')
ddi2edge = load(f)
f.close()

f = open(input_name + '_edges', 'r')
edges = load(f)
f.close()

f = open(input_name + '_nodes', 'r')
nodes = load(f)
f.close()

f = open(input_name + '_C', 'r')
C = load(f)
f.close()

f = open(input_name + '_M', 'r')
M = load(f)
f.close()

f = open(input_name + '_N', 'r')
N = load(f)
f.close()

f = open(input_name + '_P', 'r')
P = load(f)
f.close()

f = open(input_name + '_Q', 'r')
Q = load(f)
f.close()

f = open(input_name + '_R', 'r')
R = load(f)
f.close()

f = open(input_name + '_S', 'r')
S = load(f)
f.close()

f = open(input_name + '_Z', 'r')
Z = load(f)

# nodes
# Nodes, each being a protein. node_key -> [domain0, domain1, ... , domainN]
t0 = time()

if list_provided:
  f = open(dom_pair_file, 'r')
  for line in f.readlines():
    pattern = search(r'^(\S+)\s+(\S+)$', line)
    if pattern:
      d1 = pattern.group(1)
      d2 = pattern.group(2)
      pairs.append((d1,d2))
else:
  pairs = Q
  
# theta
logfile.write("Initializing theta\n")
logfile.flush()
initialize_theta(theta, P, R, S, alpha, beta)
t1 = time()
logfile.write("theta:\t%.2f\n" % (t1 - t0))
logfile.flush()

##################################################
# ERROR CHECKING FOR ALL STRUCTURES INDEXED BY i,j
##################################################
check_symmetry(Q, (R, theta, M, N))

######
# EM
######
logfile.write("Running EM\n")
logfile.flush()
L = calc_likelihood(Q, theta, M, N, Z, alpha, beta)
t1 = time()
logfile.write("Initial likelihood: %.2f\n" % (L))
logfile.flush()

for iter in range(NUM_ITERS):
  update_theta(Q, theta, M, N, Z, alpha, beta)
  update_C(C, theta)
  update_MN(C, M, N)
  L = calc_likelihood(Q, theta, M, N, Z, alpha, beta)
  logfile.write("  %d\t%f\n" % (iter, L))
  logfile.flush()
  check_symmetry(Q, (R, S, theta, M, N, Z))

t1 = time()
logfile.write("Initial EM finished:\t%.2f\n" % (t1 - t0))
logfile.flush()

######
# E
######
# Initialize
for i,j in P:
  E[(i,j)] = 0.0
  

# Now compute yet another way, rerunning EM for each excluded domain interaction
for i,j in pairs:
  if (i,j) in P:
    logfile.write("Recalculating excluding %s <-> %s\n" % (i,j))
    logfile.flush()
    # First copy the relevant data structures with i <-> j excluded
    theta_null = {}  # Re-estimated theta, theta_null[(i,j)] set to zero
    initialize_theta(theta_null, P, R, S, alpha, beta)  # Give them the standard initial values
    initialize_C(C, nodes, edges)
    populate_MN(C, M, N)
    # Now rerun EM for with this i,j eliminated
    iter = 0
    L_prev = L
    L_curr = L
    while iter < 10 or exp(L_prev - L_curr) < .99: 
    #for iter in range(10):
      update_theta(Q, theta_null, M, N, Z, alpha, beta)
      # THE KEY STEP, ALL OTHER THINGS WILL FALL INTO PLACE FROM THIS:
      theta_null[(i,j)] = p_random
      theta_null[(j,i)] = p_random
      # END KEY STEP
      update_C(C, theta_null)
      update_MN(C, M, N)
      L_prev = L_curr  # Save the old value to see how much likelihood is changing
      L_curr = calc_likelihood(Q, theta_null, M, N, Z, alpha, beta)
      logfile.write("    %d\t%f\n" % (iter, L_curr))
      logfile.flush()
      check_symmetry(Q, (theta_null, M, N))
      iter += 1
    calc_Eij_null(i, j, ddi2edge, edges, C, theta, theta_null, p_random, E)

########
# OUTPUT
########
#for i,j in Q:
    out.write("%s\t%s\t%f\t%f\t%f\n" % (i, j, float(R[(i,j)]) / S[(i,j)], theta[(i,j)], E[(i,j)]))
    out.flush()
t1 = time()
logfile.write("All parameters evaluated:\t%.2f\n" % (t1 - t0))
logfile.flush()
