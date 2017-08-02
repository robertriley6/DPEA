################################################################################
# This program measures the evidence for domain interactions in protein
# interaction networks of one or more species. It does so by:
#   1. Counting domain co-occurrence to estimate a simple probability
#      of domain i interacting with domain j, denoted P(i <-> j).
#   2. Refining the estimate of P(i <-> j) by EM to get a maximum likelihood
#      estimate.
#   3. Calculating the drop in likelihood of the observations when the domain
#      interaction i<->j is excluded from the model. Expresses this as a log
#      odds ratio and sums this over all interacting protein pairs m <-> n
#      that may interact through i<->j.
#
#      Robert Riley, November 2004, UCLA
################################################################################

import sys
from MySQLdb import *
from functions import *
from time import time
from math import exp
from re import search
from pickle import dump

TRUE = 1
FALSE = 0

# PROCESS COMMAND LINE
if len(sys.argv) != 2:
  print "ARGS: <run name>"
  sys.exit()
run_name = sys.argv[1]

# OPEN SOME OUTPUT FILES
logfile = open(run_name + '_lookup.log', 'w')

# INPUT DATA IN MYSQL TABLES
# CHANGE AS NEEDED
node_table = 'combined_nodes'
edge_table = 'combined_edges'
domain_table = 'combined_domains'

# OPEN A MYSQL CONNECTION
db = connect(host = "mysql1", user = "riley", passwd = "changeme", db = "riley")
c = db.cursor()

# DEFINE ALL DATA STRUCTURES
nodes = {}        # Nodes, each being a protein. node_key -> [domain0, domain1, ... , domainN]
edges = {}        # The edges. Each is a protein-protein interaction. edge_key -> (n1, n2) 
taxa = []         # The distinct taxa in our data. [t0, t1, ... , tN]
taxa2edge = {}    # The edges for a given taxa. t0 = [e00, e01, ... , e0N]
taxa2node = {}    # The nodes in a given taxa. t0 -> [n00, n01, ... , n0N]
D = []            # The list of all domain types in our data
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

# SOME PARAMETERS
t0 = 0.0
t1 = 0.0

# FILL UP THE DATA STRUCTURES
# nodes
# Nodes, each being a protein. node_key -> [domain0, domain1, ... , domainN]
t0 = time()
c.execute("""SELECT uid FROM %s""" % node_table)
for r in c.fetchall():
  n = r[0]
  nodes[n] = []

# node -> domain mappings
# Nodes, each being a protein. node_key -> [domain0, domain1, ... , domainN]
c.execute("""SELECT * FROM %s""" % domain_table)
for r in c.fetchall():
  n = r[0]
  d = r[1]
  nodes[n].append(d)
t1 = time()
logfile.write("nodes:\t%.2f\n" % (t1 - t0))
logfile.flush()

# edges
# The edges for a given taxa. t0 = [e00, e01, ... , e0N]
c.execute("""SELECT uid, n1, n2 FROM %s""" % edge_table)
for r in c.fetchall():
  uid = r[0]
  n1 = r[1]
  n2 = r[2]
  edges[uid] = (n1, n2)
t1 = time()
logfile.write("edges:\t%.2f\n" % (t1 - t0))
logfile.flush()

# taxa
# The distinct taxa in our data. [t0, t1, ... , tN]
c.execute("""SELECT DISTINCT organism FROM %s""" % node_table)
for r in c.fetchall():
  t = r[0]
  taxa.append(t)
  taxa2edge[t] = []
  taxa2node[t] = []

# taxa2edge and taxa2node
# The edges for a given taxa. t0 = [e00, e01, ... , e0N]
# The nodes in a given taxa. t0 -> [n00, n01, ... , n0N]
for t in taxa:
  c.execute("""SELECT e.uid FROM %s AS e, %s AS m, %s AS n 
                 WHERE e.n1 = m.uid 
                   AND e.n2 = n.uid 
                   AND m.organism = '%s' 
                   AND n.organism = '%s'""" % (edge_table, node_table, node_table, t, t))
  for r in c.fetchall():
    e = r[0]
    taxa2edge[t].append(e)
  c.execute("""SELECT uid FROM %s WHERE organism = '%s'""" % (node_table, t))
  for r in c.fetchall():
    n = r[0]
    taxa2node[t].append(n)
t1 = time()
logfile.write("taxa, taxa2edge, taxa2node:\t%.2f\n" % (t1 - t0))
logfile.flush()

# D
# The list of all domain types in our data
c.execute("""SELECT DISTINCT pf_name FROM %s""" % domain_table)
for r in c.fetchall():
  d = r[0]
  D.append(d)
t1 = time()
logfile.write("D:\t%.2f\n" % (t1 - t0))
logfile.flush()

# P
# List list of potentially interacting domain pairs. For lookup mainly.
c.execute("""SELECT DISTINCT i.pf_name, j.pf_name 
               FROM %s AS e, %s AS i, %s AS j, %s AS m, %s AS n
                 WHERE ((e.n1 = i.uid AND e.n2 = j.uid) OR (e.n1 = j.uid AND e.n2 = i.uid))
		   AND (e.n1 = m.uid AND e.n2 = n.uid)
		   AND (m.organism = n.organism)""" 
		   % (edge_table, domain_table, domain_table, node_table, node_table))
for r in c.fetchall():
  d1 = r[0]
  d2 = r[1]
  P.append((d1,d2))
t1 = time()
logfile.write("P:\t%.2f\n" % (t1 - t0))
logfile.flush()

# Q
# Nonredundant list of potentially interacting domain pairs. For iterating and for output.
# May be horribly slow if we sort out redundancy in this program.
# However, we can design a query that asks that d1 >= d2. This should result in all domain
# pairs where d1 = d2, or where d1 > d2, which means we get either (di, dj) or (dj, di) but
# not both. Ha.
c.execute("""SELECT DISTINCT i.pf_name, j.pf_name 
               FROM %s AS e, %s AS i, %s AS j, %s AS m, %s AS n
                 WHERE ((e.n1 = i.uid AND e.n2 = j.uid) OR (e.n1 = j.uid AND e.n2 = i.uid)) 
		   AND (e.n1 = m.uid AND e.n2 = n.uid)
		   AND (m.organism = n.organism)
                   AND (i.pf_name >= j.pf_name)""" 
		   % (edge_table, domain_table, domain_table, node_table, node_table))
for r in c.fetchall():
  d1 = r[0]
  d2 = r[1]
  Q.append((d1,d2))
logfile.write("%d parameters to estimate\n" % (len(Q)))
logfile.flush()
t1 = time()
logfile.write("Q:\t%.2f\n" % (t1 - t0))
logfile.flush()

#if list_provided:
#  f = open(dom_pair_file, 'r')
#  for line in f.readlines():
#    pattern = search(r'^(\S+)\s+(\S+)$', line)
#    if pattern:
#      d1 = pattern.group(1)
#      d2 = pattern.group(2)
#      pairs.append((d1,d2))
#else:
#  pairs = Q
  
#ddi2edge
# Edges dependent on a given domain interaction. (i,j) -> [e0, e1,...,eN]
for i,j in Q:
  ddi2edge[(i,j)] = []
  ddi2edge[(j,i)] = []
  c.execute("""SELECT DISTINCT e.uid FROM %s AS e, %s AS i, %s AS j 
                 WHERE (e.n1 = i.uid AND e.n2 = j.uid) 
                   AND ((i.pf_name = '%s' AND j.pf_name = '%s') 
                     OR (j.pf_name = '%s' AND i.pf_name = '%s'))""" % 
               (edge_table, domain_table, domain_table, i, j, i, j))
  for r in c.fetchall():
    ddi2edge[(i,j)].append(r[0])
    if i != j:
      ddi2edge[(j,i)].append(r[0])
t1 = time()
logfile.write("ddi2edge:\t%.2f\n" % (t1 - t0))
logfile.flush()

#R            
# Number of reported protein pairs interacting through i<->j
for i,j in P:
  c.execute("""SELECT COUNT(DISTINCT e.uid) FROM %s AS e, %s AS i, %s AS j 
                 WHERE (e.n1 = i.uid AND e.n2 = j.uid) 
                   AND ((i.pf_name = '%s' AND j.pf_name = '%s') 
                     OR (j.pf_name = '%s' AND i.pf_name = '%s'))""" % 
               (edge_table, domain_table, domain_table, i, j, i, j))
  for r in c.fetchall():
    R[(i,j)] = r[0]
t1 = time()
logfile.write("R:\t%.2f\n" % (t1 - t0))
logfile.flush()

#S
# Number of possible protein pairs interacting through i<->j
for i,j in P:
  S[(i,j)] = 0
for t in taxa:
  for i,j in P:
    c.execute("""SELECT COUNT(DISTINCT m.uid) 
                   FROM %s AS m, %s AS d
                     WHERE m.uid = d.uid
                       AND d.pf_name = '%s'
                       AND m.organism = '%s'""" % (node_table, domain_table, i, t))
    r = c.fetchone()
    x = r[0]
    c.execute("""SELECT COUNT(DISTINCT m.uid) 
                   FROM %s AS m, %s AS d
                     WHERE m.uid = d.uid
                       AND d.pf_name = '%s'
                       AND m.organism = '%s'""" % (node_table, domain_table, j, t))
    r = c.fetchone()
    y = r[0]
    S[(i,j)] += (x * y)
t1 = time()
logfile.write("S:\t%.2f\n" % (t1 - t0))
logfile.flush()

# Z
for i,j in P:
  Z[(i,j)] = S[(i,j)] - R[(i,j)]
t1 = time()
logfile.write("Z:\t%.2f\n" % (t1 - t0))
logfile.flush()

# theta
logfile.write("Initializing theta\n")
logfile.flush()
initialize_theta(theta, P, R, S, 0, 0)
t1 = time()
logfile.write("theta:\t%.2f\n" % (t1 - t0))
logfile.flush()

# C
# edge_key -> [(i0, j0) -> c0, (i1, j1) -> c1, ..., (iN, jN) -> cN]. 1 if domains interact, 0 if not
logfile.write("Initializing C\n")
logfile.flush()
initialize_C(C, nodes, edges)
t1 = time()
logfile.write("C:\t%.2f\n" % (t1 - t0))
logfile.flush()

# initialize M, N
logfile.write("Initializing M, N\n")
logfile.flush()
initialize_MN(P, M, N)
t1 = time()
logfile.write("M, N initialized:\t%.2f\n" % (t1 - t0))
logfile.flush()
#for i,j in P:
#  M[(i,j)] = 0.0
#  N[(i,j)] = 0.0

# populate M, N
logfile.write("Populating M, N\n")
logfile.flush()
populate_MN(C, M, N)
t1 = time()
logfile.write("M, N populated:\t%.2f\n" % (t1 - t0))
logfile.flush()

##################################################
# ERROR CHECKING FOR ALL STRUCTURES INDEXED BY i,j
##################################################
check_symmetry(Q, (R, theta, M, N))

#############################################################################
# The data structures we will need to run EM from another program
# nodes, edges, ddi2edge, P, Q, R, S, C, M, N, Z 
#############################################################################
f = open(run_name + '_ddi2edge', 'w')
dump(ddi2edge, f)
f.close()

f = open(run_name + '_edges', 'w')
dump(edges, f)
f.close()

f = open(run_name + '_nodes', 'w')
dump(nodes,f)
f.close()

f = open(run_name + '_C', 'w')
dump(C, f)
f.close()

f = open(run_name + '_M', 'w')
dump(M, f)
f.close()

f = open(run_name + '_N', 'w')
dump(N, f)
f.close()

f = open(run_name + '_P', 'w')
dump(P, f)
f.close()

f = open(run_name + '_Q', 'w')
dump(Q, f)
f.close()

f = open(run_name + '_R', 'w')
dump(R, f)
f.close()

f = open(run_name + '_S', 'w')
dump(S, f)
f.close()

f = open(run_name + '_Z', 'w')
dump(Z, f)
f.close()
