
##### routines that a user might actually want to call ###################

#### produce data from a triangulation ####
# triangle_splines(triangle_list, r) # returns (R, s, J) 
# cone_over_triangulation(triangle_list) # return list of triangles

#### produce data from a polyhedral tiling ####
# poly_splines(poly_list, r) # returns (R, s, J) 

#### produce random data = (R, s, J) ####
# random_RsJ(n = None, s = None)
# principal_random_RsJ(n = None, s = None)
# random_triangle_splines(n,r)    #returns (R, s, J, triangulation, r)

#### produce splines from data or spline module generators (or free basis) ####

# spline_module_generators(R, s, J) # returns list of spline module generators
# is_projective(smg) # returns true or false
# compute_free_basis(smg) # returns a free basis  of splines module generators (if the module is free)

#### Hilbert series routines from data ####

# spline_Hilbert_series(R, s, J) # returns the Hilbert series as a rational function
# spline_Hilbert_polynomial(R,s,J): # returns the Hilbert polynomial of the series above


#### ring structure routines from spline module generators (or free basis) ####

# pretty_SP_relation_ideal(smg) #returns an ideal of R[spline variables as x's]

#### true intersection ####

# sheet_intersection(smg, j, k) #returns ideal of R of true intersection

#### returns a triangulation object from an iterable of triangles ####

# type_triangulation(triangles)


###############################################################





####################### the ring OA and its ideals ##################
### the main purpose of OA is to compute the spline module generators ###

def OAey(R, s):
   n = len(R.gens())  
   OA = PolynomialRing(QQ, 'a', s+n) 
   return (OA, OA.gens()[:s], OA.gens()[s:])

def OA_to_tuple(g, R, s):
   (OA, e, y) = OAey(R, s)
   zero_list = [0 for i in range(s)]
   OA_to_R = OA.hom(zero_list+list(R.gens()), R)
   gt = tuple([OA_to_R(diff(g,e[i])) for i in range(s)])
   return gt
   
def tuple_to_OA(g):
   R = g[0].parent()
   s = len(g)
   (OA, e, y) = OAey(R, s)
   R_to_OA = R.hom(y, OA)
   ge = sum([R_to_OA(g[i]) *e[i] for i in range(s)])
   return ge   

def M_ideal(R, s):
   (OA, e, y) = OAey(R, s)
   M = ideal(e)^2
   return M

def G_ideal(R, s, J):
   (OA, e, y) = OAey(R, s)
   zero_list = [0 for i in range(s)]
   R_to_OA = R.hom(y, OA)
   #this corresponds to the projection from AsY to Y
   M = M_ideal(R, s)
   def jkfilter(j,k):
       return lambda x:((x !=e[j]) and (x !=e[k]))
   def GM(j,k):
      Jjk_in_OA = ideal([ R_to_OA(h) for h in J(j,k).gens()])
      if s<3:
         return (M + Jjk_in_OA*ideal([e[j], e[k]]) + ideal(e[j]+e[k]))
      return (M + Jjk_in_OA*ideal([e[j], e[k]]) + ideal(e[j]+e[k]) + ideal(filter(jkfilter(j,k), e)))
   def intersect(I,J):
        return I.intersection(J)
   if s < 3:
      return GM(0,1)
   G = reduce(intersect,
                 [GM(j,k) for k in range(1,s) for j in range(k)])
   return G

######################################## checking a tuple against an ideal-difference condition ########

def is_J_spline(J, g):
    s = len(g)
    for p in [(i,j) for j in range(s) for i in range(j)]:
        if not (g[p[0]] - g[p[1]] in J(p[0],p[1])):
            return False
    return True

############### generators of the spline R-module and ring.  ######################
##### also a possibly expensive routine to try to find small generating sets #############

def spline_module_generators(R, s, J):
    (OA, e, y) = OAey(R, s)
    zero_tuple = tuple([0 for i in range(s)])
    OA_to_R = OA.hom(list(zero_tuple)+list(R.gens()), R)
    G = G_ideal(R, s, J)
    Grgb = G.groebner_basis()
    long_list = [OA_to_tuple(g, R, s)  for g in Grgb]
    return [g for g in long_list if not (g == zero_tuple)]

def short_spline_module_generators(R, s, J):
   (OA, e, y) = OAey(R, s)
   zero_tuple = tuple([0 for i in range(s)])
   OA_to_R = OA.hom(list(zero_tuple)+list(R.gens()), R)
   G = G_ideal(R, s, J)
   GF = G.groebner_fan() # groebner_fan routines seem to crap-out sometimes
   rGbs = GF.reduced_groebner_bases()
   long_lists = [[OA_to_tuple(g, R, s)  for g in rgb] for rgb in rGbs]
   smg_lists = [[g for g in long_list if not (g == zero_tuple)] for long_list in long_lists]
   shortlist = min(smg_lists, key = len)
   return shortlist

def spline_ring_generators(smg):
   s = len(smg[0])
   one = tuple([1 for i in range(s)])
   srg = [generator for generator in smg if not (generator == one)]
   return srg


########################### the ring OB and its ideals #############################
#### this ring is an intermediary ring for computing the relations in the spline ring SP ####
##################### srg is the list of spline ring generators ####################

def OBexy(srg):
   R = srg[0][0].parent()
   s = len(srg[0])
   n = len(R.gens())  
   num_gens = len(srg)
   OB = PolynomialRing(QQ, s+num_gens+n, 'b')
   return (OB, OB.gens()[:s], OB.gens()[s:s+num_gens], OB.gens()[s+num_gens:])


def E_ideal(srg):
   (OB, e, x, y) = OBexy(srg)
   s = len(srg[0])
   E = ideal([1 - sum(e)] +
             [ei - ei^2 for ei in e] +
             [e[j]*e[k] for k in range(1,s) for j in range(k)])
   return E

def R_ideal(srg):
    R = srg[0][0].parent()
    (OB, e, x, y) = OBexy(srg)
    s = len(srg[0])
    R_to_OB = R.hom(y,OB)
    bigR = (ideal([x[l] - sum([R_to_OB(srg[l][j]) *e[j] for j in range(s)]) for l in range(len(srg))])+E_ideal(srg))
    R = bigR.elimination_ideal(e)
    return R




##################### the ring of splines SP #########################
######## it is computed from a list of spline ring generators ###########

#from sage.rings.polynomial.flatten import FlatteningMorphism

def SPxy(srg):
   R = srg[0][0].parent()
   py = list(R.gens())
   num_gens = len(srg)
   px = list(PolynomialRing(R, 'x', num_gens).gens()) #
   xyvars = px+py 
#   fl = FlatteningMorphism(temp_SP) #
   SP = PolynomialRing(QQ, xyvars) # 
 #  n = len(R.gens())
 #  SP = PolynomialRing(QQ, num_gens+n, 'x')
   x = SP.gens()[:num_gens]
   y = SP.gens()[num_gens:]
   return (SP, x, y)

def R_to_SP(srg):
   R = srg[0][0].parent()
   (SP, x, y) = SPxy(srg)
   out = R.hom(y, SP)
   return out

def tuple_to_SP(smg, in_tup):
   R = smg[0][0].parent()
   s = len(smg[0])
   srg = spline_ring_generators(smg)
   (SP, x, y) = SPxy(srg)
   one_tup = tuple([R(1) for i in range(s)])
   oi = None
   for i in range(len(smg)):
      if smg[i]==one_tup:
         oi = i
   if oi == None:
      oi = 0
      smg = [one_tup] + smg
   lx = list(x)
   lox = lx[:oi]+[SP(1)]+lx[oi:]
   ox = tuple(lox)
   cof = expand_spline(in_tup, smg)
   oot = sum([R_to_SP(srg)(cof[i])*ox[i] for i in range(len(ox))])
   return oot

                   
def SP_relation_ideal(srg):
   (OB, OBe, OBx, OBy) = OBexy(srg)
   (SP, x, y) = SPxy(srg)
   OB_to_SP = OB.hom([0 for ei in OBe] + list(x) + list(y), SP)
   R_gens = R_ideal(srg).gens()
   SPR = ideal(OB_to_SP(rel) for rel in R_gens)
   return SPR

def sigma_SP_to_R(SP, srg, i):
    R = srg[0][0].parent()
    y = R.gens()
    i_ev_lst = [g[i] for g in srg] + list(y)
    sigma_i = SP.hom(i_ev_lst, R)
    return sigma_i

def SP_to_tuple(SP, srg, f):
    sigma_lst = [sigma_SP_to_R(SP, srg, i) for i in range(len(srg[0]))]
    sigma_f = tuple([sig(f) for sig in sigma_lst])
    return sigma_f


 
######## Below are the "pretty" versions of the objects above #######

from sage.rings.polynomial.flatten import FlatteningMorphism

def pretty_SPxy(smg):
    srg = spline_ring_generators(smg)
    R = srg[0][0].parent()
    temp_pSP = PolynomialRing(R, 'x', len(srg))
    fl = FlatteningMorphism(temp_pSP)
    pSP = fl.codomain()
    x = temp_pSP.gens()
    y = R.gens()
    return (pSP, x, y)

def pretty_SP_relation_ideal(smg):
    srg = spline_ring_generators(smg)
    SP = SPxy(srg)[0]
    (pSP , x, y) = pretty_SPxy(srg)
    SP_to_pSP = SP.hom(list(x)+list(y), pSP)
    SPR = SP_relation_ideal(srg)
    pSPR = ideal([SP_to_pSP(rel) for rel in SPR.gens()])
    return pSPR

################# intersection of sheets #######################

def sheet_intersection(smg, j, k):
    srg = spline_ring_generators(smg)
    s = len(srg[0])
    num_gens = len(srg)
    (SP, x, y) = SPxy(srg)
    SPR = SP_relation_ideal(srg)
    R = srg[0][0].parent()
    R_to_SP = R.hom(y, SP)
    Ij = ideal([x[l] - R_to_SP(srg[l][j]) for l in range(num_gens)]) + SPR   
    Ik = ideal([x[l] - R_to_SP(srg[l][k]) for l in range(num_gens)]) + SPR   
    Ijk = Ij + Ik
    Ijk_elim = Ijk.elimination_ideal(list(x))
    SP_to_R = SP.hom([R(0) for xl in x]+list(R.gens()), R)
    out_ideal = ideal([SP_to_R(f) for f in Ijk_elim.gens()])
    return out_ideal

def contact_ideals(smg):
   s = len(smg[0])
   ideals = [[sheet_intersection(smg, j, k) for j in range(s)] for k in range(s)]
   def out_J(u,v):
      return ideals[v][u]
   return out_J

   
############# A Hilbert series routine for degree filtration #########################

def spline_Hilbert_series(R, s, J):
   G = G_ideal(R, s, J)
   Ggb = G.groebner_basis()
   LTgens = [f.lt() for f in Ggb]
   M = M_ideal(R, s)
   LT = ideal(LTgens)+M
   LThs = LT.hilbert_series()
   t = LThs.parent().gens()[0]
   Mhs = M.hilbert_series()
   shs = (Mhs - LThs)/t
   return shs

def spline_Hilbert_polynomial(R,s,J):
   G = G_ideal(R, s, J)
   Ggb = G.groebner_basis()
   LTgens = [f.lt() for f in Ggb]
   M = M_ideal(R, s)
   LT = ideal(LTgens)+M
   LTp = LT.hilbert_polynomial()
   t = LTp.parent().gens()[0]
   Mp = M.hilbert_polynomial()
   P_diff = Mp-LTp
   toss = PolynomialRing(QQ, 'd')
   d = toss.gen()
   P = P_diff(t = d+1)
   return P

########################### triangle input ###############################

def hyperplane_ideal(t1, t2):
   n = len(t1[0])
   R = PolynomialRing(QQ,'y',n)
   y = R.gens()
   if t1 == t2:
      return 0*R
   t1_cap_t2 = [p for p in t1 if p in t2]
   if t1_cap_t2 == []:
      return 1*R
   int_vects = [vector(p) for p in t1_cap_t2]
   based_vects = [v-int_vects[0] for v in int_vects[1:]]
   M = matrix( based_vects + [y])
   minor_size = min(M.nrows(), M.ncols())
   linear_terms = M.minors(minor_size)
   ideal_generators = [l - l(list(t1_cap_t2[0])) for l in linear_terms]
   h_ideal = ideal(ideal_generators)
   return h_ideal

def facet_hyperplane_ideal(t1, t2):
   n = len(t1[0])
   R = PolynomialRing(QQ,'y',n)
   y = R.gens()
   if t1 == t2:
      return 0*R
   t1_cap_t2 = [p for p in t1 if p in t2]
   if not len(t1_cap_t2) == n:
      return 1*R
   int_vects = [vector(p) for p in t1_cap_t2]
   based_vects = [v-int_vects[0] for v in int_vects[1:]]
   M = matrix( based_vects + [y])
   minor_size = min(M.nrows(), M.ncols())
   linear_terms = M.minors(minor_size)
   ideal_generators = [l - l(t1_cap_t2[0]) for l in linear_terms]
   h_ideal = ideal(ideal_generators)
   return h_ideal

def triangle_splines(triangle_list = [], r=0):
    r"""
    Returns a tuple of 
    - the coordinate ring of the affine space containing the 'triangles' (simplices), 
    - the number s of triangles(simplices), 
    - and a function on [0,..,s-1]^2 which outputs ideals of the triangle intersections.
   
    INPUT:
   
    - ``triangle_list`` -- presumbably a list of simplices. Each is itself a list of its vertices. 

    - ``r`` -- a smoothness parameter.
    
    OUTPUT: (R,s,Jr) = (a ring, a number of simplices, an (s x s)-'table' of ideals given as a function of two integers) 
    
    EXAMPLES:
 
    This example is Dalbec and Schenck counterexample to Rose's conjecture ::

        sage: (v0, v1,v2,v3, v4, v5) = ( (3,0,0), (0,3,0), (-2,-2,1), (1,0,-3), (0,0,3), (0,0,0) ) #triangle along which tetrahedra are glued is [v0, v1, v2]
        sage: (t0,t1,t2,t3,t4,t5) = ((v5,v0,v1,v4),(v5,v0,v1,v3),(v5,v2,v4,v0),(v5,v2,v4,v1),(v5,v2,v3,v0),(v5,v2,v3,v1))
        sage: (dsR, dss, dsJ) = triangle_splines((t0, t1, t2, t3, t4, t5) , 0)

    TESTS::

        sage: (v0, v1,v2,v3, v4, v5) = ( (3,0,0), (0,3,0), (-2,-2,1), (1,0,-3), (0,0,3), (0,0,0) ) #triangle along which tetrahedra are glued is [v0, v1, v2]
        sage: (t0,t1,t2,t3,t4,t5) = ((v5,v0,v1,v4),(v5,v0,v1,v3),(v5,v2,v4,v0),(v5,v2,v4,v1),(v5,v2,v3,v0),(v5,v2,v3,v1))
        sage: (dsR, dss, dsJ) = triangle_splines((t0, t1, t2, t3, t4, t5) , 0)
        sage: dsR
        Multivariate Polynomial Ring in y0, y1, y2 over Rational Field
        sage: dss
        6
        sage: dsJ(2,3)
        Ideal (-6*y0 + 6*y1) of Multivariate Polynomial Ring in y0, y1, y2 over Rational Field
    """
    t_list = list(triangle_list)
    n = len(t_list[0][0])
    s = len(t_list)
    R = PolynomialRing(QQ,'y',n)
    def Jr(j,k):
       Jjk =  hyperplane_ideal(t_list[j], t_list[k])
       return Jjk^(r+1)
    return (R,s,Jr)

def triangle_splines_facets(triangle_list, r):
    t_list = list(triangle_list)
    n = len(t_list[0][0])
    s = len(t_list)
    R = PolynomialRing(QQ,'y',n)
    def Jr(j,k):
        Jjk =  facet_hyperplane_ideal(t_list[j], t_list[k])
        return Jjk^(r+1)
    return (R,s,Jr)

def cone_over_triangle(triangle):
    triangle_as_lists = [list(v) for v in triangle]
    for v in triangle_as_lists:
        v.append(1)
    zero_list = [0 for c in triangle_as_lists[0]]
    out_as_lists = triangle_as_lists +[zero_list]
    out= [tuple(v) for v in out_as_lists]
    return tuple(out)

def cone_over_triangulation(triangle_list):
    output_list = [ cone_over_triangle(t) for t in triangle_list ]
    return tuple(output_list)

def random_triangle_list(n):
    num_points = 1+n+abs(ZZ.random_element())
    points = [random_vector(QQ,n)  for i in range(num_points)]
    PointConfiguration.set_engine('internal') 
    pc = PointConfiguration(points)
    triang = pc.triangulate()
    out= (tuple([tuple([tuple(pc[l]) for l in t]) for t in triang]), triang)
    return out

def random_triangle_splines(n,r):
    (rtl, triang) = random_triangle_list(n)
    spline_data = triangle_splines(rtl, r)
    out_list = list(spline_data) + [triang] + [r]
    return tuple(out_list)

########################### Polyhedral input ##############################

def poly_splines(poly_list, r):
    t_list = list(poly_list)
    n = len(t_list[0][0])
    s = len(t_list)
    R = PolynomialRing(QQ,'y',n)
    def J(j,k):
        outJ =  hyperplane_ideal(t_list[j], t_list[k])
        return outJ^(r+1)
    return (R,s,J)

################### two random data generators ############################

def random_RsJ(n = None, s = None):
   if n == None:
      n = 1+abs(ZZ.random_element())
   R = PolynomialRing(QQ, 'y', n)
   if s ==None:
      s = 1+abs(ZZ.random_element())
   def randomJ(j,k):
      num_gens = 1+abs(ZZ.random_element())
      out_ideal = ideal([R.random_element() for i in range(num_gens)])
      return out_idea
   J_list = [[randomJ(j,k) for k in range(j,s)] for j in range(s-1)]
   def out_J(j,k):
      if j<k:
         return J_list[j][k-j]
      else:
         return J_list[k][j-k]
   return (R, s, out_J)


def principal_random_RsJ(n = None, s = None):
   if n == None:
      n = 1+abs(ZZ.random_element())
   R = PolynomialRing(QQ, 'y', n)
   if s ==None:
      s = 1+abs(ZZ.random_element())
   def randomJ(j,k):
      out_ideal = ideal(R.random_element())
      return out_ideal
   J_list = [[randomJ(j,k) for k in range(j,s)] for j in range(s-1)]
   def out_J(j,k):
      if j<k:
         return J_list[j][k-j]
      else:
         return J_list[k][j-k]
   return (R, s, out_J)


#example:
#(rR, rs, rJ) = random_RsJ()
#rshs = spline_Hilbert_series(rR, rs, rJ)


#example:
#(prR, prs, prJ) = principal_random_RsJ()
#rshs = spline_Hilbert_series(prR, prs, prJ)


#### expands a spline tuple as in terms of smg ####

def expand_spline(g, smg):
   R = smg[0][0].parent()
   Gsmg = ideal([tuple_to_OA(generator) for generator in smg])
   expansion_coeffs= tuple_to_OA(g).lift(Gsmg)
   OA = expansion_coeffs[0].parent()
   s = len(g)
   zero_list = [0 for i in range(s)]
   OA_to_R = OA.hom(zero_list+list(R.gens()), R)
   ec_in_R = [OA_to_R(c) for c in expansion_coeffs]
   return ec_in_R




#### these routines compares the hilbert series's of the "top term" ideal and the original ideal ###
#### we just write it for prinicpal ideals ####

def top_term(f):
   f_total_degree = f.total_degree()
   f_monomials = f.monomials()
   f_top_monomials = [ m for m in f_monomials if m.total_degree() == f_total_degree]
   f_top_term = sum([f.monomial_coefficient(m) * m for m in f_top_monomials])
   return f_top_term

def principal_top_term_ideal(J,s):
   J_list = [[J(j,k) for k in range(j,s)] for j in range(s-1)]
   def out_J(j,k):
      if j<k:
         return ideal(top_term(J_list[j][k-j].gens()[0]))
      else:
         return ideal(top_term(J_list[k][j-k].gens()[0]))
   return out_J

######## the routines below require that Barwick-Stone's QuillenSuslin package is installed  #########
def is_projective(smg):
   mac2 = Macaulay2()
   mac2.eval('loadPackage "QuillenSuslin"; ')
   gen_matrix = matrix(smg).transpose()
   mGenMatrix = mac2(gen_matrix)
   mSmod = mGenMatrix.image()
   mSmodIsFree = mSmod.isProjective()
   S_mod_is_free = mSmodIsFree.to_sage()
   return S_mod_is_free

def compute_free_basis(smg):
   mac2 = Macaulay2()
   mac2.eval('loadPackage "QuillenSuslin"; ')
   gen_matrix = matrix(smg).transpose()
   mGenMatrix = mac2(gen_matrix)
   mSmod = mGenMatrix.image()
   mFreeBasis = mSmod.computeFreeBasis()
   free_basis_matrix = mFreeBasis.to_sage()
   free_basis = free_basis_matrix.columns()
   return free_basis



######### the containment graph ##########

def my_contains(I,J):
   for g in J.gens():
      if not g in I:
         return false
   return true

def spline_graph_mat_entry(J, i1, i2, j, k):
    if (i1, i2) == (j, k) or (i1, i2) == (k, j):
       return 0
    if my_contains(J(i1,i2), J(j, k)): #the containment routine is not working correctly
       return 1
    return 0

def spline_graph(s, J, i1, i2):
    adj_list = [[spline_graph_mat_entry(J, i1, i2, j, k) for k in range(s)] for j in range(s)]
    adj_matrix = matrix(adj_list)
    out_graph = Graph(adj_matrix)
    return out_graph

def are_connected(G, i1, i2):
    am =G.adjacency_matrix()
    s = G.num_verts()
    cm = am^s
    if cm[i1,i2] ==0:
       return false
    return true


######### convert a triangle list to a "triangulation" ##########

#triangulation is an iterable of triangles
def type_triangulation(triangles):
   long_points = []
   for t in triangles:
      long_points = long_points + list(t)
   short_points = set([tuple(p) for p in long_points])
   pc = PointConfiguration(short_points)
   pc_triangulations = list(pc.triangulations())
   list_triangles = [[list(p) for p in t] for t in triangles]
   for t in list_triangles:
      t.sort()
   list_triangles.sort()
   for it_triang in pc_triangulations:
      list_it_triang = [[list(pc[p]) for p in t] for t in it_triang]
      for t in list_it_triang:
         t.sort()
      list_it_triang.sort()
      if list_it_triang ==  list_triangles:
         return it_triang
   print "no match"

