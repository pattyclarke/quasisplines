

#load spline_routines.sage

#example: no induced gluing
print 'no induced gluing'

(n,s) = (1,3)
nig_R = PolynomialRing(QQ, 'y', n)

def nig_J(j,k):
    y = nig_R.gens()[0]
    if ((j,k) == (0,1)) or ((j,k) == (1,0)):
        return y*nig_R
    elif ((j,k) == (0,2)) or ((j,k) == (2,0)):
        return (y-1)^2*nig_R
    else:
        return (y-2)^3*nig_R

nig_smg = spline_module_generators(nig_R,s,nig_J)
nig_srg = spline_ring_generators(nig_smg)

def nig_I(j,k):
	return sheet_intersection(nig_srg, j,k)

nig_test = [nig_I(j,k) == nig_J(j,k) for k in range(s) for j in range(k)]



#example:  induced gluing
print 'induced gluing'

(n,s) = (1,3)
ig_R = PolynomialRing(QQ, 'y', n)

def ig_J(j,k):
    y = ig_R.gens()[0]
    if ((j,k) == (0,1)) or ((j,k) == (1,0)):
        return y*ig_R
    elif ((j,k) == (0,2)) or ((j,k) == (2,0)):
        return (y)^2*ig_R
    else:
        return (y)^3*ig_R

ig_smg = spline_module_generators(ig_R,s,ig_J)
ig_srg = spline_ring_generators(ig_smg)

def ig_I(j,k):
	return sheet_intersection(ig_srg, j,k)

ig_test = [ig_I(j,k) == ig_J(j,k) for k in range(s) for j in range(k)]


#example: hypersurface induces non-hypersurface
print 'hypersurface induces non-hypersurface'

(n,s) = (2,3)
hnh_R = PolynomialRing(QQ, 'y', n)

def hnh_J(j,k):
    y = hnh_R.gens()
    if ((j,k) == (0,1)) or ((j,k) == (1,0)):
        return y[0]^2*hnh_R
    elif ((j,k) == (0,2)) or ((j,k) == (2,0)):
        return y[1]^2*hnh_R
    else:
        return (y[0]-y[1])^2*hnh_R

hnh_smg = spline_module_generators(hnh_R,s,hnh_J)
hnh_srg = spline_ring_generators(hnh_smg)

def hnh_I(j,k):
	return sheet_intersection(hnh_srg, j,k)


hnh_test = [(len((hnh_J(j,k)).gens()), len((hnh_I(j,k)).gens()))  for k in range(s) for j in range(k)]


#example: not-free
print 'not-free'

(n,s) = (2,2)
nf_R = PolynomialRing(QQ,'y',n)

def nf_J(j,k):
    y = nf_R.gens()
    return y * nf_R

nf_smg = spline_module_generators(nf_R, s , nf_J)
nf_smod_is_free = is_projective(nf_smg)


#example: Dalbec and Schenck counterexample to Rose's conjecture
print 'Dalbec and Schenck counterexample to Rose\'s conjecture'

(v0, v1,v2,v3, v4, v5) = ( (3,0,0), (0,3,0), (-2,-2,1), (1,0,-3), (0,0,3), (0,0,0) )
#triangle along which tetrahedra are glued is [v0, v1, v2]
#the interior vertex is v5
#if f(v) = (1,1,7)*v then the triangle is at level 3
#f(v3) = -20, v(v4) = 21, f(v5) = 0
#so v5 should be an interior point of [v0,v1,v2,v3]
#(t0, t1, t2, t3, t4) = ( (v0, v1 ,v2 ,v4), (v0, v1, v2, v5), (v0, v1, v3, v5), (v0, v2, v3, v5), (v1, v2, v3, v5) )
(t0,t1,t2,t3,t4,t5) = ((v5,v0,v1,v4),(v5,v0,v1,v3),(v5,v2,v4,v0),(v5,v2,v4,v1),(v5,v2,v3,v0),(v5,v2,v3,v1))
(dsR, dss, dsJ) = triangle_splines((t0, t1, t2, t3, t4, t5) , 0)
dssmg = spline_module_generators(dsR, dss, dsJ)
dssmod_is_free = is_projective(dssmg)
dss_mod_basis = compute_free_basis(dssmg)
ds_triangulation = type_triangulation((t0, t1, t2, t3, t4, t5))


#####the last lines below are a bit slow to compute, so i've commented them out
ds_tups=[ triangle_splines((t0, t1, t2, t3, t4, t5) , r) for r in range(7)]
ds_smgs=[spline_module_generators(ds_tup[0], ds_tup[1], ds_tup[2]) for ds_tup in ds_tups]
#ds_bool = [is_projective(ds_smg) for ds_smg in ds_smgs]
#for r in range(len(ds_smgs)):
#    print r
#    print is_projective(ds_smgs[r])


#example: Morgan-Scott triangulation - symmetric version
print 'Morgan-Scott triangulation - symmetric version'

(v0, v1,v2,v3, v4, v5) = ( (4,4), (0,0), (8,0), (4,1), (5,2), (3,2) )

(t0, t1, t2, t3, t4, t5, t6) = ((v3, v4 ,v5),
				               (v0, v4 ,v5),
				               (v1, v3 ,v5),
				               (v2, v3 ,v4),
				               (v0, v1 ,v5),
				               (v1, v2 ,v3),
	                                       (v0, v2 ,v4))

(msR, mss, msJ) = triangle_splines((t0, t1, t2, t3, t4, t5, t6) , 2)
mssmg = spline_module_generators(msR, mss, msJ)
mssmod_is_free = is_projective(mssmg)
mss_mod_basis = compute_free_basis(mssmg)
ms_triangulation = type_triangulation((t0, t1, t2, t3, t4, t5, t6))
   
#example: Morgan-Scott triangulation - skew version
print 'Morgan-Scott triangulation - skew version'

(v0, v1,v2,v3, v4, v5) = ( (0,0), (4,0), (2,4), (1,1), (3,1), (2,3) )

(t0, t1, t2, t3, t4, t5, t6) = ((v3, v4 ,v5),
				               (v0, v3 ,v5),
				               (v1, v3 ,v4),
				               (v2, v4 ,v5),
				               (v0, v1 ,v3),
				               (v1, v2 ,v4),
	                                       (v0, v2 ,v5))

(ms2R, ms2s, ms2J) = triangle_splines((t0, t1, t2, t3, t4, t5, t6) , 2)
ms2smg = spline_module_generators(ms2R, ms2s, ms2J)
ms2smod_is_free = is_projective(ms2smg)
ms2s_mod_basis = compute_free_basis(ms2smg)
ms2_triangulation = type_triangulation((t0, t1, t2, t3, t4, t5, t6))


#example: non-manifold triangulation without free splines - from Billera-Rose
print 'non-manifold triangulation without free splines - from Billera-Rose'

(v0, v1,v2,v3, v4) = ( (-1,0), (0,1), (0,0), (1,0), (0,-1) )
(t0, t1) = ((v0,v1,v2),(v2,v3,v4))
(nmR, nms, nmJ) = triangle_splines((t0, t1), 0)
nmsmg = spline_module_generators(nmR, nms, nmJ)
nmsmod_is_free = is_projective(nmsmg)
#the triangulation look up here doesn't work, because the triangulation is not of the convex hull of the points
#nm_triangulation = type_triangulation((t0, t1, t2, t3, t4, t5, t6))



#example: random triangulation
print 'random triangulation' 

(R, s, J, triangulation, r) = random_triangle_splines(3,2)
smg = spline_module_generators(R, s, J)
sm_is_projective = is_projective(smg)
time fb = compute_free_basis(smg)
triangulation.plot()










#example: Billera-Rose, Modules of piecewise polynomials and their freeness
#symmetric arrangment
print 'Billera-Rose, Modules of piecewise polynomials and their freeness'

(v0, v1, v2, v3, v4, v5, v6) = ((0,0,0),(1,0,0),(0,1,0), (0,0,1),(-1,0,0),(0,-1,0),(0,0,-1))
(t0, t1, t2, t3, t4, t5, t6, t7) = (
(v0,v1,v2,v3),
(v0,v1,v2,v6),
(v0,v1,v5,v3),
(v0,v1,v5,v6),
(v0,v4,v2,v3),
(v0,v4,v2,v6),
(v0,v4,v5,v3),
(v0,v4,v5,v6),
)

(brR, brs, brJ) = triangle_splines((t0, t1, t2, t3, t4, t5, t6, t7) , 2)
brsmg = spline_module_generators(brR, brs, brJ)
brsmod_is_free = is_projective(brsmg)
brs_mod_basis = compute_free_basis(brsmg)
br_triangulation = type_triangulation((t0, t1, t2, t3, t4, t5, t6, t7))

#deformed Billera-Rose arrangement
print 'deformed Billera-Rose arrangement'



(fbr2R, fbr2s, fbr2J) = triangle_splines((t0, t1, t2, t3, t4, t5, t6, t7) , 2)



#Alfeld split
print 'Alfeld split'

def alfeld_triangles(n):
    v = tuple([tuple([0 for i in range(n)])]+[tuple([-1 for i in range(n)])]+identity_matrix(n).columns())
    t = tuple([v[:i]+ v[i+1:] for i in range(n+2)[1:]])    
    return t

def alfeld(n,r):
    (asR, ass, asJ) = triangle_splines(alfeld_triangles(n),r)    
    return (asR, ass, asJ)

as2_triang = alfeld_triangles(2)

as2_triangulation = type_triangulation(as2_triang)

(asR20, ass20, asJ20) = triangle_splines(alfeld_triangles(2),0)

(asR21, ass21, asJ21) = triangle_splines(alfeld_triangles(2),1)
as21smg = spline_module_generators(asR21, ass21, asJ21)
as21mod_is_free = is_projective(as21smg)
as21_mod_free_basis = compute_free_basis(as21smg)




