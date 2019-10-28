
### separator line ###

# For profiling: comment out the 'from' line below, move up a directory and run
#
# python -m cProfile tests/test_quasisplines.py
#

from .context import quasisplines 


import quasisplines.quasisplines as qs
import quasisplines.qspline_objects as qsob


def test_Dalbec_Schenck():
    # This example is Dalbec and Schenck counterexample to Rose's conjecture.
    # It is comprised for two tetrahedra in 3-space glued along a triangle.
    # The triangle along which tetrahedra are glued is [v0, v1, v2]. 
    (v0, v1,v2,v3, v4, v5) = ( (3,0,0), (0,3,0), (-2,-2,1), (1,0,-3), (0,0,3), (0,0,0) ) 
    (t0,t1,t2,t3,t4,t5) = ((v5,v0,v1,v4),(v5,v0,v1,v3),(v5,v2,v4,v0),(v5,v2,v4,v1),(v5,v2,v3,v0),(v5,v2,v3,v1))

    # Now, we create the ideal-difference condition.
    DSIDCob = qsob.IDC_from_TL((t0, t1, t2, t3, t4, t5) , 0)

    # A sample of relevant information extracted from the ideal-difference condition:
    R = QQ['y0, y1, y2'] 
    assert DSIDCob.coefficient_ring == R
    assert DSIDCob.sheet_number == 6
    y = R.gens()
    assert DSIDCob.difference_ideal(0,1) == R.ideal(9*y[2])

    DSQSRob = DSIDCob.get_quasi_spline_ring()

    QSR = DSQSRob.ring()

    assert len(DSQSRob.ring().args()) == 8
    assert len(DSQSRob.coefficient_variables()) == 3

    DSHp = DSQSRob.Hilbert_polynomial()
    DSHs = DSQSRob.Hilbert_series()

    d = DSHp.parent().gen()
    t = DSHs.parent().gen()
    
    assert DSHp == 3*d^2 + 2
    assert DSHs == (-t^3 - 2*t^2 - 2*t - 1)/(t^3 - 3*t^2 + 3*t - 1)

    DSQSMob = DSQSRob.module()
    DSsmg = DSQSMob.gens()
    
    contact_J_ob = DSQSRob.contact_ideals()
    contact_J = contact_J_ob.difference_ideal
    DSJ = DSIDCob.difference_ideal

    for j in range(len(DSsmg[0])):
        for i in range(len(DSsmg[0])):
            assert qs.my_contains(DSJ(i,j), contact_J(i,j))

    DSQSob = DSQSRob.random_quasi_spline()
    oDSQSob = DSQSRob.random_quasi_spline()

    sumDSQSob = DSQSob + oDSQSob
    diffDSQSob = DSQSob - oDSQSob
    prodDSQSob = DSQSob * oDSQSob

    assert [DSQSob.tuple()[i] + oDSQSob.tuple()[i] - sumDSQSob.tuple()[i] for i in range(len(sumDSQSob.tuple()))] == [0,0,0,0,0,0]

    DSQSRrel = DSQSRob.relation_ideal
    DSsp_gens = [DSQSRob.quasi_spline(f) for f in DSQSRrel.gens()]

    for f in DSsp_gens:
        assert f.tuple() == (0,0,0,0,0,0)

    tupQSob = DSQSob.tuple()
    churnQSob = DSQSRob.quasi_spline(DSQSRob.quasi_spline(tupQSob).ring_element)

    erpoly  = DSQSob.ring_element - churnQSob.ring_element
    assert (erpoly in DSQSRob.relation_ideal, DSQSRob.quasi_spline(erpoly).tuple()) == (True, (0,0,0,0,0,0))


if __name__ == "__main__":
    # print("__name__ is __main__!")
    test_Dalbec_Schenck()
    





    
