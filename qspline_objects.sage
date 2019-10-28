
import quasisplines as qs
import operator

#######################################################################
# OBJECTS
#######################################################################


# Ideal Difference Conditions Object #

class IdealDifferenceConditions:

    def __init__(self, s, J):
        self.sheet_number = s 
        self.coefficient_ring = J(0,0).ring()
        self.difference_ideal = J #function from range(s)^2 into ideals of the coefficient ring CR

    ### attributes and setters ###
        
    sheet_number = property(operator.attrgetter('_sheet_number'))

    @sheet_number.setter
    def sheet_number(self, s):
        if not (s.is_integer() and s >0): raise Exception("sheet_number must be a positive integer")
        self._sheet_number = s

    coefficient_ring = property(operator.attrgetter('_coefficient_ring'))

    @coefficient_ring.setter
    def coefficient_ring(self, R):
	if not R.is_ring(): raise Exception("The coefficients must form a ring")
        self._coefficient_ring = R 

    difference_ideal = property(operator.attrgetter('_difference_ideal'))

    @difference_ideal.setter
    def difference_ideal(self, J):
        J_dom = [(i,j) for i in range(self.sheet_number) for j in range(self.sheet_number)]
        for t in J_dom:
            if not (J*(t) == self.coefficient_ring.ideal(J(*t).gens())):
                raise Exception("The difference ideal function fails.")
        self._difference_ideal = J

    ### methods ###
        
    def is_spline(self, g):
        J = self.difference_ideal
        s = len(g)
        for p in [(i,j) for j in range(s) for i in range(j)]:
            if not (g[p[0]] - g[p[1]] in J(p[0],p[1])):
                return False
        return True

    def get_quasi_spline_ring(self):
        QSR = QSR_from_IDC(self)
        return QSR
        
# Quasi-Spline Module Object #                                                                                         
class QuasiSplineModule:

    def __init__(self, smg):

        self.coefficient_ring = smg[0][0].parent()
        
        if self._good_gens(smg):
            self._gens = smg
        else: raise Exception("Bad spline module generators.")

    ### attributes and setter ###

    coefficient_ring = property(operator.attrgetter('_coefficient_ring'))

    @coefficient_ring.setter
    def coefficient_ring(self, CR):
        if not CR.is_ring(): raise Exception("The coefficients must form a ring")
        self._coefficient_ring = CR

    ### checker ###

    def _good_gens(self, smg):
        s = len(smg[0])
        for g in smg: 
            if not len(g) == s: raise Exception("Bad module generators.")
            for f in g:
                if not f in self.coefficient_ring: raise Exception("Bad module generators.")
        return True
        
    ### methods ###
        
    def gens(self):
        return self._gens
    
    def quasi_spline_ring(self):
        return QuasiSplineRing(QSM = self)

    def is_projective(self):
        return qs.is_projective(self.gens())

    def get_free_basis(self):
        if self.is_projective():
            return qs.compute_free_basis(self.gens())
        else:
            print("Module is not free.\n")
            return self.gens()
        
        
# Quasi-Spline Ring Object #

class QuasiSplineRing:

    def __init__(self, QSM = None):

        smg = QSM.gens()
        srg = qs.spline_ring_generators(smg)
        CR = QSM.gens()[0][0].parent()

        (QSR, qsv, cv) = qs.SPxy(srg)
	RI = qs.SP_relation_ideal(srg)
 
        self.ring = QSR
        self.relation_ideal = RI
        self.generating_tuples = srg
        

        self.quasi_spline_variables = qsv #x variables in the code
        

        self._module = QSM
        self._coefficient_ring = CR
        self._coefficient_variables = cv #y variables in the code
        
    ### attributes and setters ###

    ring = property(operator.attrgetter('_ring'))

    @ring.setter
    def ring(self, QSR):
        #The 'SP' or in the code
        self._ring = QSR

    relation_ideal = property(operator.attrgetter('_relation_ideal'))

    @relation_ideal.setter
    def relation_ideal(self, RI):
        self._relation_ideal = RI

    generating_tuples = property(operator.attrgetter('_generating_tuples'))

    @generating_tuples.setter    
    def generating_tuples(self, srg):
        self._generating_tuples = srg

    quasi_spline_variables = property(operator.attrgetter('_quasi_spline_variables'))

    @quasi_spline_variables.setter
    def quasi_spline_variables(self, qsv):
        #The list of 'x' aka spline variables in the code
        self.quasi_spline_variables = qsv

    ### methods ###
    
    def module(self):
        return self._module
        
    def coefficient_ring(self):
        #The 'R' ring in the code
        return self._coefficient_ring

    def coefficient_variables(self):
        #The list of 'y' variables aka coefficient variables in the code
        return self._coefficient_variables

    def quasi_spline(self, inp):
        r'''
        This method returns a QuasiSpline object. As input, it can take either a tuple or an element of the self.ring.
        '''
        if inp.__class__ == tuple:
            tup_in_SP = qs.tuple_to_SP(self.generating_tuples, inp)
            qsob = QuasiSpline(self, tup_in_SP)
            return qsob

        R = self.ring
        return QuasiSpline(self, R(inp))

    def random_quasi_spline(self):
        DSSPrand = self.ring.random_element()
        outqs = QuasiSpline(self, DSSPrand)
        return outqs

    def contact_ideals(self):
        ctis = IDC_from_QSR(self)
        return ctis

    def Hilbert_polynomial(self):
        return Hilbert_poly_from_QSR(self)
    
    def Hilbert_series(self):
        return Hilbert_series_from_QSR(self)
    
# Quasi-Spline Object #
    
class QuasiSpline(RingElement):

    def __init__(self, QSR, qusp):
        self.parent = QSR #A quasi-spline ring object for which qs in an element of QSR.ring
        self.ring_element = qusp

    ### attributes and setters ###

    parent = property(operator.attrgetter('_parent'))

    @parent.setter
    def parent(self, QSR):
        self._parent = QSR
        
    ring_element = property(operator.attrgetter('_ring_element'))
    
    @ring_element.setter
    def ring_element(self, qusp):
        if  not (qusp in self.parent.ring): raise Exception("The quasi-spline must be in the quasi-spline ring.") 
        self._ring_element = qusp
        
    ### methods ###
    
    def tuple(self):
        DSSP = self.parent.ring
        DSsrg = self.parent.generating_tuples
        DSSPrand = self.ring_element
        tup = qs.SP_to_tuple(DSSP, DSsrg, DSSPrand)
        return tup

    def _add_(self, other):
        sum_ob = QuasiSpline(self.parent, self.ring_element + other.ring_element)
        return sum_ob
    
    def _sub_(self, other):
        diff_ob = QuasiSpline(self.parent, self.ring_element - other.ring_element)
        return diff_ob

    def _mul_(self, other):
        prod_ob = QuasiSpline(self.parent, self.ring_element * other.ring_element)
        return prod_ob

    def _div_(self, other):
        quot_ob = QuasiSpline(self.parent, self.ring_element / other.ring_element)
        return quot_ob

    def __pow__(self, n):
        pow_ob = QuasiSpline(self.parent, (self.ring_element)^n)
        return pow_ob
    
#######################################################################
# QuasiSplineModule = function(IdealDifferenceCondition)
#######################################################################

def QSM_from_IDC(IDC = None):

    #Rudimentary input check
    if IDC == None:
        return None

    J = IDC.difference_ideal
    s = IDC.sheet_number
    CR = J(0,0).ring()

    smg = qs.spline_module_generators(CR, s, J)

    QSM = QuasiSplineModule(smg)
    return QSM


#######################################################################
# QuasiSplineRing = function(QuasiSplineModule)
#######################################################################

# Just use the constructor
def QSR_from_QSM(QSM = None):

    #Rudimentary input check
    if QSM == None:
        return None

    QSR = QuasiSplineRing(QSM)
    return QSR


#######################################################################                                               
#IdealDifferenceCondition = function(QuasiSplineRing)
#######################################################################  

def IDC_from_QSR(QSR = None):

    #Rudimentary input check
    if QSR == None:
        return None

    QSM = QSR.module()

    J = qs.contact_ideals(QSM.gens())
    s = len(QSM.gens()[0])

    IDC = IdealDifferenceConditions(s,J)
    return IDC


######################### Above we did a loop ###################################
######################### Below we do the rest of the loop ######################

#######################################################################
# QuasiSplineRing = function(IdealDifferenceCondition)
#######################################################################

def QSR_from_IDC(IDC = None):

    #Rudimentary input check
    if IDC == None:
        return None

    QSM = QSM_from_IDC(IDC)
    QSR = QuasiSplineRing(QSM)
    return QSR

#######################################################################
# IdealDifferenceCondition = function(QuasiSplineModule)
#######################################################################

def IDC_from_QSM(QSM = None):

    #Rudimentary input check
    if QSM == None:
        return None

    QSR = QSR_from_QSM(QSM)
    IDC = IDC_from_QSR(QSR)
    return IDC

#######################################################################
# QuasiSplineModule = function(QuasiSplineRing)
#######################################################################

def QSM_from_QSR(QSR):
    QSM = QSR.module()
    return QSM


#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################


def IDC_from_TL(TL, r):
    (R, s, J) = qs.triangle_splines(TL, r)
    IDC = IdealDifferenceConditions(s,J)
    return IDC

def Hilbert_series_from_QSR(QSR):
    #Returns the Hilbert Series of the homogenized ring
    IDC =  IDC_from_QSR(QSR)
    
    R = QSR.coefficient_ring()
    s = IDC.sheet_number
    J = IDC.difference_ideal
    
    Hs = qs.spline_Hilbert_series(R, s, J)
    return Hs

def Hilbert_poly_from_QSR(QSR):
    #Returns the Hilbert Series of the homogenized ring
    IDC =  IDC_from_QSR(QSR)
    
    R = QSR.coefficient_ring()
    s = IDC.sheet_number
    J = IDC.difference_ideal
    
    Hp = qs.spline_Hilbert_polynomial(R, s, J)
    return Hp

