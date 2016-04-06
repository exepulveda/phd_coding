import numpy as np

lbToTon = 2204.6
ounceToTon = 1.0/31.1

def calc_dmt(
        ore_grade,
        recovery,
        concentrate_grade):
    return ore_grade * recovery / concentrate_grade / 100.0

def calc_payable_ore(
        ore_grade,
        recovery,
        concentrate_grade,
        deduction):
    
    return ore_grade * recovery * lbToTon * (concentrate_grade - deduction) / concentrate_grade

def refining_charges(
        refining_charges,
        payable_ore
        ):

    return payable_ore  * refining_charges
    
def calculate_nsr_cut():

    dmt = calc_dmt(ore_grade,recovery,concentrate_grade) #U
    
    payable_ore = calc_payable_ore(ore_grade,recovery,concentrate_grade,deduction) #W
    
    #charges
    smelter_charges_ore = smelter_charges / dmt #AM
    refining_charges_ore = refining_charges * payable_ore #AR

def calculateNSR(
        grade_cu,
        grade_au,
        recovery_cu,
        recovery_au,
        concentrate_cu,
        concentrate_au,
        price_cu, #US$/ton
        price_au, #US$/ton
        deduct_cu,
        deduct_au,
        refining_cost_cu,
        refining_cost_au,
        smelter_charges,
        refining_charges,
        freight_charges,
        penalties,
        marketing_charge,
        treatment_distribution_charge,
        verbose = False
        ):
    #main metal
    #DMT
    dmt = grade_cu * (recovery_cu / 100.0) / concentrate_cu
    
    payable_cu = grade_cu * (recovery_cu /100.0) * lbToTon * (concentrate_cu - deduct_cu)/concentrate_cu / 100.0 

    #secondary metal
    payable_pct = np.min((concentrate_au-deduct_au)/concentrate_au,0.95) * 100.0
    payable_au = grade_au*(payable_pct / 100.0)*(recovery_au/100.0)*0.03215
    
    #smelter charge
    smelter_charges = smelter_charges * dmt
    
    #refining charges
    refining_cu = (refining_cost_cu + np.max((price_cu - 0.9)*0.1,0.0))
    refining_au = refining_cost_au
    
    refining_charges = refining_cu * payable_cu + refining_au * payable_au
    
    total_treatment_charges = smelter_charges + refining_charges + penalties * dmt
    
    freight = freight_charges * dmt * 1.075
    
    #revenues
    revenue_cu = payable_cu * price_cu
    revenue_au = payable_au * price_au
    
    total_revenue = revenue_cu + revenue_au
    
    net_revenue = total_revenue * (100.0 - marketing_charge)/ 100.0
    
    net_smelter_return = net_revenue - treatment_distribution_charge

    return net_smelter_return


def calculateNSRLane(
    #cu
    grade_cu,
    recovery_cu,
    concentrate_cu,
    price_cu,
    price_to_ton_cu,
    #au
    grade_au,
    recovery_au,
    concentrate_au,
    price_au,
    price_to_ton_au,
    #cost cu
    deduct_cu,
    deduct_proportion_cu,
    refining_cost_cu,
    #cost au
    deduct_au,
    deduct_proportion_au,
    refining_cost_au,
    #common costs
    smelter_charges,
    freight_charges,
    #penalties related F
    concentrate_f = None,
    penalty_above = None,
    penalty_each = None,
    penalty_f = None,
    #others
    verbose = False):
    #concentrate gross
    concentrate_gross_cu = (concentrate_cu / 100.0) * price_cu * price_to_ton_cu #dividing 100 %units
    concentrate_gross_au = (concentrate_au * price_au * price_to_ton_au)
    total_concentrate_gross = concentrate_gross_cu + concentrate_gross_au
    
    #dmt
    dmt_cu = grade_cu * (recovery_cu/100.0) / concentrate_cu
    dmt_au = grade_cu * (recovery_cu/100.0) / concentrate_cu
    
    #calculate payable
    
    payable_cu = (concentrate_cu - deduct_cu) * deduct_proportion_cu
    payable_metal_cu = payable_cu * price_cu * price_to_ton_cu / 100.0

    payable_au = (concentrate_au - deduct_au) * deduct_proportion_au 
    payable_metal_au = payable_au * price_au * price_to_ton_au 

    total_payable = payable_metal_cu + payable_metal_au
    #print "totalPayable",totalPayable

    #penalties
    if penalty_f is not None:
        penalties_f = np.zeros_like(concentrate_f)
        '''above 300 ppm penalties of $1.50 per 100 ppm will apply'''
        indices = np.where(concentrate_f > penalty_above)
        #print len(indices),concentrate_f.shape,penalty_each,penalty_f,
        penalties_f[indices] = np.floor(concentrate_f[indices]/penalty_each) * penalty_f
    else:
        penalties_f = 0.0


    gross_payable_total = total_payable - (smelter_charges + freight_charges + penalties_f)
    
    #nsr per metal
    proportion_payable_cu = payable_metal_cu / total_payable
    gross_payable_cu = proportion_payable_cu * gross_payable_total
    refining_value_cu = refining_cost_cu * price_to_ton_cu * (payable_cu / 100.0)

    proportion_payable_au = payable_metal_au / total_payable
    gross_payable_au = proportion_payable_au * gross_payable_total
    refining_value_au = refining_cost_au * price_to_ton_au * payable_au
    
        
    #print gross_payable_cu.shape
    #print refining_value_cu.shape
    #print penalties_f.shape


    nsr_concentrate_cu = gross_payable_cu - refining_value_cu - penalties_f


    nsr_concentrate_au = gross_payable_au - refining_value_au
    
    nsr_ore_cu = nsr_concentrate_cu / concentrate_cu * (recovery_cu / 100.0) * grade_cu
    nsr_ore_au = nsr_concentrate_au / concentrate_au * (recovery_au / 100.0) * grade_au

    return nsr_ore_cu,nsr_concentrate_cu,nsr_ore_au,nsr_concentrate_au

    
def calculateNSR2(
        grade,
        gradeConvertion,
        recovery,
        concentrate,
        price,
        deduct,
        proportion,
        refiningCost,
        refiningCostConvertion,
        smelterCharges,
        freightCharges,
        verbose = False
        ):
    #concentrate gross
    concentrateGross = concentrate * price / gradeConvertion
    
    #print "concentrateGross",concentrateGross
    
    totalConcentrateGross = np.sum(concentrateGross,axis=1)
    
    #assert (totalConcentrateGross.shape[0] == concentrate.shape[0])
    
    if verbose:
        print concentrate[:5]
        print grade[:5]
        print recovery[:5]
    
    #calculate payable
    oreToConcentrate = concentrate / (grade * recovery)
    payable = (concentrate - deduct) * proportion 
    #print "(concentrate - deduct) * proportion ",concentrate,deduct,proportion

    payableMetal = payable * price / gradeConvertion


    if verbose:
        print "oreToConcentrate",oreToConcentrate[:5]
        print "payable",payable[:5]
        print "payableMetal",payableMetal[:5]

    #print "payable",payable,"payableMetal",payableMetal
    totalPayable = np.sum(payableMetal,axis=1)
    #print "totalPayable",totalPayable

    grossPayableTotal = totalPayable - (smelterCharges + freightCharges)
    
    #print "grossPayableTotal",grossPayableTotal

    #nsr per metal
    proportionPayable = np.empty_like(payableMetal)
    grossPayable = np.empty_like(proportionPayable)
    n = proportionPayable.shape[1]
    for i in xrange(n):
        proportionPayable[:,i] = payableMetal[:,i] / totalPayable
        grossPayable[:,i] = proportionPayable[:,i] * grossPayableTotal
    
    #print "grossPayable",grossPayable
    
    refiningValue = refiningCost * refiningCostConvertion

    if verbose:
        print "grossPayable",grossPayable[:5]
        print "refiningValue",refiningValue[:5]

    nsrConcentrate = grossPayable - refiningValue
    #print "nsrConcentrate",nsrConcentrate
    
    nsrOre = nsrConcentrate / oreToConcentrate

    nsrTotal = np.sum(nsrOre,axis=1)

    if verbose:
        print "nsrConcentrate",nsrConcentrate[:5]
        print "nsrOre",nsrOre[:5]
        print "nsrTotal",nsrTotal[:5]



    return nsrTotal,nsrOre,nsrConcentrate

if __name__ == "__main__":
    grade_cu = 0.17
    grade_au = 0.73
    recovery_cu = 78.0
    recovery_au = 72.5
    concentrate_cu = 25.0
    concentrate_au = 100.0
    price_cu = 1.7 #US$/lb
    price_to_ton_cu = 2204.6
    price_au = 650.0 #US$oz
    price_to_ton_au = 0.03215
    deduct_cu =1.0
    deduct_proportion_cu = 1.0
    deduct_au = 1.0
    deduct_proportion_au = 0.99
    refining_cost_cu = 0.08 #US$/lb
    refining_cost_au = 6.0 #US$/oz
    smelter_charges = 80.0
    refining_charges = 0.00 #US$/lb
    freight_charges = 35.0 #US$/ton of con
    penalties = 0
    verbose = False
    marketing_charge = 0.0 #%
    treatment_distribution_charge = 0.0

    ret = calculateNSRLane(
    #cu
    grade_cu,
    recovery_cu,
    concentrate_cu,
    price_cu,
    price_to_ton_cu,
    #au
    grade_au,
    recovery_au,
    concentrate_au,
    price_au,
    price_to_ton_au,
    #cost cu
    deduct_cu,
    deduct_proportion_cu,
    refining_cost_cu,
    #cost au
    deduct_au,
    deduct_proportion_au,
    refining_cost_au,
    #common costs
    smelter_charges,
    freight_charges,
    verbose = False)

    print ret
    quit()


    grade_cu = 0.17
    grade_au = 0.73
    recovery_cu = 78.0
    recovery_au = 72.5
    concentrate_cu = 25.0
    concentrate_au = 100.0
    price_cu =0.625 #US$/lb
    price_to_ton_cu = 2204.6
    price_au =265 #US$oz
    price_to_ton_au = 0.03215
    deduct_cu =1.0
    deduct_proportion_cu = 1.0
    deduct_au = 1.0
    deduct_proportion_au = 0.99
    refining_cost_cu = 0.10 #US$/lb
    refining_cost_au = 6.0 #US$/oz
    smelter_charges = 110.0
    refining_charges = 0.00 #US$/lb
    freight_charges = 80.0 #US$/ton of con
    penalties = 0
    verbose = False
    marketing_charge = 0.0 #%
    treatment_distribution_charge = 0.0

    ret = calculateNSRLane(
    #cu
    grade_cu,
    recovery_cu,
    concentrate_cu,
    price_cu,
    price_to_ton_cu,
    #au
    grade_au,
    recovery_au,
    concentrate_au,
    price_au,
    price_to_ton_au,
    #cost cu
    deduct_cu,
    deduct_proportion_cu,
    refining_cost_cu,
    #cost au
    deduct_au,
    deduct_proportion_au,
    refining_cost_au,
    #common costs
    smelter_charges,
    freight_charges,
    verbose = False)

    print ret
    quit()
    ret = calculateNSR(
        grade_cu,
        grade_au,
        recovery_cu,
        recovery_au,
        concentrate_cu,
        concentrate_au,
        price_cu,
        price_au,
        deduct_cu,
        deduct_au,
        refining_cost_cu,
        refining_cost_au,
        smelter_charges,
        refining_charges,
        freight_charges,
        penalties,
        marketing_charge,
        treatment_distribution_charge,
        verbose = False)

    print ret
