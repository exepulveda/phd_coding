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

def calculateNSRNewcrest(
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
    treatment_cost_cu,
    #cost au
    deduct_au,
    deduct_proportion_au,
    refining_cost_au,
    #penalties realted low copper
    min_concentrate_cu,
    penalty_min_concentrate_cu,
    #penalties related F
    concentrate_f,
    penalty_above,
    penalty_each,
    penalty_f,
    #others
    transportation_costs, # $/t
    moisture,
    marine_insurence,
    milling_costs, # $/t
    transportation_losses, # %
    #royalty
    royalty, #%
    admi_costs_royalty, # $/t
    allowable_depreciation_royalty, # $/t
    
    verbose = False):
        
    #calculate payable or smelter metal content
    
    #au_smelter_grade = grade_au*recovery_au * deduct_proportion_au * (1.0 - transportation_losses)
    #cu_smelter_grade = (grade_cu*recovery_cu * 100.0 * 1000.0) * deduct_proportion_cu * (1.0 - transportation_losses)
    
    
    
    #payable_cu = (concentrate_cu - deduct_cu) * deduct_proportion_cu
    #payable_metal_cu = payable_cu * price_cu * price_to_ton_cu / 100.0

    #payable_au = (concentrate_au - deduct_au) * deduct_proportion_au 
    #payable_metal_au = payable_au * price_au * price_to_ton_au         
        
    TON = 342.50
        
    #deductions
    au_fm = grade_au * recovery_au / 100.0 * TON * price_to_ton_au
    au_sfm = au_fm * deduct_proportion_au * (1.0 - transportation_losses)
    au_sg = au_sfm / TON / price_to_ton_au
    au_fg = au_fm / price_to_ton_au / TON
    au_ded = (au_fg - au_sg) * (price_au * price_to_ton_au)
    
    if verbose:
        print "au_fm",au_fm
        print "au_sfm",au_sfm
        print "au_sg",au_sg
        print "au_fg",au_fg
        print "au_ded",au_ded
    
    

    cu_str = (concentrate_cu - deduct_cu) / concentrate_cu
    cu_fm = grade_cu * (recovery_cu / 100.0) * TON / 100.0
    
    
    cu_sfm =  cu_fm* cu_str * (1.0-transportation_losses)
    
    cu_sg = cu_sfm*1000000/TON

    cu_fg = cu_fm* 1000000/TON
    cu_ded = (cu_fg - cu_sg) * price_cu / 453.59237

    if verbose:
        print "cu_str",cu_str
        print "cu_sfm",cu_sfm
        print "cu_sg",cu_sg
        print "cu_fm",cu_fm
        print "cu_fg",cu_fg
        print "cu_ded",cu_ded

    
    #smelter charges
    #con_dry = ((grade_cu*(recovery_cu/100.0)))/concentrate_cu/100.0
    con_dry = ((grade_cu*10000 *(recovery_cu/100.0))/10000)/(concentrate_cu)
    sme_dry = (con_dry-0.0/1000)*(1.0-transportation_losses)  
    smet_cst = sme_dry * treatment_cost_cu

    if verbose:
        print "con_dry",con_dry
        print "sme_dry",sme_dry
        print "smet_cst",smet_cst

    
    #fluorine penalties
    fpt_cst = None
    
    #low copper
    cuc_pen = (0 if concentrate_cu > min_concentrate_cu else penalty_min_concentrate_cu) * sme_dry #penalty is in $/dwt
    if verbose:
        print "cuc_pen",cuc_pen
    
    #refining charges
    au_ref =  au_sg * refining_cost_au * price_to_ton_au
    cu_ref =  cu_sg * refining_cost_cu / 453.59237
    ref_cst =  au_ref + cu_ref
    if verbose:
        print "au_ref",au_ref
        print "cu_ref",cu_ref
        print "ref_cst",ref_cst
    
    #revenues
    au_rev = au_fg * price_au * price_to_ton_au 
    cu_rev = cu_fg * price_cu / 453.59237
    rev = au_rev + cu_rev

    if verbose:
        print "au_rev",au_rev
        print "cu_rev",cu_rev
        print "rev",rev

    
    #conc_cst
    tra_miav = (au_rev + cu_rev) * marine_insurence
    ttracst = transportation_costs*(con_dry)/(1.0-moisture)
    conc_cst = tra_miav + ttracst
    
    #realisation costs
    real_cst = conc_cst + ded_cst + smet_cst + fpt_cst + cuc_pen + ref_cst # $/t

    
    #royalty costs
    roya_cst = max(0,(rev - milling_costs - real_cst - 1.0/3.0* admi_costs_royalty - allowable_depreciation_royalty) * royalty)

    #underground value
    underground_value = rev - (real_cst + roya_cst)
        

    return underground_value


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
