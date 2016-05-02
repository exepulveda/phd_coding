import numpy as np
import nsr


if __name__ == "__main__":
    grade_cu = 0.12314
    grade_au = 0.1662
    recovery_cu = 74.18866
    recovery_au = 81.45106
    concentrate_cu = 19.962952
    concentrate_au = 60.92745
    concentrate_f = 1450.2705
    
    price_cu = 3.38 #$/lb
    price_to_ton_cu = 2204.6 
    pon_to_ton = 453.59237
    price_au = 1562.50 #US$oz
    price_to_ton_au = 0.0321543408
    
    deduct_cu =1.0 #%
    deduct_proportion_cu = 0.965
    deduct_au = 1.3 #g/t conc
    deduct_proportion_au = 0.975
    
    treatment_cost_cu = 100.0
    
    refining_cost_cu = 0.1 #US$/lb
    refining_cost_au = 8.75 #US$/oz
    
    refining_charges = 0.00 #US$/lb
    
    site_costs = 1.6 #$/t
    
    penalties = 0
    verbose = False
    marketing_charge = 0.0 #%
    treatment_distribution_charge = 0.0

    transportation_costs = 100.388
    moisture = 0.085
    marine_insurence = 0.0175 / 100.0
    milling_costs = 8 # $/t
    transportation_losses = 0.008
    #royalty
    royalty = 0.04
    admi_costs_royalty = 1.6
    allowable_depreciation_royalty = 1.09

    min_concentrate_cu = 22.0
    penalty_min_concentrate_cu = 7.5

    ret = nsr.calculateNSRNewcrest(
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
        [300,700],
        [100,100],
        [1.25,1.875],
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
        True)
