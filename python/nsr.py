TONS2POUNDS = 2204.6

def nsr_calculation(
    copper_grade,
    copper_recovery,
    copper_concentrate_recovery,
    copper_deductions,
    copper_refining_charges,
    copper_price,
    gold_grade,
    gold_concentrate_recovery,
    gold_deductions,
    gold_refining_charges,
    gold_price,
    smelter_charge,
    penalties,
    freights,
    marketing_charges,
    moisture = 0.07
):

    A = copper_grade 
    B = copper_recovery
    
    M = gold_concentrate_recovery
    
    T = copper_concentrate_recovery
    
    L = 0.0 #gravity gold recovery
    
    #prices
    DU = copper_price
    DX = gold_price
  
    DZ =  1.0
    
    #DMT
    U = (A * B/T) / 100.0
    print "U",U
    
    V = copper_deductions
    #Payable cooper
    W = A*B*TONS2POUNDS*(T-V)/T / 100.0
    print "W",W

    #gold
    K = gold_grade
    AA = K #  / 0.0292 if au is onces
    #gold in dry concentrate
    AB = AA * M/U 
    print "AB",AB
    AC = gold_deductions
    #payable gold %
    if (AC > AB):
        AD = 0
    else:
        AD = min((AB-AC)/AB,0.95)
    print "AD",AD

    #payable gold oz/t of ore
    AE = AA * AD * (M) * 0.03215
    print "AE",AE
    
    #Smelter charge concentrate
    AL = smelter_charge
    #Smelter charge / t of ore
    AM = AL*U
    print "AM",AM
    
    #Refinig charges
    AN = (copper_refining_charges+ max((DU-0.9)*0.1,0))/DZ
    print "AN",AN
    AP = gold_refining_charges
    #Refining charges /t ore
    AR = (AN*W)+(AP*AE)
    print "AR",AR
    #Penalties
    AS = penalties
    
    #total treatment charges /t ore
    AT = AM + AR + (AS*U)
    AU = freights
    #freights /t ore
    AV = AU*U*(1.0 + moisture)
    
    #treatement charges total  /t ore
    DM = AT
    #freight charges total  /t ore
    DN = AV
    #total
    DO = DM + DN
    
    #Payable metals
    #copper
    DP = W
    #gold
    DS = AE # [GRAVITY + (K*L*TONS2POUNDS/2000)] # + BK + CS
    
    print "DS",DS
    
    #revenues
    EA = DP * DU/DZ #copper
    ED = DS * DX/DZ #gold
    
    EF = EA + ED #total
    
    #marketing
    EG = marketing_charges
    
    #net revenue /t ore
    EH = EF * ((100.0 - EG)/100.0)
    
    #NSR /t of ore
    EI = EH - DO
    
    return EI
        
if __name__ == "__main__":
    copper_grade = 1.259
    copper_recovery = 84.8 / 100.0
    copper_concentrate_recovery = 24.0  / 100.0
    copper_deductions = 1.06  / 100.0
    copper_refining_charges = 0.09
    copper_price = 1.21
    gold_grade = 0.788
    gold_concentrate_recovery = 33.0  / 100.0
    gold_deductions = 0.7
    gold_refining_charges = 5.83
    gold_price = 383
    smelter_charge = 81.66
    penalties = 3.50
    freights = 40.83
    marketing_charges = 2.0
    
    nsr = nsr_calculation(
        copper_grade,
        copper_recovery,
        copper_concentrate_recovery,
        copper_deductions,
        copper_refining_charges,
        copper_price,
        gold_grade,
        gold_concentrate_recovery,
        gold_deductions,
        gold_refining_charges,
        gold_price,
        smelter_charge,
        penalties,
        freights,
        marketing_charges)

    print nsr
