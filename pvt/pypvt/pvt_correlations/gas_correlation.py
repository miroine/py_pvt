import numpy as np
from utils.unit_converter import UnitConverter
import math 




class DryGas():
    """
        Dry Gas Object used to get several correlation 
    """

    
    def gas_pseudoprops(self, temp_i, press_i, sg, x_h2s= 0, x_co2= 0,x_n2= 0, unit_system = "metric",correlation ="sutton"):
        """
        Calculate Gas Pseudo-critical and Pseudo-reduced Pressure and Temperature
        * Pseudo-critical properties
            1- For range: 0.57 < sg < 1.68 (Sutton, 1985)
            2- properties based on Guo and Chalambor for H2S < 3% (2005) 
            3- Standing correlation (1981)
            4- Elsharkawy et al. (2000)
            5- Ahmed (1989)

        * Pseudo-reduced properties
            For range: x_h2s (mol%) < 0.738; x_co2 (mol%) < 0.544; 154 < p (psia) < 7026; 40 < temp (°F) < 300 (error 0.97%)
            (Wichert and Aziz, 1972)
        """
        # check the units
        if unit_system.lower() not in ["metric", "field"]:
            raise ValueError("Unknown unit system")
        
        if unit_system.lower() == "metric":
            temp = UnitConverter().convert(temp_i, "c", "f")
        else : 
            temp = temp_i

        if unit_system.lower() == "metric":
            press = UnitConverter().convert(press_i, "bar", "psi")
        else: 
            press = press_i



        


        temp_R = temp + 459.67 # convert to Rankine
        if correlation == "sutton":
            if sg > 0.57 and sg < 1.68 and x_h2s < 0.738 and x_co2 < 0.544 and press > 154 and press < 7026 and temp > 40 and temp < 300:
                # calculate pseudocritical properties (Sutton, valid for 0.57<sg<1.68)
                P_pc = 756.8 - (131.07 * sg) - (3.6 * sg**2)
                T_pc = 169.2 + (349.50 * sg) - (74 * sg**2) # in Rankine 
            elif correlation == "Guo_Chalambor":
                if x_h2s < 0.03 : 
                    P_pc = 709.604 - (58.718 * sg)
                    T_pc = 170.491 + (307.344 * sg)
                else: 
                    raise ValueError("too high H2S density, correlation not applied")
            elif correlation == "Standing":
                P_pc = 706 - (51.7 * sg) - (11.1 * sg**2)
                T_pc = 187 + (330 * sg) - (71.5 * sg**2)
            elif correlation == "Elsharkawy et al.":
                P_pc = 787.06 - (147.34 * sg) - (7.916 * sg**2)
                T_pc = 149.18 + (358.14 * sg) - (66.976 * sg**2)
            elif correlation == "Ahmed":
                P_pc = 678 - 50*(sg - 0.5)- (206.7 * x_n2) + (440 * x_co2) + (606.7 * x_h2s)
                T_pc = 326 + 315.7*(sg -0.5) - (240 * x_n2) + (83.3 * x_co2) + (133.3 * x_h2s)
            elif correlation == "Equinor":
                gasg = (sg - 0.972*x_n2 - 1.5195*x_co2 - 1.1765*x_h2s) / (1-x_n2 - x_co2 - x_h2s)
                P_pc = 677 + 15 * gasg - 37.5 * gasg**2 
                P_pc = (1- x_n2 - x_co2 - x_h2s) * P_pc + (493*x_n2 + 1071*x_co2 + 1306 *x_h2s)
                T_pc = 168 + 325 * gasg - 12.5 * gasg**2 
                T_pc = (1-x_n2 - x_co2 - x_h2s) * T_pc + (227.3 * x_n2**2 + 547.6 * x_co2**2 + 672.4 * x_h2s)
    
            else:
                P_pc, T_pc, P_pr, T_pr = np.nan, np.nan, np.nan, np.nan
  
            # calculate adjustment to pseudocritical properties for sour gas (Wiechert-Aziz, valid for x_co2<0.544 and x_h2s<0.738)
            e = (120 * (((x_h2s + x_co2)**0.9) - ((x_h2s + x_co2)**1.6))) + (15 * (x_h2s**0.5 - x_h2s**4))
  
            T_pc = T_pc - e # corrected T_pc
            P_pc = (P_pc * T_pc) / (T_pc - x_h2s * e * (1-x_h2s))  
            # calculate pseudoreduced properties
            P_pr = press / P_pc
            T_pr = temp_R / T_pc 

            return(P_pc, T_pc, P_pr, T_pr)
        

    def gas_zfactor(self,  temp_i , press_i , sg, x_h2s=0, x_co2=0,x_n2=0,P_pc =None, T_pc=None, correlation ="DAK", unit_system = "metric" ):
        """
        Calculate Gas Compressibility Factor
        1. DAK : Correlation from (Dranchuk and Aboukassem, 1975, valid for range 0.2 < P_pr < 30; 1 < T_pr < 3 )
        2. Wang Correlation (2021)
        3. Hall and Yarboroug z-factor correlation based on SPE phase behaviour monograph eq 3.42 (validation range 1< Tr < 3 and 0.2 < pr < 25-30 )
        4. Brill&Beggs z-factor correlation  Validation range 1.2 < Tr < 2 and 0.2 < pr < 15. 
        """

        # T_pr : calculated pseudoreduced temperature
        # P_pr : calculated pseudoreduced pressure   
        from scipy.optimize import fsolve # non-linear solver
        import numpy as np
        
        
        if unit_system.lower() not in ["metric", "field"]:
            raise ValueError("Unknown unit system")

        if unit_system.lower() == "metric":
            temp = UnitConverter().convert(temp_i, "c", "f")
        else : 
            temp = temp_i

        if unit_system.lower() == "metric":
                press = UnitConverter().convert(press_i, "bar", "psi")
        else:
            press = press_i



        if not P_pc or not T_pc:
            (P_pc, T_pc, P_pr, T_pr) = self.gas_pseudoprops(temp_i = temp_i, press_i = press_i, sg=sg, x_h2s=x_h2s, x_co2=x_co2,x_n2=x_n2, unit_system = unit_system)
        else:
            temp_R = temp + 459.67 # convert to Rankine 
            P_pr = press / P_pc
            T_pr = temp_R / T_pc         


        if correlation == "DAK":
            #print(T_pr, P_pr)
            if T_pr > 1 and T_pr < 3 and P_pr > 0.2 and P_pr < 30:
                a1 = 0.3265; a2 = -1.0700; a3 = -0.5339; a4 = 0.01569; a5 = -0.05165; a6 = 0.5475
                a7 = -0.7361; a8 = 0.1844; a9 = 0.1056; a10 = 0.6134; a11 = 0.7210

                def f(y):
                    rho_pr, z = y
                    c1 = a1 + (a2/T_pr) + (a3/(T_pr**3))+ (a4/(T_pr**4))+ (a5/(T_pr**5))
                    c2 = a6 + (a7/T_pr) + (a8/(T_pr**2))
                    c3 = a9*((a7/T_pr) + (a8/(T_pr**2)))
                    c4 = (a10)*(1+(a11*(rho_pr**2)))*((rho_pr**2)/(T_pr**3))*(np.exp(-a11*(rho_pr**2)))

                    f1 = z + (c3*(rho_pr**5)) - (c2*(rho_pr**2)) - (c1*(rho_pr**1)) - c4 - 1
                    f2 = rho_pr - ((0.27 * P_pr) / (z * T_pr))
                    return[f1, f2]

                _, z_factor = fsolve(f, [1, 1]) # initial guess
            else: 
                raise ValueError("Not Applicable for the Tr and Pr Ranges ")

        elif correlation == "Wang" : 
            if T_pr > 1 and T_pr < 3 and P_pr > 0.2 and P_pr < 30:
                a1 = 256.41675; a2 = 7.18202; a3 = -178.57250; a4 = 182.98704; a5 = -40.74427; a6 = 2.24427
                a7 = 47.44825; a8 = 5.28520; a9 = -0.14914; a10 = 271.50446; a11 = 16.2694; a12= -121.51728
                a13 = 167.71477; a14 = -81.73093; a15= 20.36191; a16= -2.1170; a17=124.6444; a18= -6.74331
                a19 = 0.20897; a20= -0.00314 

                c1 = a1 + a2* (1+ a3* T_pr + a4 * T_pr**2 + a5 * T_pr**3 + a6 * T_pr**4) * P_pr 
                c2 = a7* P_pr**2 + a8 * P_pr**3 + a9 * P_pr**4 
                c3 = a10 + a11 *(1+ a12* T_pr + a13* T_pr**2 + a14*  T_pr**3 + a15 * T_pr**4 + a16 * T_pr**5) * P_pr 
                c4 = a17 * P_pr**2 + a18 * P_pr**3 + a19 * P_pr**4 + a20 * P_pr**5 


                z_factor = (c1 + c2) /(c3+c4)
                
  
        elif correlation == "Equinor":
            a1 = 0.064225133 
            a2 = (0.53530771 * T_pr) - 0.61232032
            a3 = (0.31506237 * T_pr) - 1.0467099 - (0.57832729 / T_pr ** 2)
            a4 =  T_pr 
            a5 = 0.68157001 / T_pr ** 2
            a6 = 0.68446549
            a7 = 0.27 * P_pr 

            epsilon = 1E-6 
    
            rho_pr = 0.27 * P_pr / T_pr  # initial guess 
            rho_pr_old = rho_pr 
            while (rho_pr - rho_pr_old) < epsilon * rho_pr :
                rho_pr_old = rho_pr
                frho = a1 * rho_pr**6 + a2 * rho_pr**3 + a3*rho_pr**2 + a4 *rho_pr + a5 * rho_pr**3 * (1+a6*rho_pr**2) * np.exp(-a6 * rho_pr**2)- a7 
                dfrho = 6 * a1 * rho_pr ** 5 + 3 * a2 * rho_pr ** 2 + 2 * a3 * rho_pr + a4 + a5 * rho_pr ** 2 * (3 + a6 * rho_pr ** 2 * (3 - 2 * a6 * rho_pr ** 2)) * np.exp(-a6 * rho_pr ** 2)
                rho_pr = rho_pr - frho / dfrho
                if (rho_pr - rho_pr_old) >= epsilon * rho_pr:
                    break 
    
            
            z_factor = 0.27 * P_pr / rho_pr / T_pr

        elif correlation == "HallyarBorough": 
               # -------------------------------------------------------------------------
            def fy(y, alpha, pr, t):
                """
                Equation for solving reduced density param. Eq. 3-46 in Monograph

                """

                return (
                    -alpha * pr
                    + (y + y ** 2 + y ** 3 - y ** 4) / (1 - y) ** 3
                    - (14.76 * t - 9.76 * t ** 2 + 4.58 * t ** 3) * y ** 2
                    + (90.7 * t - 242.2 * t ** 2 + 42.4 * t ** 3) * y ** (2.18 + 2.82 * t)
                    )

            t = 1 / T_pr 
            alpha =  0.06125 * t * math.exp(-1.2 * (1 - t) ** 2)
            y0 = 0.001 
            y = fsolve(fy, y0, args=(alpha, P_pr, t))
            pseudo_rho = y[0]
            z_factor = alpha * P_pr / pseudo_rho
        
        elif correlation == "BrillBegs":
            if T_pr >= 1.1 and T_pr <= 2 and P_pr >= 0.2 and P_pr <= 15: 
                a0 = 1.39 * (T_pr - 0.92) ** 0.5 - 0.36 * T_pr - 0.1
                a1 = 9.0 * (T_pr - 1.0)
                a2 = ((0.62 - 0.23 * T_pr) * P_pr+ (0.066 / (T_pr - 0.86) - 0.037) * P_pr ** 2 + 0.32 * P_pr ** 2 / 10 ** a1)
                a3 = 0.132 - 0.32 * math.log(T_pr, 10)
                a4 = 0.3106 - 0.49 * T_pr + 0.1824 * T_pr ** 2
                a5 = 10 ** a4 
                z_factor = a0 + (1.0 - a0) / (math.exp(a2)) + a3 * P_pr ** a5


        else:
            z_factor = np.nan

        return(z_factor) # result z-factor
    

    def gas_density(self, temp_i, press_i, sg, z_factor=None, x_h2s=0, x_co2=0, x_n2=0, unit_system ="metric"):
        """
        Calculate Gas Density
        For range: this is not a correlation, so valid for infinite intervals
        """  
        
        if unit_system.lower() not in ["metric", "field"]:
            raise ValueError("Unknown unit system")
        
        if unit_system.lower() == "metric":
            temp = UnitConverter().convert(temp_i, "c", "f")
        else: 
            temp = temp_i

        if unit_system.lower() == "metric":
            press = UnitConverter().convert(press_i, "bar", "psi")
        else:
            press = press_i

                
        if not z_factor:
            z_factor = self.gas_zfactor( temp_i , press_i , sg, x_h2s=0, x_co2=0,x_n2=0,P_pc =None, T_pc=None, correlation ="DAK", unit_system = unit_system )


        
        R = 10.732 # gas constant in (ft3*psi)/(lb-mol*R) 

        rhogas = (28.97 * sg * press) / (z_factor * R * temp)
        if unit_system.lower() == "metric":
            print("density converted to kg/m3 ")
            rhogas = UnitConverter().convert(rhogas, "lb/ft3", "kg/m3")
        return rhogas 

    def gas_fvf(self, temp_i, press_i, sg, z_factor=None, x_h2s=0, x_co2=0, x_n2=0, unit_system ="metric"):
        """
        Calculate Gas FVF
        For range: this is not a correlation, so valid for infinite intervals
        """
 
        if unit_system.lower() == "metric":
            temp = temp_i + 273.15  # Kalvein
            t_sc = 15.56 + 273.15  # K
            p_sc = 1.013529  # bar
        else:
            temp = temp_i + 459.67  # R
            t_sc = 60 + 459.67  # R
            p_sc = 14.7  # psia
 


        
        
        if not z_factor:
            z_factor = self.gas_zfactor(temp_i , press_i , sg, x_h2s=0, x_co2=0,x_n2=0,P_pc =None, T_pc=None, correlation ="DAK", unit_system = unit_system )
        
        
        bg =  z_factor* (temp/t_sc) * (p_sc/press_i) 
        return(bg)

    def gas_viscosity(self, temp_i, press_i, sg, z_factor=None, rhogas=None, unit_system ="Field", P_pc =None , T_pc=None,x_h2s=0, x_co2=0, x_n2=0, correlation = "LBC"):
        """
        Calculate Gas Viscosity 
        For gas with CO2 and N2 composition
        1. LBC correlation: For range: 100 < temp (°F) < 340; 0.9 < x_CO2 (mol%) < 3.2; x_N2 (mol%) < 4.8 (std 2.7-9.0%) (Lee et al, 1996)
        2. Lucas gas viscosity bsaed on Whitson monograph eq 3.66  Validation range 1 < Tr < 40 and 0 < pr < 100
        """
        #if temp: 
        #    if self.unit_system.lower() == "metric":
        #        temp = UnitConverter().convert(self.temp, "c", "f")
        #if press: 
        #    if self.unit_system.lower() == "metric":
        #        press = UnitConverter().convert(self.press, "bar", "psi")
        
        if unit_system.lower() not in ["metric", "field"]:
            raise ValueError("Unknown unit system")
        if temp_i: 
            if unit_system.lower() == "metric":
                temp = UnitConverter().convert(temp_i, "c", "f")
            else:
                temp = temp_i

        if press_i: 
            if unit_system.lower() == "metric":
                press = UnitConverter().convert(press_i, "bar", "psi")
            else: 
                press = press_i
        


        if not T_pc or P_pc:
            (P_pc, T_pc, P_pr, T_pr) = self.gas_pseudoprops(temp_i = temp_i, press_i = press_i, sg=sg, x_h2s=x_h2s, x_co2=x_co2,x_n2=x_n2, correlation ="sutton",unit_system = unit_system)

        if not z_factor:
            z_factor = self.gas_zfactor(temp_i = temp_i, press_i = press_i, sg=sg, x_h2s= x_h2s, x_co2=x_co2,x_n2=x_n2,P_pc =P_pc, T_pc=T_pc, unit_system = unit_system)

        if not rhogas:
            # make sure the density is on feet 
            rhogas = self.gas_density(temp_i=temp, press_i = press, sg=sg, z_factor=z_factor, unit_system="Field")

        if correlation == "LBC":
            if temp >= 100 and temp <= 340:
                temp = temp + 459.67
                Mg = 28.97 * sg
                rhogas_lee = rhogas / 62.4 # lbm/ft3 converted to gas density unit of Lee et al (g/cm3)
                K = ((0.00094 + 2E-06*Mg)*(temp**1.5)) / (209 + 19*Mg + temp)
                x = 3.5 + (986 / temp) + (0.01 * Mg)
                y = 2.4 - 0.2*x  
                viscogas = K * np.exp(x * (rhogas_lee**y))
        
        elif correlation == "Lucas":
            if T_pr >=1 and T_pr <= 40 and P_pr >=0 and P_pr <= 100:
                Mg = 28.97 * sg
                eps = 9490 * (T_pc / (Mg ** 3 * P_pc ** 4)) ** (1.0 / 6.0)
                ugsc_eps = (0.807 * T_pr ** 0.618- 0.357 * math.exp(-0.449 * T_pr)+ 0.340 * math.exp(-4.058 * T_pr)+ 0.018)
                ugsc = ugsc_eps / eps

                # parameters
                A1 = 1.245e-3 * math.exp(5.1726 * T_pr ** -0.3286) / T_pr
                A2 = A1 * (1.6553 * T_pr - 1.2723)
                A3 = 0.4489 * math.exp(3.0578 * T_pr ** -37.7332) / T_pr
                A4 = 1.7368 * math.exp(2.231 * T_pr ** -7.6351) / T_pr
                A5 = 0.9425 * math.exp(-0.1853 * T_pr ** 0.4489)

                ug_on_ugsc = 1.0 + A1 * P_pr ** 1.3088 / (A2 * P_pr ** A5 + (1 + A3 * P_pr ** A4) ** -1)

                viscogas = ug_on_ugsc * ugsc 
                
        else:
            viscogas = np.nan
        return viscogas

    def gas_compressibility(self, temp_i, press_i, sg, z_factor=None, rhogas=None, unit_system ="metric", P_pc =None , T_pc=None,x_h2s=0, x_co2=0, x_n2=0):
        """
        Calculate Gas Isothermal Compressibility
        For range: unspecified
        (Trube, 1957; Mattar, 1975)
        """


        if unit_system.lower() not in ["metric", "field"]:
            raise ValueError("Unknown unit system")
        if temp_i: 
            if unit_system.lower() == "metric":
                temp = UnitConverter().convert(temp_i, "c", "f")
            else:
                temp = temp_i

        if press_i: 
            if unit_system.lower() == "metric":
                press = UnitConverter().convert(press_i, "bar", "psi")
            else: 
                press = press_i

        


        if not T_pc or P_pc:
            (P_pc, T_pc, P_pr, T_pr) = self.gas_pseudoprops(correlation ="sutton",temp_i = temp_i, press_i = press_i, sg=sg, x_h2s=x_h2s, x_co2=x_co2,x_n2=x_n2, unit_system = unit_system)

        if not z_factor:
            z_factor = self.gas_zfactor(temp_i = temp_i, press_i = press_i, sg=sg, x_h2s= x_h2s, x_co2=x_co2,x_n2=x_n2,P_pc =P_pc, T_pc=T_pc, unit_system = unit_system)

        if not rhogas:
            # density should be on field units
            rhogas = self.gas_density(temp_i=temp, press_i=press, sg=sg, z_factor=z_factor, unit_system="field")


        a1 = 0.3265; a2 = -1.0700; a3 = -0.5339; a4 = 0.01569; a5 = -0.05165; a6 = 0.5475
        a7 = -0.7361; a8 = 0.1844; a9 = 0.1056; a10 = 0.6134; a11 = 0.7210

        rhor = 0.27 * P_pr / (T_pr * z_factor)
        dzdrho = (a1 + (a2 / T_pr) + (a3 / T_pr**3) + (a4 / T_pr**4) + (a5 / T_pr**5) )
        dzdrho = dzdrho + (2 * rhor + (a6 + (a7/T_pr) + (a8/T_pr**2) ))
        dzdrho = dzdrho - (5 * np.power(rhor, 4) * a9 * ((a7 / T_pr) + (a8 / (T_pr * T_pr))))
        dzdrho = dzdrho + (2 * a10 * rhor / (T_pr * T_pr * T_pr)) * (1 + (a11 * rhor * rhor) - (a11 * a11 * np.power(rhor, 4))) * np.exp(-a11 * rhor * rhor)  # 2.23
        cpr = (1 / P_pr) - ((0.27 / (z_factor**2 * T_pr)) * (dzdrho / (1 + (rhor / z_factor) * dzdrho)))  # 2.22
        cg = cpr / P_pc
        
        if unit_system.lower() == "metric":
            print("density converted to 1/bar ")
            cg = UnitConverter().convert(cg, "1/psi", "1/bar")

        return(cg)

    def gas_pseudopressure(self,temp_i, press_i, sg, z_factor=None, rhogas=None, unit_system ="metrics", P_pc =None , T_pc=None,x_h2s=0, x_co2=0, x_n2=0, viscogas=None):
        """ calculation of gas pseudo pressure """

        if unit_system.lower() not in ["metric", "field"]:
            raise ValueError("Unknown unit system")        

        if unit_system.lower() == "metric":
            temp = UnitConverter().convert(temp_i, "c", "f")
        else: 
            temp = temp_i

        if unit_system.lower() == "metric":
            press = UnitConverter().convert(press, "bar", "psi")
        else:
            press = press_i


        if not z_factor:
            z_factor = self.gas_zfactor(temp_i = temp_i, press_i = press_i, sg=sg, x_h2s= x_h2s, x_co2=x_co2,x_n2=x_n2,P_pc =P_pc, T_pc=T_pc, unit_system = unit_system)


        
        if not viscogas:
            # make sure vicosity in field units
            viscogas = self.gas_viscosity(correlation = "LBC", temp_i=temp, press_i=press, sg=sg, z_factor=z_factor, rhogas=rhogas, unit_system ="field", P_pc =P_pc , T_pc=T_pc,x_h2s=x_h2s, x_co2=x_co2, x_n2=x_n2) 
        
        mp = 0
        pold = 0
        xold = 0
        pstep = press / 20
        for n in range(1,20):
            pnew = pold + pstep
            xnew = 2 * pnew / z_factor / viscogas
            mp = mp + (xold + xnew) / 2 * pstep
            pold = pnew
            xold = xnew
        pseudopress = pnew

        return pseudopress 


