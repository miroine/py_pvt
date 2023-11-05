import numpy as np
from utils.unit_converter import UnitConverter
import math
import warnings 



class Oil():
    """
        Oil Object used to get several pvt correlation
        most of the PVT equations have been derived from Reservoir Engineering handbook 3 E (Tarek Ahmed) 
    """

    def calc_api(self, calc_dens_sto):
        """
        Calculate Oil API gravity [degAPI]

        Parameters
        ----------
        calc_dens_sto : float
            Stock tank Oil (C5+) density [g/cm3]

        Returns
        -------
        calc_api : float
            Oil API gravity [degAPI]
        """
        
        calc_api = 141.5 * 0.9991 / calc_dens_sto - 131.5
        return calc_api
    
    def calc_oil_dens_api (self, sto_api):
         """ 
         calculate oil density from api 
        
        Parameters
        ----------
        calc_api : float
            Oil API gravity [degAPI]

        Returns
        -------
        calc_dens_sto : float
            Stock tank Oil (C5+) density [g/cm3]
         """

         calc_dens_sto = 141.5/(sto_api+131.5) * 998.9/1000

         return calc_dens_sto
    
    
    def oil_pbubble(self, rs, sto_api, sg_gas, temp, correlation="standing", unit_system="metric"):
        """
            calculate oil bubble pressure based on several correlations 
            Most of the correlation are part of the SPE monograph manual            

        Parameters
        ----------
            rs : float
                single-stage gas-oil solution ratio (metric: sm3/sm3, field: scf/stb)
            sto_api : float
                stock tank oil API gravity
            sg_gas : float
                specific gravity of surface gas
            temp : float
                reservoir temperature (metric: Celcius, field Farenheit)
            correlation : str 
                standing : based on Standing bubble-point correlation 
                            Application range : lighter volatile oils, low res. temperatures               
                vazquezbeggs : based on Vazquez and Beggs correlation 
                glaso  : based on glaso correlation to be applied mainly for NCS wells
                marhoun : based on Marhoun correlation, applicable mainly for middle east
                velardi : based on Velardi correlation 
                petrosky : based on Petrosky and Farshad correlation
                valko_mccain : based on Valko McCain correlation made in 2003 

                
            unit_system : str, optional
                "metric" or "field", by default "metric"

        Returns
        -------
        float
            bubble-point pressure (metric: bar, field: psi)

        Raises
        ------
            ValueError
                Unknown unit system.
        """
        correlation_list = ['standing' , "glaso","velarde","vasquezbeggs", "valkomccain", "marhoun", "petrosky"]
        if correlation.lower() not in correlation_list: 
            raise ValueError(f"correlation needs to be one of this: {correlation_list}")
        
        
        if correlation == "standing":
            """ Standing correlation : SPE phase behavior monograph eq. 3.79

            Raises:
            ValueError: use the correct unit system
            """
        
            
            
            if unit_system.lower() not in ["metric", "field"]:
                raise ValueError("Unknown unit system")

            ucnv = UnitConverter().convert

            # Correlation is in field units
            if unit_system.lower() == "metric":

                rs = ucnv(rs, "m3/m3", "ft3/bbl")
                temp = ucnv(temp, "c", "f")

            x = (rs / sg_gas) ** 0.83 * 10 ** (0.00091 * temp - 0.0125 * sto_api)
            pbub = 18.2 * (x - 1.4)

            if unit_system.lower() == "metric":
                pbub = ucnv(pbub, "psi", "bar")
        
        elif correlation == "glaso":
            """Based on Equation 2.80 from Reservoir Engineering Handbook (Tarek Ahmed)
                    to be mostly applied for the North sea fluid 

            Raises:
                ValueError: use the correct units
            """
        
            ucnv = UnitConverter().convert
            if unit_system.lower() not in ["metric", "field"]:
                raise ValueError("Unknown unit system")
                # Correlation is in field units
            if unit_system.lower() == "metric":

                rs = ucnv(rs, "m3/m3", "ft3/bbl")
                temp = ucnv(temp, "c", "f")

            A = (rs / sg_gas) ** 0.816 * (temp ** 0.172 / sto_api ** 0.989)

            pbub = 10 ** (1.7669 + 1.7447 * math.log(A, 10) - 0.30218 * (math.log(A, 10)) ** 2)

            if unit_system.lower() == "metric":
                pbub = ucnv(pbub, "psi", "bar")
        
        elif correlation == "velarde":

            ucnv = UnitConverter().convert
            if unit_system.lower() not in ["metric", "field"]:
                raise ValueError("Unknown unit system")

                # Correlation is in field units
            if unit_system.lower() == "metric":

                rs = ucnv(rs, "m3/m3", "ft3/bbl")
                temp = ucnv(temp, "c", "f")
            x = 0.013098 * temp ** 0.282372 - 8.2e-6 * sto_api ** 2.176124
            pbub = (1091.47* (rs ** 0.081465 * sg_gas ** -0.161488 * 10 ** x - 0.740152)** 5.354891)
            
            if unit_system.lower() == "metric":
                pbub = ucnv(pbub, "psi", "bar")  

        elif correlation == "vasquezbeggs":
            """based on Reservoir engineering book (Tarek Ahmed 3E)
                    equation 2-79 
                    Calculate Oil Bubble-Point Pressure
                    For range: 20 < Rsb (scf/STB) < 2,070; 0.56 < sg < 1.18; 16 < api < 58; 70 < temp (°F) < 295 (err=0.7%)
                    (Vazquez and Beggs, 1980)

            Raises:
                ValueError: use correct unit_system
            """
            
            
            ucnv = UnitConverter().convert
            if unit_system.lower() not in ["metric", "field"]:
                raise ValueError("Unknown unit system")
            
            # Correlation is in field units
            if unit_system.lower() == "metric":

                rs = ucnv(rs, "m3/m3", "ft3/bbl")
                temp = ucnv(temp, "c", "f")

           
                # c1, c2, c3 coefficient from Vazquez-Beggs

            if sto_api <= 30: 
                c1,c2,c3 =  27.624, 0.914328, 11.172
            else:
                c1,c2,c3  = 56.18, 0.84246, 10.393

            a = -c3 *sto_api / (temp + 459.67)
            #print(a)
            pbub = ((c1 * rs/sg_gas)* 10**a)**c2
            #pbub =   (rs / (c1 * sg_gas * np.exp((c3 * sto_api)/(temp + 459.67))))**(1 / c2) # convert temp to Rankine 

            if unit_system.lower() == "metric":
                pbub = ucnv(pbub, "psi", "bar")

        elif correlation == "valkomccain":
            """ applicable ranges:
                1.8 <= R_s <= 394.7 
                0.725 <= SG_o <= 1.029 
                0.555 <= SG_g <= 1.685 
                15.6 <= T, C <= 172.2 

            Raises:
                ValueError: choose correct unit_system
            """ 
            ucnv = UnitConverter().convert
            if unit_system.lower() not in ["metric", "field"]:
                raise ValueError("Unknown unit system")
            

            if unit_system.lower() == "metric":
                temp = ucnv(temp, "c", "k")
            elif unit_system.lower() == "field":
                temp = ucnv(temp, "f", "k")

            z1 = -4.814074834 + (0.7480913 * np.log(rs)) + (0.1743556 * np.log(rs)**2) + (-0.0206 * np.log(rs)**3)
            z2 = 1.27 + (-0.0449 * sto_api) + (4.36e-4 * sto_api**2) +(-4.76e-6 * sto_api**3)
            z3 = 4.51 + (-10.84 * sg_gas) + (8.39 * sg_gas**2) + (-2.34 * sg_gas**3)
            z4 = -7.2254661 + (0.043155 * temp) - (8.55e-5 * temp**2) +(6.00696e-8 * temp**3)
            z = z1 + z2 + z3 + z4
            lnpb = 2.498 + 0.713 * z + 0.0075 * z ** 2

            pbub = np.exp(lnpb) * 10 # from mega pascal to bar

            if unit_system.lower() == "field":
                pbub = ucnv(pbub, "bar", "psi")
        

        elif correlation == "marhoun": 
            """
                Based on Equation 2.82 from Reservoir Engineering Handbook (Tarek Ahmed)
                to be mostly applied for the middle east fluid

            Raises:
                ValueError: choose correct unit_system
            """
            ucnv = UnitConverter().convert
            if unit_system.lower() not in ["metric", "field"]:
                raise ValueError("Unknown unit system")
            
            # Correlation is in field units
            if unit_system.lower() == "metric":

                rs = ucnv(rs, "m3/m3", "ft3/bbl")
                temp = ucnv(temp, "c", "f")
            
            temp_R = (temp + 459.67)
            sg_oil=141.5/(sto_api+131.5) * 998.9/1000

            a = [5.38088e-3, 0.715082, -1.87784, 3.1437, 1.32657]

            pbub = a[0] * rs**a[1] * sg_gas**a[2] * sg_oil **a[3] * (temp_R)**a[4]

            if unit_system.lower() == "metric":
                pbub = ucnv(pbub, "psi", "bar")

        elif correlation == "petrosky":
            """ 
            Based on Equation 2.83 from Reservoir Engineering Handbook (Tarek Ahmed)
            to be mostly applied for the middle east fluid 

            Raises:
                ValueError: use correct unit_system
            """
            ucnv = UnitConverter().convert
            if unit_system.lower() not in ["metric", "field"]:
                raise ValueError("Unknown unit system")
            
            # Correlation is in field units
            if unit_system.lower() == "metric":

                rs = ucnv(rs, "m3/m3", "ft3/bbl")
                temp = ucnv(temp, "c", "f")
            
            x= (7.916 * (10**-4) *(sto_api)**1.5410) - (4.561 * (10**-5) * (temp)**1.3911)
            pbub = ((112.727 * rs**0.577421) / (sg_gas**0.8439 *(10)**x))-1391.051

            if unit_system.lower() == "metric":
                pbub = ucnv(pbub, "psi", "bar")
                 
               

        else:
                warnings.warn("correlation not yet implemented")
                pbub = np.nan
            
        return pbub
    

    def oil_bob (self,temp, rs, sto_api, sg_gas, unit_system="metric", correlation="standing"):
        """
            oil formation volume factor at bubble-point pressure

            Application range : xxxx

            Input:
                temp    : reservoir temperature (metric: Celcius, field Farenheit)
                rs      : single-stage gas-oil solution ratio (metric: sm3/sm3, field: scf/stb )
                sg_gas  : surface gas specific gravity (dens/dens_air)
                sto_api  : oil api
                unit_system : "metric" or "field"
                correlation: correlation to be used for know only those correlations are implemented
                             - standing : Standing correlation
                             - vasquezbeggs : Vazquez and Beggs correlation
                             - glaso : Glaso's correlation
                             - marhoun : Marhoun correlatio
                             - petrosky : Petrosky and Farshad

            Output:
                bob     : bo @ pbub

        """
        correlation_list = ["standing" , "vasquezbeggs", "glaso", "marhoun", "petrosky"]
        if correlation.lower() not in correlation_list: 
            raise ValueError(f"correlation needs to be one of this: {correlation_list}")
        
        sg_oil = self.calc_oil_dens_api (sto_api)

        if unit_system.lower() not in ["metric", "field"]:
            raise ValueError("Unknown unit system")

        ucnv = UnitConverter().convert

        if correlation == "standing":
            # Correlation is in field units
            if unit_system.lower() == "metric":
                rs = ucnv(rs, "m3/m3", "ft3/bbl")
                temp = ucnv(temp, "c", "f")

            A = rs * (sg_gas / sg_oil) ** 0.5 + 1.25 * temp

            bob = 0.9759 + 12e-5 * A ** 1.2
        
        elif correlation == "vasquezbeggs": 
            """
             equation 2.86 (Reservoir engineering handbook Tarek Ahmed)
             Required inputs:
                Rs should be in scf/stb 
                T is in R 
                sg_gas 
            """
            
            if unit_system.lower() == "metric":
                rs = ucnv(rs, "m3/m3", "ft3/bbl")
                temp = ucnv(temp, "c", "f")


            if sto_api <= 30 : 
                 c1, c2, c3 = 4.677e-4 , 1.751e-5 , -1.811e-8
            else: 
                c1 , c2, c3 = 4.67e-4 , 1.1e-5 , 1.337e-9
            
            bob = 1.0 + (c1 * rs) + ((temp + 460 - 520 ) * (sto_api / sg_gas) * (c2 + (c3 *rs)))
        
        elif correlation == "glaso":
            """ 
             Glaso equation from Reservoir engineering handbook (Tarek Ahmed)
             eq 2-88
             field units
            """
            
            if unit_system.lower() == "metric":
                rs = ucnv(rs, "m3/m3", "ft3/bbl")
                temp = ucnv(temp, "c", "f")
            
            # temperature in Rankin

            temp_R = temp + 459.66

            bob_star = rs * (sg_gas / sg_oil)**0.526 + 0.968 * (temp_R - 460) 
            A = -6.58511 + (2.91329 * np.log10(bob_star)) - (0.27683 * np.log10(bob_star)**2)

            bob = 1.0 + (10 ** A)

        elif correlation == "marhoun" : 
            """
            Marhoun equation made in 1988 from Reservoir engineering handbook (Tarek Ahmed)
            eq 2-90
            (field units)
            """
            if unit_system.lower() == "metric":
                rs = ucnv(rs, "m3/m3", "ft3/bbl")
                temp = ucnv(temp, "c", "f")

            temp_R = temp + 459.66 # temp in Rankin
            F = (rs**0.742390) * (sg_gas**0.323294) * (sg_oil ** -1.20204)

            bob = 0.497069 + (0.862963e-3 * temp_R) + (0.182594e-2 * F) + (0.318099e-5 * F**2)
        
        elif correlation == "petrosky": 
            """ 
                Based on Petrosky and Farshad (1993) from Reservoir engineering handbook (Tarek Ahmed)
                eq 2-92
            """
            if unit_system.lower() == "metric":
                rs = ucnv(rs, "m3/m3", "ft3/bbl")
                temp = ucnv(temp, "c", "f")
            temp_R = temp + 459.66   # temp in Rankin  

            a = (rs**0.3738 * (sg_gas**0.2914 / sg_oil**0.6265) + (0.24626 *(temp_R-460)**0.5371))**3.0936

            bob = 1.0113 + (7.2046e-5 * a)


        else: 
            warnings.warn("correlation not yet implemented")
            bob = np.nan 



        return bob


    def sg_evolved_gas(self,press: float, temp: float, rsb: float, sto_api: float, sg: float, unit_system="metric") -> float:
        """ 
        Returns estimated specific gravity of gas evolved from oil insitu due to depressurization below Pb
        uses McCain & Hill Correlation (1995, SPE 30773)

            press: Pressure 
            temp: Temperature 
            rsb: Oil solution GOR at Pb 
            api: Stock tank oil density (API)
            sg: Specific gravity of separator gas (relative to air)
        """ 

        if unit_system.lower() not in ["metric", "field"]:
            raise ValueError("Unknown unit system")

        ucnv = UnitConverter().convert

            # Correlation is in field units
        if unit_system.lower() == "metric":
                rs = ucnv(rs, "m3/m3", "ft3/bbl")
                press = ucnv(press, "bar", "psi")
                temp = ucnv(temp, "c", "f")

        if (press > 314.7):  # Two different sets from original 1995 paper (not reflected in Correlations book)
            a = [0,-208.0797,22885,-0.000063641,3.38346,-0.000992,-0.000081147,-0.001956,1.081956,0.394035,]
        else:
            a = [0,-214.0887,9971,-0.001303,3.12715,-0.001495,-0.000085243,-0.003667,1.47156,0.714002,]
        one_on_sgr = (a[1] / press+ a[2] / press ** 2+ a[3] * press + a[4] / temp ** 0.5 + a[5] * temp
                + a[6] * rsb + a[7] * sto_api + a[8] / sg + a[9] * sg ** 2)  # 
           
        sg_egas = max(1 / one_on_sgr , sg)

        return sg_egas


    def sg_st_gas(self, press_sp: float, rsp: float, sto_api: float, sg_sp: float, temp_sp: float , unit_system:str ="metric") -> float:
        """ 
        Estimates specific gravity of gas evolving from stock tank
        from oil API and separator gas properties & conditions
        Returns sg_st (Stock Tank gas SG relative to air).
        Correlation reproduced from Valko McCain 2003 paper Eq 4-2

        press_sp: Separator pressure 
        rsp: Separator GOR (separator )
        api: Stock tank oil density (API)
        sg_sp: Separator gas specific gravity relative to air
        temp_sp: Separator temperature 
        """ 

        if unit_system.lower() not in ["metric", "field"]:
            raise ValueError("Unknown unit system")

        ucnv = UnitConverter().convert

            # Correlation is in field units
        if unit_system.lower() == "metric":
                rsp = ucnv(rsp, "m3/m3", "ft3/bbl")
                press_sp = ucnv(press_sp, "bar", "psi")
                temp_sp = ucnv(temp_sp, "c", "f")
        
        var = [np.log(press_sp), np.log(rsp), sto_api, sg_sp, temp_sp]
        C = [
            [-17.275, -0.3354, 3.705, -155.52, 2.085],
            [7.9597, -0.3346, -0.4273, 629.61, -7.097e-2],
            [-1.1013, 0.1956, 1.818e-2, -957.38, 9.859e-4],
            [2.7735e-2, -3.4374e-2, -3.459e-4, 647.57, -6.312e-6],
            [3.2287e-3, 2.08e-3, 2.505e-6, -163.26, 1.4e-8],
            ]
        Zn = [sum([C[i][n] * var[n] ** i for i in range(5)]) for n in range(5)]
        Z = sum(Zn)
        sg_st = (1.219 + 0.198 * Z + 0.0845 * Z ** 2 + 0.03 * Z ** 3 + 0.003 * Z ** 4)

        return sg_st
    

    def rsbub (self, sto_api, temp, press, sg_gas, correlation="standing", unit_system="metric") -> float:
        """
        Rs correlations 

        Args:
            sto_api (float): oil API
            temp (float): Temperature 
            pb (float): bubble point pressure
            sg_g (float): gas gravity
            correlation (str, optional): correlation for the calculation. Defaults to "standing".
            unit_system (str, optional): unit system can be field or metric. Defaults to "metric".

        Returns:
            float: Rs for the oil
        """
        correlation_list = ["standing" , "vasquezbeggs", "glaso", "marhoun", "petrosky"]
        if correlation.lower() not in correlation_list: 
            raise ValueError(f"correlation needs to be one of this: {correlation_list}")
        

        if unit_system.lower() not in ["metric", "field"]:
            raise ValueError("Unknown unit system")

        ucnv = UnitConverter().convert

        if correlation == "standing":
            """
                based on standing correlation from Reservoir engineering handbook (Tarek Ahmed)
                Equation 2-70
            """
            # Correlation is in field units
            if unit_system.lower() == "metric":
                    press = ucnv(press, "bar", "psi")
                    temp = ucnv(temp, "c", "f")
            # temp in F
            a = -0.00091 * temp + 0.0125 * sto_api  # Eq 1.64
            

            rsbub =  sg_gas * (((press) / 18.2 + 1.4) * 10 ** a) ** (1.2048)

        elif correlation == "vasquezbeggs" :
            # Correlation is in field units
            if unit_system.lower() == "metric":
                    press = ucnv(press, "bar", "psi")
                    temp = ucnv(temp, "c", "f")

            if sto_api <= 30: 
                c1,c2,c3 = 0.0362 , 1.0937 , 25.7240 
            else: 
                c1,c2,c3 = 0.0178 , 1.1870 , 23.931 

            rsbub = c1 * sg_gas * (press)**c2 * (np.exp(c3 * (sto_api / (temp+460)))) #equation 1.73 

        elif correlation == "marhoun" :
            # Correlation is in field units
            if unit_system.lower() == "metric":
                    press = ucnv(press, "bar", "psi")
                    temp = ucnv(temp, "c", "f") 

            a0,a1,a2,a3,a4 = 185.843208 , 1.877840, -3.1437, -1.32657, 1.39844 
      
            # get the sg_oil from the api 
            sg_oil=141.5/(sto_api+131.5) * 0.9989

            rsbub = ( a0 * (sg_gas)**a1 * (sg_oil)**a2 *(temp+460)**a3 * press) ** a4 #equation 1.74
        

        elif correlation == "glaso":
            # Correlation is in field units
            if unit_system.lower() == "metric":
                    press = ucnv(press, "bar", "psi")
                    temp = ucnv(temp, "c", "f")

            temp_R = temp + 460 
            a = 2.8869 - (14.1811 - 3.3093 * np.log10(press)) ** 0.5 

            rsbub = sg_gas * (10**a  * sto_api ** 0.989 / (temp ** 0.172))**1.2255  #equation 1.75
        
        elif correlation == "petrosky" : 
            # Correlation is in field units
            if unit_system.lower() == "metric":
                    press = ucnv(press, "bar", "psi")
                    temp = ucnv(temp, "c", "f")

            a = (7.916 * 10**-4 * sto_api ** 1.541) - (4.561 * 10**-5 * temp ** 1.3911)
            rsbub = ( (press/112.727 +12.34) * sg_gas ** 0.8439 * 10 ** a) ** 1.73184 
        else: 
             rsbub = np.nan 

        if unit_system =="metric":
             rsbub = ucnv(rsbub ,"ft3/bbl", "m3/m3" )
        return rsbub    



    def oil_den (self, press: float, temp: float, rs: float = 0, sg_gas: float = 0, sg_sp: float = 0,
        pb: float = 0, sg_oil: float = 0, sto_api: float = 0, bob : float = 0,correlation: str = "petrosky", unit_system:str = "metric") -> float:

        """ Returns live oil density calculated with different correlations

            p: Pressure
            pb: Bubble point pressure . Defaults to 0, and not used for densities below Pb. A valid value is required for density calculations above Pb
            temp: Reservoir Temperature 
            rs: Oil solution gas volume (scf/stb)
            sg_gas: Weighted average specific gravity of surface gas (relative to air).
            sg_sp: Separator gas specific gravity (relative to air). If not known, an alternate nethod to estimate pseudo liquid density of surface gas will be used
            sto_api: Stock tank oil density (deg API). If undefined will calculate from sg_o. If both defined api value will prevail
            correlation
        """

        correlation_list = ["direct_calculation", "vasquezbeggs",  "petrosky"]
        if correlation.lower() not in correlation_list: 
            raise ValueError(f"correlation needs to be one of this: {correlation_list}")
        


        if unit_system.lower() not in ["metric", "field"]:
            raise ValueError("Unknown unit system")

        ucnv = UnitConverter().convert

        
        if sg_oil == 0 :
            #warnings.warn("calculating oil gravity from API ")
            sg_oil = self.calc_oil_dens_api(sto_api=sto_api)
        if bob == 0 : 
            #warnings.warn("No Bo inputs --> calculating it using standing")
            bob = self.oil_bob(temp=temp, rs =rs , sto_api=sto_api, sg_gas=sg_gas, unit_system=unit_system , correlation = correlation)
        if rs == 0: 
            #warnings.warn("No RS inputs --> calculating it using standing") 
            rs = self.rsbub(sto_api=sto_api, temp=temp , press=press, sg_gas=sg_gas, unit_system=unit_system, correlation = correlation)

        if pb == 0 : 
             #print ("pb not providing it --> calculating it")
             pb = self.oil_pbubble(sto_api=sto_api, sg_gas = sg_gas, temp=temp, rs=rs, correlation=correlation, unit_system=unit_system)

        if correlation == "direct_calculation":
            if unit_system.lower() == "metric":
                    rs = ucnv(rs , "m3/m3", "ft3/bbl")       
            rho = ((62.4 * sg_oil) + (0.0136 * rs * sg_gas)) / bob
            
        if press <= pb:

            # Correlation is in field units
            if unit_system.lower() == "metric":
                    temp = ucnv(temp, "c", "f")
                    rs = ucnv(rs , "m3/m3", "ft3/bbl")
                    
            a = (62.4 * sg_oil + 0.0136 * rs * sg_gas)
            b = 0.972 + (0.000147 * ((rs *(sg_gas/sg_oil)**0.5)+ 1.25 * temp)** 1.175)
            rho = a / b

        else :
            if correlation == "vasquezbeggs": 
                # Correlation is in field units
                    if unit_system.lower() == "metric":
                        temp = ucnv(temp, "c", "f")
                        press = ucnv(press, "bar", "psi")
                        pb = ucnv(pb, "bar", "psi")
                        rs = ucnv(rs , "m3/m3", "ft3/bbl")
                    #print ("above bubble point")
                    a = (62.4 * sg_oil + 0.0136 * rs * sg_gas)
                    b = 0.972 + (0.000147 * ((rs *(sg_gas/sg_oil)**0.5)+ 1.25 * temp)** 1.175)
                    rhob = a / b

                    rsb = self.rsbub(sto_api=sto_api, temp=temp , press=pb, sg_gas=sg_gas, unit_system=unit_system)
                    A = 10**-5 * (-1433 + 5*rsb +17.2 * temp -1180*sg_gas + 12.61 *sto_api)
                    rho = rhob * np.exp (A *np.log(press/pb))
 
            elif correlation == "petrosky":

            # Correlation is in field units
                if unit_system.lower() == "metric":
                    temp = ucnv(temp, "c", "f")
                    press = ucnv(press, "bar", "psi")
                    pb = ucnv(pb, "bar", "psi")
                    rs = ucnv(rs , "m3/m3", "ft3/bbl")
                #print ("above bubble point")
                a = (62.4 * sg_oil + 0.0136 * rs * sg_gas)
                b = 0.972 + (0.000147 * ((rs *(sg_gas/sg_oil)**0.5)+ 1.25 * temp)** 1.175)
                rhob = a / b
                rsb = self.rsbub(sto_api=sto_api, temp=temp , press=pb, sg_gas=sg_gas, unit_system=unit_system)  

                A =   4.1646 * (10 **-7) * (rsb**0.69357) * sg_gas**0.1885 * sto_api**0.3273 * temp**0.6729
                
                rho = rhob * np.exp(A * (press**0.4094 - pb**0.4094))
  
        
        if unit_system == "metric":
              #print("density generated in Kg/m3")
              rho = ucnv(rho, "lb/ft3", "kg/m3")
        
        return rho
    
    def oil_compressibility (self,temp : float, sto_api: float, press: float,sg_gas:float,rs: float = 0,pb:float = 0  ,correlation="vasquezbeggs", unit_system="metric"):
        """ return isothermal compressibility of the oil

        Args:
            temp (float): temperature
            sto_api (float): oil API
            press (float): pressure
            rs  (float) : gas solubility at pb
            sg_gas (float): gas gravity
            pb (str, optional): bubble point pressure. Defaults to zero should be input to calculate compresibility below bubble point.
            correlation (str, optional): _description_. Defaults to "vasquezbeggs".
            unit_system (str, optional): _description_. Defaults to "metric".
        """

        correlation_list = ["vasquezbeggs", "petrosky"]
        if correlation.lower() not in correlation_list: 
            raise ValueError(f"correlation needs to be one of this: {correlation_list}")      
        if unit_system.lower() not in ["metric", "field"]:
            raise ValueError("Unknown unit system")

        ucnv = UnitConverter().convert

        if correlation == "vasquezbeggs" : 
            # Correlation is in field units
            if unit_system.lower() == "metric":
                    press = ucnv(press, "bar", "psi")
                    temp = ucnv(temp, "c", "f")
                    rs = ucnv(rs , "m3/m3", "ft3/bbl")
            if rs == 0: 
                 warnings.warn ("no RS present in inputs --> calculating it")
                 rs = self.rsbub(sto_api=sto_api, temp=temp, press=press, sg_gas=sg_gas) 
                 #print(rs)
            co = (-1.433 + 5*rs + (17.2*temp) - (1.180 * sg_gas ) + (12.61*sto_api))/(10**5 *press)
            #print(co)
        elif correlation == "petrosky":
            # Correlation is in field units
            if unit_system.lower() == "metric":
                    press = ucnv(press, "bar", "psi")
                    temp = ucnv(temp, "c", "f")
                    rs = ucnv(rs , "m3/m3", "ft3/bbl") 

            co = 1.705 * 10**-7 * rs**0.69357 * sg_gas**0.1885 * sto_api**0.3272 * temp**0.6729 * press**-0.5906
            if press <= pb: 
                A = -7.573 - (1.45 * np.log(press)) - (0.383* np.log(pb)) + (1.402 * np.log(temp+460)) + (0.256 * np.log(sto_api)) + (0.449 * np.log(rs))
                co = np.exp(A)
        else: 
             warnings.warn("something wrong compressibility not generated, verify inputs")
             co = np.nan
        if unit_system == "metric":
              print("compressibility in bar")
              co = ucnv(co, "1/psi", "1/bar") 
        return co            

    def sg_gas_sep(self, sg_gas:float, sto_api:float, t_sep:float, p_sep:float, unit_system="metric"):
        """
        calculated a separator specific gravity based on test separator condition
        Correlation from Vasquez and Beggs (Eq 2-72: Reservoir engineering handbook: Tarek Ahmed)


        Args:
            sg_gas (float): gas specific gravity
            sto_api (float): oil API
            t_sep (float): separator temperature
            p_sep (float): separator pressure
            unit_system (str, optional): units could be field or metric. Defaults to "metric".
        """

        if unit_system.lower() not in ["metric", "field"]:
            raise ValueError("Unknown unit system")

        ucnv = UnitConverter().convert
        
        #field units to be used
        if unit_system.lower() == "metric":
            p_sep = ucnv(p_sep, "bar", "psi")
            t_sep = ucnv(t_sep, "c", "f")
        
        sg_gass = sg_gas * ( 1+ (5.912 * 10**-5)*sto_api * (t_sep) * np.log10(p_sep /114.7))

        return sg_gass 

    def oil_bo(self, press:float, temp:float ,sg_gas:float, pb:float = 0,rs:float=0, sto_api:float = 0,bob:float= 0 , unit_system:str ="metric", correlation:str ="vasquezbeggs"):
        """ generated undersatured oil fvf factor 

        Args:
            press (float): pressure
            temp (float): temperature
            pb (float, optional): bubble point pressure will be calculated if defaulted to zero. Defaults to 0.
            bob (float, optional): oil fvf for at bubble point pressure, will be calcualted if defaulted. Defaults to 0.
            unit_system (str, optional): unit system could be field or metric. Defaults to "metric".
            correlation (str, optional): different correlations . Defaults to "vazquezbeggs".
        """
        if unit_system.lower() not in ["metric", "field"]:
            raise ValueError("Unknown unit system")

        ucnv = UnitConverter().convert
        if rs == 0 :
            #print ("no rs inputs --> calculating it ")
            #sto_api = self.calc_api(calc_dens_sto=sg_oil)
            rs = self.rsbub(sto_api=sto_api, temp=temp , press=press, sg_gas=sg_gas, unit_system=unit_system, correlation = correlation)
        if pb == 0:
            print ("no bubble point pressure inp")
            pb = self.oil_pbubble(rs= rs, sto_api=sto_api, sg_gas=sg_gas, temp=temp,  unit_system=unit_system, correlation=correlation)


        if press < (pb): 
            warnings.warn("pressure lower than pb")
            bo = self.oil_bob(temp=temp, rs=rs, sto_api=sto_api, sg_gas=sg_gas, unit_system = unit_system, correlation = correlation)
        else:
            
            if rs !=0 : 
                rsb = rs 
            else:
                rsb =  self.rsbub(sto_api=sto_api, temp=temp , press=pb, sg_gas=sg_gas, unit_system=unit_system, correlation = correlation)
            
            if bob ==0 :
                bob = self.oil_bob(temp=temp, rs=rsb, sto_api=sto_api, sg_gas=sg_gas, unit_system = unit_system, correlation = correlation)
            
                    #field units to be used
            if unit_system.lower() == "metric":
                press = ucnv(press, "bar", "psi")
                pb = ucnv(pb, "bar", "psi")
                temp = ucnv(temp, "c", "f")
                rs = ucnv(rs, "m3/m3", "ft3/bbl")
            #a = 4.1646* (10**-7) * (rsb**0.69357) * (sg_gas**0.1885) *(sto_api**0.3272) *(temp**0.6729)

            a = 10**-5 * (-1433 + (5*rsb) + (17.2* temp) - (1180*sg_gas)+ (12.61* sto_api))
            

            bo = bob * np.exp (-a * np.log(press/pb))
            #bo = bob * np.exp(-a * (press**0.4094-pb**0.4094))
        
        return bo


    def calculate_rs (self, sg_gas:float, sto_api:float, press:float, temp:float, pb:float=0, rsb:float =0, unit_system="metric", correlation="standing"):
        
        """calculate Rs for every pressure 

            Args:
            sg_gas (float): gas gravity
            sto_api (float): oil API
            press (float): Pressure
            temp (float): Temperature
            pb (float, optional): bubble point pressure. Defaults to 0.
            rsb (float, optional):  Rsb. Defaults to 0.
            unit_system (str, optional): could be field or metric. Defaults to "metric".
            correlation (str, optional): _description_. Defaults to "standing".
        """
        
        if unit_system.lower() not in ["metric", "field"]:
            raise ValueError("Unknown unit system")

        
        if rsb ==0 and pb == 0 : 
            raise TypeError("Pb or rsb has to be defined and greated than zero")
        
        elif rsb == 0: 
            #print ("rsb not defined")        
            rsb = self.rsbub(sto_api=sto_api, temp=temp, press=pb, sg_gas=sg_gas, correlation=correlation, unit_system=unit_system) 
        
        elif pb == 0: 
            #print("calculating pb")
            sg_oil = self.calc_oil_dens_api (sto_api)
            pb = self.oil_pbubble(temp=temp, rs=rsb, sto_api= sto_api, sg_gas=sg_gas,  correlation = correlation, unit_system= unit_system)
            #print(pb)
        if press >= pb : 
            rs = rsb 
        
        else : 
            rs =  self.rsbub(sto_api=sto_api, temp=temp, press=press, sg_gas=sg_gas, correlation=correlation, unit_system=unit_system) 
        
        return rs
            
            

    
    def calculate_bt (self , temp:float , press:float, sto_api:float , sg_gas:float, rs:float = 0, unit_system:str ="metric", correlation="standing"):
        """ total volume factor calculation 

        Args:
            temp (float): temperature
            press (float): pressure
            sto_api (float): oil API
            sg_gas (float): gas gravity
            rs (float, optional): soluation gas if null than it will be calculation. Defaults to 0.
        """
        if unit_system.lower() not in ["metric", "field"]:
            raise ValueError("Unknown unit system")

        ucnv = UnitConverter().convert

        if rs == 0: 
             print ("no rs provided --> wil be calculated from correlations")
             rs = self.rsbub (sto_api=sto_api, temp=temp, press=press,sg_gas=sg_gas, correlation=correlation, unit_system=unit_system)
        #field units to be used
        if unit_system.lower() == "metric":
            press = ucnv(press, "bar", "psi")
            temp = ucnv(temp, "c", "f")
            rs = ucnv(rs, "m3/m3", "ft3/bbl")
        sg_oil = self.calc_oil_dens_api(sto_api=sto_api)
        
        if correlation == "standing":
            c = 2.9 * 10 ** (-0.00027 * rs)
            A = (np.log10( rs * (temp)**0.5 * sg_oil**c * 1/(sg_gas**0.3)) ) - (10.1 - 96.8 / (6.604 + np.log10(press)))
            bt = 10 ** (-5.223 - (47.4/(-12.22 + A)))
        
        elif correlation == "glaso":
             c = 2.9 * 10 ** (-0.00027 * rs)
             A = (rs *(temp)**0.5 * sg_oil**c/sg_gas**0.3) * press**-1.1089
             bt = 10 ** (0.080135 + 0.47257*np.log10(A) + (0.17351 * np.log10(A)**2))
        elif correlation == "marhoun":
            a0, a1, a2, a3, a4 = 0.644516 , -1.079340, 0.724874, 2.006210, -0.761910
            F = rs**a0 * sg_gas**a1 * sg_oil**a2 * (temp+460)**a3 * press**a4 
            #print(F)
            bt = 0.314693 + (0.106253* 10**-4 *F) + (0.18883 * 10**-10 * F**2)
        
        else: 
             warnings.warn("methond not defined")
             bt = np.nan 
        return bt 
    
    def oil_viscosity (self,sto_api:float, temp:float, oiltype:str="dead", pb:float = 0 ,press:float = 0,rs:float = 0,correlation="beals", unit_system:str ="metric"):
        """ Oil viscosity in cp

        Args:
            sto_api (float): Oil api
            temp (float): Temperature
            oiltype (str, optional): Oil type can be dead, saturated, undersatured. Defaults to "dead".
            correlation (str, optional): used correlation. Defaults to "beggs".
            rs (float, optional): used for saturated correlations
            unit_system (str, optional): units can be field or metric
        """
        if unit_system.lower() not in ["metric", "field"]:
            raise ValueError("Unknown unit system")

        ucnv = UnitConverter().convert

        if oiltype == "dead":
            if correlation == "beals" :
                #field units to be used
                if unit_system.lower() == "metric": 
                    temp = ucnv(temp, "c", "f")

                a = 10 ** (0.43 + 8.33/sto_api)
                mu = (0.32 + (1.8*10**7)/sto_api**4.53) * (360/(temp+200))**a
            
            elif correlation == "beggsrobinson":
                #field units to be used
                if unit_system.lower() == "metric": 
                    temp = ucnv(temp, "c", "f")

                z = 3.0324 - (0.02023 * sto_api)
                y = 10 ** z 
                x = y * (temp)** -1.163

                mu = 10**x -1

            elif correlation == "glaso":
                #field units to be used
                if unit_system.lower() == "metric": 
                    temp = ucnv(temp, "c", "f")

                a = 10.313 * (np.log10(temp)) - 36.447 
                mu = (3.141 * (10**10)) * (temp**-3.444 * np.log10(sto_api)**a)
        
        elif oiltype == "saturated": 
            if rs == 0 : 
                  raise ZeroDivisionError("rs should be defined")
            if correlation == "chewconnally":
                if unit_system.lower() == "metric": 
                    temp = ucnv(temp, "c", "f")
                    rs = ucnv(rs, "m3/m3", "ft3/bbl")
                a = 10.313 * (np.log10(temp)) - 36.447 
                mud = (3.141 * (10**10)) * (temp**-3.444 * np.log10(sto_api)**a)
                a = rs * ((2.2 * 10**-7 * rs) - 7.4 * 10**-4)
                
                c = 8.62 * 10**-5 * rs 
                d = 1.1 * 10**-3 * rs 
                e = 3.74 * 10**-3 * rs
                b = 0.68/(10**c) + 0.25 /(10**d) + 0.062 / (10**e)

                mu = 10**a * mud**b 
            elif correlation == "beggsrobinson":
                if unit_system.lower() == "metric": 
                    temp = ucnv(temp, "c", "f")
                    rs = ucnv(rs, "m3/m3", "ft3/bbl")
                z = 3.0324 - (0.02023 * sto_api)

                y = 10 ** z 
                x = y * (temp)** -1.163
                mud = 10**x -1  # deadoil viscosity 

                a = 10.715 * (rs + 100)**-0.515
                b = 5.44 * (rs + 150)**-0.338 

                mu = a * (mud)**b 
        
        elif oiltype == "undersaturated": 
            if press * rs  * pb == 0: 
                  raise ZeroDivisionError("undersaturated oil --> pressure and rs should be defined")
            
            if correlation == "vasquezbeggs": 
                if unit_system.lower() == "metric": 
                    temp = ucnv(temp, "c", "f")
                    rs = ucnv(rs, "m3/m3", "ft3/bbl")  

                z = 3.0324 - (0.02023 * sto_api)
                y = 10 ** z 
                x = y * (temp)** -1.163
                mud = 10**x -1  # deadoil viscosity 

                a = 10.715 * (rs + 100)**-0.515
                b = 5.44 * (rs + 150)**-0.338 

                mub = a * (mud)**b   # viscosity for saturated oil 

                a = (3.9 * 10**-5 *press) -5
                m = 2.6 * press**1.187 * 10**a 
                 

                mu = mub * (press /pb)**m 

        else: 
             warnings.warn ("viscosity not calculated --> verify inputs")
             mu = np.nan 
        return mu 
                    
             
                 
             

                                  

             











