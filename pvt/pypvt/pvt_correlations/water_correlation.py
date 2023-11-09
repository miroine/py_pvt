import pandas as pd 
import numpy as np 
from utils.unit_converter import UnitConverter
import math
import warnings 



class Water :
  
    ''' The water properties required for most well test application
    are FVF Bw , Viscosity mu_w and cw'''

    def waterCompressibility(self, press:float , temp:float,salt:float, rsw:float=0, unit_system:str ="metric"):
        """ Calculate water compressibility

        Args:
            press (float): Pressure
            temp (float): Temperature
            rsw(float): gas solubitity
            salt(float) salt percentage
            unit_system (str, optional): _description_. Defaults to "metric".

        Raises:
            ValueError: _description_

        Returns:
            _type_: _description_
        """

        if unit_system.lower() not in ["metric", "field"]:
                raise ValueError("Unknown unit system")

        ucnv = UnitConverter().convert

        if rsw == 0:
             rsw = self.gasSolubility(temp=temp, press=press, unit_system=unit_system)
        # Correlation is in field units
        if unit_system.lower() == "metric":
            press = ucnv(press, "bar", "psi")
            temp = ucnv(temp, "c", "f")
            rsw = ucnv(rsw , "m3/m3" , "ft3/bbl")
        

 
        
        c1 = 3.8546 - 0.000134 * press 
        c2 = -0.01052 + 4.77e-7 * press 
        c3 = 3.9267e-5 - 8.8e-10 * press 

        cw = (c1 + c2 * temp + c3 * temp**2) * 10**-6
        cw = cw * (1 + 0.0089 * rsw)   #Dissolved gas correction
        cw = cw * ((-0.052 + (0.00027 * temp) - (0.00000114 * temp ** 2) + 0.000000001121 * temp ** 3) * salt ** 0.7 + 1)


        if unit_system == "metric":
              print("compressibility in bar")
              cw = ucnv(cw, "1/psi", "1/bar")

        return cw         


    def waterFvf (self, press:float, temp:float, salt:float, unit_system:str="metric" ):
        """        water formation volume could be calculated by Eq 2-127
        Ref Reservoir engineering handbook 3E (Tarek Ahmed)

        Args:
            press (float): Pressure
            temp (float): Temperatire
            salt (float): amount of salt in weight percentage
            unit_system (str optional) : unit system

        """

        if unit_system.lower() not in ["metric", "field"]:
                raise ValueError("Unknown unit system")

        ucnv = UnitConverter().convert

        # Correlation is in field units
        if unit_system.lower() == "metric":
            press = ucnv(press, "bar", "psi")
            temp = ucnv(temp, "c", "f")
        

        a = 0.9911 + (0.0000635 * temp) + (0.00000085 * temp ** 2) 
        b = -0.000001093 - (0.000000003497 * temp) + (0.00000000000457 * temp ** 2)
        c = -0.00000000005 + (6.429e-13 * temp) - (1.43e-15 * temp ** 2)
        bw = a + (b * press) + (c * press ** 2)
        bw = bw * ((0.000000051 * press + (0.00000547 - 0.000000000195 * press) * (temp - 60) + (-0.0000000323 + 0.00000000000085 * press) * (temp - 60) ** 2) * salt + 1)              

        return bw 
    
    def waterViscosity (self,  temp:float, press:float, salt:float, unit_system:str ="metric"):
        """calculate water viscosity bsaed on Brill and Beggs correlaiton 


        Args:
            temp (float): temperature
            press (float): Pressure
            salt (float): salt amount
            unit_system (str, optional): _description_. Defaults to "metric".
        """

        if unit_system.lower() not in ["metric", "field"]:
                raise ValueError("Unknown unit system")

        ucnv = UnitConverter().convert

        # Correlation is in field units
        if unit_system.lower() == "metric":
            temp = ucnv(temp, "c", "f")


        tc = 5 / 9 * (temp - 32)
        tk = tc + 273.15

        a = -7.419242  - (0.29721 * (0.65 - 0.01 * tc)) - (0.1155286 * (0.65 - 0.01 * tc)**2) - (0.008685635 * (0.65 - 0.01 * tc)**3) + (0.001094098 * (0.65 - 0.01 * tc)**4)

        a += 0.00439993 * (0.65 - 0.01 * tc)**5 
        a += 0.002520658 * (0.65 - 0.01 * tc)**6
        a += 0.0005218684 * (0.65 - 0.01 * tc)**7

        psat = 22088 * np.exp((374.136 - tc) * a / tk)

        uw = 0.02414 * 10 ** (247.8 / (tk - 140)) * (1 + (press / 14.504 - psat) * 0.0000010467 * (tk - 305))
        uw = uw * (1 - 0.00187 * (salt ** 0.5) + 0.000218 * (salt ** 2.5) + (temp ** 0.5 - 0.0135 * temp) * (0.00276 * salt - 0.000344 * salt ** 1.5))

        return uw 
    
    def gasSolubility (self, press:float , temp:float,salt:float, unit_system:str ="metric"):
        """calculate gas solubility in water 

        Args:
            press (float): Pressure 
            temp (float): Temperature
            unit_system (str, optional): _description_. Defaults to "metric".
        """

        if unit_system.lower() not in ["metric", "field"]:
                raise ValueError("Unknown unit system")

        ucnv = UnitConverter().convert

        # Correlation is in field units
        if unit_system.lower() == "metric":
            press = ucnv(press, "bar", "psi")
            temp = ucnv(temp, "c", "f")
        
        a = 2.12 + (3.45e-3 * temp) - (3.59e-5 * temp**2)
        b = 0.0107 + (-5.26e-5 * temp) + (1.48e-7 * temp**2)
        c = 8.75e-7 + (3.9e-9 * temp) + (1.02e-11 * temp**2)

        rsw = a + b*press + (c* press**2)

        rsw = rsw * (1 - (0.0753 - 0.000173 * temp) * salt)


        if unit_system =="metric":
             rsw = ucnv(rsw ,"ft3/bbl", "m3/m3" )

        return rsw

    
      





