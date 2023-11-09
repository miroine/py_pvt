from .oil_correlation import Oil 
from .gas_correlation import DryGas
from .water_correlation import Water 
import pandas as pd 
import numpy as np 



class OilGas : 

    """
        Creates data required for Oil-Gas-Water black oil tables
        Returns dictionary of results, with index:
  
    """
    def __init__ (self, sto_api: float, sg_gas:float, temp:float, rsb:float =0, pb:float =0, pmin:float = 50, pmax:float =300, step:int =20, unit_system ="metric", \
                  correlation_pb:str = "petrosky", correlation_rs:str ="petrosky", correlation_bo:str ="petrosky", correlation_deno="petrosky", correlation_zfac:str ="DAK",
                  correlation_gasvisco:str ="LBC", salt:float=0, p_sep:float =0, t_sep:float =0, p_ref:float =0
                  ) :
        """ Class to generate black oil tables for the Oil and gas system
            Typical output would be a table or eclipse export files

        Args:
            sto_api (float): Oil Stock tank API
            sg_gas (float): gas surface gravity
            temp (float): Temperature
            rsb (float, optional): RS and standard condition. Defaults to 0.
            pb (float, optional): bubble point pressure. Defaults to 0.
            pmin (float, optional): min pressure. Defaults to 50.
            pmax (float, optional): max pressure. Defaults to 300.
            step (int, optional): pressure step. Defaults to 20.
            unit_system (str, optional): unit system could be field or metric. Defaults to "metric".

        """

        self.sto_api = sto_api

        #correct the gas saturation based on the separator conditions 
        if t_sep * p_sep != 0:
            self.sg_gas = Oil().sg_gas_sep(sg_gas=sg_gas, sto_api=sto_api, t_sep=t_sep, p_sep=p_sep, unit_system=unit_system)
        else:
            self.sg_gas = sg_gas
        self.temp = temp 
        self.rsb = rsb 
        self.pb = pb 
        self.pmin = pmin 
        self.pmax = pmax
        self.step = step 
        self.unit_system = unit_system 
        self.correlation_pb = correlation_pb 
        self.correlation_rs = correlation_rs
        self.correlation_bo = correlation_bo
        self.correlation_deno = correlation_deno
        self.correlation_zfac = correlation_zfac
        self.correlation_gasvisco = correlation_gasvisco 
        self.salt = salt
        if p_ref >0 : 
            self.p_ref = p_ref 
        else: 
            self.p_ref = (pmax - pmin)/2    # just a guess 
     
    
    def get_bo_pvt_table (self):

        if self.unit_system.lower() not in ["metric", "field"]:
            raise ValueError("Unknown unit system")
        
        results = []
        for press in range(int(self.pmin), int(self.pmax+1), self.step):
            
            if self.pb ==0 and self.rsb ==0 : 
                raise ValueError("either pb or rs needs to be provided")
            if self.pb ==0:
                # calculate pb 
                self.pb = Oil().oil_pbubble(rs = self.rsb, sto_api=self.sto_api, sg_gas=self.sg_gas, temp=self.temp, correlation=self.correlation_pb, unit_system=self.unit_system)
            
            # calculate rs 
            rs = Oil().calculate_rs(temp=self.temp, press=press, rsb=self.rsb,pb=self.pb, sg_gas=self.sg_gas, sto_api=self.sto_api, correlation=self.correlation_rs, unit_system=self.unit_system)

            # calculate bo 
            bo = Oil().oil_bo (press=press, temp=self.temp, rs=rs,  sg_gas = self.sg_gas, sto_api=self.sto_api, correlation=self.correlation_bo, unit_system=self.unit_system)

            # calculate density

            deno = Oil().oil_den (press=press, temp=self.temp,  sg_gas = self.sg_gas, sto_api=self.sto_api, rs=self.rsb, correlation=self.correlation_deno, unit_system=self.unit_system)

           # viscosity
            if press > self.pb:
                    # undersaturated

                if self.sto_api < 22:
                    # deald oil 
                    mu = Oil().oil_viscosity(sto_api=self.sto_api, temp= self.temp, press= press, rs=rs, pb=self.pb, oiltype="dead",correlation = "beals", unit_system=self.unit_system)
                else:

                    mu = Oil().oil_viscosity(sto_api=self.sto_api, temp= self.temp, press= press, rs=rs, pb=self.pb, oiltype="undersaturated",correlation = "vasquezbeggs", unit_system=self.unit_system)
            else: 
                # saturated
                #mu = Oil().oil_viscosity(sto_api=self.sto_api, temp= self.temp, press= press, rs=rs, pb=self.pb, oiltype="saturated",correlation = "chewconnally", unit_system=self.unit_system)
                mu = Oil().oil_viscosity(sto_api=self.sto_api, temp= self.temp, press= press, rs=rs, pb=self.pb, oiltype="undersaturated",correlation = "vasquezbeggs", unit_system=self.unit_system) 

            # compressibility 
            co=0
            co = Oil().oil_compressibility(temp=self.temp, sto_api=self.sto_api, press=press, sg_gas=self.sg_gas , rs=rs, pb=self.pb, unit_system=self.unit_system)

            # Gas properties 

            z = DryGas().gas_zfactor(press_i = press , temp_i=self.temp, sg=self.sg_gas, unit_system=self.unit_system, correlation=self.correlation_zfac)

            #bg 

            bg = DryGas().gas_fvf(temp_i= self.temp, press_i=press, sg=self.sg_gas, z_factor=z, unit_system=self.unit_system)

            #gas density

            deng = DryGas().gas_density(temp_i= self.temp, press_i=press, sg=self.sg_gas, z_factor=z, unit_system=self.unit_system)

            #cg 

            cg= DryGas().gas_compressibility(temp_i= self.temp, press_i=press, sg = self.sg_gas, z_factor=z,rhogas=deng, unit_system=self.unit_system)

            #ug 
            ug= 0    
            ug = DryGas().gas_viscosity(temp_i=self.temp, press_i=press, sg=self.sg_gas, z_factor=z, unit_system=self.unit_system, correlation= self.correlation_gasvisco)

            #rv 
            # calculate condensate api 
            condensate_api = Oil().calc_api(deno/1000)
            rv = DryGas().gas_rv(temp_i=self.temp , press_i=press , sto_api = condensate_api, sg_gas= self.sg_gas)
            # water 

            bw = Water().waterFvf(temp=self.temp, press=press ,salt=self.salt, unit_system=self.unit_system)

            #uw 
            uw = Water().waterViscosity(temp=self.temp, press=press, salt= self.salt, unit_system=self.unit_system)

            #dzdp
            dzpdp = DryGas().dzdp(temp_i=self.temp, press_i= press, sg_gas=self.sg_gas, z_factor=z,unit_system=self.unit_system )


            d = {"press": press , "temp": self.temp, "rs": rs , "bo":bo , "deno": deno , "rv": rv,"uo": mu , "co": co, "gas zfactor": z, "bg":bg , "deng": deng, "cg":cg, "ug":ug, "dzdp": dzpdp, "bw":bw, "uw":uw}

            results.append(d)
            #print (d)
        
        return (pd.DataFrame(results))
    

    def excport_to_eclipse_format (self, pvto =True, export_file="pvt.inc"): 
        """ function to export to eclipse format

        Args:
            pvto (bool, optional): _description_. Defaults to True.
        """

        if pvto: 
            print("pvto requested --> include gas")

            df = self.get_bo_pvt_table()

            d = list()

            # build the pressure table 
            pbs =[]
            for rsb in df['rs'] : 
                # calculate pb 
                pbs.append(Oil().oil_pbubble(sto_api=self.sto_api, sg_gas=self.sg_gas, temp= self.temp, rs =rsb, unit_system=self.unit_system))
            
            n= max(pbs)
            if self.pmax > n:
                pbs += [n + (i*(self.pmax - n)/10)  for i in range (1,10)]
            for rsb in df['rs'] : 
                # calculate pb 

                pb = Oil().oil_pbubble(sto_api=self.sto_api, sg_gas=self.sg_gas, temp= self.temp, rs =rsb, unit_system=self.unit_system)

                for p in set(sorted(pbs)):
                    if p >= pb : 
                        bo = Oil().oil_bo(press=p, temp= self.temp, sg_gas=self.sg_gas, pb=pb, sto_api=self.sto_api)
                        mu = Oil().oil_viscosity(sto_api=self.sto_api, temp=self.temp , press=p , pb= pb, rs=rsb, oiltype="undersaturated",correlation = "vasquezbeggs", unit_system=self.unit_system)

                        d.append( {"pressure": p , "pb": pb, "rsb":rsb , "bo":bo, "visco":mu})
                    

            


            df_results = pd.DataFrame(d).sort_values(['rsb', 'pressure'], ascending = [True, True]).drop_duplicates(keep='first')
            
            #export pvto data: 

            with open("pvt.inc", "w") as text_file: 

                text_file.write("DENSITY\n")
                text_file.write("-- OilDens   WaterDens    GasDens \n")
                text_file.write("-- kg/m3       kg/m3       kg/m3 \n")
                if self.unit_system == "metric":
                    text_file.write(f"  {Oil().calc_oil_dens_api(self.sto_api)*1000 :.1f}         {1000 + self.salt*100 :.1f}     {self.sg_gas*1000 :.1f}  /  \n")
                else:
                    text_file.write(f"  {Oil().calc_oil_dens_api(self.sto_api)*1000/16.0185 :.1f}         {1000/16.0185 + self.salt*100/16.0185 :.1f}     {self.sg_gas*1000/16.0185 :.1f}  /  \n")
                text_file.write(" / \n\n\n")
                text_file.write('-- Reservoir temperature \n')
                text_file.write(f"-- {self.temp}")
                text_file.write(" / \n\n\n")
                text_file.write("PVTW\n")
                bw = Water().waterFvf(press = self.p_ref, temp=self.temp, salt=self.salt , unit_system=self.unit_system)
                cw = Water().waterCompressibility(press=self.p_ref, temp=self.temp, salt=self.salt, unit_system=self.unit_system)
                vw = Water().waterViscosity(temp=self.temp, press= self.p_ref, salt=self.salt, unit_system=self.unit_system)
                text_file.write(f"   {self.p_ref:.1f}      {bw:.3f}    {cw:.3f}    {vw:.3f}    / \n")


                text_file.write('------------------------------------------------ \n')
                text_file.write("--SOLUTION  PRESSURE  OIL FVF     OIL \n")
                text_file.write("-- GOR Rs  ,  Press ,   Bo ,   Viscosity  \n")

                if self.unit_system == "metric":
                    text_file.write("-- Sm3/Sm3  ,  bar ,   m3/Sm3 ,   cP  \n")
                else:
                    text_file.write('mscf/stb        psia    rb/stb    cP \n')
                
                text_file.write("PVTO\n")
                text_file.write('------------------------------------------------ \n')
                for rs in df_results.rsb.unique():
                    temp = df_results.loc[df_results.rsb == rs]
                    r = 0
                    # make sure to get the correct unit in rs 
                    if self.unit_system == "field":
                        rs = rs/1000
                    for _, row in temp.iterrows():
                        if r ==0: 
                            text_file.write(f"{rs:.1f}      {row['pressure']:.1f}      {row['bo']:.5f}        {row['visco']:.5f}  \n")
                        else: 
                            if row['pressure'] == np.max(temp['pressure']):
                                text_file.write(f"          {row['pressure']:.1f}      {row['bo']:.5f}        {row['visco']:.5f} / \n")                            
                            else:
                                text_file.write(f"          {row['pressure']:.1f}      {row['bo']:.5f}        {row['visco']:.5f}  \n")
                        r +=1
                text_file.write(" / \n\n\n")

                # work with gas 

                text_file.write("--PRESSURE     RV    GAS FVF     GAS \n")
                text_file.write("-- Press  ,    rv ,   Bg ,   Viscosity  \n")                

                if self.unit_system == "metric":
                    text_file.write("-- bar  ,  sm3/sm3 ,   Rm3/Sm3 ,   cP  \n")
                else:
                    text_file.write('psia        mscf/stb     rb/stb    cP \n')
                
                text_file.write("PVTG\n")
                text_file.write('------------------------------------------------ \n')                
                
                # the gas is based on constant rv 
                r = 0
                for p in df['press'].unique():
                    rv_list = sorted(list(set(df['rv'].unique().tolist())), reverse=True)
                    #print(set(rv_list))
                    temp = df.loc[df.press ==p]

                    text_file.write(f"{float(temp['press']):.1f}    {float(temp['rv']):.7f}    {float(temp['bg']):.5f}        {float(temp['ug']):.4f}   \n")

                    for rv in rv_list:
                        if int(rv * 1E4) < int(float(temp['rv'])*1E4)  and int(rv*1e4)>0:
                            text_file.write(f"            {rv:.7f}    {float(temp['bg']):.5f}        {float(temp['ug']):.4f}   \n")
                    text_file.write(f"            {0:.6f}    {float(temp['bg']):.5f}        {float(temp['ug']):.4f}  / \n")
                
                text_file.write(" / \n\n\n")
                   

            return df_results


        

