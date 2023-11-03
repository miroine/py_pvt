#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
    PVT functions module
    ~~~~~~~~~~~~~~~~~~~
    This module contains the functions to calculate PVT derivated output such as API, Bo, Bg, GOR etc.
    'Calc' refers to physical relationship
    'Corr' refers to empirical correlation from litterature
"""
import sys
import os
import math
import warnings

from scipy.optimize import fsolve
from pvt_tools import constants
from pvt_tools import unit_converter


def calc_mw_c7_plus(
    m_n2,
    m_co2,
    m_h2s,
    m_c1,
    m_c2,
    m_c3,
    m_ic4,
    m_nc4,
    m_ic5,
    m_nc5,
    m_c6,
    m_c7_plus,
    m_mw_res_fluid,
):
    """
    Calculate Molecular Weight of the C7+ components (g/mole) based molar composition and reservoir fluid molecular weight measured in Lab.

    Parameters
    ----------
    m_n2 : float
        measured N2 in Lab (mole%)
    m_co2 : float
        measured CO2 in Lab (mole%)
    m_h2s : float
        measured H2S in Lab (mole%)
    m_c1 : float
        measured C1 in Lab (mole%)
    m_c2 : float
        measured C2 in Lab (mole%)
    m_c3 : float
        measured C3 in Lab (mole%)
    m_ic4 : float
        measured iC4 in Lab (mole%)
    m_nc4 : float
        measured nC4 in Lab (mole%)
    m_ic5 : float
        measured iC5 in Lab (mole%)
    m_nc5 : float
        measured nC5 in Lab (mole%)
    m_c6 : float
        measured C6 in Lab (mole%)
    m_c7_plus : float
        measured C7+ in Lab (mole%)
    m_mw_res_fluid : float
        measured Reservoir Fluid moleuclar weight (g/mole)

    Returns
    -------
    calc_mw_c7_plus : float
        Calculated C7+ molecular weight (g/mole) deuced from Lab measurments
    """
    calc_mw_c7_plus = (
        m_mw_res_fluid / 0.01
        - m_n2 * constants.N2.MOLECULAR_WEIGHT
        - m_co2 * constants.CO2.MOLECULAR_WEIGHT
        - m_h2s * constants.H2S.MOLECULAR_WEIGHT
        - m_c1 * constants.C1.MOLECULAR_WEIGHT
        - m_c2 * constants.C2.MOLECULAR_WEIGHT
        - m_c3 * constants.C3.MOLECULAR_WEIGHT
        - m_ic4 * constants.C4.MOLECULAR_WEIGHT
        - m_nc4 * constants.C4.MOLECULAR_WEIGHT
        - m_ic5 * constants.C5.MOLECULAR_WEIGHT
        - m_nc5 * constants.C5.MOLECULAR_WEIGHT
        - m_c6 * constants.C6.MOLECULAR_WEIGHT
    ) / m_c7_plus
    return calc_mw_c7_plus


def calc_mw_res_fluid(
    m_n2,
    m_co2,
    m_h2s,
    m_c1,
    m_c2,
    m_c3,
    m_ic4,
    m_nc4,
    m_ic5,
    m_nc5,
    m_c6,
    m_c7_plus,
    m_mw_c7_plus,
):
    """
    Calculate Reservoir Fluid Molecular Weight (g/mole) based on molar composition measured in Lab.

    Parameters
    ----------
    m_n2 : float
        measured N2 in Lab (mole%)
    m_co2 : float
        measured CO2 in Lab (mole%)
    m_h2s : float
        measured H2S in Lab (mole%)
    m_c1 : float
        measured C1 in Lab (mole%)
    m_c2 : float
        measured C2 in Lab (mole%)
    m_c3 : float
        measured C3 in Lab (mole%)
    m_ic4 : float
        measured iC4 in Lab (mole%)
    m_nc4 : float
        measured nC4 in Lab (mole%)
    m_ic5 : float
        measured iC5 in Lab (mole%)
    m_nc5 : float
        measured nC5 in Lab (mole%)
    m_c6 : float
        measured C6 in Lab (mole%)
    m_c7_plus : float
        measured C7+ in Lab (mole%)
    m_mw_c7_plus : float
        measured C7+ molecular weight in Lab (g/mole)

    Returns
    -------
    calc_mw_res_fluid : float
        Calculated Reservoir Fluid molecular weight (g/mole) deuced from Lab measurements
    """
    calc_mw_res_fluid = 0.01 * (
        m_n2 * constants.N2.MOLECULAR_WEIGHT
        + m_co2 * constants.CO2.MOLECULAR_WEIGHT
        + m_h2s * constants.H2S.MOLECULAR_WEIGHT
        + m_c1 * constants.C1.MOLECULAR_WEIGHT
        + m_c2 * constants.C2.MOLECULAR_WEIGHT
        + m_c3 * constants.C3.MOLECULAR_WEIGHT
        + m_ic4 * constants.C4.MOLECULAR_WEIGHT
        + m_nc4 * constants.C4.MOLECULAR_WEIGHT
        + m_ic5 * constants.C5.MOLECULAR_WEIGHT
        + m_nc5 * constants.C5.MOLECULAR_WEIGHT
        + m_c6 * constants.C6.MOLECULAR_WEIGHT
        + m_c7_plus * m_mw_c7_plus
    )
    return calc_mw_res_fluid


def calc_mw_rich_gas(m_n2, m_co2, m_h2s, m_c1, m_c2, m_c3, m_ic4, m_nc4):
    """
    Calculate Rich Gas (C4-) Molecular Weight [g/mole]

    Parameters
    ----------
    m_n2 : float
        measured N2 in Lab (mole%)
    m_co2 : float
        measured CO2 in Lab (mole%)
    m_h2s : float
        measured H2S in Lab (mole%)
    m_c1 : float
        measured C1 in Lab (mole%)
    m_c2 : float
        measured C2 in Lab (mole%)
    m_c3 : float
        measured C3 in Lab (mole%)
    m_ic4 : float
        measured iC4 in Lab (mole%)
    m_nc4 : float
        measured nC4 in Lab (mole%)

    Returns
    -------
    calc_mw_rich_gas : float
        Rich Gas molecular weight [g/mole]
    """
    calc_mw_rich_gas = (
        m_n2 * constants.N2.MOLECULAR_WEIGHT
        + m_co2 * constants.CO2.MOLECULAR_WEIGHT
        + m_h2s * constants.H2S.MOLECULAR_WEIGHT
        + m_c1 * constants.C1.MOLECULAR_WEIGHT
        + m_c2 * constants.C2.MOLECULAR_WEIGHT
        + m_c3 * constants.C3.MOLECULAR_WEIGHT
        + m_ic4 * constants.C4.MOLECULAR_WEIGHT
        + m_nc4 * constants.C4.MOLECULAR_WEIGHT
    ) / (m_n2 + m_co2 + m_h2s + m_c1 + m_c2 + m_c3 + m_ic4 + m_nc4)
    return calc_mw_rich_gas


def calc_vol_rich_gas(m_n2, m_co2, m_h2s, m_c1, m_c2, m_c3, m_ic4, m_nc4):
    """
    Calculate Rich Gas (C4-) Volume [Scm3]

    Parameters
    ----------
    m_n2 : float
        measured N2 in Lab (mole%)
    m_co2 : float
        measured CO2 in Lab (mole%)
    m_h2s : float
        measured H2S in Lab (mole%)
    m_c1 : float
        measured C1 in Lab (mole%)
    m_c2 : float
        measured C2 in Lab (mole%)
    m_c3 : float
        measured C3 in Lab (mole%)
    m_ic4 : float
        measured iC4 in Lab (mole%)
    m_nc4 : float
        measured nC4 in Lab (mole%)

    Returns
    -------
    calc_vol_rich_gas : float
        Rich Gas volume [Scm3]
    """

    calc_vol_rich_gas = (
        m_n2 + m_co2 + m_h2s + m_c1 + m_c2 + m_c3 + m_ic4 + m_nc4
    ) * constants.IDEAL_GAS.MOLAR_VOLUME
    return calc_vol_rich_gas


def calc_vol_sto(m_ic5, m_nc5, m_c6, m_c7_plus, m_mw_c7_plus, m_sg_c7_plus):
    """
    Calculate mole per unit of volume [mol/cm3] based on molar composition measured in Lab.

    Parameters
    ----------
    m_ic5 : float
        measured iC5 in Lab (mole%)
    m_nc5 : float
        measured nC5 in Lab (mole%)
    m_c6 : float
        measured C6 in Lab (mole%)
    m_c7_plus : float
        measured C7+ in Lab (mole%)
    m_mw_c7_plus : float
        measured C7+ molecular weight in Lab [g/mole]
    m_sg_c7_plus : float
        measured C7+ standard gravity in Lab [g/mole]

    Returns
    -------
    calc_vol_sto : float
        Calculated mole per unit of volume [mol/cm3] deduced from Lab measurments
    """
    calc_vol_sto = (
        m_ic5 * constants.C5.MOLECULAR_WEIGHT / (constants.IC5.STD_LIQ_DEN * 1000)
        + m_nc5 * constants.C5.MOLECULAR_WEIGHT / (constants.NC5.STD_LIQ_DEN * 1000)
        + m_c6 * constants.C6.MOLECULAR_WEIGHT / (constants.C6.STD_LIQ_DEN * 1000)
        + m_c7_plus * m_mw_c7_plus / (m_sg_c7_plus * 1000)
    )
    return calc_vol_sto


def calc_mass_rich_gas(m_n2, m_co2, m_h2s, m_c1, m_c2, m_c3, m_ic4, m_nc4):
    """
    Calculate Rich Gas (C4-) Mass [g]

    Parameters
    ----------
    m_n2 : float
        measured N2 in Lab (mole%)
    m_co2 : float
        measured CO2 in Lab (mole%)
    m_h2s : float
        measured H2S in Lab (mole%)
    m_c1 : float
        measured C1 in Lab (mole%)
    m_c2 : float
        measured C2 in Lab (mole%)
    m_c3 : float
        measured C3 in Lab (mole%)
    m_ic4 : float
        measured iC4 in Lab (mole%)
    m_nc4 : float
        measured nC4 in Lab (mole%)

    Returns
    -------
    calc_mass_rich_gas : float
        Rich Gas mass [g]
    """
    calc_mass_rich_gas = (
        m_n2 * constants.N2.MOLECULAR_WEIGHT
        + m_co2 * constants.CO2.MOLECULAR_WEIGHT
        + m_h2s * constants.H2S.MOLECULAR_WEIGHT
        + m_c1 * constants.C1.MOLECULAR_WEIGHT
        + m_c2 * constants.C2.MOLECULAR_WEIGHT
        + m_c3 * constants.C3.MOLECULAR_WEIGHT
        + m_ic4 * constants.C4.MOLECULAR_WEIGHT
        + m_nc4 * constants.C4.MOLECULAR_WEIGHT
    )
    return calc_mass_rich_gas


def calc_mass_sto(m_ic5, m_nc5, m_c6, m_c7_plus, m_mw_c7_plus):
    """
    Calculate Stock tank Oil (C5+) Mass [g]

    Parameters
    ----------
    m_ic5 : float
        measured iC5 in Lab (mole%)
    m_nc5 : float
        measured nC5 in Lab (mole%)
    m_c6 : float
        measured C6 in Lab (mole%)
    m_c7_plus : float
        measured C7+ in Lab (mole%)
    m_mw_c7_plus : float
        measured C7+ molecular weight in Lab [g/mole]
    Returns
    -------
    calc_mass_sto : float
        Rich Gas mass [g]
    """

    calc_mass_sto = (
        m_ic5 * constants.C5.MOLECULAR_WEIGHT
        + m_nc5 * constants.C5.MOLECULAR_WEIGHT
        + m_c6 * constants.C6.MOLECULAR_WEIGHT
        + m_c7_plus * m_mw_c7_plus
    )
    return calc_mass_sto


def calc_dens_rich_gas(calc_mass_rich_gas, calc_vol_rich_gas):
    """
    Calculate Rich Gas (C4-) density [g/cm3]

    Parameters
    ----------
    calc_mass_rich_gas : float
        calculated mass of rich gas from measurment in Lab [g]
    calc_vol_rich_gas : float
        calculated volume of rich gas from measurment in Lab [Scm3]
    Returns
    -------
    calc_dens_rich_gas : float
        Rich Gas density [g/cm3]
    """
    calc_dens_rich_gas = calc_mass_rich_gas / calc_vol_rich_gas / 1000
    return calc_dens_rich_gas


def calc_dens_sto(calc_mass_sto, calc_vol_sto):
    """
    Calculate Stock tank Oil (C5+) density [g/cm3]

    Parameters
    ----------
    calc_mass_sto : float
        Stock tank Oil (C5+) Mass [g]
    calc_vol_sto : float
        mole per unit of volume (mole/cm3) based on molar composition measured in Lab.
    Returns
    -------
    calc_dens_rich_gas : float
        Rich Gas density [g/cm3]
    """
    calc_dens_sto = calc_mass_sto / (calc_vol_sto * 1000)
    return calc_dens_sto


def calc_sg_rich_gas(calc_mw_rich_gas):
    """
    Calculate Rich Gas (C4-) Standard Gravity

    Parameters
    ----------
    calc_mw_rich_gas : float
        Rich Gas molecular weight [g/mole]
    Returns
    -------
    calc_sg_rich_gas : float
        Rich Gas standard gravity [g/cm3]
    """
    calc_sg_rich_gas = calc_mw_rich_gas / constants.AIR.MOLECULAR_WEIGHT
    return calc_sg_rich_gas


def calc_sst_gor(calc_vol_rich_gas, calc_vol_sto):
    """
    Calculate single stage Gas oil Ratio (Sm3/m3) from Lab measurments.

    Parameters
    ----------
    calc_vol_rich_gas : float
        Calculated Rich Gas (C4-) Volume [Scm3] from Lab measurments
    calc_vol_sto : float
        Calculated mole per unit of volume (mol/cm3) deduced from Lab measurments
    Returns
    -------
    calc_sst_gor : float
        Single stage Gas oil Ratio (Sm3/m3)
    """

    calc_sst_gor = calc_vol_rich_gas / calc_vol_sto
    return calc_sst_gor


def calc_sst_cgr(calc_sst_gor):
    """
    Calculate single stage condensate to gas Ratio (Sm3/m3) from Lab measurements.

    Parameters
    ----------
    calc_sst_gor : float
        Calculated Single stage Gas oil Ratio (Sm3/m3)

    Returns
    -------
    calc_sst_cgr : float
        Single stage Condensate to gas Ratio (Sm3/m3)
    """
    calc_sst_cgr = 1 / calc_sst_gor
    return calc_sst_cgr


def calc_ratio_mass_rich_gas(calc_mass_rich_gas, calc_mass_sto):
    """
    Calculate a ratio of the mass of rich gas vs the mass of stock tank oil to be used in bo and bg calculations.

    Parameters
    ----------
    calc_mass_rich_gas : float
        Calculated mass of rich gas (C4-) (g)
    calc_mass_sto : float
        Calculated mass of stock tank oil (C5+) (g)

    Returns
    -------
    calc_ratio_mass_rich_gas : float
        Mass ratio of rich gas
    """
    calc_ratio_mass_rich_gas = (calc_mass_rich_gas + calc_mass_sto) / calc_mass_rich_gas
    return calc_ratio_mass_rich_gas


def calc_ratio_mass_sto(calc_mass_rich_gas, calc_mass_sto):
    """
    Calculate a ratio of the mass of stock tank oil vs the mass of rich gas to be used in bo and bg calculations.

    Parameters
    ----------
    calc_mass_rich_gas : float
        Calculated mass of rich gas (C4-) (g)
    calc_mass_sto : float
        Calculated mass of stock tank oil (C5+) (g)

    Returns
    -------
    calc_ratio_mass_sto : float
        Mass ratio of stock tank oil
    """
    calc_ratio_mass_sto = (calc_mass_rich_gas + calc_mass_sto) / calc_mass_sto
    return calc_ratio_mass_sto


def calc_bo(calc_dens_sto, calc_ratio_mass_sto, density_pres):
    """
    Calculate Formation Volume Factor of Oil (bo) [Sm3/m3]
    Parameters
    ----------
    calc_dens_sto : float
        Calculate Stock tank Oil (C5+) density [g/cm3]
    calc_ratio_mass_sto : float
        Calculate a ratio of the mass of stock tank oil vs the mass of rich gas
    density_pres : float
        Density at reservoir pressure [bar]

    Returns
    -------
    calc_bo : float
        Calculated Bo
    """
    calc_bo = calc_dens_sto * calc_ratio_mass_sto / density_pres
    return calc_bo


def calc_bg(calc_dens_rich_gas, calc_ratio_mass_rich_gas, density_pres):
    """
    Calculate Formation Volume Factor of Gas (bg) [Sm3/m3]

    Parameters
    ----------
    calc_dens_rich_gas : float
        Calculated density of Rich Gas [g/cm3]
    calc_ratio_mass_rich_gas : float
        Calculated ratio of the mass of rich gas vs the mass of stock tank oil
    density_pres : float
        Density at reservoir pressure [bar]

    Returns
    -------
    calc_bg : float
        Calculated Bg [Sm3/m3]
    """
    calc_bg = calc_dens_rich_gas * calc_ratio_mass_rich_gas / density_pres
    return calc_bg


def calc_api(calc_dens_sto):
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


# ==================================
# PVT Correlations from bibliography
# ===================================

# Gas z-factor correlations


def gas_zfac_hallyarborough(
    pres, temp, sg, unit_system="metric", print_warning="always"
):
    """
    Hall and Yarboroug z-factor correlation
    Reference: SPE phase behavior monograph eq. 3.42
    Validation range 1 < Tr < 3 and 0.2 < pr < 25-30.

    Parameters
    ----------
    pres : float
        Pressure (bar/psi)
    temp : float
        Temperature (C/F)
    sg : float
        surface gas Specific Gravity (air=1)
    unit_system : str, optional
        "metric" or "field", by default "metric"
    print_warning : ['always', 'ignore', default', 'error', 'module', 'once'], optional
        'default' — display a warning the first time it is encountered
        'error' — turn the warning into an exception
        'ignore' — ignore the warning
        'always' — always show the warning, even if it was displayed before
        'module' — show the warning once per module
        'once' — show the warning only once, throughout the programe
        By default it is set to always.

    Returns
    -------
    flaot
        z_fac : gas z-factor

    Raises
    ------
    ValueError
        Unknown unit system.
    """

    warnings.simplefilter(print_warning)

    if unit_system.lower() not in ["metric", "field"]:
        raise ValueError("Unknown unit system")

    ucnv = unit_converter.UnitConverter()

    # Correlation is in field units
    if unit_system.lower() == "metric":
        temp = ucnv.convert(temp, "c", "f")
        pres = ucnv.convert(pres, "bar", "psi")
    # Sutton correlation for pseudo critical temperature (R) and pseudo critical pressure
    t_pc = 169.2 + 349.5 * sg - 74.0 * sg ** 2
    p_pc = 756.8 - 131.0 * sg - 3.6 * sg ** 2
    # pseudo reduced temperature and pressures
    t_pr = (temp + 460) / t_pc
    p_pr = pres / p_pc

    if t_pr > 3.0 or t_pr < 1.0:
        warnings.warn(
            message="WARNING: Hall and Yarboroug z-factor correlation outside range, Tr = %6.3f range is [1.0. 3]"
            % (t_pr)
        )
        # print("WARNING: Hall and Yarboroug z-factor correlation outside range")
        # print("Tr = %6.3f range is [1.0. 3]" % (t_pr))

    if p_pr > 30 or p_pr < 0.2:
        warnings.warn(
            message="WARNING: Hall and Yarboroug z-factor correlation outside range, Pr = %6.3f range is [0.2. 30]"
            % (p_pr)
        )
        # print("WARNING: Hall and Yarboroug z-factor correlation outside range")
        # print("Pr = %6.3f range is [0.2. 30]" % (p_pr))

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

    t = 1 / t_pr
    alpha = 0.06125 * t * math.exp(-1.2 * (1 - t) ** 2)
    y0 = 0.001
    y = fsolve(fy, y0, args=(alpha, p_pr, t))

    zfac = alpha * p_pr / y[0]

    return zfac


def gas_zfac_brillbeggs(pres, temp, sg, unit_system="metric", print_warning="always"):
    """
    Brill & Beggs z-factor correlation.
    Validation range 1.2 < Tr < 2 and 0.2 < pr < 15.

    Parameters
    ----------
    pres : float
        Pressure (bar/psi)
    temp : float
        Temperature (C/F)
    sg : float
        surface gas Specific Gravity (air=1)
    unit_system : str, optional
        "metric" or "field", by default "metric"
    print_warning : bool, optional
        If true prints warnings, by default True

    Returns
    -------
    float
        z_fac : gas z-factor

    Raises
    ------
    ValueError
        Unknown unit system.
    """

    warnings.simplefilter(print_warning)

    if unit_system.lower() not in ["metric", "field"]:
        raise ValueError("Unknown unit system")

    ucnv = unit_converter.UnitConverter().convert

    # Correlation is in field units
    if unit_system.lower() == "metric":
        temp = ucnv(temp, "c", "f")
        pres = ucnv(pres, "bar", "psi")

    # Sutton correlation for pseudo critical temperature (R) and
    # pseudo critical pressure
    t_pc = 169.2 + 349.5 * sg - 74.0 * sg ** 2
    p_pc = 756.8 - 131.0 * sg - 3.6 * sg ** 2

    # pseudo reduced temperature and pressures
    t_pr = (temp + 460) / t_pc
    p_pr = pres / p_pc

    if t_pr > 2.0 or t_pr < 1.5:
        warnings.warn(
            "WARNING: Brill and Beggs z-factor correlation outside range Tr = %6.3f range is [1.5, 2.0]"
            % (t_pr)
        )

    if print_warning and (p_pr > 15 or p_pr < 0.2):
        warnings.warn(
            "WARNING: Brill and Beggs z-factor correlation outside range Pr = %6.3f range is [0.2. 15]"
            % (p_pr)
        )

    a = 1.39 * (t_pr - 0.92) ** 0.5 - 0.36 * t_pr - 0.1
    e = 9.0 * (t_pr - 1.0)
    b = (
        (0.62 - 0.23 * t_pr) * p_pr
        + (0.066 / (t_pr - 0.86) - 0.037) * p_pr ** 2
        + 0.32 * p_pr ** 2 / 10 ** e
    )
    c = 0.132 - 0.32 * math.log(t_pr, 10)
    f = 0.3106 - 0.49 * t_pr + 0.1824 * t_pr ** 2
    d = 10 ** f

    zfac = a + (1.0 - a) / (math.exp(b)) + c * p_pr ** d

    return zfac

    # Gas - Bg


def gas_bg(pres, temp, z_fac, unit_system="metric"):
    """Calculates Bg based on a given z-factor
    If you dont have the gas z-factor: use the Hall Yarborough correlation

    Parameters
    ----------
    pres : float
        Pressure (bar/psi)
    temp : float
        Temperature (C/F)
    z_fac : float
        z-factor gas
    unit_system : str, optional
        "metric" or "field", by default "metric"

    Returns
    -------
    float
        bg
    """

    # Define standard condition pressure and temperature
    if unit_system == "metric":
        t = temp + 273.15  # K
        t_sc = 15.56 + 273.15  # K
        p_sc = 1.013529  # bar

    else:
        t = temp + 459.67  # R
        t_sc = 60 + 459.67  # R
        p_sc = 14.7  # psia

    bg = (t / t_sc) * (p_sc / pres) * z_fac

    return bg


# Gas viscosities


def gas_visc_lucas(
    pres, temp, mw_gas, z_fac=None, unit_system="metric", print_warning="always"
):
    """
    Lucas gas viscosity correlation
    Reference: Whitson Monograph eq. 3.66
    Validation range 1 < Tr < 40 and 0 < pr < 100

    ###TODO
    Tc and Pc is calculated from Sutton correlation -
    Should maybe use pc, tc as input?

    Parameters
    ----------
    pres : float
        Pressure (bar/psi)
    temp : float
        Temperature (C/F)
    mw_gas : float
        surface gas molecular weight
    z_fac : float, optional
        Gas z-factor (if None then calculated by correlation), by default None
    unit_system : str, optional
        "metric" or "field", by default "metric"
    print_warning : bool, optional
        If true prints warnings, by default False

    Returns
    -------
    float
        gas viscosity (cp)

    Raises
    ------
    ValueError
        Unknown unit system
    """

    warnings.simplefilter(print_warning)

    if unit_system.lower() not in ["metric", "field"]:
        raise ValueError("Unknown unit system")

    ucnv = unit_converter.UnitConverter()

    # Correlation is in field units
    if unit_system.lower() == "metric":
        temp = ucnv.convert(temp, "c", "f")
        pres = ucnv.convert(pres, "bar", "psi")

    sg = mw_gas / 28.9647

    # If gas z-factor is not given -> use Hall Yarborough
    if z_fac is None:
        z_fac = gas_zfac_hallyarborough(pres, temp, sg, "field", print_warning)

    # Sutton correlation for pseudo critical temperature (R) and
    # pseudo critical pressure (psia)
    # NOTE: correlation function should maybe be changed to allow
    # for t_pc and p_pc based on compositional data

    t_pc = 169.2 + 349.5 * sg - 74.0 * sg ** 2
    p_pc = 756.8 - 131.0 * sg - 3.6 * sg ** 2

    # pseudo reduced temperature and pressures
    t_pr = (temp + 460) / t_pc
    p_pr = pres / p_pc

    if t_pr > 40.0 or t_pr < 1.0:
        warnings.warn(
            "WARNING: Brill and Beggs z-factor correlation outside range Tr = %6.3f range is [1.5, 2.0]"
            % (t_pr)
        )
        # print("WARNING: Brill and Beggs z-factor correlation outside range")
        # print("Tr = %6.3f range is [1.5, 2.0]" % (t_pr))

    if p_pr > 100 or p_pr < 0.0:
        warnings.warn(
            "WARNING: Brill and Beggs z-factor correlation outside range Pr = %6.3f range is [0.2. 15]"
            % (p_pr)
        )
        # print("WARNING: Brill and Beggs z-factor correlation outside range")
        # print("Pr = %6.3f range is [0.2. 15]" % (p_pr))

    # Lucas correlation
    eps = 9490 * (t_pc / (mw_gas ** 3 * p_pc ** 4)) ** (1.0 / 6.0)

    ugsc_eps = (
        0.807 * t_pr ** 0.618
        - 0.357 * math.exp(-0.449 * t_pr)
        + 0.340 * math.exp(-4.058 * t_pr)
        + 0.018
    )

    ugsc = ugsc_eps / eps

    A1 = 1.245e-3 * math.exp(5.1726 * t_pr ** -0.3286) / t_pr

    A2 = A1 * (1.6553 * t_pr - 1.2723)

    A3 = 0.4489 * math.exp(3.0578 * t_pr ** -37.7332) / t_pr

    A4 = 1.7368 * math.exp(2.231 * t_pr ** -7.6351) / t_pr

    A5 = 0.9425 * math.exp(-0.1853 * t_pr ** 0.4489)

    ug_on_ugsc = 1.0 + A1 * p_pr ** 1.3088 / (
        A2 * p_pr ** A5 + (1 + A3 * p_pr ** A4) ** -1
    )

    gas_visc = ug_on_ugsc * ugsc

    return gas_visc


# Oil - bubble point pressure correlations


def oil_pbub_standing(rs, sto_api, sg, temp, unit_system="metric"):
    """
    Standing bubble-point correlation
    Reference: SPE phase behavior monograph eq. 3.79
    Application range : lighter volatile oils, low res. temperatures

    Parameters
    ----------
    rs : float
        single-stage gas-oil solution ratio (metric: sm3/sm3, field: scf/stb)
    sto_api : float
        stock tank oil API gravity
    sg : float
        specific gravity of surface gas
    temp : float
        reservoir temperature (metric: Celcius, field Farenheit)
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

    if unit_system.lower() not in ["metric", "field"]:
        raise ValueError("Unknown unit system")

    ucnv = unit_converter.UnitConverter().convert

    # Correlation is in field units
    if unit_system.lower() == "metric":

        rs = ucnv(rs, "m3/m3", "ft/bbl")
        temp = ucnv(temp, "c", "f")

    x = (rs / sg) ** 0.83 * 10 ** (0.00091 * temp - 0.0125 * sto_api)
    pbub = 18.2 * (x - 1.4)

    if unit_system.lower() == "metric":
        pbub = ucnv(pbub, "psi", "bar")

    return pbub


def oil_pbub_glaso(rs, sto_api, sg, temp, unit_system="metric"):
    """
    Glasoe bubble-point correlation
    Reference: SPE phase behavior monograph eq. 3.82
    Application range : Typical North sea oils

    Parameters
    ----------
    rs : float
        single-stage gas-oil solution ratio (metric: sm3/sm3, field: scf/stb)
    sto_api : float
        stock tank oil API gravity
    sg : float
        specific gravity of surface gas
    temp : float
        reservoir temperature (metric: Celcius, field Farenheit)
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

    if unit_system.lower() not in ["metric", "field"]:
        raise ValueError("Unknown unit system")

    ucnv = unit_converter.UnitConverter().convert

    # Correlation is in field units
    if unit_system.lower() == "metric":

        rs = ucnv(rs, "m3/m3", "ft/bbl")
        temp = ucnv(temp, "c", "f")

    A = (rs / sg) ** 0.816 * (temp ** 0.172 / sto_api ** 0.989)

    pbub = 10 ** (1.7669 + 1.7447 * math.log(A, 10) - 0.30218 * (math.log(A, 10)) ** 2)

    if unit_system.lower() == "metric":
        pbub = ucnv(pbub, "psi", "bar")

    return pbub


# Oil viscosity correlations


def oil_visc_pbub_standing(temp, sto_api, rs, unit_system="metric"):
    """
    Standing viscosity correlation (at pbub)
    Reference: SPE phase behavior monograph eq. 3.126

    Parameters
    ----------
    temp : float
        reservoir temperature (metric: Celcius, field: Farenheit)
    sto_api : float
        stock tank oil API gravity
    rs : float
        single-stage gas-oil solution ratio (metric: sm3/sm3, field: scf/stb)
    unit_system : str, optional
        "metric" or "field", by default "metric"

    Returns
    -------
    float
        Viscosity at pbub (cP)

    Raises
    ------
    ValueError
        Unknown unit system.
    """

    if unit_system.lower() not in ["metric", "field"]:
        raise ValueError("Unknown unit system")

    ucnv = unit_converter.UnitConverter().convert

    # Correlation is in field units
    if unit_system.lower() == "metric":
        rs = ucnv(rs, "m3/m3", "ft3/bbl")
        temp = ucnv(temp, "c", "f")

    # calculate dead oil viscosity
    A1 = 10 ** (0.43 + (8.33 / sto_api))
    vis_od = (0.32 + (1.8e7) / (sto_api ** 4.53)) * (360 / (temp + 200)) ** A1

    # calculate oil viscosity at pbub
    A1 = 10 ** (-rs * 7.4e-4 + (2.2e-7) * (rs ** 2))
    A2 = (
        (0.68 / (10 ** (rs * 8.62e-5)))
        + (0.25 / (10 ** (rs * 1.1e-3)))
        + (0.062 / (10 ** (rs * 3.74e-3)))
    )

    return A1 * (vis_od ** A2)


def oil_visc_pbub_bergman(temp, sto_api, rs, unit_system="metric"):
    """
    Bergman viscosity correlation (at pbub)
    Reference: SPE phase behavior monograph eq. 3.126

    Parameters
    ----------
    temp : float
        reservoir temperature (metric: Celcius, field: Farenheit)
    sto_api : float
        stock tank oil API gravity
    rs : float
        single-stage gas-oil solution ratio (metric: sm3/sm3, field: scf/stb)
    unit_system : str, optional
        "metric" or "field", by default "metric"

    Returns
    -------
    float
        Viscosity at pbub (cP)

    Raises
    ------
    ValueError
        Unknown unit system.
    """

    if unit_system.lower() not in ["metric", "field"]:
        raise ValueError("Unknown unit system")

    ucnv = unit_converter.UnitConverter().convert

    # Correlation is in field units
    if unit_system.lower() == "metric":
        rs = ucnv(rs, "m3/m3", "ft3/bbl")
        temp = ucnv(temp, "c", "f")

    # calculate dead oil viscosity
    A0 = 22.33 - 0.194 * (sto_api) + 0.00033 * (sto_api) ** 2
    A1 = -3.20 + (0.0185) * (sto_api)
    vis_od = -1 + math.exp(math.exp(A0 + A1 * math.log(temp + 310)))

    # calcuate oil viscosity at pbub
    A1 = math.exp(4.768 - 0.8359 * math.log(rs + 300))
    A2 = 0.555 + 133.5 / (rs + 300)

    return A1 * (vis_od ** A2)


# Oil Formation Volume Factor


def bob_standing(temp, rs, sg_oil, sg_gas, unit_system="metric"):
    """
    Standing oil formation volume factor at bubble-point pressure

    Reference: SPE phase behavior monograph eq.  3.111

    Application range : xxxx

    Input:
        temp    : reservoir temperature (metric: Celcius, field Farenheit)
        rs      : single-stage gas-oil solution ratio (metric: sm3/sm3, field: scf/stb )
        sg_gas  : surface gas specific gravity (dens/dens_air)
        sg_oil  : surface oil specific grabvity (dens/dens_wat)
        unit_system : "metric" or "field"

    Output:
        bob     : bo @ pbub

    """

    if unit_system.lower() not in ["metric", "field"]:
        raise ValueError("Unknown unit system")

    ucnv = unit_converter.UnitConverter().convert

    # Correlation is in field units
    if unit_system.lower() == "metric":
        rs = ucnv(rs, "m3/m3", "ft3/bbl")
        temp = ucnv(temp, "c", "f")

    A = rs * (sg_gas / sg_oil) ** 0.5 + 1.25 * temp

    bob = 0.9759 + 12e-5 * A ** 1.2

    return bob
