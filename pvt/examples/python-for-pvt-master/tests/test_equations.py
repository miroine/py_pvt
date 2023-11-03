#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
    PVT equation tests
    ~~~~~~~~~~~~~~~~~~~
    Test suite covering the PVT equations
"""

import pvt_tools.equations


class TestEquationsModule:
    def test_output_calc_mw_c7_plus(self):

        a = pvt_tools.equations.calc_mw_c7_plus(
            1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13
        )
        assert a == -215.67078333333333

    def test_output_calc_mw_res_fluid(self):

        a = pvt_tools.equations.calc_mw_res_fluid(
            1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1, 2, 3
        )
        assert a == 30.322694

    def test_output_calc_mw_rich_gas(self):

        a = pvt_tools.equations.calc_mw_rich_gas(6, 2, 3, 4, 5, 6, 3, 9)
        assert a == 40.3928

    def test_output_calc_vol_rich_gas(self):

        a = pvt_tools.equations.calc_vol_rich_gas(1, 2, 3, 4, 5, 6, 7, 8)
        assert a == 851.22

    def test_output_calc_vol_sto(self):

        a = pvt_tools.equations.calc_vol_sto(1, 2, 3, 4, 5, 6)
        assert a == 0.7338844147395526

    def test_output_calc_mass_rich_gas(self):

        a = pvt_tools.equations.calc_mass_rich_gas(1, 2, 3, 4, 5, 6, 7, 8)
        assert a == 1569.23

    def test_output_calc_mass_sto(self):

        a = pvt_tools.equations.calc_mass_sto(1, 2, 3, 4, 5)
        assert a == 494.9858

    def test_output_calc_dens_rich_gas(self):

        a = pvt_tools.equations.calc_dens_rich_gas(1, 2)
        assert a == 0.0005

    def test_output_calc_dens_sto(self):

        a = pvt_tools.equations.calc_dens_sto(1, 2)
        assert a == 0.0005

    def test_output_calc_sg_rich_gas(self):

        a = pvt_tools.equations.calc_sg_rich_gas(1)
        assert a == 0.034522042324023894

    def test_output_calc_sst_gor(self):

        a = pvt_tools.equations.calc_sst_gor(1, 2)
        assert a == 0.5

    def test_output_calc_sst_cgr(self):

        a = pvt_tools.equations.calc_sst_cgr(1)

        assert a == 1.0

    def test_output_calc_ratio_mass_rich_gas(self):

        a = pvt_tools.equations.calc_ratio_mass_rich_gas(1, 2)

        assert a == 3.0

    def test_output_calc_ratio_mass_sto(self):

        a = pvt_tools.equations.calc_ratio_mass_sto(1, 2)

        assert a == 1.5

    def test_output_calc_bo(self):

        a = pvt_tools.equations.calc_bo(1, 2, 4)

        assert a == 0.5

    def test_output_calc_bg(self):

        a = pvt_tools.equations.calc_bg(1, 2, 4)

        assert a == 0.5

    def test_output_calc_api(self):

        a = pvt_tools.equations.calc_api(1)

        assert a == 9.872649999999993
