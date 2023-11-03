# -*- coding: utf-8 -*-
"""
Petroleum unit conversions
"""


class UnitConverter:
    """
    Petroleum unit converter class
    """

    # The converter_dict is the dictionary containing the
    # conversion factors
    # The keys contain the units, and the values corresponding to each
    # unit key is a tuplet containing "type", "a1" and "a0" used for the
    # conversions, where the constants "a1" and "a0" are the conversion
    # constants: value_at_ref_unit = a1*value + a0
    # The first key for each "type" defines the reference unit, and
    # subsequent units are converted according the a1 and a0 values

    converter_dict = {
        # lengt conversions (reference is m3)
        "m": ("length", 1.0, 0.0),
        "cm": ("length", 0.01, 0.0),
        "ft": ("length", 0.3048, 0.0),
        # volume conversions (reference is m3)
        "m3": ("volume", 1.0, 0.0),
        "ft3": ("volume", 1.0 / 35.3147, 0.0),
        "litre": ("volume", 1.0e-3, 0.0),
        "bbl": ("volume", 1.0 / 6.28981, 0.0),
        # pressure conversions (reference is bar)
        "bar": ("pressure", 1.0, 0.0),
        "psi": ("pressure", 1.0 / 14.50377, 0.0),
        "pascal": ("pressure", 1.0e-5, 0.0),
        # temperature conversions (reference is c)
        "c": ("temperature", 1.0, 0.0),
        "k": ("temperature", 1.0, -273.15),
        "f": ("temperature", 1.0 / 1.8, -32.0 / 1.8),
        "r": ("temperature", 1.0 / 1.8, -491.67 / 1.8),
        # mass (reference is kg)
        "kg": ("mass", 1.0, 0.0),
        "g": ("mass", 1.0e-3, 0.0),
        "lb": ("mass", 0.45359237, 0.0),
        # density (reference is kg/m3)
        "kg/m3": ("density", 1.0, 0.0),
        "g/cm3": ("density", 1000.0, 0.0),
        "lb/ft3": ("density", 16.0185, 0.0),
        # volume ratio (reference is m3/m3)
        "m3/m3": ("vol_ratio", 1.0, 0.0),
        "ft3/bbl": ("vol_ratio", 1.0 / 5.614583, 0.0),
        # compressibility (reference is 1/bar)
        "1/bar": ("compressibility", 1.0, 0.0),
        "1/psi": ("compressibility", 14.50377, 0.0),
        # volumetric rate (reference is m3/d)
        "m3/d": ("vol_rate", 1.0, 0.0),
        "m3/h": ("vol_rate", 24.0, 0.0),
        "ft3/d": ("vol_rate", 1.0 / 35.3147, 0.0),
        "ft3/h": ("vol_rate", 24.0 / 35.3147, 0.0),
        "bbl/d": ("vol_rate", 1.0 / 6.28981, 0.0),
        "bbl/h": ("vol_rate", 24.0 / 6.28981, 0.0),
    }

    # -------------------------------------------------------------------------
    def convert(self, value, from_unit, to_unit):
        """
        Convert a value from "from_unit" to "to_unit"
        """

        try:
            value = float(value)
        except:
            raise TypeError(str(value) + " is not of float number format")

        from_unit = from_unit.lower()
        to_unit = to_unit.lower()

        if from_unit not in self.converter_dict.keys():
            raise KeyError(from_unit + " is not a defined unit in converter_dict")

        if to_unit not in self.converter_dict.keys():
            raise KeyError(to_unit + " is not a defined unit in converter_dict")

        from_type = self.converter_dict[from_unit][0]
        to_type = self.converter_dict[to_unit][0]

        if from_type != to_type:
            raise TypeError("Incompatible unit conversion ")

        to_a1 = self.converter_dict[from_unit][1]
        to_a0 = self.converter_dict[from_unit][2]

        ref_unit_value = to_a1 * value + to_a0

        from_a1 = self.converter_dict[to_unit][1]
        from_a0 = self.converter_dict[to_unit][2]

        to_unit_value = (ref_unit_value - from_a0) / from_a1

        return to_unit_value


def main():
    """
    For testing purposes
    """

    converter = UnitConverter()

    value = 1.0
    from_unit = "m3"
    to_unit = "ft3"
    to_value = converter.convert(value, from_unit, to_unit)
    print(value, from_unit, "is", to_value, to_unit)

    value = 1.0
    from_unit = "BBL"
    to_unit = "ft3"
    to_value = converter.convert(value, from_unit, to_unit)
    print(value, from_unit, "is", to_value, to_unit)

    value = 60
    from_unit = "F"
    to_unit = "C"
    to_value = converter.convert(value, from_unit, to_unit)
    print(value, from_unit, "is", to_value, to_unit)


if __name__ == "__main__":
    main()
