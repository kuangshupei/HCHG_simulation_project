from dataclasses import dataclass

@dataclass
class GasPropertyData:
    name: str
    beta2: float          # [fs^2/(cm bar)]
    beta3: float          # [fs^3/(cm bar)]
    beta4: float          # [fs^4/(cm bar)]
    n2_at_one_bar: float  # [m^2/(W bar)]
    
    def __repr__(self):
        return (f"{self.name}: beta2={self.beta2} fs^2/(cm bar), "
                f"beta3={self.beta3} fs^3/(cm bar), "
                f"beta4={self.beta4} fs^4/(cm bar), "
                f"n2_at_one_bar={self.n2_at_one_bar} m^2/(W bar)")

noble_gases = [
    GasPropertyData(name="Neon", beta2=0.0202, beta3=0.0158, beta4=0, n2_at_one_bar=0.14e-23),
    GasPropertyData(name="Argon", beta2=0.1980, beta3=0.1586, beta4=0, n2_at_one_bar=1.74e-23),
    GasPropertyData(name="Krypton", beta2=0.3996, beta3=0.3298, beta4=0, n2_at_one_bar=4.03e-23),
    GasPropertyData(name="Xenon", beta2=0.9113, beta3=0.7836, beta4=0, n2_at_one_bar=11.15e-23),
]

def loadGasParameters(gas_name, gases):
    
    gas_dict = {gas.name.lower(): gas for gas in gases}
    gas = gas_dict.get(gas_name.lower(), "Gas not found.")

    n2_at_one_bar = gas.n2_at_one_bar
    beta2 = gas.beta2 * (1e-3)**2 / 1e-5    # [ps^2/(km bar)]
    beta3 = gas.beta3 * (1e-3)**3 / 1e-5    # [ps^3/(km bar)]
    beta4 = gas.beta4 * (1e-3)**4 / 1e-5    # [ps^4/(km bar)]
    betas_at_one_bar = [beta2, beta3, beta4]
    return n2_at_one_bar, betas_at_one_bar