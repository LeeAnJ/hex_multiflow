from dataclasses import dataclass
from typing import Literal

p = Literal['kpa', 'kPa', 'KPA', 'bar', 'Bar', 'BAR']
t = Literal['k', 'K', 'c', 'C', 'oC']
m_rate = Literal['kgs', 'mols']

# aliases for Units
Kelvin = float
Bar = float
Jmol = float


@dataclass(kw_only=True)
class Parameter:
    value: float
    units: str


@dataclass(kw_only=True)
class Pressure(Parameter):
    units: Literal[p]       # see above p = Literal['kpa', 'kPa', 'KPA', 'bar', 'Bar', 'BAR']


@dataclass(kw_only=True)
class Temperature(Parameter):
    units: Literal[t]       # see above t = Literal['k', 'K', 'c', 'C', 'oC']

@dataclass(kw_only=True)
class MassRate(Parameter):
    units: Literal[m_rate]       # see above _rate = Literal['kgs', 'mols']

@dataclass(kw_only=True)
class PthDat:
    p: Bar
    t: Kelvin
    h: Jmol


