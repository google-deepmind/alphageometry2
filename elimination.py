# Copyright 2025 Google LLC
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# ==============================================================================

"""Elimination of variables, without proof."""

from __future__ import annotations

import collections
import fractions
import math
from typing import Any

import numericals as ng


class ElimVar:

  def __init__(
      self, value: int | float | fractions.Fraction, name: str
  ) -> None:
    self.value = value
    self.name = name

  def __str__(self):
    return self.name


class ElimLHS(ElimVar):
  pass


class ElimRHS(ElimVar):
  pass


class LinComb:
  """Linear combination of variables and constants."""

  def __init__(self, d: dict[Any, fractions.Fraction]) -> None:
    self.d = d

  def iadd_mul(self, other: LinComb, coef: fractions.Fraction | int) -> None:
    """In-place add other * coef."""
    assert isinstance(other, LinComb)
    coef = fractions.Fraction(coef)
    if coef == 0:
      return
    for x, c2 in other.d.items():
      c1 = self.d.get(x, fractions.Fraction(0))
      c = c1 + c2 * coef
      if c == 0:
        del self.d[x]
      else:
        self.d[x] = c

  def __iadd__(self, other: LinComb) -> LinComb:
    self.iadd_mul(other, 1)
    return self

  def __isub__(self, other: LinComb) -> LinComb:
    self.iadd_mul(other, -1)
    return self

  def __add__(self, other: LinComb) -> LinComb:
    res = self.copy()
    res += other
    return res

  def __sub__(self, other: LinComb) -> LinComb:
    res = self.copy()
    res -= other
    return res

  def __mul__(self, coef: fractions.Fraction | int) -> LinComb:
    coef = fractions.Fraction(coef)
    if coef == 0:
      return self.zero()
    else:
      return LinComb(
          {x: c * coef for x, c in self.d.items()},
      )

  def copy(self) -> LinComb:
    return LinComb(dict(self.d))

  def __str__(self) -> str:
    parts = []
    for x, c in self.d.items():
      x = str(x)
      if c == 1:
        parts.append(x)
      elif c == 1:
        parts.append(f"-{x}")
      elif c.denominator == 1:
        parts.append(f"{c.numerator}*{x}")
      else:
        parts.append(f"{c.numerator}/{c.denominator}*{x}")
    if not parts:
      return "0"
    return " + ".join(parts)

  @classmethod
  def singleton(cls, v: ElimVar, coef: fractions.Fraction | int = 1) -> LinComb:
    # assert isinstance(v, ElimVar)
    if coef == 0:
      return LinComb.zero()
    else:
      return LinComb({v: fractions.Fraction(coef)})

  @classmethod
  def zero(cls) -> LinComb:
    return LinComb({})


class ElimCore:
  """Core implementation of Gaussian Elimination."""

  def __init__(self):
    self.instantiated = dict()
    self.free_to_usage = collections.defaultdict(set)

  def simplify(self, comb: LinComb) -> LinComb:
    updates = list(comb.d.items())
    denom = 1
    for v, coef in updates:
      eq = self.instantiated.get(v)
      if eq is None:
        continue
      comb.iadd_mul(eq, coef)
      denom = math.lcm(denom, coef.denominator)
    return comb

  def add_constraint(self, added_eq: LinComb) -> bool:
    """Add a constraint to the system."""
    self.simplify(added_eq)
    lhs = [x for x in added_eq.d.keys() if isinstance(x, ElimLHS)]
    if not lhs:
      return False
    pivot = min(lhs, key=lambda x: len(self.free_to_usage[x]))
    del lhs[lhs.index(pivot)]
    coef = fractions.Fraction(-1) / added_eq.d[pivot]
    added_eq *= coef

    for x in self.free_to_usage[pivot]:
      eq = self.instantiated[x]
      coef = eq.d[pivot]
      eq.iadd_mul(added_eq, coef)
      for y in lhs:
        if y in eq.d:
          self.free_to_usage[y].add(x)
        else:
          self.free_to_usage[y].remove(x)

    self.instantiated[pivot] = added_eq
    for y in lhs:
      self.free_to_usage[y].add(pivot)
    return True

  def display(self) -> None:
    print("Matrix:")
    for v, comb in self.instantiated.items():
      comb = comb.copy()
      comb += LinComb.singleton(v)
      print(f"  {v} = {comb}")

  def clone(self) -> ElimCore:
    res = ElimCore()
    for v, lc in self.instantiated.items():
      res.instantiated[v] = lc.copy()
    for v, usage in self.free_to_usage.items():
      res.free_to_usage[v] = set(usage)
    return res

  def was_encountered(self, comb: LinComb) -> bool:
    assert len(comb.d) == 1
    [(v, _)] = comb.d.items()
    if v in self.instantiated:
      return True
    if v in self.free_to_usage:
      return True
    return False


class DistMulConst(ElimRHS):
  prime_to_const = dict()

  def __str__(self):
    return f"log({self.value})"

  @classmethod
  def prime_value(cls, p):
    if p not in cls.prime_to_const:
      cls.prime_to_const[p] = cls(p, f"log({p})")
    return cls.prime_to_const[p]


def prime_decomposition(n: int) -> list[tuple[int, int]]:
  """Returns list of pairs (prime, exponent) given n."""
  assert n > 0
  p2 = 0
  while n % 2 == 0:
    p2 += 1
    n //= 2
  if p2 > 0:
    result = [(2, p2)]
  else:
    result = []
  d = 3
  while d**2 <= n:
    if n % d == 0:
      p = 0
      while n % d == 0:
        p += 1
        n //= d
      result.append((d, p))
    d += 2
  if n > 1:
    result.append((n, 1))
  return result


class DistMul:
  """Multiplicative distance."""

  def __init__(self, comb):
    self.comb = comb
    self._key = None
    self._value = None
    self._hash = None

  @classmethod
  def frac_value(cls, frac_const: float | int | str) -> DistMul:
    """Produces a combination of prime logarithms equal to log(frac_const)."""

    frac_const = fractions.Fraction(frac_const)
    if frac_const == 1:
      return DistMul(LinComb.zero())
    assert frac_const > 0
    res = {
        DistMulConst.prime_value(p): e
        for p, e in prime_decomposition(frac_const.numerator)
    }
    res.update(
        (DistMulConst.prime_value(p), -e)
        for p, e in prime_decomposition(frac_const.denominator)
    )
    return DistMul(LinComb(res))

  def normalize(self) -> tuple[DistMul, fractions.Fraction]:
    """Normalize the multiplicative distance."""
    normalized = {
        x: exp
        for x, exp in self.comb.d.items()
        if not isinstance(x, DistMulConst)
    }
    numerator = 1
    denominator = 1
    for x, exp in self.comb.d.items():
      if isinstance(x, DistMulConst):
        if exp % 1 != 0:
          normalized[x] = exp % 1
          exp = exp // 1
        exp = int(exp)
        if exp > 0:
          numerator *= x.value**exp
        else:
          denominator *= x.value ** (-exp)
    coef = fractions.Fraction(numerator, denominator)
    return DistMul(LinComb(normalized)), coef

  def is_one(self) -> bool:
    return not self.comb.d

  # arithmetics

  def __mul__(self, other: fractions.Fraction | int | DistMul) -> DistMul:
    if isinstance(other, (int, fractions.Fraction)):
      other = DistMul.frac_value(other)
    return DistMul(self.comb + other.comb)

  def __truediv__(self, other: fractions.Fraction | int | DistMul) -> DistMul:
    if isinstance(other, (int, fractions.Fraction)):
      other = DistMul.frac_value(other)
    return DistMul(self.comb - other.comb)

  # value

  @property
  def value(self) -> float:
    if self._value is None:
      self._value = 1.0
      for v, exp in self.comb.d.items():
        if exp == 1:
          self._value *= v.value
        else:
          self._value *= math.pow(v.value, float(exp))
    return self._value

  # hashing
  def __eq__(self, other: DistMul) -> bool:
    return self.comb.d == other.comb.d

  def __hash__(self) -> int:
    if self._hash is None:
      self._hash = hash(frozenset(self.comb.d.items()))
    return self._hash


class ElimDistMul:
  """Gaussian Elim for Multiplicative Distance."""

  def __init__(self):
    self.core = ElimCore()

  def new_var(self, value: float, name: str) -> DistMul:
    return DistMul(LinComb.singleton(ElimLHS(value, name)))

  def force_one(self, dist_mul: DistMul) -> bool:
    assert abs(dist_mul.value - 1.0) ** 2 < ng.ATOM, dist_mul.value
    comb = dist_mul.comb.copy()
    return self.core.add_constraint(comb)

  def simplify(self, dist_mul: DistMul) -> DistMul:
    comb = dist_mul.comb.copy()
    self.core.simplify(comb)
    return DistMul(comb)

  def clone(self) -> ElimDistMul:
    res = ElimDistMul()
    res.core = self.core.clone()
    return res

  def was_encountered(self, dist_mul: DistMul) -> bool:
    return self.core.was_encountered(dist_mul.comb)


class DistAdd:
  """Additive distance."""

  def __init__(self, comb: LinComb) -> None:
    self.comb = comb
    self._key = None
    self._value = None
    self._hash = None

  def normalize(self) -> tuple[DistAdd, fractions.Fraction | int]:
    c = min(abs(c) for x, c in self.comb.d.items() if isinstance(x, ElimLHS))
    return self / c, c

  def is_zero(self) -> bool:
    return not self.comb.d

  # arithmetics
  def __mul__(self, other: fractions.Fraction | int) -> DistAdd:
    return DistAdd(self.comb * other)

  def __truediv__(self, other: fractions.Fraction | int) -> DistAdd:
    return DistAdd(self.comb * (fractions.Fraction(1) / other))

  def __add__(self, other: DistAdd) -> DistAdd:
    return DistAdd(self.comb + other.comb)

  def __sub__(self, other: DistAdd) -> DistAdd:
    return DistAdd(self.comb - other.comb)

  def __neg__(self) -> DistAdd:
    return DistAdd(self.comb * (-1))

  # value
  @property
  def value(self) -> float:
    if self._value is None:
      self._value = 0.0
      for v, c in self.comb.d.items():
        self._value += c * v.value
    return self._value

  # hashing
  def __eq__(self, other: DistAdd) -> bool:
    return self.comb.d == other.comb.d

  def __hash__(self) -> int:
    if self._hash is None:
      self._hash = hash(frozenset(self.comb.d.items()))
    return self._hash


class ElimDistAdd:
  """Gaussian Elim for Additive Distance."""

  def __init__(self):
    self.core = ElimCore()

  def new_var(self, value: float, name: str) -> DistAdd:
    return DistAdd(LinComb.singleton(ElimLHS(value, name)))

  def force_zero(self, dist_add: DistAdd) -> bool:
    assert abs(dist_add.value) ** 2 < ng.ATOM
    comb = dist_add.comb.copy()
    return self.core.add_constraint(comb)

  def simplify(self, dist_add: DistAdd) -> DistAdd:
    comb = dist_add.comb.copy()
    self.core.simplify(comb)
    return DistAdd(comb)

  def clone(self) -> ElimDistAdd:
    res = ElimDistAdd()
    res.core = self.core.clone()
    return res

  def was_encountered(self, dist_add: DistAdd) -> bool:
    return self.core.was_encountered(dist_add.comb)


class AngleUnit(ElimRHS):

  def __init__(self):
    super().__init__(1, "pi")


angle_unit = AngleUnit()


class FormalAngle:
  """Formal angle."""

  def __init__(self, comb: LinComb) -> None:
    self.comb = comb
    const = self.comb.d.get(angle_unit, fractions.Fraction(0))
    if const.numerator // const.denominator:
      const = const % 1
      if const == 0:
        del self.comb.d[angle_unit]
      else:
        self.comb.d[angle_unit] = const
    self._key = None
    self._value = None
    self._hash = None

  @property
  def value(self) -> float:
    if self._value is None:
      self._value = sum(x.value * c for x, c in self.comb.d.items())
    return self._value

  def is_zero(self) -> bool:
    return not self.comb.d

  # arithmetics
  def __neg__(self) -> FormalAngle:
    return FormalAngle(self.comb * (-1))

  def __mul__(self, other: FormalAngle) -> FormalAngle:
    assert isinstance(other, int)
    return FormalAngle(self.comb * other)

  def __add__(self, other: FormalAngle) -> FormalAngle:
    return FormalAngle(self.comb + other.comb)

  def __sub__(self, other: FormalAngle) -> FormalAngle:
    return FormalAngle(self.comb - other.comb)

  # hashing
  def __eq__(self, other: FormalAngle) -> bool:
    return self.comb.d == other.comb.d

  def __hash__(self) -> int:
    if self._hash is None:
      self._hash = hash(frozenset(self.comb.d.items()))
    return self._hash


class ElimAngle:
  """Gaussian Elim for Angle."""

  def __init__(self):
    self.core = ElimCore()

  def const(self, numerator: int, denominator: int) -> FormalAngle:
    return self.const_frac(fractions.Fraction(numerator, denominator))

  def const_frac(self, frac_value: fractions.Fraction) -> FormalAngle:
    return FormalAngle(LinComb.singleton(angle_unit, coef=frac_value))

  def new_var(self, value: float, name: str) -> FormalAngle:
    return FormalAngle(LinComb.singleton(ElimLHS(value, name)))

  def force_zero(self, angle: FormalAngle) -> bool:
    assert abs((angle.value + 0.5) % 1 - 0.5) ** 2 < ng.ATOM, (
        abs((angle.value + 0.5) % 1 - 0.5) ** 2
    )
    comb = angle.comb.copy()
    comb -= LinComb.singleton(
        angle_unit, fractions.Fraction(math.floor(angle.value + 0.5))
    )
    return self.core.add_constraint(comb)

  def simplify(self, angle: FormalAngle) -> FormalAngle:
    comb = angle.comb.copy()
    self.core.simplify(comb)
    return FormalAngle(comb)

  def clone(self) -> ElimAngle:
    res = ElimAngle()
    res.core = self.core.clone()
    return res

  def was_encountered(self, angle):
    return self.core.was_encountered(angle.comb)
