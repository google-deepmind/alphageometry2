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

"""Numerical implementation of euclidean geometry."""

import numpy as np

ATOM = 1e-12
NumPoint = np.ndarray


def distance(a: NumPoint, b: NumPoint) -> float:
  return np.linalg.norm(a - b)


def normalize(v: NumPoint) -> NumPoint:
  return v / np.linalg.norm(v)


def perp_rot(v: NumPoint) -> NumPoint:
  [x, y] = v
  return np.array([y, -x])


def direction(v: NumPoint) -> float:
  [x, y] = v
  return np.arctan2(y, x) / np.pi


def midpoint(a: NumPoint, b: NumPoint) -> NumPoint:
  return (a + b) / 2


def orientation(a: NumPoint, b: NumPoint, c: NumPoint) -> int:
  matrix = np.stack((b - a, c - a))
  det = np.linalg.det(matrix)
  if det > ATOM:
    return 1
  elif det < -ATOM:
    return -1
  else:
    return 0


def collinear(a: NumPoint, b: NumPoint, c: NumPoint) -> bool:
  return orientation(a, b, c) == 0


class NumLine:
  """A point x is in the line if x*n = c where n is a vector of unit length."""

  def __init__(self, n: NumPoint, c: float):
    assert np.isclose(np.linalg.norm(n), 1), n
    self.n: NumPoint = n
    self.c: float = c

  @classmethod
  def through1(cls, n: NumPoint, a: NumPoint) -> "NumLine":
    """A line passing through `a` with a normal vector `n`."""
    return cls(n, np.dot(a, n))

  @classmethod
  def through(cls, a: NumPoint, b: NumPoint) -> "NumLine":
    return cls.through1(perp_rot(normalize(b - a)), a)

  def direction(self) -> float:
    return (direction(self.n) + 0.5) % 1

  def distance(self, a: NumPoint) -> float:
    return abs(self.c - np.dot(a, self.n))

  def position(self, a: NumPoint) -> float:
    """Position of a on the line as a single number."""
    return -np.dot(perp_rot(self.n), a)


def intersect_ll(line1: NumLine, line2: NumLine):
  matrix = np.stack((line1.n, line2.n))
  b = np.array((line1.c, line2.c))
  if abs(np.linalg.det(matrix)) < ATOM:
    return None
  return np.linalg.solve(matrix, b)


def perp_bisector(a: NumPoint, b: NumPoint) -> NumLine:
  return NumLine.through1(normalize(b - a), midpoint(a, b))


class NumCircle:
  """A point x is on a circle if its distance from center = r."""

  def __init__(self, center: NumPoint, r: float):
    self.center: NumPoint = center
    self.r: float = r

  @classmethod
  def through1(cls, center: NumPoint, a: NumPoint) -> "NumCircle":
    r = distance(a, center)
    return NumCircle(center, r)

  @classmethod
  def through(cls, a: NumPoint, b: NumPoint, c: NumPoint) -> "NumCircle":
    center = intersect_ll(
        perp_bisector(a, b),
        perp_bisector(a, c),
    )
    return NumCircle.through1(center, a)

  def distance(self, a: NumPoint) -> float:
    return abs(distance(self.center, a) - self.r)
