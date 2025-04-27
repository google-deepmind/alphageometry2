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

"""Basic datastructures for AlphaGeometry, and their parsing."""

import dataclasses
import fractions
import numpy as np


Fraction = fractions.Fraction


class AGPoint:
  """Alpha geometry point, containing a name and a numerical value."""

  def __init__(self, name: str, value) -> None:
    self.name = name
    self.value = value

  def __str__(self) -> str:
    return self.name


@dataclasses.dataclass
class AGPredicate:
  """Alpha geometry predicate."""

  name: str
  points: list[str | AGPoint]
  constants: list[int | float]

  def replace_points(self, ori_to_new):
    return AGPredicate(
        name=self.name,
        points=[ori_to_new[x] for x in self.points],
        constants=self.constants,
    )

  def __str__(self):
    tokens = [self.name]
    tokens.extend(map(str, self.points))
    tokens.extend(map(str, self.constants))
    return ' '.join(tokens)

  @classmethod
  def parse(cls, line):
    """Parse an AlphaGeometry predicate string."""
    pred, *args = line.split()
    points = []
    constants = []
    for arg in args:
      if arg[0].isnumeric() or arg[0] == '-':
        if arg.isnumeric() or arg[0] == '-' and arg[1:].isnumeric():
          constants.append(int(arg))
        elif 'pi/' in arg:
          num, den = arg.split('pi/')
          num = int(num)
          den = int(den)
          constants.append(Fraction(num * 180, den))
        else:
          num, den = arg.split('/')
          num = int(num)
          den = int(den)
          constants.append(Fraction(num, den))
      else:
        points.append(arg)
    return AGPredicate(name=pred, points=points, constants=constants)


@dataclasses.dataclass
class AGProblem:
  """Structure representing a low-level geometry problem for AlphaGeometry."""
  points: list[AGPoint]
  preds: list[AGPredicate]
  goal: AGPredicate

  def __str__(self):
    points = [f'{point.name} : {point.value}' for point in self.points]
    preds = list(map(str, self.preds))
    goal = '\nGoal: ' + str(self.goal)
    return '\n'.join(points + [''] + preds + [goal])

  def replace_points(self, ori_to_new):
    return AGProblem(
        points=[ori_to_new[x] for x in self.points],
        preds=[pred.replace_points(ori_to_new) for pred in self.preds],
        goal=self.goal.replace_points(ori_to_new),
    )

  def pstring(self):
    points = '; '.join(
        f'{x.name}@{x.value[0]}_{x.value[1]} = ' for x in self.points
    )
    preds = ', '.join(map(str, self.preds))
    goal = ' ? ' + str(self.goal)
    return points + preds + goal

  @classmethod
  def parse(cls, line):
    """Parse an AlphaGeometry problem string."""
    if '?' in line:
      steps, goal = line.split('?')
    else:
      steps = line
      goal = None
    preds = []
    name_to_point = {}
    for step in steps.split(';'):
      points, constraints = step.split('=')
      new_points = []
      for point in points.split():
        if '@' not in point:
          raise ValueError(f'Need a value for point {point}')
        name, value = point.split('@')
        new_points.append(name)
        x, y = value.split('_')
        x = float(x)
        y = float(y)
        value = np.array([x, y])
        name_to_point[name] = AGPoint(name=name, value=value)
      for constraint in constraints.strip().split(','):
        constraint = constraint.strip()
        if not constraint:
          continue
        preds.append(AGPredicate.parse(constraint))

    preds = list(map(lambda pred: pred.replace_points(name_to_point), preds))
    if goal is not None:
      goal = AGPredicate.parse(goal).replace_points(name_to_point)

    return AGProblem(
        points=list(name_to_point.values()),
        preds=preds,
        goal=goal,
    )
