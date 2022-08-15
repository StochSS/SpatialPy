# SpatialPy is a Python 3 package for simulation of
# spatial deterministic/stochastic reaction-diffusion-advection problems
# Copyright (C) 2019 - 2022 SpatialPy developers.

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU GENERAL PUBLIC LICENSE Version 3 as
# published by the Free Software Foundation.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU GENERAL PUBLIC LICENSE Version 3 for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

from spatialpy.core.spatialpyerror import BoundaryConditionError


class BoundaryCondition():
    """
    Set spatial regions of the domain where a property of
    particles are held constant (updated each simulation step)

    Conditions (one or more of the following must be set):
         - xmin, xmax: (float) min or max value in the x dimension
         - ymin, ymax: (float) min or max value in the y dimension
         - zmin, zmax: (float) min or max value in the z dimension
         - type_id: type (subdomain) of the partciles
    Targets (one of the following must be set):
        property: (str), 'nu', 'rho','v'
        species: (str) name of a chemical species. \
                   Must also set deterministic=True/False flag.

    :param xmin: x-axis coordinate lower bound of **condition**
    :type xmin: float

    :param xmax: x-axis coordinate upper bound of **condition**
    :type xmax: float

    :param ymin: y-axis coordinate lower bound of **condition**
    :type ymin: float

    :param ymax: y-axis coordinate upper bound of **condition**
    :type ymax: float

    :param zmin: z-axis coordinate lower bound of **condition**
    :type zmin: float

    :param zmax: z-axis coordinate upper bound of **condition**
    :type zmax: float

    :param type_id: Set **condition** to particle type id
    :type type_id: int | str

    :param deterministic: **Must be set if target is Species.** Set True if boundary condition target is species \
    and applies to deterministic simulation. **BoundaryCondition not yet implemeneted for Stochastic Species**.
    :type deterministic: bool

    :param target: Set **target** to properties, can be 'nu' 'rho' or 'v' or species name \
    If species name, determinstic must also be set to True/False.
    :type target: str

    :param value: Value property will take in region defined by the conditions
    :type value: float or float[3]

    :param model: Target model of boundary condition
    :type model: spatialpy.core.model.Model
    """

    def __init__(self, xmin=None, xmax=None, ymin=None, ymax=None, zmin=None, zmax=None,
                 type_id=None, deterministic=True, target=None, value=None, model=None):
        from spatialpy.core.model import Model # pylint: disable=import-outside-toplevel
        from spatialpy.core.species import Species # pylint: disable=import-outside-toplevel
        if xmin is not None and not isinstance(xmin, (int, float)):
            raise BoundaryConditionError("X-min must be of type int or float.")
        if xmax is not None and not isinstance(xmax, (int, float)):
            raise BoundaryConditionError("X-max must be of type int or float.")
        if ymin is not None and not isinstance(ymin, (int, float)):
            raise BoundaryConditionError("Y-min must be of type int or float.")
        if ymax is not None and not isinstance(ymax, (int, float)):
            raise BoundaryConditionError("Y-max must be of type int or float.")
        if zmin is not None and not isinstance(zmin, (int, float)):
            raise BoundaryConditionError("Z-min must be of type int or float.")
        if zmax is not None and not isinstance(zmax, (int, float)):
            raise BoundaryConditionError("Z-max must be of type int or float.")
        if type_id is not None and not isinstance(type_id, (int, str)):
            raise BoundaryConditionError("Type-ID must be of type int.")
        elif type_id is not None:
            if "UnAssigned" in type_id:
                raise BoundaryConditionError("'UnAssigned' is not a valid type_id")
            type_id = f"type_{type_id}"
        if target is None or not (isinstance(target, (str, Species)) or
                type(target).__name__ == 'Species' or property in ('nu', 'rho', 'v')):
            raise BoundaryConditionError("Target must be of type string or SpatialPy.Species")
        if not (value is None or isinstance(value, (int, float)) or (isinstance(value, list) and len(value) == 3)):
            raise BoundaryConditionError("Value must be of type float or float[3].")
        try:
            if isinstance(value, int):
                value = float(value)
            if isinstance(value, list):
                value = [float(val) for val in value]
        except Exception as err:
            raise BoundaryConditionError("Value must be of type float or float[3].") from err
        if not (model is None or isinstance(model, Model) or type(model).__name__ == 'Model'):
            raise BoundaryConditionError("Model must be of type SpatialPy.Model.")

        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax
        self.zmin = zmin
        self.zmax = zmax
        self.type_id = type_id
        self.target = target
        self.deterministic = deterministic
        self.value = value
        self.model = model


    def expression(self):
        """
        Creates evaluable string expression of boundary condition.

        :returns: A string representation of the boundary condition.
        :rtype: str
        """
        if self.value is None:
            raise BoundaryConditionError("Must set value")
        cond=[]
        if self.xmin is not None:
            cond.append(f"(me->x[0] >= {self.xmin})")
        if self.xmax is not None:
            cond.append(f"(me->x[0] <= {self.xmax})")
        if self.ymin is not None:
            cond.append(f"(me->x[1] >= {self.ymin})")
        if self.ymax is not None:
            cond.append(f"(me->x[1] <= {self.ymax})")
        if self.zmin is not None:
            cond.append(f"(me->x[2] >= {self.zmin})")
        if self.zmax is not None:
            cond.append(f"(me->x[2] <= {self.zmax})")
        if self.type_id is not None:
            cond.append(f"(me->type == {self.type_id})")
        if len(cond)==0:
            raise BoundaryConditionError('need at least one condition on the BoundaryCondition')
        bcstr = "if(" + '&&'.join(cond) + "){"
        if self.target in self.model.listOfSpecies:
            if self.deterministic:
                s_ndx = self.model.species_map[self.model.listOfSpecies[self.target]]
                bcstr += f"me->C[{s_ndx}] = {self.value};"
            else:
                raise BoundaryConditionError(
                    "BoundaryConditions don't work for stochastic species yet."
                )
        elif self.target is not None:
            if self.target == 'v':
                for i in range(3):
                    bcstr+= f"me->v[{i}]={self.value[i]};"
            elif self.target == 'nu':
                bcstr += f"me->nu={self.value};"
            elif self.target == 'rho':
                bcstr += f"me->rho={self.value};"
            else:
                raise BoundaryConditionError(f"Unable handle boundary condition for property '{self.target}'")
        bcstr+= "}"
        return bcstr
