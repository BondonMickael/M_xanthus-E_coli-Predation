from dataclasses import dataclass
from methods.IMATConfig import IMATConfig
from methods.BasePulpVarConfig import BasePulpVarConfig
import utils.GPRMapper as GPRMapper
from cobra import Model
from typing import List, Optional, Dict
import pulp
import cobra
import pandas as pd


@dataclass
class weighted_iMAT(BasePulpVarConfig):
    gpr_mapper: GPRMapper
    lower_threshold_scaled: float
    upper_threshold_scaled: float

    def build_problem(self):
        super().build_problem()
        self._create_binary_variables()
        self._create_weight_variables()
        self._add_weighted_iMAT_constraints()
        self._add_objective_function()

    def _create_binary_variables(self):
        self.y_vars: Dict[str, tuple] = {}

        for rid in self.RH + self.RM + self.RL:
            y_tot = pulp.LpVariable(f"y_tot_{rid}", cat="Binary")
            y_f = pulp.LpVariable(f"yf_{rid}", cat="Binary")
            y_r = pulp.LpVariable(f"yr_{rid}", cat="Binary")
            self.prob += y_tot == y_f + y_r, f"y_tot_def_{rid}"
            self.y_vars[rid] = (y_tot, y_f, y_r)

    def _create_weight_variables(self):
        self.c_vars: Dict[str, float] = {}
        for rid in self.RH + self.RM:
            gpr = self.metabolicModel.reactions.get_by_id(rid).gpr
            x = self.gpr_mapper.get_reaction_expression(str(gpr))

            if rid in self.RH:
                c = (x - self.upper_threshold_scaled) + 1
            else:
                c = (x - self.lower_threshold_scaled) / (
                    self.upper_threshold_scaled - self.lower_threshold_scaled
                )
            self.c_vars[rid] = c

    def _add_weighted_iMAT_constraints(self):
        for rid in self.RH + self.RM + self.RL:
            v = self.v_vars[rid]
            lb = self.metabolicModel.reactions.get_by_id(rid).lower_bound
            ub = self.metabolicModel.reactions.get_by_id(rid).upper_bound
            y_f = self.y_vars[rid][1]
            y_r = self.y_vars[rid][2]

            if rid in (self.RH + self.RM):
                self.prob += v + y_f * (lb - self.epsilon) >= lb, f"rH_forward_{rid}"
                self.prob += v + y_r * (ub + self.epsilon) <= ub, f"rH_reverse_{rid}"
            else:
                self.prob += lb * (1 - y_f) <= v, f"rL_rM_leftHandSide_{rid}"
                self.prob += v <= (1 - y_f) * ub, f"rL_rM_rightHandSide_{rid}"

    def _add_objective_function(self):
        terms = [
            (self.y_vars[rct][1] + self.y_vars[rct][2]) * self.c_vars[rct]
            for rct in (self.RH + self.RM)
        ] + [self.y_vars[rct][1] for rct in self.RL]
        self.prob += pulp.lpSum(terms)
