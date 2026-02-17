from dataclasses import dataclass 
from methods.IMATConfig import IMATConfig
from methods.BasePulpVarConfig import BasePulpVarConfig
from cobra import Model
from typing import List, Optional, Dict
import pulp
import cobra 

@dataclass 
class iMAT(BasePulpVarConfig): 
    
    def build_problem(self): 
        super().build_problem()
        self._create_binary_variables()
        self._add_iMAT_constraints()
        self._add_objective_function()
        
    def _create_binary_variables(self):
        self.y_vars: Dict[str, tuple] = {}
        
        for rid in (self.RH + self.RL):
            #either yf or yr ==1, not both
            y_tot = pulp.LpVariable(f'y_tot_{rid}', cat = 'Binary')
            y_f = pulp.LpVariable(f'yf_{rid}', cat = 'Binary')
            y_r = pulp.LpVariable(f'yr_{rid}', cat='Binary')
            self.prob += y_tot == y_f + y_r, f'y_tot_def_{rid}'
            self.y_vars[rid] = (y_tot, y_f, y_r)

    def _add_iMAT_constraints(self):
        for rid in (self.RH + self.RL): 
            v = self.v_vars[rid]
            lb = self.metabolicModel.reactions.get_by_id(rid).lower_bound
            ub = self.metabolicModel.reactions.get_by_id(rid).upper_bound
            y_f = self.y_vars[rid][1]
            y_r = self.y_vars[rid][2]
            
            if rid in self.RH:
                self.prob += v + y_f*(lb - self.epsilon) >= lb, f'rH_forward_{rid}'
                self.prob += v + y_r*(ub + self.epsilon) <= ub, f'rH_reverse_{rid}'
            else:
                self.prob += lb * (1 - y_f) <= v, f'rL_leftHandSide_{rid}'
                self.prob += v <= (1 - y_f) * ub, f'rL_rightHandSide_{rid}'
            
    def _add_objective_function(self):
        terms = (
            [self.y_vars[rct][1] + self.y_vars[rct][2] for rct in self.RH] +
            [self.y_vars[rct][1] for rct in self.RL]
        ) 
        self.prob += pulp.lpSum(terms)