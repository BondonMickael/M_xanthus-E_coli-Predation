from dataclasses import dataclass 
import pulp
from typing import List, Optional, Dict
from cobra import Model
import numpy as np


@dataclass
class BasePulpVarConfig:
    metabolicModel: Model
    RH: List[str]
    RM: List[str]
    RL: List[str]
    epsilon: float
    oxygenLevel: Optional[float] 
    
    def build_problem(self):
        self.prob = pulp.LpProblem('iMAT', pulp.LpMaximize)
        self.v_vars = self._create_flux_variables(self.prob)
        self.__add_mass_balance(self.prob, self.v_vars)
    
    def solve(self, solver: Optional[pulp.LpSolver] = None):
        '''
        Solves the LP problem
        '''
        if solver is None: 
            solver = pulp.PULP_CBC_CMD(msg=1)
        
        self.status = self.prob.solve(solver)
        if pulp.LpStatus[self.prob.status] == 'Infeasible':
            #Save infeasible File
            raise InterruptedError(f'Problem {self.prob.name} is INFEASIBLE and will be skipped')
        else:
            self.fluxes = {rid: self.v_vars[rid].varValue for rid in self.v_vars}
            self.y_values = {
                rid: [v.varValue for v in vals[1:3]] 
                for rid, vals in self.y_vars.items()
            }
            self.c_values = self.c_vars if hasattr(self, "c_vars") else None
            return self.status, self.fluxes, self.y_values, self.c_values
        
    def __add_mass_balance(
        self,
        prob: pulp.LpProblem,
        v_vars: Dict[str, pulp.LpVariable]
        ):
        for met in self.metabolicModel.metabolites:
            constraint = pulp.lpSum(rct.metabolites.get(met,0) * v_vars[rct.id] for rct in met.reactions) == 0
            prob += constraint, f'mass_balance:{met.id}'
    
    def _create_flux_variables(
        self,
        prob: pulp.LpProblem
        ) -> Dict[str, pulp.LpVariable]:
        
        v_vars: Dict[str, pulp.LpVariable] = {}
        
        for rct in self.metabolicModel.reactions:
            rid = rct.id 
            lb = rct.lower_bound
            ub = rct.upper_bound 
            
            if self.oxygenLevel is not None and rid == 'EX_o2_e':
                v = pulp.LpVariable(f'v_{rid}', self.oxygenLevel, self.oxygenLevel)
            #CII only runs forward
            elif rid == 'CII':
                v = pulp.LpVariable(f'v_{rid}', 0.0, 1000.0, cat = 'Continuous')
            # Glucose transporter constraint
            elif rid == 'GLCt1r':
                v = pulp.LpVariable(f'v_{rid}',lb,5.0, cat='Continuous')
            # Creatine/ Phosphocreatine exchange
            elif rid =='r0942':
                v = pulp.LpVariable(f'v_{rid}',-0.05,1000.0, cat='Continuous')
            elif rid == 'r0942b_mitoMap':
                v = pulp.LpVariable(f'v_{rid}',-0.05,1000.0, cat='Continuous')

            elif lb == np.inf or ub == np.inf:  
                v = pulp.LpVariable(f'v_{rid}', cat='Continuous') 
            else:
                v = pulp.LpVariable(f'v_{rid}', lb, ub, cat='Continuous')
            v_vars[rid] = v
        # Constrain CKc and CK to run in the same direction to avoid loop
        # y_creatine = pulp.LpVariable('y_CK_sign', cat = 'Binary')
        # v_CK = v_vars['CK']
        # v_CKc = v_vars['CKc']
        # M_ck = 1000.0
        # prob += v_CK   <=  M_ck * y_creatine, f"CK_pos_upper_if_y1"
        # prob += v_CK   >= -M_ck * (1 - y_creatine), f"CK_pos_lower_if_y1"
        # prob += v_CKc  <=  M_ck * y_creatine, f"CKc_pos_upper_if_y1"
        # prob += v_CKc  >= -M_ck * (1 - y_creatine), f"CKc_pos_lower_if_y1"

        return v_vars 
        
            
            
        
    
