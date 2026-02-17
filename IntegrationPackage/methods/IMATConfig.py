import numpy as np 
import pandas as pd
import pulp
from dataclasses import dataclass, field
from pathlib import Path
from utils.Discretizer import Discretizer, DiscretizationMethod
from utils.GPRMapper import GPRMapper
from cobra import Model
from typing import Optional, List 

@dataclass
class IMATConfig:
    '''
    expression_df: RNA-seq Dataframe (rows: gene, columns: samples) with metabolic-mapped data. NaN entries allowed.
    discretization_method: string, either mean or quantile
    quantiles: List of lower and upper quantile, default [40,70]
    epsilon: Threshold of active reaction, default = 1
    outputPath: Directory where output will be saved
    oxygenLevel: Oxygen uptake constraint value.
     metabolicModel : cobra.Model
        SBML metabolic model loaded via COBRApy.
    '''
    expression_df: pd.DataFrame
    discretization_method: DiscretizationMethod 
    quantiles: Optional[List[int]] 
    epsilon: Optional[float]
    metabolicModel: Model
    
    def __post_init__(self):
        if self.epsilon is None:
            self.epsilon = 1.0
        
        # Find out whether human or mice gene IDs
        if self.expression_df.index[0].startswith('MXAN'):
            self.ignore_human = True
        elif self.expression_df.index[0].startswith('ENSG'):
            self.ignore_human = False
        else:
            raise ValueError('Expression Dataframe seems to not contain mice or human species.')
        
        if isinstance(self.discretization_method, str):
            self.discretization_method = DiscretizationMethod(
                self.discretization_method
            )
        #Quantile Check 
        if self.discretization_method.QUANTILE:
            if self.quantiles is None: 
                self.quantiles = [40, 70]
                
            if self.quantiles[0] >= self.quantiles[1]:
                raise ValueError('Lower quantile must be smaller than upper quantile.')

            
    def prepare(self):
        '''
        Runs full preprocessing pipeline:
        1. Discretize expression
        2. Map GPR to reaction
        '''
        self._apply_discretization()
        self._map_GPR_to_reaction()
        
        
    def _apply_discretization(self):
        discretizer = Discretizer(
            method = self.discretization_method, 
            quantiles = self.quantiles)
        results = discretizer.run(self.expression_df)
        self.expression_df = results.dataframe
        self.lower_threshold_scaled, self.upper_threshold_scaled = results.scaled_thresholds
        
    def _map_GPR_to_reaction(self):
        self.gpr_mapper = GPRMapper(self.metabolicModel, self.expression_df, self.ignore_human)
        output = self.gpr_mapper.create_reaction_classes()
        self.RL = output.RL
        self.RM = output.RM
        self.RH = output.RH

    