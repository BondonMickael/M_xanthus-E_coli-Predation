from dataclasses import dataclass 
from pathlib import Path
from typing import List, Optional, Dict
import pulp 
import pandas as pd
import re
from utils.Discretizer import DiscretizationMethod


@dataclass
class CreateOutput:
    output_dir: Path   #directory to output path (e.g., /output)
    method: str #either iMAT or weighted_iMAT
    prob: pulp.LpProblem
    RH: List[str]
    RM: List[str]
    RL: List[str]
    flux_distribution: Optional[Dict[str,float]] 
    y_values: Optional[Dict[str,List[int]]]
    c_values: Optional[Dict[str,float]]
    epsilon: float
    oxygenLevel: Optional[float]
    expression_df: pd.DataFrame
    discretization_method: DiscretizationMethod
    quantiles: Optional[List[float]]
    
    
    #retrieve cell type Name 
    def __post_init__(self):
        raw_name = self.expression_df.columns[0]
        self.cell_type_name = re.sub(r"\.+", "_", raw_name)
        
    
    def create_output(self, fileName: Optional[str] = None): 
        file = self._generate_fileNames(fileName=fileName)
        self._create_file_flux_classification(file = file)
    
    def _create_output_dir(self): 
        '''
        Generates the directory where output files per cell type will be created 
        parentPath/method/cell_type
        Checks whether path already exists, if not creates it
        '''
        new_dir = Path(self.output_dir) / self.method / self.cell_type_name
        new_dir.mkdir(parents=True, exist_ok=True)
        return new_dir  
    
    def _generate_fileNames(self, fileName: Optional[str] = None):
        """
        Generates the full file path, where results will be stored. 
        The filename contains: 
        - epsilon
        - quantiles (if quantiles), else mean
        - optional fileName
        - optional oxygenLevel
        """
        if fileName is not None:
            fileName = Path(fileName).stem  

        base_dir = self._create_output_dir()
        
        if self.discretization_method.QUANTILE:
            base_file =  f'epsilon_{self.epsilon}_quantiles_{self.quantiles[0]}_{self.quantiles[1]}'
        elif self.discretization_method.MEAN:
            base_file =  f'epsilon_{self.epsilon}_mean'


        # final file name
        parts = [base_file]
        if fileName:
            parts.append(fileName)
        if self.oxygenLevel is not None:
            parts.append(f'oxygenLevel_{self.oxygenLevel}')
        
        # Combine parts with underscores
        final_file_name = '_'.join(parts) + '.tsv'
        full_file_path = base_dir / final_file_name

        return full_file_path     
    
    def _create_file_flux_classification(self, file:Path):
        '''
        Writes into the file: 
        - flux distribution with
            - reaction_id
            - flux value
            - classification 
            - y_f
            - y_r
            - c_value (if weighted_iMAT)
        '''
        with open(file, "w") as f:
            # Header
            header = ["reaction_id", "flux_value", "classification", 'y_f', 'y_r']
            if self.method == "weighted_iMAT":
                header.append("c_value")
            f.write("\t".join(header) + "\n")

            for rid, flux in self.flux_distribution.items():
                if rid in self.RH: 
                    classification = 'high'
                elif rid in self.RM: 
                    classification = 'moderate'
                elif rid in self.RL: 
                    classification = 'low'
                else:
                    classification = ''
                yf, yr = ("", "")
                if rid in self.y_values:
                    yf = self.y_values[rid][0]
                    yr = self.y_values[rid][1]
                row = [rid, str(flux), classification, str(yf), str(yr)]
                
                if self.method == "weighted_iMAT" and self.c_values:
                    c_val = self.c_values.get(rid, "")
                    row.append(str(c_val))
                
                f.write("\t".join(row) + "\n")

    
    
    def _create_infeasible_file(self): 
        '''
        Generates the file, where infeasible models will be printed to
        '''
        inf_dir = Path(self.dir).parent / 'infeasible_combinations'
        inf_dir.parent.mkdir(parents=True, exist_ok=True)
        return inf_dir
    
    def handle_infeasibility(self, fileName:Optional[str] = None):
        '''
        Writes info about an infeasible model to a file.
        Includes: 
        - problem name
        - cell type
        - epsilon
        - quantiles (if quantile) or mean
        - oxygenLevel
        - file Name if given
        '''
        dir_created = False
        if dir_created == False:
            inf_file = self._create_infeasible_dir() / 'infeasible_combinations.tsv'
            dir_created = True
        if dir_created:
            row = [
                self.prob.name, 
                self.cell_type_name, 
                f'epsilon = {self.epsilon}',
            ]
            if self.discretization_method.QUANTILE: 
                row.append(f'quantiles = {self.quantiles[0]}, {self.quantiles[1]}')
            else:
                row.append('mean')
            
            if self.oxygenLevel is not None:
                row.append(self.oxygenLevel)
            
            if fileName is not None:
                row.append(fileName)
            write_header = not inf_file.exists()
            # Append to file
            with open(inf_file, "a") as f:
                if write_header:
                    header = ["problem_name", "cell_type", "epsilon", "discretization", "oxygenLevel", 'InputFile']
                    f.write("\t".join(header) + "\n")
                f.write("\t".join(row) + "\n")
            
            