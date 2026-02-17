from enum import Enum 
from typing import List
from dataclasses import dataclass
import pandas as pd
import numpy as np

class DiscretizationMethod(Enum):
    MEAN = 'mean'
    QUANTILE = 'quantile'
    

@dataclass
class DiscretizationResult:
    thresholds: List[float]
    scaled_thresholds: List[float]
    dataframe: pd.DataFrame
    
class Discretizer:
    def __init__(
        self,
        method: str,
        quantiles: List[int]         
    ):
        self.method = DiscretizationMethod(method)
        self.quantiles = quantiles

        self.lower_threshold = None
        self.upper_threshold = None
        self.lower_threshold_scaled = None
        self.upper_threshold_scaled = None
        
    def run(self, df: pd.DataFrame) -> DiscretizationResult:
        df['scaled_expression'] = (df.iloc[:,0] - df.iloc[:,0].min()) / (df.iloc[:,0].max() - df.iloc[:,0].min()) 

        if self.method is DiscretizationMethod.MEAN:
            self._mean_discretization(df)
        elif self.method is DiscretizationMethod.QUANTILE:
            self._quantile_discretization(df)
        else:
            raise ValueError(f"Unsupported discretization method: {self.method}")

        return DiscretizationResult(
            thresholds=[self.lower_threshold, self.upper_threshold],
            scaled_thresholds=[
                self.lower_threshold_scaled,
                self.upper_threshold_scaled,
            ],
            dataframe=df,
        )

    def _mean_discretization(self, df: pd.DataFrame):
        mean = df.iloc[:,0].mean()
        sd =  df.iloc[:,0].std()

        self.upper_threshold = mean + 0.5 * sd
        self.lower_threshold = mean - 0.5 * sd

        self.upper_threshold_scaled = (
            df["scaled_expression"].mean() + 0.5 * df["scaled_expression"].std()
        )
        self.lower_threshold_scaled = (
            df["scaled_expression"].mean() - 0.5 * df["scaled_expression"].std()
        )
        self._apply_thresholds(df)
        

    def _quantile_discretization(self, df: pd.DataFrame):
        if not self.quantiles or len(self.quantiles) != 2:
            raise ValueError("Quantile discretization requires two quantiles")

        q1, q2 = self.quantiles

        self.lower_threshold = df.iloc[:,0].quantile(q1 / 100)
        self.upper_threshold = df.iloc[:,0].quantile(q2 / 100)

        self.lower_threshold_scaled = df["scaled_expression"].quantile(q1 / 100)
        self.upper_threshold_scaled = df["scaled_expression"].quantile(q2 / 100)

        self._apply_thresholds(df)
        
    def _apply_thresholds(self, df: pd.DataFrame):
        rules = [
            df.iloc[:,0] > self.upper_threshold,
            df.iloc[:,0] < self.lower_threshold,
            df.iloc[:,0].between(self.lower_threshold, self.upper_threshold)
        ]
        classes = [1, -1, 0]
        df["discretization"] = np.select(rules, classes, default=np.nan)