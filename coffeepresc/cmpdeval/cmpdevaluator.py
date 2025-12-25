import numpy as np
import numpy.typing as npt
import pandas as pd

from coffeepresc.cmpdeval.query import read_all_queries
from coffeepresc.common.typing import StrPath
from coffeepresc.models.settings import CmpdEvalSettings


class CmpdEvaluator:
    def __init__(self, settings: CmpdEvalSettings) -> None:
        self.__settings = settings
        self.__done = False

    @property
    def settings(self) -> CmpdEvalSettings:
        return self.__settings

    @property
    def done(self) -> bool:
        return self.__done

    @property
    def result(self) -> pd.DataFrame:
        if not self.done:
            raise RuntimeError(
                "Evaluation not yet performed. Call run() before accessing result."
            )
        return self.__result

    @staticmethod
    def _select_rows(df: pd.DataFrame) -> pd.DataFrame:
        """extract rows with the best score for each mol_name"""
        return df.loc[df.groupby("mol_name")["score"].idxmin()].sort_values(by="score")

    def calc_fragscores(
        self,
        fragsims: np.ndarray,
        fragscores: np.ndarray,
    ) -> npt.NDArray:
        if np.any(fragsims < 0) or np.any(fragsims > 1):
            raise ValueError(
                f"Invalid fragsim. Expected 0 <= fragsim <= 1, but got {fragsims}"
            )
        return fragscores + np.maximum(0, -self.settings.penalty_coef * (fragsims - 1))

    def run(self, matched_path: StrPath, query_path: StrPath) -> None:
        df: pd.DataFrame = pd.read_csv(matched_path)
        query_list = read_all_queries(query_path)

        qobj = df["query_id"].apply(lambda i: query_list[i])  # type: ignore
        df = df.assign(
            query_obj=qobj,
            score_1=self.calc_fragscores(
                df["sim_1"].to_numpy(),
                np.array(qobj.apply(lambda q: q.score1)).astype(float),
            ),
            score_2=self.calc_fragscores(
                df["sim_2"].to_numpy(),
                np.array(qobj.apply(lambda q: q.score2)).astype(float),
            ),
        )
        df = df.assign(score=df["score_1"] + df["score_2"])

        best = self._select_rows(df)
        best = best.assign(query=best["query_obj"].astype(str))
        self.__result = best[["mol_name", "score", "query"]]
        self.__done = True

    def save(self, output_path: StrPath, query_row: bool = True) -> None:
        if not self.done:
            raise RuntimeError("Evaluation not yet performed. Call run() before save().")
        header = ["mol_name", "score"] if not query_row else ["mol_name", "score", "query"]
        self.__result[header].to_csv(output_path, index=False, header=None)
