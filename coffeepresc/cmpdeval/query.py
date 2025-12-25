from dataclasses import dataclass
from pathlib import Path

import pandas as pd


@dataclass(frozen=True, slots=True)
class Query:
    f1: str
    f2: str
    dmin: float
    dmax: float
    pair_score: float
    rank1: int
    rank2: int
    score1: float
    score2: float

    def __str__(self) -> str:
        return (
            f"{self.f1} {self.f2} "
            f"{self.dmin:.1f} {self.dmax:.1f} "
            f"{self.rank1} {self.rank2} "
            f"{self.score1:.3f} {self.score2:.3f} "
            f"{self.pair_score:.3f}"
        )

def read_all_queries(file: Path | str) -> list[Query]:
    if not Path(file).exists():
        raise FileNotFoundError(f"Query file not found at {file}")
    df = pd.read_csv(file)
    required_cols = ["f_1", "f_2", "start", "end", "pair_score", "rank_1", "rank_2", "score_1", "score_2"]
    missing = [col for col in required_cols if col not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")
    if df[required_cols].isnull().any().any():
        raise ValueError("Missing values detected in required columns.")
    queries = []
    for _, row in df.iterrows():
        try:
            queries.append(Query(
                f1=row["f_1"],
                f2=row["f_2"],
                dmin=float(row["start"]),
                dmax=float(row["end"]),
                pair_score=float(row["pair_score"]),
                rank1=int(row["rank_1"]),
                rank2=int(row["rank_2"]),
                score1=float(row["score_1"]),
                score2=float(row["score_2"]),
            ))
        except Exception as e:
            raise ValueError(f"Row parse error: {e}\nRow: {row}")
    return queries