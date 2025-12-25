import tempfile
import tomllib
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Self

import numpy as np
import numpy.typing as npt

from coffeepresc.common.typing import StrPath


class FragQuerySettings:
    # docking region
    inner_width: npt.NDArray[np.int32]
    outer_width: npt.NDArray[np.int32]
    center: npt.NDArray[np.float64]
    score_pitch: npt.NDArray[np.float64]
    # parameters
    candidate_poses: int
    cluster_size: float
    distance_width: float
    no_oberrorlog: bool
    log: Path | None

    def __init__(
        self,
        docking_config: StrPath,
        log: StrPath | None = None,
        candidate_poses: int = 40,
        cluster_size: float = 1.0,
        distance_width: float = 0.1,
        no_oberrorlog: bool = False,
        **kwargs,
    ) -> None:
        self._parse_docking_config(docking_config)
        self.log = Path(log) if log is not None else None
        self.candidate_poses = candidate_poses
        self.cluster_size = cluster_size
        self.distance_width = distance_width
        self.no_oberrorlog = no_oberrorlog

    def validate(self) -> None:
        if self.candidate_poses <= 0:
            raise ValueError(
                "The number of promising poses must be a positive integer."
            )
        if self.cluster_size <= 0:
            raise ValueError("Cluster size must be a positive value.")
        if self.distance_width <= 0:
            raise ValueError("Distance width must be a positive value.")

    def _parse_docking_config(self, config_path: StrPath) -> None:
        with open(config_path, "r") as f:
            lines = f.readlines()
            for line in lines:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                k, v = line.split(" ", 1)
                if k == "INNERBOX":
                    self.inner_width = np.array(
                        list(map(int, v.split(","))), dtype=np.int32
                    )
                elif k == "OUTERBOX":
                    self.outer_width = np.array(
                        list(map(int, v.split(","))), dtype=np.int32
                    )
                elif k == "BOX_CENTER":
                    self.center = np.array(
                        list(map(float, v.split(","))), dtype=np.float64
                    )
                elif k == "SCORING_PITCH":
                    self.score_pitch = np.array(
                        list(map(float, v.split(","))), dtype=np.float64
                    )
                else:
                    warnings.warn(
                        f"Skipping key '{k}' in the config file.", RuntimeWarning
                    )

    def dumps_docking_config(self) -> str:
        self.validate()
        lines = []
        lines.append(f"INNERBOX {', '.join(map(str, self.inner_width))}")
        lines.append(f"OUTERBOX {', '.join(map(str, self.outer_width))}")
        lines.append(f"BOX_CENTER {', '.join(map(str, self.center))}")
        lines.append(f"SCORING_PITCH {', '.join(map(str, self.score_pitch))}")
        if self.log is not None:
            lines.append(f"LOG {self.log}")
        lines.append(f"PROMISING_POSE {self.candidate_poses}")
        lines.append(f"CLUSTER_SIZE {self.cluster_size}")
        lines.append(f"DISTANCE_WIDTH {self.distance_width}")
        lines.append(f"NO_OBERRORLOG {str(self.no_oberrorlog).lower()}")
        return "\n".join(lines)

    def dump_docking_config(self, path: StrPath) -> None:
        with open(path, "w") as f:
            f.write(self.dumps_docking_config())


class RetrievalSettings:
    storage_path: Path
    similarity_th: float

    def __init__(
        self,
        storage: StrPath = "data.h5",
        similarity_th: float = 0.45,
        **kwargs,
    ) -> None:
        self.storage_path = Path(storage)
        self.similarity_th = similarity_th

    def validate(self) -> None:
        # path
        if not self.storage_path.exists():
            raise FileNotFoundError(
                f"Storage file '{self.storage_path}' does not exist."
            )
        # values
        if not (0.0 <= self.similarity_th <= 1.0):
            raise ValueError("Similarity threshold must be between 0.0 and 1.0.")


class CmpdEvalSettings:
    penalty_coef: float

    def __init__(
        self,
        penalty_coef: float = 8.0,
        **kwargs,
    ) -> None:
        self.penalty_coef = penalty_coef

    def validate(self) -> None:
        # values
        if self.penalty_coef < 0:
            raise ValueError("Score slope must be a non-negative value.")
       

class RepclusSettings:
    clustering_molecules: Path
    clustering_output: Path
    n_clusters: int

    def __init__(
        self,
        clustering_molecules: StrPath,
        clustering_output: StrPath,
        n_clusters: int = 5,
        **kwargs,
    ) -> None:
        if clustering_molecules is None:
            raise ValueError("clustering_molecules must be specified.")
        if clustering_output is None:
            raise ValueError("clustering_output must be specified.")
        if n_clusters is None:
            raise ValueError("n_clusters must be specified.")

        self.clustering_molecules = Path(clustering_molecules)
        self.clustering_output = Path(clustering_output)
        self.n_clusters = n_clusters

    def validate(self) -> None:
        # path
        if not self.clustering_molecules.exists():
            raise FileNotFoundError(
                f"Molecules file '{self.clustering_molecules}' does not exist."
            )
        # values
        if self.n_clusters <= 0:
            raise ValueError("Number of clusters must be a positive integer.")
       
        
class DecomposeSettings:
    ligand_path: Path
    annotated_path: Path

    def __init__(
        self,
        ligand_path: StrPath,
        annotated_path: StrPath,
        **kwargs,
    ) -> None:
        self.ligand_path = Path(ligand_path)
        self.annotated_path = Path(annotated_path)

    def validate(self) -> None:
        # path
        if not self.ligand_path.exists():
            raise FileNotFoundError(
                f"Ligand file '{self.ligand_path}' does not exist."
            )
        if not self.annotated_path.parent.exists():
            raise FileNotFoundError(
                f"Annotated path '{self.annotated_path.parent}' does not exist."
            )


@dataclass(slots=True)
class CoffeePrescSettings:
    # each module settings
    fragquery: FragQuerySettings
    retrieval: RetrievalSettings
    scoring: CmpdEvalSettings

    # input/output/intermediate files
    receptor_path: Path
    fragments_path: Path
    output_path: Path
    query_path: Path
    matched_path: Path
    grid_folder: Path
    log_path: Path | None

    verbosity: int

    def validate(self) -> None:
        self.fragquery.validate()
        self.retrieval.validate()
        self.scoring.validate()

        # input files
        if not self.receptor_path.exists():
            raise FileNotFoundError(
                f"Receptor file '{self.receptor_path}' does not exist."
            )
        if not self.fragments_path.exists():
            raise FileNotFoundError(
                f"Fragments file '{self.fragments_path}' does not exist."
            )

        if not self.verbosity in [0, 1, 2]:
            raise ValueError(
                f"Verbosity must be 0 (quiet), 1 (normal), or 2 (verbose): got {self.verbosity}."
            )

    @classmethod
    def parse(
        cls,
        docking_config: StrPath | None = None,
        receptor: StrPath | None = None,
        fragments: StrPath | None = None,
        output: StrPath | None = None,
        query: StrPath | None = None,
        grid: StrPath | None = None,
        matched: StrPath | None = None,
        log: StrPath | None = None,
        verbosity: int = 1,
        setting: StrPath | None = None,
        **kwargs,
    ) -> Self:
        data: dict[str, Any] = {}
        # load from toml file
        if setting is not None:
            with open(setting, "rb") as f:
                data = tomllib.load(f)
        # override with kwargs
        for k, v in kwargs.items():
            if v is not None:
                data[k] = v

        # merge explicit args and values from `data` with concise logic
        params: dict[str, Any] = {
            "docking_config": docking_config,
            "receptor": receptor,
            "fragments": fragments,
            "output": output,
            "query": query,
            "grid": grid,
            "matched": matched,
            "log": log,
            "verbosity": verbosity,
        }

        # if a value is missing in params but present in data, take it from data.
        for key in list(params.keys()):
            if params[key] is None and key in data:
                params[key] = data.pop(key)

        # required fields
        required = [
            "docking_config",
            "receptor",
            "fragments",
            "output",
        ]
        for key in required:
            if params.get(key) is None:
                raise ValueError(
                    f"{key} must be specified either in setting file or as argument."
                )

        # tempfile
        if params.get("query") is None:
            tmp_query = tempfile.NamedTemporaryFile(suffix=".csv", delete=False)
            params["query"] = tmp_query.name
            tmp_query.close()
        if params.get("matched") is None:
            tmp_matched = tempfile.NamedTemporaryFile(suffix=".csv", delete=False)
            params["matched"] = tmp_matched.name
            tmp_matched.close()
        # temp folder
        if params.get("grid") is None:
            tmp_grid = tempfile.TemporaryDirectory()
            params["grid"] = tmp_grid.name
        # default name
        if params.get("log") is None:
            params["log"] = Path("coffeepresc.log")
        if params.get("verbosity") is None:
            params["verbosity"] = 1

        return cls(
            fragquery=FragQuerySettings(
                docking_config=params["docking_config"], **data
            ),
            retrieval=RetrievalSettings(**data),
            scoring=CmpdEvalSettings(**data),
            receptor_path=Path(params["receptor"]),
            fragments_path=Path(params["fragments"]),
            output_path=Path(params["output"]),
            query_path=Path(params["query"]),
            matched_path=Path(params["matched"]),
            grid_folder=Path(params["grid"]),
            log_path=Path(params["log"]),
            verbosity=params["verbosity"],
        )
