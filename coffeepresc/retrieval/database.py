import os
import subprocess
from pathlib import Path

from coffeepresc.common.typing import StrPath

_PATH_FBDB_CLI = f"{os.path.dirname(__file__)}/FBDBv2/cli.py"


class Database:
    def __init__(
        self,
        storage_path: StrPath | None = None,
        verbosity: int = 1,
        log_path: StrPath | None = None,
    ):
        if storage_path is None:
            storage_path = "data.h5"
        self.__storage_path = Path(storage_path)
        self.set_log(log_path)
        self.set_verbosity(verbosity)
        # TODO: version check

    @property
    def storage_path(self) -> Path:
        return self.__storage_path

    def set_log(self, log_path: StrPath | None) -> None:
        if log_path is not None:
            self.__cmd_l = ["--log_path", str(log_path)]

    def set_verbosity(self, verbosity: int) -> None:
        if verbosity == 0:
            self.__cmd_v = ["--quiet"]
        elif verbosity == 1:
            self.__cmd_v = []
        elif verbosity == 2:
            self.__cmd_v = ["--verbose"]
        else:
            raise ValueError(
                f"Verbosity must be 0 (quiet), 1 (normal), or 2 (verbose):",
                f"got {verbosity}.",
            )

    def _run(self, cmd: list[str]) -> None:
        subprocess.run(
            cmd + self.__cmd_v + self.__cmd_l,
            check=True,
        )

    def create(
        self,
        sdf_path: StrPath,
        fragment_path: StrPath,
        chunk_size: int = 100000,
        dropped_path: StrPath | None = None,
    ) -> None:
        cmd = [
            "python",
            _PATH_FBDB_CLI,
            "create",
            str(sdf_path),
            str(fragment_path),
            "--storage_path",
            str(self.storage_path),
            "--chunk_size",
            str(chunk_size),
        ]
        if dropped_path is not None:
            cmd.extend(["--dropped_path", str(dropped_path)])
        self._run(cmd)

    def search(
        self,
        query_path: StrPath,
        output_path: StrPath,
        similarity_th: float = 0.45,
    ) -> None:
        cmd = [
            "python",
            _PATH_FBDB_CLI,
            "search",
            str(query_path),
            str(output_path),
            "--storage_path",
            str(self.storage_path),
            "--similarity_th",
            str(similarity_th),
        ]
        self._run(cmd)
