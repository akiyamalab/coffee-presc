import os
import subprocess
import tempfile
from pathlib import Path

from coffeepresc.common.typing import StrPath
from coffeepresc.models import FragQuerySettings

__PATH_BIN_DIR = f"{os.path.dirname(__file__)}/bin"
__PATH_ATOMGRID_GEN = f"{__PATH_BIN_DIR}/atomgrid-gen"
__PATH_FRAGMENT_QUERY = f"{__PATH_BIN_DIR}/fragment-query"


def _generate_atomgrid(
    settings: FragQuerySettings,
    receptor_path: StrPath,
    grid_folder: StrPath,
    log_path: StrPath | None = None,
    verbosity: int = 1,
) -> None:
    """
    execute atomgrid-gen based on the given config file.
    """
    with tempfile.NamedTemporaryFile(mode="w", delete=False) as fout:
        fout.write(settings.dumps_docking_config())
        fout.flush()
        config_path = Path(fout.name)
        if log_path is None:
            log_path = "/dev/null"
        cmd = [
            __PATH_ATOMGRID_GEN,
            str(config_path),
            "-r", str(receptor_path),
            "-g", str(grid_folder),
            "--log", str(log_path),
            "--verbosity", str(verbosity),
        ]
        subprocess.run(
            cmd,
            check=True,
        )


def _generate_query(
    settings: FragQuerySettings,
    receptor_path: StrPath,
    fragments_path: StrPath,
    output_path: StrPath,
    grid_folder: StrPath,
    log_path: StrPath | None = None,
    verbosity: int = 1,
) -> None:
    """
    execute fragment-query based on the given config file.
    """
    with tempfile.NamedTemporaryFile(mode="w", delete=False) as fout:
        fout.write(settings.dumps_docking_config())
        fout.flush()
        config_path = Path(fout.name)
        if log_path is None:
            log_path = "/dev/null"
        cmd = [
            __PATH_FRAGMENT_QUERY,
            str(config_path),
            "-o", str(output_path),
            "-r", str(receptor_path),
            "-g", str(grid_folder),
            "--log", str(log_path),
            "-f", str(fragments_path),
            "--verbosity", str(verbosity),
        ]
        subprocess.run(
            cmd,
            check=True,
        )


class FragQuery:
    def __init__(self, settings: FragQuerySettings):
        self.__settings = settings

    @property
    def settings(self) -> FragQuerySettings:
        return self.__settings

    def run(
        self,
        receptor_path: StrPath,
        fragments_path: StrPath,
        output_path: StrPath,
        grid_folder: StrPath | None = None,
        log_path: StrPath | None = None,
        verbosity: int = 1,
    ) -> None:
        temp_folder = None
        if grid_folder is None:
            temp_folder = tempfile.TemporaryDirectory()

        # grid_folder may be a TemporaryDirectory object or a path-like string
        folder_path = temp_folder.name if temp_folder is not None else str(grid_folder)

        try:
            _generate_atomgrid(
                self.settings,
                receptor_path,
                folder_path,
                log_path,
                verbosity=verbosity,
            )
            _generate_query(
                self.settings,
                receptor_path,
                fragments_path,
                output_path,
                folder_path,
                log_path,
                verbosity=verbosity,
            )
        except Exception as e:
            raise e
        finally:
            if temp_folder is not None:
                try:
                    temp_folder.cleanup()
                except Exception:
                    pass
