# /// script
# requires-python = ">=3.12"
# dependencies = [
#   "wandb[media]",
#   "loguru",
# ]
# ///

import argparse
import json
import subprocess
import tempfile
import tomllib
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from loguru import logger

import wandb

ENTITY = "jspaezp"
PROJECT = "timsseek"


@dataclass
class TimsseekRunner:
    fasta_file_location: Path
    speclib_location: Path
    raw_file_location: Path
    config_dict: dict[str, Any] | None = None

    def build_speclib(self):
        if self.speclib_location.exists():
            logger.info("Skipping speclib build bc already exists")
            return

        logger.info("Building speclib")
        args = [
            "uv",
            "run",
            "speclib_build_fasta",
            "--fasta_file",
            str(self.fasta_file_location),
            "--decoy_strategy",
            "REVERSE",
            "--max_ions",
            "10",
            "--outfile",
            str(self.speclib_location),
            "--model",
            "onnx",
        ]
        res = subprocess.run(args, check=True)
        return res

    def setup_run(self):
        if self.config_dict is None:
            self.config_dict = self.default_timsseek_config()

        logger.info("Building release versions")
        subprocess.run(["cargo", "b", "--release"], check=True)

    def loggable_config_dict(self) -> dict[str, Any]:
        out = {}
        out.update(self.config_dict)
        out["raw_file"] = self.raw_file_location.name
        out["speclib"] = self.speclib_location.name
        return out

    def run(
        self,
        output_loc: Path | None = None,
        wandb_kwargs: dict | None = None,
    ):
        self.setup_run()

        with tempfile.TemporaryDirectory() as temp_dir:
            tmpdir = Path(temp_dir)
            if output_loc is None:
                output_loc = tmpdir

            results_path = output_loc / "res"
            summary_dir = output_loc / "summ"
            results_path.mkdir(exist_ok=True, parents=True)
            summary_dir.mkdir(exist_ok=True, parents=True)
            config_path = results_path / "config.json"

            with open(config_path, "w") as f:
                f.write(json.dumps(self.config_dict))

            with wandb_context(
                self.loggable_config_dict(),
                wandb_kwargs=wandb_kwargs,
            ) as wandb_experiment:
                self._run(
                    config_path=config_path,
                    speclib_path=self.speclib_location,
                    output_path=results_path,
                    raw_file=self.raw_file_location,
                )
                # Results are now in a subdirectory named after the raw file
                raw_file_stem = self.raw_file_location.stem
                file_results_dir = results_path / raw_file_stem
                self._rescore(results_dir=file_results_dir, summary_dir=summary_dir)
                self.log_results(wandb_experiment, output_loc, raw_file_stem)

    @staticmethod
    def _run(config_path, speclib_path, output_path, raw_file):
        if not raw_file.exists():
            raise FileNotFoundError(f"Raw file {raw_file} does not exist")
        if not speclib_path.exists():
            raise FileNotFoundError(f"Speclib file {speclib_path} does not exist")
        if not config_path.exists():
            raise FileNotFoundError(f"Config file {config_path} does not exist")

        logger.info(f"Running timsseek on {raw_file.name}")
        args = [
            "cargo",
            "run",
            "--release",
            "--bin",
            "timsseek",
            "--",
            "--overwrite",
            "--config",
            str(config_path),
            "--speclib-file",
            str(speclib_path),
            "--output-dir",
            str(output_path),
            "--dotd-files",
            str(raw_file),
        ]
        logger.info(f"Running command: {' '.join(args)}")
        stdout_file = output_path / "timsseek_stdout.log"
        stderr_file = output_path / "timsseek_stderr.log"

        logger.info(f"Starting timsseek, logging to {stdout_file} and {stderr_file}")
        try:
            res = subprocess.run(
                args,
                stdout=open(stdout_file, "w"),
                stderr=open(stderr_file, "w"),
                check=True,
            )
        finally:
            # Log stdout and stderr
            logger.info(stdout_file.read_text())
            logger.error(stderr_file.read_text())
        logger.info(f"Timsseek completed with return code {res.returncode}")
        return res

    @staticmethod
    def _rescore(results_dir: Path, summary_dir: Path):
        args = [
            "uv",
            "run",
            "python",
            "-m",
            "timsseek_rescore",
            "--results_dir",
            str(results_dir),
            "--output_dir",
            str(summary_dir),
        ]
        stdout_file = summary_dir / "timsseek_stdout.log"
        stderr_file = summary_dir / "timsseek_stderr.log"

        logger.info(f"Starting rescoring, logging to {stdout_file} and {stderr_file}")
        try:
            res = subprocess.run(
                args,
                stdout=open(stdout_file, "w"),
                stderr=open(stderr_file, "w"),
                check=True,
            )
        finally:
            # Log stdout and stderr
            logger.info(stdout_file.read_text())
            logger.error(stderr_file.read_text())

        logger.info(f"Rescoring completed with return code {res.returncode}")
        return res

    def log_results(self, wandb_experiment, results_loc, raw_file_stem):
        metrics = self.crunch_metrics(results_loc, raw_file_stem)
        with open("latest_metrics.json", "w") as f:
            serializable_metrics = {
                k: v
                for k, v in metrics.items()
                if isinstance(v, (int, float, str, bool, list, dict))
            }
            assert serializable_metrics
            json.dump(serializable_metrics, f, indent=4)
        wandb_experiment.log(metrics)

    @staticmethod
    def default_timsseek_config():
        config = {
            "analysis": {
                "chunk_size": 20000,
                "tolerance": {
                    "ms": {"ppm": [15.0, 15.0]},
                    "mobility": {"percent": [10.0, 10.0]},
                    "quad": {"absolute": [0.1, 0.1]},
                    "rt": {"minutes": [5, 5]},
                },
            }
        }
        return config

    @staticmethod
    def crunch_metrics(output_dir: Path, raw_file_stem: str) -> dict[str, Any]:
        metrics = {}
        xgboost_images = [
            ("variable_importance_plot", "importances.png"),
            ("mass_error_plot", "mass_error_rt_1pct.png"),
            ("mobility_error_plot", "mobility_error_rt_1pct.png"),
            ("mz_mobility_plot", "mz_mobility_1pct.png"),
            ("predicted_rt_obs_plot", "predicted_rt_obs_rt_1pct.png"),
            ("mass_error_mz_1pct_plot", "mass_error_mz_1pct.png"),
        ]

        for metric_name, file_name in xgboost_images:
            img_path = output_dir / "summ" / "xgboost" / file_name
            if img_path.exists():
                metrics[metric_name] = wandb.Image(img_path)
            else:
                logger.warning(f"Image {img_path} does not exist, skipping")

        report_toml = output_dir / "summ" / "xgboost" / "report.toml"
        with open(report_toml, "rb") as f:
            report = tomllib.load(f)
        metrics.update(report["report"])

        performance_report_path = (
            output_dir / "res" / raw_file_stem / "performance_report.json"
        )
        if performance_report_path.exists():
            with open(performance_report_path, "r") as f:
                performance_report = json.load(f)
            metrics.update(performance_report)
        else:
            logger.warning(
                f"Performance report {performance_report_path} does not exist"
            )
        return metrics


@contextmanager
def wandb_context(config_dict: dict[str, Any], wandb_kwargs=None):
    if wandb_kwargs is None:
        wandb_kwargs = {}
    # Start a new wandb run to track this script.
    run = wandb.init(
        # Set the wandb entity where your project will be logged (generally your team name).
        entity=ENTITY,
        # Set the wandb project where this run will be logged.
        project=PROJECT,
        # Track hyperparameters and run metadata.
        config=config_dict,
        **wandb_kwargs,
    )
    try:
        yield run
    except KeyboardInterrupt as e:
        logger.warning("Keyboard interrupt, finishing wandb run")
        run.finish(1)
        raise e
    finally:
        run.finish()


def main(wandb_kwargs: dict | None = None):

    fasta_file = Path.home() / "fasta/hela_gt20peps.fasta"
    speclib_path = Path.home() / "fasta/asdad.msgpack.zstd"

    prefix = Path.home() / "data/decompressed_timstof/"
    dotd_files = [
        # prefix / "MSR28858_EXP80_Plate3_G08_DMSO_DIA_S5-G8_1_7079.d",
        # prefix / "MSR28893_EXP80_Plate4_B07_DMSO_DIA_S6-B7_1_7115.d",
        # prefix / "250225_Desnaux_200ng_Hela_ICC_on_DIA.d",
        prefix / "250225_Desnaux_200ng_Hela_ICC_off_DIA.d",
    ]

    for file in dotd_files:
        runner = TimsseekRunner(
            fasta_file_location=fasta_file,
            speclib_location=speclib_path,
            raw_file_location=file,
        )
        runner.build_speclib()
        runner.run(wandb_kwargs=wandb_kwargs)


def build_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--notes",
        type=str,
        help="The notes to add to the wandb run",
    )
    return parser


if __name__ == "__main__":
    parser = build_parser()
    args, unkargs = parser.parse_known_args()
    if unkargs:
        raise ValueError(f"Unknown arguments: {unkargs}")

    wandb_kwargs = {}

    if args.notes is not None:
        wandb_kwargs["notes"] = args.notes

    main(wandb_kwargs=wandb_kwargs)
