"""Load and validate BRB pipeline configuration."""

import csv
from pathlib import Path

try:
    import yaml
except ImportError:
    raise ImportError("PyYAML is required. Install with: pip install pyyaml")


class Sample:
    """Represents a single sample from the mapping file."""

    def __init__(self, sample_name, rt_barcode, well_name=None):
        self.sample_name = sample_name
        self.rt_barcode = rt_barcode
        self.well_name = well_name

    def __repr__(self):
        well_str = f", well={self.well_name}" if self.well_name else ""
        return f"Sample({self.sample_name}, barcode={self.rt_barcode}{well_str})"


class Config:
    """Pipeline configuration object."""

    def __init__(self, yaml_dict, mapping_path):
        """Initialize config from YAML dict and mapping file."""
        self.yaml_dict = yaml_dict
        self.mapping_path = Path(mapping_path)
        self.samples = []
        self._validate()
        self._load_mapping()

    def _validate(self):
        """Validate required YAML fields."""
        required_sections = ["project", "reference", "processing", "reads", "multiqc"]
        for section in required_sections:
            if section not in self.yaml_dict:
                raise ValueError(f"Missing required section: {section}")

        # Validate project
        project = self.yaml_dict["project"]
        if "name" not in project or "output_dir" not in project:
            raise ValueError("project section must have 'name' and 'output_dir'")

        # Validate reference
        reference = self.yaml_dict["reference"]
        required_ref = ["species", "star_index", "gtf", "bed"]
        for field in required_ref:
            if field not in reference:
                raise ValueError(f"reference section must have '{field}'")

        # Validate processing
        processing = self.yaml_dict["processing"]
        if "demultiplex" not in processing or "remove_intermediate" not in processing:
            raise ValueError(
                "processing section must have 'demultiplex' and 'remove_intermediate'"
            )

        # Validate reads
        reads = self.yaml_dict["reads"]
        if "read1" not in reads or "read2" not in reads:
            raise ValueError("reads section must have 'read1' and 'read2'")

        # Validate multiqc
        multiqc = self.yaml_dict["multiqc"]
        if "library_name" not in multiqc:
            raise ValueError("multiqc section must have 'library_name'")

    def _load_mapping(self):
        """Load and parse the TSV mapping file."""
        if not self.mapping_path.exists():
            raise FileNotFoundError(f"Mapping file not found: {self.mapping_path}")

        with open(self.mapping_path, "r") as f:
            reader = csv.DictReader(f, delimiter="\t")
            required_fields = {"sample_name", "rt_barcode"}

            if reader.fieldnames is None:
                raise ValueError("Mapping file is empty")

            missing_fields = required_fields - set(reader.fieldnames)
            if missing_fields:
                raise ValueError(
                    f"Mapping file must have columns: sample_name, rt_barcode. "
                    f"Missing: {', '.join(missing_fields)}"
                )

            for row in reader:
                sample = Sample(
                    sample_name=row["sample_name"],
                    rt_barcode=row["rt_barcode"],
                    well_name=row.get("well_name"),
                )
                self.samples.append(sample)

        if not self.samples:
            raise ValueError("Mapping file contains no samples")

    # Convenience properties
    @property
    def project_name(self):
        return self.yaml_dict["project"]["name"]

    @property
    def output_dir(self):
        return Path(self.yaml_dict["project"]["output_dir"])

    @property
    def species(self):
        return self.yaml_dict["reference"]["species"]

    @property
    def star_index(self):
        return self.yaml_dict["reference"]["star_index"]

    @property
    def gtf(self):
        return self.yaml_dict["reference"]["gtf"]

    @property
    def bed(self):
        return self.yaml_dict["reference"]["bed"]

    @property
    def demultiplex(self):
        return self.yaml_dict["processing"]["demultiplex"]

    @property
    def remove_intermediate(self):
        return self.yaml_dict["processing"]["remove_intermediate"]

    @property
    def read1(self):
        return self.yaml_dict["reads"]["read1"]

    @property
    def read2(self):
        return self.yaml_dict["reads"]["read2"]

    @property
    def library_name(self):
        return self.yaml_dict["multiqc"]["library_name"]

    @property
    def multiqc_template(self):
        template = self.yaml_dict["multiqc"].get("config_template")
        return template if template else None

    @property
    def cpus_per_task(self):
        slurm = self.yaml_dict.get("slurm", {})
        return slurm.get("cpus_per_task", 4)

    @property
    def mem(self):
        slurm = self.yaml_dict.get("slurm", {})
        return slurm.get("mem", 75000)


def load_config(config_path):
    """
    Load and validate BRB pipeline configuration.

    Parameters
    ----------
    config_path : str or Path
        Path to the YAML config file.

    Returns
    -------
    Config
        A validated Config object.

    Raises
    ------
    FileNotFoundError
        If config file or mapping file does not exist.
    ValueError
        If config is missing required fields or invalid.
    """
    config_path = Path(config_path)
    if not config_path.exists():
        raise FileNotFoundError(f"Config file not found: {config_path}")

    with open(config_path, "r") as f:
        yaml_dict = yaml.safe_load(f)

    if yaml_dict is None:
        raise ValueError("Config file is empty or not valid YAML")

    # Get the mapping file path (can be relative to config file directory)
    if "mapping_file" not in yaml_dict:
        raise ValueError("Config must specify 'mapping_file'")

    mapping_file = yaml_dict["mapping_file"]
    mapping_path = Path(mapping_file)

    # If mapping file is relative, resolve it relative to the current working directory first.
    # If that path does not exist, fall back to resolving relative to the config file location.
    if not mapping_path.is_absolute():
        cwd_path = Path.cwd() / mapping_path
        if cwd_path.exists():
            mapping_path = cwd_path
        else:
            mapping_path = config_path.parent / mapping_path

    return Config(yaml_dict, mapping_path)
