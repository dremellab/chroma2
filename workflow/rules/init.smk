###################################################################################
# Function definitions
###################################################################################

def _is_true(variable):
    if variable == True or variable == "True" or variable == "TRUE":
        return True
    else:
        return False

###################################################################################

def _convert_to_int(variable):
    if variable:
        return 1  # True
    if not variable:
        return 0  # False
    return -1  # Unknown

###################################################################################

def append_files_in_list(flist, ofile):
    if not os.path.exists(ofile):
        print("FILE %s does not exist! Creating it!"%(ofile),flush=True)
        with open(ofile, "w") as outfile:
            for fname in flist:
                with open(fname) as infile:
                    for line in infile:
                        outfile.write(line)
    return True

###################################################################################

def _logdir(rule_name):
    return join(RESULTSDIR, "logs", rule_name)

###################################################################################

def get_peorse(wildcards):
    peorse = SAMPLESDF.loc[SAMPLESDF['sampleName'] == wildcards.sample, 'PEorSE'].values[0]
    return peorse

###################################################################################

def get_fastqs(wildcards):
    d = dict()
    peorse = SAMPLESDF.loc[SAMPLESDF['sampleName'] == wildcards.sample, 'PEorSE'].values[0]
    # print(f"peorse: {peorse}")
    if peorse == "PE":
        d["R1"] = SAMPLESDF.loc[SAMPLESDF['sampleName'] == wildcards.sample, 'path_to_R1_fastq'].values[0]
        d["R2"] = SAMPLESDF.loc[SAMPLESDF['sampleName'] == wildcards.sample, 'path_to_R2_fastq'].values[0]
    else:
        d["R1"] = SAMPLESDF.loc[SAMPLESDF['sampleName'] == wildcards.sample, 'path_to_R1_fastq'].values[0]
        d["R2"] = DUMMYFILE

    # print(f"sample: {wildcards.sample}")
    # print(f"R1: ##{d['R1']}##")
    # print(f"R2: ##{d['R2']}##")
    return d

###################################################################################

def _get_threads(rule_name, profile_config):
    """
    Return threads for a rule from profile_config.
    Falls back to default if not defined.
    """
    if (
        "set-resources" in profile_config
        and rule_name in profile_config["set-resources"]
        and "threads" in profile_config["set-resources"][rule_name]
    ):
        return profile_config["set-resources"][rule_name]["threads"]
    return profile_config["default-resources"]["threads"]

###################################################################################

def _get_mem_mb(rule_name, profile_config):
    """
    Return mem_mb for a rule from profile_config.
    Falls back to default if not defined.
    """
    if (
        "set-resources" in profile_config
        and rule_name in profile_config["set-resources"]
        and "mem_mb" in profile_config["set-resources"][rule_name]
    ):
        return profile_config["set-resources"][rule_name]["mem_mb"]
    return profile_config["default-resources"].get("mem_mb", 4096)


def getmemG(rule_name):
    # Return memory in MB for tools that expect MB (despite name).
    return int(_get_mem_mb(rule_name, profile_config))

###################################################################################
###################################################################################

# import yaml
# from pathlib import Path

# Locate your cluster profile (relative to workdir or absolute path)
profile_path = join(Path(os.environ["PROFILE"]) , "config.yaml")
with open(profile_path) as f:
    profile_config = yaml.safe_load(f)

# Now profile_config is a normal dict
# pprint(profile_config)
# sys.exit(1)


# print("printing config...")
WORKDIR = os.getcwd()
print(WORKDIR)
configfilepath = join(WORKDIR, "config.yaml") # this is workflow config .. not to be confused with snakemake cluster profile above
try:
    with open(configfilepath, "r") as f:
        config = yaml.safe_load(f)
except Exception as e:
    print(f"❌ File does not exist: {configfilepath}")
    print(f"❌ Error opening config file: {e}")
    sys.exit(1)
_is_unlock = "--unlock" in sys.argv
if not _is_unlock:
    print("Snakemake working directory:", WORKDIR)
# print(config)
# print("end of config")

# resource absolute path
# WORKDIR = config["workdir"]
TEMPDIR = config["tempdir"]
SCRIPTS_DIR = config["scriptsdir"]
RESOURCES_DIR = config["resourcesdir"]
FASTAS_GTFS_DIR = config["fastas_gtfs_dir"]

REF_DIR = join(WORKDIR, "ref")
if not os.path.exists(REF_DIR):
    os.mkdir(REF_DIR)

# strip trailing slashes if any
for varname in [
    "WORKDIR", "SCRIPTS_DIR", "RESOURCES_DIR", "FASTAS_GTFS_DIR",
    "REF_DIR", "TEMPDIR"
]:
    globals()[varname] = globals()[varname].rstrip(r"\/")

HOST = config["host"].strip()  # hg38 or mm39
ADDITIVES = config["additives"].strip()  # ERCC and/or BAC16Insert
ADDITIVES = ADDITIVES.replace(" ", "")
VIRUSES = config["viruses"].strip()
VIRUSES = VIRUSES.replace(" ", "")

if HOST != "":
    if ADDITIVES != "":
        HOST_ADDITIVES = HOST + "," + ADDITIVES
    else:
        HOST_ADDITIVES = HOST

    if VIRUSES != "":
        HOST_ADDITIVES_VIRUSES = HOST_ADDITIVES + "," + VIRUSES
        HOST_VIRUSES = HOST + "," + VIRUSES
    else:
        HOST_VIRUSES = HOST
        HOST_ADDITIVES_VIRUSES = HOST_ADDITIVES
else:
    if ADDITIVES != "":
        HOST_ADDITIVES = ADDITIVES
        HOST_ADDITIVES_VIRUSES = ADDITIVES
        HOST_VIRUSES = ""
    else:
        HOST_ADDITIVES = ""
        HOST_ADDITIVES_VIRUSES = ""
        HOST_VIRUSES = ""
    if VIRUSES != "":
        HOST_ADDITIVES_VIRUSES = VIRUSES
        HOST_VIRUSES = VIRUSES
    else:
        raise ValueError("Both host and viruses are not set. Please set at least one of them.")

REFERENCE_GTF_BY_HOST = config.get("reference_gtf", {})
TRNAS_GTF_BY_HOST = config.get("trnas_gtf", {})
CHRR_GTF_BY_HOST = config.get("chrr_gtf", {})
INCLUDE_TRNAS_GTF_IN_REF = _is_true(config.get("include_trnas_gtf_in_ref", True))


def _host_annotation_aliases(host_name):
    aliases = []
    if host_name:
        aliases.append(host_name)
        if "_" in host_name:
            aliases.append(host_name.split("_", 1)[0])
    return aliases


def _resolve_optional_host_annotation(mapping, host_name, fallback_templates):
    if host_name == "":
        return ""

    for alias in _host_annotation_aliases(host_name):
        mapped_name = mapping.get(alias)
        if mapped_name:
            mapped_path = (
                mapped_name if os.path.isabs(mapped_name) else join(FASTAS_GTFS_DIR, mapped_name)
            )
            if os.path.exists(mapped_path):
                return mapped_path

    for alias in _host_annotation_aliases(host_name):
        for template in fallback_templates:
            candidate_name = template.format(host=alias)
            candidate_path = (
                candidate_name
                if os.path.isabs(candidate_name)
                else join(FASTAS_GTFS_DIR, candidate_name)
            )
            if os.path.exists(candidate_path):
                return candidate_path

    return ""


if HOST != "":
    TRNAS_GTF = _resolve_optional_host_annotation(
        TRNAS_GTF_BY_HOST,
        HOST,
        ["{host}.tRNAs.{host}chroms.gtf"],
    )
    CHRR_GTF = _resolve_optional_host_annotation(
        CHRR_GTF_BY_HOST,
        HOST,
        ["{host}.chrR.gtf"],
    )
else:
    TRNAS_GTF = ""
    CHRR_GTF = ""

# print (f"Using TRNAS_GTF: {TRNAS_GTF}")
# exit()

HOST_ADDITIVES_VIRUSES = HOST_ADDITIVES_VIRUSES.split(",")
HOST_VIRUSES = HOST_VIRUSES.split(",")
HOST_LIST = [h for h in HOST.split(",") if h] if HOST != "" else []
ADDITIVE_LIST = [a for a in ADDITIVES.split(",") if a] if ADDITIVES != "" else []
VIRUS_LIST = [v for v in VIRUSES.split(",") if v] if VIRUSES != "" else []
FASTAS = [join(FASTAS_GTFS_DIR, f + ".fa") for f in HOST_ADDITIVES_VIRUSES]
REGIONS = [join(FASTAS_GTFS_DIR, f + ".fa.regions") for f in HOST_ADDITIVES_VIRUSES]
if HOST != "":
    REGIONS_HOST = [join(FASTAS_GTFS_DIR, f + ".fa.regions") for f in HOST.split(",")]
else:
    REGIONS_HOST = []
if VIRUSES != "":
    REGIONS_VIRUSES = [join(FASTAS_GTFS_DIR, f + ".fa.regions") for f in VIRUSES.split(",")]
else:
    REGIONS_VIRUSES = []
def _get_genome_gtf(genome_name):
    # First check if explicitly defined in reference_gtf config
    if genome_name in REFERENCE_GTF_BY_HOST:
        gtf_file = REFERENCE_GTF_BY_HOST[genome_name]
        return gtf_file if os.path.isabs(gtf_file) else join(FASTAS_GTFS_DIR, gtf_file)
    # Fall back to default {genome}.gtf
    return join(FASTAS_GTFS_DIR, genome_name + ".gtf")

GTFS = [_get_genome_gtf(f) for f in HOST_ADDITIVES_VIRUSES]
if CHRR_GTF != "":
    GTFS.append(CHRR_GTF)
if INCLUDE_TRNAS_GTF_IN_REF and TRNAS_GTF != "":
    GTFS.append(TRNAS_GTF)
FASTAS_REGIONS_GTFS = FASTAS.copy()
FASTAS_REGIONS_GTFS.extend(REGIONS)
FASTAS_REGIONS_GTFS.extend(GTFS)
EGS = join(FASTAS_GTFS_DIR, "effectiveGenomeSizes.tsv")

if not _is_unlock:
    print("FASTAS_REGIONS_GTFS:")
    for item in FASTAS_REGIONS_GTFS:
        print(f"  - {item}")

REF_FA = join(REF_DIR, "ref.fa")
REF_REGIONS = join(REF_DIR, "ref.fa.regions")
REF_REGIONS_HOST = join(REF_DIR, "ref.fa.regions.host")
REF_REGIONS_VIRUSES = join(REF_DIR, "ref.fa.regions.viruses")
REF_REGIONS_HOST_VIRUSES = join(REF_DIR, "ref.fa.regions.host_viruses")
REF_GTF = join(REF_DIR, "ref.gtf")
append_files_in_list(FASTAS, REF_FA)
append_files_in_list(REGIONS, REF_REGIONS)

###################################################################################################
# check if sequence IDs are unique for unique genome names

if not _is_unlock:
    print("Validating ref.regions file for unique genome names and sequence IDs...")
    input_file = REF_REGIONS

    seqid_to_genomes = defaultdict(set)
    seen_genomes = set()
    duplicate_genomes = set()

    with open(input_file) as f:
        for line_number, line in enumerate(f, 1):
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            genome = parts[0]
            seq_ids = parts[1:]

            # Check for repeated genome names
            if genome in seen_genomes:
                duplicate_genomes.add(genome)
            seen_genomes.add(genome)

            # Track which genomes each sequence ID appears under
            for seq in seq_ids:
                seqid_to_genomes[seq].add(genome)

    # Detect conflicts in sequence IDs
    conflicts = {seq: genomes for seq, genomes in seqid_to_genomes.items() if len(genomes) > 1}

    # Report issues
    if duplicate_genomes:
        print("\n❌ The following genome names are repeated:")
        for g in sorted(duplicate_genomes):
            print(f"  {g}")

    if conflicts:
        print("\n❌ The following sequence IDs are assigned to multiple genomes:")
        for seq, genomes in sorted(conflicts.items()):
            print(f"  {seq}: {', '.join(sorted(genomes))}")

    if not duplicate_genomes and not conflicts:
        print("\n✅ All genome names are unique and sequence IDs are uniquely assigned.")
    else:
        sys.exit(1)

###################################################################################################


append_files_in_list(REGIONS_HOST, REF_REGIONS_HOST)
append_files_in_list(REGIONS_HOST + REGIONS_VIRUSES, REF_REGIONS_HOST_VIRUSES)
append_files_in_list(REGIONS_VIRUSES, REF_REGIONS_VIRUSES)

if not os.path.exists(REF_GTF):

    ###################################################################################################
    # check if gene_id and gene_name are unique across GTF files

    print("Validating GTF files for unique gene_id and gene_name...")
    # Extract gene_id and gene_name from attribute string
    def parse_attributes(attr_string):
        gene_id = gene_name = None
        matches = re.findall(r'(\S+)\s+"([^"]+)"', attr_string)
        for key, val in matches:
            if key == "gene_id":
                gene_id = val
            elif key == "gene_name":
                gene_name = val
        return gene_id, gene_name

    # Dicts to track where each ID was found
    gene_id_to_file = {}
    gene_name_to_file = {}

    # Read GTF files from command line arguments
    gtf_files = GTFS

    for gtf_file in gtf_files:
        with open(gtf_file) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                cols = line.strip().split("\t")
                if len(cols) < 9 or cols[2] != "gene":
                    continue
                attr_string = cols[8]
                gene_id, gene_name = parse_attributes(attr_string)

                if gene_id:
                    if gene_id in gene_id_to_file and gene_id_to_file[gene_id] != gtf_file:
                        print(f"❌ gene_id '{gene_id}' found in both '{gene_id_to_file[gene_id]}' and '{gtf_file}'")
                    gene_id_to_file[gene_id] = gtf_file

                if gene_name:
                    if gene_name in gene_name_to_file and gene_name_to_file[gene_name] != gtf_file:
                        print(f"❌ gene_name '{gene_name}' found in both '{gene_name_to_file[gene_name]}' and '{gtf_file}'")
                    gene_name_to_file[gene_name] = gtf_file

    print("✅ Done checking gene_id and gene_name uniqueness across GTF files.")

###################################################################################################

append_files_in_list(GTFS, REF_GTF)

# read in the samplesheet
# Step 1: Read the tab-delimited file with headers
MANIFEST_FILE = config.get("samples")
SAMPLESDF = pd.read_csv(MANIFEST_FILE, sep="\t", dtype=str).fillna("")

required_columns = [
    "sampleName",
    "groupName",
    "batch",
    "path_to_R1_fastq",
    "path_to_R2_fastq",
    "role",
    "target",
    "host_input_pool",
    "virus_input_pool",
]
# Check if all required columns are present
missing_columns = [col for col in required_columns if col not in SAMPLESDF.columns]
if missing_columns:
    print("Headers in the samplesheet:", [f'"{header}"' for header in SAMPLESDF.columns])
    raise ValueError(f"Missing required columns: {', '.join(missing_columns)}")

# Step 2: Confirm that sampleNames are unique
dup_mask = SAMPLESDF['sampleName'].duplicated(keep=False)
if dup_mask.any():
    dup_names = sorted(SAMPLESDF.loc[dup_mask, 'sampleName'].astype(str).tolist())
    raise ValueError(f"Duplicate sampleNames found: {dup_names}")

# SAMPLESDF.set_index(["sampleName"], inplace=False)
SAMPLES = list(SAMPLESDF["sampleName"])

# Step 3: Ensure each sampleName has a non-empty groupName
if (SAMPLESDF['groupName'].str.strip() == "").any():
    raise ValueError("Some sampleNames have empty groupName!")

# Step 3b: Normalize batch column (empty -> single batch)
if (SAMPLESDF['batch'].str.strip() == "").any():
    SAMPLESDF['batch'] = SAMPLESDF['batch'].replace("", "batch1")

# Step 3c: Validate role and target columns
valid_roles = {"case", "control"}
valid_targets = {"host", "virus", "both"}
if not SAMPLESDF['role'].isin(valid_roles).all():
    bad = sorted(SAMPLESDF.loc[~SAMPLESDF['role'].isin(valid_roles), 'role'].unique())
    raise ValueError(f"Invalid role values found: {bad}. Allowed: {sorted(valid_roles)}")
if not SAMPLESDF['target'].isin(valid_targets).all():
    bad = sorted(SAMPLESDF.loc[~SAMPLESDF['target'].isin(valid_targets), 'target'].unique())
    raise ValueError(f"Invalid target values found: {bad}. Allowed: {sorted(valid_targets)}")

# Step 4: Check if files in R1 and R2 paths exist and are readable
def check_file(path):
    return path != "" and os.path.isfile(path) and os.access(path, os.R_OK)

SAMPLESDF['R1_exists'] = SAMPLESDF['path_to_R1_fastq'].apply(check_file)
SAMPLESDF['R2_exists'] = SAMPLESDF['path_to_R2_fastq'].apply(check_file)

if not SAMPLESDF['R1_exists'].all():
    raise FileNotFoundError("Some R1 files are missing or not readable.")

# R2 may be missing for single-end, so we'll handle that below

# Step 5: Determine paired-end vs single-end
SAMPLESDF['PEorSE'] = SAMPLESDF['path_to_R2_fastq'].apply(lambda x: x.strip() != "")
SAMPLESDF['PEorSE'] = SAMPLESDF['PEorSE'].apply(lambda x: "PE" if x else "SE")

# Step 6: Create SAMPLENAME2GROUPNAME
SAMPLENAME2GROUPNAME = dict(zip(SAMPLESDF['sampleName'], SAMPLESDF['groupName']))

# Step 7: Create GROUPNAME2SAMPLENAME
from collections import defaultdict
GROUPNAME2SAMPLENAME = defaultdict(list)
for sample, group in zip(SAMPLESDF['sampleName'], SAMPLESDF['groupName']):
    GROUPNAME2SAMPLENAME[group].append(sample)

# Step 8: Create SAMPLENAMEISPE
SAMPLENAMEISPE = dict(zip(SAMPLESDF['sampleName'], SAMPLESDF['PEorSE']))

# Step 9: Case/control splits and input pools
CASE_SAMPLES = SAMPLESDF.loc[SAMPLESDF['role'] == "case", 'sampleName'].tolist()
CONTROL_SAMPLES = SAMPLESDF.loc[SAMPLESDF['role'] == "control", 'sampleName'].tolist()

HOST_INPUT_POOL_BY_SAMPLE = dict(zip(SAMPLESDF['sampleName'], SAMPLESDF['host_input_pool']))
VIRUS_INPUT_POOL_BY_SAMPLE = dict(zip(SAMPLESDF['sampleName'], SAMPLESDF['virus_input_pool']))

HOST_POOL_CONTROLS = defaultdict(list)
VIRUS_POOL_CONTROLS = defaultdict(list)
for _, row in SAMPLESDF.iterrows():
    if row['role'] != "control":
        continue
    if row['target'] in ["host", "both"] and row['host_input_pool'].strip() != "":
        HOST_POOL_CONTROLS[row['host_input_pool']].append(row['sampleName'])
    if row['target'] in ["virus", "both"] and row['virus_input_pool'].strip() != "":
        VIRUS_POOL_CONTROLS[row['virus_input_pool']].append(row['sampleName'])

HOST_INPUT_POOLS = sorted(HOST_POOL_CONTROLS.keys())
VIRUS_INPUT_POOLS = sorted(VIRUS_POOL_CONTROLS.keys())


def _opt_input(path):
    return [] if path == "" else path


def get_host_control_qname_bam(wildcards):
    pool = HOST_INPUT_POOL_BY_SAMPLE.get(wildcards.sample, "")
    if pool == "" or pool not in HOST_POOL_CONTROLS:
        return ""
    return join(RESULTSDIR, "inputs", "host", f"{pool}.qname.bam")


def get_host_control_filtered_bam(wildcards):
    pool = HOST_INPUT_POOL_BY_SAMPLE.get(wildcards.sample, "")
    if pool == "" or pool not in HOST_POOL_CONTROLS:
        return ""
    return join(RESULTSDIR, "inputs", "host", f"{pool}.filtered.bam")


def get_virus_control_qname_bam(wildcards):
    pool = VIRUS_INPUT_POOL_BY_SAMPLE.get(wildcards.sample, "")
    if pool == "" or pool not in VIRUS_POOL_CONTROLS:
        return ""
    return join(RESULTSDIR, "inputs", "virus", pool, f"{wildcards.virus}.qname.bam")


def get_virus_control_filtered_bam(wildcards):
    pool = VIRUS_INPUT_POOL_BY_SAMPLE.get(wildcards.sample, "")
    if pool == "" or pool not in VIRUS_POOL_CONTROLS:
        return ""
    return join(RESULTSDIR, "inputs", "virus", pool, f"{wildcards.virus}.filtered.bam")

DUMMYFILE = join(RESOURCES_DIR, "dummy")
RESULTSDIR = join(WORKDIR, "results")
if not os.path.exists(RESULTSDIR):
    os.mkdir(RESULTSDIR)
LOGSDIR = join(WORKDIR, "logs")
if not os.path.exists(LOGSDIR):
    os.mkdir(LOGSDIR)
TMPDIR = join(WORKDIR, "tmp")
if not os.path.exists(TMPDIR):
    os.mkdir(TMPDIR)

def _resolve_workdir_path(path_value):
    path_str = str(path_value or "").strip()
    if path_str == "":
        return ""
    if os.path.isabs(path_str):
        return path_str
    return join(WORKDIR, path_str)


def _slugify_label(value):
    slug = re.sub(r"[^A-Za-z0-9._-]+", "_", str(value).strip())
    slug = re.sub(r"_+", "_", slug).strip("._-")
    return slug or "group"


def _parse_bool_config(section_name, key, default):
    if key not in section_name:
        return default
    value = section_name.get(key)
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        lowered = value.strip().lower()
        if lowered in {"true", "1", "yes"}:
            return True
        if lowered in {"false", "0", "no"}:
            return False
    raise ValueError(f"deseq2.{key} must be a boolean. Found: {value!r}")


def _parse_int_config(section_name, key, default, minimum=None):
    try:
        parsed = int(section_name.get(key, default))
    except (TypeError, ValueError):
        raise ValueError(f"deseq2.{key} must be an integer. Found: {section_name.get(key)!r}")
    if minimum is not None and parsed < minimum:
        raise ValueError(f"deseq2.{key} must be >= {minimum}. Found: {parsed}")
    return parsed


def _parse_float_config(section_name, key, default, minimum=None, maximum=None, inclusive_max=True):
    try:
        parsed = float(section_name.get(key, default))
    except (TypeError, ValueError):
        raise ValueError(f"deseq2.{key} must be numeric. Found: {section_name.get(key)!r}")
    if minimum is not None and parsed < minimum:
        raise ValueError(f"deseq2.{key} must be >= {minimum}. Found: {parsed}")
    if maximum is not None:
        if inclusive_max and parsed > maximum:
            raise ValueError(f"deseq2.{key} must be <= {maximum}. Found: {parsed}")
        if not inclusive_max and parsed >= maximum:
            raise ValueError(f"deseq2.{key} must be < {maximum}. Found: {parsed}")
    return parsed


DESEQ2_CONFIG = config.get("deseq2", {}) or {}
if not isinstance(DESEQ2_CONFIG, dict):
    raise ValueError("config['deseq2'] must be a mapping if provided.")

DESEQ2_ENABLED = _parse_bool_config(DESEQ2_CONFIG, "enabled", False)
DESEQ2_CONTRASTS_FILE = _resolve_workdir_path(DESEQ2_CONFIG.get("contrasts_tsv", ""))
DESEQ2_OUTDIR = _resolve_workdir_path(
    DESEQ2_CONFIG.get("outdir", join("results", "deseq2"))
)
DESEQ2_ALPHA = _parse_float_config(DESEQ2_CONFIG, "alpha", 0.05, minimum=0.0, maximum=1.0)
if DESEQ2_ALPHA == 0:
    raise ValueError("deseq2.alpha must be > 0.")
DESEQ2_LFC_THRESHOLD = _parse_float_config(
    DESEQ2_CONFIG, "lfc_threshold", 1.0, minimum=0.0
)
DESEQ2_MIN_REPLICATES_PER_GROUP = _parse_int_config(
    DESEQ2_CONFIG, "min_replicates_per_group", 2, minimum=1
)
DESEQ2_MIN_TOTAL_COUNT = _parse_int_config(
    DESEQ2_CONFIG, "min_total_count", 1, minimum=0
)
DESEQ2_REPORT_TOP_N = _parse_int_config(
    DESEQ2_CONFIG, "report_top_n", 50, minimum=1
)
DESEQ2_REPORT_LABEL_TOP_N = _parse_int_config(
    DESEQ2_CONFIG, "report_label_top_n", 15, minimum=0
)
DESEQ2_REPORT_MAX_TABLE_ROWS = _parse_int_config(
    DESEQ2_CONFIG, "report_max_table_rows", 50000, minimum=1
)
DESEQ2_SIZE_FACTOR_TYPE = str(
    DESEQ2_CONFIG.get("size_factor_type", "poscounts")
).strip()
if DESEQ2_SIZE_FACTOR_TYPE not in {"ratio", "poscounts", "iterate"}:
    raise ValueError(
        "deseq2.size_factor_type must be one of ratio, poscounts, iterate."
    )
DESEQ2_FIT_TYPE = str(DESEQ2_CONFIG.get("fit_type", "parametric")).strip()
if DESEQ2_FIT_TYPE not in {"parametric", "local", "mean", "glmGamPoi"}:
    raise ValueError(
        "deseq2.fit_type must be one of parametric, local, mean, glmGamPoi."
    )
DESEQ2_SHRINK_TYPE = str(DESEQ2_CONFIG.get("shrink_type", "apeglm")).strip()
if DESEQ2_SHRINK_TYPE not in {"apeglm", "ashr", "normal"}:
    raise ValueError("deseq2.shrink_type must be one of apeglm, ashr, normal.")
DESEQ2_P_ADJUST_METHOD = str(
    DESEQ2_CONFIG.get("p_adjust_method", "BH")
).strip()
if DESEQ2_P_ADJUST_METHOD not in {
    "holm",
    "hochberg",
    "hommel",
    "bonferroni",
    "BH",
    "BY",
    "fdr",
    "none",
}:
    raise ValueError(
        "deseq2.p_adjust_method must be one of holm, hochberg, hommel, bonferroni, BH, BY, fdr, none."
    )
DESEQ2_COOKS_CUTOFF = _parse_bool_config(
    DESEQ2_CONFIG, "cooks_cutoff", True
)
DESEQ2_INDEPENDENT_FILTERING = _parse_bool_config(
    DESEQ2_CONFIG, "independent_filtering", True
)
DESEQ2_VST_BLIND = _parse_bool_config(DESEQ2_CONFIG, "vst_blind", True)
DESEQ2_HTML_SELF_CONTAINED = _parse_bool_config(
    DESEQ2_CONFIG, "html_self_contained", True
)
DESEQ2_DESIGN_FACTORS = DESEQ2_CONFIG.get("design_factors", ["group"])
if isinstance(DESEQ2_DESIGN_FACTORS, str):
    DESEQ2_DESIGN_FACTORS = [DESEQ2_DESIGN_FACTORS]
if not isinstance(DESEQ2_DESIGN_FACTORS, list) or not all(
    isinstance(item, str) for item in DESEQ2_DESIGN_FACTORS
):
    raise ValueError("deseq2.design_factors must be a list of strings.")
DESEQ2_DESIGN_FACTORS = [item.strip() for item in DESEQ2_DESIGN_FACTORS if item.strip()]
if DESEQ2_ENABLED and DESEQ2_DESIGN_FACTORS != ["group"]:
    raise ValueError(
        "Future DESeq2 covariates are not implemented yet. Set deseq2.design_factors to ['group']."
    )

DESEQ2_HOST_GENE_FLANK_SIZE = int(
    config.get("tn5_motif", {}).get("host_gene_flank_size", 250)
)
DESEQ2_VIRUS_BIN_SIZE = int(config.get("tn5_motif", {}).get("virus_bin_size", 100))
DESEQ2_TRNA_GENE_FLANK_SIZE = int(config.get("tn5_motif", {}).get("trna_gene_flank_size", 100))
DESEQ2_CONTRASTS = []
DESEQ2_CONTRASTS_BY_COMPARISON = {}

if DESEQ2_ENABLED:
    if DESEQ2_CONTRASTS_FILE == "":
        raise ValueError(
            "deseq2.enabled is true, but deseq2.contrasts_tsv is empty."
        )
    if not os.path.isfile(DESEQ2_CONTRASTS_FILE) or not os.access(
        DESEQ2_CONTRASTS_FILE, os.R_OK
    ):
        raise FileNotFoundError(
            f"deseq2.contrasts_tsv is not readable: {DESEQ2_CONTRASTS_FILE}"
        )
    deseq2_container = str(
        config.get("containers", {}).get("deseq2_report", "")
    ).strip()
    if deseq2_container == "":
        raise ValueError(
            "deseq2.enabled is true, but containers.deseq2_report is not set."
        )

    try:
        DESEQ2_CONTRASTS_DF = pd.read_csv(
            DESEQ2_CONTRASTS_FILE, sep="\t", dtype=str
        ).fillna("")
    except Exception as exc:
        raise ValueError(
            f"Failed to read deseq2.contrasts_tsv '{DESEQ2_CONTRASTS_FILE}': {exc}"
        )

    if list(DESEQ2_CONTRASTS_DF.columns) != ["group1", "group2"]:
        raise ValueError(
            "deseq2.contrasts_tsv must have exactly two tab-delimited headers: group1 and group2."
        )
    if DESEQ2_CONTRASTS_DF.empty:
        raise ValueError("deseq2.contrasts_tsv does not contain any contrasts.")

    known_groups = set(GROUPNAME2SAMPLENAME.keys())
    for row_idx, row in enumerate(
        DESEQ2_CONTRASTS_DF.itertuples(index=False), start=2
    ):
        group1 = str(row.group1).strip()
        group2 = str(row.group2).strip()
        if group1 == "" or group2 == "":
            raise ValueError(
                f"deseq2.contrasts_tsv line {row_idx} contains an empty group."
            )
        if group1 == group2:
            raise ValueError(
                f"deseq2.contrasts_tsv line {row_idx} compares the same group '{group1}' to itself."
            )
        if group1 not in known_groups or group2 not in known_groups:
            raise ValueError(
                f"deseq2.contrasts_tsv line {row_idx} references unknown groups: {group1!r}, {group2!r}."
            )
        if len(GROUPNAME2SAMPLENAME[group1]) < DESEQ2_MIN_REPLICATES_PER_GROUP:
            raise ValueError(
                f"Group '{group1}' has fewer than deseq2.min_replicates_per_group={DESEQ2_MIN_REPLICATES_PER_GROUP} samples."
            )
        if len(GROUPNAME2SAMPLENAME[group2]) < DESEQ2_MIN_REPLICATES_PER_GROUP:
            raise ValueError(
                f"Group '{group2}' has fewer than deseq2.min_replicates_per_group={DESEQ2_MIN_REPLICATES_PER_GROUP} samples."
            )
        comparison = f"{_slugify_label(group1)}_vs_{_slugify_label(group2)}"
        if comparison in DESEQ2_CONTRASTS_BY_COMPARISON:
            raise ValueError(
                f"deseq2.contrasts_tsv produces duplicate comparison slug '{comparison}'."
            )
        contrast = {
            "comparison": comparison,
            "group1": group1,
            "group2": group2,
        }
        DESEQ2_CONTRASTS.append(contrast)
        DESEQ2_CONTRASTS_BY_COMPARISON[comparison] = contrast

# Optional: print or return results
if "workflow" in globals() and getattr(workflow, "dryrun", False):
    print("SAMPLENAME2GROUPNAME:", SAMPLENAME2GROUPNAME)
    print("GROUPNAME2SAMPLENAME:", dict(GROUPNAME2SAMPLENAME))
    print("SAMPLENAMEISPE:", SAMPLENAMEISPE)
    print("SAMPLESDF:\n", SAMPLESDF)
    print("SAMPLES:\n", SAMPLES)
    print("Initialization complete.")
###################################################################################
