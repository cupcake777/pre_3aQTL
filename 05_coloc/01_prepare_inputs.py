#!/usr/bin/env python3
import csv
import gzip
import os
import re
import subprocess
import threading
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed


QTL_HEADER = [
    "chr",
    "pos",
    "ref",
    "alt",
    "phenotype_id",
    "variant_id",
    "start_distance",
    "end_distance",
    "af",
    "ma_samples",
    "ma_count",
    "pval_nominal",
    "slope",
    "slope_se",
]

GWAS_MA_HEADER = ["SNP", "A1", "A2", "freq", "b", "se", "p", "n"]
ROOT_DIR = Path(__file__).resolve().parent.parent
SIG_DIR = ROOT_DIR / "03_QTL"
STAGES = ["stage1", "stage2", "stage3", "stage4"]
_UNION_LOCK = threading.Lock()

LEGACY_EASYCOLOC_ROOTS = [
    Path("/mnt/share_group_folder/work/EasyColoc"),
    Path("/home/lyc/share_group_folder/work/EasyColoc"),
]


def chrom_sort_key(chrom: str) -> tuple[int, str]:
    value = chrom.lower().removeprefix("chr")
    mapping = {"x": 23, "y": 24, "m": 25, "mt": 25}
    if value in mapping:
        return (mapping[value], value)
    if value.isdigit():
        return (int(value), value)
    return (999, value)


def parse_variant_id(variant_id: str) -> tuple[str, str, str, str] | None:
    match = re.match(r"^(?:chr)?([^:]+):([0-9]+):([ACGT]+):([ACGT]+)$", variant_id)
    if not match:
        return None
    chrom, pos, ref, alt = match.groups()
    return (f"chr{chrom}", pos, ref, alt)


def load_gwas_meta(config_path: Path) -> dict[str, dict[str, str]]:
    meta: dict[str, dict[str, str]] = {}
    current_id: str | None = None
    with config_path.open() as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            stripped = line.strip()
            if stripped.startswith("- id:"):
                current_id = stripped.split(":", 1)[1].strip().strip('"')
                meta[current_id] = {}
                continue
            if not current_id or not stripped or stripped.startswith("#"):
                continue
            if ":" not in stripped:
                continue
            key, value = stripped.split(":", 1)
            key = key.strip()
            value = value.strip().strip('"')
            if key in {"build", "pop", "type", "sample_size_n", "prop", "name"}:
                meta[current_id][key] = value
    return meta


def load_stage_variant_map(
    stage_for_smr_path: Path,
) -> dict[str, tuple[str, str, str, str]]:
    variant_map: dict[str, tuple[str, str, str, str]] = {}
    with stage_for_smr_path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            snp = str(row.get("SNP", "")).strip()
            chrom = str(row.get("Chr", "")).strip()
            pos = str(row.get("BP", "")).strip()
            ref = str(row.get("A1", "")).strip()
            alt = str(row.get("A2", "")).strip()
            if not snp or not chrom or not pos or not ref or not alt:
                continue
            variant_map.setdefault(
                snp, (f"chr{chrom.removeprefix('chr')}", pos, ref, alt)
            )
    return variant_map


def count_samples(sample_file: Path) -> int:
    if not sample_file.exists():
        return 0
    with sample_file.open() as handle:
        return sum(1 for line in handle if line.strip())


def convert_qtl_record(
    record: dict[str, str], variant_map: dict[str, tuple[str, str, str, str]]
) -> list[str] | None:
    variant_id = str(record.get("variant_id", "")).strip()
    parsed = parse_variant_id(variant_id)
    if parsed is None:
        parsed = variant_map.get(variant_id)
    if parsed is None:
        return None
    chrom, pos, ref, alt = parsed

    # Handle stage4 (interaction) which has different column names
    pval = (
        record.get("pval_nominal")
        or record.get("interaction_pval")
        or record.get("pval_g", "")
    )
    slope = record.get("slope", "")
    slope_se = record.get("slope_se", "")

    return [
        chrom,
        pos,
        ref,
        alt,
        str(record.get("phenotype_id", "")),
        variant_id,
        str(record.get("start_distance", "")),
        str(record.get("end_distance", "")),
        str(record.get("af", "")),
        str(record.get("ma_samples", "")),
        str(record.get("ma_count", "")),
        pval,
        slope,
        slope_se,
    ]


def load_qtl_rows(
    qtl_path: Path,
    variant_map: dict[str, tuple[str, str, str, str]],
) -> tuple[list[list[str]], int, set[str]]:
    rows: list[list[str]] = []
    phenotypes: set[str] = set()
    skipped = 0
    with gzip.open(qtl_path, "rt") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for record in reader:
            output_row = convert_qtl_record(record, variant_map)
            if output_row is None:
                skipped += 1
                continue
            rows.append(output_row)
            phenotypes.add(output_row[4])
    rows.sort(key=lambda x: (chrom_sort_key(x[0]), int(x[1]), x[4], x[5]))
    return rows, skipped, phenotypes


def write_qtl_tables(
    stage: str,
    sig_path: Path,
    nominal_path: Path,
    stage_for_smr_path: Path,
    input_dir: Path,
    variant_union: set[str],
    rsid_union: set[str],
) -> dict[str, str] | None:
    variant_map = load_stage_variant_map(stage_for_smr_path)
    if not variant_map:
        raise RuntimeError(
            f"Could not build SNP coordinate map from {stage_for_smr_path}"
        )

    all_rows, all_skipped, all_phenotypes = load_qtl_rows(nominal_path, variant_map)
    sig_rows, sig_skipped, sig_phenotypes = load_qtl_rows(sig_path, variant_map)

    if not all_rows or not sig_rows:
        return None

    sig_plain = input_dir / f"{stage}.signif_pairs.txt"
    all_plain = input_dir / f"{stage}.all_pairs.txt"
    for output_path, rows in ((sig_plain, sig_rows), (all_plain, all_rows)):
        with output_path.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(QTL_HEADER)
            writer.writerows(rows)

    with _UNION_LOCK:
        for row in sig_rows:
            chrom, pos, ref, alt, _, variant_id = row[:6]
            variant_union.add(f"{chrom.removeprefix('chr')}:{pos}:{ref}:{alt}")
            variant_union.add(f"{chrom}:{pos}:{ref}:{alt}")
            variant_union.add(variant_id)
            rsid_union.add(variant_id)

    return {
        "stage": stage,
        "sig_plain": str(sig_plain),
        "all_plain": str(all_plain),
        "n_sig_rows": str(len(sig_rows)),
        "n_all_rows": str(len(all_rows)),
        "n_sig_phenotypes": str(len(sig_phenotypes)),
        "n_all_phenotypes": str(len(all_phenotypes)),
        "sig_skipped": str(sig_skipped),
        "all_skipped": str(all_skipped),
    }


def compress_and_index(tsv_path: Path) -> Path:
    gz_path = tsv_path.with_suffix(tsv_path.suffix + ".gz")
    subprocess.run(["bgzip", "-f", str(tsv_path)], check=True)
    subprocess.run(
        ["tabix", "-f", "-s", "1", "-b", "2", "-e", "2", "-S", "1", str(gz_path)],
        check=True,
    )
    return gz_path


def dataset_id_from_name(path: Path) -> str:
    name = path.name
    name = re.sub(r"_b(?:19|38)to38_harmonized\.tsv$", "", name)
    return name


def resolve_easycoloc_root(root_dir: Path) -> Path:
    env_root = os.environ.get("EASYCOLOC_ROOT")
    if env_root:
        candidate = Path(env_root).expanduser()
        if candidate.exists():
            return candidate
        raise FileNotFoundError(f"EASYCOLOC_ROOT does not exist: {candidate}")

    sibling_clone = root_dir.parent / "EasyColoc"
    if sibling_clone.exists():
        return sibling_clone

    for candidate in LEGACY_EASYCOLOC_ROOTS:
        if candidate.exists():
            return candidate

    searched = [
        str(root_dir.parent / "EasyColoc"),
        *[str(path) for path in LEGACY_EASYCOLOC_ROOTS],
    ]
    raise FileNotFoundError(
        "Could not locate EasyColoc. Set EASYCOLOC_ROOT to a cloned public EasyColoc repository. "
        f"Searched: {', '.join(searched)}"
    )


def write_subset_gwas(
    harmony_file: Path,
    dataset_id: str,
    gwas_subset_dir: Path,
    gwas_ma_dir: Path,
    variant_union: set[str],
    rsid_union: set[str],
    gwas_meta: dict[str, dict[str, str]],
) -> tuple[int, int]:
    subset_path = gwas_subset_dir / f"{dataset_id}.harmonized.tsv"
    ma_path = gwas_ma_dir / f"{dataset_id}.ma"
    meta = gwas_meta.get(dataset_id, {})
    default_n = meta.get("sample_size_n", "")
    kept = 0
    ma_kept = 0
    seen_ma: set[str] = set()

    with (
        harmony_file.open() as src,
        subset_path.open("w", newline="") as subset_out,
        ma_path.open("w", newline="") as ma_out,
    ):
        reader = csv.DictReader(src, delimiter="\t")
        subset_writer = csv.DictWriter(
            subset_out, fieldnames=reader.fieldnames, delimiter="\t"
        )
        ma_writer = csv.writer(ma_out, delimiter="\t")
        subset_writer.writeheader()
        ma_writer.writerow(GWAS_MA_HEADER)

        for row in reader:
            snpid = str(row.get("SNPID", "")).strip()
            rsid = str(row.get("rsID", "")).strip()
            if snpid not in variant_union and rsid not in rsid_union:
                continue

            kept += 1
            subset_writer.writerow(row)

            smr_snp = rsid if rsid and rsid != "." else snpid
            if not smr_snp or smr_snp in seen_ma:
                continue
            seen_ma.add(smr_snp)

            freq = row.get("EAF", "") or row.get("RAF", "")
            n_value = row.get("N", "") or default_n
            ma_writer.writerow(
                [
                    smr_snp,
                    row.get("EA", ""),
                    row.get("NEA", ""),
                    freq,
                    row.get("BETA", ""),
                    row.get("SE", ""),
                    row.get("P", ""),
                    n_value,
                ]
            )
            ma_kept += 1

    if kept == 0:
        subset_path.unlink(missing_ok=True)
        ma_path.unlink(missing_ok=True)
    return kept, ma_kept


def write_qtl_summary(
    stage_infos: list[dict[str, str]], qtl_summary_path: Path, root_dir: Path
) -> None:
    sample_map = {
        "stage1": count_samples(root_dir / "01_pre" / "sample.Stage1.final"),
        "stage2": count_samples(root_dir / "01_pre" / "sample.Stage2.final"),
        "stage3": count_samples(root_dir / "01_pre" / "sample.Stage3.final"),
        # stage4 is the pooled interaction dataset (Stage1 + Stage2 + Stage3 samples)
        "stage4": count_samples(
            root_dir / "03_QTL" / "interaction" / "pre_input" / "samples.interaction.txt"
        ),
    }
    with qtl_summary_path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "Type",
                "NumberRNASeqandGTSamples",
                "NumberRNASeqSamples",
                "Number_of_sigpheno",
                "allPairsTabixFilename",
                "sigPairsTabixFilename",
            ]
        )
        for info in stage_infos:
            sample_n = sample_map.get(info["stage"], 0)
            writer.writerow(
                [
                    info["stage"],
                    sample_n,
                    sample_n,
                    info["n_sig_phenotypes"],
                    info["all_gz"],
                    info["sig_gz"],
                ]
            )


def write_qtl_yaml(qtl_yaml_path: Path, qtl_summary_path: Path) -> None:
    qtl_yaml_path.write_text(
        "\n".join(
            [
                "qtl_info:",
                '  build: "hg38"',
                f'  file: "{qtl_summary_path}"',
                "  columns:",
                '    id: "Type"',
                '    sample_size: "NumberRNASeqandGTSamples"',
                '    all_filename: "allPairsTabixFilename"',
                '    sig_filename: "sigPairsTabixFilename"',
                "",
                "QTL_all_header:",
                '  - "chr"',
                '  - "pos"',
                '  - "ref"',
                '  - "alt"',
                '  - "phenotype_id"',
                '  - "variant_id"',
                '  - "start_distance"',
                '  - "end_distance"',
                '  - "af"',
                '  - "ma_samples"',
                '  - "ma_count"',
                '  - "pval_nominal"',
                '  - "slope"',
                '  - "slope_se"',
                "",
                "QTL_cols:",
                '  phenotype: "phenotype_id"',
                '  chrom: "chr"',
                '  pos: "pos"',
                '  pval: "pval_nominal"',
                '  beta: "slope"',
                '  se: "slope_se"',
                "",
            ]
        )
    )


def write_gwas_yaml(
    gwas_yaml_path: Path,
    subset_dir: Path,
    selected_ids: list[str],
    gwas_meta: dict[str, dict[str, str]],
) -> None:
    lines = ["datasets:"]
    for dataset_id in selected_ids:
        meta = gwas_meta.get(dataset_id, {})
        subset_file = subset_dir / f"{dataset_id}.harmonized.tsv"
        if not subset_file.exists():
            continue
        pop = meta.get("pop", dataset_id.split("_", 1)[0])
        gwas_type = meta.get("type", "cc")
        sample_n = meta.get("sample_size_n", "0")
        prop = meta.get("prop")
        lines.extend(
            [
                f'  - id: "{dataset_id}"',
                f'    name: "{dataset_id}"',
                f'    file: "{subset_file}"',
                '    build: "hg38"',
                f'    pop: "{pop}"',
                f'    type: "{gwas_type}"',
                f"    sample_size_n: {sample_n}",
            ]
        )
        if prop:
            lines.append(f"    prop: {prop}")
        lines.extend(
            [
                "    columns:",
                '      snp: "SNPID"',
                '      chrom: "CHR"',
                '      pos: "POS"',
                '      a1: "EA"',
                '      a2: "NEA"',
                '      pval: "P"',
                '      n: "N"',
                '      af: "EAF"',
                '      beta: "BETA"',
                '      se: "SE"',
            ]
        )
    gwas_yaml_path.write_text("\n".join(lines) + "\n")


def write_global_yaml(
    global_yaml_path: Path,
    easycoloc_root: Path,
    work_root: Path,
    results_dir: Path,
    gwas_subset_dir: Path,
) -> None:
    template = easycoloc_root / "config" / "local" / "global.yaml"
    text = template.read_text()
    replacements = {
        'project_root: "../.."': f'project_root: "{work_root}"',
        'output_dir: "../apa/04_coloc/results"': f'output_dir: "{results_dir}"',
        'temp_dir: "temp"': f'temp_dir: "{work_root / "05_coloc" / "temp"}"',
        'harmonize_dir: "harmony"': f'harmonize_dir: "{gwas_subset_dir}"',
    }
    for old, new in replacements.items():
        text = text.replace(old, new)
    global_yaml_path.write_text(text)


def write_run_script(
    run_script_path: Path, easycoloc_root: Path, config_dir: Path, results_dir: Path
) -> None:
    script = "\n".join(
        [
            "#!/usr/bin/env bash",
            "set -euo pipefail",
            f'repo_root="{easycoloc_root}"',
            f'results_dir="{results_dir}"',
            'run_log="${results_dir}/run_easycoloc.log"',
            'mkdir -p "${results_dir}"',
            'cd "${repo_root}"',
            f'export EASYCOLOC_GLOBAL_CONFIG="{config_dir / "global.yaml"}"',
            f'export EASYCOLOC_GWAS_CONFIG="{config_dir / "gwas.yaml"}"',
            f'export EASYCOLOC_QTL_CONFIG="{config_dir / "qtl.yaml"}"',
            'export EASYCOLOC_OUTPUT_DIR="${results_dir}"',
            'export EASYCOLOC_LOG_FILE="${run_log}"',
            'export EASYCOLOC_RUN_LABEL="pre3aqtl_05_coloc"',
            'exec Rscript run_coloc.r > "${run_log}" 2>&1',
            "",
        ]
    )
    run_script_path.write_text(script)
    run_script_path.chmod(0o755)


def main() -> int:
    root_dir = ROOT_DIR
    easycoloc_root = resolve_easycoloc_root(root_dir)
    coloc_root = root_dir / "05_coloc"
    input_dir = coloc_root / "input"
    gwas_subset_dir = coloc_root / "gwas_subset"
    gwas_ma_dir = coloc_root / "gwas_ma"
    config_dir = coloc_root / "config"
    results_dir = coloc_root / "results"

    for path in [input_dir, gwas_subset_dir, gwas_ma_dir, config_dir, results_dir]:
        path.mkdir(parents=True, exist_ok=True)

    sig_dir = SIG_DIR
    if not sig_dir.exists():
        raise FileNotFoundError(f"Cannot find required QTL input directory: {sig_dir}")

    stage_infos: list[dict[str, str]] = []
    variant_union: set[str] = set()
    rsid_union: set[str] = set()

    def process_stage(stage: str):
        if stage == "stage4":
            sig_path = sig_dir / "interaction" / "interaction_sig_QTL.txt.gz"
            # Use the full all-pairs nominal output generated by 04_interaction.sh
            # (interaction_all_pairs.txt.gz, produced without --best_only).
            # coloc.abf requires complete summary statistics for all variants in a
            # locus; the top-association-only file (cis_qtl_top_assoc.txt.gz) is
            # unsuitable and will produce unreliable posterior probabilities.
            nominal_path = sig_dir / "interaction" / "interaction_all_pairs.txt.gz"
            stage_for_smr_path = root_dir / "06_smr" / f"{stage}_for_smr.txt"
        else:
            sig_path = sig_dir / f"{stage}_sig_QTL.txt.gz"
            nominal_path = sig_dir / "nominal" / f"{stage}.cis_qtl_pairs.txt.gz"
            stage_for_smr_path = root_dir / "06_smr" / f"{stage}_for_smr.txt"

        if not sig_path.exists():
            print(f"[WARN] Missing {sig_path}, skip.")
            return None
        if not nominal_path.exists():
            print(f"[WARN] Missing {nominal_path}, skip stage.")
            return None
        if not stage_for_smr_path.exists():
            raise FileNotFoundError(f"Missing stage_for_smr file: {stage_for_smr_path}")

        info = write_qtl_tables(
            stage,
            sig_path,
            nominal_path,
            stage_for_smr_path,
            input_dir,
            variant_union,
            rsid_union,
        )
        if info is None:
            print(
                f"[WARN] {stage} had no usable chr:pos:ref:alt variants after parsing."
            )
            return None
        for key in ("sig_plain", "all_plain"):
            gz_path = compress_and_index(Path(info[key]))
            info["sig_gz" if key == "sig_plain" else "all_gz"] = str(gz_path)
        print(
            f"[INFO] {stage}: sig_rows={info['n_sig_rows']} all_rows={info['n_all_rows']} sig_phenotypes={info['n_sig_phenotypes']} all_phenotypes={info['n_all_phenotypes']}."
        )
        return info

    # Process stages in parallel
    with ThreadPoolExecutor(max_workers=4) as executor:
        futures = {executor.submit(process_stage, stage): stage for stage in STAGES}
        for future in as_completed(futures):
            stage = futures[future]
            try:
                info = future.result()
                if info:
                    stage_infos.append(info)
            except Exception as e:
                print(f"[ERROR] Stage {stage} failed: {e}")

    if not stage_infos:
        raise RuntimeError("No stage inputs were prepared.")

    qtl_summary_path = config_dir / "QTL_summary.csv"
    qtl_yaml_path = config_dir / "qtl.yaml"
    write_qtl_summary(stage_infos, qtl_summary_path, root_dir)
    write_qtl_yaml(qtl_yaml_path, qtl_summary_path)

    gwas_meta = load_gwas_meta(easycoloc_root / "config" / "local" / "gwas.yaml")
    harmony_dir = easycoloc_root / "harmony"
    harmony_files = sorted(harmony_dir.glob("*_harmonized.tsv"))

    selected_ids: list[str] = []
    for harmony_file in harmony_files:
        dataset_id = dataset_id_from_name(harmony_file)
        kept, ma_kept = write_subset_gwas(
            harmony_file,
            dataset_id,
            gwas_subset_dir,
            gwas_ma_dir,
            variant_union,
            rsid_union,
            gwas_meta,
        )
        if kept > 0:
            selected_ids.append(dataset_id)
            print(
                f"[INFO] {dataset_id}: kept {kept} harmonized rows, wrote {ma_kept} GWAS.ma rows."
            )

    if not selected_ids:
        raise RuntimeError("No GWAS rows overlapped the significant QTL variant set.")

    write_gwas_yaml(config_dir / "gwas.yaml", gwas_subset_dir, selected_ids, gwas_meta)
    write_global_yaml(
        config_dir / "global.yaml",
        easycoloc_root,
        root_dir,
        results_dir,
        gwas_subset_dir,
    )
    write_run_script(
        coloc_root / "02_run_easycoloc.sh", easycoloc_root, config_dir, results_dir
    )

    print(f"[DONE] Prepared 05_coloc inputs under {coloc_root}")
    print(f"[DONE] GWAS.ma files are in {gwas_ma_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
