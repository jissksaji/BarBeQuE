#!/usr/bin/env python3
"""Compare tea-list-linked amplicons with tea-list database sequences.

Amplicons are unique consensus accessions whose accession_taxid OR assigned_taxid
matches a taxid in teeliste.tsv. Database sequences are summed from the selected
database's taxid-count file for the same unique tea-list taxids.
"""

import argparse
import csv
import math
from pathlib import Path

from PIL import Image, ImageDraw

from plot_five_primer_rank_resolution import font, load_tea_taxids
from plot_selected_primers_horizontal_resolution import PRIMERS, load_barcodes


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", type=Path, required=True)
    parser.add_argument("--tea-list", type=Path, required=True)
    parser.add_argument("--primer-catalog", type=Path, required=True)
    parser.add_argument("--database", default="euphyllophyta")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--output-tsv", type=Path, required=True)
    return parser.parse_args()


def database_sequence_count(path: Path, tea_taxids: set[str]) -> int:
    """Sum reference sequences whose exact taxid occurs in the tea list."""
    total = 0
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if len(row) < 2 or row[0].strip() not in tea_taxids:
                continue
            try:
                total += int(float(row[1]))
            except ValueError:
                continue
    return total


def amplified_accession_count(path: Path, tea_taxids: set[str]) -> int:
    """Count each amplified accession once when either taxid field matches."""
    accessions = set()
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if len(row) < 6:
                continue
            accession_taxid = row[2].strip()
            assigned_taxid = row[5].strip()
            if accession_taxid in tea_taxids or assigned_taxid in tea_taxids:
                accessions.add(row[1].strip())
    return len(accessions)


def write_counts(path: Path, rows: list[dict[str, object]]) -> None:
    """Export exact plotted values for inspection and reuse."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            fieldnames=[
                "primer",
                "barcode",
                "amplified_amplicons",
                "database_sequences",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    tea_taxids = set(load_tea_taxids(args.tea_list))
    barcodes = load_barcodes(args.primer_catalog)
    db_path = (
        args.results_dir
        / "build_db_taxids"
        / f"{args.database}.db_taxids_counts.tsv"
    )
    db_sequences = database_sequence_count(db_path, tea_taxids)

    rows = []
    for primer in PRIMERS:
        consensus_path = (
            args.results_dir
            / "consensus"
            / f"{primer}_{args.database}.cluster_consensus.tsv"
        )
        if not consensus_path.is_file():
            continue
        rows.append(
            {
                "primer": primer,
                "barcode": barcodes.get(primer, "Unknown"),
                "amplified_amplicons": amplified_accession_count(
                    consensus_path, tea_taxids
                ),
                "database_sequences": db_sequences,
            }
        )

    write_counts(args.output_tsv, rows)

    # A log10 scale is necessary because database and amplicon counts differ by
    # several orders of magnitude. Exact values remain printed beside each bar.
    width = 3800
    row_height = 100
    top = 245
    bottom_margin = 230
    height = top + len(rows) * row_height + bottom_margin
    image = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(image)
    left, right = 1540, 3400
    maximum_log = math.ceil(math.log10(db_sequences + 1))

    def pixel_x(value: int) -> float:
        return left + (math.log10(value + 1) / maximum_log) * (right - left)

    ticks = [
        (0, "0"),
        (10, "10"),
        (100, "100"),
        (1_000, "1k"),
        (10_000, "10k"),
        (100_000, "100k"),
        (1_000_000, "1m"),
    ]
    for value, label in ticks:
        if math.log10(value + 1) > maximum_log:
            continue
        x = pixel_x(value)
        draw.line((x, top - 25, x, top + len(rows) * row_height), fill="#E1E5E9", width=2)
        draw.text((x, top - 40), label, font=font(28), fill="#59636E", anchor="ms")

    amplified_color = "#4C8CC8"
    database_color = "#9A928A"

    for index, row in enumerate(rows):
        center_y = top + index * row_height + row_height / 2
        primer = str(row["primer"])
        barcode = str(row["barcode"])
        amplicons = int(row["amplified_amplicons"])
        database = int(row["database_sequences"])

        draw.text(
            (left - 28, center_y),
            f"{primer}  [{barcode}]",
            font=font(32),
            fill="#30363D",
            anchor="rm",
        )

        # Two aligned bars per primer: amplicons above, database sequences below.
        draw.rectangle(
            (left, center_y - 31, pixel_x(amplicons), center_y - 3),
            fill=amplified_color,
        )
        draw.rectangle(
            (left, center_y + 5, pixel_x(database), center_y + 33),
            fill=database_color,
        )
        draw.text(
            (pixel_x(amplicons) + 18, center_y - 17),
            f"{amplicons:,}",
            font=font(24, bold=True),
            fill="#245B8A",
            anchor="lm",
        )
        draw.text(
            (pixel_x(database) + 18, center_y + 19),
            f"{database:,}",
            font=font(24, bold=True),
            fill="#5F5954",
            anchor="lm",
        )

    draw.text(
        (55, 35),
        "Tea-list amplicons vs reference-database sequences",
        font=font(58, bold=True),
        fill="#20242A",
    )
    draw.text(
        (55, 112),
        (
            "Unique amplicons match accession or assigned taxid · "
            "Database sequences use exact tea-list taxids · log10 scale"
        ),
        font=font(29),
        fill="#59636E",
    )

    legend_y = top + len(rows) * row_height + 45
    draw.rectangle((55, legend_y, 90, legend_y + 35), fill=amplified_color)
    draw.text((108, legend_y - 3), "Amplified amplicons", font=font(29), fill="#30363D")
    draw.rectangle((500, legend_y, 535, legend_y + 35), fill=database_color)
    draw.text((553, legend_y - 3), "Database sequences", font=font(29), fill="#30363D")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    image.save(args.output, format="PNG", optimize=True)


if __name__ == "__main__":
    main()
