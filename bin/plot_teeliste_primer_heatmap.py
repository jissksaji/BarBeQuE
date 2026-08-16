#!/usr/bin/env python3
"""Plot tea-list amplicon counts by primer beside database sequence counts."""

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path

from PIL import Image, ImageDraw

from plot_five_primer_rank_resolution import font
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


def load_targets(path: Path) -> list[dict[str, str]]:
    """Load every tea-list entry, retaining duplicate taxids as separate rows."""
    targets = []
    with path.open(newline="", encoding="utf-8") as handle:
        rows = csv.reader(handle, delimiter="\t")
        next(rows, None)
        for row in rows:
            if len(row) >= 4:
                targets.append(
                    {
                        "german_name": row[0].strip(),
                        "latin_name": row[1].strip(),
                        "taxid": row[2].strip(),
                        "rank": row[3].strip().lower() or "unknown",
                    }
                )
    return targets


def load_database_counts(path: Path) -> dict[str, int]:
    """Aggregate reference sequence counts by exact taxid."""
    counts: dict[str, int] = defaultdict(int)
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if len(row) < 2:
                continue
            try:
                counts[row[0].strip()] += int(float(row[1]))
            except ValueError:
                continue
    return dict(counts)


def load_primer_counts(path: Path, wanted_taxids: set[str]) -> dict[str, int]:
    """Count unique accessions per target using accession OR assigned taxid."""
    accessions_by_taxid: dict[str, set[str]] = defaultdict(set)
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if len(row) < 6:
                continue
            accession = row[1].strip()
            accession_taxid = row[2].strip()
            assigned_taxid = row[5].strip()
            if accession_taxid in wanted_taxids:
                accessions_by_taxid[accession_taxid].add(accession)
            if assigned_taxid in wanted_taxids:
                accessions_by_taxid[assigned_taxid].add(accession)
    return {taxid: len(accessions) for taxid, accessions in accessions_by_taxid.items()}


def blend(low: tuple[int, int, int], high: tuple[int, int, int], amount: float) -> str:
    """Interpolate between two RGB colors and return a Pillow-compatible hex."""
    amount = max(0.0, min(1.0, amount))
    rgb = tuple(round(a + (b - a) * amount) for a, b in zip(low, high))
    return "#{:02x}{:02x}{:02x}".format(*rgb)


def red_white_blue(amount: float) -> str:
    """Map low values to red, middle values to white, and high values to blue."""
    red = (178, 24, 43)
    white = (247, 247, 247)
    blue = (33, 102, 172)
    if amount <= 0.5:
        return blend(red, white, amount * 2)
    return blend(white, blue, (amount - 0.5) * 2)


def compact_count(value: int) -> str:
    if value >= 100_000:
        return f"{value / 1000:.0f}k"
    if value >= 10_000:
        return f"{value / 1000:.1f}k"
    return f"{value:,}"


def write_tsv(
    path: Path,
    targets: list[dict[str, object]],
    primer_counts: dict[str, dict[str, int]],
) -> None:
    """Export the full taxon-by-primer matrix with exact values."""
    fields = [
        "german_name",
        "latin_name",
        "taxid",
        "rank",
        "database_sequences",
        *PRIMERS,
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        for target in targets:
            row = dict(target)
            taxid = str(target["taxid"])
            for primer in PRIMERS:
                row[primer] = primer_counts.get(primer, {}).get(taxid, 0)
            writer.writerow(row)


def main() -> None:
    args = parse_args()
    targets = load_targets(args.tea_list)
    wanted_taxids = {target["taxid"] for target in targets}
    barcodes = load_barcodes(args.primer_catalog)
    database_path = (
        args.results_dir
        / "build_db_taxids"
        / f"{args.database}.db_taxids_counts.tsv"
    )
    database_counts = load_database_counts(database_path)

    primer_counts = {}
    for primer in PRIMERS:
        consensus_path = (
            args.results_dir
            / "consensus"
            / f"{primer}_{args.database}.cluster_consensus.tsv"
        )
        if consensus_path.is_file():
            primer_counts[primer] = load_primer_counts(
                consensus_path, wanted_taxids
            )

    for target in targets:
        target["database_sequences"] = database_counts.get(target["taxid"], 0)
    targets.sort(
        key=lambda target: (
            -int(target["database_sequences"]),
            str(target["latin_name"]).lower(),
        )
    )
    write_tsv(args.output_tsv, targets, primer_counts)

    cell_width = 145
    cell_height = 54
    left = 820
    top = 780
    columns = [*PRIMERS, "Database"]
    width = left + len(columns) * cell_width + 100
    height = top + len(targets) * cell_height + 260
    image = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(image)

    maximum = max(int(target["database_sequences"]) for target in targets)
    maximum_log = math.log10(maximum + 1)
    primer_low, primer_high = (240, 246, 251), (36, 91, 138)
    database_low, database_high = (252, 248, 235), (181, 125, 0)

    # Rotate long column labels to preserve space for all 21 primers.
    for column_index, column in enumerate(columns):
        if column == "Database":
            label = "Database sequences"
        else:
            label = f"{column} [{barcodes.get(column, 'Unknown')}]"
        label_layer = Image.new("RGBA", (650, cell_width), (255, 255, 255, 0))
        label_draw = ImageDraw.Draw(label_layer)
        label_draw.text(
            (10, cell_width / 2),
            label,
            font=font(22),
            fill="#30363D",
            anchor="lm",
        )
        rotated = label_layer.rotate(90, expand=True)
        x = left + column_index * cell_width + (cell_width - rotated.width) / 2
        image.paste(rotated, (round(x), top - rotated.height - 18), rotated)

    for row_index, target in enumerate(targets):
        y = top + row_index * cell_height
        taxid = str(target["taxid"])
        draw.text(
            (left - 18, y + cell_height / 2),
            f"{target['latin_name']} [{target['rank']}]",
            font=font(21),
            fill="#30363D",
            anchor="rm",
        )

        values = [
            primer_counts.get(primer, {}).get(taxid, 0) for primer in PRIMERS
        ]
        values.append(int(target["database_sequences"]))

        for column_index, value in enumerate(values):
            x = left + column_index * cell_width
            intensity = math.log10(value + 1) / maximum_log if value else 0
            if column_index == len(columns) - 1:
                color = blend(database_low, database_high, intensity)
            else:
                color = blend(primer_low, primer_high, intensity)
            draw.rectangle(
                (x, y, x + cell_width - 2, y + cell_height - 2),
                fill=color,
            )
            text_color = "white" if intensity >= 0.58 else "#30363D"
            draw.text(
                (x + cell_width / 2, y + cell_height / 2),
                compact_count(value),
                font=font(17, bold=True),
                fill=text_color,
                anchor="mm",
            )

    draw.text(
        (55, 38),
        "Tea-list amplicons recovered by primer",
        font=font(54, bold=True),
        fill="#20242A",
    )
    draw.text(
        (55, 112),
        (
            "Rows sorted by database representation · cells show unique amplicons "
            "matching accession or assigned taxid · log10 color scale"
        ),
        font=font(27),
        fill="#59636E",
    )

    # Compact gradient legend.
    legend_y = height - 145
    draw.text((55, legend_y), "Cell color (log10 count + 1)", font=font(26, bold=True), fill="#30363D")
    gradient_x = 420
    for index in range(240):
        amount = index / 239
        color = blend(primer_low, primer_high, amount)
        draw.rectangle(
            (gradient_x + index * 4, legend_y, gradient_x + index * 4 + 4, legend_y + 30),
            fill=color,
        )
    draw.text((gradient_x, legend_y + 42), "0", font=font(20), fill="#59636E", anchor="ma")
    draw.text(
        (gradient_x + 960, legend_y + 42),
        f"{maximum:,}",
        font=font(20),
        fill="#59636E",
        anchor="ma",
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    image.save(args.output, format="PNG", optimize=True)


if __name__ == "__main__":
    main()
