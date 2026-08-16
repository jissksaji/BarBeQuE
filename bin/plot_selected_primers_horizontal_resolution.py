#!/usr/bin/env python3
"""Create a horizontal tea-list amplification and resolution comparison."""

import argparse
import csv
from pathlib import Path

from PIL import Image, ImageDraw

from plot_five_primer_rank_resolution import (
    CATEGORIES,
    font,
    load_target_resolutions,
    load_tea_taxids,
)


# Exact run names, in the order requested by the user.
PRIMERS = [
    "ITS1_ITS2_White1990",
    "ITS2_collapsed",
    "Pr33",
    "Pr56",
    "Pr62",
    "atpFatpH_Lahaye2008",
    "PCR_psbA-trnH",
    "psbKpsbI_Lahaye2008",
    "Pr26",
    "Pr28",
    "Pr31",
    "Pr6",
    "Rubisco-ribulose-1,5-bisphosphate-carboxylase_oxygenase-large-subunit-180",
    "rpoB_CBOL2009",
    "rpoC1_CBOL2009",
    "PCR_trnL",
    "Pr15",
    "Pr50",
    "intron-region-of-a-transfer-RNA-gene-133-tnrl",
    "intron-region-of-a-transfer-RNA-gene-200-1000",
    "trnLtrnF_ef_Taberlet1991",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", type=Path, required=True)
    parser.add_argument("--tea-list", type=Path, required=True)
    parser.add_argument("--primer-catalog", type=Path, required=True)
    parser.add_argument("--database", default="euphyllophyta")
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def load_barcodes(path: Path) -> dict[str, str]:
    """Map each exact primer name to its barcode/locus family."""
    mapping = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            mapping[row["primer"].strip()] = row["barcode"].strip()
    return mapping


def main() -> None:
    args = parse_args()
    tea_taxids = load_tea_taxids(args.tea_list)
    barcodes = load_barcodes(args.primer_catalog)
    data = []
    missing = []

    for primer in PRIMERS:
        consensus_path = (
            args.results_dir
            / "consensus"
            / f"{primer}_{args.database}.cluster_consensus.tsv"
        )
        if not consensus_path.is_file():
            missing.append(primer)
            continue
        counts = load_target_resolutions(consensus_path, tea_taxids)
        data.append((primer, barcodes.get(primer, "Unknown"), counts))

    if not data:
        raise SystemExit("None of the requested consensus runs was found.")

    total_targets = len(tea_taxids)
    width = 3800
    row_height = 86
    top = 215
    bottom_margin = 250
    height = top + len(data) * row_height + bottom_margin
    image = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(image)

    # Reserve substantial left space for long primer names plus barcode labels.
    left, right = 1550, 3480
    plot_width = right - left
    bar_height = 58

    def pixel_x(value: int) -> float:
        return left + (value / total_targets) * plot_width

    # Vertical count grid and x-axis labels.
    for tick in range(0, total_targets + 1, 10):
        x = pixel_x(tick)
        draw.line((x, top - 20, x, top + len(data) * row_height), fill="#E2E6EA", width=2)
        draw.text((x, top - 38), str(tick), font=font(30), fill="#59636E", anchor="ms")

    for index, (primer, barcode, counts) in enumerate(data):
        center_y = top + index * row_height + row_height / 2
        upper_y = center_y - bar_height / 2
        lower_y = center_y + bar_height / 2

        # Barcode is appended to the primer name as requested.
        draw.text(
            (left - 24, center_y),
            f"{primer}  [{barcode}]",
            font=font(34),
            fill="#30363D",
            anchor="rm",
        )

        cumulative = 0
        for category, color in CATEGORIES:
            value = counts.get(category, 0)
            if value == 0:
                continue
            x0 = pixel_x(cumulative)
            cumulative += value
            x1 = pixel_x(cumulative)
            draw.rectangle((x0, upper_y, x1, lower_y), fill=color, outline="white", width=2)

            # Only label segments wide enough to remain legible.
            if value >= 6:
                label_color = "#30363D" if category == "Not amplified" else "white"
                draw.text(
                    ((x0 + x1) / 2, center_y),
                    str(value),
                    font=font(27, bold=True),
                    fill=label_color,
                    anchor="mm",
                )

        amplified = total_targets - counts.get("Not amplified", 0)
        draw.text(
            (right + 22, center_y),
            f"{amplified}/{total_targets}",
            font=font(30, bold=True),
            fill="#30363D",
            anchor="lm",
        )

    draw.text(
        (55, 35),
        "Tea-list amplification by taxonomic resolution",
        font=font(60, bold=True),
        fill="#20242A",
    )
    draw.text(
        (right + 30, top - 38),
        "Amplified",
        font=font(28, bold=True),
        fill="#59636E",
        anchor="ms",
    )

    # Compact two-row legend beneath the bars.
    legend_y = top + len(data) * row_height + 40
    draw.text((55, legend_y), "Assigned resolution", font=font(34, bold=True), fill="#30363D")
    for index, (category, color) in enumerate(CATEGORIES):
        column = index % 4
        row = index // 4
        x = 55 + column * 900
        y = legend_y + 58 + row * 65
        draw.rectangle((x, y, x + 32, y + 32), fill=color)
        draw.text((x + 48, y - 5), category, font=font(29), fill="#30363D")

    if missing:
        draw.text(
            (width - 70, height - 35),
            "Missing: " + ", ".join(missing),
            font=font(18),
            fill="#8A3B2D",
            anchor="ra",
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    image.save(args.output, format="PNG", optimize=True)


if __name__ == "__main__":
    main()
