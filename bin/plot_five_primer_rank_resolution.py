#!/usr/bin/env python3
"""Plot tea-list amplification split by assigned taxonomic resolution.

Each teeliste.tsv entry is counted once. A target matches when either its
original accession_taxid or the cluster's assigned_taxid equals the tea-list
taxid. The most specific assigned_rank reached by matching clusters determines
its category. Targets absent from every amplified cluster are Not amplified.
"""

import argparse
import csv
from collections import Counter
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


PRIMERS = [
    "ITS1_ITS2_White1990",
    "ITS2_collapsed",
    "Pr33",
    "Pr56",
    "Pr62",
]

# Detailed ranks are consolidated into reader-friendly resolution categories.
BELOW_SPECIES = {"varietas", "variety", "subspecies", "forma", "form", "genotype"}
CLASS_AND_ABOVE = {
    "class",
    "subclass",
    "phylum",
    "division",
    "domain",
    "kingdom",
}

CATEGORIES = [
    ("Variety / below species", "#245B8A"),
    ("Species", "#4C8CC8"),
    ("Genus-level", "#D5A021"),
    ("Family-level", "#D66A3A"),
    ("Order and above", "#8A5144"),
    ("Clade / other", "#8A8178"),
    ("Unresolved", "#C8CDD3"),
    ("Not amplified", "#EEF0F2"),
]

# Lower numbers represent more specific, preferred assignments.
CATEGORY_PRIORITY = {
    category: priority
    for priority, (category, _) in enumerate(CATEGORIES)
    if category != "Not amplified"
}


def font(size: int, bold: bool = False) -> ImageFont.FreeTypeFont:
    name = "DejaVuSans-Bold.ttf" if bold else "DejaVuSans.ttf"
    return ImageFont.truetype(f"/usr/share/fonts/truetype/dejavu/{name}", size)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", type=Path, required=True)
    parser.add_argument("--tea-list", type=Path, required=True)
    parser.add_argument("--database", default="euphyllophyta")
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def resolution_category(rank: str) -> str:
    """Map the exact assigned rank to one of the displayed resolution levels."""
    normalized = rank.strip().lower()
    if normalized in BELOW_SPECIES:
        return "Variety / below species"
    if normalized == "species":
        return "Species"
    if normalized in {"subsection", "section", "subgenus", "genus"}:
        return "Genus-level"
    if normalized in {"subtribe", "tribe", "subfamily", "family"}:
        return "Family-level"
    if normalized in {"order", "suborder"} or normalized in CLASS_AND_ABOVE:
        return "Order and above"
    if not normalized or normalized in {"no rank", "unknown", "unclassified", "na"}:
        return "Unresolved"
    return "Clade / other"


def load_tea_taxids(path: Path) -> list[str]:
    """Read one taxid for every entry in teeliste.tsv."""
    taxids = []
    with path.open(newline="", encoding="utf-8") as handle:
        rows = csv.reader(handle, delimiter="\t")
        next(rows, None)
        for row in rows:
            if len(row) >= 3 and row[2].strip():
                taxids.append(row[2].strip())
    return taxids


def load_target_resolutions(path: Path, tea_taxids: list[str]) -> Counter:
    """Categorize targets matched by either accession or assigned taxid."""
    wanted_taxids = set(tea_taxids)
    ranks_by_taxid: dict[str, set[str]] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if len(row) < 7:
                continue
            accession_taxid = row[2].strip()
            assigned_taxid = row[5].strip()
            assigned_rank = row[6].strip()

            # A set prevents a target from being counted twice when both taxids
            # are identical or when many accessions belong to the same cluster.
            if accession_taxid in wanted_taxids:
                ranks_by_taxid.setdefault(accession_taxid, set()).add(assigned_rank)
            if assigned_taxid in wanted_taxids:
                ranks_by_taxid.setdefault(assigned_taxid, set()).add(assigned_rank)

    counts: Counter = Counter()
    for taxid in tea_taxids:
        ranks = ranks_by_taxid.get(taxid)
        if not ranks:
            counts["Not amplified"] += 1
            continue
        categories = [resolution_category(rank) for rank in ranks]
        best_category = min(categories, key=CATEGORY_PRIORITY.get)
        counts[best_category] += 1
    return counts


def rounded_maximum(value: int) -> int:
    """Round the y-axis maximum upward to a convenient 5,000-cluster boundary."""
    interval = 5000
    return ((value + interval - 1) // interval) * interval


def main() -> None:
    args = parse_args()
    tea_taxids = load_tea_taxids(args.tea_list)
    data = {}
    missing = []

    for primer in PRIMERS:
        path = (
            args.results_dir
            / "consensus"
            / f"{primer}_{args.database}.cluster_consensus.tsv"
        )
        if path.is_file():
            data[primer] = load_target_resolutions(path, tea_taxids)
        else:
            missing.append(primer)

    if not data:
        raise SystemExit("None of the requested consensus files was found.")

    totals = {primer: sum(counts.values()) for primer, counts in data.items()}
    y_max = max(totals.values())

    width, height = 2200, 1650
    image = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(image)
    left, right, top, bottom = 210, 2110, 250, 1160

    def pixel_y(value: int) -> float:
        return bottom - (value / y_max) * (bottom - top)

    # Draw a restrained count grid behind the bars.
    tick_interval = 10
    for tick in range(0, y_max + 1, tick_interval):
        y = pixel_y(tick)
        draw.line((left, y, right, y), fill="#E1E5E9", width=2)
        draw.text(
            (left - 25, y),
            f"{tick:,}",
            font=font(23),
            fill="#59636E",
            anchor="rm",
        )

    draw.line((left, top, left, bottom), fill="#7B858F", width=3)
    draw.line((left, bottom, right, bottom), fill="#7B858F", width=3)

    primers = list(data)
    slot_width = (right - left) / len(primers)
    bar_width = 210

    for primer_index, primer in enumerate(primers):
        center_x = left + slot_width * (primer_index + 0.5)
        x0, x1 = center_x - bar_width / 2, center_x + bar_width / 2
        cumulative = 0
        total = totals[primer]
        amplified = total - data[primer].get("Not amplified", 0)

        for category, color in CATEGORIES:
            value = data[primer].get(category, 0)
            if value == 0:
                continue
            lower_y = pixel_y(cumulative)
            cumulative += value
            upper_y = pixel_y(cumulative)
            draw.rectangle(
                (x0, upper_y, x1, lower_y),
                fill=color,
                outline="white",
                width=2,
            )

            # Label only sufficiently large segments to prevent collisions.
            percentage = 100 * value / total
            if percentage >= 4:
                label_color = "white" if category != "Unresolved" else "#30363D"
                draw.text(
                    (center_x, (upper_y + lower_y) / 2),
                    f"{percentage:.0f}%",
                    font=font(22, bold=True),
                    fill=label_color,
                    anchor="mm",
                )

        draw.text(
            (center_x, pixel_y(total) - 18),
            f"{amplified}/{total} amplified",
            font=font(24, bold=True),
            fill="#30363D",
            anchor="mb",
        )
        draw.text(
            (center_x, bottom + 30),
            primer,
            font=font(23),
            fill="#30363D",
            anchor="ma",
        )

    draw.text(
        (left, 60),
        "Tea-list amplification by taxonomic resolution",
        font=font(42, bold=True),
        fill="#20242A",
    )
    # Vertical y-axis label.
    label_layer = Image.new("RGBA", (500, 60), (255, 255, 255, 0))
    label_draw = ImageDraw.Draw(label_layer)
    label_draw.text(
        (250, 30),
        "Tea-list targets",
        font=font(28),
        fill="#30363D",
        anchor="mm",
    )
    rotated = label_layer.rotate(90, expand=True)
    image.paste(rotated, (45, int((top + bottom - rotated.height) / 2)), rotated)

    # Legend follows the same fine-to-coarse order as the stacked bars.
    legend_x, legend_y = left, 1260
    draw.text((legend_x, legend_y), "Assigned resolution", font=font(25, bold=True), fill="#30363D")
    for index, (category, color) in enumerate(CATEGORIES):
        column = index % 4
        row = index // 4
        cursor_x = legend_x + column * 440
        cursor_y = legend_y + 52 + row * 58
        draw.rectangle((cursor_x, cursor_y, cursor_x + 24, cursor_y + 24), fill=color)
        draw.text(
            (cursor_x + 36, cursor_y - 4),
            category,
            font=font(20),
            fill="#30363D",
        )

    if missing:
        draw.text(
            (right, 1425),
            "Missing runs: " + ", ".join(missing),
            font=font(19),
            fill="#8A3B2D",
            anchor="ra",
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    image.save(args.output, format="PNG", optimize=True)


if __name__ == "__main__":
    main()
