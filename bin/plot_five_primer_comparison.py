#!/usr/bin/env python3
"""Create a static database-representation versus amplification comparison.

The script plots the five requested primers together. Each point represents one
target taxon from teeliste.tsv for one primer. Taxa are matched using the
original accession_taxid, which is the same rule used by pages/analysis.py.
"""

import argparse
import csv
import math
from collections import Counter, defaultdict
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


PRIMERS = [
    "ITS1_ITS2_White1990",
    "ITS2_collapsed",
    "Pr33",
    "Pr56",
    "Pr62",
]

# Fixed colors and markers keep primer identities consistent and distinguishable.
STYLES = {
    "ITS1_ITS2_White1990": ("#2563A6", "o"),
    "ITS2_collapsed": ("#D28E00", "s"),
    "Pr33": ("#C95D3A", "^"),
    "Pr56": ("#6B8E23", "D"),
    "Pr62": ("#9A5FA3", "P"),
}


def font(size: int, bold: bool = False) -> ImageFont.FreeTypeFont:
    """Load a readable system font for the exported image."""
    name = "DejaVuSans-Bold.ttf" if bold else "DejaVuSans.ttf"
    return ImageFont.truetype(f"/usr/share/fonts/truetype/dejavu/{name}", size)


def draw_marker(
    draw: ImageDraw.ImageDraw,
    x: float,
    y: float,
    color: str,
    marker: str,
    radius: int = 5,
) -> None:
    """Draw distinct marker shapes so primer identity does not rely on color."""
    box = (x - radius, y - radius, x + radius, y + radius)
    if marker == "o":
        draw.ellipse(box, fill=color, outline="white", width=1)
    elif marker == "s":
        draw.rectangle(box, fill=color, outline="white", width=1)
    elif marker == "^":
        draw.polygon(
            [(x, y - radius), (x - radius, y + radius), (x + radius, y + radius)],
            fill=color,
            outline="white",
        )
    elif marker == "D":
        draw.polygon(
            [(x, y - radius), (x - radius, y), (x, y + radius), (x + radius, y)],
            fill=color,
            outline="white",
        )
    else:
        draw.line((x - radius, y, x + radius, y), fill=color, width=3)
        draw.line((x, y - radius, x, y + radius), fill=color, width=3)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", type=Path, required=True)
    parser.add_argument("--tea-list", type=Path, required=True)
    parser.add_argument("--database", default="euphyllophyta")
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def load_targets(path: Path) -> list[dict[str, str]]:
    """Load the target names and taxids that define the plotted observations."""
    targets = []
    with path.open(newline="", encoding="utf-8") as handle:
        rows = csv.reader(handle, delimiter="\t")
        next(rows, None)
        for row in rows:
            if len(row) < 5:
                continue
            targets.append(
                {
                    "german_name": row[0].strip(),
                    "latin_name": row[1].strip(),
                    "taxid": row[2].strip(),
                    "rank": row[3].strip(),
                    "source": row[4].strip(),
                }
            )
    return targets


def load_database_counts(path: Path) -> dict[str, int]:
    """Return the number of reference sequences associated with each taxid."""
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


def load_amplified_counts(path: Path) -> Counter:
    """Count amplified accessions using their original reference taxid."""
    counts: Counter = Counter()
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if len(row) >= 3 and row[2].strip():
                counts[row[2].strip()] += 1
    return counts


def half_step_ticks(maximum: float) -> list[float]:
    """Return axis ticks at 0.5-unit intervals."""
    final_tick = math.floor(maximum * 2) / 2
    return [index / 2 for index in range(int(final_tick * 2) + 1)]


def main() -> None:
    args = parse_args()
    targets = load_targets(args.tea_list)
    database_path = (
        args.results_dir
        / "build_db_taxids"
        / f"{args.database}.db_taxids_counts.tsv"
    )
    database_counts = load_database_counts(database_path)

    # The database x-coordinate is shared across primers for the same target.
    database_values = [database_counts.get(row["taxid"], 0) for row in targets]
    x_values = [math.log10(value + 1) for value in database_values]

    # PIL is used so the plot can be exported without a GUI or plotting package.
    width, height = 2200, 1500
    image = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(image)
    left, right, top, bottom = 175, 2120, 220, 1160
    summary_lines = []
    primer_points = []

    for primer in PRIMERS:
        consensus_path = (
            args.results_dir
            / "consensus"
            / f"{primer}_{args.database}.cluster_consensus.tsv"
        )
        if not consensus_path.is_file():
            summary_lines.append(f"{primer}: run unavailable")
            continue

        amplified_counts = load_amplified_counts(consensus_path)
        amplified_values = [
            amplified_counts.get(row["taxid"], 0) for row in targets
        ]
        y_values = [math.log10(value + 1) for value in amplified_values]
        color, marker = STYLES[primer]

        primer_points.append((primer, color, marker, x_values, y_values))
        amplified_taxa = sum(value > 0 for value in amplified_values)
        summary_lines.append(
            f"{primer}: {amplified_taxa:,}/{len(targets):,} target taxa amplified"
        )

    all_x = [value for points in primer_points for value in points[3]]
    all_y = [value for points in primer_points for value in points[4]]
    x_max = max(1.0, max(all_x) * 1.04)
    y_max = max(1.0, max(all_y) * 1.08)

    def pixel_x(value: float) -> float:
        return left + (value / x_max) * (right - left)

    def pixel_y(value: float) -> float:
        return bottom - (value / y_max) * (bottom - top)

    # Quiet gridlines and exact tick labels preserve scale readability.
    for tick in half_step_ticks(x_max):
        x = pixel_x(float(tick))
        draw.line((x, top, x, bottom), fill="#E2E6EA", width=2)
        label = f"{tick:g}"
        draw.text((x, bottom + 16), label, font=font(24), fill="#59636E", anchor="ma")
    for tick in half_step_ticks(y_max):
        y = pixel_y(float(tick))
        draw.line((left, y, right, y), fill="#E2E6EA", width=2)
        draw.text((left - 20, y), f"{tick:g}", font=font(24), fill="#59636E", anchor="rm")

    draw.line((left, top, left, bottom), fill="#7B858F", width=3)
    draw.line((left, bottom, right, bottom), fill="#7B858F", width=3)

    # Plot each primer's taxa after the grid so every marker remains visible.
    for primer, color, marker, xs, ys in primer_points:
        for x_value, y_value in zip(xs, ys):
            draw_marker(draw, pixel_x(float(x_value)), pixel_y(float(y_value)), color, marker)

    draw.text(
        (left, 70),
        "Database representation vs amplification for five primers",
        font=font(42, bold=True),
        fill="#20242A",
    )
    draw.text(
        (left, 132),
        f"Database: {args.database} · Counts shown as log10(count + 1)",
        font=font(25),
        fill="#59636E",
    )
    draw.text(
        ((left + right) / 2, bottom + 78),
        "Database sequences (log10(count + 1))",
        font=font(28),
        fill="#30363D",
        anchor="ma",
    )

    # Draw the vertical y-axis label on a transparent layer and rotate it.
    label_layer = Image.new("RGBA", (900, 60), (255, 255, 255, 0))
    label_draw = ImageDraw.Draw(label_layer)
    label_draw.text(
        (450, 30),
        "Primer-amplified sequences (log10(count + 1))",
        font=font(28),
        fill="#30363D",
        anchor="mm",
    )
    rotated_label = label_layer.rotate(90, expand=True)
    image.paste(rotated_label, (35, int((top + bottom - rotated_label.height) / 2)), rotated_label)

    # Legend uses both color and shape, supporting grayscale interpretation.
    legend_x, legend_y = left + 25, top + 25
    draw.text((legend_x, legend_y), "Primer", font=font(26, bold=True), fill="#30363D")
    for index, (primer, color, marker, _, _) in enumerate(primer_points, start=1):
        y = legend_y + index * 43
        draw_marker(draw, legend_x + 8, y + 13, color, marker, radius=7)
        draw.text((legend_x + 30, y), primer, font=font(23), fill="#30363D")

    # A compact summary below the chart provides exact target-coverage counts.
    summary_y = 1300
    draw.text((left, summary_y - 45), "Coverage summary", font=font(25, bold=True), fill="#30363D")
    for index, line in enumerate(summary_lines):
        draw.text((left, summary_y + index * 34), line, font=font(22), fill="#4B5560")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    image.save(args.output, format="PNG", optimize=True)


if __name__ == "__main__":
    main()
