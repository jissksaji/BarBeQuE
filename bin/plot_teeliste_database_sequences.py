#!/usr/bin/env python3
"""Plot reference-database sequence counts for all teeliste.tsv entries."""

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path

from PIL import Image, ImageDraw

from plot_five_primer_rank_resolution import font


RANK_COLORS = {
    "species": "#4C8CC8",
    "genus": "#D5A021",
    "varietas": "#286A91",
    "variety": "#286A91",
    "subspecies": "#286A91",
}
OTHER_COLOR = "#8A8178"
ZERO_COLOR = "#C8CDD3"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tea-list", type=Path, required=True)
    parser.add_argument("--database-counts", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--output-tsv", type=Path, required=True)
    return parser.parse_args()


def load_database_counts(path: Path) -> dict[str, int]:
    """Aggregate the database sequence count for each exact taxid."""
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


def load_targets(path: Path, database_counts: dict[str, int]) -> list[dict[str, object]]:
    """Load all 83 entries and attach their exact-taxid database counts."""
    targets = []
    with path.open(newline="", encoding="utf-8") as handle:
        rows = csv.reader(handle, delimiter="\t")
        next(rows, None)
        for row in rows:
            if len(row) < 5:
                continue
            taxid = row[2].strip()
            targets.append(
                {
                    "german_name": row[0].strip(),
                    "latin_name": row[1].strip(),
                    "taxid": taxid,
                    "rank": row[3].strip().lower() or "unknown",
                    "database_sequences": database_counts.get(taxid, 0),
                }
            )
    return sorted(
        targets,
        key=lambda row: (
            -int(row["database_sequences"]),
            str(row["latin_name"]).lower(),
        ),
    )


def write_tsv(path: Path, targets: list[dict[str, object]]) -> None:
    """Export the sorted values underlying the figure."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            fieldnames=[
                "german_name",
                "latin_name",
                "taxid",
                "rank",
                "database_sequences",
            ],
        )
        writer.writeheader()
        writer.writerows(targets)


def main() -> None:
    args = parse_args()
    database_counts = load_database_counts(args.database_counts)
    targets = load_targets(args.tea_list, database_counts)
    write_tsv(args.output_tsv, targets)

    # Split the ranked list into two columns to keep 83 labels readable.
    split = math.ceil(len(targets) / 2)
    panels = [targets[:split], targets[split:]]
    width, height = 4400, 2940
    image = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(image)
    top = 285
    row_height = 58
    panel_width = 2110
    panel_starts = [55, 2245]
    label_width = 780
    bar_width = 1130

    maximum = max(int(row["database_sequences"]) for row in targets)
    maximum_log = math.ceil(math.log10(maximum + 1))

    def bar_length(value: int) -> float:
        if value <= 0:
            return 0
        return (math.log10(value + 1) / maximum_log) * bar_width

    tick_values = [
        (0, "0"),
        (10, "10"),
        (100, "100"),
        (1_000, "1k"),
        (10_000, "10k"),
        (100_000, "100k"),
        (1_000_000, "1m"),
    ]

    for panel_index, panel in enumerate(panels):
        panel_x = panel_starts[panel_index]
        plot_x = panel_x + label_width

        draw.text(
            (panel_x, 205),
            f"Ranks {panel_index * split + 1}–{panel_index * split + len(panel)}",
            font=font(27, bold=True),
            fill="#59636E",
        )

        for tick, label in tick_values:
            if math.log10(tick + 1) > maximum_log:
                continue
            x = plot_x + bar_length(tick)
            draw.line(
                (x, top - 18, x, top + len(panel) * row_height),
                fill="#E1E5E9",
                width=2,
            )
            draw.text(
                (x, top - 28),
                label,
                font=font(21),
                fill="#59636E",
                anchor="ms",
            )

        for index, target in enumerate(panel):
            center_y = top + index * row_height + row_height / 2
            latin_name = str(target["latin_name"])
            rank = str(target["rank"])
            count = int(target["database_sequences"])
            color = RANK_COLORS.get(rank, OTHER_COLOR)

            draw.text(
                (plot_x - 20, center_y),
                f"{latin_name}  [{rank}]",
                font=font(21),
                fill="#30363D",
                anchor="rm",
            )

            if count > 0:
                end_x = plot_x + bar_length(count)
                draw.rectangle(
                    (plot_x, center_y - 16, end_x, center_y + 16),
                    fill=color,
                )
                label_x = min(end_x + 12, panel_x + panel_width - 8)
            else:
                draw.ellipse(
                    (plot_x - 6, center_y - 6, plot_x + 6, center_y + 6),
                    fill=ZERO_COLOR,
                )
                label_x = plot_x + 14

            draw.text(
                (label_x, center_y),
                f"{count:,}",
                font=font(20, bold=True),
                fill="#30363D",
                anchor="lm",
            )

    total_sequences = sum(int(row["database_sequences"]) for row in targets)
    represented = sum(int(row["database_sequences"]) > 0 for row in targets)
    draw.text(
        (55, 38),
        "Reference-database sequences for 83 tea-list entries",
        font=font(52, bold=True),
        fill="#20242A",
    )
    draw.text(
        (55, 108),
        (
            f"{represented}/{len(targets)} entries represented · "
            f"{total_sequences:,} sequences across entries · exact taxid matching · "
            "log10 scale"
        ),
        font=font(27),
        fill="#59636E",
    )

    legend_y = 2785
    legend_items = [
        ("Species", RANK_COLORS["species"]),
        ("Genus", RANK_COLORS["genus"]),
        ("Variety / below species", RANK_COLORS["varietas"]),
        ("Other rank", OTHER_COLOR),
        ("No database sequences", ZERO_COLOR),
    ]
    draw.text((55, legend_y), "Taxonomic rank", font=font(28, bold=True), fill="#30363D")
    for index, (label, color) in enumerate(legend_items):
        x = 55 + index * 780
        y = legend_y + 48
        draw.rectangle((x, y, x + 30, y + 30), fill=color)
        draw.text((x + 43, y - 3), label, font=font(24), fill="#30363D")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    image.save(args.output, format="PNG", optimize=True)


if __name__ == "__main__":
    main()
