#!/usr/bin/env python3
"""Create a two-panel dot plot of database sequences for tea-list entries."""

import argparse
import csv
import math
from pathlib import Path

from PIL import Image, ImageDraw

from plot_five_primer_rank_resolution import font


DOT_COLORS = {
    "species": "#356FA3",
    "genus": "#D39A18",
    "varietas": "#287F7B",
    "variety": "#287F7B",
    "subspecies": "#287F7B",
}
OTHER_COLOR = "#7D756E"
ZERO_COLOR = "#B94A3D"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-tsv", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    with args.input_tsv.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    rows.sort(
        key=lambda row: (
            -int(row["database_sequences"]),
            row["latin_name"].lower(),
        )
    )

    split = math.ceil(len(rows) / 2)
    panels = [rows[:split], rows[split:]]
    width, height = 4600, 3050
    image = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(image)
    panel_starts = [70, 2350]
    panel_width = 2180
    label_width = 820
    plot_width = 1190
    top = 330
    row_height = 61

    maximum = max(int(row["database_sequences"]) for row in rows)
    maximum_log = math.ceil(math.log10(maximum + 1))

    def offset(value: int) -> float:
        if value <= 0:
            return 0
        return math.log10(value + 1) / maximum_log * plot_width

    ticks = [
        (0, "0"),
        (10, "10"),
        (100, "100"),
        (1_000, "1k"),
        (10_000, "10k"),
        (100_000, "100k"),
    ]

    for panel_index, panel in enumerate(panels):
        panel_x = panel_starts[panel_index]
        plot_x = panel_x + label_width
        first_rank = panel_index * split + 1
        last_rank = panel_index * split + len(panel)
        draw.text(
            (panel_x, 235),
            f"Ranks {first_rank}–{last_rank}",
            font=font(28, bold=True),
            fill="#30363D",
        )

        for tick, label in ticks:
            x = plot_x + offset(tick)
            draw.line(
                (x, top - 18, x, top + len(panel) * row_height),
                fill="#DCE1E6",
                width=2,
            )
            draw.text(
                (x, top - 30),
                label,
                font=font(21),
                fill="#59636E",
                anchor="ms",
            )

        for index, row in enumerate(panel):
            center_y = top + index * row_height + row_height / 2
            count = int(row["database_sequences"])
            rank = row["rank"].lower()
            color = DOT_COLORS.get(rank, OTHER_COLOR)

            draw.text(
                (plot_x - 20, center_y),
                f"{row['latin_name']}  [{rank}]",
                font=font(21),
                fill="#30363D",
                anchor="rm",
            )

            if count:
                x = plot_x + offset(count)
                # A quiet horizontal guide helps readers connect labels to dots
                # without turning the figure back into a bar chart.
                draw.line(
                    (plot_x, center_y, x, center_y),
                    fill="#D6DCE2",
                    width=2,
                )
                draw.ellipse(
                    (x - 11, center_y - 11, x + 11, center_y + 11),
                    fill=color,
                    outline="white",
                    width=2,
                )
                value_x = min(x + 17, panel_x + panel_width - 10)
            else:
                x = plot_x
                draw.line(
                    (x - 9, center_y - 9, x + 9, center_y + 9),
                    fill=ZERO_COLOR,
                    width=4,
                )
                draw.line(
                    (x - 9, center_y + 9, x + 9, center_y - 9),
                    fill=ZERO_COLOR,
                    width=4,
                )
                value_x = x + 17

            draw.text(
                (value_x, center_y),
                f"{count:,}",
                font=font(20, bold=True),
                fill="#30363D" if count else ZERO_COLOR,
                anchor="lm",
            )

    represented = sum(int(row["database_sequences"]) > 0 for row in rows)
    draw.text(
        (70, 42),
        "Reference-database sequence counts for tea-list entries",
        font=font(48, bold=True),
        fill="#20242A",
    )
    draw.text(
        (70, 110),
        f"{represented}/{len(rows)} entries represented · logarithmic sequence-count axis",
        font=font(26),
        fill="#59636E",
    )

    legend_y = 2900
    legend = [
        ("Species", DOT_COLORS["species"]),
        ("Genus", DOT_COLORS["genus"]),
        ("Variety / below species", DOT_COLORS["varietas"]),
        ("Other rank", OTHER_COLOR),
    ]
    draw.text((70, legend_y), "Taxonomic rank", font=font(27, bold=True), fill="#30363D")
    for index, (label, color) in enumerate(legend):
        x = 70 + index * 780
        draw.ellipse((x, legend_y + 47, x + 24, legend_y + 71), fill=color)
        draw.text((x + 38, legend_y + 43), label, font=font(23), fill="#30363D")
    zero_x = 3300
    draw.line(
        (zero_x, legend_y + 47, zero_x + 22, legend_y + 69),
        fill=ZERO_COLOR,
        width=4,
    )
    draw.line(
        (zero_x, legend_y + 69, zero_x + 22, legend_y + 47),
        fill=ZERO_COLOR,
        width=4,
    )
    draw.text(
        (zero_x + 38, legend_y + 43),
        "No exact-taxid sequences",
        font=font(23),
        fill="#30363D",
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    image.save(args.output, format="PNG", optimize=True, dpi=(300, 300))


if __name__ == "__main__":
    main()
