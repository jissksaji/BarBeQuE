#!/usr/bin/env python3
"""Split the combined database dot plot into two tightly cropped panels."""

import argparse
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


def font(size: int, bold: bool = False) -> ImageFont.FreeTypeFont:
    name = "DejaVuSans-Bold.ttf" if bold else "DejaVuSans.ttf"
    return ImageFont.truetype(f"/usr/share/fonts/truetype/dejavu/{name}", size)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output-first", type=Path, required=True)
    parser.add_argument("--output-second", type=Path, required=True)
    return parser.parse_args()


def make_panel(
    source: Image.Image,
    crop_box: tuple[int, int, int, int],
    title: str,
    subtitle: str,
) -> Image.Image:
    """Crop one panel and add a compact standalone header."""
    cropped = source.crop(crop_box)
    header_height = 115
    output = Image.new("RGB", (cropped.width, cropped.height + header_height), "white")
    output.paste(cropped, (0, header_height))
    draw = ImageDraw.Draw(output)
    # Remove the cropped fragment of the original shared panel heading; the new
    # standalone subtitle already identifies the rank range.
    draw.rectangle((0, header_height, 210, header_height + 95), fill="white")
    draw.text((25, 18), title, font=font(34, bold=True), fill="#20242A")
    draw.text((25, 65), subtitle, font=font(21), fill="#59636E")
    return output


def main() -> None:
    args = parse_args()
    source = Image.open(args.input).convert("RGB")

    # Crop away the shared title, legend, central gap, and outer margins.
    first = make_panel(
        source,
        (180, 175, 2260, 2815),
        "Reference-database sequence counts",
        "Tea-list ranks 1–42 · logarithmic sequence-count axis",
    )
    second = make_panel(
        source,
        (2460, 175, 4585, 2815),
        "Reference-database sequence counts",
        "Tea-list ranks 43–83 · logarithmic sequence-count axis",
    )

    args.output_first.parent.mkdir(parents=True, exist_ok=True)
    args.output_second.parent.mkdir(parents=True, exist_ok=True)
    first.save(args.output_first, format="PNG", optimize=True, dpi=(300, 300))
    second.save(args.output_second, format="PNG", optimize=True, dpi=(300, 300))


if __name__ == "__main__":
    main()
