#!/usr/bin/env python3
"""
Regenerates the lecture download-links table inside README.md.

Crawls the repo root for folders named `xx_<slug>` (xx = lecture number),
reads the lecture title out of `xx_<slug>_meta.tex` (\title{...}) inside
each folder, and writes a Markdown table with raw-GitHub download links for
the slide deck and article PDF. The table is written between the
<!-- LECTURES:START --> / <!-- LECTURES:END --> markers in README.md;
everything else in the file is left untouched.

Run from the repo root:
    python update_readme.py
"""

import re
from pathlib import Path

# --- Repo identity, used to build raw.githubusercontent.com download links ---
GITHUB_USERNAME = "ami-null"
GITHUB_REPO = "intro-regression"

BRANCH = "main"
REPO_ROOT = Path(__file__).resolve().parent
README_PATH = REPO_ROOT / "README.md"

FOLDER_PATTERN = re.compile(r"^(\d+)_(.+)$")
TITLE_PATTERN = re.compile(r"\\title\{(.*?)\}", re.DOTALL)

TABLE_START_MARKER = "<!-- LECTURES:START -->"
TABLE_END_MARKER = "<!-- LECTURES:END -->"


def raw_url(relative_path: str) -> str:
    """Build a raw.githubusercontent.com download link for a file in the repo."""
    return f"https://raw.githubusercontent.com/{GITHUB_USERNAME}/{GITHUB_REPO}/{BRANCH}/{relative_path}"


def find_meta_title(folder: Path, folder_name: str) -> str:
    """Read the lecture title out of the folder's *_meta.tex file."""
    meta_candidates = sorted(folder.glob("*_meta.tex"))
    if not meta_candidates:
        raise FileNotFoundError(f"No *_meta.tex file found in {folder_name}")
    if len(meta_candidates) > 1:
        print(f"Warning: multiple *_meta.tex files in {folder_name}, using {meta_candidates[0].name}")

    meta_text = meta_candidates[0].read_text(encoding="utf-8")
    match = TITLE_PATTERN.search(meta_text)
    if not match:
        raise ValueError(f"No \\title{{...}} found in {meta_candidates[0]}")

    # Collapse any internal newlines/extra whitespace from a multi-line \title{...}
    return " ".join(match.group(1).split())


def collect_lectures():
    """Find every xx_<slug> folder in the repo root and gather its data."""
    lectures = []

    for entry in REPO_ROOT.iterdir():
        if not entry.is_dir():
            continue

        match = FOLDER_PATTERN.match(entry.name)
        if not match:
            continue  # doesn't match xx_*** naming, skip

        number_str, _slug = match.groups()
        folder_name = entry.name

        try:
            title = find_meta_title(entry, folder_name)
        except (FileNotFoundError, ValueError) as e:
            print(f"Skipping {folder_name}: {e}")
            continue

        presentation_rel = f"{folder_name}/{folder_name}_presentation.pdf"
        article_rel = f"{folder_name}/{folder_name}_article.pdf"

        lectures.append({
            "number": int(number_str),
            "number_str": number_str,
            "title": title,
            "slides_url": raw_url(presentation_rel),
            "article_url": raw_url(article_rel),
        })

    lectures.sort(key=lambda lec: lec["number"])
    return lectures


def build_table(lectures) -> str:
    header = "| # | Title | Slides | Article |\n|---|---|---|---|"
    rows = [
        f"| {lec['number_str']} | {lec['title']} | [Download]({lec['slides_url']}) | [Download]({lec['article_url']}) |"
        for lec in lectures
    ]
    return "\n".join([header, *rows])


def update_readme(table_markdown: str):
    if not README_PATH.exists():
        raise FileNotFoundError(f"README.md not found at {README_PATH}")

    readme_text = README_PATH.read_text(encoding="utf-8")

    if TABLE_START_MARKER not in readme_text or TABLE_END_MARKER not in readme_text:
        raise ValueError(
            f"README.md is missing the {TABLE_START_MARKER} / {TABLE_END_MARKER} markers"
        )

    pattern = re.compile(
        re.escape(TABLE_START_MARKER) + r".*?" + re.escape(TABLE_END_MARKER),
        re.DOTALL,
    )
    replacement = f"{TABLE_START_MARKER}\n{table_markdown}\n{TABLE_END_MARKER}"
    new_readme_text = pattern.sub(replacement, readme_text)

    README_PATH.write_text(new_readme_text, encoding="utf-8")


def main():
    lectures = collect_lectures()
    if not lectures:
        print("No lecture folders found — README lecture table left empty.")
    table_markdown = build_table(lectures)
    update_readme(table_markdown)
    print(f"Updated {README_PATH} with {len(lectures)} lecture(s).")


if __name__ == "__main__":
    main()
