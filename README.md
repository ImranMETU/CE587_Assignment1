# CE587 Assignment 1 — Structural Optimization

This repository contains the final report and Excel Solver workbook for Questions 1–3, along with calculation code and supporting verification files.

## Main deliverables

- `CE587_HW1_Final.pdf` — final report
- `CE587_HW1_Solver.xlsx` — workbook with a Solver model for each question

## Project contents

- `CE587_HW1_Final.tex`, `report/` — report source, style, and Excel evidence
- `python/`, `Q2/`, `Q3/` — calculation and verification code
- `figures/submission/` — figures used in the final report
- `verification/` — independent checks and Solver run records
- `CE587 Homework 1 (2026-2027 Fall Semeter).pdf` — assignment source

## Rebuild the report

From the project root, with a LaTeX installation available:

```powershell
pdflatex -interaction=nonstopmode -halt-on-error CE587_HW1_Final.tex
pdflatex -interaction=nonstopmode -halt-on-error CE587_HW1_Final.tex
```

The Python calculation scripts use the packages listed in `python/requirements.txt`.
