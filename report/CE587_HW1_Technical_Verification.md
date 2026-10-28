# CE587 Homework 1: technical verification record

This file accompanies the concise academic report. It retains the numerical and execution evidence that is useful for verification but does not belong in the submitted solution. The prior mathematical report remains available as `CE587_HW1_Report.pdf` and `report/CE587_HW1_Report.tex`; this redesign does not replace either file or modify the Solver workbook or validated model implementations.

## Q1

- Analytical minimum: 4.680000 kg, attained for `Rt = 1/pi cm^2` and `2 <= R <= sqrt(20/pi) cm`. The reported representative is the right endpoint, `R = 2.523132522020 cm`, `t = 0.126156626101 cm`.
- Python response: area `2.000000000000e-4 m^2`, inertia `6.366197723676e-8 m^4`, mass `4.680000000000 kg`, stress `250.000000 MPa`, deflection `1.206077 mm`, `R/t = 20.000000`.
- Excel Solver ran with return code 0. It returned `R = 2.523132521836 cm`, `t = 0.126156626092 cm`, mass `4.679999999317 kg`, and stress `250000000.036496 Pa`. The normalized stress residual `1.46e-10` is within the saved `1e-9` Solver precision and the `1e-7` feasibility tolerance. The tiny sub-bound mass is finite-precision slack.
- Seven Q1 tests passed. The 25-start SLSQP set included feasible and infeasible initial points; all returned candidates passed the independent feasibility and mass criteria. The full analytical optimal interval was also checked at 1001 points. The existing JSON files retain every run and residual.

## Q2

- Analytical result: `Ro = 0.400000 m`, `Ri = 0.398301201450 m`, wall thickness `1.698799 mm`, mass `334.447684923 kg`.
- Independent Python result: mass `334.447684930 kg`, bending `164.999999998 MPa`, shear `6.572009505 MPa`.
- Excel Solver ran with return code 0 and returned mass `334.447687562 kg`, bending `164.999998704 MPa`, and shear `6.572009453 MPa`. Its mass is about `2.63e-6 kg` above Python, consistent with Solver tolerance.
- The independent proof establishes shear redundancy and monotonic mass reduction along the bending boundary as `Ro` increases. The wall guard is `1e-8 m`; it is numerical only and far below the solved thickness. The preserved Q2 verification files include all starts, including the heavier feasible stationary result.

## Q3

- Analytical result: `H = 50.000000 cm`, `D = 3.423228131285 cm`, mass `6.605448065566 kg`.
- Python independently checked a three-dimensional tripod equilibrium in SI units; its calculated responses agree with the N-cm model. The height-bound proof shows mass increases with `H` along the buckling boundary throughout the assigned interval.
- Excel Solver ran with return code 0 and returned `H = 50 cm`, `D = 3.423228131300 cm`, mass `6.605448065623 kg`, a difference of about `5.70e-11 kg` from Python.
- Six Q3 tests passed. The saved JSON results retain all 25 starts and constraint residuals.

## Solver evidence

The Excel workbook `CE587_HW1_Solver.xlsx` contains the actual executed solutions and native Solver model blocks for all three questions. The submitted report reproduces their settings and final cell results in compact tables. It does not include a screenshot of the live Solver Parameters dialog; no authentic dialog capture was available for this edit. No synthetic screenshot was made.


## Final submission correction pass (2026-10-17)

The final report source and PDF were revised without changing the engineering model, verified results, workbook formulas, or Solver runs. The previous final PDF, source, and figure files are preserved in `checkpoints/final_submission_corrections_before/`.

- The submitted figures now pair the full assigned design domain with an optimum close-up for each question. Q1's full panel uses exactly 2 to 20 cm radius and 0.1 to 1 cm thickness; its zoom displays the complete nonunique optimum segment. Q2 uses the 0 to 40 cm radius domain and a close-up, with the 334.448 kg contour label separated from the optimum marker. Q3 uses the 50 to 200 cm height and 0.5 to 10 cm diameter domain, with the close-up and legend below the boundaries.
- Each formulation states its design vector, objective, and normalized inequalities. Q2 retains the strict mathematical condition `Ri < Ro` separately from the workbook's `1e-8 m` guard.
- Solver tables now separate setting, workbook cell reference, configuration, and result. Workbook formulas and setup sheets were read to verify the objective, changing cells, constraints, bounds, and GRG Nonlinear method.
- The actual workbook remains unchanged. The current desktop inventory had no Excel window available for capture, so neither worksheet nor Solver Parameters screenshots are included; the report identifies this evidence gap and records the saved native Solver settings.
- The revised PDF compiles to 8 A4 pages. The final LaTeX log has no warnings, unresolved references, or overfull boxes; all eight pages were rendered and visually inspected.

## Excel worksheet appendix (2026-10-17)

The final report now includes nine rendered views of the saved workbook: calculation, Solver setup, and graphical worksheets for Q1, Q2 and Q3. The workbook itself was not modified. The updated report compiles to 17 A4 pages, with the appendix occupying nine pages; the added pages were rendered and visually inspected. The sheet views are workbook renders rather than live Excel screenshots.

