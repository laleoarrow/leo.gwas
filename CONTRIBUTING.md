# leo.gwas Development Specification

This document defines practical development rules for `leo.gwas`.
The goal is to keep code, docs, and CI behavior consistent and reproducible.

## 1. Scope

- Applies to all code and docs in this repository.
- Preferred language:
  - Code/comments: English.
  - Commit messages/issues/PR discussion: English or Chinese.

## 2. Repository Structure

- `R/`: implementation files.
- `man/`: generated `.Rd` files from roxygen2 (do not hand-edit).
- `docs/`: generated pkgdown site output.
- `vignettes/`: long-form tutorials.
- `.github/workflows/`: CI workflows (`R-CMD-check`, `pkgdown`, `test-coverage`).
- `_pkgdown.yml`: pkgdown navigation and reference index.

## 3. Function and File Conventions

- Function naming:
  - Exported functions: `snake_case` or established project naming.
  - Internal helpers: prefix with `.`.
- Keep one coherent topic per file.
- Prefer explicit namespace calls (for example `dplyr::mutate`) inside functions.
- Use `stop(..., call. = FALSE)` for hard failures.
- For progress/status messages, prefer `leo.basic::leo_log()` in long-running pipelines.

## 4. roxygen2 and Documentation Rules

- Every exported function must include:
  - `@param` for all arguments.
  - `@return`.
  - `@examples` (use `\\dontrun{}` for heavy/network/system commands).
  - `@export`.
- Generate docs only via roxygen2:
  - `Rscript -e "devtools::document('.')"`
- Never manually edit `man/*.Rd` content.

## 5. pkgdown Rules (Critical)

- Every exported topic must appear in `_pkgdown.yml` `reference.contents`,
  or be marked with `@keywords internal`.
- If a new export is added but not indexed, CI will fail with:
  - `topic missing from index`.
- Before push, run:
  - `Rscript -e "pkgdown::check_pkgdown('.')"`
  - `Rscript -e \"pkgdown::build_site_github_pages(new_process = FALSE, install = FALSE)\"`

## 6. Dependency and Namespace Rules

- Add required runtime packages to `DESCRIPTION` `Imports`.
- Add optional packages to `Suggests`.
- Update `NAMESPACE` through roxygen2 (not manual edits when avoidable).
- For GitHub-only dependencies, document rationale in `DESCRIPTION` `Remotes`.

## 7. Local Validation Before Push

Run the minimal validation set:

1. `Rscript -e "devtools::document('.')"`
2. `Rscript -e "pkgdown::check_pkgdown('.')"`
3. `Rscript -e \"pkgdown::build_site_github_pages(new_process = FALSE, install = FALSE)\"`
4. `Rscript -e "rcmdcheck::rcmdcheck(args = c('--no-manual'), error_on = 'warning')"`

If step 4 is too heavy locally, at least ensure steps 1-3 pass before push.

## 8. CI Contract

- `R-CMD-check.yaml` must stay green on `main/master`.
- `pkgdown.yaml` must build and deploy without reference index errors.
- Do not merge changes that introduce deterministic CI failures.

## 9. Release and Changelog

- User-visible changes should update `NEWS.md`.
- Keep versioning in `DESCRIPTION` consistent with release intent.

## 10. Pull Request Checklist

- [ ] Code follows this spec.
- [ ] roxygen docs regenerated.
- [ ] `_pkgdown.yml` updated for new exported functions.
- [ ] Local validation commands passed (or failure reason documented).
- [ ] `NEWS.md` updated when behavior changes.
