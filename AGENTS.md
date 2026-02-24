# AGENTS.md

This file defines repository-level instructions for coding agents
working in `leo.gwas`.

## 1. Primary Rule

- Always read and follow `CONTRIBUTING.md` before making changes.
- Treat `CONTRIBUTING.md` as the development contract for this
  repository.

## 2. Language and Style

- Use Chinese for user-facing conversation by default.
- Keep code, comments, and technical identifiers in English.
- Keep changes minimal, focused, and reproducible.

## 3. Documentation and Exports

- Use roxygen2 for documentation updates.
- Do not hand-edit files under `man/`.
- If you add or export a new function, ensure it is included in
  `_pkgdown.yml` `reference.contents` (or mark as internal).

## 4. Minimum Validation Before Push

Run these checks from repo root:

1.  `Rscript -e "devtools::document('.')"`
2.  `Rscript -e "pkgdown::check_pkgdown('.')"`
3.  `Rscript -e \"pkgdown::build_site_github_pages(new_process = FALSE, install = FALSE)\"`

If a check is skipped, explain why in the final report.

## 5. CI and Safety

- Do not push changes that introduce deterministic CI failures.
- Prefer small, auditable commits with clear commit messages.
