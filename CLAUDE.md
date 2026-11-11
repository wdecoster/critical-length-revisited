# critical-length-revisited

Revisits read length vs structural variant discovery (and, as a follow-up, de novo assembly) with current ONT data
and tools, for blog posts. **This repository is public**: aggregate statistics only; no sample identifiers of
non-public data, no server names or paths, no cohort details.

## Documentation

Keep current state apart from history, so the docs do not drift:
- `docs/plan.md` and `docs/assembly_plan.md` describe the current plan and methods. When a method or decision
  changes, rewrite the relevant part to describe the current state and add a dated line to the plan's
  **Decisions** section (what, why) or its **Progress log** (runs and results, aggregate only). Do not leave
  superseded plans in the main text.
- `README.md`: what the workflow does now and how to run it.
- Blog figures come from `scripts/blog_figures.py`; results that end up in a post are reproducible from the workflow.
- Scripts take server paths from config (`config/server.yaml`, not committed), never hard-coded.
- Errors are not silenced: tolerate only failures our own code causes on purpose; anything else fails the job.
- Never sync to the server while a run of this repository is active.
