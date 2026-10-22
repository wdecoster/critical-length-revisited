# critical-length-revisited

This repo revisits [De Coster et al. 2020, *Critical length in long-read resequencing*](https://doi.org/10.1093/nargab/lqz027) (NAR Genom Bioinform 2(1):lqz027) with current data and tools, for a follow-up blog post. The first focus is structural variant discovery; phasing comes later.

- [docs/assessment.md](docs/assessment.md): a critical assessment of the original SV methods.
- [docs/plan.md](docs/plan.md): the replication plan. It uses the HG002 Q100 truth, ONT R10 reads (real reads cut to length plus simulation), Sniffles2 and cuteSV, and stratified truvari evaluation.

The blog post will appear on [Gigabase or gigabyte](https://wdecoster.github.io/gigabaseorgigabyte/).

The original code is at https://github.com/wdecoster/read_length_SV_discovery.
