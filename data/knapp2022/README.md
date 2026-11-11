# Knapp et al. 2022: Dusy Basin 2012 swabs

- **File**: `dusy-2012-swabs.csv`, 371 skin swabs from subadult *Rana sierrae* in Dusy Basin,
  California, July to September 2012. Unmodified.
- **Source**: https://github.com/SNARL1/bd-mitigation-report, `data/dusy-2012-swabs.csv`,
  commit 949cdc0039 (2022-02-24); also deposited at Dryad,
  https://doi.org/10.25349/D9D90C. Retrieved 2026-11-10.
- **Paper**: Knapp RA, Joseph MB, Smith TC, Hegeman EE, et al. (2022) Effectiveness of
  antifungal treatments during chytridiomycosis epizootics in populations of an endangered
  frog. *PeerJ* 10:e12712. https://doi.org/10.7717/peerj.12712
- **Licence**: CC BY 4.0 (data; the repository's code is MIT). Attribution: Knapp et al. 2022.

## Fields

From the source repository's `data/README.md`:

| Field | Meaning |
|:------|:--------|
| `site_id` | 5-digit site code |
| `visit_date` | date of the visit (YYYY-MM-DD) |
| `group_date` | identical to `visit_date` here |
| `capture_stage` | life stage (subadult) |
| `toe_clip1` | group mark: L2 control, R2 treated |
| `swab_id` | swab code |
| `bd_load` | Bd ITS copies on the swab |
| `jliv_ge` | *Janthinobacterium lividum* genomic equivalents on the swab |
| `expt_or_wild` | in the experiment (`expt`) or not (`wild`) |
| `expt_trt` | `treated` or `control` (`NA` for wild frogs) |
| `trt_stage` | life stage treated |

## Design (checked against the paper)

Treated subadults received the antifungal itraconazole daily on 12–18 July 2012, were swabbed
on 19 July just before the first of two baths in *J. lividum* culture (19 and 20 July), and
were released into the pond. Frogs were group-marked, not individually tagged, so swabs
cannot be linked across visits. Control frogs were recaptured too rarely after 31 July to
compare. The *J. lividum* qPCR is species-level; the applied isolate came from a frog in the
same basin.
