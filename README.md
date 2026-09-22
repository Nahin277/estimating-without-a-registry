# Estimating a National Disease Burden With No National Survey

**How do you estimate how many women under 30 are diagnosed with breast cancer in Bangladesh, when the country has no population-based cancer registry, no national survey, and only a handful of single-hospital case series?**

This repository is the full, reproducible analysis behind the Medium post
**"She Was 28. Her Country Has No Number for That."**

It is a worked example of a problem that shows up constantly in global health and almost never gets handled honestly: **you have real data, and none of it is representative.**

The short version of the answer:

> You can estimate it. But only by turning the selection bias into a named, explicit parameter — and the honest output is a **range**, not a point estimate.

**Base R only. No packages. Runs in about four minutes.**

---

## The problem in one number

The WHO's Global Cancer Observatory publishes a fact sheet for Bangladesh. GLOBOCAN 2024: **9,480 new female breast cancer cases.** Specific, authoritative, no error bar.

Then you read the footnote:

> *"Rates based on incidence estimates or registry data from neighbouring countries."*

And the one attached to the prevalence figure:

> *"Computed using sex-, site- and age-specific incidence to 1-, 3- and 5-year prevalence ratios from **Nordic countries** for the period (2011-2020), and scaled using **Human Development Index (HDI) ratios**."*

This is not a scandal — IARC states it in public, every time, and borrowing is the responsible thing to do when you must publish a global number today. But it means that when you ask how many Bangladeshi women under 30 are diagnosed each year, the honest first answer is: **nobody has ever looked.**

What exists instead is a scatter of single-hospital studies. A few hundred women, every one of them a real observation, none of them a sample of the country.

---

## Quickstart

```bash
git clone https://github.com/<you>/<repo>.git
cd <repo>/R
Rscript run_all.R
```

No packages to install. Outputs land in `figures/` and `results/`.
To force a refit after editing inputs, delete `results/fit.rds` first.

Optional — rebuild the tables document (the only part that needs Python):

```bash
pip install -r requirements.txt
python3 make_tables_docx.py
```

---

## What the analysis shows

### 1. The shortcut gives a confident, wrong answer

Pool the six hospital case series with a proper random-effects meta-analysis and you get **5.2% (95% CI 3.9–6.9%)**, implying **492 women a year**. Publishable. Cited method. An interval.

But every single study sits on one side of the band where regional registries live. When independent measurements all miss in the same direction, the problem is not noise.

![Six hospital case series](figures/fig1_what_the_studies_say.png)

A hospital is not a country. A woman in a published Bangladeshi case series has already passed through a long chain of filters — noticing, being believed, reaching a facility that can biopsy, being referred onward, being treated somewhere that publishes — and **every one of those filters has an age gradient.**

### 2. More data does not fix it — it makes it worse

Simulate it. True national share 2%, every hospital series tilted the same modest amount. Pool 2 studies, 10, 500.

![Precision is not accuracy](figures/fig2_precision_is_not_accuracy.png)

| Studies pooled | Chance the 95% interval contains the truth |
|---:|---:|
| 2 | 32% |
| 5 | 7% |
| 10 | 1% |
| 20+ | **0%** |

The interval got *narrower* and *more confidently wrong* at the same time, and nothing in the output would warn you. Heterogeneity statistics won't catch it — the studies agree beautifully, because they share a bias.

> **Sampling error shrinks with more data. Selection bias does not. It just gets more precise.**

This is also why no amount of gradient boosting or deep learning rescues you. The women who never reach a hospital are not in the data under any column. There is nothing there for an algorithm to learn from. It is an identification problem, not a computational one.

### 3. The only honest move: name the bias

```
logit(what hospitals see) = logit(the national truth) + δ
```

**δ** is the referral tilt — how much more likely a young patient is to end up in a readable case series. `exp(δ)` is just "young patients are this many times over-represented."

One equation, two unknowns. The data pin down the left side and say nothing about how to split the right. Any claim about the national truth is a claim about δ in disguise.

So instead of picking a δ and hiding it, trace the whole curve:

![The identification curve](figures/fig3_the_honest_answer.png)

Every point on that line fits the hospital data exactly as well as every other point.

δ can be learned without collecting anything new: some places have **both** a population registry **and** hospital series from the same catchment, so you can measure the gap directly and transport it with its uncertainty attached. That assumption might be wrong — but it is stated, estimable, and attackable with evidence, which is the entire game.

---

## Results

| Model | Share under 30 | 95% interval | National count |
|---|---:|---:|---:|
| Naive pooling of the studies | 5.2% | 3.9–6.9% | 492 |
| + hierarchical structure | 4.7% | 3.0–6.6% | 448 |
| + regional registry prior | 3.7% | 2.7–5.2% | 352 |
| **+ referral-tilt correction (full)** | **2.0%** | **1.2–3.3%** | **192** |

![Posterior distributions](figures/fig4_posterior.png)

**The headline, stated the way it deserves to be stated:**

> Under any defensible assumption about referral patterns — young patients somewhere between **1.5× and 3×** over-represented in published case series — between **1.8% and 3.5%** of Bangladeshi breast cancers occur before age 30. That is roughly **170 to 334 women a year**. Carrying the small-sample noise in the underlying studies as well widens it to **1.3%–4.7%**, or 126 to 447 women.
>
> Among women aged 20–29, that is an annual incidence on the order of **1 to 3 per 100,000**.
>
> This is a model-based estimate, not an observation. It moves by a factor of four across the plausible range of one assumption.

A range that wide is not a failure. It is the correct width. Anyone handing you a tighter number for this quantity has not done better statistics — they have made a stronger assumption and not told you what it was.

**Rare and hundreds at the same time.** For an individual woman aged 20–29, that rate is roughly a **1 in 50,000** chance this year — rare, by any ordinary use of the word. The same estimate says **170–334 women** will hear this diagnosis before their 30th birthday this year. *Rare* is a rate. *Hundreds* is a count. Confusing the two is how a 28-year-old with a lump gets reassured instead of imaged.

**Diagnostics:** R̂ = 1.001 on θ across three chains · largest leave-one-study-out shift 0.15 pp · posterior referral tilt OR 2.26 (1.50–3.42).

---

## The model

```
x_j        ~ Binomial(n_j, q_j)                 # Bangladeshi hospital series
logit(q_j) = logit(θ) + δ + u_j                 # referral tilt + study noise
u_j        ~ Normal(0, τ²)

logit(θ)   ~ Normal(m_regional, s_regional²)    # what regional registries see
δ          ~ Normal(μ_δ, σ_δ²)                  # tilt learned where BOTH a
                                                #   registry and hospital
                                                #   series exist
τ          ~ HalfNormal(0.5)

C_<30      = θ × C_all,   C_all ~ LogNormal around GLOBOCAN's 9,480
```

`θ` is the national quantity. The data identify `logit(θ) + δ` and nothing finer; the split comes entirely from the priors. **That is the point, not a defect** — and it is why the headline is a range.

The last line matters too. GLOBOCAN's 9,480 is a borrowed estimate with no published interval. Multiplying a posterior share by it as though it were exact manufactures precision out of nothing, so it gets a distribution and the count is computed draw by draw.

The sampler is a forty-line random-walk Metropolis written out in base R rather than a call to Stan. Deliberately: when the whole argument is *be transparent about where the answer comes from*, it helps to be able to read the machine that produces it.

---

## Repository layout

```
.
├── R/
│   ├── 00_inputs.R            every number in one place, each tagged
│   ├── 01_naive_pooling.R     pooling + the precision-is-not-accuracy simulation
│   ├── 02_bounds.R            partial identification: the θ-vs-δ curve
│   ├── 03_bayes_synthesis.R   Metropolis sampler, model ladder, LOO, sensitivity
│   ├── 04_figures.R           four figures + all result CSVs
│   └── run_all.R              runs the lot
├── figures/                   four PNGs
├── results/                   eight CSVs, every number the post quotes
├── make_tables_docx.py        builds tables_for_post.docx from the CSVs
├── tables_for_post.docx       nine formatted tables
├── medium_post.md             the post
└── README.md
```

---

## Honesty about the inputs

Every input in `00_inputs.R` carries a tag: `[VERIFIED]`, `[APPROX]`, `[ILLUSTRATIVE]` or `[ASSUMPTION]`. The most important one:

> **The published Bangladeshi hospital studies report mean ages and `<40` / `<50` cutoffs. None of them publishes a count under 30.**

A mean age of 46.2 is compatible with an under-30 share of 1% or of 10% — infinitely many age distributions share a mean. So the under-30 counts in `00_inputs.R` are **constructed** to be consistent with what those papers do report. They are tagged `[ILLUSTRATIVE]`, and every downstream number in this repository inherits that tag.

| | |
|---|---|
| **Verified** | GLOBOCAN 2024 Bangladesh figures (9,480 new cases; 4,211 deaths; 23,625 five-year prevalence; female population 88,218,713) and their footnotes; study sizes for three of the six series |
| **Approximate, marked** | Population denominator for women 20–29; comparator registry shares; the referral-tilt prior |
| **Illustrative** | The under-30 counts |

**Swap in real counts and every number moves. None of the conclusions do.** The argument is about structure — that pooling non-representative sources produces confident wrong answers, that the bias parameter is unidentified from the data alone, that the honest output is a range, and that the fix is a registry rather than a cleverer model.

---

## What would actually fix this

Everything here is a bridge, not a road. The thing that ends this exercise is a **population-based cancer registry** — defined catchment, mandatory notification, active case-finding, age published as counts in five-year bands.

Three things would help immediately, and two cost nothing:

1. **Publish counts, not means.** Every case series that reported "mean age 46.2" instead of an age table destroyed most of its own information on the way to the printer.
2. **Describe your catchment.** Say who could have reached your hospital and who could not. That sentence is what lets a future analyst estimate δ instead of guessing it.
3. **Report the assumption, not just the interval.** If the estimate moves by a factor of four across defensible values of one parameter, that fact *is* the finding.

---

## The line that generalises

> **A model can borrow information and quantify assumptions. It cannot manufacture representativeness.**

When you have a non-random sample and no external anchor, the honest output is not a point estimate with a small standard error. It is a range, plus an explicit statement of what you would have to believe to land anywhere in it. The best estimate is not the one with the tightest interval — it is the one whose target, observation process, transport assumptions and uncertainty are all written down where somebody can disagree with them.

---

## Methodological background

- Bayesian multiparameter evidence synthesis — Sweeting et al., hepatitis C prevalence in England and Wales
- Multilevel regression and poststratification — Gelman & Little
- Doubly robust inference from non-probability samples — Chen, Li & Wu
- Bayesian multistate disease-burden modelling — Jackson et al.

**Primary source:** Global Cancer Observatory (IARC / WHO), Bangladesh fact sheet, GLOBOCAN 2024 — <https://gco.iarc.who.int>

---

## Citation

```bibtex
@misc{nahin2026registry,
  author = {Kazi Sabbir Ahmad Nahin},
  title  = {Estimating a National Disease Burden With No National Survey},
  year   = {2026},
  url    = {https://github.com/<you>/<repo>}
}
```

---

## License

MIT. See [LICENSE](LICENSE).
