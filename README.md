<div align="center">

# 🧠 Longitudinal Modeling of Motor Progression in Parkinson's Disease and Serum NfL Using Bayesian Multilevel Models

**Does baseline serum neurofilament light chain (NfL) predict how fast motor symptoms worsen in Parkinson's disease?**

![Python](https://img.shields.io/badge/Python-3.14-3776AB?logo=python&logoColor=white)
![PyMC](https://img.shields.io/badge/PyMC-Bayesian%20Modeling-blue)
![Sampler](https://img.shields.io/badge/Sampler-NUTS%20(HMC)-orange)
![Data](https://img.shields.io/badge/Data-PPMI-green)
![Convergence](https://img.shields.io/badge/R--hat-%E2%89%A4%201.01-brightgreen)
![Paper](https://img.shields.io/badge/Paper-Interspeech%20template-lightgrey)

*Ricardo Vidal Cuadrado · Universidad San Pablo CEU, Madrid, Spain*

</div>

---

## 📑 Table of Contents

1. [TL ; DR](#-tldr)
2. [Background](#-background)
3. [Data](#-data)
4. [Model](#-model)
5. [Sampling Setup](#-sampling-setup)
6. [Results](#-results)
7. [Interpretation](#-interpretation)
8. [Repository Structure](#-repository-structure)
9. [Getting Started](#-getting-started)
10. [Reproducibility](#-reproducibility)
11. [Limitations](#-limitations)
12. [Citation](#-citation)
13. [Acknowledgments](#-acknowledgments)
14. [Generative AI Use Disclosure](#-generative-ai-use-disclosure)

---

## ⚡ TL;DR

- 🎯 **Goal:** quantify whether **baseline serum NfL** is associated with faster motor deterioration in Parkinson's disease (PD), measured with the **MDS-UPDRS Part III** scale.
- 📊 **Data:** the **Parkinson's Progression Marker Initiative (PPMI)** cohort, restricted to patients with **≥ 3 clinic visits** (**13,120 longitudinal observations**).
- 🧮 **Method:** a **Bayesian linear mixed-effects model** (random intercept + random slope, correlated via a Cholesky-parameterized LKJ prior) fitted with **PyMC** and the **NUTS** sampler.
- ✅ **Result:** the **NfL × time interaction** ($\beta_5$) is **positive** with a 95% HDI that excludes zero: **0.308 [0.097, 0.515]**. Patients with higher baseline NfL show **faster motor progression**.
- 🔁 **Convergence:** $\hat{R}$ equal to 1.00 or 1.01 for all fixed effects, with high ESS.

---

## 🔬 Background

Parkinson's disease is a clinically heterogeneous neurodegenerative disorder with **motor** (tremor, bradykinesia, postural instability) and **non-motor** symptoms (loss of smell, constipation, depression, blood-pressure changes). Motor severity is usually assessed with the **MDS-UPDRS Part III**, a clinician-rated scale that is inherently subjective, which complicates **prognosis** and **clinical trial design**.

**Neurofilament light chain (NfL)** is a structural protein of large-caliber myelinated axons. It is released into the CSF and eventually reaches the bloodstream when axonal damage occurs, which makes **blood NfL** a minimally invasive marker of neurodegeneration. NfL is *not* specific to PD (it also rises in Alzheimer's disease, FTD, and ALS), but it is promising as a **prognostic** marker of disease progression.

Because age, sex and genetics can also shift NfL levels, the model adjusts for **age** and **sex**, and the data pipeline also considers **APOE ε4** carrier status and **GBA** mutations.

<img width="695" height="300" alt="nfl_csf_plasma" src="https://github.com/user-attachments/assets/ee5dd82d-90c4-4dca-9b29-f71aa72a813f" />
<p align="center">
  <img src="docs/figures/nfl_csf_plasma.png"  width="70%">
  <br>
  <em>Figure 1. Comparison of plasma and CSF NfL concentrations across experimental groups, illustrating their parallel increase and strong correlation.</em>
</p>

---

## 🗂 Data

| Item | Description |
|---|---|
| **Source** | [Parkinson's Progression Marker Initiative (PPMI)](https://www.ppmi-info.org/) |
| **Design** | Longitudinal, observational, multi-visit cohort |
| **Outcome** | MDS-UPDRS Part III motor score |
| **Key predictor** | Baseline serum NfL (log₂-transformed, mean-centered) |
| **Covariates** | Age (mean-centered), sex, APOE ε4 carrier status (binary), GBA mutation status |
| **Inclusion criterion** | ≥ 3 clinic visits per patient |
| **Final sample** | **13,120** longitudinal observations |

### Preprocessing

- **Log₂ transform** of NfL, because its distribution is strongly right-skewed.
- **Mean-centering** of continuous variables (age and log₂ baseline NfL), which makes the units comparable and improves MCMC sampling efficiency.
- **Binarization** of APOE ε4 status (carrier vs. non-carrier).
- **Time** $t_{ij}$ is measured in **years since the first visit**.

> ⚠️ **PPMI data is not redistributed in this repository.** Access must be requested directly from PPMI under their data use agreement.

<!--
📍 FIGURE 2 — INSERT HERE
File in paper: spaghetti_plot_hito1_real.png
Caption: Spaghetti plot of individual MDS-UPDRS Part III trajectories against age
at clinical visit in the PPMI cohort.
-->
<p align="center">
  <img src="docs/figures/spaghetti_plot_hito1_real.png" alt="Spaghetti plot of MDS-UPDRS III trajectories" width="70%">
  <br>
  <em>Figure 2. Spaghetti plot of individual MDS-UPDRS Part III trajectories against age at clinical visit in the PPMI cohort.</em>
</p>

The plot shows how noisy and heterogeneous patient trajectories are, which is exactly why a **random-effects** formulation is appropriate.

---

## 🧮 Model

A **Bayesian hierarchical (mixed-effects) linear model**. For patient *i* at visit *j*:

$$
y_{ij} \sim \mathcal{N}(\mu_{ij}, \sigma^2)
$$

$$
\mu_{ij} = (\beta_0 + u_{0i}) + (\beta_1 + u_{1i})\,t_{ij} + \beta_2\,\text{Age}_i + \beta_3\,\text{Sex}_i + \beta_4\,\text{NfL}_{i,0} + \beta_5\,(t_{ij} \times \text{NfL}_{i,0})
$$

| Parameter | Meaning |
|---|---|
| $\beta_0$ | Population intercept (baseline motor score) |
| $\beta_1$ | Population-level progression rate (UPDRS points / year) |
| $\beta_2$ | Effect of age |
| $\beta_3$ | Effect of sex |
| $\beta_4$ | Effect of baseline NfL (log₂, centered) on baseline motor score |
| **$\beta_5$** | **NfL × time interaction: how much baseline NfL accelerates motor decline** ⭐ |
| $u_{0i}, u_{1i}$ | Patient-specific random intercept and random slope (correlated) |

### Priors and regularization

- **Fixed effects** $\beta_1 \dots \beta_5$: weakly informative $\mathcal{N}(0, 5)$.
- **Intercept** $\beta_0$: centered at **20.0**.
- **Random effects covariance:** **Cholesky parameterization** with an **LKJ prior** ($\eta = 2.0$), which regularizes the covariance matrix and helps avoid convergence problems.

---

## ⚙️ Sampling Setup

| Setting | Value |
|---|---|
| Library | **PyMC** |
| Sampler | **NUTS** (No-U-Turn Sampler, a Hamiltonian Monte Carlo variant) |
| Chains | 4 independent chains |
| Warm-up (tuning) | 1,500 iterations per chain |
| Sampling draws | 2,000 per chain |
| Total posterior samples | **8,000** |
| Target acceptance rate | **0.95** (to avoid divergent transitions) |

Supporting libraries: **Pandas** and **NumPy** for data handling and transformation.

---

## 📈 Results

### ✅ Convergence diagnostics

All fixed-effect parameters converged: $\hat{R}$ values of **1.00–1.01**, high **bulk and tail ESS**, and low **MCSE**.

<!--
📍 TABLE 1 — INSERT HERE
Table in paper: "MCMC Diagnostics for Fixed Effects Parameters"
(Reproduced below in Markdown. Alternatively, replace with a screenshot, e.g. docs/figures/table1_mcmc_diagnostics.png)
-->

**Table 1. MCMC diagnostics for fixed-effects parameters**

| Parameter | | Mean | SD | 95% HDI | ESS | $\hat{R}$ |
|---|---|---:|---:|---:|---:|---:|
| $\beta_0$ | Intercept | 10.78 | 0.44 | [10.0, 11.0] | 640 | 1.01 |
| $\beta_1$ | Time | 1.136 | 0.093 | [0.99, 1.30] | 1141 | 1.00 |
| $\beta_2$ | Age | −0.066 | 0.039 | [−0.13, −0.005] | 373 | 1.01 |
| $\beta_3$ | Sex | 1.67 | 0.57 | [0.76, 2.60] | 619 | 1.00 |
| $\beta_4$ | NfL₀ | 1.32 | 0.48 | [0.55, 2.10] | 480 | 1.01 |
| **$\beta_5$** | **Interaction (NfL × time)** | **0.308** | **0.107** | **[0.097, 0.515]** | **1147** | **1.00** |

### 📉 Motor progression trajectories

<!--
📍 FIGURE 3 — INSERT HERE
File in paper: bayesian_stratification_panel.png
Caption: Posterior predictive motor progression trajectories by baseline NfL
percentiles (A) and 95% HDI forest plot of fixed-effect parameters (B).
-->
<p align="center">
  <img src="docs/figures/bayesian_stratification_panel.png" alt="Bayesian stratification panel" width="95%">
  <br>
  <em>Figure 3. Posterior predictive motor progression trajectories by baseline NfL percentiles (A) and 95% HDI forest plot of fixed-effect parameters (B).</em>
</p>

### 🎯 The key parameter: $\beta_5$

<!--
📍 FIGURE 4 — INSERT HERE
File in paper: beta5_posterior_real.png
Caption: Posterior distribution of the interaction parameter (β5) confirming the
acceleration of motor decline driven by baseline NfL.
-->
<p align="center">
  <img src="docs/figures/beta5_posterior_real.png" alt="Posterior distribution of beta5" width="70%">
  <br>
  <em>Figure 4. Posterior distribution of the interaction parameter (β₅) confirming the acceleration of motor decline driven by baseline NfL.</em>
</p>

| Quantity | Estimate |
|---|---|
| Posterior mean of $\beta_5$ | **0.308** |
| 95% HDI | **[0.097, 0.515]** (entirely positive, excludes 0) |
| Mean progression rate $\beta_1$ | **1.136** UPDRS-III points per year |

---

## 💡 Interpretation

- The **95% HDI of $\beta_5$ is entirely above zero**, so the data support that a **higher baseline NfL is associated with a steeper increase in MDS-UPDRS Part III scores over time**.
- Patients with baseline NfL above the cohort mean follow a **steeper motor decline curve** than the average patient, over and above the average yearly progression ($\beta_1$).
- The model **adjusts for age and sex**, so the NfL effect is not simply reflecting those factors.
- **Clinical relevance:**
  - 🧑‍⚕️ **Patient stratification:** identifying who is likely to progress faster, enabling better planning for clinicians and patients.
  - 🧪 **Clinical trials:** better separation of treatment and control arms, reducing the masking effect of PD heterogeneity.
  - 🩸 **Practicality:** blood sampling is far less invasive than lumbar puncture, and blood and CSF NfL are strongly correlated.

---

## 📁 Repository Structure

> 🛠️ *Adjust this tree to match the actual contents of the repository.*

```text
nfl_biomarker_mcmc/
├── data/                  # PPMI data goes here (NOT included, see "Data")
├── notebooks/             # Data loading, filtering, modeling and analysis
├── src/                   # (optional) reusable Python modules
├── docs/
│   └── figures/           # Figures used in the paper and in this README
│       ├── nfl_csf_plasma.png
│       ├── spaghetti_plot_hito1_real.png
│       ├── bayesian_stratification_panel.png
│       └── beta5_posterior_real.png
├── requirements.txt
└── README.md
```

---

## 🚀 Getting Started

```bash
# 1. Clone the repository
git clone https://github.com/ricvidaal/nfl_biomarker_mcmc.git
cd nfl_biomarker_mcmc

# 2. (Recommended) create a virtual environment
python -m venv .venv
source .venv/bin/activate        # Windows: .venv\Scripts\activate

# 3. Install dependencies
pip install -r requirements.txt
```

**Main dependencies:** `pymc`, `pandas`, `numpy`, and the usual plotting / diagnostics stack (e.g. `arviz`, `matplotlib`).

**Run the analysis:**

1. Request and download the data from [PPMI](https://www.ppmi-info.org/) and place it in `data/`.
2. Run the data loading and filtering step (≥ 3 visits per patient).
3. Fit the Bayesian model with PyMC (NUTS, 4 chains, 1,500 warm-up + 2,000 draws, `target_accept=0.95`).
4. Inspect diagnostics ($\hat{R}$, ESS, MCSE) and generate the figures.

---

## 🔁 Reproducibility

- Fixed model specification, priors and sampler configuration are documented above.
- Convergence is verified with $\hat{R}$, bulk/tail ESS and MCSE.
- Setting a **random seed** in the sampling call is recommended to reproduce exact draws.
- The raw PPMI data must be obtained independently, since it cannot be redistributed.

---

## ⚠️ Limitations

- 🔗 **Association, not causation:** the model shows that baseline NfL is associated with faster motor decline, not that it causes it.
- 🧬 NfL is **not disease-specific**: it also rises in Alzheimer's disease, FTD and ALS, so it is not a standalone diagnostic tool for PD.
- 📏 The MDS-UPDRS Part III outcome keeps some **inter-rater subjectivity**, even though NfL itself is an objective measurement.
- 📐 The model is a **linear** progression model; non-linear trajectories are not captured.

---

## 📚 Citation

If you use this work, please cite it (update with the final reference once published):

```bibtex
@misc{vidalcuadrado_nfl_pd,
  author       = {Vidal Cuadrado, Ricardo},
  title        = {Longitudinal Modeling of Motor Progression in Parkinson's Disease
                  and Serum NfL Using Bayesian Multilevel Models},
  howpublished = {\url{https://github.com/ricvidaal/nfl_biomarker_mcmc}},
  note         = {Universidad San Pablo CEU, Madrid, Spain}
}
```

---

## 🙏 Acknowledgments

- To the **Parkinson's Progression Marker Initiative (PPMI)** for providing access to the dataset, as well as information about the initiative, how the data were collected, and their scope.
- To my project advisor, **Juan Balsa**, for his availability and valuable assistance.

---

## 🤖 Generative AI Use Disclosure

- **Open Evidence:** used to obtain summarized and accessible information from the medical literature and to answer specific questions about parts of the project.
- **Gemini:** used to translate specific phrases into English and to paraphrase certain parts of the paper to improve grammar and spelling.
- **Claude:** used to generate some parts of the code and to obtain information about certain PyMC functions.

---

<div align="center">

⭐ If you find this project useful, consider starring the repository.

</div>
