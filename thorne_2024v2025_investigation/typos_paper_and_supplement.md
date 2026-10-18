# Typos and inconsistencies: main paper and supplement (17 Oct 2026 drafts)

## Main paper

### Numbers and internal consistency (most important)

- **Kal_flexLin weight, Section 5.1.** It is given as 7.9 % in the weights list, but as 5.4 % in "(5.4%, 7.4%, and 10.5% respectively)". The 5.4 % is that method's crossing probability, not its weight.
- **Uncertainty reduction.** The text says "approximately 3-fold" (Section 5.1) but "approximately 5-fold" (Section 6, point 3).
- **Footnote 39.** "RMSE of 0.489, 0.418, 0.450" should read 0.0489, 0.0418, 0.0450.
- **Table 3, butterworth row.** The SSP2-4.5 RMSE is "0.00597", presumably 0.0597.
- **Table 3 header vs caption.** The header says "RMSE 2025-2100"; the caption says "(2025-2090)".
- **"as of the end of 2024, only 4 of the 56 evaluated methods exceed 1.5°C" (Section 4.6.1).** The year should probably be 2025 for this update. The method count is 71 analysed / 40 retained elsewhere.
- **Year ranges that don't match:**
  - The text says "historical period (1850-2024)", but the Table 2 caption says 1850–2025.
  - The Figure 27 caption says "observational syntheses from 1850-2024".
- **Rate range.** "0.32 [0.21 - 0.43] °C/decade": the README gives 0.212677–0.424970, which rounds to [0.21–0.42].
- **Section numbering in Tables 2–4.** The tables label "4.4.2: Constrained Emulators" and "4.4.3: Energy Balance Kalman Filter", but the text has 4.4.2 = EBM-KF and 4.4.3 = emulators.
- **Figure 22 caption.** It refers to both "ensemble member 20" and "The 11th ensemble member is thicker". Check which member is highlighted.
- **Figure 28 caption.** "distinguished from their counterparts in Figure 26" should be Figure 27.
- **Section 5.2 cross-references.** "(Section 4.6.3)" for volcanic eruptions should be 4.6.2. "(Section 4.6.2)" for reaching/exceeding/returning should also be checked.

### Spelling and wording

- "1850-1900 to 2013-2012" should be "2003–2012".
- Affiliation ³² appears twice; the first, Stripe Inc., should be ³¹.
- "EMBKF_ta4" should be "EBMKF_ta4" (Section 5.1 weights list).
- "Foster et al., 2025" and "Foster et al 2025" should be "Forster" (Figure 4 caption, twice).
- "wee chose" should be "we chose" (end of Section 4.7).
- "has B* tag because it hasare highlighted" (Table 3 caption) should be "have B* tags because they are highlighted".
- "Mastraneda et al., 2010" should be "Mastrandrea".
- "Klimaseniorennen" should be "KlimaSeniorinnen" (the footnote spells it correctly).
- "Hasselman (1979)" should be "Hasselmann".
- "multimodel ensembleover" should be "ensemble over" (Figure 18 caption).
- "ensemble membersto" should be "members to"; "monthly upates" should be "updates" (Section 5 summary).
- "exceednace" should be "exceedance" (Figure 29 caption).
- "RSME" should be "RMSE" (Section 4.6.3).
- "(see Supplement 3).." has a doubled full stop.
- "Smith and Reynolds, 2005[3]": check the footnote placement.

## Supplement

- **Affiliation 3:** "Geneva, SwitzerlandIndependent researcher, France" is a merged or leftover line.
- **Figure S11 caption** says "addressed in v5.2", but the paper uses HadCRUT5.1 throughout. Check the version.
- **Table S2:**
  - "Milllar" should be "Millar"; "Vissier 2018" should be "Visser 2018".
  - "FaIR_antho" / "FaIR_antho_unB" should be "anthro".
  - "GWI_antrho_CGWL" and "eROF_antrho" should be "anthro".
  - "except for for attributable" has a doubled word.
  - Section labels such as "4.3/2_ERF_FaIR" and "4.3/3_Kalman" don't match the main-text sections (4.4.x).
  - "$(Fig 14)" / "$(Fig 15)" / "$(Fig 17)" figure numbers look stale.
- **Section 2.3.1:**
  - "NorESM Voc ensemble" (twice) should be "Volc".
  - "GWI_tot_anthro" should probably be "GWI_anthro_orig".
  - "4.4.6: Human Induced" doesn't exist; it should be 4.4.4.
- **Supplement text:**
  - "entire 1850-20100 window" should be 2100.
  - "Kullback-Liebler" should be "Kullback-Leibler" (Table S4 header).
  - "von Shuckmann" should be "von Schuckmann".
  - "none of the datasets space the whole" should be "span".
  - "EMB-KF-ta4" should be "EBM-KF-ta4".
  - "Figure S12.Tthe" should be "Figure S12. The".
  - "observations are sparse .It" should be "sparse. It".
- **Section 2.2:** "The main paper in Figure 20 illustrated … 1℃" should be Figure 21. "Figure S1 and Figure 20" should be Figure 21.
- **Table S3 caption:** "The red shading" refers to shading lost in plain export; check the formatting.
- **Table S4 caption:** it describes a comparison with GWI_tot_CGWL and EBMKF_ta4, but those rows aren't in the table.
- **References:**
  - "Scalfe" should be "Scaife" (Betts et al.).
  - "Zscheisler" should be "Zscheischler".
  - Otto et al. 2015: "Myles, R." should be "Allen, M. R.".
  - Samset et al. 2023: "Zelinka, D." should be "Zelinka, M. D.".
