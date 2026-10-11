# Kryptos K4: crib-constrained eliminations not (fully) covered elsewhere

Author: Thanos Zekios. The computations were carried out with the help of Claude (Anthropic).

This note documents a set of exact searches on K4 that are **absent or only partly covered** in the three public elimination efforts I compared against:
- jcolinpatrick/kryptos (KryptosBot)
- nitsuah/kryptos
- keithadler/kryptos

Every section was checked against those three repositories. Tests they already cover with the same or larger scope have been removed. Where they cover a family partly, a **Prior work** line says what exists and what is new here. All scripts are in `scripts/`.

## Conventions

- **Ciphertext.** K4 = `OBKRUOXOGHULBSOLIFBBWFLRVQQPRNGKSSOTWTQSJQSSEKZZWATJKLUDIAWINFBNYPVTTMZFPKWGDKZXTJCDIGKUHUAUEKCAR` (97 letters).
- **Cribs.** 0-based positions 21–33 = `EASTNORTHEAST` and 63–73 = `BERLINCLOCK` (24 letters). In 1-based terms these are 22–34 and 64–74.
- **Alphabets.** `ST` = A–Z; `KA` = `KRYPTOSABCDEFGHIJLMNQUVWXZ`. Rules tested are Vigenère (c = p + k), Beaufort (c = k − p) and Variant (c = p − k), in index arithmetic of the alphabet.
- **Every test has a positive control.** A planted ciphertext is built with the cribs at the same positions, and the search must find it.
- **Most tests also have a null control.** The same search is run on shuffled K4, to calibrate the false-positive rate.
- **What "excluded" means.** It means "no configuration in the stated space is consistent with all 24 crib letters". It does not exclude neighbouring spaces.

## 1. The five W's as segment separators (new)

K4 has W at 1-based positions 21, 37, 49, 59, 75. Of these:
- position **21** sits immediately before `EASTNORTHEAST` (22–34);
- position **75** sits immediately after `BERLINCLOCK` (64–74).

Segments between W's have lengths 20 · 15 · 11 · 9 · 15 · 22. `EASTNORTHEAST` starts its segment, so 2 unknown letters follow it before the next W. `BERLINCLOCK` ends its segment, so 4 unknown letters precede it.

**Tests** (`wseg.py`, `wseg2.py`, `nulls.py`, `depth.py`, `guess.py`)

| Model | Result |
|---|---|
| Key restarts after every W, any key (even non-repeating), ST/KA Vig/Beau. Crib letters at equal in-segment offsets must share a key value (9 overlaps). | **Excluded:** agreement 1–2/9. |
| Same, arbitrary alphabet per in-segment offset ("segments in depth"). | **Not excluded.** 0 contradictions, 1 strong check. |
| W's are inserted nulls (removed, 92 letters): periodic key 1–24 with arbitrary alphabets, and columnar widths 2–8 (all orders) plus a key of length ≤12. (ST/KA without transposition is covered by nitsuah's null-gap check.) | Excluded (one threshold-level "any-alphabet" hit, at the same rate as the null). |
| Every other letter treated as a null class (direct test). | Excluded. |
| Segment lengths used as row lengths of a transposition (read by columns) plus a periodic key. | Excluded. |
| 23 two-letter words after `EASTNORTHEAST` × ~560 four-letter fillers before `BERLINCLOCK`, added as extra cribs under the W-restart key (jcolinpatrick's filler tests use one global key). | No second coincidence. |

**Observation, not a result.** Under W-segmentation, the R of NORTHEAST (pos 28) and the R of BERLIN (pos 66) fall on the **same in-segment offset (7th letter)**, and both encrypt to **P**. This is the known k[27] = k[65] "Bean equality", seen from a new angle.
- With W positions fixed and the other letters shuffled, the event "0 contradictions and ≥ 1 strong match" occurs with p ≈ 0.026 (20,000 shuffles).
- W was chosen post hoc, so this is not evidence.
- Two consequences under this hypothesis: a 4-letter word before `BERLIN` cannot be `EAST` (FLRV vs INFB), and the letters after `EASTNORTHEAST` cannot start with C or have K second.
- Related anomaly reported by KryptosBot: the width-21 bigram anomaly disappears when the W's are deleted.

## 2. Light / mask over the tableau or panel selecting the KEY (new)

The hypothesis: a mask or the light passing through holes selects key letters from the Vigenère tableau, or from the K1–K3 panel, as carved. The selected letters are the key. The K4 letters themselves stay in place, as the position-aligned cribs require.

| Mask family (`mask.py`, `mask2.py`, `yar.py`, `sinemask.py`) | Space | Result |
|---|---|---|
| n holes per row (n = 4–32) at fixed columns, all start rows, offsets, straight or serpentine; tableau and panel, top-down and bottom-up | exact per-slot column feasibility | 0 (null 0) |
| Mask sliding by α columns per row over the tableau (α = 0–25) | as above | 0 |
| "YAR-shaped" holes (relative columns 0,1,3), tiled every q = 4–32 columns | exact | best 7/24 (null 8) |
| Sinusoidal (S-curve) path: amplitude 1–6 rows, wavelength 8–97, all phases, row drift, start rows | exact | best 9/24 (null 9) |

Regular masks over the tableau reduce to odometer/Alberti keys (section 6), and those are excluded too. An irregular, hand-made mask is not testable. **Prior work:** jcolinpatrick (E-ANTIPODES-06, E-BESPOKE-11) and nitsuah (tableau keystream sweeps) score single-path readings of the tableau; keithadler reads the tableau and panel row by row as a running key. The multi-hole, YAR-shaped and sinusoidal masks with an exact crib test are new.

## 3. Hill with a periodic key or a transposition (extension)

`hillx.py`, `hillT.py`, `hillchk.py`. Exact linear algebra mod 2 and mod 13 over **all** matrices and block alignments, ST and KA indexing.

**Prior work.** Direct Hill 2×2–4×4 and affine Hill are excluded by keithadler, nitsuah and jcolinpatrick (DOT-STANDING-001/002/010/014); Hill 2×2 followed by keyed columnar widths 2–10 is jcolinpatrick's DOT-STANDING-030; columnar followed by Hill is partly covered by jcolinpatrick E-S-41 (widths 5–8), E-S-151 (width 9) and E-ANTIPODES-01. New here:
- **"Periodic key then Hill"** (block-class offsets, ≤ 6 classes): inconsistent for 2×2 and 3×3.
- **Columnar → Hill at widths 2–4**, and KA indexing at all widths 2–9: 0 hits.
- **Hill 3×3 → columnar** at widths 2–9: 0 hits at widths 2–8. Width 9 gives hits with only 3 equations for 2 unknowns; most are non-invertible matrices, and all 276 invertible ones decrypt to gibberish (quadgram ≈ −7.5/char).

## 4. Paradigm-style three-layer keyword pipelines (new)

`layers3.py`, `layers3b.py`, `big3.py`, `biglist.py`. These are exact crib checks for the PK puzzle structures:
- T(w1)→Q(w2)→Q(w3), as in PK4;
- T(w1)→T(w2)→Q(w3), as in PK6;
- Q→Q→T;
- Q→Q→Q, as in PK8.

Here T is keyword columnar (and its inverse) and Q is Quagmire III, each step ±. Both KA and ST were run.

- **Controls:** with only the 24 crib positions, the search recovers PK4 (`UNDERLAY`→`OCHRE`→`VERDIGRIS`) and PK6 (`HANDIWORK`→`SMITHWORK`→`PORTAL`), including inside the 363-word list.
- **K4, 89-word and 363-word thematic lists** (Sanborn works, CIA/Cold War, Berlin, Egypt, crypto history, magnetism/geodesy, light/time, materials, K0–K3 words, all PK keys; see `biglist.py`): **0 hits in every structure.**

## 5. Any-alphabet periodic key with transposition (extension)

`anyalpha.py`, `gen.py`. This is the most general periodic substitution: P completely independent alphabets, so it covers Quagmire I–IV with any two alphabets. Within a key class, same plaintext letter ⇔ same ciphertext letter.

- **Direct, fully independent alphabets per key class:** periods 1–7, 9, 10, 14, 15, 17, 21 and 25 are inconsistent. All other periods give no check at all and are undecidable this way. (With one shared mixed alphabet, i.e. general Quagmire I–IV, keithadler's exact solver covers periods up to 26; at 19 Quagmire IV stays satisfiable because of the R→P repeat 38 = 2×19 apart.)
- **Per crib segment alone:** `EASTNORTHEAST` excludes periods 1–5 and 9; `BERLINCLOCK` excludes 1–4 and 7.
- **With ~1.1M transpositions, periods 1–12, key indexed by plaintext or ciphertext position:** 0 hits (null 0). The transpositions are:
  - columnar widths 2–9, all orders, both directions;
  - "24-face cylinder" readings (24/12 columns, straight or helical, any start/step/direction);
  - opposite-direction rotations (two halves, even/odd, grid rows/columns as gears).
- **K3-style double rotation,** all widths 2–96 × 2–96, CW/CCW, both directions (72,200 permutations): 1 threshold hit (null also 1). The implementation reproduces K3 exactly (336/336).

**Prior work.** keithadler (columnar 2–14, K3 double rotation), nitsuah (columnar 2–14, 21,096 rotation layouts, 7,680 geometric maps) and jcolinpatrick (route bank) test these transpositions with ST/KA tableau keys. New here is the arbitrary-alphabet version at periods ≤ 12, and the cylinder / opposite-rotation families.

## 6. Mechanisms not found elsewhere

The rotor is new. The odometer extends jcolinpatrick's additive grid key (E-S-10, ST only); the Alberti stepping with N ≠ p extends the progressive keys of nitsuah and keithadler (N = p). Trifid is covered by keithadler and jcolinpatrick and is not repeated here.

| Mechanism (script) | Space | Result |
|---|---|---|
| Single stepping rotor (Hebern type), **unknown wiring**, step every N = 1–48, both directions, ST/KA (`rotor.py`) | exact injectivity check | 0 (null 0/100) |
| Odometer key k = a[i mod p] + b[⌊i/p⌋ mod q] (`odometer.py`) | p + q ≤ 20 | 0 |
| Alberti-disk stepping k = a[i mod p] + d·⌊i/N⌋ (`odometer.py`) | p ≤ 18, N = 2–48 | 0 |
| 6×6 Bifid (`frac.py`) | periods 2–97, all offsets for ≤ 24 (jcolinpatrick E-S-09 covers 2–50; 51–97 are new) | infeasible everywhere |

## 7. Error tolerance under transposition (extension)

`err1.py`. The models are: one crib letter wrong (drop each anchor in turn), one extra letter (97 deletions) and one missing letter (98 insertions). Tested:
- periodic key 1–24 with arbitrary alphabets (keithadler covers ST/KA with 1 wrong letter or 1–3 inserted/dropped letters);
- **columnar widths 2–8 with periods ≤ 12, drop-one-anchor model, arbitrary alphabets** (nitsuah has a statistical version with ST/KA families).

Result: 0 hits. Controls recovered a planted typo and a planted extra letter.

## 8. Numeric "masking" keys (extension)

`numkeys.py`, `numsrc.py`, `coordcol.py`. Sequences used as a non-repeating shift, at every start offset, ×1/−1/2/3, add or subtract, ST/KA:
- digits of π, e, √2, √3 and φ as 2-digit blocks and sliding 2-digit values (single digits are covered by keithadler's any-0–9-keystream search);
- squares, cubes, triangular numbers, primes and prime gaps, powers of 2/3/5/7, Lucas;
- K2 coordinates and dates (9 Nov 1989; 3 Nov 1990);
- Morse code of the sculpture text (code lengths, dots, dashes, binary values, 5-bit blocks);
- Weltzeituhr: UTC offset per city or per letter, and name lengths;
- W segment lengths and positions.

About 120,000 configurations gave best 6–7/24, at the null level. **Prior work:** jcolinpatrick E-S-50 (Fibonacci, powers, π/e at one alignment, ST) and E-WTZ-00 (UTC offsets, sampled). Coordinate and date digits used as column orders are already covered by the all-orders columnar searches of keithadler and nitsuah.

## 9. W-segment restart of keyword alphabets

Keyed alphabets from `WILLIAMWEBSTER`, `WILLIAMHWEBSTER`, `WILLIAMHEDGCOCKWEBSTER`, `JUDGEWEBSTER`, `ORDINATE`, `COORDINATE(S)`, `APPLICATE`, `YAXIS`, `ABSCISSAORDINATE`, the misspelling letters (`QUA`, `LOE`, …) and `YAR`/`RAY`, with the key **restarting at every W segment**: all at chance level (≤ 6/24). The global (non-restarting) versions are covered by keithadler's any-alphabet Quagmire solver and nitsuah's dictionary-alphabet search.

## Caveats

- Everything above assumes the published cribs are exact and position-aligned. Section 7 relaxes this by one character.
- The arbitrary-alphabet tests decide only the periods listed in section 5 when there is no transposition, and nothing above 24.
- Three-layer results only exclude the listed words.

## Running

The scripts need Python 3 and numpy. `pk_controls.json` holds the PK4 and PK6 plaintexts, used only as controls. `quad.pkl` contains English quadgram log-probabilities. Each script prints its control result first, then the K4 result.

## 10. Period-19 Quagmire IV: see section 35

An earlier annealer failed its planted control, so its runs carried no information and it has been removed. Section 35 replaces it with a working search (planted control recovered).

## 11. Weltzeituhr (Alexanderplatz) face transcription (2017 state) and the 1990 problem

The per-face, per-line city order is listed as a blocker in nitsuah's open list ("per-letter Weltzeituhr lookup keys"). Below is my transcription from user photos (Google Maps contributions, Feb/Jul 2017, plus later screenshots).

Notes on the format:
- Order is top to bottom as engraved.
- "top" = above the rotating hour ring, "bottom" = below it.
- `(blank)` = map only, no city names.
- `?` = not readable.

| UTC | Top | Bottom |
|---|---|---|
| −12 / −11 / −10 / −9 | ? (small text, not readable; includes Alaska/Hawaii names such as NOME, FAIRBANKS, ANCHORAGE, HONOLULU per other sources) | ? |
| −8 | VANCOUVER · DAWSON · SAN FRANCISCO · LOS ANGELES | (blank) |
| −7 | EDMONTON · DENVER | (blank) |
| −6 | NEW ORLEANS · MEXIKO-STADT | GUATEMALA-STADT · MANAGUA · GALÁPAGOS I. |
| −5 | MONTREAL · WASHINGTON · NEW YORK · HAVANNA | PANAMA · SANTAFÉ DE BOGOTÁ · QUITO · LIMA |
| −4 | HALIFAX | CARACAS · LA PAZ · ASUNCIÓN · SANTIAGO DE CHILE |
| −3 | WESTGRÖNLAND | BRASILIA · RIO DE JANEIRO · SÃO PAULO · MONTEVIDEO · BUENOS AIRES |
| −2 | OSTGRÖNLAND | ? |
| −1 | AZOREN | KAP VERDE |
| 0 | REYKJAVIK · DUBLIN · LONDON · LISSABON · MADEIRA · BISSAU | CASABLANCA · CONAKRY · DAKAR · BAMAKO · ACCRA |
| +1 | AMSTERDAM · BERLIN · BRÜSSEL · BUDAPEST · MADRID · PARIS · PRAG · STOCKHOLM · WARSCHAU | OSLO · KOPENHAGEN · WIEN · BERN · PRESSBURG · BELGRAD · ROM · TUNIS · KINSHASA |
| +2 | HELSINKI · RIGA · TALLINN · WILNA · MINSK · KIEW · BUKAREST · SOFIA · NIKOSIA | ANKARA · ISTANBUL · ATHEN · TEL AVIV · JERUSALEM · BEIRUT · DAMASKUS · KAIRO · KAPSTADT |
| +3 | MURMANSK · ST. PETERSBURG · MOSKAU | TEHERAN +30′ · BAGDAD · ADEN · SANAA · ADDIS ABEBA · MOGADISCHU · DARESSALAM · ANTANANARIVO · KUWAIT |
| +4 | NISCHNIJ NOWGOROD · WOLGOGRAD · BAKU · TIFLIS · ERIWAN | KABUL +30′ · MAURITIUS |
| +5 | JEKATERINBURG · ASCHGABAT · BISCHKEK · DUSCHANBE | NEW DELHI +30′ · KARACHI · COLOMBO +30′ |
| +6 | OMSK · ALMATY · TASCHKENT · NOWOSIBIRSK | RANGUN +30′ · DHAKA |
| +7 | KRASNOJARSK | HANOI · BANGKOK · PHNOM PENH · JAKARTA |
| +8 | IRKUTSK · ULAN-BATOR | PEKING · SHANGHAI · MANILA · PERTH · HONGKONG · KUALA LUMPUR · SINGAPUR |
| +9 | JAKUTSK | PJÖNGJANG · TOKYO · SEOUL |
| +10 | CHABAROWSK · WLADIWOSTOK | SYDNEY · CANBERRA · MELBOURNE |
| +11 | MAGADAN · SACHALIN | (blank) |
| +12 | KAMTSCHATKA | DATUMSGRENZE · WELLINGTON |

Confidence varies by face:
- **Highest** for −6…+5, from clear photos.
- **Lower** for the order within +6…+10, which comes from less sharp images.

**The 1990 problem.** The clock was renovated in 1997. Leningrad became St. Petersburg, Alma Ata became Almaty and Bratislava became Pressburg, about 20 cities were added (e.g. Tel Aviv, Jerusalem), and some zone assignments changed (e.g. Kiew). Sources: de.wikipedia "Weltzeituhr (Alexanderplatz)"; the-passenger.de.

Kryptos was dedicated in 1990. Any "BERLIN CLOCK → city list" hypothesis therefore needs the **pre-1997** faces, which I could not find in readable form online. Good leads:
- Bundesarchiv / Wikimedia Commons GDR photos;
- 1970s–80s postcards (eBay, Delcampe);
- Berlinische Galerie.

Faces showing LENINGRAD or ALMA ATA identify the old state.

**Prior work.** jcolinpatrick's `e_s_134_weltzeituhr_ddr_era.py` already holds per-face city lists from Hungarian Wikipedia (including the −12…−9 faces) and an inferred DDR-era list (Leningrad, Alma-Ata, Gorki, Bogota; Kiew at +3); nitsuah holds 130 photo-read names. New here: the top/bottom split, the order within each side, the photo provenance and the +30′ marks. City names as keyed alphabets with repeating keys are covered by keithadler's any-alphabet solver; only multi-word face strings and Quagmire I/II at lengths keithadler leaves open were run here (`cities.py`, `cities_eu.py`): no hits.

## 12. Sub → transposition → sub with two independent periodic keys; width-21 routes

**Model** (`stackT.py`):
- c[m(i)] − e·p[i] = K1[i mod a] + K2[m(i) mod b];
- the first key is indexed by plaintext position and the second by ciphertext position, after transposition m;
- e = ±1 covers all Vigenère/Beaufort combinations;
- the check is exact, by weighted union-find mod 26.

This is a fully exhaustive version, within these ranges, of the "substitution + transposition + substitution" class. **Prior work:** jcolinpatrick's E-FRAC-52 covers widths 6/8/9 with Bean-filtered orders, ST, periods ≤ 12, and is marked disputed in its own `docs/elimination_tiers.md`. New here: widths 2–5 and 7 with all orders, KA, no Bean pre-filter.

**Results:**
- Columnar widths 2–7, all column orders and both directions, a = 2–8, b = 2–8, ST and KA: **0 hits** (null 0).
- Control: a = 8 (ABSCISSA), width 5, b = 7 (PALIMPS) is found uniquely.

**Width 21** (`w21.py`). KryptosBot reports 11 repeated vertical bigrams at width 21 (p ≈ 1.6e-4; row-major specific). Their audit shows the effect is length-fragile rather than W-specific. I tested 6,312 width-21 transpositions:
- row-major with every cyclic column rotation;
- column reads (down, up, right-to-left);
- serpentine by rows and by columns;
- both diagonals;
- spiral and reverse spiral;
- keyword columnar from all 21-letter one- or two-word phrases of the 363-word list.

Each was combined with any periodic key (period ≤ 21, any alphabet, ST/KA) and with two stacked periodic keys (a, b ≤ 8). **0 hits** (null 0). A planted spiral + period-7 control was recovered.

## 13. Running key under transposition: exact search (Carter's book and an 11M-letter corpus)

**Prior work.** keithadler tests Carter vols. 1–3 as running keys without transposition (149M offsets); jcolinpatrick's `e_s_31` tests Carter vol. 1 with columnar widths 5–10, all orders, A–Z, Vigenère/Beaufort; jcolinpatrick's FRAC-49/50 used a Bean filter that its own `docs/research_questions.md` now discounts; nitsuah has a statistical (fragment-scored) version for any English text with columnar 2–8. New here: widths 2–4, the KA alphabet, the Variant rule, reversed text, both key-index models, the K3-rotation / cylinder / opposite-rotation families, and an exact 24/24 test on a named 11M-letter corpus. Check (`rk_trans.py`):
- for every key offset in the text and every transposition, all 24 crib letters must match;
- the key is indexed by ciphertext position (substitution after transposition) or by plaintext position (substitution before);
- rules: Vigenère, Beaufort or Variant; alphabets KA and ST.

**Howard Carter, *The Tomb of Tut-ankh-Amen* Vol. I** (the source of K3; archive.org full text, 287,513 letters, forward and reversed; `carter_run.py`):

| Transposition family | Size | KA | ST |
|---|---|---|---|
| K3-style double rotation, all widths 2–96 | 72,200 | 0 | 0 |
| 24/12-face cylinder readings | 88,992 | 0 | 0 |
| opposite-direction rotations | 196,126 | 0 | 0 |
| columnar widths 2–9, all orders, both directions | 818,496 | 0 | 0 |

Control: a key taken from the "SLOWLY DESPERATELY SLOWLY" passage plus a K3 rotation was recovered at the exact offset.

**Corpus of 11.16M letters:** NLTK Gutenberg (KJV Bible, Moby Dick, Shakespeare, Milton, Whitman, Blake, Austen, Chesterton, Carroll…), all US inaugural addresses, State of the Union addresses 1945–2006, and K1–K3. Result: columnar widths 2–8, all orders and both directions, KA and ST: **0 hits**. Control: a Moby Dick offset plus width-6 columnar was recovered.

## 14. PK9-style pipelines with spiral routes; chemistry and music keys

**Q(word) → spiral → T(word)** (`spiral3.py`):
- spirals on grids of width 2–48 (ragged last row), 4 start corners × 2 rotation senses, used as either read or write order;
- combined with keyword columnar (363 words, plus inverse and none) in both orders, and a Quagmire III keyword (363 words, ±);
- KA and ST.

Control: PK9, using only 24 plaintext letters at the K4 crib positions, recovers CLEPSYDRA → spiral(12, TR, D-L-U-R) → BEAMWORK. K4: **0 hits**.

Spirals (752) were also tested with any periodic key (period ≤ 24, any alphabet) and with two stacked keys (a, b ≤ 10): 0 hits. KryptosBot's RQ-3 lists spiral readings as untested.

**Chemistry and music** (`chem_music.py`) were used as numeric keys, at every offset, ×±1, ST/KA, Vig/Beau:
- atomic masses, atomic numbers and groups;
- piano-key frequencies (A4 = 440 Hz), C-major MIDI numbers, Westminster chime, the BACH motif, chromatic semitones.

Best match 6/24 (null 5). Element-symbol and element-name streams as running keys, direct or with columnar widths 2–8: 0.

Curiosity: BERLINCLOCK is spellable in element symbols (B·Er·Li·N·Cl·O·C·K). EASTNORTHEAST and K4 are not, since K4 contains Q and J.

## 15. Dictionary-scale PK-style search; periodic-table layouts

**PK9 structure with a 370,105-word dictionary** (`dictsearch.py`; dwyl/english-words, which contains every PK keyword including CLEPSYDRA and BEAMWORK):
- Q(any dictionary word, any length, ±) applied first, so the key is plaintext-indexed (exact hash lookup);
- then none, or any of 752 spirals (widths 2–48, 4 corners × 2 senses, read or write);
- then none, or columnar widths 2–8 with **all** column orders, both directions, in either order relative to the spiral;
- KA and ST.

Control: PK9 is recovered from 24 plaintext letters (CLEPSYDRA, spiral 12 TR, width-8 order). **K4: 0 hits.**

**PK4 structure, any keys** (`stack2ct.py`): transposition, then Q(key a), then Q(key b), with both keys ciphertext-indexed. Exact union-find for any keys of lengths a < b ≤ 12 (b not a multiple of a), columnar widths 1–7, all orders, ±, ST/KA.

Control: T5 → OCHRE → VERDIGRIS is recovered. **K4: 0 hits.** At width 7 the shuffled null gave 6 threshold-level hits; K4 gave 0.

**Periodic table as layout** (`ptable.py`). In 1988 IUPAC recommended the 1–18 group numbering, published in the 1990 Red Book, the year Kryptos was dedicated.
- Routes: K4 letters placed on elements Z = s … s+96 (s = 1–22) in the 18-column and 32-column layouts, then read by rows or columns, forward or reverse (352 permutations). Tested with any periodic key (period ≤ 24, any alphabet), two stacked keys (a, b ≤ 10), and dictionary Q + columnar 2–5: **0**.
- Masks: the table shape laid over the tableau or the K1–K3 panel at every offset, with key letters read in atomic-number order (59,466 configurations): best 5/24, at chance level.

## 16. Gear trains, safe-dial transpositions and clock-hand keys

All of these are exact linear tests mod 2 and mod 13 (`mech.py`, `mech_run.py`, `mech_gen.py`, `mech_err.py`).

**Gear trains as key.** The key is the sum of two gears, k_i = g1[i mod p1] + g2[i mod p2], for every pair 2 ≤ p1 < p2 ≤ 30 (406 pairs). Gears with tooth counts taken from the raised YAR and from the misspelling letters (Q/U/A/E, L/O/E/I and their shifts, in ST and KA, 0- and 1-based) are included as 2–4-gear trains. Vigenère/Beaufort, ST/KA.
- Direct: 2 consistent cases at redundancy 1, against 4.8 expected by chance; the shuffled null gave 4. **Nothing.** Most pairs with p1+p2 > 25 have too many unknowns to be tested from 24 cribs. (Pairs up to 26 are covered by keithadler §4F and nitsuah's two-key check; new here are the YAR/misspelling multi-gear trains.)
- Control: a planted (4,15) gear train is recovered.

**Clock hands as key.** k_i = a·i + b·⌊(i+s)/r⌋ + c, i.e. a minute hand plus an hour hand stepping every r letters. All r = 2–60 and every phase s; also a three-hand version (r1, r2) ∈ {12×12, 24×12, 12×5, 60×12, 5×12}. This covers a progressive "dial sweep" key (b = 0). For r ≤ 26 the model is a progressive key, already excluded by nitsuah and keithadler; new here are r = 27–60, the three-hand version and the one-error tolerance.
- Each model has redundancy 21.
- Direct: **0**.
- One-error tolerance (any one of the 24 crib letters wrong, for all models with redundancy ≥ 5): **0**.
- Control: a planted r = 12, s = 5 case is recovered.

**Safe-dial transposition.** The 97 letters sit on a dial; a pointer alternates direction (L, R, L, …), moving n1, n2, n3, … positions among the remaining letters and removing the letter it lands on (an alternating-direction Josephus). The removal order is used as the transposition, in both directions.
- Combinations: 1 number (1–97), 2 numbers (1–97 each), 3 numbers (1–26 each), 4 numbers (1–12 each), first turn L or R. That is 191,272 permutations including inverses.
- A literal-combination set adds 3,072 permutations: YAR and misspelling-letter values in every order, with 0–3 extra full turns (+26 each).
- On top of each permutation: clock-hand models (redundancy ≥ 6, 535 models), gear trains, and any periodic key with period ≤ 12, ST/KA or any alphabet. Columnar widths 2–8 (all orders) were also run against the clock/gear models.
- **Results: 0** for clock/gear models (0.019 expected by chance) and **0** for periodic keys on the 191k set.
- The literal-combination set gave one any-alphabet hit at the minimum threshold: combination (86, 82, 89, 92) = L-O-E-I + 3 turns, ct-indexed period 8, 4 coincidences. The measured chance rate is 6 hits per 1.53M permutation tests on shuffled K4, so over the ~194k permutations run here about 0.8 hits were expected. The hit survives Quagmire IV (too unconstrained to matter) but not Quagmire III or ST/KA alphabets. **We read it as chance.**
- Controls: a dial (7,3,11) + KA Vigenère period 7 is recovered; a dial (5,9,3) + clock key is recovered.

## 17. KRYPTOS as a safe combination (`kdial.py`, `kdial_c.py`)

The word KRYPTOS turned into dial numbers in six ways:
- ST 0-based, ST 1-based, KA 1-based;
- ST distance between consecutive letters, signed and shortest;
- tableau header;

each forward and reversed, with 0–3 uniform extra turns or "safe" decreasing turns.

Three uses were tested:
- **(a) Dial transposition** (step-wise and absolute stops, both directions, 432 permutations), with periodic/any-alphabet keys and clock/gear models on top: **0**.
- **(b) Dial-sweep keystream.** A pointer clicks around a KA or ST ring to each letter of KRYPTOS (or reversed), alternating direction, with uniform or decreasing turn counts, every start and offset. Each click gives one key letter. Vigenère/Beaufort/variant, KA/ST: 299,520 keystreams. Best 6/24, at chance.
- **(c) Null mask.** Cumulative dial stops on the 97-letter ring mark 1–14 letters as nulls; only masks that leave the cribs intact are kept. Each reduced text is tested with any periodic key ≤ 24, ST/KA or any alphabet. 740 masks: **0**. A planted 7-null + KRYPTOS-period-7 control is recovered.

## 18. Simplified Hagelin M-209 / C-38 (`hagelin.py`, `hagelin2.py`)

**Why test it.** Scheidt's CIA years coincide with the CIA's covert ownership of Crypto AG. The M-209 pin wheels have 26, 25, 23, 21, 19 and 17 positions, and 19 is the length at which keithadler's general Quagmire IV solver stays satisfiable.

**Model.** k_i = c0 + Σ_j w_j·pin_j[(i+off_j) mod p_j], with lug weights and no lug overlap. ST alphabet; Beaufort (the M-209 rule, up to the constant c0) and Vigenère. The full M-209, with 131 free pins, 27 bars and lug overlap, cannot be decided from 24 cribs.

**Part 1 — any two wheels, every pin pattern free, any weights 1–25, any c0.** All 15 wheel pairs, solved exactly with 2-SAT. **K4: 0**, shuffled null: 0. A planted (19,17) case is recovered.

**Part 2 — pins fixed by a letter set on the real wheel lettering, offsets enumerated, weights solved linearly mod 2 and mod 13.**
- Letter sets: YAR (the raised letters), QUAE, LOEI, QUAE+LOEI, the misspelling pairs, KRYPTOS; each also complemented.
- Wheels: every 1- and 2-wheel subset, and 3-wheel subsets containing 19 and 17.
- 910,464 configurations. **K4: 0**, shuffled null: 0. A planted YAR-pinned (23,19) case is recovered.

**Part 3 — the triple 25 + 21 + 17** (`hagelin3.py`). These are the 1-based values of Y, U and Q, which come from YAR and from the misspellings, and are also M-209 wheel sizes. Pins fixed by the same letter sets, plus complements; all 8,925 offset combinations; Beaufort and Vigenère: 214,200 configurations. **K4: 0**, shuffled null: 0. A planted QUAE-pinned case is recovered.

## 19. "Three turns of the dial": K1, K2 and K3's methods chained (`chain3.py`, `chain3b.py`)

**Prior work.** jcolinpatrick's `e_s_143_progressive_onion` chains K1–K3-style layers but models K3 as a width-7 columnar. New here: the true K3 double rotation at all widths 2–96, key phase continued from K1/K2, and both application directions.

**Hypothesis.** Each earlier section is one turn of the dial, applied in its own way. K4 = some chain of:
- K1's method: KA Vigenère, PALIMPSEST;
- K2's method: KA Vigenère, ABSCISSA;
- K3's method: double rotation, any widths 2–96, both senses, forward or inverse.

**Search.** All 6 orders; each substitution in either direction (left/right); key phase 0 or continued from where K1 and K2 end (63 mod 10, 369 mod 8). In total 72,200 transpositions × 96 chains, scored by crib matches.
- Control: a planted ABSCISSA → rotation(13,7) → PALIMPSEST cipher scores 24/24.
- **K4: best 8/24**, the same as the shuffled null (8/24): chance level.

**Section sizes as the combination.** The dial numbers are taken from the sections:
- lengths (63, 369, 336), also with the '?' (372, 337) and with 97 appended;
- key lengths (10, 8), alone or with the K3 grid widths;
- raw, mod 97 and mod 26, in every order, as an alternating-direction dial transposition (488 permutations).

Tested with any periodic key ≤ 24 (ST/KA/any alphabet) and the clock/gear models: **0**.

## 20. Greek hypotheses (`greek.py`, `greek_astro.py`)

**Words.** 236 transliterated Greek words used as keys. Topics: hidden/truth (KRYPTOS, ALETHEIA, APOKALYPSIS), light/shadow (PHOS, SKIA), clocks (GNOMON, HOROLOGION, KLEPSYDRA), mechanisms (ANTIKYTHERA, MECHANE, TROCHOS), ciphers (SKYTALE, POLYBIOS), magnetism (MAGNES, MAGNESIA), Egypt/tombs, the eight winds of the Tower of the Winds in Athens (a clock with compass directions, incl. KAIKIAS = NE), letter names, numbers and gods.

**Base alphabets.**
- ST, KA;
- Greek-keyboard order ABGDEZHUIKLMNJOPRSTYFXCVQW;
- an alternative transliteration order;
- KRYPTOS-keyed Greek order.

**Pipelines (validated on PK9, PK4, PK6).**
- Q(word) → none or 752 spirals → columnar 2–7, all orders;
- T+Q+Q, Q+Q+T, Q+Q+Q, T+T+Q with keyword columnar.

Result: **0** in every alphabet.

**Long periods 27–60 and period 38 (= 2×19, also an Antikythera gear count; periods 25–26 are covered by keithadler and nitsuah).** Not consistent for ST/KA/Greek alphabets with ≥ 3 checks. With arbitrary alphabets period 38 has a single check, the R→P pair at 28/66, so it remains undecidable.

**Calendar-gear keys.** k_i = a·i + b·⌊(p·i+s)/q⌋ + c, every phase s, five alphabets, Vigenère/Beaufort (22,760 configurations). Ratios used:
- Metonic 235/19 and 19/235, sidereal 254/19;
- Saros 223/18, draconic 242/223, anomalistic 239/223;
- Callippic 940/76, exeligmos 669/54, Olympiad 1/4;
- six Antikythera gear pairs.

Result: **0**. A planted Metonic key is recovered.

## 21. Clock-state masks and clock vocabulary (`clocks.py`, `clockwords.py`)

**Prior art.** The Mengenlehreuhr is already in jcolinpatrick/kryptos as band routes, a Cardan-grille reading order and keyword tests (E-S-122, E-MENGEN-GRILLE-01, `k4_clock_cipher.py`), and in nitsuah/kryptos as bearings. What follows adds exact tests over every displayed time with free keys.

**Masks (11,613 states).**
- Mengenlehreuhr: all 1,440 minutes × seconds lamp on/off, as 24-lamp vectors (1+4+4+11+4), in both reading directions.
- "Fill-bar" level clocks (water-clock style, as in the Europa-Center *Uhr der fließenden Zeit*): bars of 12, 24 or 60 vessels at every fill level.
- DCF77 59-bit minute frames for every minute of 3 Nov 1990, 9 Nov 1989, 3 Oct 1990 and 3 Nov 1988.

**Results.**
- **(1) Lamps as key:** k = c + Σ_row w_row·lamp, with any weight per lamp row and every phase; ST/KA × Vigenère/Beaufort. For the 24-lamp Mengenlehreuhr and 12/24 bars this is a special case of the period-24/12 keys that keithadler already excludes; the new part is the 59-bit DCF77 frames and the 60-vessel bars. **K4: 0**, null: 0. A planted control is recovered.
- **(2) Lamps as grille:** read lit cells first, then unlit (per block and globally, plus inverses); 45,903 distinct permutations. Tested with any periodic key ≤ 24 (ST/KA/any alphabet) and the clock/gear models: **0**. (jcolinpatrick's E-S-122 turns the 1,440 lamp patterns into permutations with fixed keyword keys; new here are free periodic keys, arbitrary alphabets and the DCF77/fill-bar grilles.)
- **(3) Lamps as null mask:** lit or unlit cells, tiled at every phase, removed; only the 503 distinct masks that leave the cribs intact are kept. The reduced text is tested with any periodic key ≤ 24: **0**.

**Clock vocabulary as keywords** (PK9's CLEPSYDRA is a clock word). 207 English and German horology words: escapement, fusee, verge, Hemmung, Unruh, Zifferblatt, Mengenlehreuhr, Weltzeituhr, Wasseruhr, fließende Zeit, DCF77, etc. Same pipelines and five base alphabets as §20 (Q → spiral → columnar 2–7; T+Q+Q, Q+Q+T, Q+Q+Q, T+T+Q): **0**.

## 22. Safe wheel-pack keystream (`safe.py`)

**Mechanism.** Dialling a 3-wheel combination lock moves all 3 wheels together in the first phase, then 2, then 1, alternating direction. If each click enciphers one letter, the key is piecewise linear in position with slopes −3, +2, −1 (or the mirror +3, −2, +1).

**Test.**
- Generalised to any 3-phase slope sequence from {±1, ±2, ±3} (216 sequences) and any two phase boundaries.
- Key indexed by ciphertext position, on top of 355,937 transpositions: identity, columnar 2–8, alternating dial transpositions, K3-style double rotations (both directions).
- ST/KA × Vigenère/Beaufort.
- Exact: each phase must have a constant intercept, giving 21 redundant checks.

Results: **K4: 0**, shuffled null: 0. A planted dial(15,25,2) + (−3, +2, −1) keystream is recovered.

**Literal lock simulation** (`safe2.py`).
- **Model.** The dial is ringed with ST or KA letters. Four phases: wheels {1,2,3} turn T1 full turns then stop at a; {2,3} reverse direction, T2 turns, stop at b; {3} reverse, T3 turns, stop at c; then the cam-only opening stroke.
- **Keystream.** One letter per click; the key is the sum of the wheel positions (or the dial position), ± sign, with a free constant.
- **Settings.** Turn patterns (3,2,1), (4,3,2), (2,1,0), (1,1,1), (0,0,0), (1,2,3); both starting directions; every window offset of the 97 letters into the click sequence; ST/KA × Vigenère/Beaufort.
- **Exhaustive, anchors decide:** all 26³ stop combinations, dial starting at 0. 611M window tests: **K4 0**.
- **Stops at the earlier sections' misspellings, read K3 → K2 → K1:** wrong letters A,U,Q; correct letters E,O,L; the K1 → K2 → K3 order; the misspelling positions in the plaintexts (10, 114, 56); shift sizes. Every dial start 0–25; 7.2M window tests: **K4 0**, null 0.
- **Stops at YAR / RAY** (letters and 1-based numbers), every dial start: 3.6M window tests: **K4 0**, null 0.
- A planted A-U-Q (3,2,1) lock keystream is recovered.

## 23. The clue words themselves, and guessed Berlin words as extra cribs (`cribwords.py`, `cribwords2.py`, `gapwords.py`)

**(a) Clue words as keys.** EASTNORTHEAST, BERLINCLOCK, BERLIN, CLOCK, EAST, NORTHEAST, ENE and their concatenations, plus 70 Berlin/ENE-bearing words, run through the PK pipelines in five base alphabets: **0**.

**(b) Berlin words as extra cribs at every free position.** LENIN, LENINPLATZ, LENINALLEE, MOLLSTRASSE, FERNSEHTURM, WELTZEITUHR, HAUSDESREISENS, SBAHN, WALL, CHECKPOINT and others; 3,558 placements.
- Against a period-19 key with arbitrary alphabets they are consistent at the same rate as random strings (0.2%). The extra letters do not make period 19 decidable.

**(c) Train/rail/station/direction words, with or without THE/A/AN/DER/DIE/DAS, in the gap 34–62** (224 phrases, 4,994 placements): no difference from random strings for period 19.

**Prior work.** jcolinpatrick's `e_s_139`/`e_s_153` place Berlin words at every position against periods 2–26 (A–Z); keithadler §5 slides 52 theme words against long keys. New here: the period-19 arbitrary-alphabet test calibrated against random strings, and these word lists.


## 24. Weltzeituhr mechanism (`wcity.py`, `wring.py`)

**Source.** Sanborn's November 2025 open letter (Scientific American): "The BERLINCLOCK mentioned in previous K4 clues refers to the World Clock in Berlin."

**(A) The city drum as running key.** The face texts of §11 are read as the drum turns:
- east or west, from any starting face;
- top+bottom, bottom+top, top only or bottom only;
- umlauts stripped or written as OE/UE/AE;
- with or without guessed Pacific faces;
- the 2017 names and a 1990 variant (LENINGRAD, ALMAATA, BRATISLAVA, GORKI; Tel Aviv and Jerusalem removed);
- each also letter-reversed.

That gives 128 streams of about 500 letters, at every offset, ST/KA × Vigenère/Beaufort/variant.
- Direct: **best 8/24** (null 6/24).
- With columnar 2–7, all orders, key on ciphertext positions: **best 10/24**; the shuffled null also gives 10/24, about what chance predicts.
- A planted drum-text key is recovered at 24/24. A 1990-correct transcription remains the main caveat.

**Prior work.** jcolinpatrick's `e_s_128`/`e_s_134` use concatenated city lists (modern and DDR-era) as running keys, with columnar width 7. New here: the 128 stream variants, widths 2–6, KA, the Variant rule and both key-index models. Hour-ring and city-count keys are special cases of period-24 keys and are not repeated.


## 25. Coupling the two Berlin clocks (`twoclocks.py`)

**Mechanism coupling.** Both clocks have 24 elements: 24 lamps and 24 faces. For each Mengenlehreuhr state (1,440 minutes × seconds lamp), lamp j is mapped onto a Weltzeituhr face:
- any starting face, either direction;
- either the lit or the unlit lamps select faces.

The city text of the selected faces, in drum order (top+bottom, top or bottom; 2017 or 1990 names), is then used as a running key: every offset, ST/KA × Vigenère/Beaufort/variant.
- 3.3 × 10⁹ window tests: **K4 best 8/24**, shuffled null 8/24, i.e. chance.
- A planted 14:37 coupling key is recovered at 24/24.

## 26. PK9 structure with the two Berlin clocks as the route (`pk9clock.py`, `pk9clock8.py`)

**Prior work.** jcolinpatrick's `e_berlin_clock_perm.py` and E-S-127 read clock-shaped routes with one repeating or keyword key. New here: the three-layer stack with all 370k dictionary Q keys, helical decimations mod 97 and all 120×32 band orders.

**Idea.** PK9 = Q(CLEPSYDRA, a clock word) → spiral route → T(BEAMWORK). Here the spiral is replaced by clock geometry.

**Routes (4,750 distinct, read and write).**
- Weltzeituhr cylinder: rows of 24 faces, read face by face from any start face, east or west, up or down.
- Helical readings: every decimation k mod 97; all starts for k = 23, 24, 25.
- Mengenlehreuhr band routes on 24-lamp blocks: all 120 band orders × 32 band reversals.

**Pipeline.** Q(any of the 370,105 dictionary words, plaintext-indexed) → clock route → columnar transposition with widths 2–8, all orders (width 8 covers every 8-letter keyword such as BEAMWORK). Route and columnar in either order; KA and ST.

**Result.** A planted Q(CLEPSYDRA) → cylinder → width-6 cipher is recovered. **K4: 0** for widths 2–7 and for width 8, in both alphabets.

## 27. Every ST/KA-only family re-run with 270 sculpture-keyed alphabets (`sculptalpha.py`, `sculptalpha2.py`, `allalpha.py`)

**Alphabets.** Keyword-mixed alphabets built from every word (≥ 3 letters) of the K1–K3 plaintexts, the Morse panel, the keys (PALIMPSEST, ABSCISSA, KRYPTOS), the clues (BERLIN, CLOCK, EASTNORTHEAST, …), the misspelled words (IQLUSION, UNDERGRUUND, DESPARATLY, DIGETAL) and phrases such as LAYERTWO, IDBYROWS, WHATSTHEPOINT. Each word is followed by the rest of the alphabet in normal or reverse order; 268 distinct alphabets, plus ST and KA.

**Not re-run.** Repeating keys with a shared mixed alphabet (general Quagmire I–IV) are covered by keithadler's exact solver for any alphabet, and with transpositions by our arbitrary-alphabet tests at periods ≤ 12.

**Re-run, direct.** Each family has a planted PALIMPSEST-alphabet control that is recovered.
- Q(370,105 dictionary words) → none / 752 spirals, and Q → columnar 2–7 (all orders): **0**.
- Linear models (3,779): gear pairs 2–30, clock hands r = 2–60, three-hand clocks, Greek/astronomical ratios, Weltzeituhr hour-ring alignment: **0**.
- Safe wheel-pack piecewise keys (216 slope sequences, any phase boundaries): **0**.
- Clock-state key bits (Mengenlehreuhr, fill bars, DCF77; 482k configurations × 540 alphabet/rule variants): **0**.
- Hill and affine Hill 2×2 to 5×5: **0**.
- Plaintext and ciphertext autokey, primers 1–40: **0**.
- Exact running key over 11.66M letters (Gutenberg, inaugural, State of the Union, Carter vol. I, K1–K3, 128 Weltzeituhr drum streams), forward and reversed, key letters read in the mixed or the standard alphabet: **0**.
- Shuffled-K4 nulls: 0 everywhere.

**PK-style three layers** (T+Q+Q, Q+Q+T, Q+Q+Q, T+T+Q with the 363-word list) on 16 base alphabets: PALIMPSEST, ABSCISSA, IQLUSION, UNDERGRUUND, DESPARATLY, BERLINCLOCK, EASTNORTHEAST, SHADOW, LUCID, DIGETAL, INTERPRETATU, LANGLEY, MAGNETIC, CANDLE, LAYERTWO, IDBYROWS. Result: **0** for all 16.

## 28. Checking two published proposals (numberworld.blog, parts 4 and 5) (`scott.py`, `scott2.py`, `xor5.py`)

**"Scott's method"** (part 4): Q3(KRYPTOS / EMUFPHZLRFA, the first 11 K1 letters), then the position map y = 77 + 38x (mod 97), with numbers taken from K2's coordinates.
- Re-run independently, both directions and both signs.
- BERLIN does appear, but at 1-based position **41**, not 64, and it is followed by **PEEKR…, not CLOCK**.
- EASTNORTHEAST is absent; the best crib score is 2/24.
- The post's three-step HXDKZEK extension scores 0–1/24.
- Over all 855 eleven-letter windows of the sculpture ciphertext × both signs × all 9,312 affine maps, "BERLIN somewhere" occurs exactly once, at about 0.05 expected. That is mildly lucky but post-hoc, and the method is falsified by the 2014 CLOCK clue.

**5-bit / XOR proposal** (part 5; OBKR as key, CIAW/CIAX message key). The post itself reports that it fails. We tested 5-bit XOR with any periodic key ≤ 24:
- ST/KA, 0- or 1-based codes, with or without the b1↔b2 bit swap;
- **0** of 192 configurations.

## 29. Teleprinter / computer-era XOR family (`xor5b.py`)

Encodings: 5-bit A=0/A=1 in ST or KA order, ITA2 (Baudot–Murray teleprinter code, as used on Vernam and Lorenz tape machines), and bit-reversed ITA2.
- **Periodic XOR key ≤ 24:** **0**.
- **Plaintext/ciphertext autokey XOR, primers 1–40:** **0**.
- **Running-key XOR** (key = letters of a text): impossible for every encoding. Some crib pairs XOR to codes that are not letters; for example the S→S self-encryption gives 0, which is not a letter in A=1 or ITA2.
- **LFSR keystream (Vernam with a shift-register generator).** Berlekamp–Massey linear complexity of the crib keystream bits:
  - K4: 26–35 for the 65-bit and 55-bit blocks;
  - shuffled null: 25–33;
  - a random sequence of n bits has complexity about n/2.
  - So no LFSR of degree below ~25 generates the key. A planted degree-17 LFSR is recovered with complexity 17.

## 30. Punched-card grille, 3×3×3 fractionation with transposition, period 19 with keyword alphabets (`hollerith.py`, `frac3*.py`, `p19keyed.py`)

**Hollerith punched card as a grille.** (jcolinpatrick's `team/e_hollerith_*` scripts use the card to remove nulls with named keywords; here it is a reading order with free keys.)
- Mask: 572 words (sculpture vocabulary, clock words, Greek words, dates) punched in IBM/Hollerith code (12 rows: 12, 11, 0–9), tiled at every shift.
- Layout: K4 laid down 12-row card columns or across 9-wide rows. Letters under the holes are read first, then the rest (and the inverse): 15,964 permutations.
- Tested with any periodic key ≤ 24 (ST/KA/any alphabet) and the clock/gear models: **0**, null 0. A planted KRYPTOS-card + KA period-7 control is recovered.

**ADFGVX-style fractionation on a 3×3×3 cube** (the only "square" that can output exactly 26 letters plus one spare cell): letters → trits → columnar transposition (widths 2–7, all orders, both directions) → regroup.
- Cribs force coordinate equalities.
- Read direction: 5,299 of 5,912 orders are eliminated by a forced collision (null: 3,741).
- Write direction: only 2 are eliminated; the rest cannot be decided from 24 cribs.
- Random cubes give letter-only output only ~7% of the time (the spare cell appears otherwise), so the family is a priori unlikely.
- In 20 planted letter-only cases the true order always survives.

**Period 19 and 38 with keyword-mixed alphabets** (direct): Quagmire I/II over 370,091 dictionary + sculpture words (reversed-ST base and the ~138k words outside nitsuah's 231,933-word list) and Quagmire IV over all 268 × 268 sculpture alphabets (reversed tails, non-dictionary words such as IQLUSION and LAYERTWO), plus period 38: **0**. Quagmire III at period 19 is already impossible for any alphabet (keithadler `results/14b`).

## 31. What has NOT been tried yet (cross-checked against jcolinpatrick, nitsuah, keithadler and this package)

**Ground rule.** Anything below can only be confirmed by the cribs if it is constrained enough. Several items are "undecidable from 24 letters" rather than "untried"; they need either a named key rule or K5.

**A. Transposition families never searched exhaustively with a key**
1. **Turning grilles (Fleissner).** A physical light/mask device of exactly the kind Scheidt described ("masking"); a 10×10 grille covers 97 letters + 3 nulls. jcolinpatrick (E-S-18, E-S-70, `e_fleissner_*`) and §32 searched heuristically; §36 adds an exact search that has so far covered about 42% of one of eight variants (0 survivors). The rest is still open.
2. **Myszkowski, disrupted and incomplete columnar**, and keyed routes. Now covered exactly in §36 for Myszkowski widths 3–9, VIC-disrupted widths 5–9 (10 sampled) and irregular columnar widths 2–8. Still open: Myszkowski ≥ 10, disrupted ≥ 10 exhaustively, keyed routes.
3. **Columnar widths ≥ 15 with arbitrary column orders.** Keyword-ordered widths 12–30 (keithadler) and clock-keyword widths 15–26 (nitsuah) are done. Any transposition + period-19 arbitrary alphabets is undecidable without a named alphabet.

**B. Masking that changes length**
4. Letters inserted or dropped inside words, mirroring the K1–K3 misspellings (IQLUSION, UNDERGRUUND, DESPARATLY). Tested so far: nulls between crib blocks, many 24-null masks (jcolinpatrick), clock-derived null masks, a crib-alignment diagnostic (nitsuah), and in §36 a shift of ±1…10 between the two crib blocks under periodic, transposed, running and linear keys. Still open: edits combined with ciphertext-indexed keys under transposition, and shifts larger than 10.
5. Carving errors **combined with a transposition**, or two or more errors. Without transposition, keithadler already allows 1 wrong crib letter or 1–3 inserted/dropped letters.

**C. A method that changes part-way.** Sanborn on the German submission: "the first part had started to look right, but then the rest didn't."
6. Two different transpositions or two different key systems for two segments, with a split point anywhere. Only fixed splits (at the W's, or with two periodic keys) have been tested.

**D. Key sources tied to the two events Sanborn named (Nov 2025: the 1986 Egypt trip and the 1989 fall of the Wall)**
7. As running keys **under transposition**: Carter vols. II–III (keithadler tested them without transposition only); Sanborn's own writings about the Egypt trip; Schabowski's 9 Nov 1989 press-conference transcript. (Reagan 1987 and Kennedy 1963 have already been run under several transposition families by jcolinpatrick.) Without transposition, any English running key is already excluded statistically.
8. The K0 Morse panel and the Cyrillic Projector as running sources under transposition (listed open by nitsuah).
9. **Per-letter lookups on the 1990 Weltzeituhr** (city, zone or hour per letter). This needs the pre-1997 plate order, which no one has photographed completely.

**E. Larger algebraic systems**
10. Hill 6×6 and above (nitsuah has only a diagnostic), and Hill ≥ 3×3 with keyed columnar ≥ 10.
11. Fractionation + transposition on the 3×3×3 cube in the write direction (§30). It is undecidable from 24 cribs.

**F. Physical / on-site**
12. The measured bearing of the CIA compass rose and the lodestone deflection; the layout of the site ("buried out there somewhere", K5). This needs a person on site or a FOIA request.

**G. What would change everything**
13. **K5.** Sanborn says it uses the same system and shares coded words with K4; Paradigm plans to release it. Twice the ciphertext under one method would turn most of the "undecidable" items above (period 19, two free alphabets, long keys, layered systems) into decidable ones. Every script in this package can be rerun on K5 the day it is published.

## 32. Turning grille, two-segment split, keyboard layouts (`fleissner.py`, `split2.py`, `keyboard.py`)

**Turning grille (Fleissner) 10×10 + periodic key.** (jcolinpatrick's searches are also heuristic and inconclusive; new here are KA, the ciphertext-indexed key, both rotation directions and the shuffled-null calibration.)
- Setup: 25 holes, 4 rotations (clockwise or counter-clockwise), 97 letters plus 3 blank cells.
- Search: simulated annealing over the 4²⁵ grilles; key ST/KA × Vigenère/Beaufort, periods 1–12, indexed by plaintext or ciphertext position.
- Results:
  - planted control: best 20/24;
  - K4: best 20/24;
  - shuffled nulls: 19/24.
- **Inconclusive.** The annealer does not reach 24/24 even on the planted case, so it cannot separate K4 from noise. An exact (constraint-propagation) grille search is still needed.

**Two segments with independent systems.** (jcolinpatrick's E-SEGM-01 / `e_segmented_tape_01` uses two repeating keys without transposition; new here is a separate transposition per segment and arbitrary alphabets.)
- Split point s = 34…63: segment A holds EASTNORTHEAST, segment B holds BERLINCLOCK.
- Each segment gets its own identity or columnar transposition (widths 2–7, all orders, both directions) plus a periodic key: ST/KA ≤ 12 or any alphabet, plaintext- or ciphertext-indexed.
- Results:
  - The planted split (s = 50) is recovered.
  - K4: segment A consistent at 1 split (2 configurations), segment B at none, so **no split works for both**.
  - Shuffled nulls: A 0–1, B 0.

**Computer-keyboard layouts.**
- Alphabets: QWERTY, QWERTZ, AZERTY, Dvorak and the column-wise QWERTY/QWERTZ orders, each also reversed, run through every §27 family:
  - linear/clock/gear models, safe keys, clock-state masks;
  - Hill, autokey, running key over 11.66M letters;
  - Q(370k dictionary) → spirals / columnar.
- Keyboard geometry as a route: K4 written in staggered rows of 10/9/7, read by columns; tested with any periodic key ≤ 24, arbitrary alphabets. (Aligned rows are plain columnar widths 7/9/10, already covered.)
- Result: **0**, null 0.
- A key-neighbour "shift" (each letter replaced by the adjacent key) is a monoalphabetic substitution and is already excluded.

## 33. Nefertiti: the Egypt–Berlin link (`nefertiti.py`)

**Why.** Sanborn's November 2025 letter names a 1986 Egypt trip and the 1989 fall of the Wall. The Nefertiti bust (found at Amarna in 1912 by Borchardt) has been in Berlin since 1913. From 1967 to 2005 it stood in the eastern Stülerbau opposite Charlottenburg Palace in West Berlin, and in the early 1990s it got its own dimly lit ground-floor room (SMB).

**Geography.** Bearing from the 1990 site to the Weltzeituhr: 88.8° (due east, not ENE), 7.9 km. To the Mengenlehreuhr's 1990 site: 133.7°. No ENE alignment.

**Vocabulary (68 words).** NEFERTITI, NOFRETETE, NEFERNEFERUATEN, AKHENATEN/ECHNATON, ATEN, AMARNA, THUTMOSE, BORCHARDT, JAMESSIMON, CHARLOTTENBURG, STUELERBAU, TUTANKHAMUN, CARTER, CARNARVON, NILE, ANKH, etc. That gives 204 keyed alphabets (ST, reversed-ST and KA bases).

**Results.**
- §27 families with these alphabets (linear/clock/gear models, safe keys, clock-state masks, Hill, autokey, running key over 11.66M letters): **0**, null 0.
- The words as Q keys → 752 spirals → columnar 2–7 (all orders), ST/KA: **0**.
- The 370k-word dictionary as Q key with each Nefertiti alphabet (spirals, or columnar 2–7): **0**.
- Three layers (T+Q+Q, Q+Q+T, Q+Q+Q, T+T+Q; 363 words + Nefertiti words) with ST, KA, NEFERTITI- and AKHENATEN-keyed alphabets: **0**.

**Nefertiti words as extra cribs** (`nefcrib.py`). NEFERTITI, NOFRETETE, QUEENNEFERTITI, BUSTOFNEFERTITI, AKHENATEN, AMARNA, THUTMOSE, BORCHARDT and EGYPT, each placed at every free position (31–61 positions per word):
- 0 consistent with any ST/KA periodic key ≤ 48;
- 0 consistent with any-alphabet periodic keys ≤ 48 at the ≥ 3-coincidence level, the same as random strings of equal length.

Extra cribs cannot revive excluded families. In the still-open ones, a correct word and a random word behave identically, so the cribs alone cannot confirm or rule out Nefertiti in the plaintext.

## 34. Further structural ideas (`mirror.py`, `new_frame.py`, `new_panel31.py`, `new_mixsplit.py`, `new_autocol.py`, `new_doublecol.py`)

Each test has a planted positive control (recovered) and shuffled-K4 nulls (0 or at the K4 rate).
- **Reflection autokey** (key = plaintext/ciphertext at the mirror position c−i, every centre, 4 rules, ST/KA; crib positions 23–33 mirror onto 63–73 so 11 pairs are fully known): 1,320 models, **0**.
- **Frame stripping** (jcolinpatrick's `e_cable_format_hypotheses` tests OBKRU as an indicator with keyword ciphers). OBKR and/or KCAR treated as indicators, or a '?' cell prepended (lengths 93, 89, 98, 94), columnar 2–8 + periodic ≤ 12 (ST/KA/any alphabet): 462k permutation tests, **0**.
- **Carved 4×31 block read non-row-wise** (columns, serpentine, diagonals, spirals; widths 30–33; with or without K3 context cells; 882 + 72×28 permutations) + periodic ≤ 24: **0**.
- **Mixed families by segment** (one half K3-style substitution-free double rotation, the other a periodic key, any split 34–63): **0**, and provably impossible: K4 has E only at 44 and 92 and L only at 11, 15, 22, 53, while EASTNORTHEAST needs two E's and BERLINCLOCK two L's. No substitution-free transposition of either half can exist (the whole-text version is in jcolinpatrick's `docs/elimination_tiers.md`).
- **Autokey under columnar 2–7** (plaintext autokey, reflected keys, and plaintext lags 26–52; the ciphertext-autokey part is already excluded by nitsuah): 11,824 perms × all configurations, **0**; the null floor is 4 checks, a real instance gives 10–24.
- **Double columnar transposition** (jcolinpatrick's E-FRAC-46 covers only the 6×6 pair, scored; here widths 2–6 × 2–6, all orders, forward composites and inverses; 1,520,252 distinct permutations) + periodic ≤ 12: **0** ST/KA; any-alphabet hits at the null rate.

## 35. Period 19 with arbitrary (non-keyword) alphabets, English-scored; letters as numbers (`p19sa.py`, `numcode.py`)

**Period-19 Quagmire IV with arbitrary alphabets.** (Quagmire III at period 19 is already impossible for any alphabet by keithadler's exact solver.) Simulated annealing over the unknown alphabet(s), with the 19-letter key forced by the cribs. Score = English quadgram log-probability (Gutenberg corpus) − 6 × crib conflicts.
- Planted controls (English text, random alphabets) are recovered almost letter-perfect: Quagmire IV at 0 conflicts and score −456, 1 of 12 runs at 150k iterations.
- **Quagmire IV**, 12 runs × 150k iterations:
  - K4 best: −467/−480 with 4 conflicts;
  - null: −477 with 5 conflicts.
- **K4 never reaches a crib-consistent English decryption and is indistinguishable from shuffled text.**
- Caveat: at the measured control success rate, a Quagmire IV solution could be missed with roughly 35% probability in 12 runs.

**Letters as numbers.** A=0/1 or KRYPTOS order, concatenated into a digit stream, then searched for dates and coordinates (1990, 1989, 1986, 1912, 3857, 7708, 1109, …), runs of equal differences, and row sums of the carved 31-wide rows.
- K4 shows 5 target substrings against 3.3 expected by chance.
- "1912" (the Nefertiti find year) appears in 53% of shuffled texts.
- "1109" appears in 8% of shuffled texts.
- The ST 1-based letter sum 1326 = 26 × 51 (a 1-in-26 coincidence that is invariant under shuffling).
- Nothing exceeds chance.

## 36. Exhaustive follow-ups on the historically motivated gaps (`new_lenmask.py`, `new_fleissner_exact.py` + `fl_exact7.c`, `new_tradecraft.py`, `new_subtext_4layer.py`, `cp_test.py`, `cp_qii.py`, `qii_more.py`)

Each test has a planted positive control (recovered) and shuffled-K4 nulls.

**Length-changing masking** (letters added or removed in the gap, so the BERLINCLOCK block's key index is shifted by Δ = ±1…±10 relative to EASTNORTHEAST):
- Direct periodic key, periods 1–26, ST/KA and arbitrary alphabets: **0** (null 0.03–0.08 per shuffle).
- Plaintext-indexed periodic key ≤ 12 under 164,664 transpositions (columnar 2–8 all orders, K3 double rotations): **0**.
- Running key over 11.66M letters × 270 alphabets: **0**; neither crib key fragment occurs anywhere in the corpus, so this holds for any Δ.
- 2,360 linear clock/gear models: **0**.

**Turning grille, exact.** A depth-first search over the 4²⁵ grilles that prunes as soon as no key model survives (ST/KA × Vigenère/Beaufort, periods 1–12). Plaintext-indexed key; 8 variants (rotation direction, row/column read, blanks at the start or end). The planted grille is found in all 8 variants. On K4, the clockwise / row-read / blanks-at-end variant has covered about 42% of the grille space with **0 survivors** (null: 17% covered, 0 survivors). The other seven variants and the rest of this one would need many more CPU hours. Arbitrary alphabets per class are not decisive for grilles: even random text matches 10⁵–10¹⁴ grilles.

**Tradecraft transpositions** with periodic key ST/KA ≤ 24 and arbitrary alphabets ≤ 12, both directions (18.6M permutations):
- Myszkowski, every key pattern for widths 3–9;
- VIC-style disrupted columnar, widths 5–9 with every key order and width 10 sampled at 10%, in two triangle variants;
- irregular columnar, widths 2–8.

Result: **0** keyed hits; 5 weak arbitrary-alphabet hits at the minimum threshold, at the null rate.

**Only part of K4 enciphered** (the method changes part-way). 11 substrings bounded by the W's or the ends that contain both cribs; identity, columnar 2–8 and K3 rotations, periodic key ≤ 12 (ST/KA/any alphabet): **0** keyed hits (one weak any-alphabet hit, null two). A planted cipher on positions 21–96 after random junk is recovered.

**Four-layer PK-style stacks** (T+Q+Q+Q, Q+Q+Q+T, T+T+Q+Q, Q+T+Q+T; 89 thematic words incl. all PK keys; KA and ST; exact, meet-in-the-middle): **0**. All 8 planted stacks are recovered.

**Antipodes and the Cyrillic Projector.**
- Antipodes' K4 matches Kryptos K4 letter for letter (in both copies on Antipodes, per jcolinpatrick's transcription).
- The Cyrillic Projector cipher is a Quagmire II: plaintext in the normal Russian alphabet, ciphertext in the ТЕНЬ-keyed alphabet, repeating key МЕДУЗА. This was re-derived and verified on 72 letters.
- The same form on K4 (normal/KA alphabets in all four pairings): no period 1–26 is consistent directly. With a periodic key ≤ 12 it also gives **0** under columnar 2–8, all 72,200 K3 rotations, 191,272 dial permutations and 752 spirals. Fixed keys (MEDUSA, SHADOW, TEN, KRYPTOS, PALIMPSEST, ABSCISSA, BERLIN, CLOCK, ANTIPODES) reach < 4/24.
- **Correction to earlier sections.** The arbitrary-alphabet check used in this package (`gen.run`) needs repeated plaintext letters in a key class. It did not by itself cover Quagmire I/II with *fixed* different alphabets for plaintext and ciphertext (e.g. normal/KA). That form has now been run separately over the transposition families listed above (`cp_qii.py`, `qii_more.py`); other families in this package used ST/ST and KA/KA keys only.

## 37. Flower, lotus and paper fortune teller (opening/closing)

**Idea.** The text is laid out like a flower and read as it opens and closes: phyllotaxis spirals (sunflower/lotus seed heads), Fibonacci parastichies, rings of petals, or the folded paper "fortune teller" (cootie catcher) that opens and closes in two directions.

*Prior work:* jcolinpatrick's `e_clock_reading_order.py` tried a 137.508° angular step around three centres on the sculpture grid. That is a geometric sweep of the physical layout. The Vogel-index orders, parastichies, petal rings and paper-fold orders below are not in the community repos we checked. Folding was mentioned at PhreakNIC 26 ("You fold the paper into origami…"), but by the speaker, not by Sanborn.

**Phyllotaxis (`flower.py`)**
- **Vogel spiral.** Seed *n* sits at angle *n*·θ and radius √*n*, with θ = golden 137.508°, 222.492°, Lucas 99.502° and noble 151.137°.
- **Reading orders:**
  - by angle (both directions, and outward/inward);
  - parastichy families 3/5/8/13/21/34/55/89, straight or alternating ("open/close");
  - petal rings 4/6/8/12/24, outward or alternating;
  - all decimations *k* = 2–96 from every start.
- 9,343 unique permutations, each used as encryption and as decryption.
- **Tests:** a periodic key ≤ 24 (ST, KA and arbitrary alphabet per class), and Quagmire I/II with normal/KA pairs ≤ 24.
- **Results:**
  - Planted control is recovered (golden parastichy-13 + KA period 6).
  - **K4: 0. Shuffled-K4 null: 0.**

**Paper fortune teller (`cootie.py`)**
- **Fold model.** A 10×10 sheet holds 97 letters plus 3 pads (at the end, the start, or 3 corners). The letters are written in rows, columns, boustrophedon or a spiral. The sheet is then:
  - blintzed (four corners to the centre),
  - optionally flipped,
  - optionally blintzed again,
  - optionally folded in half.
- Each cell's final position and layer is tracked by geometric reflection.
- **Reading orders:**
  - position-then-layer (up/down, rows/cols);
  - layer-then-position;
  - by angle around the centre with layer.
- 1,904 unique permutations.
- **Results:**
  - Control is recovered (KA period 7 through a fold order).
  - **K4: 0** for a periodic key ≤ 24 (ST/KA/arbitrary) and for Quagmire I/II normal/KA. Null: 0.

**Flower and fortune-teller words as keyed alphabets**
- **Words (52):** LOTUS, BLUELOTUS, NYMPHAEA(CAERULEA), SESHEN, NELUMBO, PAPYRUS, SUNFLOWER, HELIANTHUS, HOROLOGIUMFLORAE, LINNAEUS, FLOWERCLOCK, BLUMENUHR, MORNINGGLORY, MIRABILIS, PHYLLOTAXIS, FIBONACCI, ORIGAMI, COOTIECATCHER, FORTUNETELLER, CHATTERBOX, LOULOUDI and others (full list in `cootie.py`).
- **Alphabets.** Each word gives a keyword-mixed alphabet, used as Quagmire III (W/W), I (W/ST), II (ST/W), W/KA and KA/W. That makes 260 alphabet pairs.
- **Tests:** Vigenère and Beaufort sign, periods 1–26, key indexed by pt or ct position, under the identity, all fold permutations and all flower permutations (11,248).
- **Results:**
  - Control (LOTUS QIII, period 5) is recovered.
  - **K4: 0. Null: 0.**
- A periodic Vigenère key that is itself a flower word needs no separate run: it is already covered by the arbitrary periodic key ≤ 24 above.

**Conclusion.** Under the forms tested, neither a flower-like reading (spiral, petals, opening/closing) nor a paper-fortune-teller fold, nor flower vocabulary as alphabets, produces K4 with a periodic or Quagmire key.

## 38. Hand-made masks from features visible at home; raised YAR as the mask's top-left corner (`handmask.py`, `handmask2.py`, `nullshift.py`)

**Idea.** Sanborn says K4 can be solved without visiting the site. So a hand-cut mask would have to be built from features visible in public photos and transcriptions:
- the raised letters YAR,
- the four `?`,
- the misspellings,
- the extra E's in the entrance Morse code,
- the tableau,
- the K1–K3 texts.

Thanos Zekios's specific proposal is that the raised YAR marks where the mask's top-left corner goes.

**Layout used.** We take the line breaks of the cipher side as published (Wikipedia transcription, 28 lines). With that layout:
- YAR sits on line 15, columns 4, 5 and 7 (1-based). Line 15 is the first line of the right-hand panel (K3+K4), so YAR is near that panel's top-left corner.
- K4 occupies lines 25–28 of that panel.
- So a mask anchored at YAR must be at least 11 rows tall to reach K4. The entrance Morse lines (7 rows, or 12 if split by plate) barely reach K4.
- K3's plaintext reads "…A TINY BREACH IN THE UPPER LEFT HAND CORNER…". We record the coincidence; it is not evidence.

**Mask anchors (13).** The mask's top-left corner was placed at each of:
- Y, A, R;
- the panel corner and the sheet corner, for comparison;
- each of the four `?`;
- the cipher cells of the three misspellings: IQLUSION Q, UNDERGRUUND U, and DESPARATLY A (mapped through K3's double rotation, verified 336/336).

**Hole patterns** (all also mirrored):
- the Morse filler E's as holes, or their complement (lines whole or split by plate);
- the left panel's `?` and misspelling cells;
- the YAR shape;
- regular grids h×w (1–8);
- diagonal lattices;
- tableau coincidences: holes where the tableau letter, laid at the anchor straight or mirrored, equals the K4 letter;
- the left panel (K1–K2 ciphertext or plaintext) laid at Y/A/R with holes where letters coincide (`handmask2.py`);
- the K1–K3 plaintext laid over K4 by index or on the panel, with holes at letter sets (KRYPTOS, YAR, QUA, LOE, X, E and each single letter).

**Josephus skip readings.** Step sequences tested:
- Morse E counts 2,1,6,3,2,2,5,3,1,1 (these sum to 26), also without the last E;
- the K2 coordinate digits;
- YAR (24,0,17 or 25,1,18);
- the misspelling letters (QUA / LOE);
- KRYPTOS.

Each runs from every start, with step or step+1.

**Uses tested:**
- **Holes first, then the rest (or the reverse), as a reading order**, plus the inverse. That gives 7,494 unique permutations, tested with a periodic key ≤ 24 (ST, KA and arbitrary alphabet) and with Quagmire I/II normal/KA ≤ 24.
  - Control: a Y-anchored 3×5 grid plus KA period 7 is recovered.
  - **K4: 0. Null: 0.**
- **Tableau laid at every offset (straight or mirrored) as a running key** (Vigenère, Beaufort and variant Beaufort, ST/KA): best **5/24**, the same as the null (5/24).
- **Holes as nulls.** For a periodic ST/KA key on the reduced text, only two numbers matter for the cribs: the nulls before position 21 and the nulls between the two cribs. So *every* null mask, of any size and in any position, can be tested exhaustively (`nullshift.py`):
  - **ST/KA Vigenère and Beaufort, periods 1–26: 0** (control recovered, null 0).
  - With an arbitrary alphabet per class the test is undecidable from 24 cribs: K4 and the shuffled null give the same count.
  - *Prior work:* nitsuah covers null masks with ST/KA keys; jcolinpatrick's `e_k0_morse_null_mask.py` derives Morse null masks; ours subsumes both for periodic ST/KA keys.
- *Prior work* on overlays as keys: our `board.py` (K1–K3 panels as key at offsets) and jcolinpatrick's E-GRILLE-01 / installation sweep (fold-as-Cardan).

**Conclusion.** These forms all give 0: a mask from the public anomalies, anchored at YAR or elsewhere, used as a reading order, a null mask or a key overlay with periodic or Quagmire keys. Still open:
- an irregular mask whose holes are not one of the patterns listed;
- the combination of such a mask with a non-periodic key.

## 39. Reading along light rays of a compass bearing ("look through the hole toward ENE") (`raylight.py`)

**Idea.** The raised YAR is read as "look through here", and the light ray follows a compass bearing: perhaps EAST-NORTH-EAST (22.5° above east), perhaps another.

**Parallel light is already covered.** A parallel beam through a flat sheet onto a parallel surface only shifts the image by a constant offset. The constant-offset overlays (tableau and K1–K3 panels as key at every offset, §38 and `board.py`) already cover this case: best 5/24, the same as the null.

**What this section tests: rays as a reading order.**
- **Layouts:** K4 laid out in its carved 4-line layout, in rows of width 4–48, and with every starting offset (the paper's corner not at K4's first letter).
- **Directions (16):** compass bearings E, NE, N, NW, SE and the half-points, with ENE also approximated as 1:2, 2:5, 5:12 and 1:3.
- **Reading:** ray by ray, in both orders across the rays and both senses along each ray, plus inverses.
- **Total:** 73,836 unique permutations, tested with a periodic key ≤ 24 (ST/KA/arbitrary alphabet) and Quagmire I/II normal/KA.
- **Results:**
  - Control (carved layout, ENE 1:2, KA period 8) is recovered.
  - **K4: 0. Null: 0.**

**Not testable from home:** projection through the curved S-shaped copper onto another part of the screen. That needs the exact 3D shape and orientation of the sculpture, which are not published to the required precision.

## 40. Folding the sheet: K1–K3 (solved or unsolved) and the tableau projected onto K4 (`foldproj.py`); sunlight simulator (`kryptos_sun.html`)

**Idea (Thanos Zekios).** The masking is done by folding the sheet, as in the paper fortune teller, or simply in two or in four. Whatever lands on K4 is the mask or key. What lands there could be K1–K3 already solved (the plaintext), or K1–K3 unsolved (the ciphertext), with or without the misspellings and with or without the raised YAR.

**Sources (11 variants, plus the tableau):**
- **Solved panel** (K1–K3 plaintext laid on their carved cells):
  - as carved;
  - misspelled letters corrected in place (ILLUSION, UNDERGROUND, DESPERATLY);
  - fully corrected (DESPERATELY, which shifts the rest of K3 by one cell).
- **Unsolved panel:**
  - as carved;
  - with the three misspelled letters correctly re-encrypted (K1 with KA + PALIMPSEST, K2 with ABSCISSA, K3 letter in place).
- **Raised YAR:** each of the five above is used with the YAR letters kept, or blanked.
- **The tableau,** straight or mirrored, with row shifts −4…+4.

**Folds:**
- in half horizontally (row *r* onto 27−*r*);
- in half vertically (sheet width 33 or each line's own length);
- in quarters;
- the K1–K2 panel stacked on the K3–K4 panel, straight or flipped, with column shifts −8…+8;
- the top-left quarter laid on K4's rows.

That gives 427 combinations.

**Uses tested:**
- **(a) Landing letter as the key** (Vigenère, Beaufort, variant Beaufort; ST/KA).
  - **K4 best 6/22**; the shuffled null reaches 5/20, so this is chance.
  - A planted control scores 24/24.
- **(b) Landing letter + a periodic key** (periods 1–26, consistency required at all available cribs).
  - **K4: 0. Null: 0.** A planted period-5 control is recovered.
- **(c) Landing letter as a mask** (holes where it equals the K4 letter, or where nothing lands) → reading orders. That is 521 permutations, tested with a periodic key ≤ 24 (ST/KA/arbitrary alphabet) and Quagmire I/II.
  - **K4: 0. Null: 0.**
- As null masks these are covered by `nullshift.py` (§38): 0 for any null set with a periodic ST/KA key.

**Sunlight simulator.** `kryptos_sun.html` is a stand-alone page:
- **What it shows:** the cipher side (and optionally the tableau) cut out of A4 paper, flat or bent into a C or an S shape, lit by the sun at any azimuth and elevation.
- **Where the light lands:** a wall, the floor, or the sheet itself.
- **On the sheet itself:** it reports which letters the light passing through each K4 hole lands on, and scores them as a key against the cribs.
- **Sweep:** it can sweep every sun position. With the sculpture-like S bent to 67.5°, the best is 3/16, i.e. chance.
- **Rotation is covered:** self-projection depends only on the sun's angle relative to the sheet, so the azimuth sweep also covers every orientation of the sheet.
