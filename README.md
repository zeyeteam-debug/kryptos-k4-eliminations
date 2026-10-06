# Kryptos K4: crib-constrained eliminations not (fully) covered elsewhere

Author: Thanos Zekios, October 2026. The computations were carried out with the help of Claude (Anthropic).

This note documents a set of exact searches on K4 that, as of early October 2026, are **absent or only partly covered** in the three public elimination efforts I compared against:
- jcolinpatrick/kryptos (KryptosBot)
- nitsuah/kryptos
- keithadler/kryptos

I left out everything those repositories already cover. Where they cover a family partly, I describe only the extension. All scripts are in `scripts/`.

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
| W's are inserted nulls (removed, 92 letters): periodic key 1–24 (any alphabet), and columnar widths 2–8 plus a key of length ≤12. | Excluded (one threshold-level "any-alphabet" hit, at the same rate as the null). |
| Every other letter treated as a null class (direct test). | Excluded. |
| Segment lengths used as row lengths of a transposition (read by columns) plus a periodic key. | Excluded. |
| 23 two-letter words after `EASTNORTHEAST` × ~560 four-letter fillers before `BERLINCLOCK`, added as extra cribs (W-restart hypothesis). | No second coincidence. |

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
| "YAR-shaped" holes (relative columns 0,1,3), tiled every q = 4–32 columns; plus running keys starting exactly at the raised Y/A/R, and the tableau read from row N (the extra-L row) | exact | best 7/24 (null 8) |
| Sinusoidal (S-curve) path: amplitude 1–6 rows, wavelength 8–97, all phases, row drift, start rows | exact | best 9/24 (null 9) |

Regular masks over the tableau reduce to odometer/Alberti keys (section 6), and those are excluded too. An irregular, hand-made mask is not testable.

## 3. Hill: any matrix, plus transposition (extension)

`hillx.py`, `hillT.py`, `hillchk.py`. Exact linear algebra mod 2 and mod 13 over **all** matrices, all block alignments, ST and KA indexing.
- **Direct:** 2×2, 3×3 and 4×4 are inconsistent.
- **Affine** (c = Mp + b) and **"periodic key then Hill"** (block-class offsets, ≤ 6 classes): inconsistent for 2×2 and 3×3. 4×4 affine and 5×5 and above are untestable (0 redundancy).
- **Hill (2×2, 3×3) with columnar transposition at widths 2–9, both orders** (transpose→Hill and Hill→transpose):
  - widths 2–8: 0 hits;
  - width 9: 684 hits, all with only 3 equations for 2 unknowns. Most solutions are **non-invertible** matrices. All 276 invertible ones decrypt to gibberish (quadgram ≈ −7.5/char).
  - 3×3 transpose-first was run at widths ≤ 8.

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

- **Direct:** inconsistent for every period 1–24 except **19**, which has a single strong check and is untestable. This matches keithadler's finding that Quagmire IV at length 19 stays satisfiable because of the R→P repeat 38 = 2×19 apart.
- **Per crib segment alone:** `EASTNORTHEAST` excludes periods 1–5 and 9; `BERLINCLOCK` excludes 1–4 and 7.
- **With ~1.1M transpositions, periods 1–12, key indexed by plaintext or ciphertext position:** 0 hits (null 0). The transpositions are:
  - columnar widths 2–9, all orders, both directions;
  - "24-face cylinder" readings (24/12 columns, straight or helical, any start/step/direction);
  - opposite-direction rotations (two halves, even/odd, grid rows/columns as gears).
- **K3-style double rotation,** all widths 2–96 × 2–96, CW/CCW, both directions (72,200 permutations): 1 threshold hit (null also 1). The implementation reproduces K3 exactly (336/336).

## 6. Mechanisms not found elsewhere

| Mechanism (script) | Space | Result |
|---|---|---|
| Single stepping rotor (Hebern type), **unknown wiring**, step every N = 1–48, both directions, ST/KA (`rotor.py`) | exact injectivity check | 0 (null 0/100) |
| Odometer key k = a[i mod p] + b[⌊i/p⌋ mod q] (`odometer.py`) | p + q ≤ 20 | 0 |
| Alberti-disk stepping k = a[i mod p] + d·⌊i/N⌋ (`odometer.py`) | p ≤ 18, N = 2–48 | 0 |
| 6×6 Bifid (`frac.py`) | periods 2–97 (all offsets for ≤ 24) | infeasible everywhere |
| Trifid (`frac.py`, `trisa.py`) | periods 2–97; the 4 feasible cases (13/5, 16/0, 16/7, 16/8) checked by SA | SA cannot fit the cribs (15–17/24); control fully recovered |

## 7. Error tolerance under transposition (extension)

`err1.py`. The models are: one crib letter wrong (drop each anchor in turn), one extra letter (97 deletions) and one missing letter (98 insertions). Tested:
- periodic key 1–24 with ST/KA and any alphabet;
- **plus columnar widths 2–8 with periods ≤ 12, drop-one-anchor model.**

Result: 0 hits. Controls recovered a planted typo and a planted extra letter.

## 8. Numeric "masking" keys (extension)

`numkeys.py`, `numsrc.py`, `coordcol.py`. Sequences used as a non-repeating shift, at every start offset, ×1/−1/2/3, add or subtract, ST/KA:
- digits of π, e, √2, √3 and φ (1-digit, 2-digit blocks, sliding 2-digit);
- squares, cubes, triangular numbers, primes and prime gaps, powers of 2/3/5/7, Lucas;
- K2 coordinates and dates (9 Nov 1989; 3 Nov 1990);
- Morse code of the sculpture text (code lengths, dots, dashes, binary values, 5-bit blocks);
- Weltzeituhr: UTC offset per city or per letter, and name lengths;
- W segment lengths and positions.

About 120,000 configurations gave best 6–7/24, at the null level. Coordinate and date digits were also used as columnar orders (widths 3–13, ascending, descending or plain) with any periodic key ≤ 12: 0.

## 9. Small explicit checks (for completeness)

The following were explicitly checked, as Quagmire keys and keyed alphabets in all combinations, globally and restarting per W segment. All are at chance level (≤ 6/24).
- Keyword sets around WW / William Webster: `WILLIAMWEBSTER`, `WILLIAMHWEBSTER` (15 letters, equal to the two crib segments), `WILLIAMHEDGCOCKWEBSTER` (22, equal to the last segment), `JUDGEWEBSTER`.
- `ORDINATE`, `COORDINATE(S)`, `APPLICATE`, `YAXIS`, `ABSCISSAORDINATE` and others. Any key of length 13–16 with columnar widths 2–8 also gives 0.
- Misspelling letters (`QUA`, `LOE`, …) and `YAR`/`RAY`.

## Caveats

- Everything above assumes the published cribs are exact and position-aligned. Section 7 relaxes this by one character.
- The "any-alphabet" tests cannot reach periods above 24, and they are weak at 19.
- Three-layer results only exclude the listed words.

## Running

The scripts need Python 3 and numpy. `pk_controls.json` holds the PK4 and PK6 plaintexts, used only as controls. `quad.pkl` contains English quadgram log-probabilities. Each script prints its control result first, then the K4 result.

## 10. Period-19 Quagmire IV with free alphabets: a capability limit

The any-alphabet test (section 5) leaves period 19 open, as keithadler also reports. The reason is the R→P repeat 38 = 2×19 apart.

I attacked it directly with simulated annealing (`q4p19.py`):
- unknowns: the plaintext alphabet, the ciphertext alphabet and the 19 shifts;
- crib consistency imposed as a penalty;
- score: English quadgrams.

**Positive control fails.** A planted Quagmire IV ciphertext with period 19 was built from 97 letters of English with the same cribs and random alphabets. Ten runs of up to 3M iterations each did **not** recover it:
- accuracy was about 33%, essentially the cribs only;
- the best score was −4.81/char, against −4.58 for the true plaintext.

So a negative K4 run with this method would mean nothing. Period-19 Quagmire IV with free alphabets cannot be decided from the 97 letters with this approach. Keyword-generated alphabets are already covered by nitsuah's dictionary × dictionary Quagmire IV search.

## 11. Weltzeituhr (Alexanderplatz) face transcription (2017 state) and the 1990 problem

The per-face, per-line city order is listed as a blocker in nitsuah's open list ("per-letter Weltzeituhr lookup keys"). Below is my transcription from user photos (Google Maps contributions, Feb/Jul 2017, plus 2026 screenshots).

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

Checks already run with the 2017 names as keyed alphabets or keys (Quagmire I–III, periods 1–24), including WILNA, NIKOSIA, the full +2 face and the old names LENINGRAD, ALMAATA and BRATISLAVA: no hits. Script: `cities.py`, `cities_eu.py`.

## 12. Sub → transposition → sub with two independent periodic keys; width-21 routes

**Model** (`stackT.py`):
- c[m(i)] − e·p[i] = K1[i mod a] + K2[m(i) mod b];
- the first key is indexed by plaintext position and the second by ciphertext position, after transposition m;
- e = ±1 covers all Vigenère/Beaufort combinations;
- the check is exact, by weighted union-find mod 26.

This is a fully exhaustive version, within these ranges, of the "substitution + transposition + substitution" class that KryptosBot's E-FRAC-52 sampled at 1.32%.

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

KryptosBot's running-key-plus-columnar search relied on a Bean filter that is invalid across a transposition, so it is not exhaustive. Here the check is exact (`rk_trans.py`):
- for every key offset in the text and every transposition, all 24 crib letters must match;
- the key is indexed by ciphertext position (substitution after transposition) or by plaintext position (substitution before);
- rules: Vigenère, Beaufort or Variant; alphabets KA and ST.

**Howard Carter, *The Tomb of Tut-ankh-Amen* Vol. I** (the source of K3; archive.org full text, 287,513 letters, forward and reversed; `carter_run.py`):

| Transposition family | Size | KA | ST |
|---|---|---|---|
| none | 1 | 0 | 0 |
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
