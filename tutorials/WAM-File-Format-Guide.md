# Reading & Writing OpenWAM `.WAM` Files — a topology-first guide

> Audience: yourself (Claude), other AI agents, and humans who need to **reconstruct an engine's
> architecture from a `.WAM`**, or **hand-write / edit** one. Companion to `AllFields.md` (which covers
> the per-field UI meanings). This doc is about **connectivity**: how the blocks reference each other to
> encode a gas-dynamic network. Worked example throughout: `Test5(4S-G-T)/Mitsubishi 3B21T engine.WAM`.

---

## 0. First principles

- A `.WAM` is **XML-tagged but read with `fscanf`**. The `<BLOQUE_… >`, `<TUBO id='n'>`, `<NODO id='n'>`
  tags are **landmarks for humans**; the parser consumes a **flat stream of whitespace-separated numbers**
  in a fixed order. **Position matters, tags don't.** When hand-editing, never change the *count* of
  numbers on a line without knowing what each one is — everything after shifts.
- **Units are SI**: metres, seconds, degrees, kg, K internally — but the file stores **bar** for pressure,
  **°C** for most temperatures, **metres** for lengths/diameters (the UI shows mm — see
  `[[openwam-file-units-meters-not-mm]]`). Angles in **degrees**.
- IDs are **1-based** in cross-references and resolved as `Array[ID-1]`.

---

## 1. The connectivity model (the one thing to internalise)

There are **two independent ID spaces**, and the whole topology is the mapping between them:

| Space | Block | Role in the graph |
|---|---|---|
| **Ducts / `TUBO`** | BLOQUE_IV | **Edges** — 1-D gas-dynamic pipes |
| **Nodes / `NODO`** | BLOQUE_VIII | **Vertices** — every place a duct *end* attaches |

**Rule:** every `TUBO`'s first two integers are the **node id at its left end and its right end**. Every
`NODO` declares **what kind of vertex it is** and, if it touches a 0-D element, **which element**.

**0-D elements are NOT nodes.** Cylinders (BLOQUE_III), deposits (BLOQUE_VI), the compressor (BLOQUE_VII)
are separate ID spaces. A duct never lists a cylinder or deposit directly — it lists a **node**, and that
node carries the reference to the cylinder/deposit/compressor. This is the WAM encoding of the UI rule
"a duct must reach a 0-D element through a connection element" (valve / throttle / node-to-deposit /
stator-rotor).

So the graph is: **duct → node → {another duct | a 0-D element via a connector | a boundary condition}.**

### Reconstruction algorithm (read a WAM → topology)
1. **Edges:** for each `<TUBO id='t'>`, record `t: (nodeL, nodeR)` from the first two integers. Grab the
   duct-type digit (line 2) and the `Ø_in length Ø_out` triple for labelling.
2. **Vertices:** for each `<NODO id='n'>`, read the **type code** (first integer) and its payload
   (§4). Classify as boundary / tube-union / ramification / valve-to-cylinder / union-to-deposit /
   compressor.
3. **Resolve refs:** for valve/union nodes, resolve the payload to the cylinder (BLOQUE_III), the
   deposit (BLOQUE_VI) + the connecting valve (BLOQUE_V), or the compressor (BLOQUE_VII).
4. **Walk it:** start at a boundary node (atmosphere / reservoir) and follow duct↔node links. Deposits
   and cylinders are "hubs" reached through their connector nodes.

---

## 2. Block map (order is fixed)

| Block | Contents | Topology relevance |
|---|---|---|
| BLOQUE_0 | version | none |
| BLOQUE_I | execution (dt out, ambient P/T, mode) | none |
| BLOQUE_II | atmospheric composition | none |
| BLOQUE_III | **engine + cylinders + combustion (Wiebe) + firing order** | cylinders referenced by valve nodes |
| BLOQUE_IV | **`TUBO` ducts** | **edges** |
| BLOQUE_V | **`VALVULA` valves & unions** (connectors) | referenced by valve/union nodes |
| BLOQUE_VI | **`DEPOSITO` 0-D volumes** (plenums, turbine, chambers) | referenced by union-to-deposit nodes |
| BLOQUE_VII | **`COMPRESOR`** (+ map) | referenced by the compressor node |
| BLOQUE_VIII | **`NODO` nodes** | **vertices** |
| BLOQUE_IX | **`TURBOGRUPO`** (shaft) | couples compressor ↔ turbine deposit |
| BLOQUE_S | **`SENSOR`s** | control feedback / inputs |
| BLOQUE_C | **`CONTROLADOR`s + `TABLA1D`** | control logic; setpoint chaining |
| BLOQUE_X / XI / XII | mean / instantaneous / spatial results requested | none (output selection) |
| BLOQUE_XIII | **compressor ↔ suction-deposit link** | compressor suction chamber |
| BLOQUE_XIV / XV | external-calc vars / free text | none |
| COORDENADAS | UI canvas x/y for every element | none (layout only) |

---

## 3. Ducts — BLOQUE_IV `TUBO`

```
<TUBO id='1'>
17 18 1 12 1.000000e-02 4.000000e+01 8.000000e+01 1.100000e+00 5.000000e+00      <- line A
 1 1.000000e+00 1.000000e+00                                                     <- line B: duct-type, HT-adj, friction-adj
0.000000e+00 1.000000e+00                                                        <- initial state (velocity, ...)
5.000000e-03 1                                                                   <- mesh size (m), ...
2 8.000000e-01                                                                   <- solver/limiter params
4.000000e-02 3.000000e-01 4.000000e-02                                           <- Ø_in(m)  length(m)  Ø_out(m)
1.000000e+00 8.000000e-01 0                                                      <- (variable-wall extras / IC)
1                                                                                <- n wall layers
1 0 2.700000e+03 9.000000e+02 2.000000e+02 2.000000e+00 0 0                      <- layer: main/fluid, ρ, cp, k, thickness(mm), emis...
</TUBO>
```
- **Line A, fields 1–2 = the node endpoints** (`17 18`). *This is the topology.* Field 4 = **N of
  intercooler ducts** (`12` here → this is the intercooler; `1` for a normal duct). Remaining line-A
  fields are mesh/initial-temperature/Courant-type numbers (see `AllFields.md`).
- **Line B, field 1 = duct type** (`TTubo.cpp` `switch(TipTC)`): **1 = intake pipe (`tubo`), 2 = exhaust
  pipe, 3 = exhaust port (`pipa`), 4 = intake port; 0 = plain pipe (no engine).** Remember
  `[[duct-type-pipa-port-tubo-pipe]]`: **`pipa` = PORT, `tubo` = PIPE**; ports pin the wall to engine
  coolant and hide the Heat-Transfer tab.
- The `Ø_in length Ø_out` triple gives geometry (all in **metres**). Equal in/out Ø = straight pipe.
- Wall-layer block appears only for **variable wall-temperature** ducts (e.g., the intercooler).

Only lines A(1–2) and B(1) matter for topology; the rest is geometry/physics — defer to `AllFields.md`.

---

## 4. Nodes — BLOQUE_VIII `NODO` (the vertex dictionary)

The block opens with the node count and a short header of category counts — **don't rely on the header;
trust the per-`<NODO>` tags.** Each node's **first integer is a type code**:

| Code | `type=` | Meaning | Payload after the code | Resolves to |
|---:|---|---|---|---|
| 0 | EXTREMO ATMOSFERA | atmosphere outlet BC | `1.0` (ref) | boundary |
| 1 | EXTREMO REMANSO | reservoir/stagnation inlet BC | `P_bar T_C … ` + composition | boundary (a fixed-state inlet) |
| 6 | UNION ENTRE TUBOS | plain 2-duct join | `0.0 0.0` | connects the 2 ducts whose ends = this node |
| 7 | VALVULA DE ADMISION | intake valve | `self  cyl  valveDef` | **cylinder** `cyl` via **valve** `valveDef` |
| 8 | VALVULA DE ESCAPE | exhaust valve | `self  cyl  valveDef` | **cylinder** `cyl` via **valve** `valveDef` |
| 9 | PERDIDA DE PRESION LAMINAR | linear pressure-loss connector | `FK` | joins the 2 ducts on this node; loss ∝ `FK·\|U\|` (`FK>0`) |
| 10 | PERDIDA DE PRESION TURBULENTA | quadratic pressure-loss connector | `FK` | joins the 2 ducts on this node; loss ∝ `FK·U²` (`FK>0`) |
| 11 | UNION A DEPOSITO | duct ↔ 0-D volume | `self  deposit  valveDef` | **deposit** `deposit` via **valve/connector** `valveDef` |
| 12 | RAMIFICACION | multi-duct branch | *(none)* | all ducts whose ends = this node |
| 13 | COMPRESOR DE TORNILLO | screw compressor / supercharger | `cvId  ctrl  spd/ratio  Pin Tin  C1 C2 C3  D1 D2 D3  E1 E2 E3 E4 E5 E6` + composition | boundary (prescribed inflow into the one attached duct); see §4a |
| 17 | COMPRESOR | compressor connection | `compId` | **compressor** `compId` (discharge duct is the tube on this node) |

Notes:
- `self` (2nd number on valve/union nodes) equals the node's own id here — a redundant
  boundary-condition instance index; leave it equal to the node id when writing.
- **Which ducts touch a node?** Not stored in the node — found by scanning BLOQUE_IV for tubes whose
  endpoint = that node id. A RAMIFICACION with three tubes pointing at it is a 3-way junction; the WAM
  models manifolds as chains of these (no plenum needed).
- A `type=17` compressor node names the compressor; the compressor's **suction** side is *not* here —
  it's a deposit named in **BLOQUE_XIII** (`compId  depositId`). Its **discharge** is simply the `TUBO`
  whose end is this node.

### 4a. Screw compressor node (type 13) payload — exact layout

`TCCCompresorVolumetrico`. A `<NODO type='COMPRESOR DE TORNILLO'>` reads, **in this order** (whitespace
/ newlines are interchangeable — the real sample splits it across lines):

```
13            node type code
cvId          screw-compressor index. MUST be present but the value is now IGNORED for linking
              (populated sequentially like every other BC array). The UI writes it 0-based; the engine
              used to require 1-based → segfault (fixed 2026-08-22, TOpenWAM.cpp:1523).
ctrl          speed-control mode: 0 = Proper (own fixed rpm), 1 = Engine (slaved to engine rpm)
spd|ratio     if ctrl==0: compressor Speed (rpm); if ctrl==1: Speed ratio (comp-rpm / engine-rpm)
Pin  Tin      Inlet (suction) Pressure (bar), Inlet Temperature (°C)  — set intake density only
C1 C2 C3      Flow-equation coefficients      Q[l/s] = C1 + C2·P + C3·(SR·N)
D1 D2 D3      Temperature-equation coeffs      T_out[°C] = D1·P² + D2·P + D3
E1 E2 E3      Power-equation coeffs (part 1)   Power[W] = E1·P³ + E2·P² + E3·P + E4 + E5·exp(E6·SR·N)
E4 E5 E6      Power-equation coeffs (part 2)
Y0 … Y(n-2)   species mass fractions (n-1 values; the last species is the 1-complement filler)
```

where **P** = discharge (duct-node) pressure in **bar**, **N** = compressor rpm, **SR** = 1 (Proper) or
the entered ratio (Engine). Mass flow = `Q · ρ_in/1000`, `ρ_in = Pin/(R·Tin)`. Flow is **always into**
the single attached duct (the suction is an implicit atmosphere at Pin/Tin — no second pipe). Power is a
reported output only. Reduced example (Proper, N=5000, verified by flow-bench): `13  1  0  5000  1.013 25
10 -2 2e-4  0 8 30  0 0 1000 500 0 0  0.0 1.0`. **Note the UI↔file order swap:** the UI shows
P/T/coeffs first and the speed controls last, but the file stores `ctrl` and `spd/ratio` **before** Pin.
Outputs available: Power (W) / Mass Flow (kg/s) / Output Pressure (bar).

---

## 5. Connectors — BLOQUE_V `VALVULA`

Each `<VALVULA id='n' type='…'>` is a **connection element** that a UNION-A-DEPOSITO or valve node points
at. Types seen:

| `type=` | Meaning | First-line fields |
|---|---|---|
| Estator turbina | turbine **stator** nozzle | `objType Cd_in Cd_out Ø_ref(m)` |
| Rotor turbina | turbine **rotor** nozzle | `objType Cd_in Cd_out Ø_ref(m)` |
| Valvula leva | **cam valve** (intake/exhaust) | `… Ø_ref  Npts  mult  openingAngle(deg)  Ø…` then the lift table + a `Cd_in/Cd_out/Torb` flow table |
| Valvula mariposa | **throttle / butterfly** | `… 0 Ø_ref(m)` then lift (`1.0`=open) + controlled-flag |
| Valvula waste-gate | **mechanical spring wastegate** | Like all valves, line starts with its **type code `6`**, then **13 data values**: `ctrlDuct dist P_capsule C1 C2 mass damp k preload A_diaf A_plate Ø_ref mode` (internal units: bar, m, m²) — **14 numbers total**. `ctrlDuct` is a duct **index** (integer) = the post-throttle intake. See the AllFields "Wastegate Valve" entry for meanings/values. **Cracks at ≈ P_capsule + preload/A_diaf.** |
| Nodo a deposito | **node-to-deposit** (lossless connector) | `0 1.0 1.0 0` |

- The **cam valve's opening angle** (e.g., `3.5e2`=350° intake, `1.3e2`=130° exhaust) is the crank angle
  it opens — see `[[openwam-4stroke-angle-convention]]`.
- A **throttle** whose 2nd-from-last block is `1 \n 0 1` is **controlled** (a controller actuates it); a
  trailing `0` means fixed. In the example, valve 6 (WG throttle) is controlled, valve 5 (intake
  throttle) is fixed-open.
- **Node-to-deposit** is the neutral connector used to attach a duct to a 0-D volume when no valve/throttle
  is physically there (e.g., airbox inlet, catalyst in/out, throttle-chamber outlet).

---

## 6. 0-D volumes — BLOQUE_VI `DEPOSITO`

Opens with the deposit count + a few category counts. The `type=` string carries a leading **deposit-type
code** (first integer of the body): `Volumen constante`=0, **single-entry turbine**=2, **twin-scroll
turbine**=3, **Venturi**=4, **Union direccional**=5. (Both turbine kinds use the same `type='Turbina 1'`
string — the 2 vs 3 code distinguishes them.) `type=` seen:
- **`Volumen constante`** — a plenum/chamber. Body: `0 / 0.0 1.0 / volume_m3  P_bar  T_C`.
- **`Venturi N`** (code **4**) — a plenum + venturi throat model. Body after the plenum init
  (`4 / 1 / 0.0 1.0 / volume_m3 P_bar T_C`) is a venturi line:
  `numid  nodeIn nodeOut nodeLateral  areaRatio  efficiency  heatLoss` (`TVenturi.cpp:92-94`). `numid` is
  a legacy WAMer id. The three nodes are **main inlet, main outlet, lateral(suction) inlet** (each a
  UNION-A-DEPOSITO node) → "2 inputs + 1 output". **`areaRatio` is the UI-hidden section ratio and the UI
  writes it `0.0`, which makes the venturi INERT** (throat calc gated `>0`); set it >0 by hand to engage
  the throat. `efficiency` 0–1; `heatLoss` small (see §10). Example: `1 5 7 6 0.0 0.3 0.4`.
- **`Union direccional N`** (code **5**) — a plenum that merges **two inlets → one outlet** with an
  anti-backflow rule. Body after plenum init (`5 / 0.0 1.0 / volume_m3 P_bar T_C`) is a junction line:
  `numid  in0 in1 out  CDmax0 cut0 end0  CDmax1 cut1 end1` (`TUnionDireccional.cpp:107-112`). The three
  nodes are inlet-0, inlet-1, outlet (UNION-A-DEPOSITO nodes). Per inlet: max reverse-flow discharge Cd
  (0–1), cut speed and end speed (m/s) — the reverse Cd ramps `CDmax→0` between cut and end of the *other*
  inlet's velocity. **cut ≠ end** (slope divides by `cut−end`). Example: `2 8 9 10 0.1 0.2 0.3 0.4 0.5 0.6`.
- **`Turbina 1`** (code **2** single-entry, **3** twin-scroll) — the turbine's volume + init P/T, then a
  **turbine block** whose first digit is the model (**0 Fixed / 1 2-Nozzle / 2 Map**). A 2-Nozzle turbine
  carries a polynomial efficiency (max-efficiency BSR, max BSR, max η) + a ref diameter; a **Map** turbine
  embeds a full swallowing map instead (and hides the efficiency tab). See **§7** for the full
  field-by-field layout of both. The stator/rotor are the valves that UNION-A-DEPOSITO nodes attach — a
  **single-entry** turbine (code 2) has **one stator + one rotor**; a **twin-scroll** turbine (code 3,
  `TTurbinaTwin`) has **two stators (one per scroll) + one rotor**. The twin uses the same single map for
  both scrolls (each stator gets half the effective area) and sums their work; see the AllFields
  "Twin-Scroll Turbine" entry. Init the free turbo shaft speed near its equilibrium to avoid a startup
  transient.

A deposit is reached **only** through nodes that name it. Example deposits: 1 = turbine, 2 = airbox
(= compressor suction, named in BLOQUE_XIII), 3 = waste-gate chamber, 4 = throttle body, 5 = catalyst.
> Field-order caveat: const-vol init reads as `volume P_bar T_C` for deposits 3/4/5, but deposit 2 reads
> `1.5e-3 60 1.0` which is only sensible as `volume T_C P_bar` — i.e. its **P/T look swapped**. Worth a
> check when editing (init state; the compressor drives the suction anyway).

### 6a. EGR loop — no dedicated component, just plumbing
There is **no EGR object**. EGR is (1) a **species flag** + (2) an ordinary **pipe/valve loop**:
- **Flag:** engine "Working Conditions → EGR = yes" sets `ThereIsEGR` (header line 8, token 4 → `1`). This
  only adds the burnt-gas "EGR" species so recirculated exhaust is tracked; it is physically inert alone.
- **Loop (built by hand from normal elements):** tap an **exhaust-manifold plenum** → a (preferably
  **cooled**, i.e. heat-transfer-enabled) **pipe** → a **throttle** (`Valvula mariposa` = the EGR valve,
  needs a ≥2-pt Cd table) → an **intake-manifold plenum**. Each duct↔plenum join is a **Node-to-Deposit**
  (§4/§5). A convenient pattern (validated on `FILIPI_ENGINE_VGT_EGR.WAM`): insert a small const-vol
  plenum as an **exhaust collector** between the manifolds and the turbine, and tap the EGR pipe off it —
  the main flow still reaches the turbine while a fraction bleeds to EGR. Exhaust (≈2.3 bar) sits above
  intake (≈2.1 bar), so the natural Δp drives flow; the throttle opening sets the rate.
- **Behaviour (validated):** mass conserved per species to machine precision across the GJM junctions;
  EGR displaces fresh air (AFR/O₂ down, power down, BSFC up); with **ACT** combustion the inducted EGR
  lowers flame temp → **exhaust NOx decreases** (small at full-load/low-EGR, larger at part load). No
  measured EGR/NOx data in the Filipi papers, so this is a **capability demonstration**, not a validation.
> Control is optional: the papers' strategy is "VGT holds Δp_exh−int, EGR valve modulates" — reproducing
> that needs a Δp signal (a sensor per manifold + a differencing controller) or a PID on the EGR valve to a
> target EGR%. A fixed throttle opening is fine for a demonstration.

---

## 7. Turbo — BLOQUE_VII `COMPRESOR`, the turbine deposit map, BLOQUE_IX `TURBOGRUPO`, BLOQUE_XIII

### Compressor (BLOQUE_VII)
- Header: `0 1 2` (flags/instance) · a `format` digit (**0 = legacy WAMer map, 1 = SAE map**) · a flag ·
  ref `P_bar T_C` · the `massMult CRMult effMult` multipliers · `Npoints` · then the map rows.
- **Map rows: `Speed(rpm)  MassFlow(kg/s)  PressureRatio  Efficiency`** — 4 columns, grouped into speed
  lines (rows sharing the same first value). Verified row: `4.610000e+04 8.012000e-02 1.211130e+00
  6.500000e-01` = 46 100 rpm, 0.0801 kg/s, PR 1.211, η 0.65.
- **SAE-map requirements** (`[[prefer-sae-compressor-map-format]]`): each speed line must reach a **PR = 1**
  choke point, and the header **mass multiplier must be non-zero**.
- **Discharge** = the tube on its `type=17` node; **suction** = the deposit named in BLOQUE_XIII.

### Turbine map — a `DEPOSITO type='Turbina 1'` with model = **Map**
The turbine is a 0-D deposit: after the usual deposit species line + `volume_m3 P_bar T_C` init line comes
the **turbine block**. Its first digit is the **turbine-model** (`TTurbina.cpp` `tipoturb`):
**0 = Fixed, 1 = Variable-geometry / 2-Nozzle, 2 = Map.** For **Map** the block reads, in order:
```
2                          <- tipoturb = Map
6.000000e-02               <- wheel diameter (m)                    [FDiametroRodete]
5.2e-2 1.6e-2 4.5e-2 60    <- rotorOut  hub(nut)  inlet diam (m) + critical vane angle (deg)
1                          <- NumPositions          (== 1  ⇒  FIXED turbine, TTurbineMap.cpp:80)
28 70 5.300000e+01         <- Nrows   position(INTEGER!)   vane_angle(deg)      (per position)
   <Nrows data rows>:  Speed   ExpansionRatio   MassFlow   Efficiency
0 7.000000e+01             <- numctrl (=0) then FRack (fixed rack, since uncontrolled)
1.000000e+00               <- FAjustRendTurb (efficiency multiplier, applied on top of the map η)
```
- **Row columns = `Speed  ExpansionRatio  MassFlow  Efficiency`** (`TTurbPosition::ReadTurbinPosition:52`).
  **Rows MUST be grouped by ascending Speed, and ER strictly increasing within each speed group.** A row
  whose Speed < the previous is **silently dropped**; ER not increasing within a group is dropped too.
  A mis-grouped/short map → *"Hermite interpolation table has 1 x-point"* crash in `LeeRendimientoTurbina`.
- **MassFlow is the corrected/reduced flow `ṁ·√T₀₀ / p₀₀ × 10⁶`** (the code multiplies the column by
  `1e-6` in `TIsoSpeedLine::EffectiveSection`). Small-diesel range ≈ 5–35; at a ~0.22 kg/s, 890 K,
  1.95 bar operating point ≈ 34. **Smaller MassFlow at a given ER = a smaller / more-restrictive turbine =
  MORE boost** — this is *the* swallowing/boost knob (efficiency does not set swallowing).
- **Efficiency lives IN the map** for Map turbines: the UI's separate **Turbine→Efficiency tab disappears**
  when model = Map (see `AllFields.md`). The 2-Nozzle models instead carry a polynomial
  `rcOpt rcMax etaMax` + multiplier (the `0.7 1.4 0.72` / `1.0` block on a non-map turbine deposit).
- **`.TMP` import file** = the same layout minus the deposit wrapper: `NumPositions` / per-position
  `Nrows pos angle` / rows. Sample: `Notes/SampleFiles/SampleTurbineMap…TMP`. An `Adiab`/measure-temp line
  precedes `NumPositions` **only if built with `tchtm`** — absent in this build (file starts with NumPositions).
- **⚠ The rack position must be an INTEGER.** The importer reads it with `%d`; writing it as a float
  (`70.0`) leaves a stray `.0` that is misread as the vane angle and **shifts every subsequent number by one
  column**, corrupting the whole table. Symptom: a phantom first row `<angle> <speed> <ER> <MF>` (e.g.
  `53 40 1.05 6.888`) and the header vane-angle forced to `0`. Always emit the position as a bare integer.

### Shaft + suction link
- **BLOQUE_XIII** `2 2` → *"compressor 1's boundary is deposit 2"* — the suction chamber. A compressor's
  inlet must be a const-vol deposit.
- **`TURBOGRUPO`** (BLOQUE_IX): `speed0  variacion  inertia  flags…` then the coupled compressor id.
  `variacion` 0 = variable-speed (power balance solved, inertia used), 1 = imposed speed. Couples the
  compressor to the turbine deposit → the mechanical shaft.

---

## 8. Control — BLOQUE_S `SENSOR`, BLOQUE_C `CONTROLADOR`/`TABLA1D`

This is where the "engine speed → table → PID setpoint" chain lives. Verified from `TSensor.cpp`,
`TPIDController.cpp`, `TTable1D.cpp`, `TOpenWAM.cpp::ReadControllers`.

### Sensors — `obj prm [objID] [dist] delay gain`
| `obj` | object | `prm` values |
|---:|---|---|
| 0 | execution | 0 = time |
| 1 | **pipe** (reads `objID`=tube, then `dist`=position 0–1) | 1 = pressure, 2 = temperature, 3 = mass flow |
| 2 | **deposit** (reads `objID`) | 1 = pressure, 2 = temperature |
| 3 | **engine** | 4 = fuel, 5 = **engine speed** |

Examples: `3 5` = engine-speed sensor; `1 1 2 0.01` = pressure in **pipe 2** at x=0.01.

### Controllers — count, then per controller a **type digit** (`1 PID, 2 Table1D, 3 Decisor, 4 Gain`)
- **PID** (`TPIDController::LeeController`): `Kp+ Ki+ Kd+  Kp- Ki- Kd-  out out0 max min  period delay
  gain  setPoint  setPointCtrlID  sensorID`.
  - `setPointCtrlID`: **0 = fixed setpoint**, else = **the controller whose Output() is the setpoint**
    (this is how a table drives a PID).
  - `sensorID` = the **feedback / process-variable** sensor.
- **Table1D** (`TTable1D::LeeController`): `fromfile  [xnum + (x,y)…]  period  interpType(0 lin/1
  Hermite/2 steps)  sensorID`. `Output()` returns `interp(Sensor[sensorID-1].Output())` → **`sensorID`
  is the table's input**.

All ids are **1-based** (`Sensor[ID-1]`, `Controller[ID-1]`).

### The example's boost-control loop
```
[4T engine] --speed--> Sensor 1 (obj3/prm5)        [Duct 2] --pressure--> Sensor 2 (obj1/prm1, pipe2)
                          │  (intended table input)                          │
                          ▼                                                   ├── feedback ─┐
                     Table 1 (ctrl 2, rpm→boost) ── setpoint ─▶ PID (ctrl 1) ◀─────────────┘
                                                                    │
                                                                    ▼  actuates
                                                             WG throttle (valve 6)
```
> **Wiring gotcha this file demonstrates:** the sensors were re-ordered (speed→id 1, pressure→id 2). The
> PID's `sensorID` was updated 1→2 (correct: feedback = pressure), **but the Table's `sensorID` was left
> at 2** — so the table now interpolates against the *pressure* sensor, not engine speed. **Fix: Table's
> `sensorID` 2 → 1.** Lesson: when you renumber sensors, grep **every** controller/table `sensorID` and
> `setPointCtrlID`, because they are bare 1-based indices with no symbolic tie to the sensor's meaning.

---

## 8a. Space-time results — BLOQUE_XII (`RESULTADOS ESPACIO-TEMPORALES`)
The x–t wave-map request. Format (`TOutputResults::ReadSpaceTimeResults`):
```
Nmag                     first int = number of magnitudes; if 0, the block is empty (nothing read/output)
Nc Nd Nt                 counts: cylinders, deposits(plenums), tubes  (only present when Nmag>0)
<cylIDs...>              Nc ids   (1-based; blank line if Nc=0)
<depIDs...>              Nd ids   (blank line if Nd=0)
<tubeIDs...>             Nt ids
<paramIDs...>            Nmag ids: 0=pressure 1=temperature 2=velocity 3=mass-flow 4=species
```
Example (`0 0 2 / (blank) (blank) / 1 22 / 0 1 2 3`) = 4 magnitudes on pipes 1 and 22. **Both the
element list and ≥1 magnitude are required** — the UI writes `Nmag=0` (empty) if you tick only the
per-element "space-time" boxes without choosing magnitudes in the Execution node, and then nothing is
output. Selecting a **pipe** is what gives the spatial axis (one column per cell); a cylinder/plenum
contributes a single scalar column.

**Output:** one file per magnitude beside the INS file — `<name>_pre.DAT` / `_tem.DAT` / `_vel.DAT` /
`_air.DAT` (+ species `_YGQ`/`_YAF`/…). **Tab-separated** (like AVG/INS, so it opens into columns): line 1
is a **labelled header row** (`Angle(deg) <qty>_Cyl<n> … <qty>_Pipe<n>_cell<k> …`, e.g. `P_bar_Pipe1_cell0`),
then data rows `angle │ [cyls] │ [plenums] │ [pipe cells]`, one row per crank-angle step of the **last**
cycle. (The old numeric metadata block — `Nc Nd Nt`, id lists, cell counts, row count — was dropped
2026-08-25; the header row supersedes it.)
> **Gating fix (2026-08-24):** the engine write-gate used `Run.CycleDuration` (=720°, angle-per-cycle) as
> if it were the cycle count, so it only triggered on cycle 719 — effectively never. Now uses the run's
> `SimulationDuration` (cycle count) and writes the last cycle (`TOpenWAM.cpp:3829`,
> `TOutputResults::WriteSpaceTime`). See AllFields "Space-time results" for the UI steps.

---

## 9. Writing a WAM from an engine architecture (procedure)

1. **Sketch the graph**: list every duct (edge) and every junction / connector / boundary / 0-D volume
   (vertex). Decide cylinder count and firing order.
2. **Assign node ids** to every vertex; **assign tube ids** to every edge. Keep a table.
3. **BLOQUE_III**: engine, cylinders, Wiebe law, firing order.
4. **BLOQUE_IV**: for each duct write `nodeL nodeR …`, the duct type digit (pipe vs port), and
   `Ø_in length Ø_out` in **metres**. Intercooler → set N-intercooler-ducts > 1 + variable wall + layers.
5. **BLOQUE_V**: define every connector — stator/rotor for a turbine, cam valves for intake/exhaust,
   throttles, node-to-deposits. Mark controlled throttles.
6. **BLOQUE_VI**: define deposits (plenums, turbine, chambers). Give each a volume + init P/T.
7. **BLOQUE_VII/IX/XIII**: compressor (+map, SAE preferred — `[[prefer-sae-compressor-map-format]]`),
   turbogroup (couple to turbine deposit), and the compressor→suction-deposit link.
8. **BLOQUE_VIII**: for every vertex emit a `<NODO>` with the right **type code** (§4) and payload
   (`self deposit valveDef` / `self cyl valveDef` / `compId` / boundary state). **Cross-check**: each
   node id you used as a tube endpoint must exist here, and each `deposit`/`cyl`/`valveDef`/`compId` must
   exist in its block.
9. **BLOQUE_S/C**: sensors + controllers; wire `sensorID` and `setPointCtrlID` (1-based). Double-check
   the indices point at the intended sensors (see §8 gotcha).
10. **COORDENADAS**: x/y per element for the UI (cosmetic; can be rough).

**Validation pass (do this every time):**
- Every `TUBO` endpoint appears as a `NODO`. Every `NODO` element ref (deposit/cyl/valve/comp) exists.
- Each 0-D volume is reached through the right connector count (a compressor suction needs a
  node-to-deposit; a throttle needs its own chamber + an outlet connector; a turbine needs stator in +
  rotor out).
- Sensor/controller indices resolve to the intended objects.
- Counts at the top of each block equal the number of entries you wrote.

---

## 10. Gotchas checklist

- **Tags are landmarks, numbers are data** — never add/remove a number without knowing the field.
- **`pipa` = PORT, `tubo` = PIPE** (`[[duct-type-pipa-port-tubo-pipe]]`).
- **Lengths/diameters in metres** in-file, mm in the UI (`[[openwam-file-units-meters-not-mm]]`).
- **0-D elements are reached only through nodes** — a duct never names a deposit/cylinder directly.
- **Manifolds may be pure ramifications** (no plenum) — the example has none.
- **Compressor suction = a deposit** named in BLOQUE_XIII; discharge = the tube on its node.
- **Control indices are bare 1-based integers** — renumbering sensors silently breaks tables/PIDs unless
  you update every `sensorID`/`setPointCtrlID`.
- **SAE compressor maps** need each speed line to reach **PR = 1** and use a **non-zero mass multiplier**
  in the header (`tchtm` off → 6-field header `refP refT massMult CRMult effMult Npoints`). Rows are
  `Speed(rpm) MassFlow(kg/s) PR Eff` (§7).
- **Turbine maps**: rows are `Speed ER MassFlow Eff` **grouped by ascending Speed, ER increasing within a
  group** — mis-order and rows are silently dropped → Hermite crash. `NumPositions==1` ⇒ fixed turbine.
  MassFlow is reduced flow `ṁ√T₀₀/p₀₀·1e6` and is the swallowing/boost knob; efficiency is **in** the map
  (no efficiency tab for Map turbines). See §7.
- **Turbine-map rack position is an INTEGER** — a float like `70.0` shifts every value by one column and
  corrupts the whole map (phantom `<angle> <spd> <ER> <MF>` first row). Same class of bug as the control-
  index gotcha: one wrong token cascades through the fixed-position stream.
- **Pressure-loss `FK` must be > 0** (node types 9/10) — the solver brackets with `1/√FK`, so `FK=0`
  divides by zero. Order-unity values are physical (loss ∝ `FK·|U|` laminar, `FK·U²` turbulent).
- **Venturi area-ratio `0.0` = inert** (UI writes it) — the throat model is gated `ratio>0`; set >0 to
  engage. When active, keep `heatLoss` small: throat velocity uses `√(eff·V²−2·heatLoss)`.
- **Directional-junction `cut ≠ end`** — the reverse-Cd slope divides by `(cut−end)`; equal ⇒ inf. Use
  cut/end in m/s bracketing the real flow; in steady co-flow the rule is dormant anyway.
- Cross-reference `AllFields.md` for exhaustive per-field UI meanings; this doc is the wiring map.

---

## 11. Known crash conditions (silent segfaults during setup)

OpenWAM often **segfaults (0xC0000005 / SIGSEGV) with no error message** on malformed input, and stdout is
block-buffered so the last prints are lost — the crash looks like it's "right after the last line printed"
but is usually later. To localize: the reads run in the order `Pipes → Valves → Plenums → Compressors →
Connections(nodes) → TurboAxis → Sensors → Controllers`; the pipe INFO lines flush (they use `endl`), the
rest don't. A quick way to pin it down is to add flushed `std::cout << … << std::endl;` markers before each
read call in `TOpenWAM::ReadInputData` (and inside `ReadConnections`' node loop), rebuild, and read the last
marker. (gdb on the release build is useless — MSVC build is stripped, `?? ()`.)

| Symptom | Cause | Fix |
|---|---|---|
| Segfault in `ReadConnections` at the first **throttle→deposit** union | A **throttle (`Valvula mariposa`) with an empty Cd-vs-lift table** (`FNumLev = 0`, i.e. valve data reads `10 0 <Ø> …`). `TMariposa`'s copy ctor builds `Hermite_interp` on empty vectors (`TMariposa.cpp:64`). | Give every throttle a **Cd-vs-lift table with ≥2 points** (e.g. lift 0→Cd 0, lift 1→Cd 0.9). A working throttle reads `10 2 <Ø>` + rows `lift CDin CDout`. |
| Segfault after setup when a **screw-compressor** node requests any output | The UI writes the screw-compressor index (`cvId`) **0-based**, but the lookup-array linking required 1-based, leaving `VolumetricCompressor[0]` null; the output setup then dereferenced it (`TOutputResults.cpp:671`). | **Fixed 2026-08-22** (`TOpenWAM.cpp:1523` now links sequentially). Older builds: hand-edit `cvId` from 0 to 1 in the node. |
| Screw-compressor run **hangs** (timestep collapses ~1 step in, progress frozen) | `stComprVol::operator()` shadowed its `entropia` member with a local, so the boundary wrote **NaN entropy** into the duct on step 1. | **Fixed 2026-08-22** (`BoundaryFunctions.h`: assign the member, no local). |
| Segfault on the first step with a **mechanical wastegate** (`Valvula waste-gate`) | `TWasteGate::CalculoNodoAdm` (which wires the control-duct pointer `FTuboAdm`) was never called, so `GetCDin/GetCDout` dereferenced a null `FTuboAdm`. | **Fixed 2026-08-22** (`TOpenWAM.cpp:1597`: wired when the wastegate BC list is built). |
| Run dies mid-simulation with `ERROR : plenum n. N too small` on a **mechanical wastegate** model | A **zero wastegate mass** (`FMasa`, the 6th data value): `FddX = (F − k·x − c·ẋ)/FMasa` divides by zero → NaN valve → exhaust plenum empties. WAMer **zeroes `FMasa` when it regenerates a parametric-study case `.WAM`**, so every swept case of a wastegate model hits this. | **Guarded 2026-08-26** — `TWasteGate::LeeDatosIniciales` now throws `ERROR: WasteGate #n has zero/negative mass` at read time. Restore the wastegate line in the generated `.WAM`s, or sweep a non-wastegate model. |

**Build/debug crib:** the parallel exe is MSVC + Ninja (x86). Rebuild after a source edit with the VS env
sourced: `call ".../VC/Auxiliary/Build/vcvars32.bat" && cmake --build build-parallel` (Ninja does an
incremental relink; ~30 s for one changed file). Run a diagnostic copy from a **space-free path** (the
importer writes `tmp.wam` in the cwd) with `OMP_NUM_THREADS=1` for a clean single-threaded trace.
```
