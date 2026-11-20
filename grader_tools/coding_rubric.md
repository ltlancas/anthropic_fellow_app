# Code Simplicity Rubric

This rubric grades a code base on **simplicity**: whether the code does what the task requires with the least machinery, and whether a competent reader (e.g. a computational astrophysicist with Python experience) can understand, check, and change it quickly.

Simplicity does **not** mean "fewest lines at any cost". Dense one-liners, clever tricks, and hidden magic count against simplicity. The goal is code where the path from the physics or specification to the implementation is short and obvious.

Correctness is graded elsewhere. Only grade correctness here when complexity is hiding it (e.g. code too convoluted to check).

---

## Scoring Overview

Score each of the 8 criteria below from **0 to 4**, then add them up for a total out of **32**.

| Score | Meaning |
| --- | --- |
| 4 | Excellent: nothing meaningful to simplify |
| 3 | Good: minor, isolated issues |
| 2 | Adequate: several issues; the code works, but a reader has to fight it in places |
| 1 | Poor: complexity is widespread and gets in the way of understanding |
| 0 | Unacceptable: the code base is mostly incidental complexity |

| Total | Overall Grade |
| --- | --- |
| 28–32 | A: simple, clear, minimal |
| 22–27 | B: mostly simple, a few unnecessary complications |
| 15–21 | C: works, but noticeably over-engineered or disorganized |
| 8–14 | D: hard to follow; significant rewrite needed for clarity |
| 0–7 | F: complexity dominates |

---

## Criteria

### 1. Scope: Only What Was Asked (0–4)

Does the code implement the required functionality and no more?

- **4:** Every module, function, and option serves a stated requirement. No speculative features.
- **2:** Some unused options, dead code paths, or "future-proofing" (e.g. support for configurations never exercised by the task).
- **0:** Large parts of the code base are unused, speculative, or unrelated to the task.

Red flags: unused functions or parameters; commented-out code blocks; plugin/registry systems for one implementation; configuration flags that are never toggled; support for geometries, equations of state, or solvers the task never asked for.

### 2. Structure and Abstraction (0–4)

Are abstractions earned, and does the file and module layout match the size of the problem?

- **4:** Abstractions match real concepts in the problem (e.g. grid, state, Riemann solver, time integrator). Plain functions are used where a class adds nothing. Layout is flat and easy to navigate.
- **2:** Some unnecessary layers: classes that wrap a single function, base classes with one subclass, getters and setters for plain data, many tiny files.
- **0:** Deep inheritance, factories, managers, or indirection that obscures what actually runs.

Ask: *If an abstraction were inlined, would the code get easier to read?* If yes, the abstraction is not earned.

### 3. Readability and Naming (0–4)

Can the code be read top to bottom and understood?

- **4:** Names are descriptive and match the physics/specification (`rho`, `pressure`, `r_fb`, `Mdot_w`). Functions are short and do one thing. Control flow is linear and obvious.
- **2:** Some cryptic names, long functions mixing several responsibilities, or nested logic that takes effort to trace.
- **0:** Opaque names (`tmp2`, `arr`, `f3`), very long functions, deeply nested branches, or "clever" constructs that need careful decoding.

### 4. Direct Mapping from Math to Code (0–4)

For scientific code: is it easy to see which equation each block implements?

- **4:** Key equations are implemented in a form that visibly matches their mathematical statement, with a short comment or reference (equation number, paper) where it helps. Units are handled in one obvious place and are consistent.
- **2:** The mapping is recoverable but requires effort; unit conversions are scattered or implicit; some unexplained magic numbers.
- **0:** Equations are fragmented, rearranged beyond recognition, or full of unexplained constants; unit handling is ad hoc.

Red flags: hard-coded numerical constants with no name or comment; the same physical constant defined in several places; mixed unit systems inside a single function.

### 5. No Duplication (0–4)

Is each piece of logic written once?

- **4:** Shared logic (e.g. the flux calculation, unit conversions, output writing) lives in one place. First- and second-order variants share everything they reasonably can and differ only where the method actually differs.
- **2:** Some copy-pasted blocks with small edits; parallel functions that drift apart.
- **0:** Large-scale duplication (e.g. separate near-identical solvers or scripts per experiment).

Note: do not reward premature de-duplication. Two short, similar blocks can be simpler than a generic helper with many flags. Penalize duplication only when it makes the code harder to maintain or check.

### 6. Dependencies and Tooling (0–4)

Is the dependency footprint as small as the task allows?

- **4:** Only the standard library plus what is needed (here, `numpy`/`scipy`). Environment setup is a single short, reproducible step.
- **2:** One or two unnecessary dependencies, or a setup process with extra steps.
- **0:** Heavy or unnecessary frameworks, complex build systems, or an environment that is hard to reproduce.

### 7. Ease of Running (0–4)

Can a new user run the code and reproduce the outputs with minimal effort?

- **4:** One obvious entry point (or one per task) with a clear command. Parameters are explicit and in one place. Outputs land in the expected locations without manual steps.
- **2:** Running requires reading the source to find the right script, editing files by hand, or running steps in an undocumented order.
- **0:** Unclear how to run at all, or reproduction requires significant manual intervention.

### 8. Tests and Documentation Proportionate to the Code (0–4)

Are tests and docs short, focused, and useful, rather than absent or bloated?

- **4:** A small number of tests that check meaningful behaviour (e.g. conservation, a known analytic solution, a shock tube) and run quickly. Documentation explains structure and usage concisely. Comments explain *why*, not *what*.
- **2:** Tests are trivial or redundant, or documentation is padded and restates the code. Comments narrate obvious lines.
- **0:** No usable tests or docs, or tests and docs so sprawling that they are themselves a maintenance burden.

---

## Grading Procedure

1. **Skim the layout first.** List the files and their sizes. Note how many modules, classes, and entry points there are relative to the size of the task.
2. **Trace one end-to-end path.** For example, follow one experiment from input parameters to output file. Note every layer of indirection you pass through.
3. **Score each criterion.** For each, give the score and 1–3 sentences of justification citing specific files and line numbers (`path/to/file.py:42`).
4. **List the top simplifications.** Give up to five concrete changes that would most reduce complexity, ordered by impact (e.g. "Merge `FirstOrderSolver` and `SecondOrderSolver` into one function with an `order` argument; they share ~80% of their code").
5. **Report the total and overall grade.**

## Output Template

```markdown
## Simplicity Grade: <letter> (<total>/32)

| Criterion | Score | Justification |
| --- | --- | --- |
| 1. Scope | x/4 | ... |
| 2. Structure and Abstraction | x/4 | ... |
| 3. Readability and Naming | x/4 | ... |
| 4. Math-to-Code Mapping | x/4 | ... |
| 5. No Duplication | x/4 | ... |
| 6. Dependencies and Tooling | x/4 | ... |
| 7. Ease of Running | x/4 | ... |
| 8. Tests and Docs | x/4 | ... |

### Top Simplifications
1. ...
2. ...

### Summary
<2–4 sentences on the overall character of the code base.>
```
