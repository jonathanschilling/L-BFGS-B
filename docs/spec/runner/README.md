# Conformance test runner

`conformance.py` reads JSON test vectors from `docs/spec/data/` and
runs each one through the Python reference implementation under
`docs/spec/reference_impl/core/`. Reports pass / fail with diagnostics
on mismatch.

## Usage (Python reference impl)

```
# Strict mode (bit-for-bit, requires reference BLAS).
python3 docs/spec/runner/conformance.py --strict

# Tolerance mode (per-subroutine numerical slack).
python3 docs/spec/runner/conformance.py --tolerance

# Filter to a single subroutine.
python3 docs/spec/runner/conformance.py --filter projgr
```

Exit status: `0` on full pass, `1` on any failure.

## Adapting the runner to another port

The runner has a per-subroutine `HANDLERS` dispatch table. Each
handler:

1. Takes the JSON `inputs` dict.
2. Constructs the subroutine call (build arrays, set parameters).
3. Calls the implementation under test.
4. Returns the actual outputs as a dict matching the `expected` dict
   in the JSON.

To validate a port in another language (C, C++, Java, ...), there are
two approaches:

### Option 1: language-native test runner

Each port writes its own runner that mirrors `conformance.py`. The
JSON format is language-neutral; the runner just needs to parse
inputs, call the port's implementation, and compare outputs.

This is the recommended approach for production-quality ports.
Reference: the Python `conformance.py` is ~250 LOC and has one
handler per subroutine; the same structure transfers.

### Option 2: subprocess wrapper

The port exposes a CLI: read JSON test-vector input from `stdin`,
write JSON output to `stdout`. The Python runner can be extended to
spawn this CLI per JSON case and parse the output:

```
echo '{"inputs": {...}}' | port-binary --subroutine projgr
# stdout: {"sbgnrm": 3.0}
```

This is useful for quick validation during initial port development.

## Tolerances

Per-subroutine tolerances for `--tolerance` mode are defined in
`conformance.py:TOLERANCES`, matching the table in
`docs/spec/07_conformance.md`. Override locally if your port needs
looser tolerances; document the override in the port's docs.

## Handler signature

Each handler is a `Callable[[dict], dict]` registered with the
`@handler(name)` decorator. The dict in is the JSON's `inputs`
field; the dict out should match the `expected` field's keys.
