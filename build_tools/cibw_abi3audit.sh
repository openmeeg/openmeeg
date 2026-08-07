#!/bin/bash -ef

# abi3audit only applies to the stable-ABI wheel. The free-threaded build has no
# limited API at all (CPython gh-111506), so it ships a version-specific wheel that
# abi3audit rejects outright rather than passing trivially.

WHEEL="$1"

case "$(basename "$WHEEL")" in
    *-abi3-*)
        pipx run abi3audit --strict --report "$WHEEL"
        ;;
    *)
        echo "Skipping abi3audit for non-abi3 wheel: $(basename "$WHEEL")"
        ;;
esac
