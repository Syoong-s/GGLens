#!/bin/sh
set -eu

mode=$1
python_prefix=$2
patch_tool=$3
binary=${4-}
soname=libstdc++.so.6

unavailable() {
    if [ "$mode" = auto ]; then
        if [ -n "$binary" ]; then
            printf 'C++ runtime kept from the original link. %s\n' "$1" >&2
        fi
        exit 0
    fi
    printf 'C++ runtime selection failed. %s\n' "$1" >&2
    exit 1
}

case "$mode" in
    system) exit 0 ;;
    auto|/*) ;;
    *) unavailable 'Set CXX_RUNTIME to auto, system, or an absolute path to libstdc++.so.6.' ;;
esac

[ "$(uname -s)" = Linux ] || unavailable 'Selecting an explicit C++ runtime requires Linux ELF.'

library=$mode
if [ "$mode" = auto ]; then
    case "$python_prefix" in
        /*) ;;
        *) unavailable 'PYTHON_CONFIG did not return an absolute Python prefix.' ;;
    esac
    for directory in "$python_prefix/lib" "$python_prefix/lib64"; do
        library="$directory/$soname"
        if [ -f "$library" ] && [ -r "$library" ]; then break; fi
    done
    [ -f "$library" ] && [ -r "$library" ] || unavailable "No readable $soname under $python_prefix/lib or lib64."
fi
[ -f "$library" ] && [ -r "$library" ] || unavailable "No readable C++ runtime at $library."

patcher=$(command -v "$patch_tool" || :)
loader=$(command -v ldd || :)

# Without a binary, record inputs that can change the outcome on the next build.
if [ -z "$binary" ]; then
    printf '%s\n' "$library" "$patcher" "$loader"
    stat -Lc '%d %i %s %y' -- "$library" || :
    exit 0
fi

[ -n "$patcher" ] || unavailable "Install patchelf or set PATCHELF to use $library."
[ -n "$loader" ] || unavailable 'ldd is required to check the linked dependencies.'
[ "$("$patcher" --print-soname "$library" 2>/dev/null)" = "$soname" ] ||
    unavailable "$library does not declare $soname as its ELF SONAME."
needed=$("$patcher" --print-needed "$binary") || unavailable 'Could not inspect the linked ELF dependencies.'
printf '%s\n' "$needed" | grep -Fxq "$soname" || unavailable "The linked ELF has no direct $soname dependency."

dependency_report() (
    unset LD_PRELOAD
    LC_ALL=C "$loader" "$@" 2>&1
)

# Extract resolved paths from the loader's address lines, including absolute DT_NEEDED entries.
resolved_paths() {
    awk '/[[:space:]]+\(0x[[:xdigit:]]+\)[[:space:]]*$/ {
        sub(/^[[:space:]]+/, "")
        sub(/[[:space:]]+\(0x[[:xdigit:]]+\)[[:space:]]*$/, "")
        sub(/^[^/]* => /, "")
        print
    }'
}

mpi_and_python_paths() {
    grep -E '/lib(mpi[^/]*|open-(pal|rte)|python[^/]*)\.so(\.|$)' | LC_ALL=C sort
}

original=$(dependency_report "$binary") || unavailable 'Could not resolve the original ELF dependencies.'
patched=$(mktemp "$binary.cxx.XXXXXX") || unavailable 'Could not create a temporary ELF for validation.'
trap 'rm -f -- "$patched"' EXIT
trap 'exit 1' HUP INT TERM
cp -p -- "$binary" "$patched" || unavailable 'Could not copy the linked ELF for validation.'
"$patcher" --replace-needed "$soname" "$library" "$patched" || unavailable 'patchelf could not replace the C++ runtime dependency.'
if ! report=$(dependency_report -r -v "$patched"); then
    printf '%s\n' "$report" | sed -n '1,8p' >&2
    unavailable "The loader could not inspect $library."
fi
errors=$(printf '%s\n' "$report" | grep -E '=> not found$|: version .* not found \(required by |: no version information available \(required by |^[[:space:]]*undefined symbol:' || :)
if [ -n "$errors" ]; then
    printf '%s\n' "$errors" | awk '!seen[$0]++' | sed -n '1,8p' >&2
    unavailable "$library did not satisfy the linked library requirements."
fi

paths=$(printf '%s\n' "$report" | resolved_paths)
if ! printf '%s\n' "$paths" | SELECTED_CXX_LIBRARY="$library" awk '
    $0 == ENVIRON["SELECTED_CXX_LIBRARY"] { ++selected; next }
    /(^|\/)libstdc[+][+]\.so\.6(\.|$)/ { ++other }
    END { exit !(selected == 1 && other == 0) }
'; then
    unavailable "The loader did not select $library as the only C++ runtime."
fi
original_paths=$(printf '%s\n' "$original" | resolved_paths)
if [ "$(printf '%s\n' "$paths" | mpi_and_python_paths)" != "$(printf '%s\n' "$original_paths" | mpi_and_python_paths)" ]; then
    unavailable 'The candidate changes the resolved MPI or libpython libraries.'
fi

mv -f -- "$patched" "$binary"
printf 'C++ runtime pinned to %s\n' "$library"
