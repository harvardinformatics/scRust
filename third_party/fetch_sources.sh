#!/usr/bin/env bash
# Download the exact upstream sources that scRust ports from (docs/PLAN.md A1 and B6, milestone
# M1b) into third_party/src/ (gitignored), check each tarball against the sha256 recorded in
# third_party/SHA256SUMS, and unpack each one into its own directory:
#   third_party/src/tarballs/<file>           the tarball as downloaded
#   third_party/src/<name>-<version>/         its contents, without the top-level directory
#
# These are reading references for the ports. Nothing in them is built, sourced or run.
#
# Usage, from any directory:
#   third_party/fetch_sources.sh            download what is missing and verify everything
#   third_party/fetch_sources.sh --record   also record the sha256 of sources that have none yet
#
# --record never replaces a recorded sha256. If an upstream file really changed, delete its line
# from SHA256SUMS on purpose, then record it again. The checksums can also be checked directly:
#   (cd third_party/src/tarballs && sha256sum -c ../../SHA256SUMS)
set -euo pipefail

here=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
src="$here/src"
tarballs="$src/tarballs"
sums="$here/SHA256SUMS"

record=0
case "${1:-}" in
    "") ;;
    --record) record=1 ;;
    -h|--help) sed -n '2,19p' "${BASH_SOURCE[0]}"; exit 0 ;;
    *) echo "fetch_sources.sh: unknown argument: $1" >&2; exit 2 ;;
esac

CRAN=https://cran.r-project.org
BIOC=https://bioconductor.org/packages/3.22/bioc

# name, version, tarball, where to get it. "cran" and "bioc" try the current package index first
# and then the archive, where CRAN and Bioconductor move a version once a newer one is released;
# the sha256 makes sure every URL yields the same bytes. The versions are those of the production
# R env (envs/r-reference/README.md), except where a package's compiled code includes another
# package's headers: then it is the version the binary was built against (milestone M4:
# glmGamPoi 1.22.0 and sctransform 0.4.3 were built with RcppArmadillo 15.2.3-1 and Rcpp 1.1.1,
# per their conda recipes' host requirements; the env has RcppArmadillo 15.6.2-1 installed).
SOURCES='
R              4.5.3       R-4.5.3.tar.gz                  https://cran.r-project.org/src/base/R-4/R-4.5.3.tar.gz
sctransform    0.4.3       sctransform_0.4.3.tar.gz        cran
glmGamPoi      1.22.0      glmGamPoi_1.22.0.tar.gz         bioc
matrixStats    1.5.0       matrixStats_1.5.0.tar.gz        cran
Seurat         5.5.1       Seurat_5.5.1.tar.gz             cran
SeuratObject   5.4.0       SeuratObject_5.4.0.tar.gz       cran
irlba          2.4.1       irlba_2.4.1.tar.gz              cran
KernSmooth     2.23-27     KernSmooth_2.23-27.tar.gz       cran
fields         17.3        fields_17.3.tar.gz              cran
DoubletFinder  1b244d8f0d54b4b1cb4365639931bbb16f01e1cd  DoubletFinder-1b244d8f0d54b4b1cb4365639931bbb16f01e1cd.tar.gz  https://github.com/chris-mcginnis-ucsf/DoubletFinder/archive/1b244d8f0d54b4b1cb4365639931bbb16f01e1cd.tar.gz
Matrix         1.7-6       Matrix_1.7-6.tar.gz             cran
RcppArmadillo  15.2.3-1    RcppArmadillo_15.2.3-1.tar.gz   cran
Rcpp           1.1.1       Rcpp_1.1.1.tar.gz               cran
future         1.82.0      future_1.82.0.tar.gz            cran
future.apply   1.21.0      future.apply_1.21.0.tar.gz      cran
parallelly     1.50.0      parallelly_1.50.0.tar.gz        cran
MatrixGenerics 1.22.0      MatrixGenerics_1.22.0.tar.gz    bioc
DelayedMatrixStats 1.32.0  DelayedMatrixStats_1.32.0.tar.gz  bioc
'

die() { echo "fetch_sources.sh: $*" >&2; exit 1; }

urls_for() {  # name file where
    case "$3" in
        cran) echo "$CRAN/src/contrib/$2 $CRAN/src/contrib/Archive/$1/$2" ;;
        bioc) echo "$BIOC/src/contrib/$2 $BIOC/src/contrib/Archive/$1/$2" ;;
        https://*) echo "$3" ;;
        *) die "bad source location for $1: $3" ;;
    esac
}

recorded_sha() {  # file -> recorded sha256, or nothing
    [ -f "$sums" ] || return 0
    awk -v f="$1" '$2 == f || $2 == "*" f { print $1; exit }' "$sums"
}

download() {  # name file where dest
    local url
    for url in $(urls_for "$1" "$2" "$3"); do
        echo "  trying $url"
        if curl -fsL --proto '=https' --tlsv1.2 --retry 3 --connect-timeout 30 -o "$4" "$url" < /dev/null; then
            return 0
        fi
        echo "  not available there"
    done
    rm -f "$4"
    die "could not download $2"
}

# Unpack a tarball into a new, empty directory, dropping its single top-level directory. Refuse
# member names that are absolute, contain "..", or lie outside that top-level directory.
unpack() {  # tarball dest
    local top tmp
    top=$(tar -tzf "$1" | awk 'NR == 1')
    top=${top%%/*}
    [ -n "$top" ] && [ "$top" != . ] && [ "$top" != .. ] || die "$1: no top-level directory"
    if tar -tzf "$1" | awk -v top="$top" '
            $0 ~ /^\// || $0 ~ /(^|\/)\.\.(\/|$)/ { bad = 1 }
            $0 != top && $0 != top "/" && index($0, top "/") != 1 { bad = 1 }
            END { exit !bad }'; then
        die "$1: unsafe or unexpected member names; not unpacking"
    fi
    tmp="$2.unpacking"
    rm -rf "$tmp"
    mkdir "$tmp"
    tar -xzf "$1" -C "$tmp" --strip-components=1 --no-same-owner --no-same-permissions
    rm -rf "$2"
    mv "$tmp" "$2"
}

mkdir -p "$tarballs"
new_sums=()
status=0

while read -r name version file where; do
    [ -n "$name" ] || continue
    echo "$name $version"
    tarball="$tarballs/$file"
    want=$(recorded_sha "$file")

    if [ -f "$tarball" ] && [ -n "$want" ] && [ "$(sha256sum < "$tarball" | cut -d' ' -f1)" = "$want" ]; then
        echo "  present: $file"
    else
        download "$name" "$file" "$where" "$tarball.part"
        got=$(sha256sum < "$tarball.part" | cut -d' ' -f1)
        if [ -n "$want" ] && [ "$got" != "$want" ]; then
            rm -f "$tarball.part"
            echo "  sha256 MISMATCH for $file: recorded $want, downloaded $got" >&2
            status=1
            continue
        fi
        if [ -z "$want" ]; then
            if [ "$record" -eq 1 ]; then
                echo "  recording sha256 $got"
                new_sums+=("$got  $file")
            else
                rm -f "$tarball.part"
                echo "  no sha256 recorded for $file in SHA256SUMS; rerun with --record" >&2
                status=1
                continue
            fi
        fi
        mv "$tarball.part" "$tarball"
        echo "  downloaded: $file ($(stat -c %s "$tarball") bytes)"
    fi
    sha=$(sha256sum < "$tarball" | cut -d' ' -f1)

    # A GitHub archive records its commit in the tar's pax global header ("52 comment=<sha>").
    if [ "$name" = DoubletFinder ]; then
        # (pipefail off: gzip may die of SIGPIPE once head has its 1024 bytes)
        commit=$(set +o pipefail
                 gzip -dc "$tarball" | head -c 1024 | tail -c 512 | tr -d '\0' | sed -n 's/^52 comment=//p')
        [ "$commit" = "$version" ] || die "$file: archive is of commit '$commit', not $version"
        echo "  commit verified: $commit"
    fi

    dest="$src/$name-$version"
    stamp="$src/.$name-$version.sha256"
    if [ -d "$dest" ] && [ -f "$stamp" ] && [ "$(cat "$stamp")" = "$sha" ]; then
        echo "  unpacked: ${dest#"$here"/}"
    else
        unpack "$tarball" "$dest"
        echo "$sha" > "$stamp"
        echo "  unpacked into ${dest#"$here"/}"
    fi
done <<< "$SOURCES"

if [ "${#new_sums[@]}" -gt 0 ]; then
    printf '%s\n' "${new_sums[@]}" >> "$sums"
    # Keep SHA256SUMS in the order of the source list above.
    order=$(awk 'NF { print $3 }' <<< "$SOURCES")
    awk 'NR == FNR { rank[$1] = NR; next } { print rank[$2] "\t" $0 }' <(echo "$order") "$sums" \
        | sort -n -k1,1 | cut -f2- > "$sums.tmp"
    mv "$sums.tmp" "$sums"
    echo "recorded ${#new_sums[@]} new sha256 in ${sums#"$here"/}"
fi

[ "$status" -eq 0 ] || die "some sources failed; see above"
echo "all sources present and verified in ${src#"$here"/}/"
