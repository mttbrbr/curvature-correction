#!/bin/sh
# SPDX-License-Identifier: GPL-3.0-or-later
set -eu

: "${WM_PROJECT_DIR:?runner must preload OpenFOAM.com v2412}"
: "${CURVATURE_VALIDATION_ARCHIVE:?configure the read-only validation archive mount}"

campaign_id=${1:?usage: run_private_validation.sh CAMPAIGN_ID}
case "$campaign_id" in
    *[!A-Za-z0-9._-]*|'')
        echo "ERROR: invalid campaign identifier: $campaign_id" >&2
        exit 2
        ;;
esac

repo_dir=$(CDPATH='' cd -- "$(dirname "$0")/.." && pwd)
archive_dir=$CURVATURE_VALIDATION_ARCHIVE
manifest="$archive_dir/manifest.json"
public_manifest="$archive_dir/public-manifest.json"
runner="$archive_dir/bin/run-campaign"
output_dir="$repo_dir/reports/curvature-validation/$campaign_id"
work_dir=$(mktemp -d /tmp/curvature-validation.XXXXXX)

cleanup()
{
    rm -rf -- "$work_dir"
}
trap cleanup EXIT HUP INT TERM

test -r "$manifest" || { echo "ERROR: missing readable $manifest" >&2; exit 2; }
test -r "$public_manifest" || { echo "ERROR: missing readable $public_manifest" >&2; exit 2; }
test -x "$runner" || { echo "ERROR: missing executable $runner" >&2; exit 2; }

(cd "$repo_dir/src/curvatureCorrectionModels" && wmake libso)
"$runner" --campaign "$campaign_id" --repository "$repo_dir" \
    --manifest "$manifest" --output "$work_dir/validation-input.json"

validation_status=0
python3 "$repo_dir/tools/generate_curvature_validation_report.py" \
    "$work_dir/validation-input.json" --output-dir "$output_dir" \
    --strict-acceptance || validation_status=$?

cp "$public_manifest" "$output_dir/archive-manifest.json"
manifest_digest=$(sha256sum "$manifest")
manifest_digest=${manifest_digest%% *}
printf '%s  manifest.json\n' "$manifest_digest" \
    > "$output_dir/archive-manifest.sha256"

exit "$validation_status"
