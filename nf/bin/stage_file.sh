#!/bin/bash
# stage_file.sh - Shared file staging utility for OSS and local paths
# Usage: source this file in your Nextflow process script
#   source ${projectDir}/bin/stage_file.sh
#   R1_LOCAL=$(stage_file "${exp_r1}")

stage_file() {
    local src="$1"
    local name=$(basename "$src")
    if [[ "$src" == oss://* ]]; then
        echo "Downloading $src" >&2
        local script_dir
        script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
        cfg="${STAGE_OSSUTILCONFIG:-${script_dir}/.ossutilconfig}"
        if [ -f "$cfg" ]; then
            ossutil cp --sign-version v4 --region cn-beijing -c "$cfg" "$src" "$name" 1>&2
        else
            ossutil cp --sign-version v4 --region cn-beijing -e oss-cn-beijing-internal.aliyuncs.com -i "$AccessKeyId" -k "$AccessKeySecret" "$src" "$name" 1>&2
        fi
        if [ $? -ne 0 ]; then
            echo "ERROR: ossutil failed for $src" >&2
            exit 1
        fi
    else
        if [ ! -e "$src" ]; then
            echo "ERROR: Source file not found: $src" >&2
            exit 1
        fi
        if [ -f "$src" ]; then
            ln -sf "$src" "$name"
        else
            cp "$src" "$name" || { echo "ERROR: Failed to copy $src" >&2; exit 1; }
        fi
    fi
    echo "$name"
}
