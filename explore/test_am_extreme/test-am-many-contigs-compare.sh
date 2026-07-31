#!/usr/bin/bash
set -ueo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
NLINES="${NLINES:-400000}"
NJOBS="${NJOBS:-32}"
FCOV="${FCOV:-20}"
BATCH_SIZE="${BATCH_SIZE:-8192}"
MEM_LIMIT_KIB="${MEM_LIMIT_KIB:-$((10 * 1024 * 1024))}"
POLL_SECONDS="${POLL_SECONDS:-1}"
GENERATED="${HERE}/generated"
FASTA="${GENERATED}/many_contigs_${NLINES}.fa"

mkdir -p "${GENERATED}"
if [ ! -f "${FASTA}" ]; then
    "${HERE}/c/bin/generate_many_contigs" 2>/dev/null | head -n "${NLINES}" >"${FASTA}" || true
fi
printf 'input: %s (%s contigs), parallel=%s batch=%s mem_limit=%s MiB\n' \
    "$(basename "${FASTA}")" "$(grep -c '^>' "${FASTA}")" "${NJOBS}" "${BATCH_SIZE}" "$((MEM_LIMIT_KIB / 1024))"

watch_memory() {
    local pid="$1" peak_file="$2" kill_file="$3" peak=0 rss
    while kill -0 "${pid}" 2>/dev/null; do
        rss="$(awk '/^VmRSS:/{print $2}' "/proc/${pid}/status" 2>/dev/null || true)"
        rss="${rss:-0}"
        if [ "${rss}" -gt "${peak}" ]; then
            peak="${rss}"
        fi
        if [ "${rss}" -gt "${MEM_LIMIT_KIB}" ]; then
            echo "${peak}" >"${peak_file}"
            echo 1 >"${kill_file}"
            kill -9 "${pid}" 2>/dev/null || true
            return 0
        fi
        sleep "${POLL_SECONDS}"
    done
    echo "${peak}" >"${peak_file}"
}

run_backend() {
    local name="$1"
    local exe="${HERE}/opt/am_${name}/bin/art_modern"
    local log="${GENERATED}/compare_${name}.log"
    local peak_file="${GENERATED}/compare_${name}.peak"
    local kill_file="${GENERATED}/compare_${name}.killed"
    local rc=0 pid watchdog peak killed

    rm -f "${peak_file}" "${kill_file}"
    "${exe}" \
        --i-file "${FASTA}" \
        --i-type fasta \
        --mode template \
        --lc se \
        --i-parser stream \
        --i-fcov "${FCOV}" \
        --i-batch_size "${BATCH_SIZE}" \
        --parallel "${NJOBS}" \
        --o-hl_sam /dev/null \
        --reporting_interval-job_executor 10 \
        --reporting_interval-job_pool 50 \
        >"${log}" 2>&1 &
    pid=$!
    watch_memory "${pid}" "${peak_file}" "${kill_file}" &
    watchdog=$!
    wait "${pid}" || rc=$?
    wait "${watchdog}" 2>/dev/null || true

    peak="$(cat "${peak_file}" 2>/dev/null || echo 0)"
    killed="$(cat "${kill_file}" 2>/dev/null || echo 0)"
    if [ "${killed}" = "1" ]; then
        printf '%-5s OOM-KILLED at %s MiB\n' "${name}" "$((peak / 1024))"
    elif [ "${rc}" -gt 128 ]; then
        printf '%-5s exit=%-4s SIG%-8s peak=%s MiB\n' "${name}" "${rc}" "$(kill -l "$((rc - 128))")" "$((peak / 1024))"
    else
        printf '%-5s exit=%-4s %-12s peak=%s MiB\n' "${name}" "${rc}" "" "$((peak / 1024))"
    fi
    grep -iE "fatal|not 0 when exiting" "${log}" | tail -3 || true
    return "${rc}"
}

status=0
run_backend bs || status=1
run_backend asio || status=1
exit "${status}"
