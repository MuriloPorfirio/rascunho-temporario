#!/usr/bin/env bash
set -Eeuo pipefail

BASE="/raid/biomafia/murilo.aguiar"
REFERENCE_ROOT="${BASE}/GENOMAS-DE-REFERENCIA-E-KNOWN-SITES/grch37_b37"
CALL_ROOT="${BASE}/FASTQS-MELANOMA/teste-piloto-5-pacientes/CHAMADA-VARIANTES-PARABRICKS-HG19-TESTE-PILOTO"
REFERENCE="${REFERENCE_ROOT}/parabricks-bwa64/Homo_sapiens_assembly19.fasta"

GATK_IMAGE="broadinstitute/gatk:4.6.2.0"
BCFTOOLS_IMAGE="staphb/bcftools:1.19"
CPUS="8"
MEMORY="16g"
JAVA_MEMORY="12g"

PATIENTS=("ML_05" "ML_07" "ML_09" "ML_10")
declare -A TUMORS=(
    [ML_05]="ML_05_T_169823"
    [ML_07]="ML_07_T_919084"
    [ML_09]="ML_09_T_102617"
    [ML_10]="ML_10_T_367417"
)
declare -A NORMALS=(
    [ML_05]="ML_05_S_642203"
    [ML_07]="ML_07_S_535706"
    [ML_09]="ML_09_S_658074"
    [ML_10]="ML_10_S_295764"
)
declare -A RUNS=(
    [ML_05]="WGS-ML_05_20260821-144112"
    [ML_07]="WGS-ML_07_20260821-154526"
    [ML_09]="WGS-ML_09_20260821-170440"
    [ML_10]="WGS-ML_10_20260821-181042"
)

STAMP="$(date +'%Y%m%d-%H%M%S')"
CONTROL_DIR="${CALL_ROOT}/CONTROLE-FILTER-MUTECT-CALLS/FILTER-LOTE_${STAMP}"
BATCH_LOG="${CONTROL_DIR}/filtermutectcalls_lote.log"
SUMMARY="${CONTROL_DIR}/RESUMO-FILTERMUTECTCALLS.tsv"

section() {
    printf '\n============================================================\n%s\n============================================================\n' "$1"
}

die() {
    printf 'ERRO: %s\n' "$*" >&2
    exit 1
}

require_file() {
    [[ -s "$1" ]] || die "arquivo ausente ou vazio: $1"
    printf 'OK: %s\n' "$1"
}

run_bcftools_run() {
    local run_dir="$1"
    shift
    docker run --rm \
        --cpus=4 --memory=8g --memory-swap=8g --pids-limit=256 \
        -v "${run_dir}:/run:ro" \
        "$BCFTOOLS_IMAGE" \
        bcftools "$@"
}

run_bcftools_tmp() {
    local tmp_dir="$1"
    shift
    docker run --rm \
        --user "$(id -u):$(id -g)" \
        --cpus=4 --memory=8g --memory-swap=8g --pids-limit=256 \
        -v "${tmp_dir}:/work:rw" \
        "$BCFTOOLS_IMAGE" \
        bcftools "$@"
}

validate_completed_patient() {
    local out_dir="$1"
    local checksums="$2"
    (
        cd "$out_dir"
        sha256sum -c "$checksums" >/dev/null 2>&1
    )
}

mkdir -p "$CONTROL_DIR"
touch "$BATCH_LOG"
exec > >(tee -a "$BATCH_LOG") 2>&1

section "VERIFICAÇÃO DAS ENTRADAS DO LOTE"
require_file "$REFERENCE"
require_file "${REFERENCE}.fai"
require_file "${REFERENCE_ROOT}/parabricks-bwa64/Homo_sapiens_assembly19.dict"

docker image inspect "$GATK_IMAGE" >/dev/null 2>&1 || \
    die "imagem Docker não encontrada: ${GATK_IMAGE}"
docker image inspect "$BCFTOOLS_IMAGE" >/dev/null 2>&1 || \
    die "imagem Docker não encontrada: ${BCFTOOLS_IMAGE}"
printf 'OK: imagem %s\n' "$GATK_IMAGE"
printf 'OK: imagem %s\n' "$BCFTOOLS_IMAGE"

for patient in "${PATIENTS[@]}"; do
    tumor="${TUMORS[$patient]}"
    normal="${NORMALS[$patient]}"
    pair="${tumor}_vs_${normal}"
    run_dir="${CALL_ROOT}/AMOSTRAS/${patient}/${RUNS[$patient]}"

    [[ -d "$run_dir" ]] || die "corrida WGS não encontrada: $run_dir"
    require_file "${run_dir}/01-POSTPON/${pair}.mutectcaller.pon-annotated.vcf.gz"
    require_file "${run_dir}/01-POSTPON/${pair}.mutectcaller.pon-annotated.vcf.gz.tbi"
    require_file "${run_dir}/01-POSTPON/${pair}.mutectcaller.pon-annotated.vcf.gz.stats"
    require_file "${run_dir}/02-READ-ORIENTATION-MODEL/${pair}.artifact-priors.tar.gz"
    require_file "${run_dir}/03-CONTAMINATION/${tumor}.contamination.table"
    require_file "${run_dir}/03-CONTAMINATION/${tumor}.segments.table"
done

section "PLANO DO LOTE FILTERMUTECTCALLS"
printf 'PACIENTES:      ML_05, ML_07, ML_09 e ML_10\n'
printf 'ORDEM:          sequencial, um paciente por vez\n'
printf 'GATK:           4.6.2.0\n'
printf 'CPUs MÁXIMAS:   %s\n' "$CPUS"
printf 'RAM MÁXIMA:     %s\n' "$MEMORY"
printf 'GPUs:           não utilizadas\n'
printf 'CONTROLE:       %s\n' "$CONTROL_DIR"
printf 'LOG DO LOTE:    %s\n' "$BATCH_LOG"

printf 'PACIENTE\tSTATUS\tTOTAL\tPASS\tFILTRADAS\tPERCENTUAL_PASS\tTHRESHOLD\tFDR_ESTIMADO\tSENSIBILIDADE_ESTIMADA\tTEMPO_SEGUNDOS\tDIRETORIO\n' > "$SUMMARY"

process_patient() {
    local patient="$1"
    local tumor="${TUMORS[$patient]}"
    local normal="${NORMALS[$patient]}"
    local pair="${tumor}_vs_${normal}"
    local run_name="${RUNS[$patient]}"
    local run_dir="${CALL_ROOT}/AMOSTRAS/${patient}/${run_name}"

    local input_vcf="${run_dir}/01-POSTPON/${pair}.mutectcaller.pon-annotated.vcf.gz"
    local input_stats="${run_dir}/01-POSTPON/${pair}.mutectcaller.pon-annotated.vcf.gz.stats"
    local priors="${run_dir}/02-READ-ORIENTATION-MODEL/${pair}.artifact-priors.tar.gz"
    local contamination="${run_dir}/03-CONTAMINATION/${tumor}.contamination.table"
    local segments="${run_dir}/03-CONTAMINATION/${tumor}.segments.table"

    local out_dir="${run_dir}/04-FILTER-MUTECT-CALLS"
    local patient_stamp
    patient_stamp="$(date +'%Y%m%d-%H%M%S')"
    local tmp_dir="${BASE}/TMP_parabricks/MELANOMA-PILOTO-HG19-MUTECTCALLER/AMOSTRAS/${patient}/${run_name}/FILTERMUTECTCALLS_${patient_stamp}"

    local output_vcf="${pair}.mutectcaller.filtered.vcf.gz"
    local filtering_stats="${pair}.mutectcaller.filtering-stats.tsv"
    local distribution="${pair}.mutectcaller.filter-distribution.tsv"
    local manifest="MANIFESTO-FILTERMUTECTCALLS.txt"
    local patient_log="filtermutectcalls_${patient,,}.log"
    local checksums="SHA256SUMS.txt"

    section "PACIENTE: ${patient}"
    printf 'TUMOR:       %s\n' "$tumor"
    printf 'NORMAL:      %s\n' "$normal"
    printf 'CORRIDA:     %s\n' "$run_name"
    printf 'SAÍDA:       %s\n' "$out_dir"

    mkdir -p "$out_dir"

    local expected=(
        "$output_vcf"
        "${output_vcf}.tbi"
        "$filtering_stats"
        "$distribution"
        "$manifest"
        "$patient_log"
        "$checksums"
    )
    local existing=0
    local name
    for name in "${expected[@]}"; do
        [[ -e "${out_dir}/${name}" ]] && existing=$((existing + 1))
    done

    if [[ "$existing" -eq "${#expected[@]}" ]]; then
        validate_completed_patient "$out_dir" "$checksums" || \
            die "arquivos existentes de ${patient} falharam nos checksums"

        local total pass filtered pass_percent threshold fdr sensitivity
        total="$(run_bcftools_run "$run_dir" index --nrecords "/run/04-FILTER-MUTECT-CALLS/${output_vcf}")"
        pass="$(awk -F '\t' '$2=="PASS"{print $1; exit}' "${out_dir}/${distribution}")"
        pass="${pass:-0}"
        filtered="$((total - pass))"
        pass_percent="$(awk -v p="$pass" -v t="$total" 'BEGIN{printf "%.4f", 100*p/t}')"
        threshold="$(awk -F= '/^#<METADATA>threshold=/{gsub(/>$/, "", $2); print $2; exit}' "${out_dir}/${filtering_stats}")"
        fdr="$(awk -F= '/^#<METADATA>fdr=/{gsub(/>$/, "", $2); print $2; exit}' "${out_dir}/${filtering_stats}")"
        sensitivity="$(awk -F= '/^#<METADATA>sensitivity=/{gsub(/>$/, "", $2); print $2; exit}' "${out_dir}/${filtering_stats}")"

        printf 'STATUS: já concluído e validado; paciente ignorado.\n'
        printf '%s\tJA_CONCLUIDO\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t0\t%s\n' \
            "$patient" "$total" "$pass" "$filtered" "$pass_percent" \
            "$threshold" "$fdr" "$sensitivity" "$out_dir" >> "$SUMMARY"
        return 0
    fi

    [[ "$existing" -eq 0 ]] || \
        die "há saídas parciais em ${out_dir}; nenhuma foi substituída"

    mkdir -p "$tmp_dir"
    touch "${tmp_dir}/${patient_log}"
    printf 'TEMPORÁRIO:  %s\n' "$tmp_dir"
    printf 'LOG:         %s\n' "${out_dir}/${patient_log}"

    local input_inside="/run/01-POSTPON/$(basename "$input_vcf")"
    local input_records input_index_records input_samples input_variant_hash input_pon

    section "${patient} — VALIDAÇÃO DO VCF DE ENTRADA"
    input_records="$(run_bcftools_run "$run_dir" view -H "$input_inside" | wc -l)"
    input_index_records="$(run_bcftools_run "$run_dir" index --nrecords "$input_inside")"
    input_samples="$(run_bcftools_run "$run_dir" query -l "$input_inside" | paste -sd, -)"
    input_variant_hash="$(run_bcftools_run "$run_dir" query -f '%CHROM\t%POS\t%REF\t%ALT\n' "$input_inside" | sha256sum | awk '{print $1}')"
    input_pon="$(run_bcftools_run "$run_dir" query -f '%INFO/PON\n' "$input_inside" | awk '$1==1{n++} END{print n+0}')"

    [[ "$input_records" -gt 0 ]] || die "VCF de entrada sem registros em ${patient}"
    [[ "$input_records" -eq "$input_index_records" ]] || die "índice de entrada inconsistente em ${patient}"
    [[ "$input_samples" == "${normal},${tumor}" ]] || die "amostras inesperadas em ${patient}: ${input_samples}"
    [[ "$input_pon" -gt 0 ]] || die "VCF de entrada sem flags INFO/PON em ${patient}"

    local contamination_value segments_count
    contamination_value="$(awk '!/^#/ && $1!="sample" && NF{print $2; exit}' "$contamination")"
    segments_count="$(awk 'BEGIN{n=0} !/^#/ && $1!="contig" && NF{n++} END{print n}' "$segments")"

    printf 'REGISTROS:             %s\n' "$input_records"
    printf 'REGISTROS PELO ÍNDICE: %s\n' "$input_index_records"
    printf 'AMOSTRAS:              %s\n' "$input_samples"
    printf 'VARIANTES INFO/PON:    %s\n' "$input_pon"
    printf 'CONTAMINAÇÃO:          %s\n' "$contamination_value"
    printf 'SEGMENTOS:             %s\n' "$segments_count"
    printf 'HASH CHROM/POS/REF/ALT: %s\n' "$input_variant_hash"

    section "${patient} — EXECUTANDO FILTERMUTECTCALLS"
    local start_epoch end_epoch elapsed
    start_epoch="$(date +%s)"

    docker run --rm \
        --user "$(id -u):$(id -g)" \
        --cpus="$CPUS" \
        --memory="$MEMORY" \
        --memory-swap="$MEMORY" \
        --pids-limit=512 \
        --tmpfs /tmp:rw,exec,size=2g \
        -v "${REFERENCE_ROOT}:/reference:ro" \
        -v "${run_dir}:/run:ro" \
        -v "${tmp_dir}:/work:rw" \
        --entrypoint /gatk/gatk \
        "$GATK_IMAGE" \
        --java-options "-Xms1g -Xmx${JAVA_MEMORY} -Djava.io.tmpdir=/tmp" \
        FilterMutectCalls \
        -R /reference/parabricks-bwa64/Homo_sapiens_assembly19.fasta \
        -V "$input_inside" \
        --stats "/run/01-POSTPON/$(basename "$input_stats")" \
        --orientation-bias-artifact-priors "/run/02-READ-ORIENTATION-MODEL/$(basename "$priors")" \
        --contamination-table "/run/03-CONTAMINATION/$(basename "$contamination")" \
        --tumor-segmentation "/run/03-CONTAMINATION/$(basename "$segments")" \
        --filtering-stats "/work/${filtering_stats}" \
        -O "/work/${output_vcf}" \
        2>&1 | tee -a "${tmp_dir}/${patient_log}" || \
        die "FilterMutectCalls falhou em ${patient}; temporário preservado: ${tmp_dir}"

    end_epoch="$(date +%s)"
    elapsed="$((end_epoch - start_epoch))"
    printf 'TEMPO DO PACIENTE: %s segundos\n' "$elapsed"

    section "${patient} — RECRIAÇÃO PADRONIZADA DO ÍNDICE TBI"
    require_file "${tmp_dir}/${output_vcf}"
    require_file "${tmp_dir}/${filtering_stats}"
    run_bcftools_tmp "$tmp_dir" index --tbi --force "/work/${output_vcf}"
    require_file "${tmp_dir}/${output_vcf}.tbi"

    section "${patient} — VALIDAÇÃO DO RESULTADO"
    local output_records output_index_records output_samples output_variant_hash output_pon
    local contigs indexed_records pass filtered pass_percent filter_headers
    output_records="$(run_bcftools_tmp "$tmp_dir" view -H "/work/${output_vcf}" | wc -l)"
    output_index_records="$(run_bcftools_tmp "$tmp_dir" index --nrecords "/work/${output_vcf}")"
    output_samples="$(run_bcftools_tmp "$tmp_dir" query -l "/work/${output_vcf}" | paste -sd, -)"
    output_variant_hash="$(run_bcftools_tmp "$tmp_dir" query -f '%CHROM\t%POS\t%REF\t%ALT\n' "/work/${output_vcf}" | sha256sum | awk '{print $1}')"
    output_pon="$(run_bcftools_tmp "$tmp_dir" query -f '%INFO/PON\n' "/work/${output_vcf}" | awk '$1==1{n++} END{print n+0}')"
    contigs="$(run_bcftools_tmp "$tmp_dir" index -s "/work/${output_vcf}" | awk '$3>0{print $1}' | paste -sd, -)"
    [[ -n "$contigs" ]] || die "novo índice de ${patient} sem contigs com registros"
    indexed_records="$(run_bcftools_tmp "$tmp_dir" view -H -r "$contigs" "/work/${output_vcf}" | wc -l)"

    [[ "$output_records" -eq "$input_records" ]] || die "quantidade de variantes modificada em ${patient}"
    [[ "$output_index_records" -eq "$output_records" ]] || die "metadado do novo TBI inconsistente em ${patient}"
    [[ "$indexed_records" -eq "$output_records" ]] || die "consulta indexada omitiu registros em ${patient}"
    [[ "$output_samples" == "$input_samples" ]] || die "amostras modificadas em ${patient}"
    [[ "$output_variant_hash" == "$input_variant_hash" ]] || die "CHROM/POS/REF/ALT modificados em ${patient}"
    [[ "$output_pon" -eq "$input_pon" ]] || die "flags INFO/PON modificadas em ${patient}"

    pass="$(run_bcftools_tmp "$tmp_dir" query -f '%FILTER\n' "/work/${output_vcf}" | awk '$1=="PASS"{n++} END{print n+0}')"
    filtered="$((output_records - pass))"
    pass_percent="$(awk -v p="$pass" -v t="$output_records" 'BEGIN{printf "%.4f", 100*p/t}')"
    filter_headers="$(run_bcftools_tmp "$tmp_dir" view -h "/work/${output_vcf}" | grep -c '^##FILTER=<ID=' || true)"
    [[ "$filter_headers" -gt 0 ]] || die "VCF de ${patient} sem definições FILTER"

    {
        printf 'COUNT\tFILTER\n'
        run_bcftools_tmp "$tmp_dir" query -f '%FILTER\n' "/work/${output_vcf}" |
            LC_ALL=C sort |
            uniq -c |
            awk '{$1=$1; count=$1; $1=""; sub(/^ /,""); print count "\t" $0}' |
            sort -nr
    } > "${tmp_dir}/${distribution}"

    local threshold fdr sensitivity
    threshold="$(awk -F= '/^#<METADATA>threshold=/{gsub(/>$/, "", $2); print $2; exit}' "${tmp_dir}/${filtering_stats}")"
    fdr="$(awk -F= '/^#<METADATA>fdr=/{gsub(/>$/, "", $2); print $2; exit}' "${tmp_dir}/${filtering_stats}")"
    sensitivity="$(awk -F= '/^#<METADATA>sensitivity=/{gsub(/>$/, "", $2); print $2; exit}' "${tmp_dir}/${filtering_stats}")"
    [[ -n "$threshold" && -n "$fdr" && -n "$sensitivity" ]] || \
        die "estatísticas globais ausentes em ${patient}"

    printf 'REGISTROS ANTES:        %s\n' "$input_records"
    printf 'REGISTROS DEPOIS:       %s\n' "$output_records"
    printf 'METADADO DO ÍNDICE:     %s\n' "$output_index_records"
    printf 'RECUPERADOS PELO TBI:   %s\n' "$indexed_records"
    printf 'VARIANTES PASS:         %s (%s%%)\n' "$pass" "$pass_percent"
    printf 'VARIANTES FILTRADAS:    %s\n' "$filtered"
    printf 'THRESHOLD:              %s\n' "$threshold"
    printf 'FDR ESTIMADO:           %s\n' "$fdr"
    printf 'SENSIBILIDADE ESTIMADA: %s\n' "$sensitivity"
    printf 'INFO/PON PRESERVADAS:   %s\n' "$output_pon"
    printf 'REGISTROS/ALELOS/AMOSTRAS/TBI: OK\n'

    cat > "${tmp_dir}/${manifest}" <<EOF
ETAPA: FilterMutectCalls
PACIENTE: ${patient}
TUMOR: ${tumor}
NORMAL: ${normal}
GATK: 4.6.2.0
IMAGEM_GATK: ${GATK_IMAGE}
IMAGEM_BCFTOOLS: ${BCFTOOLS_IMAGE}
REFERENCIA: ${REFERENCE}
VCF_ENTRADA: ${input_vcf}
STATS: ${input_stats}
ORIENTATION_PRIORS: ${priors}
CONTAMINATION_TABLE: ${contamination}
TUMOR_SEGMENTATION: ${segments}
CONTAMINATION_VALUE: ${contamination_value}
REGISTROS_ENTRADA: ${input_records}
REGISTROS_SAIDA: ${output_records}
VARIANTES_PASS: ${pass}
VARIANTES_FILTRADAS: ${filtered}
PERCENTUAL_PASS: ${pass_percent}
THRESHOLD: ${threshold}
FDR_ESTIMADO: ${fdr}
SENSIBILIDADE_ESTIMADA: ${sensitivity}
VARIANTES_INFO_PON: ${output_pon}
HASH_CHROM_POS_REF_ALT: ${output_variant_hash}
CPUS_MAXIMAS: ${CPUS}
RAM_MAXIMA: ${MEMORY}
GPUS: nao utilizadas
TEMPO_SEGUNDOS: ${elapsed}
EOF

    mv "${tmp_dir}/${output_vcf}" "${out_dir}/${output_vcf}"
    mv "${tmp_dir}/${output_vcf}.tbi" "${out_dir}/${output_vcf}.tbi"
    mv "${tmp_dir}/${filtering_stats}" "${out_dir}/${filtering_stats}"
    mv "${tmp_dir}/${distribution}" "${out_dir}/${distribution}"
    mv "${tmp_dir}/${manifest}" "${out_dir}/${manifest}"
    mv "${tmp_dir}/${patient_log}" "${out_dir}/${patient_log}"

    (
        cd "$out_dir"
        sha256sum \
            "$output_vcf" \
            "${output_vcf}.tbi" \
            "$filtering_stats" \
            "$distribution" \
            "$manifest" \
            "$patient_log" \
            > "$checksums"
        sha256sum -c "$checksums"
    ) || die "falha na validação final de ${patient}"

    rmdir "$tmp_dir" 2>/dev/null || true

    printf '%s\tCONCLUIDO\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$patient" "$output_records" "$pass" "$filtered" "$pass_percent" \
        "$threshold" "$fdr" "$sensitivity" "$elapsed" "$out_dir" >> "$SUMMARY"

    section "${patient} CONCLUÍDO E VALIDADO"
    printf 'VARIANTES TOTAIS:    %s\n' "$output_records"
    printf 'VARIANTES PASS:      %s (%s%%)\n' "$pass" "$pass_percent"
    printf 'VARIANTES FILTRADAS: %s\n' "$filtered"
    printf 'TEMPO:               %s segundos\n' "$elapsed"
    printf 'SAÍDA:               %s\n' "$out_dir"
}

batch_start="$(date +%s)"
completed=0

for patient in "${PATIENTS[@]}"; do
    process_patient "$patient"
    completed=$((completed + 1))
done

batch_end="$(date +%s)"
batch_elapsed="$((batch_end - batch_start))"

section "LOTE FILTERMUTECTCALLS CONCLUÍDO"
printf 'PACIENTES PROCESSADOS OU JÁ VALIDADOS: %s\n' "$completed"
printf 'PACIENTES ESPERADOS:                  %s\n' "${#PATIENTS[@]}"
printf 'TEMPO TOTAL DO LOTE:                  %s segundos\n' "$batch_elapsed"
printf 'RESULTADO:                            TODOS CONCLUÍDOS E VALIDADOS\n'
printf 'RESUMO:                               %s\n' "$SUMMARY"
printf 'LOG:                                  %s\n' "$BATCH_LOG"
