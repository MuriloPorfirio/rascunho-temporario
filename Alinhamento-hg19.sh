#!/bin/bash

set -Eeuo pipefail
umask 027

COMPOSE="/raid/biomafia/murilo.aguiar/dockers_compose/docker-compose.yml"

TRIM_DIR="/raid/biomafia/murilo.aguiar/FASTQS-MELANOMA/teste-piloto-5-pacientes/FASTQS-TRIMADOS-TESTE-PILOTO"
ALIGN_ROOT="/raid/biomafia/murilo.aguiar/FASTQS-MELANOMA/teste-piloto-5-pacientes/ALINHAMENTOS-PARABRICKS-HG19-TESTE-PILOTO"
RESULTS_ROOT="${ALIGN_ROOT}/AMOSTRAS"
MANIFEST="${ALIGN_ROOT}/manifesto_alinhamento_hg19.csv"
LOG_DIR="${ALIGN_ROOT}/logs"
LOCK="${ALIGN_ROOT}/alinhamento_parabricks_hg19.lock"

TMP_ROOT="/raid/biomafia/murilo.aguiar/TMP_parabricks/MELANOMA-PILOTO-HG19/AMOSTRAS"

REFERENCE_ROOT="/raid/biomafia/murilo.aguiar/GENOMAS-DE-REFERENCIA-E-KNOWN-SITES/grch37_b37"
REFERENCE_DIR="${REFERENCE_ROOT}/parabricks-bwa64"
KNOWN_DIR="${REFERENCE_ROOT}/known-sites"
REFERENCE_HOST="${REFERENCE_DIR}/Homo_sapiens_assembly19.fasta"

MILLS_HOST="${KNOWN_DIR}/Mills_and_1000G_gold_standard.indels.b37.vcf.gz"
INDELS_HOST="${KNOWN_DIR}/1000G_phase1.indels.b37.vcf.gz"
SNPS_HOST="${KNOWN_DIR}/1000G_phase1.snps.high_confidence.b37.vcf.gz"

RUN_ID="$(date '+%Y%m%d_%H%M%S')"
RUN_LOG="${LOG_DIR}/alinhamento_parabricks_hg19_${RUN_ID}.log"

GPUS=4
THREADS_PER_GPU=16
TOTAL_CPU_THREADS=64
BWA_STREAMS=4
MEMORY_LIMIT_GB=192

mkdir -p \
    "$RESULTS_ROOT" \
    "$LOG_DIR" \
    "$TMP_ROOT"

chmod 750 \
    "$RESULTS_ROOT" \
    "$LOG_DIR" \
    "$TMP_ROOT"

exec > >(tee -a "$RUN_LOG") 2>&1

exec 9>"$LOCK"

if ! flock -n 9; then
    echo "ERRO: já existe outro alinhamento deste projeto em andamento."
    exit 1
fi

for CMD in docker flock nvidia-smi awk sed grep find tee wc chmod df hostname date; do
    if ! command -v "$CMD" >/dev/null 2>&1; then
        echo "ERRO: comando obrigatório não encontrado: $CMD"
        exit 1
    fi
done

docker compose \
    -f "$COMPOSE" \
    config --quiet

REQUIRED_REFERENCE_FILES=(
    "$REFERENCE_HOST"
    "${REFERENCE_HOST}.fai"
    "${REFERENCE_HOST}.amb"
    "${REFERENCE_HOST}.ann"
    "${REFERENCE_HOST}.bwt"
    "${REFERENCE_HOST}.pac"
    "${REFERENCE_HOST}.sa"
    "${REFERENCE_DIR}/Homo_sapiens_assembly19.dict"
    "$MILLS_HOST"
    "${MILLS_HOST}.tbi"
    "$INDELS_HOST"
    "${INDELS_HOST}.tbi"
    "$SNPS_HOST"
    "${SNPS_HOST}.tbi"
)

for FILE in "${REQUIRED_REFERENCE_FILES[@]}"; do
    if [ ! -s "$FILE" ]; then
        echo "ERRO: referência ou known-site ausente ou vazio:"
        echo "$FILE"
        exit 1
    fi
done

if [ ! -s "$MANIFEST" ]; then
    echo "ERRO: manifesto ausente ou vazio:"
    echo "$MANIFEST"
    exit 1
fi

if grep -q $'\r' "$MANIFEST"; then
    echo "ERRO: o manifesto possui finais de linha do Windows (CRLF)."
    exit 1
fi

EXPECTED_HEADER="sample,patient,type,r1,r2,rg_id,rg_lb,rg_pl,rg_pu"
ACTUAL_HEADER="$(sed -n '1p' "$MANIFEST")"

if [ "$ACTUAL_HEADER" != "$EXPECTED_HEADER" ]; then
    echo "ERRO: cabeçalho inesperado no manifesto."
    echo "Esperado: $EXPECTED_HEADER"
    echo "Obtido:  $ACTUAL_HEADER"
    exit 1
fi

mapfile -t MANIFEST_LINES < <(sed '1d' "$MANIFEST")

if [ "${#MANIFEST_LINES[@]}" -ne 10 ]; then
    echo "ERRO: o manifesto possui ${#MANIFEST_LINES[@]} amostras; eram esperadas 10."
    exit 1
fi

if ! awk -F ',' '
    NR > 1 && NF != 9 {
        printf "ERRO: linha %d possui %d campos; eram esperados 9.\n", NR, NF
        found_error=1
    }
    END {
        exit found_error ? 1 : 0
    }
' "$MANIFEST"; then
    exit 1
fi

declare -A SEEN_SAMPLES=()
declare -A SEEN_RG_IDS=()

for LINE in "${MANIFEST_LINES[@]}"; do
    IFS=',' read -r \
        SAMPLE PATIENT TYPE R1 R2 RG_ID RG_LB RG_PL RG_PU \
        <<< "$LINE"

    if [ -z "$SAMPLE" ] ||
       [ -z "$PATIENT" ] ||
       [ -z "$TYPE" ] ||
       [ -z "$R1" ] ||
       [ -z "$R2" ] ||
       [ -z "$RG_ID" ] ||
       [ -z "$RG_LB" ] ||
       [ -z "$RG_PL" ] ||
       [ -z "$RG_PU" ]; then
        echo "ERRO: campo vazio no manifesto:"
        echo "$LINE"
        exit 1
    fi

    if [ -n "${SEEN_SAMPLES[$SAMPLE]:-}" ]; then
        echo "ERRO: amostra duplicada no manifesto: $SAMPLE"
        exit 1
    fi

    if [ -n "${SEEN_RG_IDS[$RG_ID]:-}" ]; then
        echo "ERRO: RG ID duplicado no manifesto: $RG_ID"
        exit 1
    fi

    SEEN_SAMPLES[$SAMPLE]=1
    SEEN_RG_IDS[$RG_ID]=1

    if [ "$R1" != "${SAMPLE}_R1_val_1.fq.gz" ] ||
       [ "$R2" != "${SAMPLE}_R2_val_2.fq.gz" ]; then
        echo "ERRO: nomes de FASTQ incompatíveis com a amostra $SAMPLE."
        exit 1
    fi

    if [ "$RG_LB" != "$SAMPLE" ] || [ "$RG_PL" != "ILLUMINA" ]; then
        echo "ERRO: LB ou PL inválido para $SAMPLE."
        exit 1
    fi

    case "$TYPE" in
        NORMAL)
            if [[ "$SAMPLE" != "${PATIENT}_S_"* ]]; then
                echo "ERRO: amostra NORMAL incompatível: $SAMPLE"
                exit 1
            fi
            ;;
        TUMOR)
            if [[ "$SAMPLE" != "${PATIENT}_T_"* ]]; then
                echo "ERRO: amostra TUMOR incompatível: $SAMPLE"
                exit 1
            fi
            ;;
        *)
            echo "ERRO: tipo inválido para $SAMPLE: $TYPE"
            exit 1
            ;;
    esac

    for FASTQ in "$TRIM_DIR/$R1" "$TRIM_DIR/$R2"; do
        if [ ! -s "$FASTQ" ]; then
            echo "ERRO: FASTQ ausente ou vazio:"
            echo "$FASTQ"
            exit 1
        fi
    done
done

check_selected_gpus() {
    local GPU_STATUS
    local GPU_INDEX
    local GPU_MEMORY
    local GPU_UTIL
    local GPU_ROWS=0
    local BUSY=0

    if ! GPU_STATUS="$(
        nvidia-smi \
            -i 4,5,6,7 \
            --query-gpu=index,memory.used,utilization.gpu \
            --format=csv,noheader,nounits
    )"; then
        echo "ERRO: não foi possível consultar as GPUs 4,5,6,7."
        return 1
    fi

    while IFS=',' read -r GPU_INDEX GPU_MEMORY GPU_UTIL; do
        GPU_INDEX="${GPU_INDEX//[[:space:]]/}"
        GPU_MEMORY="${GPU_MEMORY//[[:space:]]/}"
        GPU_UTIL="${GPU_UTIL//[[:space:]]/}"

        GPU_ROWS=$((GPU_ROWS + 1))

        printf 'GPU %s: %s MiB; utilização %s%%\n' \
            "$GPU_INDEX" \
            "$GPU_MEMORY" \
            "$GPU_UTIL"

        if [ "$GPU_MEMORY" -gt 100 ] || [ "$GPU_UTIL" -gt 0 ]; then
            BUSY=1
        fi
    done <<< "$GPU_STATUS"

    if [ "$GPU_ROWS" -ne 4 ]; then
        echo "ERRO: eram esperadas informações de quatro GPUs."
        return 1
    fi

    if [ "$BUSY" -ne 0 ]; then
        echo "ERRO: pelo menos uma das GPUs 4,5,6,7 parece estar ocupada."
        return 1
    fi
}

echo
echo "=========================================================="
echo "ALINHAMENTO PARABRICKS FQ2BAM — HG19/B37"
echo "Data:                 $(date '+%Y-%m-%d %H:%M:%S')"
echo "Host:                 $(hostname)"
echo "Amostras:             10"
echo "Execução:             estritamente sequencial"
echo "GPUs:                 4,5,6,7"
echo "Threads por GPU:      $THREADS_PER_GPU"
echo "CPUs máximas:         $TOTAL_CPU_THREADS"
echo "Memória Docker:       256 GiB"
echo "Memória sort/postsort: ${MEMORY_LIMIT_GB} GB"
echo "Streams por GPU:      $BWA_STREAMS"
echo "Referência:           Homo_sapiens_assembly19.fasta"
echo "Known-sites:          3"
echo "Resultados:           $RESULTS_ROOT"
echo "Temporários:          $TMP_ROOT"
echo "Manifesto:            $MANIFEST"
echo "Log geral:            $RUN_LOG"
echo "=========================================================="
echo

df -h "$ALIGN_ROOT"

echo
echo "Estado inicial das GPUs selecionadas:"
check_selected_gpus

PROCESSADOS=0
PULADOS=0
ATUAL=0

for LINE in "${MANIFEST_LINES[@]}"; do
    IFS=',' read -r \
        SAMPLE PATIENT TYPE R1 R2 RG_ID RG_LB RG_PL RG_PU \
        <<< "$LINE"

    ATUAL=$((ATUAL + 1))

    SAMPLE_DIR="${RESULTS_ROOT}/${SAMPLE}"
    QC_DIR="${SAMPLE_DIR}/qc_metrics"
    TMP_SAMPLE="${TMP_ROOT}/${SAMPLE}"

    BAM_HOST="${SAMPLE_DIR}/${SAMPLE}.bam"
    BAI_HOST="${BAM_HOST}.bai"
    CHRS_HOST="${SAMPLE_DIR}/${SAMPLE}_chrs.txt"
    RECAL_HOST="${SAMPLE_DIR}/${SAMPLE}.recal.table"
    DUP_HOST="${SAMPLE_DIR}/${SAMPLE}.dup_metrics.txt"
    SAMPLE_LOG="${SAMPLE_DIR}/${SAMPLE}.fq2bam.log"
    VALIDATED_MARKER="${SAMPLE_DIR}/${SAMPLE}.fq2bam.validated.ok"

    echo
    echo "----------------------------------------------------------"
    echo "[$ATUAL/10] $SAMPLE"
    echo "Paciente: $PATIENT"
    echo "Tipo:     $TYPE"
    echo "R1:       $R1"
    echo "R2:       $R2"
    echo "RG ID:    $RG_ID"
    echo "Início:   $(date '+%Y-%m-%d %H:%M:%S')"

    if [ -s "$VALIDATED_MARKER" ]; then
        EXISTING_OUTPUTS=(
            "$BAM_HOST"
            "$BAI_HOST"
            "$CHRS_HOST"
            "$RECAL_HOST"
            "$DUP_HOST"
            "$SAMPLE_LOG"
        )

        for FILE in "${EXISTING_OUTPUTS[@]}"; do
            if [ ! -s "$FILE" ]; then
                echo "ERRO: marcador de validação existe, mas falta:"
                echo "$FILE"
                exit 1
            fi
        done

        EXISTING_QC_COUNT="$(
            find "$QC_DIR" \
                -maxdepth 1 \
                -type f \
                -size +0c |
                wc -l
        )"

        if [ "$EXISTING_QC_COUNT" -lt 20 ]; then
            echo "ERRO: marcador de validação existe, mas o QC está incompleto."
            exit 1
        fi

        echo "Já estava concluída e validada. Pulando."
        PULADOS=$((PULADOS + 1))
        continue
    fi

    if [ -e "$SAMPLE_DIR" ]; then
        echo "ERRO: existe diretório parcial sem marcador de validação:"
        echo "$SAMPLE_DIR"
        echo "Nada será apagado ou sobrescrito."
        exit 1
    fi

    if [ -e "$TMP_SAMPLE" ]; then
        echo "ERRO: existe diretório temporário anterior:"
        echo "$TMP_SAMPLE"
        echo "Nada será apagado ou sobrescrito."
        exit 1
    fi

    echo
    echo "Verificando as GPUs antes de iniciar $SAMPLE:"
    check_selected_gpus

    mkdir -p \
        "$SAMPLE_DIR" \
        "$QC_DIR" \
        "$TMP_SAMPLE"

    chmod 750 \
        "$SAMPLE_DIR" \
        "$QC_DIR" \
        "$TMP_SAMPLE"

    READ_GROUP="@RG\\tID:${RG_ID}\\tLB:${RG_LB}\\tPL:${RG_PL}\\tSM:${SAMPLE}\\tPU:${RG_PU}"

    echo
    echo "READ GROUP:"
    printf '%s\n' "$READ_GROUP"
    echo

    set +e

    docker compose \
        -f "$COMPOSE" \
        run \
        --rm \
        --no-deps \
        --name "murilo-melanoma-fq2bam-${SAMPLE}" \
        parabricks \
        fq2bam \
        --ref /reference/parabricks-bwa64/Homo_sapiens_assembly19.fasta \
        --in-fq \
            "/input/$R1" \
            "/input/$R2" \
            "$READ_GROUP" \
        --knownSites /reference/known-sites/Mills_and_1000G_gold_standard.indels.b37.vcf.gz \
        --knownSites /reference/known-sites/1000G_phase1.indels.b37.vcf.gz \
        --knownSites /reference/known-sites/1000G_phase1.snps.high_confidence.b37.vcf.gz \
        --out-bam "/output/AMOSTRAS/${SAMPLE}/${SAMPLE}.bam" \
        --out-recal-file "/output/AMOSTRAS/${SAMPLE}/${SAMPLE}.recal.table" \
        --out-duplicate-metrics "/output/AMOSTRAS/${SAMPLE}/${SAMPLE}.dup_metrics.txt" \
        --out-qc-metrics-dir "/output/AMOSTRAS/${SAMPLE}/qc_metrics" \
        --tmp-dir "/parabricks-tmp/AMOSTRAS/${SAMPLE}" \
        --num-gpus "$GPUS" \
        --num-cpu-threads-per-stage "$THREADS_PER_GPU" \
        --bwa-nstreams "$BWA_STREAMS" \
        --bwa-options="-K 10000000" \
        --memory-limit "$MEMORY_LIMIT_GB" \
        --gpusort \
        --gpuwrite \
        --monitor-usage \
        --x3 \
        2>&1 |
    tee "$SAMPLE_LOG"

    STATUS_PARABRICKS="${PIPESTATUS[0]}"

    set -e

    echo
    printf 'CÓDIGO DE SAÍDA DO PARABRICKS: %s\n' "$STATUS_PARABRICKS"

    if [ "$STATUS_PARABRICKS" -ne 0 ]; then
        echo "ERRO: o fq2bam falhou para $SAMPLE."
        echo "Resultados parciais e temporários foram preservados."
        exit "$STATUS_PARABRICKS"
    fi

    EXPECTED_OUTPUTS=(
        "$BAM_HOST"
        "$BAI_HOST"
        "$CHRS_HOST"
        "$RECAL_HOST"
        "$DUP_HOST"
        "$SAMPLE_LOG"
    )

    for FILE in "${EXPECTED_OUTPUTS[@]}"; do
        if [ ! -s "$FILE" ]; then
            echo "ERRO: resultado esperado ausente ou vazio:"
            echo "$FILE"
            exit 1
        fi
    done

    QC_COUNT="$(
        find "$QC_DIR" \
            -maxdepth 1 \
            -type f \
            -size +0c |
            wc -l
    )"

    if [ "$QC_COUNT" -lt 20 ]; then
        echo "ERRO: somente $QC_COUNT arquivos de QC não vazios foram encontrados."
        exit 1
    fi

    echo
    echo "Validando BAM, BAI, cabeçalho, BQSR e métricas..."

    docker run \
        --rm \
        --network lakehouse \
        --cpus=2 \
        --memory=4g \
        --memory-swap=4g \
        --pids-limit=128 \
        --read-only \
        --tmpfs /tmp:size=256m \
        --cap-drop ALL \
        --security-opt no-new-privileges:true \
        --user 1006:1002 \
        --volume "${SAMPLE_DIR}:/data:ro" \
        --entrypoint sh \
        staphb/samtools:1.19 \
        -lc '
set -eu

SAMPLE="$1"
RG_ID="$2"
RG_LB="$3"
RG_PL="$4"
RG_PU="$5"

BAM="/data/${SAMPLE}.bam"
BAI="${BAM}.bai"
RECAL="/data/${SAMPLE}.recal.table"
DUP="/data/${SAMPLE}.dup_metrics.txt"

for FILE in "$BAM" "$BAI" "$RECAL" "$DUP"; do
    if [ ! -s "$FILE" ]; then
        echo "ERRO: arquivo ausente ou vazio dentro do validador: $FILE"
        exit 1
    fi
done

samtools quickcheck -v "$BAM"
samtools idxstats "$BAM" > /tmp/idxstats.txt

IDX_LINES="$(wc -l < /tmp/idxstats.txt)"

if [ "$IDX_LINES" -le 1 ]; then
    echo "ERRO: o índice BAI não produziu estatísticas."
    exit 1
fi

HD_LINE="$(samtools view -H "$BAM" | grep -m1 "^@HD")"
SQ_LINE="$(samtools view -H "$BAM" | grep -m1 "^@SQ")"
RG_LINE="$(samtools view -H "$BAM" | grep -m1 "^@RG")"

printf "%s\n" "$HD_LINE" | grep -F "SO:coordinate" >/dev/null
printf "%s\n" "$SQ_LINE" | grep -F "SN:1" >/dev/null
printf "%s\n" "$SQ_LINE" | grep -F "LN:249250621" >/dev/null

for TAG in \
    "ID:${RG_ID}" \
    "LB:${RG_LB}" \
    "PL:${RG_PL}" \
    "SM:${SAMPLE}" \
    "PU:${RG_PU}"
do
    printf "%s\n" "$RG_LINE" | grep -F "$TAG" >/dev/null
done

grep -q "^#:GATKReport" "$RECAL"
grep -q "RecalTable0" "$RECAL"
grep -q "RecalTable1" "$RECAL"
grep -q "RecalTable2" "$RECAL"

grep -q "^## METRICS CLASS" "$DUP"
grep -F "$SAMPLE" "$DUP" >/dev/null

echo "VALIDAÇÃO ESTRUTURAL: OK"
echo "Linhas do idxstats: $IDX_LINES"
echo "$HD_LINE"
echo "$RG_LINE"
' sh \
        "$SAMPLE" \
        "$RG_ID" \
        "$RG_LB" \
        "$RG_PL" \
        "$RG_PU"

    chmod 640 \
        "$BAM_HOST" \
        "$BAI_HOST" \
        "$CHRS_HOST" \
        "$RECAL_HOST" \
        "$DUP_HOST" \
        "$SAMPLE_LOG"

    find "$QC_DIR" \
        -maxdepth 1 \
        -type f \
        -exec chmod 640 {} +

    {
        printf 'sample=%s\n' "$SAMPLE"
        printf 'patient=%s\n' "$PATIENT"
        printf 'type=%s\n' "$TYPE"
        printf 'reference=Homo_sapiens_assembly19.fasta\n'
        printf 'parabricks=4.5.1-1\n'
        printf 'gpus=4,5,6,7\n'
        printf 'threads_per_gpu=%s\n' "$THREADS_PER_GPU"
        printf 'completed=%s\n' "$(date '+%Y-%m-%d %H:%M:%S')"
    } > "$VALIDATED_MARKER"

    chmod 640 "$VALIDATED_MARKER"

    PROCESSADOS=$((PROCESSADOS + 1))

    echo
    echo "CONCLUÍDO E VALIDADO: $SAMPLE"
    echo "BAM:       $BAM_HOST"
    echo "BAI:       $BAI_HOST"
    echo "BQSR:      $RECAL_HOST"
    echo "Duplicata: $DUP_HOST"
    echo "QC:        $QC_COUNT arquivos"
    echo "Fim:       $(date '+%Y-%m-%d %H:%M:%S')"
done

chmod 640 "$RUN_LOG"

echo
echo "=========================================================="
echo "ALINHAMENTOS HG19 FINALIZADOS"
echo "Amostras selecionadas: 10"
echo "Processadas agora:     $PROCESSADOS"
echo "Já validadas:          $PULADOS"
echo "Resultados:            $RESULTS_ROOT"
echo "Log geral:             $RUN_LOG"
echo "=========================================================="
