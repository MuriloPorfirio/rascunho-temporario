#!/usr/bin/env bash

set -Eeuo pipefail
umask 027

COMPOSE_DIR="/raid/biomafia/murilo.aguiar/dockers_compose"

ALIGN_ROOT="/raid/biomafia/murilo.aguiar/FASTQS-MELANOMA/teste-piloto-5-pacientes/ALINHAMENTOS-PARABRICKS-HG19-TESTE-PILOTO"

OUTPUT_ROOT="/raid/biomafia/murilo.aguiar/FASTQS-MELANOMA/teste-piloto-5-pacientes/CHAMADA-VARIANTES-PARABRICKS-HG19-TESTE-PILOTO"

TMP_ROOT="/raid/biomafia/murilo.aguiar/TMP_parabricks/MELANOMA-PILOTO-HG19-MUTECTCALLER"

REFERENCE_ROOT="/raid/biomafia/murilo.aguiar/GENOMAS-DE-REFERENCIA-E-KNOWN-SITES/grch37_b37"
REFERENCE_DIR="${REFERENCE_ROOT}/parabricks-bwa64"
RESOURCE_DIR="${REFERENCE_ROOT}/somatic-resources"

REFERENCE_NAME="Homo_sapiens_assembly19.fasta"
REFERENCE_FILE="${REFERENCE_DIR}/${REFERENCE_NAME}"
REFERENCE_FAI="${REFERENCE_FILE}.fai"
REFERENCE_DICT="${REFERENCE_DIR}/Homo_sapiens_assembly19.dict"
EXPECTED_REFERENCE_MD5="886ba1559393f75872c1cf459eb57f2d"

GERMLINE_NAME="af-only-gnomad.raw.sites.b37-parabricks-compatible.vcf.gz"
PON_NAME="Mutect2-WGS-panel-b37.vcf.gz"

GERMLINE_FILE="${RESOURCE_DIR}/${GERMLINE_NAME}"
PON_FILE="${RESOURCE_DIR}/${PON_NAME}"

SAMTOOLS_IMAGE="staphb/samtools:1.19"
BCFTOOLS_IMAGE="staphb/bcftools:1.19"
PARABRICKS_SERVICE="parabricks-mutectcaller"

MIN_FREE_BYTES=$((500 * 1024 * 1024 * 1024))

PATIENT="${1:-}"

if [[ ! "$PATIENT" =~ ^ML_(04|05|07|09|10)$ ]]; then
    echo "USO: $0 ML_04"
    echo
    echo "PACIENTES PERMITIDOS: ML_04 ML_05 ML_07 ML_09 ML_10"
    exit 2
fi

SAMPLES_DIR="${ALIGN_ROOT}/AMOSTRAS"

echo "============================================================"
echo "DESCOBERTA AUTOMÁTICA DO PAR TUMOR-NORMAL"
echo "============================================================"

mapfile -t TUMOR_DIRS < <(
    find "$SAMPLES_DIR" \
        -mindepth 1 \
        -maxdepth 1 \
        -type d \
        -name "${PATIENT}_T_*" \
        -print |
    sort
)

mapfile -t NORMAL_DIRS < <(
    find "$SAMPLES_DIR" \
        -mindepth 1 \
        -maxdepth 1 \
        -type d \
        -name "${PATIENT}_S_*" \
        -print |
    sort
)

if [ "${#TUMOR_DIRS[@]}" -ne 1 ]; then
    echo "ERRO: era esperado exatamente um diretório tumoral para ${PATIENT}."
    echo "Encontrados: ${#TUMOR_DIRS[@]}"
    printf '%s\n' "${TUMOR_DIRS[@]:-}"
    exit 1
fi

if [ "${#NORMAL_DIRS[@]}" -ne 1 ]; then
    echo "ERRO: era esperado exatamente um diretório normal para ${PATIENT}."
    echo "Encontrados: ${#NORMAL_DIRS[@]}"
    printf '%s\n' "${NORMAL_DIRS[@]:-}"
    exit 1
fi

TUMOR_DIR="${TUMOR_DIRS[0]}"
NORMAL_DIR="${NORMAL_DIRS[0]}"

TUMOR="$(basename "$TUMOR_DIR")"
NORMAL="$(basename "$NORMAL_DIR")"

TUMOR_BAM="${TUMOR_DIR}/${TUMOR}.bam"
TUMOR_BAI="${TUMOR_BAM}.bai"
TUMOR_RECAL="${TUMOR_DIR}/${TUMOR}.recal.table"

NORMAL_BAM="${NORMAL_DIR}/${NORMAL}.bam"
NORMAL_BAI="${NORMAL_BAM}.bai"
NORMAL_RECAL="${NORMAL_DIR}/${NORMAL}.recal.table"

echo "PACIENTE: $PATIENT"
echo "TUMOR:    $TUMOR"
echo "NORMAL:   $NORMAL"

echo
echo "============================================================"
echo "VERIFICAÇÃO DOS ARQUIVOS"
echo "============================================================"

for REQUIRED_FILE in \
    "$TUMOR_BAM" \
    "$TUMOR_BAI" \
    "$TUMOR_RECAL" \
    "$NORMAL_BAM" \
    "$NORMAL_BAI" \
    "$NORMAL_RECAL" \
    "$REFERENCE_FILE" \
    "$REFERENCE_FAI" \
    "$REFERENCE_DICT" \
    "$GERMLINE_FILE" \
    "${GERMLINE_FILE}.tbi" \
    "$PON_FILE" \
    "${PON_FILE}.tbi" \
    "${PON_FILE}.pon"
do
    if [ ! -s "$REQUIRED_FILE" ]; then
        echo "ERRO: arquivo ausente ou vazio:"
        echo "$REQUIRED_FILE"
        exit 1
    fi

    echo "OK: $REQUIRED_FILE"
done

for IMAGE in "$SAMTOOLS_IMAGE" "$BCFTOOLS_IMAGE"
do
    if ! docker image inspect "$IMAGE" >/dev/null 2>&1; then
        echo "ERRO: imagem Docker ausente: $IMAGE"
        exit 1
    fi

    echo "OK: imagem $IMAGE"
done

if ! docker compose \
    --project-directory "$COMPOSE_DIR" \
    config --services 2>/dev/null |
    grep -Fxq "$PARABRICKS_SERVICE"
then
    echo "ERRO: serviço Compose não encontrado: $PARABRICKS_SERVICE"
    exit 1
fi

echo "OK: serviço Compose $PARABRICKS_SERVICE"

echo
echo "============================================================"
echo "IDENTIDADE DA REFERÊNCIA"
echo "============================================================"

REFERENCE_MD5="$(md5sum "$REFERENCE_FILE" | awk '{print $1}')"

echo "MD5 esperado: $EXPECTED_REFERENCE_MD5"
echo "MD5 obtido:   $REFERENCE_MD5"

if [ "$REFERENCE_MD5" != "$EXPECTED_REFERENCE_MD5" ]; then
    echo "ERRO: a referência não é a b37 oficial validada."
    exit 1
fi

echo "REFERÊNCIA B37 OFICIAL: OK"

samtools_header() {
    local sample_dir="$1"
    local bam_name="$2"

    docker run \
        --rm \
        --user "$(id -u):$(id -g)" \
        --network none \
        --cpus=2 \
        --memory=2g \
        --memory-swap=2g \
        --pids-limit=128 \
        --read-only \
        --tmpfs /tmp:rw,nosuid,nodev,size=256m \
        --cap-drop ALL \
        --security-opt no-new-privileges:true \
        --volume "${sample_dir}:/sample:ro" \
        --entrypoint samtools \
        "$SAMTOOLS_IMAGE" \
        view -H "/sample/${bam_name}"
}

check_bam_identity() {
    local sample_dir="$1"
    local bam_name="$2"
    local expected_sm="$3"
    local sample_label="$4"
    local observed_sms

    observed_sms="$(
        samtools_header "$sample_dir" "$bam_name" |
        awk -F '\t' '
            $1 == "@RG" {
                for (i = 2; i <= NF; i++) {
                    if ($i ~ /^SM:/) {
                        sub(/^SM:/, "", $i)
                        print $i
                    }
                }
            }
        ' |
        sort -u
    )"

    echo "$sample_label — SM observado: ${observed_sms:-AUSENTE}"

    if [ "$observed_sms" != "$expected_sm" ]; then
        echo "ERRO: SM incompatível no BAM $bam_name"
        echo "Esperado: $expected_sm"
        echo "Obtido:   ${observed_sms:-AUSENTE}"
        exit 1
    fi

    if ! diff -q \
        <(cut -f1,2 "$REFERENCE_FAI") \
        <(
            samtools_header "$sample_dir" "$bam_name" |
            awk -F '\t' '
                $1 == "@SQ" {
                    nome=""
                    comprimento=""

                    for (i = 2; i <= NF; i++) {
                        if ($i ~ /^SN:/) {
                            nome=$i
                            sub(/^SN:/, "", nome)
                        }

                        if ($i ~ /^LN:/) {
                            comprimento=$i
                            sub(/^LN:/, "", comprimento)
                        }
                    }

                    print nome "\t" comprimento
                }
            '
        ) >/dev/null
    then
        echo "ERRO: ordem, nomes ou comprimentos dos contigs do BAM são incompatíveis."
        exit 1
    fi

    echo "$sample_label — contigs compatíveis: OK"
}

echo
echo "============================================================"
echo "IDENTIDADE DOS BAMS"
echo "============================================================"

check_bam_identity "$TUMOR_DIR" "${TUMOR}.bam" "$TUMOR" "TUMOR"
check_bam_identity "$NORMAL_DIR" "${NORMAL}.bam" "$NORMAL" "NORMAL"

run_bcftools_resource() {
    docker run \
        --rm \
        --user "$(id -u):$(id -g)" \
        --network none \
        --cpus=2 \
        --memory=2g \
        --memory-swap=2g \
        --pids-limit=128 \
        --read-only \
        --tmpfs /tmp:rw,nosuid,nodev,size=256m \
        --cap-drop ALL \
        --security-opt no-new-privileges:true \
        --volume "${REFERENCE_ROOT}:/reference:ro" \
        --entrypoint bcftools \
        "$BCFTOOLS_IMAGE" \
        "$@"
}

vcf_contigs() {
    local container_vcf="$1"

    run_bcftools_resource view -h "$container_vcf" |
    awk '
        /^##contig=<ID=/ {
            linha=$0

            id=linha
            sub(/^##contig=<ID=/, "", id)
            sub(/,.*/, "", id)

            comprimento=linha
            sub(/^.*length=/, "", comprimento)
            sub(/[,>].*$/, "", comprimento)

            print id "\t" comprimento
        }
    '
}

echo
echo "============================================================"
echo "COMPATIBILIDADE DOS RECURSOS VCF"
echo "============================================================"

for RESOURCE_VCF in \
    "/reference/somatic-resources/${GERMLINE_NAME}" \
    "/reference/somatic-resources/${PON_NAME}"
do
    if ! diff -q \
        <(cut -f1,2 "$REFERENCE_FAI") \
        <(vcf_contigs "$RESOURCE_VCF") >/dev/null
    then
        echo "ERRO: cabeçalho incompatível com a referência:"
        echo "$RESOURCE_VCF"
        exit 1
    fi

    echo "ORDEM, NOMES E COMPRIMENTOS: OK — $RESOURCE_VCF"
done

echo
echo "============================================================"
echo "ESPAÇO DISPONÍVEL"
echo "============================================================"

mkdir -p -- "$OUTPUT_ROOT" "$TMP_ROOT"

OUTPUT_AVAILABLE="$(df --output=avail -B1 "$OUTPUT_ROOT" | tail -n 1 | tr -d ' ')"
TMP_AVAILABLE="$(df --output=avail -B1 "$TMP_ROOT" | tail -n 1 | tr -d ' ')"

df -h "$OUTPUT_ROOT" "$TMP_ROOT"

if [ "$OUTPUT_AVAILABLE" -lt "$MIN_FREE_BYTES" ]; then
    echo "ERRO: menos de 500 GiB disponíveis no destino das saídas."
    exit 1
fi

if [ "$TMP_AVAILABLE" -lt "$MIN_FREE_BYTES" ]; then
    echo "ERRO: menos de 500 GiB disponíveis no destino temporário."
    exit 1
fi

echo "ESPAÇO MÍNIMO DE 500 GiB: OK"

echo
echo "============================================================"
echo "SITUAÇÃO DAS GPUs FÍSICAS 6 E 7"
echo "============================================================"

nvidia-smi \
    -i 6,7 \
    --query-gpu=index,name,memory.used,utilization.gpu \
    --format=csv

GPU_PROCESSES="$(
    nvidia-smi \
        -i 6,7 \
        --query-compute-apps=pid,process_name,used_gpu_memory \
        --format=csv,noheader,nounits 2>/dev/null |
    sed '/^[[:space:]]*$/d' || true
)"

if [ -n "$GPU_PROCESSES" ]; then
    echo
    echo "ERRO: há processos utilizando as GPUs físicas 6 ou 7:"
    echo "$GPU_PROCESSES"
    exit 1
fi

echo "GPUs físicas 6 e 7 sem processos: OK"

RUN_ID="$(date +'%Y%m%d-%H%M%S')"
RUN_DIR_NAME="WGS-${PATIENT}_${RUN_ID}"

OUTPUT_RELATIVE="AMOSTRAS/${PATIENT}/${RUN_DIR_NAME}"
TMP_RELATIVE="AMOSTRAS/${PATIENT}/${RUN_DIR_NAME}"

OUTPUT_DIR="${OUTPUT_ROOT}/${OUTPUT_RELATIVE}"
TMP_DIR="${TMP_ROOT}/${TMP_RELATIVE}"

if [ -e "$OUTPUT_DIR" ] || [ -e "$TMP_DIR" ]; then
    echo "ERRO: diretório da corrida já existe. Nada será sobrescrito."
    echo "$OUTPUT_DIR"
    echo "$TMP_DIR"
    exit 1
fi

mkdir -p -- "$OUTPUT_DIR" "$TMP_DIR"
chmod 750 -- "$OUTPUT_DIR" "$TMP_DIR"

PREFIX="${TUMOR}_vs_${NORMAL}"
VCF_NAME="${PREFIX}.mutectcaller.unfiltered.vcf.gz"
F1R2_NAME="${PREFIX}.f1r2.tar.gz"
LOG_NAME="mutectcaller_wgs.log"
CHECKSUM_NAME="SHA256SUMS.txt"
MANIFEST_NAME="MANIFESTO-DA-CORRIDA.txt"

VCF_FILE="${OUTPUT_DIR}/${VCF_NAME}"
VCF_INDEX="${VCF_FILE}.tbi"
STATS_FILE="${VCF_FILE}.stats"
F1R2_FILE="${OUTPUT_DIR}/${F1R2_NAME}"
LOG_FILE="${OUTPUT_DIR}/${LOG_NAME}"

echo
echo "============================================================"
echo "PLANO DA CORRIDA WGS"
echo "============================================================"
echo "PACIENTE:       $PATIENT"
echo "TUMOR:          $TUMOR"
echo "NORMAL:         $NORMAL"
echo "SAÍDA:          $OUTPUT_DIR"
echo "TEMPORÁRIO:     $TMP_DIR"
echo "VCF BRUTO:      $VCF_NAME"
echo "F1R2:           $F1R2_NAME"
echo "GPUs FÍSICAS:   6 e 7"
echo "GPUs CONTAINER: 0 e 1"
echo "THREADS:        24"

echo
echo "============================================================"
echo "EXECUTANDO MUTECTCALLER WGS"
echo "============================================================"

cd "$COMPOSE_DIR"

TIMEFORMAT=$'\nTEMPO TOTAL DO MUTECTCALLER WGS: %3R segundos'

{
    time docker compose run \
        --rm \
        --no-deps \
        --no-TTY \
        "$PARABRICKS_SERVICE" \
        mutectcaller \
        --ref "/reference/parabricks-bwa64/${REFERENCE_NAME}" \
        --in-tumor-bam "/input/AMOSTRAS/${TUMOR}/${TUMOR}.bam" \
        --in-normal-bam "/input/AMOSTRAS/${NORMAL}/${NORMAL}.bam" \
        --in-tumor-recal-file "/input/AMOSTRAS/${TUMOR}/${TUMOR}.recal.table" \
        --in-normal-recal-file "/input/AMOSTRAS/${NORMAL}/${NORMAL}.recal.table" \
        --tumor-name "$TUMOR" \
        --normal-name "$NORMAL" \
        --mutect-germline-resource "/reference/somatic-resources/${GERMLINE_NAME}" \
        --pon "/reference/somatic-resources/${PON_NAME}" \
        --out-vcf "/output/${OUTPUT_RELATIVE}/${VCF_NAME}" \
        --mutect-f1r2-tar-gz "/output/${OUTPUT_RELATIVE}/${F1R2_NAME}" \
        --num-gpus 2 \
        --num-htvc-threads 24 \
        --tmp-dir "/parabricks-tmp/${TMP_RELATIVE}"
} 2>&1 |
tee "$LOG_FILE"

echo
echo "============================================================"
echo "VALIDAÇÃO DAS SAÍDAS"
echo "============================================================"

for PRODUCED_FILE in \
    "$VCF_FILE" \
    "$VCF_INDEX" \
    "$STATS_FILE" \
    "$F1R2_FILE" \
    "$LOG_FILE"
do
    if [ ! -s "$PRODUCED_FILE" ]; then
        echo "ERRO: saída ausente ou vazia:"
        echo "$PRODUCED_FILE"
        exit 1
    fi

    echo "OK: $PRODUCED_FILE"
done

docker run \
    --rm \
    --user "$(id -u):$(id -g)" \
    --network none \
    --cpus=2 \
    --memory=2g \
    --memory-swap=2g \
    --pids-limit=128 \
    --read-only \
    --tmpfs /tmp:rw,nosuid,nodev,size=256m \
    --cap-drop ALL \
    --security-opt no-new-privileges:true \
    --volume "${OUTPUT_DIR}:/data:ro" \
    --entrypoint bcftools \
    "$BCFTOOLS_IMAGE" \
    view -h "/data/${VCF_NAME}" >/dev/null

VARIANT_COUNT="$(
    docker run \
        --rm \
        --user "$(id -u):$(id -g)" \
        --network none \
        --cpus=2 \
        --memory=2g \
        --memory-swap=2g \
        --pids-limit=128 \
        --read-only \
        --tmpfs /tmp:rw,nosuid,nodev,size=256m \
        --cap-drop ALL \
        --security-opt no-new-privileges:true \
        --volume "${OUTPUT_DIR}:/data:ro" \
        --entrypoint bcftools \
        "$BCFTOOLS_IMAGE" \
        view -H "/data/${VCF_NAME}" |
    wc -l
)"

VCF_SAMPLES="$(
    docker run \
        --rm \
        --user "$(id -u):$(id -g)" \
        --network none \
        --cpus=2 \
        --memory=2g \
        --memory-swap=2g \
        --pids-limit=128 \
        --read-only \
        --tmpfs /tmp:rw,nosuid,nodev,size=256m \
        --cap-drop ALL \
        --security-opt no-new-privileges:true \
        --volume "${OUTPUT_DIR}:/data:ro" \
        --entrypoint bcftools \
        "$BCFTOOLS_IMAGE" \
        query -l "/data/${VCF_NAME}" |
    paste -sd ',' -
)"

if ! tar -tzf "$F1R2_FILE" >/dev/null; then
    echo "ERRO: o arquivo F1R2 não é um tar.gz legível."
    exit 1
fi

echo "LEITURA DO VCF: OK"
echo "AMOSTRAS NO VCF: $VCF_SAMPLES"
echo "VARIANTES BRUTAS: $VARIANT_COUNT"
echo "LEITURA DO F1R2: OK"

echo
echo "============================================================"
echo "CHECKSUMS E MANIFESTO"
echo "============================================================"

(
    cd "$OUTPUT_DIR"

    sha256sum \
        "$VCF_NAME" \
        "${VCF_NAME}.tbi" \
        "${VCF_NAME}.stats" \
        "$F1R2_NAME" \
        "$LOG_NAME" \
        > "$CHECKSUM_NAME"
)

{
    printf 'PACIENTE\t%s\n' "$PATIENT"
    printf 'TUMOR\t%s\n' "$TUMOR"
    printf 'NORMAL\t%s\n' "$NORMAL"
    printf 'RUN_ID\t%s\n' "$RUN_ID"
    printf 'PARABRICKS\t4.5.1-1\n'
    printf 'REFERENCIA\t%s\n' "$REFERENCE_FILE"
    printf 'REFERENCIA_MD5\t%s\n' "$REFERENCE_MD5"
    printf 'GERMLINE_RESOURCE\t%s\n' "$GERMLINE_FILE"
    printf 'PON\t%s\n' "$PON_FILE"
    printf 'VCF\t%s\n' "$VCF_FILE"
    printf 'F1R2\t%s\n' "$F1R2_FILE"
    printf 'VARIANTES_BRUTAS\t%s\n' "$VARIANT_COUNT"
    printf 'AMOSTRAS_NO_VCF\t%s\n' "$VCF_SAMPLES"
    printf 'GPUS_FISICAS\t6,7\n'
    printf 'NUM_GPUS\t2\n'
    printf 'THREADS\t24\n'
} > "${OUTPUT_DIR}/${MANIFEST_NAME}"

echo "OK: ${OUTPUT_DIR}/${CHECKSUM_NAME}"
echo "OK: ${OUTPUT_DIR}/${MANIFEST_NAME}"

echo
echo "ARQUIVOS PRODUZIDOS:"

find "$OUTPUT_DIR" \
    -maxdepth 1 \
    -type f \
    -printf '%f\t%s bytes\n' |
sort

echo
echo "============================================================"
echo "MUTECTCALLER WGS CONCLUÍDO COM SUCESSO — $PATIENT"
echo "============================================================"
echo "SAÍDA: $OUTPUT_DIR"
echo
echo "O VCF É BRUTO. NÃO É UM VCF FINAL PARA INTERPRETAÇÃO."
echo "O F1R2 FOI PRESERVADO PARA REFINAMENTO FUTURO."
