(
    set -euo pipefail
    umask 007

    test "$(findmnt -rn --mountpoint /mnt/apollo -o FSTYPE)" = "nfs4"

    hg38_repo='/raid/biomafia/murilo.aguiar/BRISANCE-pipeline'
    hg38_projeto='/mnt/apollo/murilo.aguiar/MELANOMA-PILOTO-HG38'

    bash "$hg38_repo/interface/run_brisance.sh" \
        "$hg38_repo/ML_04_HG38.local.ini" \
        --run 2>&1 |
        tee "$hg38_projeto/controle/execucao_ML04_$(date -u +%Y%m%dT%H%M%SZ).log"
)
