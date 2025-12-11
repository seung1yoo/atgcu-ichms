# atgcu-ichms

Markers 추출 및 LIS API 전송 파이프라인 CLI.

## 주요 기능
- `extract`: APT/Imputed VCF에서 마커를 추출해 OA(txt)와 JSON을 생성
- No-call을 Imputed VCF로 보완, 없으면 NORMAL_ALLELE로 채움
- 샘플 시트(`--sample-sheet`) 사용: SAMPLE_NAME ↔ SAMPLE_ID_IN_A_PLATE 매핑 후 출력 ID로 사용
- 샘플/마커 순서 유지 정렬, allele 정규화(`N`, `*`, `-`, `None` → `.`)
- 상세 QC 리포트 + RS-level QC TSV 생성
- `api-to-lis`: 생성된 JSON을 conf.json에 지정된 API로 POST(타임아웃/재시도 지원)

## 요구사항
- Python 3.8+
- pandas, requests
- pysam (선택: imputed 모드 `tabix/auto`에서 tabix 인덱스 사용 시 필요)

## 설정 예시 (`conf.json`)
```json
{
  "dtc_62": {
    "apt_vcf": "data/Analysis_For_All/08.make_Genotype_format/GenotypeCall_for_imput.vcf.gz",
    "imputed_vcf": "data/Analysis_For_All/15.final_vcf/imputed_all.vcf.gz",
    "rsmarker_list": "src/dtc_62_rs_markers.GRCh38.tsv",
    "output_dir": "extract_dtc_62",
    "qc_report": "extract_dtc_62/qc.report.txt"
  }
}
```

## 사용법
진입점: `python -m atgcu_ichms.cli ...`

### 1) 추출 (OA + JSON 생성)
```bash
python -m atgcu_ichms.cli extract \
  --product dtc_62 \
  --sample-sheet data/Analysis_For_All.sample_sheet.tsv \
  --path-anl-for-all /abs/path/Analysis_For_All \
  --path-wkdir /abs/path/workdir \
  --imputed-mode auto \
  --conf-json conf.json \
  --log-level INFO
```
- `--path-anl-for-all`: conf의 `apt_vcf`, `imputed_vcf`가 상대경로일 때 붙여 절대경로로 사용
- `--path-wkdir`: `output_dir`, `qc_report`를 이 경로 하위로 생성
- `--imputed-mode`: `auto`(기본: tabix 시도 후 미스 시 stream), `tabix`(탐색 실패해도 stream 미수행), `stream`(tabix 생략)
- 출력:
  - OA: `<output_dir>/<product>.oa.txt`
  - JSON: `<output_dir>/<product>.oa.json`
  - QC 리포트: `<output_dir>/<product>.oa.txt.qc.report.txt`
  - RS-level QC TSV: `<output_dir>/<product>.qc_rs.tsv`

### 2) JSON을 LIS API로 전송
```bash
python -m atgcu_ichms.cli api-to-lis \
  --product dtc_62 \
  --json extract_dtc_62/dtc_62.oa.json \
  --conf-json conf.json \
  --timeout 10 \
  --retries 3
```
- `conf.json`의 `product.api-url` 또는 별도 API 설정을 사용

## QC 리포트/로그 주요 항목
- 실행 인자: product, sample-sheet, path-anl-for-all, path-wkdir, imputed-mode
- Imputed 조회 상세: tabix/stream 히트 수, 고유 위치 수, indexed_variant_count
- RS 상태 요약: IMPUTED/NORMAL_FILLED/APT/MISSING 개수
- RS-level TSV(`*.qc_rs.tsv`): RS별 APT 존재 여부, 좌표, imputed tabix/stream 히트 여부, 최종 소스(IMPUTED/NORMAL_FILLED/APT/MISSING)

## 참고
- Imputed VCF는 bgzip + tabix 인덱스가 있으면 성능 향상
  - 예: `bgzip imputed_all.vcf.gz` 후 `tabix -p vcf imputed_all.vcf.gz`
  - `.tbi`/`.csi`가 있으면 `pysam`으로 랜덤 액세스
- `pysam` 미설치 시 `--imputed-mode tabix/auto`에서 tabix를 쓸 수 없으며, `auto`는 stream fallback, `tabix`는 실패할 수 있음
- 마커 TSV에 `CHROM`, `POS`가 비어 있으면 `utils/rsid_coordi_fetch.py`로 좌표 채우기
  - 예:
    ```bash
    python utils/rsid_coordi_fetch.py \
      --input src/dtc_53_rs_markers.txt \
      --output src/dtc_53_rs_markers.filled.txt \
      --assembly GRCh38 \
      --sleep 0.2
    ```
  - 입력 TSV 컬럼: `RS_ID, NORMAL_ALLELE, RISK_ALLELE, CHROM, POS`
  - CHROM/POS가 비어 있는 RS만 dbSNP API에서 조회해 채움

